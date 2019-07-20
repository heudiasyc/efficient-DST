#ifndef OW_BFT_MASS_AGGREGATE_HPP
#define OW_BFT_MASS_AGGREGATE_HPP

#include <aggregate.hpp>
#include <mass.hpp>

namespace ow_bft{

	/*
	 * mass_aggregate is a bft_function that results in the sum of images of a mass function
	 * i.e. commonality, belief, implicability, plausibility, conjunctive weights and disjunctive weights
	 */
	template <typename T = double>
	class mass_aggregate : public aggregate<T> {
	protected:

		const mass<T> mass_equivalent;

		const FOD *fod;
		/*
		 * special_elements are the only elements (and images through this aggregate) needed to allow for the reconstruction
		 * of the mass images of focal elements (could be only focal elements)
		 */
		powerset_btree<T> special_elements;

		void display_message_no_mass_function_provided(){
			std::clog << "No mass function provided for this mass function aggregate." << std::endl
					<< "Initializing with vacuous case." << std::endl;
		}

	public:

		mass_aggregate(const mass<T>& m) :
			mass_equivalent(m),
			fod(&(mass_equivalent.get_FOD())),
			special_elements(*fod, this->block_size)
		{
			if(!this->mass_equivalent.is_valid()){
				std::cerr << "\nError: The sum of all images of the mass function provided to this mass aggregate isn't equal to 1.";
				exit(1);
			}
		}

		mass_aggregate(const powerset_btree<T>& m_focal_elements) :
			mass_equivalent(m_focal_elements),
			fod(&(mass_equivalent.get_FOD())),
			special_elements(*fod, this->block_size)
		{
			if(!this->mass_equivalent.is_valid()){
				std::cerr << "\nError: The sum of all images of the mass function provided to this mass aggregate isn't equal to 1.";
				exit(1);
			}
		}

		mass_aggregate(const powerset_btree<T>& m_focal_elements, const powerset_btree<T>& _special_elements) : // @suppress("Class members should be properly initialized")
			mass_aggregate(m_focal_elements)
		{
			this->special_elements.copy(_special_elements);
		}

		mass_aggregate(const FOD& fod) : mass_aggregate(fod, vacuous)
		{
			display_message_no_mass_function_provided();
		}

		mass_aggregate(const FOD& _fod, const Special_case s_case) :
			mass_equivalent(_fod, s_case),
			fod(&(mass_equivalent.get_FOD())),
			special_elements(*fod, this->block_size)
		{}

		virtual ~mass_aggregate(){}


		const FOD& get_FOD() const {
			return *(this->fod);
		}

		const powerset_btree<T>& get_special_elements() const {
			return this->special_elements;
		}

		const mass<T>& get_mass_equivalent() const {
			return this->mass_equivalent;
		}

		T at_emptyset() const {
			set_N_value<T>* set_value = this->special_elements.sub_fod_of_size(0);
			if(set_value){
				return set_value->value;
			}
			return this->compute_aggregation_at_emptyset();
		}

		T at_fod() const {
			set_N_value<T>* set_value = this->special_elements.sub_fod_of_size(this->fod->size());
			if(set_value){
				return set_value->value;
			}
			return this->compute_aggregation_at_fod();
		}

		T operator[](const std::vector<std::string>& labels) const {
			return find(this->fod->to_set(labels));
		}

		T find(const boost::dynamic_bitset<>& set) const {
			set_N_value<T>* set_value = this->special_elements[set];
			if(set_value){
				return set_value->value;
			}
			return this->compute_aggregation(set);
		}

/*
		void set_special_elements(const powerset_btree<T>& f_elements, bool update_mass_equivalent=true) {
			aggregate<T>::set_special_elements(f_elements);

			if(update_mass_equivalent){
				set_mass_equivalent(to_mass(f_elements), false);
				// check if the new mass function is valid (i.e. if the sum of all its images is 1)
				if(!this->mass_equivalent.is_valid()){
					std::cerr << "\nInvalid BBA given : once converted to a mass function, the sum of its images isn't equal to 1 or some image is negative.\n";
					exit(1);
				}
			}
		}
*/
/*
		void set_mass_equivalent(const mass<T>& m, bool update_focal_elements=true){
			this->mass_equivalent = m;

			if(update_focal_elements){
				//this->fod->erase_powerset(this->focal_elements);
				//this->fod = &(this->mass_equivalent.get_FOD());
				//this->fod->push_back_powerset(this->focal_elements);
				if(&(this->mass_equivalent.get_FOD()) != this->fod){
					std::cerr << "\nCan't replace mass equivalent with a mass function using another FOD.\n";
					return;
				}
				this->special_elements.nullify();
				set_values_for_special_elements();
			}
		}
*/

/*
	public:
		bool has_been_changed(){
			if(this->changes != nullptr)
				return this->changes.size() != 0;
			else
				return false;
		}

		bft_function(bft_function<T>& f){
			if(f.origin == nullptr){
				if(f.type == mass_t && this->type == mass_t){
					this->fod = f.fod;
					this->origin = f.origin;
					this->before_changes = f.before_changes;
					this->changes = f.changes;
				}else{
					std::cerr << "Function created ignoring predecessor in parameter :"
							  << "If this function is not of type mass, it must have a mass function as predecessor"
							  << std::endl;
				}
			}else {
				bft_function(f.fod);
				if(f.origin.type == mass_t){
					// example of possible case : this->type = commonality, old_f.type = mass
					this->origin = f;
					this->before_changes = ow_bft::powerset_btree<T>(&this->fod);
					this->fod.push_back_powerset(&this->before_changes);
					std::vector<set_N_value<T>* > elements = this->origin.changes.elements();
					// insertion of value of type this->type at each focal element
					for (int i = 0; i < elements.size(); ++i) {
						this->before_changes.insert(elements[i]->fod_elements, original_accessor(elements[i]->fod_elements));
					}
				}else{
					std::cerr << "Function created ignoring predecessor in parameter :"
							  << "The predecessor must be a mass function"
							  << std::endl;
				}
			}

			////////////////////////////
			bool return_to_original_representation = false;
			bft_function<T> old_f = f;
			while(!return_to_original_representation && !old_f.has_been_changed() && old_f.origin != nullptr){
				if(old_f.origin.type == this->type){
					return_to_original_representation = true;
				}
				old_f = old_f.origin;
			}
			if((old_f.origin == nullptr && !old_f.has_been_changed()) || return_to_original_representation){
				this->fod = old_f.fod;
				this->origin = old_f.origin;
				this->before_changes = old_f.before_changes;
				this->changes = old_f.changes;
			}else{
				// allows for the change of bft_function interface and for modifications in one and only one other space
				if(old_f.origin != nullptr || old_f.type != mass_t){
					std::cerr << "Cannot have more than one origin and that origin must be a mass function" << std::endl
							  << "example : this->type = commonality, old_f.type = mass"
							  << std::endl;
					exit(EXIT_FAILURE);
				}
				// example of possible case : this->type = commonality, old_f.type = mass
				bft_function(old_f.fod);
				this->origin = old_f;
				this->before_changes = ow_bft::powerset_btree<T>(&this->fod);
				this->fod.push_back_powerset(&this->before_changes);
			}
		}


		// =============================================================================

		void nullify(std::vector<std::string> labels){
			std::vector<fod_element*> fod_elements = this->fod.to_elements(labels);

			if(this->mass_equivalent != nullptr){
				this->before_changes.nullify(this->before_changes[fod_elements]);
			}
			this->changes.nullify(this->changes[fod_elements]);
		}

		void set_values(std::unordered_map<std::vector<std::string>, T> values){
			for (std::pair<std::vector<std::string>, T> labels_U_value : values){
				set_value(labels_U_value.first, labels_U_value.second);
			}
		}

		void set_value(std::vector<std::string> labels, T value){
			std::vector<fod_element*> fod_elements = this->fod.to_elements(labels);

			if(value < 0)
				value = 0;

			if(this->origin == nullptr){
				this->changes.insert(fod_elements, value);
				return;
			}

			set_N_value<T>* s_N_v = this->before_changes[fod_elements];

			if(s_N_v != nullptr){	// i.e. if s_N_v exists in before_changes and has a non null value
				// if this set value has already be changed, check if that value isn't canceling the previous change
				if(s_N_v->value != value){
					this->changes.insert(fod_elements, value);
				}else{
					// if so, nullify the previous change
					this->before_changes.nullify(s_N_v);
					this->changes.nullify(this->changes[fod_elements]);
				}
			}else{
				// if this set value hadn't been changed before, check if this value is different from the original value
				T old_value = original_accessor(fod_elements);

				if(value != old_value){
					this->before_changes.insert(fod_elements, old_value);
					this->changes.insert(fod_elements, value);
				}
			}
		}

		void set_emptyset_value(T value){
			if(value < 0)
				value = 0;

			if(this->origin == nullptr){
				this->changes.set_value_of_sub_fod_of_size(0, value);
				return;
			}

			set_N_value<T>* s_N_v = this->before_changes.sub_fod_of_size(0);

			if(s_N_v != nullptr){	// i.e. if s_N_v exists in before_changes and has a non null value
				// if this set value has already be changed, check if that value isn't canceling the previous change
				if(s_N_v->value != value){
					this->changes.set_value_of_sub_fod_of_size(0, value);
				}else{
					// if so, nullify the previous change
					this->before_changes.nullify(s_N_v);
					this->changes.nullify(this->changes.sub_fod_of_size(0));
				}
			}else{
				// if this set value hadn't been changed before, check if this value is different from the original value
				T old_value = original_at_emptyset();

				if(value != old_value){
					this->before_changes.set_value_of_sub_fod_of_size(0, value);
					this->changes.set_value_of_sub_fod_of_size(0, value);
				}
			}
		}

		void set_fod_value(T value){
			if(value < 0)
				value = 0;

			if(this->origin == nullptr){
				this->changes.set_value_of_sub_fod_of_size(this->fod.size(), value);
				return;
			}

			set_N_value<T>* s_N_v = this->before_changes.sub_fod_of_size(this->fod.size());

			if(s_N_v != nullptr){	// i.e. if s_N_v exists in before_changes and has a non null value
				// if this set value has already be changed, check if that value isn't canceling the previous change
				if(s_N_v->value != value){
					this->changes.set_value_of_sub_fod_of_size(this->fod.size(), value);
				}else{
					// if so, nullify the previous change
					this->before_changes.nullify(s_N_v);
					this->changes.nullify(this->changes.sub_fod_of_size(this->fod.size()));
				}
			}else{
				// if this set value hadn't been changed before, check if this value is different from the original value
				T old_value = original_at_fod();

				if(value != old_value){
					this->before_changes.set_value_of_sub_fod_of_size(this->fod.size(), value);
					this->changes.set_value_of_sub_fod_of_size(this->fod.size(), value);
				}
			}
		}
*/
	};
}		// namespace ow_bft

#endif // OW_BFT_MASS_AGGREGATE_HPP
