#ifndef OW_BFT_MASS_HPP
#define OW_BFT_MASS_HPP

#include <bft_function.hpp>
#include <powerset_function.hpp>

namespace ow_bft{

	template <typename T = double>
	class mass : public bft_function<T>{
	protected:
		FOD fod;
		/*
		 * only focal elements (as defined in the mass space) and respective images
		 */
		powerset_btree<T> focal_elements;

	public:

		mass(const mass<T>& m) : mass(m.get_focal_elements())
		{}

		mass(const powerset_btree<T>& _focal_elements) : mass(*_focal_elements.fod)
		{
			this->focal_elements.copy(_focal_elements);
		}

		mass() :
			focal_elements(fod, this->block_size)
		{}

		mass(const FOD& _fod) :
			fod(_fod),
			focal_elements(fod, this->block_size)
		{
			//this->fod->push_back_powerset_function(*this);
		}

		mass(const FOD& fod, const Special_case s_case) : mass(fod)
		{
			switch(s_case){
				// create a mass function with all mass attributed to the empty set
			    case degenerate  : set_emptyset_value(1);	break;
				// create a mass function with all mass attributed to the FOD set
			    case vacuous	:  set_fod_value(1);	break;
			}
		}

		~mass(){
			//this->fod->erase_powerset_function(*this);
		}

		const FOD& get_FOD() const {
			return this->fod;
		}

		const powerset_btree<T>& get_focal_elements() const {
			return this->focal_elements;
		}

		void set_focal_elements(const powerset_btree<T>& f_elements){
			/*
			if(&(f_elements.get_FOD()) != &(this->fod)){
				std::cerr << "\nCan't replace focal elements with elements using another FOD.\n";
				std::cerr << "You have to create another mass function.\n";
				return;
			}*/
			this->focal_elements.nullify();
			this->focal_elements.copy(f_elements);
		}

		void erase_elements_containing_fod_element(const std::string& element_label){
			size_t position = this->fod.to_element(element_label)->position_in_fod;
			this->focal_elements.erase_elements_containing_fod_element(position);
			this->fod.erase(position);
			//powerset_function::erase_elements_containing_fod_element(position);
		}

		template <class fusion_rule>
		mass<T> apply(const fusion_rule fusion, const mass<T>& m2) const {
			return fusion(*this, m2);
		}

		void nullify(const std::vector<std::string>& labels) {
			this->focal_elements.nullify(this->focal_elements[labels]);
		}

		void set_values(const std::unordered_map<std::vector<std::string>, T>& values) {
			for (std::pair<std::vector<std::string>, T> labels_U_value : values){
				set_value(labels_U_value.first, labels_U_value.second);
			}
		}

		void set_value(const std::vector<std::string>& labels, T value) {
			this->focal_elements.insert(labels, value);
		}

		void set_emptyset_value(const T& value) {
			this->focal_elements.set_value_of_sub_fod_of_size(0, value);
		}

		void set_fod_value(const T& value) {
			this->focal_elements.set_value_of_sub_fod_of_size(this->fod.size(), value);
		}

		T at_emptyset() const {
			set_N_value<T>* set_value = this->focal_elements.sub_fod_of_size(0);
			if(set_value != nullptr)
				return set_value->value;
			else
				return 0;
		}

		T at_fod() const {
			set_N_value<T>* set_value = this->focal_elements.sub_fod_of_size(this->fod.size());
			if(set_value != nullptr)
				return set_value->value;
			else
				return 0;
		}

		T operator[](const std::vector<std::string>& labels) const {
			return find(this->fod.to_set(labels));
		}

		T find(const boost::dynamic_bitset<>& key) const {
			set_N_value<T>* set_value = this->focal_elements[key];
			if(set_value != nullptr)
				return set_value->value;
			else
				return 0;
		}

		///////////////////////////////////////////////////////////////////

		/// Check if this mass function is valid (i.e. if the sum of all its images is 1)
		bool is_valid() const {
			const std::vector<set_N_value<T>* >& f_elements = this->focal_elements.elements();
			T sum = 0;
			for (size_t i = 0; i < f_elements.size(); ++i) {
				if(f_elements[i]->value < 0 && !this->is_equivalent_to_zero(f_elements[i]->value))
					return false;
				sum += f_elements[i]->value;
			}

			return this->is_equivalent_to_zero(1-sum);
		}

		/// A focal set has a non-zero mass attributed.
		bool is_focal(const boost::dynamic_bitset<>& set) const {
			return !this->is_equivalent_to_zero(find(set));
		}

		bool is_focal(const std::vector<std::string>& labels) const {
			return !this->is_equivalent_to_zero(this->operator[](labels));
		}

		/// Normal mass function has no mass at empty set (Conflict).
		bool is_normal() const {
			return this->is_equivalent_to_zero(this->at_emptyset());
		}

		/// Regular mass function has no mass at empty set (Conflict).
		bool is_regular() const {
			return is_normal();
		}

		/// Subnormal mass function has some mass at empty set (Conflict).
		bool is_subnormal() const {
			return !is_normal();
		}

		/// Dogmatic mass function has no mass at FOD.
		bool is_dogmatic() const {
			return this->is_equivalent_to_zero(this->at_fod());
		}

		/// Non-dogmatic mass function has some mass at FOD.
		bool is_nondogmatic() const {
			return !is_dogmatic();
		}

		/// Vacuous mass function has all mass at FOD.
		bool is_vacuous() const {
			return this->is_equivalent_to_zero(1 - this->at_fod());
		}

		/// Degenerate mass function has all mass at empty set (Conflict).
		bool is_degenerate() const {
			return this->is_equivalent_to_zero(1 - this->at_emptyset());
		}

		/// Categorical mass function has only one focal set.
		bool is_categorical() const {
			return this->focal_elements.elements().size() == 1;
		}

		/// Simple mass function has at most two focal sets, and if it has two,
		/// Omega is one of them.
		bool is_simple() const {
			const std::vector<set_N_value<T>* >& elements = this->focal_elements.elements();

			if(elements.size() == 2){
				return at_fod() != nullptr;
			}else if(elements.size() < 2){
				return true;
			}else{
				return false;
			}
		}

		/// In bayesian mass function, all focal sets are singletons.
		bool is_bayesian() const {
			T belief = 0;
			const std::vector<set_N_value<T>* >& singletons = this->focal_elements.singletons();

			for (size_t i = 0; i < singletons.size(); ++i) {
				belief += singletons[i]->value;
			}
			return this->is_equivalent_to_zero(1-belief);
		}

		/// Mass function is consonant when all focal sets are nested.
		///
		/// Consonant evidence can be represented as a nested structure of subsets
		/// where the
		/// elements of the smallest set are included in the next larger set, all of
		/// whose elements are
		/// included in the next larger set and so on.
		bool is_consonant() const {

			const std::vector<std::vector<set_N_value<T>* > >& elements_by_cardinality = this->focal_elements.elements_by_set_cardinality();

			for (size_t i = 0; i < elements_by_cardinality.size(); ++i) {
				if(elements_by_cardinality[i].size() > 1)
					// if there is more than one element of same cardinality, the structure isn't nested
					return false;
			}
			// at this point, the structure has at most one element per cardinality

			// get smallest set
			size_t i = 0;
			while(i < elements_by_cardinality.size() && elements_by_cardinality[i].size() == 0){
				++i;
			}
			// get next smallest set
			size_t j = i + 1;
			while(j < elements_by_cardinality.size() && elements_by_cardinality[j].size() == 0){
				++j;
			}
			while (j < elements_by_cardinality.size()) {
				if (!this->fod.is_or_is_subset_of(
						elements_by_cardinality[i][0]->fod_elements,
						elements_by_cardinality[j][0]->fod_elements)) {
					return false;
				}
				i = j;
				// get next smallest set
				++j;
				while(j < elements_by_cardinality.size() && elements_by_cardinality[j].size() == 0){
					++j;
				}
			}
			return true;
		}

		/// Mass function is consistent when there is at least one set that is
		/// common to all
		/// focal sets.
		bool is_consistent() const {

			const std::vector<set_N_value<T>* >& elements_by_cardinality = this->focal_elements.elements_by_set_cardinality();
			size_t i = 0;
			// get smallest cardinality
			while(i < elements_by_cardinality.size() && elements_by_cardinality[i].size() == 0){
				++i;
			}
			if(i == elements_by_cardinality.size())
				// if this is true, then there is no focal set...
				// ...so it's inconsistent but in another way
				return false;
			if(elements_by_cardinality[i].size() > 1)
				// if there is more than one element of smallest cardinality, it is
				// inconsistent since two sets of same cardinality cannot contain each other
				return false;

			// at this point, there is only one set that has the smallest cardinality

			for (size_t j = i + 1; j < elements_by_cardinality.size(); ++j) {
				for (size_t k = 0; k < elements_by_cardinality.size(); ++k) {
					if (!this->fod.is_or_is_subset_of(
							elements_by_cardinality[i][0]->fod_elements,
							elements_by_cardinality[j][k]->fod_elements)) {
						return false;
					}
				}
			}
			return true;
		}

		/// Arbitrary evidence corresponds to the situation where there is no
		/// element common to
		/// all focal sets, though some sets may have elements in common.
		/// In other words, the intersection of all focal elements is empty
		bool is_arbitrary() const {

			const std::vector<set_N_value<T>* >& elements = this->focal_elements.elements();
			const boost::dynamic_bitset<>& seti = this->fod.to_set(elements[0]->fod_elements);

			for (size_t j = 1; j < elements.size(); ++j) {
				const boost::dynamic_bitset<>& setj = this->fod.to_set(elements[j]->fod_elements);
				seti = this->fod.set_intersection(seti, setj);
				if (this->fod.is_emptyset(seti)){
					// if the intersection of the j+1 first focal elements is already empty
					// then the intersection with the rest of the focal elements will be empty
					return true;
				}
			}
			return false;
		}

		/// Disjoint evidence implies that any two focal elements have no elements in
		/// common with any
		/// other focal element.
		/// In other words, the intersection of each focal element with the union of other elements is empty
		bool is_disjoint() const {

			const std::vector<set_N_value<T>* >& elements = this->focal_elements.elements();
			const boost::dynamic_bitset<>& focals_union = this->fod.to_set(elements[0]->fod_elements);

			for (size_t j = 1; j < elements.size(); ++j) {
				const boost::dynamic_bitset<>& setj = this->fod.to_set(elements[j]->fod_elements);
				if (!this->fod.are_disjoint(setj, focals_union)) {
					return false;
				}
				focals_union = this->fod.set_union(setj, focals_union);
			}
			return true;
		}

		/// Focal sets of a partitioned mass function do not intersect and their
		/// union is equal to the frame of discernment \f$ \Omega \f$.
		/// In other words, the union of all focal elements is omega,
		/// and the intersection of each focal element with the union of other elements is empty
		bool is_partitioned() const {
			const std::vector<set_N_value<T>* >& elements = this->focal_elements.elements();
			const boost::dynamic_bitset<>& focals_union = this->fod.to_set(elements[0]->fod_elements);

			for (size_t j = 1; j < elements.size(); ++j) {
				const boost::dynamic_bitset<>& setj = this->fod.to_set(elements[j]->fod_elements);
				if (!this->fod.are_disjoint(setj, focals_union)) {
					return false;
				}
				focals_union = this->fod.set_union(setj, focals_union);
			}
			return focals_union.count() == this->fod.size();
		}

		/// Mass function without internal conflict is a one where all pairs of
		/// focal elements have a non-empty intersection.
		/// \f[
		/// \forall A, B \subseteq \Omega, m(A) > 0, m(B) > 0 : A \cap B \neq
		/// \emptyset
		/// \f]
		/// Therefore, there is internal conflict when there is at least one empty intersection between pairs,
		/// i.e. when the intersection of every focal elements is empty
		bool has_internal_conflict() const {
			return is_arbitrary();
		}
	};

} // namespace ow_bft

#endif // OW_BFT_MASS_HPP