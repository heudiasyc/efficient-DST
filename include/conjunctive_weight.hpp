#ifndef EFFICIENT_DST_CONJUNCTIVE_WEIGHT_HPP
#define EFFICIENT_DST_CONJUNCTIVE_WEIGHT_HPP

#include <decomposition_weight.hpp>
#include <powerset_function.hpp>
#include <zeta_transform.hpp>

namespace efficient_DST{

	template <typename T = double>
	class conjunctive_weight : public decomposition_weight<T> {
	public:

		conjunctive_weight(const conjunctive_weight<T>& w) : decomposition_weight<T>(w.get_definition())
		{}

		conjunctive_weight(const powerset_btree<T>& focal_log_sets_values) : decomposition_weight<T>(focal_log_sets_values)
		{}

		conjunctive_weight(FOD& fod) : decomposition_weight<T>(fod)
		{}

		conjunctive_weight(const zeta_transform<T>& q) : decomposition_weight<T>(q)
		{
			if (q.order_relation != order_relation_t::superset) {
				std::cerr << "The given Möbius aggregate is not the commonality function and thus can only be the implicability one. "
						<< "\nDoing so, you got the disjunctive weight function instead." << std::endl;
			}
			invert_values(this->definition);
		}


		/*
		 * The conjunctive weight function is the inverse of the multiplicative Möbius transform of the commonality function.
		 * So, invert its values before computing it.
		 */
		powerset_btree<T> inverted_definition() const {
			powerset_btree<T> inverted_definition(this->definition);
			invert_values(inverted_definition);
			return inverted_definition;
		}

		void set_value(const std::vector<std::string>& labels, const T& value) {
			set_value_directly(this->definition.get_FOD()->to_set(labels), value);
		}

		void set_value_directly(const boost::dynamic_bitset<>& set, const T& value) {
			this->definition.insert(set, value);
			// The following part ensures that w(fod) is defined (as it is not when directly building w).
			// Indeed, w(fod) may be required for the computation of the commonality function.
			// w(fod)^{-1} is equal to the product of all other weights.
			set_N_value<T>* fod = this->definition.sub_fod_of_size(this->definition.get_FOD_size());
			if(fod){
				fod->value /= value;
			}else{
				this->definition.set_value_of_sub_fod_of_size(this->definition.get_FOD_size(), 1/value);
			}
		}

		void set_emptyset_value(const T& value) {
			set_value_directly(boost::dynamic_bitset<>(this->definition.get_FOD_size()), value);
		}

		void nullify(const std::vector<std::string>& labels) {
			if(this->definition.get_FOD()->to_set(labels).count() == 1)
				return;
			set_N_value<T>* A = this->definition[labels];
			if(A){
				set_N_value<T>* fod = this->definition.sub_fod_of_size(this->definition.get_FOD_size());
				if(A != fod){
					if(fod){
						fod->value *= A->value;
					}else{
						this->definition.set_value_of_sub_fod_of_size(this->definition.get_FOD_size(), A->value);
					}
				}
				this->definition.nullify(A);
			}
		}

		void compute_fod_value_from_definition(){
			compute_fod_value_from_definition(this->definition);
		}

		static void compute_fod_value_from_definition(powerset_btree<T>& definition){
			boost::dynamic_bitset<> fod(definition.get_FOD_size());
			fod.set();
			const std::vector<set_N_value<T>* >& focal_log_elements = definition.strict_subsets_of(fod);
			T val = 1;
			for (size_t i = 0; i < focal_log_elements.size(); ++i){
				val /= focal_log_elements[i]->value;
			}
			set_N_value<T>* fod_set_N_value = definition.sub_fod_of_size(definition.get_FOD_size());
			if(fod_set_N_value){
				fod_set_N_value->value = val;
			}else{
				definition.set_value_of_sub_fod_of_size(definition.get_FOD_size(), val);
			}
		}

		template <class fusion_rule>
		conjunctive_weight<T> apply(const conjunctive_weight<T>& w2) const {
			const fusion_rule fusion;
			return fusion(*this, w2);
		}

	protected:

		static void invert_values(powerset_btree<T>& values) {
			std::vector<set_N_value<T>* > elements = values.elements();
			for (size_t i = 0; i < elements.size(); ++i){
				elements[i]->value = 1 / elements[i]->value;
			}
		}
	};

} // namespace efficient_DST

#endif // EFFICIENT_DST_CONJUNCTIVE_WEIGHT_HPP
