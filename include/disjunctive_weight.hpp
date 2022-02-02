#ifndef EFFICIENT_DST_DISJUNCTIVE_WEIGHT_HPP
#define EFFICIENT_DST_DISJUNCTIVE_WEIGHT_HPP

#include <decomposition_weight.hpp>
#include <powerset_function.hpp>
#include <zeta_transform.hpp>

namespace efficient_DST{

	template <typename T, size_t N>
	class disjunctive_weight : public decomposition_weight<T, N> {
	public:

		disjunctive_weight(const disjunctive_weight<T, N>& v) : decomposition_weight<T, N>(v.get_definition())
		{}

		disjunctive_weight(const powerset_btree<T, N>& support) : decomposition_weight<T, N>(support)
		{}

		disjunctive_weight(FOD<N>& fod) : decomposition_weight<T, N>(fod)
		{}

		disjunctive_weight(const zeta_transform<T, N, down_inclusion<T, N> >& b) : decomposition_weight<T, N>(b.inversion(operation_type_t::multiplication))
		{
			invert_values(this->definition);
		}


		/*
		 * The disjunctive weight function is the inverse of the multiplicative Möbius transform of the implicability function.
		 * So, invert its values before computing it.
		 */
		powerset_btree<T, N> inverted_definition() const {
			powerset_btree<T, N> inverted_definition(this->definition);
			invert_values(inverted_definition);
			return inverted_definition;
		}

		void set_value(const std::vector<std::string>& labels, const T& value) {
			set_value_directly(this->definition.get_FOD()->to_set(labels), value);
		}

		void set_value_directly(const std::bitset<N>& set, const T& value) {
			this->definition.insert(set, value);
			// The following part ensures that v(emptyset) is defined (as it is not when directly building v).
			// Indeed, v(emptyset) may be required for the computation of the implicability function.
			// v(emptyset)^{-1} is equal to the product of all other weights.
			set_N_value<T, N>* emptyset = this->definition.sub_fod_of_size(0);
			if(emptyset){
				emptyset->value /= value;
			}else{
				this->definition.insert(std::bitset<N>(0), 1/value);
			}
		}

		void set_fod_value(const T& value) {
			set_value_directly(~std::bitset<N>(0), value);
		}

		void nullify(const std::vector<std::string>& labels) {
			//if(this->definition.get_FOD()->to_set(labels).count() == 0)
			//	return;
			set_N_value<T, N>* A = this->definition[labels];
			if(A){
				set_N_value<T, N>* emptyset = this->definition.sub_fod_of_size(0);
				if(A != emptyset){
					if(emptyset){
						emptyset->value *= A->value;
					}else{
						// if there is no emptyset in definition, this means that its associated value is 1.
						// Thus, as above, nullifying A implies multiplying 1 with A->value,
						// which translates here to inserting A->value at emptyset.
						this->definition.insert(std::bitset<N>(0), A->value);
					}
				}
				this->definition.nullify(A);
			}
		}

		void compute_emptyset_value_from_definition(){
			compute_emptyset_value_from_definition(this->definition);
		}

		static void compute_emptyset_value_from_definition(powerset_btree<T, N>& definition){
			std::bitset<N> emptyset = 0;
			const std::vector<set_N_value<T, N>* >& focal_log_elements = definition.elements();
			T val = 1;
			for (size_t i = 0; i < focal_log_elements.size(); ++i){
				if (focal_log_elements[i]->set != emptyset)
					val /= focal_log_elements[i]->value;
			}
			definition.update_or_insert(emptyset, val);
		}

		template <class fusion_rule>
		disjunctive_weight<T, N> fuse_with(const disjunctive_weight<T, N>& v2) const {
			const fusion_rule fusion;
			return fusion(*this, v2);
		}

	protected:

		static void invert_values(powerset_btree<T, N>& values) {
			std::vector<set_N_value<T, N>* > elements = values.elements();
			for (size_t i = 0; i < elements.size(); ++i){
				elements[i]->value = 1 / elements[i]->value;
			}
		}
	};

} // namespace efficient_DST

#endif // EFFICIENT_DST_DISJUNCTIVE_WEIGHT_HPP
