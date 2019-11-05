#ifndef EFFICIENT_DST_CONJUNCTIVE_WEIGHT_HPP
#define EFFICIENT_DST_CONJUNCTIVE_WEIGHT_HPP

#include <decomposition_weight.hpp>
#include <powerset_function.hpp>
#include <zeta_transform.hpp>

namespace efficient_DST{

	template <typename T, size_t N>
	class conjunctive_weight : public decomposition_weight<T, N> {
	public:

		conjunctive_weight(const conjunctive_weight<T, N>& w) : decomposition_weight<T, N>(w.get_definition())
		{}

		conjunctive_weight(const powerset_btree<T, N>& focal_log_sets_values) : decomposition_weight<T, N>(focal_log_sets_values)
		{}

		conjunctive_weight(FOD<N>& fod) : decomposition_weight<T, N>(fod)
		{}

		conjunctive_weight(const zeta_transform<T, N>& q) : decomposition_weight<T, N>(q)
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
			// The following part ensures that w(fod) is defined (as it is not when directly building w).
			// Indeed, w(fod) may be required for the computation of the commonality function.
			// w(fod)^{-1} is equal to the product of all other weights.
			set_N_value<T, N>* fod = this->definition.sub_fod_of_size(N);
			if(fod){
				fod->value /= value;
			}else{
				this->definition.insert(~std::bitset<N>(0), 1/value);
			}
		}

		void set_emptyset_value(const T& value) {
			set_value_directly(std::bitset<N>(N), value);
		}

		void nullify(const std::vector<std::string>& labels) {
			//if(this->definition.get_FOD()->to_set(labels).count() == 1)
			//	return;
			set_N_value<T, N>* A = this->definition[labels];
			if(A){
				set_N_value<T, N>* fod = this->definition.sub_fod_of_size(N);
				if(A != fod){
					if(fod){
						fod->value *= A->value;
					}else{
						this->definition.insert(~std::bitset<N>(0), A->value);
					}
				}
				this->definition.nullify(A);
			}
		}

		void compute_fod_value_from_definition(){
			compute_fod_value_from_definition(this->definition);
		}

		static void compute_fod_value_from_definition(powerset_btree<T, N>& definition){
			std::bitset<N> fod_set = 0;
			fod_set.set();
			const std::vector<set_N_value<T, N>* >& focal_log_elements = definition.strict_subsets_of(fod_set);
			T val = 1;
			for (size_t i = 0; i < focal_log_elements.size(); ++i){
				val /= focal_log_elements[i]->value;
			}
			set_N_value<T, N>* fod_set_N_value = definition.sub_fod_of_size(N);
			if(fod_set_N_value){
				fod_set_N_value->value = val;
			}else{
				definition.insert(fod_set, val);
			}
		}

		template <class fusion_rule>
		conjunctive_weight<T, N> apply(const conjunctive_weight<T, N>& w2) const {
			const fusion_rule fusion;
			return fusion(*this, w2);
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

#endif // EFFICIENT_DST_CONJUNCTIVE_WEIGHT_HPP
