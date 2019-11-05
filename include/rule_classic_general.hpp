#ifndef EFFICIENT_DST_RULE_CLASSIC_GENERAL_HPP
#define EFFICIENT_DST_RULE_CLASSIC_GENERAL_HPP

#include <commonality.hpp>
#include <conjunctive_weight.hpp>
#include <mass.hpp>

namespace efficient_DST{

	/*
	 * Virtual class.
	 * Used to define the conjunctive and disjunctive rules.
	 */
	template <typename T, size_t N>
	class rule_classic_general {
	public:

		std::string to_string() const {
			return "General classic rule";
		}

		virtual ~rule_classic_general(){}

		mass<T, N> operator()(
				const mass<T, N>& m1,
				const mass<T, N>& m2,
				std::function<const std::bitset<N>(
						const std::bitset<N>&,
						const std::bitset<N>&
						)> set_operator
				) const {
			const std::vector<set_N_value<T, N>* >& focal_sets_1 = m1.get_definition().elements();
			const std::vector<set_N_value<T, N>* >& focal_sets_2 = m2.get_definition().elements();
			powerset_btree<T, N> focal_sets_12(m1.get_definition().get_FOD(), m1.get_definition().get_block_size());

			for (size_t i1 = 0; i1 < focal_sets_1.size(); ++i1){
				for (size_t i2 = 0; i2 < focal_sets_2.size(); ++i2){
					const std::bitset<N>& set = set_operator(focal_sets_1[i1]->set, focal_sets_2[i2]->set);
					set_N_value<T, N>* node = focal_sets_12[set];
					if (node){
						node->value += focal_sets_1[i1]->value * focal_sets_2[i2]->value;
					}else{
						focal_sets_12.insert(set, focal_sets_1[i1]->value * focal_sets_2[i2]->value);
					}
				}
			}
			mass<T, N> m12(focal_sets_12);
			m12.remove_negligible_values();
			m12.normalize();
			return m12;
		}

	protected:

		static T multiply(const T& val1, const T& val2){
			return val1 * val2;
		}

		powerset_btree<T, N> weight_fusion(const powerset_btree<T, N>& w1_definition, const powerset_btree<T, N>& w2_definition) const {
			powerset_btree<T, N> w12_definition(w1_definition.get_FOD(), w1_definition.get_block_size());
			w12_definition.fill_with_union_of_powersets(w1_definition, w2_definition, multiply, 1);
			decomposition_weight<T, N>::remove_negligible_values(w12_definition);
			decomposition_weight<T, N>::normalize(w12_definition);
			return w12_definition;
		}
	};

} // namespace efficient_DST

#endif // EFFICIENT_DST_RULE_CLASSIC_GENERAL_HPP
