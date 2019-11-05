#ifndef EFFICIENT_DST_RULE_DUBOIS_PRADE_HPP
#define EFFICIENT_DST_RULE_DUBOIS_PRADE_HPP

#include <mass.hpp>

namespace efficient_DST{

	template <typename T, size_t N>
	class rule_Dubois_Prade {
	public:

		std::string to_string() const {
			return "Dubois-Prade adaptative rule";
		}

		mass<T, N> operator()(const mass<T, N>& m1, const mass<T, N>& m2) const {
			const std::vector<set_N_value<T, N>* >& focal_sets_1 = m1.get_definition().elements();
			const std::vector<set_N_value<T, N>* >& focal_sets_2 = m1.get_definition().elements();
			powerset_btree<T, N> focal_sets_12(m1.get_definition().get_FOD(), m1.get_definition().get_block_size());

			for (size_t i1 = 0; i1 < focal_sets_1.size(); ++i1){
				for (size_t i2 = 0; i2 < focal_sets_2.size(); ++i2){
					std::bitset<N> set = focal_sets_1[i1]->set & focal_sets_2[i2]->set;

					if(FOD<N>::is_emptyset(set)){
						set = focal_sets_1[i1]->set | focal_sets_2[i2]->set;
					}

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
	};

} // namespace efficient_DST

#endif // EFFICIENT_DST_RULE_DUBOIS_PRADE_HPP
