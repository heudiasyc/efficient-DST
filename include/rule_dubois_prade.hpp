#ifndef EFFICIENT_DST_RULE_DUBOIS_PRADE_HPP
#define EFFICIENT_DST_RULE_DUBOIS_PRADE_HPP

#include <mass.hpp>

namespace efficient_DST{

	template <typename T = double>
	class rule_Dubois_Prade {
		std::string to_string() const {
			return "Dubois-Prade adaptative rule";
		}

	public:

		mass<T> operator()(const mass<T>& m1, const mass<T>& m2) const {
			const std::vector<set_N_value<T>* >& focal_sets_1 = m1.get_definition().elements();
			const std::vector<set_N_value<T>* >& focal_sets_2 = m1.get_definition().elements();
			powerset_btree<T> focal_sets_12(m1.get_definition().get_FOD(), m1.get_definition().get_block_size());

			for (size_t i1 = 0; i1 < focal_sets_1.size(); ++i1){
				for (size_t i2 = 0; i2 < focal_sets_2.size(); ++i2){
					boost::dynamic_bitset<> set = FOD::set_intersection(focal_sets_1[i1]->set, focal_sets_2[i2]->set);

					if(FOD::is_emptyset(set)){
						set = FOD::set_union(focal_sets_1[i1]->set, focal_sets_2[i2]->set);
					}

					set_N_value<T>* node = focal_sets_12[set];
					if (node){
						node->value += focal_sets_1[i1]->value * focal_sets_2[i2]->value;
					}else{
						focal_sets_12.insert(set, focal_sets_1[i1]->value * focal_sets_2[i2]->value);
					}
				}
			}
			mass<T> m12(focal_sets_12);
			m12.remove_negligible_values();
			m12.normalize();
			return m12;
		}
	};

} // namespace efficient_DST

#endif // EFFICIENT_DST_RULE_DUBOIS_PRADE_HPP
