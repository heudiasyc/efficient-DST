#ifndef EFFICIENT_DST_RULE_CONJUNCTIVE_CAUTIOUS_HPP
#define EFFICIENT_DST_RULE_CONJUNCTIVE_CAUTIOUS_HPP

#include <conjunctive_weight.hpp>
#include <mass.hpp>
#include <rule_conjunctive.hpp>
#include <detail/is_small.hpp>

namespace efficient_DST{

	template <typename T>
	static T multiply(T val1, T val2){
		if(val1 > 0 && val2 > 0)
			return val1 * val2;
		else
			return 0;
	}

	struct rule_conjunctive_cautious
	{
		std::string to_string() const
		{
			return "cautious conjunctive rule";
		}

		template <typename T>
		static conjunctive_weight<T> operator()(const conjunctive_weight<T>& w1, const conjunctive_weight<T>& w2) {


			// initialization
			mass<T> m12(q1.fod);
			// store the conjunction of q1 and q2 (i.e. q12(A) = q1(A).q2(A)) in m12 and retrieve these added nodes
			// for transformation into masses
			std::vector<std::vector<set_N_value* > > m12_sets = m12.focal_elements
					.fill_with_union_of_elements_and_return_by_set_cardinality(
							w1.commonality_equivalent->focal_elements,
							w2.commonality_equivalent->focal_elements,
							multiply<T>
			);

			// copy of the conjunction of q1 and q2 in q12
			q12.set_focal_elements(powerset_btree<T>(m12));
			// computation and storage of masses in m12
			for (int i = m12_sets.size()-1; i >= 0; --i) {
				for (int j = 0; j < m12_sets[i].size(); ++j) {
					std::vector<set_N_value* > subsets = m12.focal_elements.strict_subsets_of(m12_sets[i][j]->fod_elements);
					for (int k = 0; k < subsets.size(); ++k) {
						subsets[k]->value -= m12_sets[i][j]->value;
					}
				}
			}
			q12.mass_equivalent = m12.focal_elements;
			return q12;
		}
	};

} // namespace efficient_DST

#endif // EFFICIENT_DST_RULE_CONJUNCTIVE_CAUTIOUS_HPP
