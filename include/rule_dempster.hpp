#ifndef OW_BFT_RULE_DEMPSTER_HPP
#define OW_BFT_RULE_DEMPSTER_HPP

#include <commonality.hpp>
#include <conjunctive_decomposition.hpp>
#include <mass.hpp>
#include <converter_to_mass_focal_elements.hpp>
#include <rule_conjunctive.hpp>

namespace ow_bft{

	struct rule_Dempster {
		std::string to_string() const {
			return "conjunctive rule";
		}

	public:

		template <typename T>
		powerset_btree<T> fusion_of_commonalities_to_mass_focal_elements(const commonality<T>& q1, const commonality<T>& q2, T precision) const {

			// initialization
			powerset_btree<T> m12_focal_elements(q1.get_FOD(), q1.block_size);

			// store the conjunction of q1 and q2 (i.e. q12(A) = q1(A).q2(A)) in m12

			const std::vector<set_N_value<T>* >& elements1 = q1.get_special_elements().elements();
			powerset_btree<T> q2_special_elements_copy(q2.get_special_elements());

			for (size_t i = 0; i < elements1.size(); ++i) {
				const std::vector<std::string>& labels = q1.get_FOD().to_labels(elements1[i]->fod_elements);
				set_N_value<T>* q2_matching_element = q2_special_elements_copy[labels];
				T q2_element_value;
				if(q2_matching_element){
					q2_element_value = q2_matching_element->value;
					q2_special_elements_copy.nullify(q2_matching_element);
				}else{
					q2_element_value = q2[labels];
				}

				m12_focal_elements.insert(elements1[i]->fod_elements, elements1[i]->value * q2_element_value);
			}

			const std::vector<set_N_value<T>* >& elements2 = q2_special_elements_copy.elements();

			for (size_t i = 0; i < elements2.size(); ++i) {
				const std::vector<std::string>& labels = q2.get_FOD().to_labels(elements2[i]->fod_elements);

				m12_focal_elements.insert(labels, q1[labels] * elements2[i]->value);
			}

			// transform m12_focal_elements (which is q12 at this point) into focal elements of a mass function, i.e. into the actual m12 focal elements
			commonality<T>::to_mass_focal_elements_without_initialization(m12_focal_elements);

			rule_conjunctive::fusion_of_commonalities_to_mass_focal_elements(q1, q2, precision)

			return m12_focal_elements;
		}
/*
		template <typename T>
		commonality<T> operator()(conjunctive_weights<T> const& w1, conjunctive_weights<T> const& w2) const {
			commonality<T> q12(q1.fod);

			// initialization
			mass<T> m12(q1.fod);
			// store the conjunction of q1 and q2 (i.e. q12(A) = q1(A).q2(A)) in m12 and retrieve these added nodes
			// for transformation into masses
			std::vector<std::vector<set_N_value* > > m12_sets = m12.focal_elements
					.fill_with_union_of_elements_and_return_by_set_cardinality(q1.before_changes, q2.before_changes, multiply<T>);

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
		}*/

		template <typename T>
		mass<T> operator()(const mass<T>& m1, const mass<T>& m2) const {
			commonality<T> q1(m1);
			commonality<T> q2(m2);

			const powerset_btree<T>& m12_focal_elements = fusion_of_commonalities_to_mass_focal_elements(q1, q2, m1.precision);
			return mass<T>(m12_focal_elements);
		}

		template <typename T>
		commonality<T> operator()(const commonality<T>& q1, const commonality<T>& q2) const {
			const powerset_btree<T>& m12_focal_elements = fusion_of_commonalities_to_mass_focal_elements(q1, q2, q1.precision);
			return commonality<T>(m12_focal_elements);
		}

		template <typename T>
		conjunctive_decomposition<T> operator()(const conjunctive_decomposition<T>& w1, const conjunctive_decomposition<T>& w2) const {
			const powerset_btree<T>& m12_focal_elements = fusion_of_commonalities_to_mass_focal_elements(
					w1.get_commonality_equivalent(),
					w2.get_commonality_equivalent(),
					w1.precision);
			return conjunctive_decomposition<T>(m12_focal_elements);
		}
	};

} // namespace ow_bft

#endif // OW_BFT_RULE_DEMPSTER_HPP
