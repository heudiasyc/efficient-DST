#ifndef OW_BFT_RULE_CONJUNCTIVE_HPP
#define OW_BFT_RULE_CONJUNCTIVE_HPP

#include <commonality.hpp>
#include <conjunctive_decomposition.hpp>
#include <mass.hpp>
#include <converter_to_mass_focal_elements.hpp>

namespace ow_bft{

	template <typename T>
	static T multiply(const powerset_btree<T>& powerset, T val1, T val2){
		if(val1 == powerset.null || val2 == powerset.null)
			return powerset.null;
		else
			return val1 * val2;
	}

	struct rule_conjunctive {
		std::string to_string() const {
			return "conjunctive rule";
		}

	public:

		template <typename T>
		powerset_btree<T> fusion_of_commonalities_to_mass_focal_elements(const commonality<T>& q1, const commonality<T>& q2, T precision) const {

			const FOD& fod = q1.get_FOD();

			// initialization
			powerset_btree<T> m12_focal_elements(fod, q1.block_size);

			// store the conjunction of q1 and q2 (i.e. q12(A) = q1(A).q2(A)) in m12

			const std::vector<set_N_value<T>* >& elements1 = q1.get_special_elements().elements();
			const std::vector<set_N_value<T>* >& elements2 = q2.get_special_elements().elements();

			for (size_t i = 0; i < elements1.size(); ++i) {
				for (size_t ii = 0; ii < elements2.size(); ++ii) {
					const boost::dynamic_bitset<>& m12_focal_element = fod.set_intersection(
							fod.to_set(elements1[i]->fod_elements),
							fod.to_set(q2.get_FOD().to_labels(elements2[i]->fod_elements))
						);
					const set_N_value<T>* m12_focal_element_node = m12_focal_elements[m12_focal_element];
					if(!m12_focal_element_node)
						m12_focal_elements.insert(m12_focal_element, q1.find(m12_focal_element) * q2[fod.to_labels(fod.to_elements(m12_focal_element))]);
				}
			}

			// transform m12_focal_elements (which is q12 at this point) into focal elements of a mass function, i.e. into the actual m12 focal elements
			commonality<T>::to_mass_focal_elements_without_initialization(m12_focal_elements);

			/*
			 * NE DEVRAIT PAS ÊTRE NÉCESSAIRE CAR S'IL Y A CONFLIT, L'INTERSECTION EMPTYSET A DÛ ÊTRE GÉNÉRÉE ET DONC AFFECTÉE DÉJÀ À EMPTYSET
			const std::vector<set_N_value<T>* >& m12_elements = m12_focal_elements.elements();
			T sum = 0;
			for (size_t i = 0; i < m12_elements.size(); ++i) {
				sum += m12_elements[i]->value;
			}
			set_N_value<T>* emptyset = m12_focal_elements.sub_fod_of_size(0);
			if(emptyset)
				sum -= emptyset->value;

			m12_focal_elements.set_value_of_sub_fod_of_size(0, 1 - sum);
			*/
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

#endif // OW_BFT_RULE_CONJUNCTIVE_HPP
