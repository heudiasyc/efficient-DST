#ifndef EFFICIENT_DST_PLAUSIBILITY_HPP
#define EFFICIENT_DST_PLAUSIBILITY_HPP

#include <implicability.hpp>
#include <mobius_aggregate.hpp>

namespace efficient_DST{

	template <typename T = double>
	class plausibility : public mobius_aggregate<T> {
	protected:

		static T compute_aggregation_at_emptyset(const powerset_btree<T>& m_focal_elements) {
			return 0;
		}

		static T compute_aggregation_at_fod(const powerset_btree<T>& m_focal_elements) {
			return implicability<T>::compute_aggregation_at_fod(m_focal_elements);
		}

		static T compute_aggregation(const powerset_btree<T>& m_focal_elements, const boost::dynamic_bitset<>& set) {
			return 1 - implicability<T>::compute_aggregation(m_focal_elements, m_focal_elements.fod->set_negate(set));
		}

		static T compute_aggregation(const powerset_btree<T>& m_focal_elements, const std::vector<fod_element*>& fod_elements) {
			return compute_aggregation(m_focal_elements, m_focal_elements.fod->to_set(fod_elements));
		}

		T compute_aggregation_at_emptyset() const {
			return compute_aggregation_at_emptyset(this->mass_equivalent.get_focal_elements());
		}

		T compute_aggregation_at_fod() const {
			return compute_aggregation_at_fod(this->mass_equivalent.get_focal_elements());
		}

		T compute_aggregation(const boost::dynamic_bitset<>& set) const {
			return compute_aggregation(this->mass_equivalent.get_focal_elements(), set);
		}

	public:

		plausibility(const mass<T>& m) : mass_aggregate<T>(m)
		{
			compute_values_for_negation_of_mass_focal_elements(this->mass_equivalent.get_focal_elements(), this->special_elements);
		}

		plausibility(const mobius_aggregate<T>& ma) : plausibility(ma.get_mass_equivalent())
		{}

		plausibility(const plausibility<T>& p) : mass_aggregate<T>(p.get_mass_equivalent())
		{
			this->special_elements.copy(p.get_special_elements());
		}

		plausibility(const FOD& fod) : mass_aggregate<T>(fod)
		{
			compute_values_for_negation_of_mass_focal_elements(this->mass_equivalent.get_focal_elements(), this->special_elements);
		}

		plausibility(const FOD& fod, const Special_case s_case) : mass_aggregate<T>(fod, s_case)
		{
			compute_values_for_negation_of_mass_focal_elements(this->mass_equivalent.get_focal_elements(), this->special_elements);
		}

		template <class fusion_rule>
		plausibility<T> apply(const fusion_rule fusion, const plausibility<T>& p2) const {
			return fusion(*this, p2);
		}

		static void compute_values_for_negation_of_mass_focal_elements(const powerset_btree<T>& m_focal_elements, powerset_btree<T>& special_elements) {
			const std::vector<set_N_value<T>* >& elements = m_focal_elements.elements();
			// pre-calculation for all focal elements
			for (size_t i = 0; i < elements.size(); ++i) {
				const boost::dynamic_bitset<>& focal_element_negation = m_focal_elements.fod->set_negate(elements[i]->set);
				special_elements.insert(focal_element_negation, 1 - implicability<T>::compute_aggregation(m_focal_elements, elements[i]->set));
			}
		}

		/*
		 * pl_tree must feature the value associated to the NEGATION of every focal element
		 */
		static void to_mass_focal_elements(const powerset_btree<T>& pl_tree, powerset_btree<T>& m_tree) {
			const std::vector<set_N_value<T>* >& elements = pl_tree.elements();
			powerset_btree<T> bel_tree(pl_tree.fod, pl_tree.block_size);

			// pre-calculation for all focal elements
			for (size_t i = 0; i < elements.size(); ++i) {
				const boost::dynamic_bitset<>& focal_element = pl_tree.fod->set_negate(elements[i]->set);
				bel_tree.insert(focal_element, 1 - elements[i]->value);
			}
			implicability<T>::to_mass_focal_elements(bel_tree, m_tree);
		}
	};

} // namespace efficient_DST

#endif // EFFICIENT_DST_PLAUSIBILITY_HPP

