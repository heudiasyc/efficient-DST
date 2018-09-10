#ifndef OW_BFT_BELIEF_HPP
#define OW_BFT_BELIEF_HPP

#include <mass_aggregate.hpp>

namespace ow_bft{

	template <typename T = double>
	class belief : public mass_aggregate<T> {
	protected:

		static T compute_aggregation_at_emptyset(const powerset_btree<T>& m_focal_elements) {
			return 0;
		}

		static T compute_aggregation_at_fod(const powerset_btree<T>& m_focal_elements) {
			return 1;
		}

		static T compute_aggregation(const powerset_btree<T>& m_focal_elements, const boost::dynamic_bitset<>& key) {
			T sum = 0;
			const std::vector<set_N_value<T>* >& subsets = m_focal_elements.subsets_of(key);

			if(subsets.size() == 1 && subsets[0] == m_focal_elements.sub_fod_of_size(0))
				return compute_aggregation_at_emptyset(m_focal_elements);

			for (size_t i = 0; i < subsets.size(); ++i) {
				sum += subsets[i]->value;
			}
			return sum;
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

		T compute_aggregation(const boost::dynamic_bitset<>& key) const {
			return compute_aggregation(this->mass_equivalent.get_focal_elements(), key);
		}

		T compute_aggregation(const std::vector<fod_element*>& fod_elements) const {
			return compute_aggregation(this->mass_equivalent.get_focal_elements(), fod_elements);
		}

	public:

		belief(const mass<T>& m) : mass_aggregate<T>(m)
		{
			compute_values_for_mass_focal_elements(this->mass_equivalent.get_focal_elements(), this->special_elements);
		}

		belief(const powerset_btree<T>& m_focal_elements) : mass_aggregate<T>(m_focal_elements)
		{
			compute_values_for_mass_focal_elements(this->mass_equivalent.get_focal_elements(), this->special_elements);
		}

		belief(const belief<T>& bel) : belief(bel.get_mass_equivalent().get_focal_elements(), bel.get_special_elements())
		{}

		belief(const powerset_btree<T>& m_focal_elements, const powerset_btree<T>& _special_elements) :
			mass_aggregate<T>(m_focal_elements, _special_elements)
		{}

		belief(const mass_aggregate<T>& ma) : belief(ma.get_mass_equivalent())
		{}

		belief(const FOD& fod) : mass_aggregate<T>(fod)
		{
			compute_values_for_mass_focal_elements(this->mass_equivalent.get_focal_elements(), this->special_elements);
		}

		belief(const FOD& fod, const Special_case s_case) : mass_aggregate<T>(fod, s_case)
		{
			compute_values_for_mass_focal_elements(this->mass_equivalent.get_focal_elements(), this->special_elements);
		}

		virtual ~belief()
		{}

		template <class fusion_rule>
		const belief<T> apply(const fusion_rule fusion, const belief<T>& b2) const {
			return fusion(*this, b2);
		}

		static void compute_values_for_mass_focal_elements(const powerset_btree<T>& m_focal_elements, powerset_btree<T>& special_elements) {
			const std::vector<set_N_value<T>* >& elements = m_focal_elements.elements();
			// pre-calculation for all focal elements
			for (size_t i = 0; i < elements.size(); ++i) {
				special_elements.insert(elements[i]->fod_elements, compute_aggregation(m_focal_elements, elements[i]->fod_elements));
			}
		}

		static void to_mass_focal_elements(const powerset_btree<T>& bel_tree, powerset_btree<T>& m_tree) {

			m_tree.copy(bel_tree);

			const std::vector<std::vector<set_N_value<T>* > >& elements_by_cardinality = m_tree().elements_by_set_cardinality();

			// computation based on f_elements
			for (size_t c = 0; c < elements_by_cardinality.size()-1; ++c) {
				for (size_t i = 0; i < elements_by_cardinality[c].size(); ++i) {
					const std::vector<set_N_value<T>* >& supersets = m_tree.strict_supersets_of(elements_by_cardinality[c][i]->fod_elements);
					for (size_t k = 0; k < supersets.size(); ++k) {
						supersets[k]->value -= elements_by_cardinality[c][i]->value;
					}
				}
			}
		}
	};

} // namespace ow_bft

#endif // OW_BFT_BELIEF_HPP
