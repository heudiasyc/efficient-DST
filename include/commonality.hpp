#ifndef OW_BFT_COMMONALITY_HPP
#define OW_BFT_COMMONALITY_HPP

#include <mass_aggregate.hpp>
#include <mass.hpp>

namespace ow_bft{

	template <typename T = double>
	class commonality : public mass_aggregate<T>{
	protected:

		static T compute_aggregation_at_emptyset(const powerset_btree<T>& m_focal_elements) {
			return 1;
		}

		static T compute_aggregation_at_fod(const powerset_btree<T>& m_focal_elements) {
			return m_focal_elements.sub_fod_of_size(m_focal_elements.fod->size())->value;
		}

		static T compute_aggregation(const powerset_btree<T>& m_focal_elements, const boost::dynamic_bitset<>& set) {
			T sum = 0;
			const std::vector<set_N_value<T>* >& supersets = m_focal_elements.supersets_of(set);

			if(supersets.size() == m_focal_elements.size())
				return 1;

			for (size_t i = 0; i < supersets.size(); ++i) {
				sum += supersets[i]->value;
			}
			return sum;
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

		commonality(const mass<T>& m) : mass_aggregate<T>(m)
		{
			compute_values_for_mass_focal_elements(this->mass_equivalent.get_focal_elements(), this->special_elements);
		}

		commonality(const powerset_btree<T>& m_focal_elements) : mass_aggregate<T>(m_focal_elements)
		{
			compute_values_for_mass_focal_elements(this->mass_equivalent.get_focal_elements(), this->special_elements);
		}

		commonality(const commonality<T>& q) : commonality(q.get_mass_equivalent().get_focal_elements(), q.get_special_elements())
		{}

		commonality(const powerset_btree<T>& m_focal_elements, const powerset_btree<T>& _special_elements) :
			mass_aggregate<T>(m_focal_elements, _special_elements)
		{}

		commonality(const mass_aggregate<T>& ma) : commonality(ma.get_mass_equivalent())
		{}

		commonality(const FOD& fod) : mass_aggregate<T>(fod)
		{
			compute_values_for_mass_focal_elements(this->mass_equivalent.get_focal_elements(), this->special_elements);
		}

		commonality(const FOD& fod, const Special_case s_case) : mass_aggregate<T>(fod, s_case)
		{
			compute_values_for_mass_focal_elements(this->mass_equivalent.get_focal_elements(), this->special_elements);
		}

		template <class fusion_rule>
		commonality<T> apply(const fusion_rule fusion, const commonality<T>& q2) const {
			return fusion(*this, q2);
		}

		static void compute_values_for_mass_focal_elements(const powerset_btree<T>& m_focal_elements, powerset_btree<T>& special_elements) {
			const std::vector<set_N_value<T>* >& elements = m_focal_elements.elements();
			// pre-calculation for all focal elements
			for (size_t i = 0; i < elements.size(); ++i) {
				special_elements.insert(elements[i]->set, compute_aggregation(m_focal_elements, elements[i]->set));
			}
		}

		static void to_mass_focal_elements(const powerset_btree<T>& q_tree, powerset_btree<T>& m_tree) {

			m_tree.copy(q_tree);
			to_mass_focal_elements_without_initialization(m_tree);
		}

		static void to_mass_focal_elements_without_initialization(powerset_btree<T>& m_tree) {

			std::unordered_map<size_t, std::vector<set_N_value<T>* > > elements_by_cardinality = m_tree.elements_by_set_cardinality();
			const std::vector<std::vector<set_N_value<T>* >* >& ordered_vector = m_tree.get_vector_of_vectors_ordered_by_cardinality(elements_by_cardinality);

			// computation based on f_elements
			for (size_t c = ordered_vector.size()-1; c > 0; --c) {
				for (size_t i = 0; i < ordered_vector[c]->size(); ++i) {
					const std::vector<set_N_value<T>* >& subsets = m_tree.strict_subsets_of((*ordered_vector[c])[i]->set);
					for (size_t k = 0; k < subsets.size(); ++k) {
						subsets[k]->value -= (*ordered_vector[c])[i]->value;
					}
				}
			}
		}
	};

} // namespace ow_bft

#endif // OW_BFT_COMMONALITY_HPP
