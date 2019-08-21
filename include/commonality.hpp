#ifndef EFFICIENT_DST_COMMONALITY_HPP
#define EFFICIENT_DST_COMMONALITY_HPP

#include <EMT_aggregation.hpp>
#include <mass.hpp>
#include <mobius_aggregate.hpp>

namespace efficient_DST{

	template <typename T = double>
	class commonality : public mobius_aggregate<T>{
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
			return compute_aggregation_at_emptyset(this->mass_equivalent.get_focal_sets_values());
		}

		T compute_aggregation_at_fod() const {
			return compute_aggregation_at_fod(this->mass_equivalent.get_focal_sets_values());
		}

		T compute_aggregation(const boost::dynamic_bitset<>& set) const {
			return compute_aggregation(this->mass_equivalent.get_focal_sets_values(), set);
		}

	public:

		commonality(const mass<T>& m) : mobius_aggregate<T>(m)
		{
			compute_values_for_mass_focal_elements(this->mass_equivalent.get_focal_sets_values(), this->special_elements);
		}

		commonality(const powerset_btree<T>& m_focal_elements) : mobius_aggregate<T>(m_focal_elements)
		{
			compute_values_for_mass_focal_elements(this->mass_equivalent.get_focal_sets_values(), this->special_elements);
		}

		commonality(const commonality<T>& q) : commonality(q.get_mass_equivalent().get_focal_sets_values(), q.get_special_elements())
		{}

		commonality(const powerset_btree<T>& m_focal_elements, const powerset_btree<T>& _special_elements) :
			mobius_aggregate<T>(m_focal_elements, _special_elements)
		{}

		commonality(const mobius_aggregate<T>& ma) : commonality(ma.get_mass_equivalent())
		{}

		commonality(const FOD& fod) : mobius_aggregate<T>(fod)
		{
			compute_values_for_mass_focal_elements(this->mass_equivalent.get_focal_sets_values(), this->special_elements);
		}

		commonality(const FOD& fod, const Special_case s_case) : mobius_aggregate<T>(fod, s_case)
		{
			compute_values_for_mass_focal_elements(this->mass_equivalent.get_focal_sets_values(), this->special_elements);
		}

		template <class fusion_rule>
		commonality<T> apply(const fusion_rule fusion, const commonality<T>& q2) const {
			return fusion(*this, q2);
		}

		static void compute_values_for_mass_focal_elements(const powerset_btree<T>& m_focal_elements, powerset_btree<T>& special_elements) {
			EMT_aggregation<T> mobius_computation_scheme(m_focal_elements, order_relation_t::superset);
			const std::vector<set_N_value<T>* >& elements = mobius_computation_scheme.structure.elements();

			switch (mobius_computation_scheme.scheme_type) {
				case scheme_type_t::other:
					std::clog << "other" << std::endl;
					for (size_t i = 0; i < elements.size(); ++i) {
						special_elements.insert(elements[i]->set, compute_aggregation(m_focal_elements, elements[i]->set));
					}
					break;
				case scheme_type_t::consonant:
					std::clog << "consonant" << std::endl;
					for (size_t c = 1; c < ordered_vector.size(); ++c) {
						val = (*ordered_vector[c-1])[0]->value / (*ordered_vector[c])[0]->value;
						special_elements.insert(
							(*ordered_vector[c])[0]->set,
							val
						);
					}
					break;
				case scheme_type_t::semilattice:
					std::clog << "semilattice" << std::endl;
					break;
				case scheme_type_t::lattice:
					std::clog << "lattice" << std::endl;
					break;
				default:
					break;
			}
		}

		static void to_mass_focal_elements(const powerset_btree<T>& q_tree, powerset_btree<T>& m_tree) {

			m_tree.copy(q_tree);
			to_mass_focal_elements_without_initialization(m_tree);
		}

		static void to_mass_focal_elements_without_initialization(powerset_btree<T>& m_tree) {

			std::unordered_map<size_t, std::vector<set_N_value<T>* > > elements_by_cardinality = m_tree.elements_by_set_cardinality();
			const std::vector<size_t>& ordered_vector = m_tree.get_sorted_cardinalities(elements_by_cardinality);

			// computation based on f_elements
			for (size_t c = ordered_vector.size()-1; c > 0; --c) {
				for (size_t i = 0; i < elements_by_cardinality[ordered_vector[c]].size(); ++i) {
					const std::vector<set_N_value<T>* >& subsets = m_tree.strict_subsets_of(elements_by_cardinality[ordered_vector[c]][i]->set);
					for (size_t k = 0; k < subsets.size(); ++k) {
						subsets[k]->value -= elements_by_cardinality[ordered_vector[c]][i]->value;
					}
				}
			}
		}
	};

} // namespace efficient_DST

#endif // EFFICIENT_DST_COMMONALITY_HPP
