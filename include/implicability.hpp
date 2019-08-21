#ifndef EFFICIENT_DST_IMPLICABILITY_HPP
#define EFFICIENT_DST_IMPLICABILITY_HPP

#include <mobius_aggregate.hpp>

namespace efficient_DST{

	template <typename T = double>
	class implicability : public mobius_aggregate<T> {
	protected:

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

		implicability(const mass<T>& m) : mobius_aggregate<T>(m)
		{
			compute_values_for_mass_focal_elements(this->mass_equivalent.get_focal_sets_values(), this->special_elements);
		}

		implicability(const powerset_btree<T>& m_focal_elements) : mobius_aggregate<T>(m_focal_elements)
		{
			compute_values_for_mass_focal_elements(this->mass_equivalent.get_focal_sets_values(), this->special_elements);
		}

		implicability(const implicability<T>& b) : implicability(b.get_mass_equivalent().get_focal_sets_values(), b.get_special_elements())
		{}

		implicability(const powerset_btree<T>& m_focal_elements, const powerset_btree<T>& _special_elements) :
			mobius_aggregate<T>(m_focal_elements, _special_elements)
		{}

		implicability(const mobius_aggregate<T>& ma) : implicability(ma.get_mass_equivalent())
		{}

		implicability(const FOD& fod) : mobius_aggregate<T>(fod)
		{
			compute_values_for_mass_focal_elements(this->mass_equivalent.get_focal_sets_values(), this->special_elements);
		}

		implicability(const FOD& fod, const Special_case s_case) : mobius_aggregate<T>(fod, s_case)
		{
			compute_values_for_mass_focal_elements(this->mass_equivalent.get_focal_sets_values(), this->special_elements);
		}


		static T compute_aggregation_at_emptyset(const powerset_btree<T>& m_focal_elements) {
			return m_focal_elements.sub_fod_of_size(0)->value;
		}

		static T compute_aggregation_at_fod(const powerset_btree<T>& m_focal_elements) {
			T sum = 0;
			const std::vector<set_N_value<T>* >& subsets = m_focal_elements.elements();

			for (size_t i = 0; i < subsets.size(); ++i) {
					sum += subsets[i]->value;
			}
			return sum;
		}

		static T compute_aggregation(const powerset_btree<T>& m_focal_elements, const boost::dynamic_bitset<>& set) {
			T sum = 0;
			const std::vector<set_N_value<T>* >& subsets = m_focal_elements.subsets_of(set);

			for (size_t i = 0; i < subsets.size(); ++i) {
					sum += subsets[i]->value;
			}
			return sum;
		}

		static T compute_aggregation(const powerset_btree<T>& m_focal_elements, const std::vector<fod_element*>& fod_elements) {
			return compute_aggregation(m_focal_elements, m_focal_elements.fod->to_set(fod_elements));
		}


		template <class fusion_rule>
		implicability<T> apply(const fusion_rule fusion, const implicability<T>& i2) const {
			return fusion(*this, i2);
		}

		static void compute_values_for_mass_focal_elements(const powerset_btree<T>& m_focal_elements, powerset_btree<T>& special_elements) {
			const std::vector<set_N_value<T>* >& elements = m_focal_elements.elements();
			// pre-calculation for all focal elements
			for (size_t i = 0; i < elements.size(); ++i) {
				special_elements.insert(elements[i]->set, compute_aggregation(m_focal_elements, elements[i]->set));
			}
		}

		static void to_mass_focal_elements(const powerset_btree<T>& bel_tree, powerset_btree<T>& m_tree) {

			m_tree.copy(bel_tree);

			std::unordered_map<size_t, std::vector<set_N_value<T>* > > elements_by_cardinality = m_tree.elements_by_set_cardinality();
			const std::vector<size_t>& ordered_vector = m_tree.get_sorted_cardinalities(elements_by_cardinality);

			// computation based on f_elements
			for (size_t c = 0; c < ordered_vector.size()-1; ++c) {
				for (size_t i = 0; i < elements_by_cardinality[ordered_vector[c]].size(); ++i) {
					const std::vector<set_N_value<T>* >& supersets = m_tree.strict_supersets_of(elements_by_cardinality[ordered_vector[c]][i]->set);
					for (size_t k = 0; k < supersets.size(); ++k) {
						supersets[k]->value -= elements_by_cardinality[ordered_vector[c]][i]->value;
					}
				}
			}
		}
	};

} // namespace efficient_DST

#endif // EFFICIENT_DST_IMPLICABILITY_HPP
