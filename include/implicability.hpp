#ifndef OW_BFT_IMPLICABILITY_HPP
#define OW_BFT_IMPLICABILITY_HPP

#include <belief.hpp>

namespace ow_bft{

	template <typename T = double>
	class implicability : public belief<T> {
	protected:

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

		static T compute_aggregation(const powerset_btree<T>& m_focal_elements, const boost::dynamic_bitset<>& key) {
			T sum = 0;
			const std::vector<set_N_value<T>* >& subsets = m_focal_elements.subsets_of(key);

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

		implicability(const mass<T>& m) : mass_aggregate<T>(m)
		{
			compute_values_for_special_elements(this->mass_equivalent.get_focal_elements(), this->special_elements);
		}

		implicability(const mass_aggregate<T>& ma) : implicability(ma.get_mass_equivalent())
		{}

		implicability(const implicability<T>& b) : mass_aggregate<T>(b.get_mass_equivalent())
		{
			this->special_elements.copy(b.get_special_elements());
		}

		implicability(const FOD& fod) : mass_aggregate<T>(fod)
		{
			compute_values_for_special_elements(this->mass_equivalent.get_focal_elements(), this->special_elements);
		}

		implicability(const FOD& fod, const Special_case s_case) : mass_aggregate<T>(fod, s_case)
		{
			compute_values_for_special_elements(this->mass_equivalent.get_focal_elements(), this->special_elements);
		}

		template <class fusion_rule>
		implicability<T> apply(const fusion_rule fusion, const implicability<T>& i2) const {
			return fusion(*this, i2);
		}

		static void compute_values_for_mass_focal_elements(const powerset_btree<T>& m_focal_elements, powerset_btree<T>& special_elements) {
			const std::vector<set_N_value<T>* >& elements = m_focal_elements.elements();
			// pre-calculation for all focal elements
			for (size_t i = 0; i < elements.size(); ++i) {
				special_elements.insert(elements[i]->fod_elements, compute_aggregation(m_focal_elements, elements[i]->fod_elements));
			}
		}

		static void to_mass_focal_elements(const powerset_btree<T>& bel_tree, powerset_btree<T>& m_tree) {
			belief<T>::to_mass_focal_elements(bel_tree, m_tree);
		}
	};

} // namespace ow_bft

#endif // OW_BFT_IMPLICABILITY_HPP
