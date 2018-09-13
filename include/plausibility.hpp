#ifndef OW_BFT_PLAUSIBILITY_HPP
#define OW_BFT_PLAUSIBILITY_HPP

#include <mass_aggregate.hpp>
#include <implicability.hpp>

namespace ow_bft{

	template <typename T = double>
	class plausibility : public mass_aggregate<T> {
	protected:

		static T compute_aggregation_at_emptyset(const powerset_btree<T>& m_focal_elements) {
			return 0;
		}

		static T compute_aggregation_at_fod(const powerset_btree<T>& m_focal_elements) {
			return belief<T>::compute_aggregation_at_fod(m_focal_elements);
		}

		static T compute_aggregation(const powerset_btree<T>& m_focal_elements, const boost::dynamic_bitset<>& key) {
			return 1 - implicability<T>::compute_aggregation(m_focal_elements, m_focal_elements.fod->set_negate(key));
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

		plausibility(const mass<T>& m) : mass_aggregate<T>(m)
		{
			//this->set_values_for_special_elements();
		}

		plausibility(const mass_aggregate<T>& ma) : plausibility(ma.get_mass_equivalent())
		{}

		plausibility(const plausibility<T>& p) : mass_aggregate<T>(p.get_mass_equivalent())
		{
			this->special_elements.copy(p.get_special_elements());
		}

		plausibility(const FOD& fod) : mass_aggregate<T>(fod)
		{
			//this->set_values_for_special_elements();
		}

		plausibility(const FOD& fod, const Special_case s_case) : mass_aggregate<T>(fod, s_case)
		{
			//this->set_values_for_special_elements();
		}

		template <class fusion_rule>
		plausibility<T> apply(const fusion_rule fusion, const plausibility<T>& p2) const {
			return fusion(*this, p2);
		}

		static void compute_values_for_mass_focal_elements(const powerset_btree<T>& m_focal_elements, powerset_btree<T>& special_elements) {
			const std::vector<set_N_value<T>* >& elements = m_focal_elements.elements();
			// pre-calculation for all focal elements
			for (size_t i = 0; i < elements.size(); ++i) {
				special_elements.insert(elements[i]->fod_elements, compute_aggregation(m_focal_elements, elements[i]->fod_elements));
			}
		}

		static void to_mass_focal_elements(const powerset_btree<T>& pl_tree, powerset_btree<T>& m_tree) {
			std::cerr << "\nUnimplemented method to_mass_focal_elements of plausibility. Impossible to infer masses from plausibility on focal elements only.\n"
						<< "But you can create an implicability function based on this plausibility function that will allow you\n"
						<< "to retrieve the corresponding mass function.\n";
		}
	};

} // namespace ow_bft

#endif // OW_BFT_PLAUSIBILITY_HPP

