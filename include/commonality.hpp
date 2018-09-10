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

		static T compute_aggregation(const powerset_btree<T>& m_focal_elements, const boost::dynamic_bitset<>& key) {
			T sum = 0;
			const std::vector<set_N_value<T>* >& supersets = m_focal_elements.supersets_of(key);

			if(supersets.size() == m_focal_elements.size())
				return 1;

			for (size_t i = 0; i < supersets.size(); ++i) {
				sum += supersets[i]->value;
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

		/*static powerset_btree<T> to_mass_focal_elements(const powerset_btree<T>& q_tree){
			return powerset_btree<T>(q_tree);
		}*/

		static void compute_values_for_mass_focal_elements(const powerset_btree<T>& m_focal_elements, powerset_btree<T>& special_elements) {
			const std::vector<set_N_value<T>* >& elements = m_focal_elements.elements();
			// pre-calculation for all focal elements
			for (size_t i = 0; i < elements.size(); ++i) {
				special_elements.insert(elements[i]->fod_elements, compute_aggregation(m_focal_elements, elements[i]->fod_elements));
			}
		}

		static void to_mass_focal_elements(const powerset_btree<T>& q_tree, powerset_btree<T>& m_tree) {

			m_tree.copy(q_tree);
			to_mass_focal_elements_without_initialization(m_tree);
		}

		static void to_mass_focal_elements_without_initialization(powerset_btree<T>& m_tree) {

			const std::vector<std::vector<set_N_value<T>* > >& elements_by_cardinality = m_tree.elements_by_set_cardinality();

			// computation based on f_elements
			for (size_t c = elements_by_cardinality.size()-1; c > 0; --c) {
				for (size_t i = 0; i < elements_by_cardinality[c].size(); ++i) {
					const std::vector<set_N_value<T>* >& subsets = m_tree.strict_subsets_of(elements_by_cardinality[c][i]->fod_elements);
					for (size_t k = 0; k < subsets.size(); ++k) {
						subsets[k]->value -= elements_by_cardinality[c][i]->value;
					}
				}
			}
		}

/*
	public:
		void apply_changes(){
			ow_bft::powerset_btree buffer = new ow_bft::powerset_btree(&this->fod);
			std::vector<std::vector<set_N_value<T>* > > changed_q = this->changes.nodes_by_depth;
			std::unordered_map<set_N_value<T>*, std::vector<set_N_value<T>* > > affected_masses;

			for (int i = 0; i < changed_q.size(); ++i) {
				for (int j = 0; j < changed_q[i].size(); ++j) {

					boost::dynamic_bitset key = this->fod.to_set(changed_q[i][j]->fod_elements);
					std::vector<set_N_value<T>* > masses_affected_by_q = this->fod.powerset.non_null_values_N_elements_of_subsets_of(key);

					for (int k = 0; k < masses_affected_by_q.size(); ++k) {
						std::vector<set_N_value<T>* > m = affected_masses[masses_affected_by_q[k]];
						m.push_back(changed_q[i][j]);
					}
				}
			}
			for (std::pair<set_N_value<T>*, std::vector<set_N_value<T>* > > mass_N_qs : affected_masses){
				T sum = 0;
				for (int i = 0; i < mass_N_qs.second.size(); ++i) {
					int card_diff = this->fod.set_minus(mass_N_qs.second[i]->fod_elements, mass_N_qs.first->fod_elements).size();
					if(card_diff % 2 == 0){
						sum += mass_N_qs.second[i]->value - this->before_changes[mass_N_qs.second[i]->fod_elements];
					}else{
						sum -= mass_N_qs.second[i]->value - this->before_changes[mass_N_qs.second[i]->fod_elements];
					}
				}
				buffer.insert(mass_N_qs.first->fod_elements, mass_N_qs.first->value + sum);
			}

			this->fod.powerset = buffer;
		}
*/
	};

} // namespace ow_bft

#endif // OW_BFT_COMMONALITY_HPP
