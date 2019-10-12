#ifndef EFFICIENT_DST_RULE_CONJUNCTIVE_CAUTIOUS_HPP
#define EFFICIENT_DST_RULE_CONJUNCTIVE_CAUTIOUS_HPP

#include <commonality.hpp>
#include <conjunctive_weight.hpp>
#include <mass.hpp>

namespace efficient_DST{

	template <typename T = double>
	class rule_conjunctive_cautious {
	public:

		std::string to_string() const {
			return "Cautious conjunctive rule";
		}

		mass<T> operator()(const mass<T>& m1, const mass<T>& m2) const {
			zeta_transform<T> q1(m1.get_definition(), order_relation_t::superset, operation_t::addition);
			zeta_transform<T> q2(m2.get_definition(), order_relation_t::superset, operation_t::addition);
			commonality<T> q12 = operator()(*(commonality<T>*) &q1, *(commonality<T>*) &q2);
			return mass<T>(q12);
		}


		commonality<T> operator()(const commonality<T>& q1, const commonality<T>& q2) const {
			const conjunctive_weight<T>& w1(q1);
			const conjunctive_weight<T>& w2(q2);
			return commonality<T>(operator ()(w1, w2));
		}


		conjunctive_weight<T> operator()(const conjunctive_weight<T>& w1, const conjunctive_weight<T>& w2) const {
			return conjunctive_weight<T>(weight_fusion(w1.get_definition(), w2.get_definition()));
		}


	protected:

		static T min(const T& val1, const T& val2){
			return std::min(val1, val2);
		}

		powerset_btree<T> weight_fusion(const powerset_btree<T>& w1_definition, const powerset_btree<T>& w2_definition) const {
			powerset_btree<T> w12_definition(w1_definition.get_FOD(), w1_definition.get_block_size());
			w12_definition.fill_with_union_of_powersets(w1_definition, w2_definition, min, 1);
			decomposition_weight<T>::remove_negligible_values(w12_definition);
			conjunctive_weight<T>::compute_fod_value_from_definition(w12_definition);
			return w12_definition;
		}
	};

} // namespace efficient_DST

#endif // EFFICIENT_DST_RULE_CONJUNCTIVE_CAUTIOUS_HPP
