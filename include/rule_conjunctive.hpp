#ifndef EFFICIENT_DST_RULE_CONJUNCTIVE_HPP
#define EFFICIENT_DST_RULE_CONJUNCTIVE_HPP

#include <rule_classic_general.hpp>
#include <commonality.hpp>
#include <conjunctive_weight.hpp>
#include <mass.hpp>

namespace efficient_DST{

	template <typename T = double>
	class rule_conjunctive : public rule_classic_general<T>{
	public:

		std::string to_string() const {
			return "Conjunctive rule";
		}

		mass<T> operator()(const mass<T>& m1, const mass<T>& m2) const {
			const powerset_btree<T>& m1_definition = m1.get_definition();
			const powerset_btree<T>& m2_definition = m2.get_definition();

			if (m1_definition.size() * m2_definition.size() < 0.5 * m1_definition.get_FOD_size() * pow(2, m1_definition.get_FOD_size())){
				return rule_classic_general<T>::operator ()(m1, m2, FOD::set_intersection);
			}else{
				zeta_transform<T> q1(m1_definition, order_relation_t::superset, operation_t::addition);
				zeta_transform<T> q2(m2_definition, order_relation_t::superset, operation_t::addition);
				commonality<T> q12 = operator()(*(commonality<T>*) &q1, *(commonality<T>*) &q2);
				return mass<T>(q12);
			}
		}


		commonality<T> operator()(const commonality<T>& q1, const commonality<T>& q2) const {
			const powerset_btree<T>& w1_inverted_definition = q1.inversion(operation_t::multiplication);
			const powerset_btree<T>& w2_inverted_definition = q2.inversion(operation_t::multiplication);
			zeta_transform<T> q12(
				rule_classic_general<T>::weight_fusion(w1_inverted_definition, w2_inverted_definition),
				order_relation_t::superset,
				operation_t::multiplication
			);
			return *(commonality<T>*) &q12;
		}


		conjunctive_weight<T> operator()(const conjunctive_weight<T>& w1, const conjunctive_weight<T>& w2) const {
			return conjunctive_weight<T>(rule_classic_general<T>::weight_fusion(w1.get_definition(), w2.get_definition()));
		}
	};

} // namespace efficient_DST

#endif // EFFICIENT_DST_RULE_CONJUNCTIVE_HPP
