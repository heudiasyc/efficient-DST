#ifndef EFFICIENT_DST_RULE_CONJUNCTIVE_HPP
#define EFFICIENT_DST_RULE_CONJUNCTIVE_HPP

#include <rule_classic_general.hpp>
#include <commonality.hpp>
#include <conjunctive_weight.hpp>
#include <mass.hpp>

namespace efficient_DST{

	template <typename T, size_t N>
	class rule_conjunctive : public rule_classic_general<T, N>{
	public:

		std::string to_string() const {
			return "Conjunctive rule";
		}

		mass<T, N> operator()(const mass<T, N>& m1, const mass<T, N>& m2) const {
			const powerset_btree<T, N>& m1_definition = m1.get_definition();
			const powerset_btree<T, N>& m2_definition = m2.get_definition();

			if (m1_definition.size() * m2_definition.size() < 0.5 * m1_definition.get_FOD_size() * pow(2, m1_definition.get_FOD_size())){
				return rule_classic_general<T, N>::operator ()(m1, m2, FOD<N>::set_intersection);
			}else{
				zeta_transform<T, N, up_inclusion<T, N> > q1(m1_definition, operation_type_t::addition);
				zeta_transform<T, N, up_inclusion<T, N> > q2(m2_definition, operation_type_t::addition);
				commonality<T, N> q12 = operator()(*(commonality<T, N>*) &q1, *(commonality<T, N>*) &q2);
				return mass<T, N>(q12);
			}
		}


		commonality<T, N> operator()(const commonality<T, N>& q1, const commonality<T, N>& q2) const {
			const powerset_btree<T, N>& w1_inverted_definition = q1.inversion(operation_type_t::multiplication);
			const powerset_btree<T, N>& w2_inverted_definition = q2.inversion(operation_type_t::multiplication);
			zeta_transform<T, N, up_inclusion<T, N> > q12(
				rule_classic_general<T, N>::weight_fusion(w1_inverted_definition, w2_inverted_definition),
				operation_type_t::multiplication
			);
			return *(commonality<T, N>*) &q12;
		}


		conjunctive_weight<T, N> operator()(const conjunctive_weight<T, N>& w1, const conjunctive_weight<T, N>& w2) const {
			return conjunctive_weight<T, N>(rule_classic_general<T, N>::weight_fusion(w1.get_definition(), w2.get_definition()));
		}
	};

} // namespace efficient_DST

#endif // EFFICIENT_DST_RULE_CONJUNCTIVE_HPP
