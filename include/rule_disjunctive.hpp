#ifndef EFFICIENT_DST_RULE_DISJUNCTIVE_HPP
#define EFFICIENT_DST_RULE_DISJUNCTIVE_HPP

#include <rule_classic_general.hpp>
#include <implicability.hpp>
#include <disjunctive_weight.hpp>
#include <mass.hpp>

namespace efficient_DST{

	template <typename T = double>
	class rule_disjunctive : public rule_classic_general<T> {
	public:

		std::string to_string() const {
			return "Disjunctive rule";
		}

		mass<T> operator()(const mass<T>& m1, const mass<T>& m2) const {
			const powerset_btree<T>& m1_definition = m1.get_definition();
			const powerset_btree<T>& m2_definition = m2.get_definition();

			if (m1_definition.size() * m2_definition.size() < 0.5 * m1_definition.get_FOD_size() * pow(2, m1_definition.get_FOD_size())){
				return rule_classic_general<T>::operator ()(m1, m2, FOD::set_union);
			}else{
				mobius_aggregate<T> b1(m1_definition, order_relation_t::subset, mobius_transformation_form_t::additive);
				mobius_aggregate<T> b2(m2_definition, order_relation_t::subset, mobius_transformation_form_t::additive);
				implicability<T> b12 = operator()(*(implicability<T>*) &b1, *(implicability<T>*) &b2);
				return mass<T>(b12);
			}
		}


		implicability<T> operator()(const implicability<T>& b1, const implicability<T>& b2) const {
			const powerset_btree<T>& v1_inverted_definition = b1.inversion(mobius_transformation_form_t::multiplicative);
			const powerset_btree<T>& v2_inverted_definition = b2.inversion(mobius_transformation_form_t::multiplicative);
			mobius_aggregate<T> b12(
				rule_classic_general<T>::weight_fusion(v1_inverted_definition, v2_inverted_definition),
				order_relation_t::subset,
				mobius_transformation_form_t::multiplicative
			);
			return *(implicability<T>*) &b12;
		}


		disjunctive_weight<T> operator()(const disjunctive_weight<T>& v1, const disjunctive_weight<T>& v2) const {
			return disjunctive_weight<T>(rule_classic_general<T>::weight_fusion(v1.get_definition(), v2.get_definition()));
		}
	};

} // namespace efficient_DST

#endif // EFFICIENT_DST_RULE_DISJUNCTIVE_HPP
