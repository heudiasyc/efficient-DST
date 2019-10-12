#ifndef EFFICIENT_DST_RULE_DISJUNCTIVE_BOLD_HPP
#define EFFICIENT_DST_RULE_DISJUNCTIVE_BOLD_HPP

#include <implicability.hpp>
#include <disjunctive_weight.hpp>
#include <mass.hpp>

namespace efficient_DST{

	template <typename T = double>
	class rule_disjunctive_bold {
	public:

		std::string to_string() const {
			return "Bold disjunctive rule";
		}

		mass<T> operator()(const mass<T>& m1, const mass<T>& m2) const {
			zeta_transform<T> b1(m1.get_definition(), order_relation_t::subset, operation_t::addition);
			zeta_transform<T> b2(m2.get_definition(), order_relation_t::subset, operation_t::addition);
			implicability<T> b12 = operator()(*(implicability<T>*) &b1, *(implicability<T>*) &b2);
			return mass<T>(b12);
		}


		implicability<T> operator()(const implicability<T>& b1, const implicability<T>& b2) const {
			const disjunctive_weight<T>& v1(b1);
			const disjunctive_weight<T>& v2(b2);
			return implicability<T>(operator ()(v1, v2));
		}


		disjunctive_weight<T> operator()(const disjunctive_weight<T>& v1, const disjunctive_weight<T>& v2) const {
			return disjunctive_weight<T>(weight_fusion(v1.get_definition(), v2.get_definition()));
		}


	protected:

		static T max(const T& val1, const T& val2){
			return std::max(val1, val2);
		}

		powerset_btree<T> weight_fusion(const powerset_btree<T>& v1_definition, const powerset_btree<T>& v2_definition) const {
			powerset_btree<T> v12_definition(v1_definition.get_FOD(), v1_definition.get_block_size());
			v12_definition.fill_with_union_of_powersets(v1_definition, v2_definition, max, 1);
			decomposition_weight<T>::remove_negligible_values(v12_definition);
			disjunctive_weight<T>::compute_emptyset_value_from_definition(v12_definition);
			return v12_definition;
		}
	};

} // namespace efficient_DST

#endif // EFFICIENT_DST_RULE_DISJUNCTIVE_BOLD_HPP
