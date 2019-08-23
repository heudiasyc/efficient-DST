#ifndef EFFICIENT_DST_IMPLICABILITY_HPP
#define EFFICIENT_DST_IMPLICABILITY_HPP

#include <mass.hpp>
#include <disjunctive_weight.hpp>
#include <mobius_aggregate.hpp>

namespace efficient_DST{

	template <typename T = double>
	class implicability : public mobius_aggregate<T> {
	public:

		implicability(const mass<T>& m) : mobius_aggregate<T>(m.get_definition(), order_relation_t::subset, mobius_transformation_form_t::additive)
		{}

		implicability(const disjunctive_weight<T>& v) : mobius_aggregate<T>(v.inverted_definition(), order_relation_t::subset, mobius_transformation_form_t::multiplicative)
		{}

		implicability(const implicability<T>& b) : mobius_aggregate<T>(b)
		{}

		implicability(const powerset_btree<T>& focal_points_values) : mobius_aggregate<T>(focal_points_values, order_relation_t::subset)
		{}

		implicability(const std::vector<T>& powerset_values, const FOD& fod) : mobius_aggregate<T>(powerset_values, fod, order_relation_t::subset)
		{}


		template <class fusion_rule>
		implicability<T> apply(const fusion_rule fusion, const implicability<T>& b2) const {
			return fusion(*this, b2);
		}
	};

} // namespace efficient_DST

#endif // EFFICIENT_DST_IMPLICABILITY_HPP
