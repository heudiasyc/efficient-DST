#ifndef EFFICIENT_DST_COMMONALITY_HPP
#define EFFICIENT_DST_COMMONALITY_HPP

#include <mass.hpp>
#include <conjunctive_weight.hpp>
#include <mobius_aggregate.hpp>

namespace efficient_DST{

	template <typename T = double>
	class commonality : public mobius_aggregate<T> {
	public:

		commonality(const mass<T>& m) : mobius_aggregate<T>(m.get_definition(), order_relation_t::superset, mobius_transformation_form_t::additive)
		{}

		commonality(const conjunctive_weight<T>& w) : mobius_aggregate<T>(w.inverted_definition(), order_relation_t::superset, mobius_transformation_form_t::multiplicative)
		{}

		commonality(const commonality<T>& q) : mobius_aggregate<T>(q)
		{}

		commonality(const powerset_btree<T>& focal_points_values) : mobius_aggregate<T>(focal_points_values, order_relation_t::superset)
		{}

		commonality(const std::vector<T>& powerset_values, const FOD& fod) : mobius_aggregate<T>(powerset_values, fod, order_relation_t::superset)
		{}


		template <class fusion_rule>
		commonality<T> apply(const fusion_rule fusion, const commonality<T>& q2) const {
			return fusion(*this, q2);
		}
	};

} // namespace efficient_DST

#endif // EFFICIENT_DST_COMMONALITY_HPP
