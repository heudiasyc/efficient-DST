#ifndef EFFICIENT_DST_COMMONALITY_HPP
#define EFFICIENT_DST_COMMONALITY_HPP

#include <mass.hpp>
#include <conjunctive_weight.hpp>
#include <zeta_transform.hpp>

namespace efficient_DST{

	template <typename T = double>
	class commonality : public zeta_transform<T> {
	public:

		commonality(const mass<T>& m) : zeta_transform<T>(m.get_definition(), order_relation_t::superset, operation_t::addition)
		{}

		commonality(const conjunctive_weight<T>& w) : zeta_transform<T>(w.inverted_definition(), order_relation_t::superset, operation_t::multiplication)
		{}

		commonality(const commonality<T>& q) : zeta_transform<T>(q)
		{}

		commonality(const powerset_btree<T>& focal_points_values) : zeta_transform<T>(focal_points_values, order_relation_t::superset)
		{}

		commonality(const std::vector<T>& powerset_values, FOD& fod) : zeta_transform<T>(powerset_values, fod, order_relation_t::superset)
		{}


		template <class fusion_rule>
		commonality<T> apply(const commonality<T>& q2) const {
			const fusion_rule fusion;
			return fusion(*this, q2);
		}
	};

} // namespace efficient_DST

#endif // EFFICIENT_DST_COMMONALITY_HPP
