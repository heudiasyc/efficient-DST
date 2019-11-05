#ifndef EFFICIENT_DST_COMMONALITY_HPP
#define EFFICIENT_DST_COMMONALITY_HPP

#include <mass.hpp>
#include <conjunctive_weight.hpp>
#include <zeta_transform.hpp>

namespace efficient_DST{

	template <typename T, size_t N>
	class commonality : public zeta_transform<T, N> {
	public:

		commonality(const mass<T, N>& m) : zeta_transform<T, N>(m.get_definition(), order_relation_t::superset, operation_t::addition)
		{}

		commonality(const conjunctive_weight<T, N>& w) : zeta_transform<T, N>(w.inverted_definition(), order_relation_t::superset, operation_t::multiplication)
		{}

		commonality(const commonality<T, N>& q) : zeta_transform<T, N>(q)
		{}

		commonality(const powerset_btree<T, N>& focal_points_values) : zeta_transform<T, N>(focal_points_values, order_relation_t::superset)
		{}

		commonality(const std::vector<T>& powerset_values, FOD<N>& fod) : zeta_transform<T, N>(powerset_values, fod, order_relation_t::superset)
		{}


		template <class fusion_rule>
		commonality<T, N> apply(const commonality<T, N>& q2) const {
			const fusion_rule fusion;
			return fusion(*this, q2);
		}
	};

} // namespace efficient_DST

#endif // EFFICIENT_DST_COMMONALITY_HPP
