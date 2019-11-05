#ifndef EFFICIENT_DST_IMPLICABILITY_HPP
#define EFFICIENT_DST_IMPLICABILITY_HPP

#include <mass.hpp>
#include <disjunctive_weight.hpp>
#include <zeta_transform.hpp>

namespace efficient_DST{

	template <typename T, size_t N>
	class implicability : public zeta_transform<T, N> {
	public:

		implicability(const mass<T, N>& m) : zeta_transform<T, N>(m.get_definition(), order_relation_t::subset, operation_t::addition)
		{}

		implicability(const disjunctive_weight<T, N>& v) : zeta_transform<T, N>(v.inverted_definition(), order_relation_t::subset, operation_t::multiplication)
		{}

		implicability(const implicability<T, N>& b) : zeta_transform<T, N>(b)
		{}

		implicability(const powerset_btree<T, N>& focal_points_values) : zeta_transform<T, N>(focal_points_values, order_relation_t::subset)
		{}

		implicability(const std::vector<T>& powerset_values, FOD<N>& fod) : zeta_transform<T, N>(powerset_values, fod, order_relation_t::subset)
		{}


		template <class fusion_rule>
		implicability<T, N> apply(const implicability<T, N>& b2) const {
			const fusion_rule fusion;
			return fusion(*this, b2);
		}
	};

} // namespace efficient_DST

#endif // EFFICIENT_DST_IMPLICABILITY_HPP
