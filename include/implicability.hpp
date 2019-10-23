#ifndef EFFICIENT_DST_IMPLICABILITY_HPP
#define EFFICIENT_DST_IMPLICABILITY_HPP

#include <mass.hpp>
#include <disjunctive_weight.hpp>
#include <zeta_transform.hpp>

namespace efficient_DST{

	template <typename T = double>
	class implicability : public zeta_transform<T> {
	public:

		implicability(const mass<T>& m) : zeta_transform<T>(m.get_definition(), order_relation_t::subset, operation_t::addition)
		{}

		implicability(const disjunctive_weight<T>& v) : zeta_transform<T>(v.inverted_definition(), order_relation_t::subset, operation_t::multiplication)
		{}

		implicability(const implicability<T>& b) : zeta_transform<T>(b)
		{}

		implicability(const powerset_btree<T>& focal_points_values) : zeta_transform<T>(focal_points_values, order_relation_t::subset)
		{}

		implicability(const std::vector<T>& powerset_values, FOD& fod) : zeta_transform<T>(powerset_values, fod, order_relation_t::subset)
		{}


		template <class fusion_rule>
		implicability<T> apply(const implicability<T>& b2) const {
			const fusion_rule fusion;
			return fusion(*this, b2);
		}
	};

} // namespace efficient_DST

#endif // EFFICIENT_DST_IMPLICABILITY_HPP
