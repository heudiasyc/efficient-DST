#ifndef EFFICIENT_DST_IMPLICABILITY_HPP
#define EFFICIENT_DST_IMPLICABILITY_HPP

#include <mass.hpp>
#include <disjunctive_weight.hpp>
#include <zeta_transform.hpp>

namespace efficient_DST{

	template <typename T, size_t N>
	class implicability : public zeta_transform<T, N, down_inclusion<T, N> > {
	public:

		implicability(const mass<T, N>& m) : zeta_transform<T, N, down_inclusion<T, N> >(m.get_definition(), operation_type_t::addition)
		{}

		implicability(const disjunctive_weight<T, N>& v) : zeta_transform<T, N, down_inclusion<T, N> >(v.inverted_definition(), operation_type_t::multiplication)
		{}

		implicability(const implicability<T, N>& b) : zeta_transform<T, N, down_inclusion<T, N> >(b)
		{}

		implicability(
			const powerset_btree<T, N>& focal_points_values,
			const scheme_type_t& scheme_type,
			const std::vector<std::bitset<N> >& iota_sequence,
			const T& neutral_value
		) : zeta_transform<T, N, down_inclusion<T, N> >(
				focal_points_values,
				scheme_type,
				iota_sequence,
				neutral_value
			)
		{}


		template <class fusion_rule>
		implicability<T, N> apply(const implicability<T, N>& b2) const {
			const fusion_rule fusion;
			return fusion(*this, b2);
		}
	};

} // namespace efficient_DST

#endif // EFFICIENT_DST_IMPLICABILITY_HPP
