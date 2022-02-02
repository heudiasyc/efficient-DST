#ifndef EFFICIENT_DST_COMMONALITY_HPP
#define EFFICIENT_DST_COMMONALITY_HPP

#include <mass.hpp>
#include <conjunctive_weight.hpp>
#include <zeta_transform.hpp>

namespace efficient_DST{

	template <typename T, size_t N>
	class commonality : public zeta_transform<T, N, up_inclusion<T, N> > {
	public:

		commonality(const mass<T, N>& m) : zeta_transform<T, N, up_inclusion<T, N> >(m.get_definition(), operation_type_t::addition)
		{}

		commonality(const conjunctive_weight<T, N>& w) : zeta_transform<T, N, up_inclusion<T, N> >(w.inverted_definition(), operation_type_t::multiplication)
		{}

		commonality(const commonality<T, N>& q) : zeta_transform<T, N, up_inclusion<T, N> >(q)
		{}

		commonality(
			const powerset_btree<T, N>& focal_points_values,
			const scheme_type_t& scheme_type,
			const std::vector<std::bitset<N> >& iota_sequence,
			const T& neutral_value
		) : zeta_transform<T, N, up_inclusion<T, N> >(
				focal_points_values,
				scheme_type,
				iota_sequence,
				neutral_value
			)
		{}


		template <class fusion_rule>
		commonality<T, N> apply(const commonality<T, N>& q2) const {
			const fusion_rule fusion;
			return fusion(*this, q2);
		}
	};

} // namespace efficient_DST

#endif // EFFICIENT_DST_COMMONALITY_HPP
