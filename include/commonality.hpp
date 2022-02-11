#ifndef EFFICIENT_DST_COMMONALITY_HPP
#define EFFICIENT_DST_COMMONALITY_HPP

#include <mass.hpp>
#include <conjunctive_weight.hpp>
#include <zeta_transform.hpp>

namespace efficient_DST{

	template <size_t N, typename T = float>
	class commonality : public zeta_transform<up_inclusion<N, T>, N, T> {
	public:

		commonality(
			const mass<N, T>& m
		) : zeta_transform<up_inclusion<N, T>, N, T>(m.get_sample_space(), m.get_definition(), m.get_default_value(), operation_type_t::addition)
		{}

		commonality(
			const conjunctive_weight<N, T>& w
		) : zeta_transform<up_inclusion<N, T>, N, T>(w.get_sample_space(), w.get_definition(), w.get_default_value(), operation_type_t::multiplication)
		{}

		commonality(
			const commonality<N, T>& q
		) : zeta_transform<up_inclusion<N, T>, N, T>(q)
		{}

//		commonality(
//			const powerset_btree<N, T>& focal_points_values,
//			const scheme_type_t& scheme_type,
//			const std::vector<subset >& iota_sequence,
//			const T& neutral_value
//		) : zeta_transform<up_inclusion<N, T>, N, T>(
//				focal_points_values,
//				scheme_type,
//				iota_sequence,
//				neutral_value
//			)
//		{}


		template <class fusion_rule>
		commonality<N, T> fuse_with(const commonality<N, T>& q2) const {
			const fusion_rule fusion;
			return fusion(*this, q2);
		}
	};

} // namespace efficient_DST

#endif // EFFICIENT_DST_COMMONALITY_HPP
