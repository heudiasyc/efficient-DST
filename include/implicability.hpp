#ifndef EFFICIENT_DST_IMPLICABILITY_HPP
#define EFFICIENT_DST_IMPLICABILITY_HPP

#include <mass.hpp>
#include <disjunctive_weight.hpp>
#include <zeta_transform.hpp>

namespace efficient_DST{

	template <size_t N, typename T = float>
	class implicability : public zeta_transform<down_inclusion<N, T>, N, T> {
	public:

		implicability(
			const mass<N, T>& m
		) : zeta_transform<down_inclusion<N, T>, N, T>(m.get_sample_space(), m.get_definition(), m.get_default_value(), operation_type_t::addition)
		{}

		implicability(
			const disjunctive_weight<N, T>& v
		) : zeta_transform<down_inclusion<N, T>, N, T>(v.get_sample_space(), v.get_definition(), v.get_default_value(), operation_type_t::multiplication)
		{}

		implicability(
			const implicability<N, T>& b
		) : zeta_transform<down_inclusion<N, T>, N, T>(b)
		{}

//		implicability(
//			const powerset_btree<N, T>& focal_points_values,
//			const scheme_type_t& scheme_type,
//			const std::vector<std::bitset<N> >& iota_sequence,
//			const T& neutral_value
//		) : zeta_transform<T, N, down_inclusion<N, T> >(
//				focal_points_values,
//				scheme_type,
//				iota_sequence,
//				neutral_value
//			)
//		{}


		template <class fusion_rule>
		implicability<N, T> fuse_with(const implicability<N, T>& b2) const {
			const fusion_rule fusion;
			return fusion(*this, b2);
		}
	};

} // namespace efficient_DST

#endif // EFFICIENT_DST_IMPLICABILITY_HPP
