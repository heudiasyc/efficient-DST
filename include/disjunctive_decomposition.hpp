#ifndef EFFICIENT_DST_DISJUNCTIVE_DECOMPOSITION_HPP
#define EFFICIENT_DST_DISJUNCTIVE_DECOMPOSITION_HPP

#include <decomposition.hpp>
#include <zeta_transform.hpp>

namespace efficient_DST{

	template <size_t N, typename T = float>
	class disjunctive_decomposition : public decomposition<down_inclusion<N, T>, N, T> {
	public:
		using typename powerset_function<N, T>::subset;
		using powerset_function<N, T>::emptyset;
		using powerset_function<N, T>::fullset;

		disjunctive_decomposition(
			const weight_function<N, T>& w,
			const bool& adaptive_uncertainty = true
		) : decomposition<down_inclusion<N, T>, N, T>(w, adaptive_uncertainty)
		{}

		disjunctive_decomposition(
			const sample_space<N>& outcomes,
			const bool& adaptive_uncertainty = true
		) : decomposition<down_inclusion<N, T>, N, T>(outcomes, adaptive_uncertainty)
		{}

		disjunctive_decomposition(
			const zeta_transform<down_inclusion<N, T>, N, T>& b,
			const bool& adaptive_uncertainty = true
		) : decomposition<down_inclusion<N, T>, N, T>(b, adaptive_uncertainty)
		{}

		template <class fusion_rule>
		disjunctive_decomposition<N, T> fuse_with(const disjunctive_decomposition<N, T>& v2) const {
			const fusion_rule fusion;
			return fusion(*this, v2);
		}
	};

} // namespace efficient_DST

#endif // EFFICIENT_DST_DISJUNCTIVE_DECOMPOSITION_HPP
