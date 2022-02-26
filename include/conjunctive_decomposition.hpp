#ifndef EFFICIENT_DST_CONJUNCTIVE_DECOMPOSITION_HPP
#define EFFICIENT_DST_CONJUNCTIVE_DECOMPOSITION_HPP

#include <decomposition.hpp>
#include <zeta_transform.hpp>

namespace efficient_DST{

	template <size_t N, typename T = float, bool adaptive_uncertainty = true>
	class conjunctive_decomposition : public decomposition<up_inclusion<N, T>, N, T, adaptive_uncertainty> {
	public:
		using typename powerset_function<N, T>::subset;
		using powerset_function<N, T>::emptyset;
		using powerset_function<N, T>::fullset;

		conjunctive_decomposition(const weight_function<N, T>& w) : decomposition<up_inclusion<N, T>, N, T, adaptive_uncertainty>(w)
		{}

		conjunctive_decomposition(const sample_space<N>& outcomes) : decomposition<up_inclusion<N, T>, N, T, adaptive_uncertainty>(outcomes)
		{}

		conjunctive_decomposition(
			const zeta_transform<up_inclusion<N, T>, N, T>& q
		) : decomposition<up_inclusion<N, T>, N, T, adaptive_uncertainty>(q)
		{}

		template <class fusion_rule>
		conjunctive_decomposition<N, T, adaptive_uncertainty> fuse_with(const conjunctive_decomposition<N, T, adaptive_uncertainty>& w2) const {
			const fusion_rule fusion;
			return fusion(*this, w2);
		}
	};

} // namespace efficient_DST

#endif // EFFICIENT_DST_CONJUNCTIVE_DECOMPOSITION_HPP
