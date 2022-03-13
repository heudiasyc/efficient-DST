#ifndef EFFICIENT_DST_CONJUNCTIVE_DECOMPOSITION_HPP
#define EFFICIENT_DST_CONJUNCTIVE_DECOMPOSITION_HPP

#include <decomposition.hpp>
#include <zeta_transform.hpp>

namespace efficient_DST{

	template <size_t N, typename T = float>
	class conjunctive_decomposition : public decomposition<up_inclusion<N, T>, N, T> {
	public:
		using typename powerset_function<N, T>::subset;
		using powerset_function<N, T>::emptyset;
		using powerset_function<N, T>::fullset;

		conjunctive_decomposition(
			const weight_function<N, T>& w,
			const bool& adaptive_uncertainty = true
		) : decomposition<up_inclusion<N, T>, N, T>(w, adaptive_uncertainty)
		{}

		conjunctive_decomposition(
			const sample_space<N>& outcomes,
			const bool& adaptive_uncertainty = true
		) : decomposition<up_inclusion<N, T>, N, T>(outcomes, adaptive_uncertainty)
		{}

		conjunctive_decomposition(
			const zeta_transform<up_inclusion<N, T>, N, T>& q,
			const bool& adaptive_uncertainty = true
		) : decomposition<up_inclusion<N, T>, N, T>(q, adaptive_uncertainty)
		{}

		template <class fusion_rule>
		conjunctive_decomposition<N, T> fuse_with(const conjunctive_decomposition<N, T>& w2) const {
			const fusion_rule fusion;
			return fusion(*this, w2);
		}
	};

} // namespace efficient_DST

#endif // EFFICIENT_DST_CONJUNCTIVE_DECOMPOSITION_HPP
