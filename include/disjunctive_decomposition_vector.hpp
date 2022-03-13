#ifndef EFFICIENT_DST_DISJUNCTIVE_DECOMPOSITION_VECTOR_HPP
#define EFFICIENT_DST_DISJUNCTIVE_DECOMPOSITION_VECTOR_HPP

#include <decomposition_vector.hpp>
#include <zeta_transform.hpp>

namespace efficient_DST{

	template <size_t N, typename T = float>
	class disjunctive_decomposition_vector : public decomposition_vector<down_inclusion<N, T>, N, T> {
	public:

		disjunctive_decomposition_vector(
			const weight_vector<N, T>& w,
			const bool& adaptive_uncertainty = true
		) : decomposition_vector<down_inclusion<N, T>, N, T>(w, adaptive_uncertainty)
		{}

		disjunctive_decomposition_vector(
			const sample_space<N>& outcomes,
			const bool& adaptive_uncertainty = true
		) : decomposition_vector<down_inclusion<N, T>, N, T>(outcomes, adaptive_uncertainty)
		{}

		disjunctive_decomposition_vector(
			const zeta_transform<down_inclusion<N, T>, N, T>& b,
			const bool& adaptive_uncertainty = true
		) : decomposition_vector<down_inclusion<N, T>, N, T>(b, (1 << N)-1, b.augmented_zero(b.get_definition()), adaptive_uncertainty)
		{}


		template <class fusion_rule>
		disjunctive_decomposition_vector<N, T> fuse_with(const disjunctive_decomposition_vector<N, T>& v2) const {
			const fusion_rule fusion;
			return fusion(*this, v2);
		}
	};

} // namespace efficient_DST

#endif // EFFICIENT_DST_DISJUNCTIVE_DECOMPOSITION_VECTOR_HPP
