#ifndef EFFICIENT_DST_CONJUNCTIVE_DECOMPOSITION_VECTOR_HPP
#define EFFICIENT_DST_CONJUNCTIVE_DECOMPOSITION_VECTOR_HPP

#include <decomposition_vector.hpp>
#include <zeta_transform_vector.hpp>

namespace efficient_DST{

	template <size_t N, typename T = float>
	class conjunctive_decomposition_vector : public decomposition_vector<up_inclusion<N, T>, N, T> {
	public:

		conjunctive_decomposition_vector(
			const weight_vector<N, T>& w,
			const bool& adaptive_uncertainty = true
		) : decomposition_vector<up_inclusion<N, T>, N, T>(w, adaptive_uncertainty)
		{}

		conjunctive_decomposition_vector(
			const sample_space<N>& outcomes,
			const bool& adaptive_uncertainty = true
		) : decomposition_vector<up_inclusion<N, T>, N, T>(outcomes, adaptive_uncertainty)
		{}

		conjunctive_decomposition_vector(
			const zeta_transform_vector<up_inclusion<N, T>, N, T>& q,
			const bool& adaptive_uncertainty = true
		) : decomposition_vector<up_inclusion<N, T>, N, T>(q, q.reduced_core(q.get_definition()), 0, adaptive_uncertainty)
		{}


		template <class fusion_rule>
		conjunctive_decomposition_vector<N, T> fuse_with(const conjunctive_decomposition_vector<N, T>& w2) const {
			const fusion_rule fusion;
			return fusion(*this, w2);
		}
	};

} // namespace efficient_DST

#endif // EFFICIENT_DST_CONJUNCTIVE_DECOMPOSITION_VECTOR_HPP
