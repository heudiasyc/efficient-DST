#ifndef EFFICIENT_DST_COMMONALITY_VECTOR_HPP
#define EFFICIENT_DST_COMMONALITY_VECTOR_HPP

#include <mass_vector.hpp>
#include <conjunctive_decomposition_vector.hpp>
#include <zeta_transform_vector.hpp>

namespace efficient_DST{

	template <size_t N, typename T = float>
	class commonality_vector : public zeta_transform_vector<up_inclusion<N, T>, N, T> {
	public:

		commonality_vector(
			const mass_vector<N, T>& m,
			const bool& core_reduced = true
		) : zeta_transform_vector<up_inclusion<N, T>, N, T>(m.get_sample_space(), m.get_definition(), m.get_default_value(), operation_type_t::addition, core_reduced)
		{}

		commonality_vector(
			const weight_vector<N, T>& w,
			const bool& core_reduced = true
		) : zeta_transform_vector<up_inclusion<N, T>, N, T>(w.get_sample_space(), w.get_definition(), w.get_default_value(), operation_type_t::multiplication, core_reduced)
		{}

		commonality_vector(
			const conjunctive_decomposition_vector<N, T>& w,
			const bool& core_reduced = true
		) : zeta_transform_vector<up_inclusion<N, T>, N, T>(
				w.get_sample_space(), w.get_definition(), w.get_default_value(),
				operation_type_t::multiplication, w.has_adaptive_uncertainty() ? false : core_reduced)
		{
//			std::cout << "FMT vector commonality:\n";
//			for (size_t i = 1; i < this->definition.size(); ++i){
//				std::cout << i << " : " << this->definition[i] << "\n";
//			}
			this->normalize(this->definition.size()-1, w.normalizing_value);
			this->core = w.normalizing_set;
		}

		commonality_vector(
			const commonality_vector<N, T>& q
		) : zeta_transform_vector<up_inclusion<N, T>, N, T>(q)
		{}


		template <class fusion_rule>
		commonality_vector<N, T> fuse_with(const commonality_vector<N, T>& q2) const {
			const fusion_rule fusion;
			return fusion(*this, q2);
		}
	};

} // namespace efficient_DST

#endif // EFFICIENT_DST_COMMONALITY_VECTOR_HPP
