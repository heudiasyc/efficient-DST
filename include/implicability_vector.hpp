#ifndef EFFICIENT_DST_IMPLICABILITY_VECTOR_HPP
#define EFFICIENT_DST_IMPLICABILITY_VECTOR_HPP

#include <mass_vector.hpp>
#include <disjunctive_decomposition_vector.hpp>
#include <zeta_transform_vector.hpp>

namespace efficient_DST{

	template <size_t N, typename T = float>
	class implicability_vector : public zeta_transform_vector<down_inclusion<N, T>, N, T> {
	public:

		implicability_vector(
			const mass_vector<N, T>& m,
			const bool& core_reduced = true
		) : zeta_transform_vector<down_inclusion<N, T>, N, T>(m.get_sample_space(), m.get_definition(), m.get_default_value(), operation_type_t::addition, core_reduced)
		{}

		implicability_vector(
			const weight_vector<N, T>& v,
			const bool& core_reduced = true
		) : zeta_transform_vector<down_inclusion<N, T>, N, T>(v.get_sample_space(), v.get_definition(), v.get_default_value(), operation_type_t::multiplication, core_reduced)
		{}

		implicability_vector(
			const disjunctive_decomposition_vector<N, T>& v,
			const bool& core_reduced = true
		) : zeta_transform_vector<down_inclusion<N, T>, N, T>(
				v.get_sample_space(), v.get_definition(), v.get_default_value(),
				operation_type_t::multiplication, v.has_adaptive_uncertainty() ? false : core_reduced)
		{
			this->normalize(0, v.normalizing_value);
			this->zero = v.normalizing_set;
		}

		implicability_vector(
			const implicability_vector<N, T>& b
		) : zeta_transform_vector<down_inclusion<N, T>, N, T>(b)
		{}


		template <class fusion_rule>
		implicability_vector<N, T> fuse_with(const implicability_vector<N, T>& b2) const {
			const fusion_rule fusion;
			return fusion(*this, b2);
		}

		bool is_a_belief_function(){
			return powerset_function<N, T>::is_equivalent_to_zero(this->at_emptyset());
		}
	};

} // namespace efficient_DST

#endif // EFFICIENT_DST_IMPLICABILITY_VECTOR_HPP
