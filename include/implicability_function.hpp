/*
 * Copyright (C) 2019-2023  Maxime Chaveroche (maxime.chaveroche@gmail.com)
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the CeCILL License, either version 2.1 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * CeCILL License for more details.
 * 
 * You should have received a copy of the CeCILL License
 * along with this program. If not, see <http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html>.
 */
 
#ifndef EFFICIENT_DST_IMPLICABILITY_FUNCTION_HPP
#define EFFICIENT_DST_IMPLICABILITY_FUNCTION_HPP

#include <mass_function.hpp>
#include <disjunctive_decomposition.hpp>
#include <zeta_transform.hpp>

namespace efficient_DST{

	template <size_t N, typename T = float>
	class implicability_function : public zeta_transform<down_inclusion<N, T>, N, T> {
	public:

		implicability_function(
			const mass_function<N, T>& m
		) : zeta_transform<down_inclusion<N, T>, N, T>(m.get_sample_space(), m.get_definition(), m.get_default_value(), operation_type_t::addition)
		{}

		implicability_function(
			const weight_function<N, T>& v
		) : zeta_transform<down_inclusion<N, T>, N, T>(v.get_sample_space(), v.get_definition(), v.get_default_value(), operation_type_t::multiplication)
		{}

		implicability_function(
			const disjunctive_decomposition<N, T>& v
		) : zeta_transform<down_inclusion<N, T>, N, T>(v.get_sample_space(), v.get_definition(), v.get_default_value(), operation_type_t::multiplication)
		{
			this->normalize(v.normalizing_set, v.normalizing_value);
		}

		implicability_function(
			const mass_function<N, T>& m,
			scheme_type_t scheme_type
		) : zeta_transform<down_inclusion<N, T>, N, T>(m.get_sample_space(), m.get_definition(), m.get_default_value(), operation_type_t::addition, scheme_type)
		{}

		implicability_function(
			const weight_function<N, T>& v,
			scheme_type_t scheme_type
		) : zeta_transform<down_inclusion<N, T>, N, T>(v.get_sample_space(), v.get_definition(), v.get_default_value(), operation_type_t::multiplication, scheme_type)
		{}

		implicability_function(
			const disjunctive_decomposition<N, T>& v,
			scheme_type_t scheme_type
		) : zeta_transform<down_inclusion<N, T>, N, T>(v.get_sample_space(), v.get_definition(), v.get_default_value(), operation_type_t::multiplication, scheme_type)
		{
			this->normalize(v.normalizing_set, v.normalizing_value);
		}

		implicability_function(
			const implicability_function<N, T>& b
		) : zeta_transform<down_inclusion<N, T>, N, T>(b)
		{}


		template <class fusion_rule>
		implicability_function<N, T> fuse_with(const implicability_function<N, T>& b2) const {
			const fusion_rule fusion;
			return fusion(*this, b2);
		}

		bool is_a_belief_function(){
			return powerset_function<N, T>::is_equivalent_to_zero(this->at_emptyset());
		}
	};

} // namespace efficient_DST

#endif // EFFICIENT_DST_IMPLICABILITY_FUNCTION_HPP
