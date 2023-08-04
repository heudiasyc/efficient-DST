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
 
#ifndef EFFICIENT_DST_WEIGHT_VECTOR_HPP
#define EFFICIENT_DST_WEIGHT_VECTOR_HPP

#include <mobius_transform_vector.hpp>
#include <zeta_transform_vector.hpp>

namespace efficient_DST{

	template <size_t N, typename T = float>
	class weight_vector : public mobius_transform_vector<N, T>{
	public:

		weight_vector(const weight_vector<N, T>& w) : mobius_transform_vector<N, T>(w.outcomes, w.definition, 1)
		{}

		weight_vector(
			const sample_space<N>& outcomes,
			const std::vector<T>& log_focal_sets
		) : mobius_transform_vector<N, T>(outcomes, log_focal_sets, 1)
		{
			this->remove_negligible_values();
			this->normalize();
		}

		weight_vector(
			const sample_space<N>& outcomes
		) : mobius_transform_vector<N, T>(outcomes, 1)
		{}

		weight_vector(const zeta_transform_vector<up_inclusion<N, T>, N, T >& q) : weight_vector<N, T>(q.get_sample_space(), q.inversion(operation_type_t::multiplication))
		{}

		weight_vector(const zeta_transform_vector<down_inclusion<N, T>, N, T >& b) : weight_vector<N, T>(b.get_sample_space(), b.inversion(operation_type_t::multiplication))
		{}


		void regularize() {
			this->definition[0] = 1;
			this->normalize();
		}

		void normalize() {
			normalize(this->definition);
		}

		static void normalize(std::vector<T>& definition) {
			T prod = 1;
			for (size_t i = 0; i < definition.size(); ++i) {
				prod *= definition[i];
			}
			if(prod == 0){
				std::cerr << "\nProduct of weight values equal to 0."
						<< "\nThis means that at least one of these weights is null, which is not allowed in a weight function." << std::endl;
				exit(1);
			}
			if(prod != 1){
				// normalize
				T factor = pow(prod, 1/definition.size());
				for (size_t i = 0; i < definition.size(); ++i) {
					definition[i] /= factor;
				}
			}
		}
	};

} // namespace efficient_DST

#endif // EFFICIENT_DST_WEIGHT_VECTOR_HPP
