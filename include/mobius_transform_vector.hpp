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
 
#ifndef EFFICIENT_DST_MOBIUS_TRANSFORM_VECTOR_HPP
#define EFFICIENT_DST_MOBIUS_TRANSFORM_VECTOR_HPP

#include <powerset_vector.hpp>


namespace efficient_DST{

	template <size_t N, typename T = float>
	class mobius_transform_vector : public powerset_vector<N, T>{
	public:

		mobius_transform_vector(
			const sample_space<N>& outcomes,
			const std::vector<T>& support_values,
			const T& default_value
		) : powerset_vector<N, T>(outcomes, support_values, default_value)
		{
			this->definition = support_values;
		}

		mobius_transform_vector(const sample_space<N>& outcomes, const T& default_value) : powerset_vector<N, T>(outcomes, default_value)
		{
			this->definition = std::vector<T>(1 << N, default_value);
		}

//		mobius_transform(
//			const zeta_transform<T, N, up_inclusion<N, T>, mobius_additive_operation<T> >& z
//		) : powerset_function<N, T>(z.additive_inversion())
//		{}

		void clear() {
			for (size_t i = 0; i < this->definition.size(); ++i){
				this->definition[i] = this->default_value;
			}
		}

		void nullify(const std::vector<std::string>& labels) {
			this->nullify(this->outcomes.get_subset_index(labels));
		}

		void nullify(const size_t& index) {
			this->definition[index] = this->default_value;
		}

		static void remove_negligible_values(std::vector<T>& definition, const T& neutral_value) {
			for (size_t i = 0; i < definition.size(); ++i){
				if(powerset_function<N, T>::is_equivalent_to_zero(definition[i] - neutral_value)){
					definition[i] = neutral_value;
				}
			}
		}

		void remove_negligible_values() {
			remove_negligible_values(this->definition, this->default_value);
		}

		void assign(const std::unordered_map<size_t, T>& values) {
			for (const auto& labels_U_value : values) {
				this->assign(labels_U_value.first, labels_U_value.second);
			}
		}

		void assign(const std::unordered_map<std::vector<std::string>, T>& values) {
			for (const auto& labels_U_value : values) {
				this->assign(labels_U_value.first, labels_U_value.second);
			}
		}

		void assign(const size_t& index, const T& value) {
			this->definition[index] = value;
		}

		void assign(const std::vector<std::string>& labels, const T& value) {
			this->assign(this->outcomes.get_subset_index(labels), value);
		}

		void assign_emptyset(const T& value) {
			this->assign(0, value);
		}

		void assign_fullset(const T& value) {
			this->assign((1 << N)-1, value);
		}
	};

} // namespace efficient_DST

#endif // EFFICIENT_DST_MOBIUS_TRANSFORM_VECTOR_HPP
