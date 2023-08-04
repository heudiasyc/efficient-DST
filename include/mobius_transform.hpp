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

#ifndef EFFICIENT_DST_MOBIUS_TRANSFORM_HPP
#define EFFICIENT_DST_MOBIUS_TRANSFORM_HPP

#include <powerset_function.hpp>


namespace efficient_DST{

	template <size_t N, typename T = float>
	class mobius_transform : public powerset_function<N, T>{
	public:
		using typename powerset_function<N, T>::subset;
		using powerset_function<N, T>::emptyset;
		using powerset_function<N, T>::fullset;
//		using powerset_function<N, T>::operator[];

		mobius_transform(
			const sample_space<N>& outcomes,
			const powerset_btree<N, T>& support_values,
			const T& default_value
		) : powerset_function<N, T>(outcomes, support_values, default_value)
		{}

		mobius_transform(
			const sample_space<N>& outcomes,
			const std::vector<T>& support_values,
			const T& default_value
		) : powerset_function<N, T>(outcomes, default_value)
		{
			for (size_t i = 0; i < support_values.size(); ++i){
				if(!powerset_function<N, T>::is_equivalent_to_zero(support_values[i] - default_value)){
					this->assign((subset) i, support_values[i]);
				}
			}
		}

		mobius_transform(const sample_space<N>& outcomes, const T& default_value) : powerset_function<N, T>(outcomes, default_value)
		{}

//		mobius_transform(
//			const zeta_transform<T, N, up_inclusion<N, T>, mobius_additive_operation<T> >& z
//		) : powerset_function<N, T>(z.additive_inversion())
//		{}

		powerset_btree<N, T> definition_(){
			return this->definition;
		}

		void clear() {
			this->definition.nullify();
		}

		void nullify(const std::vector<std::string>& labels) {
			this->nullify(this->outcomes.get_subset(labels));
		}

		void nullify(const subset& set) {
			this->definition.nullify(this->definition[set]);
		}

		static void remove_negligible_values(powerset_btree<N, T>& definition, const T& neutral_value) {
			const std::vector<size_t>& indices = definition.elements_indices();
			for (size_t i = 0; i < indices.size(); ++i) {
				if(powerset_function<N, T>::is_equivalent_to_zero(definition.get_node(indices[i]).value - neutral_value)){
					definition.nullify(indices[i]);
				}
			}
		}

		void remove_negligible_values() {
			remove_negligible_values(this->definition, this->default_value);
		}

		void assign(const std::unordered_map<subset, T>& values) {
			for (const auto& labels_U_value : values) {
				this->assign(labels_U_value.first, labels_U_value.second);
			}
		}

		void assign(const std::unordered_map<std::vector<std::string>, T>& values) {
			for (const auto& labels_U_value : values) {
				this->assign(labels_U_value.first, labels_U_value.second);
			}
		}

		void assign(const subset& set, const T& value) {
			this->definition.update_or_insert(set, value);
		}

		void assign(const std::vector<std::string>& labels, const T& value) {
			this->assign(this->outcomes.get_subset(labels), value);
		}

		void assign_emptyset(const T& value) {
			this->assign(emptyset, value);
		}

		void assign_fullset(const T& value) {
			this->assign(fullset, value);
		}
	};

} // namespace efficient_DST

#endif // EFFICIENT_DST_MOBIUS_TRANSFORM_HPP
