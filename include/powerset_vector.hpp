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

#ifndef EFFICIENT_DST_POWERSET_VECTOR_HPP
#define EFFICIENT_DST_POWERSET_VECTOR_HPP

#include <vector>
#include <iostream>
#include <iomanip>

#include <sample_space.hpp>
#include <powerset_btree.hpp>
#include <powerset_function.hpp>


namespace efficient_DST{

	template <size_t N, typename T = float>
	class powerset_vector {
	protected:
		sample_space<N> outcomes;
		std::vector<T> definition;
		const T default_value;

	public:
		powerset_vector (const sample_space<N>& outcomes, const std::vector<T>& definition, const T& default_value) :
			outcomes(outcomes),
			definition(definition),
			default_value(default_value)
		{}

		powerset_vector (const sample_space<N>& outcomes, const T& default_value) :
			outcomes(outcomes),
//			definition(),
			default_value(default_value)
		{}

		virtual ~powerset_vector()
		{}

		// =============================================================================

		const std::vector<T>& get_definition() const {
			return this->definition;
		}

		const sample_space<N>& get_sample_space() const {
			return this->outcomes;
		}

		const T& get_default_value() const {
			return this->default_value;
		}

		T at_emptyset() const {
			return this->definition[0];
		}

		T at_fullset() const {
			return this->definition[(1 << N)-1];
		}

		T operator[](const std::vector<std::string>& labels) const {
			return (*this)[this->outcomes.get_subset_index(labels)];
		}

		T operator[](const size_t& index) const {
			return this->definition[index];
		}

		std::ostream& print(const bool& including_null = false) const {
			std::cout << std::endl;
			for (size_t i = 0; i < this->definition.size(); ++i) {
				if (including_null || !powerset_function<N, T>::is_equivalent_to_zero(this->definition[i] - this->default_value)){
					std::cout << set_N_value<N, T>::to_string(this->definition[i]) + "\t <- " + outcomes.to_string(i) << std::endl;
				}
			}
			return std::cout;
		}
	};
}		// namespace efficient_DST

#endif // EFFICIENT_DST_POWERSET_VECTOR_HPP
