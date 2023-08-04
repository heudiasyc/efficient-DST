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

#ifndef EFFICIENT_DST_PIGNISTIC_PROBABILITY_FUNCTION_HPP
#define EFFICIENT_DST_PIGNISTIC_PROBABILITY_FUNCTION_HPP

#include <mass_function.hpp>

namespace efficient_DST{

	template <size_t N, typename T = float>
	class pignistic_probability_function {
	protected:
		const mass_function<N, T> source_mass;

		T compute_aggregation(const std::bitset<N>& setA, const bool& singleton) const {
			T sum = 0;

			if (singleton){
				const std::vector<set_N_value<N, T>* >& supersets = this->source_mass.get_definition().supersets_of(setA);

				for (size_t i = 0; i < supersets.size(); ++i) {
					sum += supersets[i]->value / supersets[i]->cardinality;
				}
			}else{
				const std::map<size_t, std::vector<set_N_value<N, T>* >, std::less<size_t> >& cardinality_map = this->source_mass.get_definition().elements_by_ascending_cardinality();

				for (auto kv : cardinality_map) {
					if (kv.first != 0) {
						for (size_t i = 0; i < kv.second.size(); ++i) {
							const set_N_value<N, T>* B = kv.second[i];
							const std::bitset<N>& intersection = setA & B->set;

							sum += B->value * intersection.count() / kv.first;
						}
					}
				}
			}
			return sum;
		}

	public:

		pignistic_probability_function(const mass_function<N, T>& m) : source_mass(m)
		{}

		pignistic_probability_function(const pignistic_probability_function<N, T>& bet_p) : source_mass(bet_p.source_mass)
		{}

		std::vector<T> get_contour() {
			std::vector<T> contour;
			contour.reserve(N);
			std::bitset<N> singleton = 0;

			for(size_t i = 0; i < N; ++i){
				singleton[i] = true;
				contour.emplace_back(compute_aggregation(singleton, true));
				singleton[i] = false;
			}
			return contour;
		}

		T at_emptyset() const {
			return 0;
		}

		T at_fod() const {
			std::bitset<N> fod = 0;
			fod.set();
			return compute_aggregation(fod, false);
		}

		T operator[](const std::vector<std::string>& labels) const {
			return find(this->source_mass.get_FOD().to_set(labels));
		}

		T find(const std::bitset<N>& set) const {
			bool is_singleton = false;
			if (set.count() == 1){
				is_singleton = true;
			}
			return compute_aggregation(set, is_singleton);
		}
	};

} // namespace efficient_DST

#endif // EFFICIENT_DST_PIGNISTIC_PROBABILITY_FUNCTION_HPP
