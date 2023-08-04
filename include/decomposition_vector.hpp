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

#ifndef EFFICIENT_DST_DECOMPOSITION_VECTOR_HPP
#define EFFICIENT_DST_DECOMPOSITION_VECTOR_HPP

#include <weight_vector.hpp>
#include <zeta_transform_vector.hpp>

namespace efficient_DST{

	template <class inclusion, size_t N, typename T = float>
	class decomposition_vector : public weight_vector<N, T> {
	public:
		size_t normalizing_set;
		T normalizing_value = 1;
		const bool adaptive_uncertainty;


		decomposition_vector(
			const weight_vector<N, T>& w,
			const bool& adaptive_uncertainty = true
		) :
			weight_vector<N, T>(w.get_sample_space(), w.get_definition()),
			adaptive_uncertainty(adaptive_uncertainty)
		{
			this->remove_negligible_values();
			compute_normalizing_assignment();
		}

		decomposition_vector(
			const sample_space<N>& outcomes,
			const bool& adaptive_uncertainty = true
		) :
			weight_vector<N, T>(outcomes),
			adaptive_uncertainty(adaptive_uncertainty)
		{
			compute_normalizing_set();
		}

		decomposition_vector(
			const zeta_transform_vector<inclusion, N, T>& q,
			const size_t& core,
			const size_t& zero,
			const bool& adaptive_uncertainty = true
		) :
			weight_vector<N, T>(
				q.get_sample_space(),
				q.expansion_from_core(
					q.inversion(operation_type_t::multiplication, adaptive_uncertainty),
					1,
					core,
					zero
				)
			),
			adaptive_uncertainty(adaptive_uncertainty)
		{
			compute_normalizing_assignment();
		}


		bool has_adaptive_uncertainty() const {
			return this->adaptive_uncertainty;
		}

		T at_emptyset() const {
			return (*this)[0];
		}

		T at_fullset() const {
			return (*this)[(1 << N)-1];
		}

		T operator[](const std::vector<std::string>& labels) const {
			return (*this)[this->outcomes.get_subset(labels)];
		}

		T operator[](const size_t& set) const {
			return 1 - 1/weight_vector<N, T>::operator[](set);
		}

		std::ostream& print() const {
//			sample_space<N> reduced_outcomes(this->outcomes.get_labels(this->normalizing_set));
			std::cout << std::endl;
			for (size_t i = 0; i < this->definition.size(); ++i) {
				if (inclusion::set_operation(i, this->normalizing_set) == i && i != this->normalizing_set && this->definition[i] != 1){
					std::cout << set_N_value<N, T>::to_string((*this)[i]) + "\t <- " + this->outcomes.to_string(i) << "\t,\t "
					<< set_N_value<N, T>::to_string(1/this->definition[i]) + "\t <- " + this->outcomes.to_string(this->normalizing_set) << std::endl;
				}
			}
//			this->definition.print_layout();
			return std::cout;
		}

		void assign(const size_t& set, const T& mass) {
			T weight = 1-mass;
			if(weight == 0){
				std::cout << "Cannot have any categorical simple belief assignment in an uncertainty decomposition. Ignoring command.\n";
				return;
			}
			this->definition[set] = 1/weight;
			this->normalizing_value *= weight;
			if(this->adaptive_uncertainty){
				this->normalizing_set = inclusion::set_dual_operation(this->normalizing_set, set);
			}
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

		void assign(const std::vector<std::string>& labels, const T& value) {
			this->assign(this->outcomes.get_subset_index(labels), value);
		}

		void assign_emptyset(const T& value) {
			this->assign(0, value);
		}

		void assign_fullset(const T& value) {
			this->assign((1 << N)-1, value);
		}

		void nullify(const std::vector<std::string>& labels) {
			this->nullify(this->outcomes.get_subset_index(labels));
		}

		void nullify(const size_t& set) {
//			if (set != fullset){
			if (this->definition[set] != 1){
				this->definition[set] = 1;
				this->compute_normalizing_set_assignment();
			}
//			}else{
//				std::cout << "Cannot directly assign the normalizing value of a decomposition. Ignoring command.\n";
//			}
		}

//		set_N_value<N, T> compute_encompassing_set_assignment(){
//			compute_encompassing_set_assignment(this->definition);
//		}

		void compute_normalizing_set(){
			if(this->adaptive_uncertainty){
				this->normalizing_set = inclusion::set_operation(0, (1 << N)-1);
				for (size_t i = 0; i < (1 << N); ++i){
					if (!powerset_function<N, T>::is_equivalent_to_zero(this->definition[i] - 1))
						this->normalizing_set = inclusion::set_dual_operation(this->normalizing_set, i);
				}
			}else{
				this->normalizing_set = inclusion::set_dual_operation(0, (1 << N)-1);
			}
		}

		void compute_normalizing_value(){
			this->normalizing_value = 1;
			for (size_t i = 0; i < this->definition.size(); ++i){
				if (i != this->normalizing_set)
					this->normalizing_value /= this->definition[i];
			}
		}

		void compute_normalizing_assignment(){
			compute_normalizing_set();
			compute_normalizing_value();
		}

		std::bitset<N> get_core() const {
			return inclusion::FMT_target(~std::bitset<N>(this->normalizing_set));
		}

//		const sample_space<N> get_sample_space() const {
//			return sample_space<N>(this->outcomes.get_labels(get_core()));
//		}

		std::vector<T> get_definition() const {
			if (adaptive_uncertainty){
				size_t n = get_core().count();
				std::vector<T> def(1 << n, this->default_value);
				size_t j = 0;
				for (size_t i = 0; i < (1 << N); ++i){
					if (inclusion::set_operation(i, this->normalizing_set) == i){
						def[j] = this->definition[i];
						++j;
					}
				}
				return def;
			}else{
				return this->definition;
			}
		}
	};

} // namespace efficient_DST

#endif // EFFICIENT_DST_DECOMPOSITION_VECTOR_HPP
