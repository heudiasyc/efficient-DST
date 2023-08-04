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
 
#ifndef EFFICIENT_DST_DECOMPOSITION_HPP
#define EFFICIENT_DST_DECOMPOSITION_HPP

#include <weight_function.hpp>
#include <zeta_transform.hpp>

namespace efficient_DST{

	template <class inclusion, size_t N, typename T = float>
	class decomposition : public weight_function<N, T> {
	public:
		using typename powerset_function<N, T>::subset;
		using powerset_function<N, T>::emptyset;
		using powerset_function<N, T>::fullset;

		subset normalizing_set;
		T normalizing_value = 1;
		const bool adaptive_uncertainty;


		decomposition(
			const weight_function<N, T>& w,
			const bool& adaptive_uncertainty = true
		) :
			weight_function<N, T>(w.get_sample_space(), w.get_definition()),
			adaptive_uncertainty(adaptive_uncertainty)
		{
			this->remove_negligible_values();
			compute_normalizing_assignment();
		}

		decomposition(
			const sample_space<N>& outcomes,
			const bool& adaptive_uncertainty = true
		) :
			weight_function<N, T>(outcomes),
			adaptive_uncertainty(adaptive_uncertainty)
		{
			compute_normalizing_set();
		}

		decomposition(
			const zeta_transform<inclusion, N, T>& q,
			const bool& adaptive_uncertainty = true
		) :
			weight_function<N, T>(q),
			adaptive_uncertainty(adaptive_uncertainty)
		{
			this->remove_negligible_values();
			compute_normalizing_assignment();
		}


		bool has_adaptive_uncertainty() const {
			return this->adaptive_uncertainty;
		}

		T at_emptyset() const {
			return (*this)[emptyset];
		}

		T at_fullset() const {
			return (*this)[fullset];
		}

		T operator[](const std::vector<std::string>& labels) const {
			return (*this)[this->outcomes.get_subset(labels)];
		}

		T operator[](const subset& set) const {
			return 1 - 1/weight_function<N, T>::operator[](set);
		}

		std::ostream& print(const bool& including_null = false) {
			const std::vector<set_N_value<N, T> const * >& values = this->definition.elements(including_null);
			std::cout << std::endl;
			for (size_t i = 0; i < values.size(); ++i) {
				if(values[i]->set != this->normalizing_set){
					if (values[i]->is_null)
						std::cout << "NULL\t <- " + this->outcomes.to_string(values[i]->set) << std::endl;
					else
						std::cout << set_N_value<N, T>::to_string(1-1/values[i]->value) + "\t <- " + this->outcomes.to_string(values[i]->set) << "\t,\t "
						<< set_N_value<N, T>::to_string(1/values[i]->value) + "\t <- " + this->outcomes.to_string(this->normalizing_set) << std::endl;
				}
			}
//			this->definition.print_layout();
			return std::cout;
		}

		void assign(const subset& set, const T& mass) {
			T weight = 1-mass;
			if(weight == 0){
				std::cout << "Cannot have any categorical simple belief assignment in an uncertainty decomposition. Ignoring command.\n";
				return;
			}
			this->definition.update_or_insert(set, 1/weight);
			this->normalizing_value *= weight;
			if(adaptive_uncertainty){
				this->normalizing_set = inclusion::set_dual_operation(this->normalizing_set, set);
			}
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

		void assign(const std::vector<std::string>& labels, const T& value) {
			this->assign(this->outcomes.get_subset(labels), value);
		}

		void assign_emptyset(const T& value) {
			this->assign(emptyset, value);
		}

		void assign_fullset(const T& value) {
			this->assign(fullset, value);
		}

		void nullify(const std::vector<std::string>& labels) {
			this->nullify(this->outcomes.get_subset(labels));
		}

		void nullify(const subset& set) {
//			if (set != fullset){
			size_t index = this->definition[set];
			if(index < this->definition.number_of_nodes()){
				this->definition.nullify(index);
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
			if(adaptive_uncertainty){
				this->normalizing_set = inclusion::absorbing_set_for_operation();
				const std::vector<set_N_value<N, T> const * >& focal_log_elements = this->definition.elements();
				for (size_t i = 0; i < focal_log_elements.size(); ++i){
					this->normalizing_set = inclusion::set_dual_operation(this->normalizing_set, focal_log_elements[i]->set);
				}
			}else{
				this->normalizing_set = inclusion::absorbing_set_for_dual_operation();
			}
		}

		void compute_normalizing_value(){
			this->normalizing_value = 1;
			const std::vector<set_N_value<N, T>* >& focal_log_elements = this->definition.elements();
			for (size_t i = 0; i < focal_log_elements.size(); ++i){
				if (focal_log_elements[i]->set != this->normalizing_set)
					this->normalizing_value /= focal_log_elements[i]->value;
			}
		}

		void compute_normalizing_assignment(){
			const std::vector<set_N_value<N, T> const * >& focal_log_elements = this->definition.elements();
			this->normalizing_value = 1;
			if(adaptive_uncertainty){
				this->normalizing_set = inclusion::absorbing_set_for_operation();
				for (size_t i = 0; i < focal_log_elements.size(); ++i){
					this->normalizing_set = inclusion::set_dual_operation(this->normalizing_set, focal_log_elements[i]->set);
				}
			}else{
				this->normalizing_set = inclusion::absorbing_set_for_dual_operation();
			}
			for (size_t i = 0; i < focal_log_elements.size(); ++i){
				if (focal_log_elements[i]->set != this->normalizing_set)
					this->normalizing_value /= focal_log_elements[i]->value;
			}
		}
	};

} // namespace efficient_DST

#endif // EFFICIENT_DST_DECOMPOSITION_HPP
