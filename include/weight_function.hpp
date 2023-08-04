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
 
#ifndef EFFICIENT_DST_WEIGHT_FUNCTION_HPP
#define EFFICIENT_DST_WEIGHT_FUNCTION_HPP

#include <mobius_transform.hpp>
#include <zeta_transform.hpp>
#include <weight_vector.hpp>

namespace efficient_DST{

	template <size_t N, typename T = float>
	class weight_function : public mobius_transform<N, T>{
	public:
		using typename powerset_function<N, T>::subset;
		using powerset_function<N, T>::emptyset;
		using powerset_function<N, T>::fullset;

		weight_function(const weight_function<N, T>& w) : mobius_transform<N, T>(w.outcomes, w.definition, 1)
		{}

		weight_function(const weight_vector<N, T>& w) : mobius_transform<N, T>(w.outcomes, w.definition, 1)
		{}

		weight_function(
			const sample_space<N>& outcomes,
			const powerset_btree<N, T>& log_focal_sets
		) : mobius_transform<N, T>(outcomes, log_focal_sets, 1)
		{
			this->remove_negligible_values();
			this->normalize();
		}

		weight_function(
			const sample_space<N>& outcomes
		) : mobius_transform<N, T>(outcomes, 1)
		{}

		weight_function(const zeta_transform<up_inclusion<N, T>, N, T >& q) : weight_function<N, T>(q.get_sample_space(), q.inversion(operation_type_t::multiplication))
		{}

		weight_function(const zeta_transform<down_inclusion<N, T>, N, T >& b) : weight_function<N, T>(b.get_sample_space(), b.inversion(operation_type_t::multiplication))
		{}


		void regularize() {
			this->definition.nullify(this->definition[emptyset]);
			this->normalize();
		}

		void normalize() {
			normalize(this->definition);
		}

		static void normalize(powerset_btree<N, T>& definition) {
			T prod = 1;
			const std::vector<set_N_value<N, T>* >& elements = definition._elements();
			for (size_t i = 0; i < elements.size(); ++i) {
				prod *= elements[i]->value;
			}
			if(prod == 0){
				definition.print(true);
				std::cerr << "\nProduct of weight values equal to 0."
						<< "\nThis means that at least one of these weights is null, which is not allowed in a weight function." << std::endl;
				exit(1);
			}
			if(prod != 1){
				// normalize
				T factor = pow(prod, 1/elements.size());
				for (size_t i = 0; i < elements.size(); ++i) {
					elements[i]->value /= factor;
				}
			}
		}

		template<class inclusion>
		weight_function<N, T> natural_fusion_with(const weight_function<N, T>& w2) const {
			size_t max_def_size = std::max(this->definition.size(), w2.definition.size());
			weight_function<N, T> w12;
			subset core;
			if (w2.definition.size() == max_def_size){
				w12 = weight_function<N, T>(w2);
				natural_fusion<inclusion>(
					this->definition,
					w2.definition,
					w12,
					core
				);
			} else {
				w12 = weight_function<N, T>(*this);
				natural_fusion<inclusion>(
					w2.definition,
					this->definition,
					w12,
					core
				);
			}
			return w12;
		}

		weight_function<N, T> conjunctive_fusion_with(const weight_function<N, T>& w2) const {
			return natural_fusion_with<up_inclusion<N, T>>(w2);
		}

		weight_function<N, T> disjunctive_fusion_with(const weight_function<N, T>& w2) const {
			return natural_fusion_with<down_inclusion<N, T>>(w2);
		}

//	protected:

		template<class inclusion>
		static std::vector<subset> natural_fusion(
			const powerset_btree<N, T>& w1_definition,
			const powerset_btree<N, T>& w2_definition,
			powerset_btree<N, T>& w12_definition,
			subset& core
		){
			const std::vector<set_N_value<N, T> const * >& focal_points_w1 = w1_definition.elements();
			const std::vector<set_N_value<N, T> const * >& focal_points_w2 = w2_definition.elements();
			core = inclusion::absorbing_set_for_operation();
			subset core2 = inclusion::absorbing_set_for_operation();
			for (size_t i = 0; i < focal_points_w1.size(); ++i) {
				core = inclusion::set_dual_operation(core, focal_points_w1[i]->set);
			}
			for (size_t i = 0; i < focal_points_w2.size(); ++i) {
				core2 = inclusion::set_dual_operation(core2, focal_points_w2[i]->set);
			}
			core = inclusion::set_operation(core, core2);
			bool all_included = inclusion::has_any_element_related(w2_definition, core);
			std::unordered_set<subset> manifest;
			std::vector<subset> new_sets;
			T val;
			for (size_t i = 0; i < focal_points_w1.size(); ++i) {
				subset& set = focal_points_w1[i]->set;
				if (inclusion::set_dual_operation(core, set) == core){
					auto occurrence = manifest.find(set);
					if (occurrence == manifest.end()){
						manifest.emplace(set);
						val = focal_points_w1[i]->value;
						const size_t& assignment_index = w2_definition[set];
						if (assignment_index < w2_definition.number_of_nodes()){
							val *= w2_definition.get_node(assignment_index).value;
						} else if (all_included || inclusion::has_any_element_related(w2_definition, set)){
							new_sets.emplace_back(set);
						} else {
							continue;
						}
//						if (!powerset_function<N, T>::is_equivalent_to_zero(val - 1))
						w12_definition.update_or_insert(set, val);
					}
				}
			}
			all_included = inclusion::has_any_element_related(w1_definition, core);
			for (size_t i = 0; i < focal_points_w2.size(); ++i) {
				subset& set = focal_points_w2[i]->set;
				if (inclusion::set_dual_operation(core, set) == core){
					auto occurrence = manifest.find(set);
					if (occurrence == manifest.end()){
						manifest.emplace(set);
						if (!all_included && !inclusion::has_any_element_related(w1_definition, set)){
							w12_definition.nullify(w12_definition[set]);
						}
					}
				}
			}
			return new_sets;
		}
	};

} // namespace efficient_DST

#endif // EFFICIENT_DST_WEIGHT_FUNCTION_HPP
