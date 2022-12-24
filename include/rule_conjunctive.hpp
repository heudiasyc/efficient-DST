#ifndef EFFICIENT_DST_RULE_CONJUNCTIVE_HPP
#define EFFICIENT_DST_RULE_CONJUNCTIVE_HPP

#include <unordered_set>

#include <rule_classic_general.hpp>
#include <commonality_function.hpp>
#include <conjunctive_decomposition.hpp>
#include <mass_function.hpp>
#include <mobius_inversion.hpp>


namespace efficient_DST{

	template <class inclusion, size_t N, typename T = float>
	class rule_conjunctive : public rule_classic_general<up_inclusion<N, T>, N, T>{
	public:
		using typename rule_classic_general<up_inclusion<N, T>, N, T>::subset;

		std::string to_string() const {
			return "Conjunctive rule";
		}

		mass_function<N, T> operator()(const mass_function<N, T>& m1, const mass_function<N, T>& m2) const {
			const powerset_btree<N, T>& m1_definition = m1.get_definition();
			const powerset_btree<N, T>& m2_definition = m2.get_definition();
			size_t prod = m1_definition.size() * m2_definition.size();
			mass_function<N, T> m12(m1.get_sample_space());
			if (prod == 0)
				return m12;
			if (log2(prod / N) < N + 2){
				m1.template natural_fusion_with<up_inclusion<N, T>>(m2, m12);
				m12.remove_negligible_values();
				m12.normalize();
				return m12;
			}else{
				return mass_function<N, T>((*this)(commonality_function<N, T>(m1), commonality_function<N, T>(m2)));
			}
		}


		commonality_function<N, T> operator()(const commonality_function<N, T>& q1, const commonality_function<N, T>& q2) const {
			return commonality_function<N, T>((*this)(weight_function<N, T>(q1), weight_function<N, T>(q2)));
		}


		weight_function<N, T> operator()(const weight_function<N, T>& w1, const weight_function<N, T>& w2) const {
			const powerset_btree<N, T>& w1_definition = w1.get_definition();
			const powerset_btree<N, T>& w2_definition = w2.get_definition();
			size_t max_def_size = std::max(w1_definition.size(), w2_definition.size());
			powerset_btree<N, T> w12_definition;//(max_def_size);
			if (w2_definition.size() == max_def_size){
				w12_definition = w2_definition;
				fuse_weight_functions(
					w1_definition,
					w2_definition,
					w12_definition
				);
			} else {
				w12_definition = w1_definition;
				fuse_weight_functions(
					w2_definition,
					w1_definition,
					w12_definition
				);
			}
			weight_function<N, T> w12(w1.get_sample_space(), w12_definition);
//			powerset_btree<N, T>& w12_definition = w12.definition_();
//			set_N_value<N, T>* node;
//			T val;
//			for (const auto& set : manifest) {
////				w12.assign(set, w1[set] * w2[set]);
//				val = 1;
//				node = w1_definition[set];
//				if (node){
//					val = node->value;
//				}
//				node = w2_definition[set];
//				if (node){
//					val *= node->value;
//				}
//				if (!powerset_function<N, T>::is_equivalent_to_zero(val - 1))
//					w12_definition.insert(set, val);
//			}
			return w12;
		}

		static void fuse_weight_functions(
			const powerset_btree<N, T>& w1_definition,
			const powerset_btree<N, T>& w2_definition,
			powerset_btree<N, T>& w12_definition
		){
			const std::vector<set_N_value<N, T> const * >& focal_points_w1 = w1_definition.elements();
			const std::vector<set_N_value<N, T> const * >& focal_points_w2 = w2_definition.elements();
			subset core1 = inclusion::absorbing_set_for_operation();
			subset core2 = inclusion::absorbing_set_for_operation();
			for (size_t i = 0; i < focal_points_w1.size(); ++i) {
				core1 = inclusion::set_dual_operation(core1, focal_points_w1[i]->set);
			}
			for (size_t i = 0; i < focal_points_w2.size(); ++i) {
				core2 = inclusion::set_dual_operation(core2, focal_points_w2[i]->set);
			}
			core1 = inclusion::set_operation(core1, core2);
			bool all_included = inclusion::has_any_element_related(w2_definition, core1);
			std::unordered_set<subset> manifest;
			T val;
			for (size_t i = 0; i < focal_points_w1.size(); ++i) {
				subset& set = focal_points_w1[i]->set;
				if (inclusion::set_dual_operation(core1, set) == core1){
					auto occurrence = manifest.find(set);
					if (occurrence == manifest.end()){
						manifest.emplace(set);
						val = focal_points_w1[i]->value;
						const size_t& assignment_index = w2_definition[set];
						if (assignment_index < w2_definition.number_of_nodes()){
							val *= w2_definition.get_node(assignment_index).value;
						} else if (!all_included && !inclusion::has_any_element_related(w2_definition, set)){
							continue;
						}
//						if (!powerset_function<N, T>::is_equivalent_to_zero(val - 1))
						w12_definition.update_or_insert(set, val);
					}
				}
			}
			all_included = inclusion::has_any_element_related(w1_definition, core1);
			for (size_t i = 0; i < focal_points_w2.size(); ++i) {
				subset& set = focal_points_w2[i]->set;
				if (inclusion::set_dual_operation(core1, set) == core1){
					auto occurrence = manifest.find(set);
					if (occurrence == manifest.end()){
						manifest.emplace(set);
						if (all_included || inclusion::has_any_element_related(w1_definition, set)){
							w12_definition.update_or_insert(set, focal_points_w2[i]->value);
						} else {
							continue;
						}
					}
				}
			}
		}
	};

} // namespace efficient_DST

#endif // EFFICIENT_DST_RULE_CONJUNCTIVE_HPP
