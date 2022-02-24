#ifndef EFFICIENT_DST_RULE_CONJUNCTIVE_HPP
#define EFFICIENT_DST_RULE_CONJUNCTIVE_HPP

#include <unordered_set>

#include <rule_classic_general.hpp>
#include <commonality_function.hpp>
#include <conjunctive_decomposition.hpp>
#include <mass_function.hpp>

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
			if (prod == 0)
				return mass_function<N, T>(m1.get_sample_space());
			if (log2(prod / N) < N + 2){
				return rule_classic_general<up_inclusion<N, T>, N, T>::operator()(m1, m2);
			}else{
				return mass_function<N, T>((*this)(commonality_function<N, T>(m1), commonality_function<N, T>(m2)));
			}
		}


		commonality_function<N, T> operator()(const commonality_function<N, T>& q1, const commonality_function<N, T>& q2) const {
			return commonality_function<N, T>((*this)(weight_function<N, T>(q1), weight_function<N, T>(q2)));
		}


		weight_function<N, T> operator()(const weight_function<N, T>& w1, const weight_function<N, T>& w2) const {
			std::unordered_set<subset > manifest;
			const powerset_btree<N, T>& w1_definition = w1.get_definition();
			const powerset_btree<N, T>& w2_definition = w2.get_definition();
			const std::vector<set_N_value<N, T>* >& focal_points_w1 = w1_definition.elements();
			const std::vector<set_N_value<N, T>* >& focal_points_w2 = w2_definition.elements();
			for (size_t i = 0; i < focal_points_w1.size(); ++i) {
				manifest.emplace(focal_points_w1[i]->set);
			}
			for (size_t i = 0; i < focal_points_w2.size(); ++i) {
				manifest.emplace(focal_points_w2[i]->set);
			}
			weight_function<N, T> w12(w1.get_sample_space());
			powerset_btree<N, T>& w12_definition = w12.definition_();
			set_N_value<N, T>* node;
			T val;
			for (const auto& set : manifest) {
//				w12.assign(set, w1[set] * w2[set]);
				val = 1;
				node = w1_definition[set];
				if (node){
					val = node->value;
				}
				node = w2_definition[set];
				if (node){
					val *= node->value;
				}
				if (!powerset_function<N, T>::is_equivalent_to_zero(val - 1))
					w12_definition.insert(set, val);
			}
			return w12;
		}
	};

} // namespace efficient_DST

#endif // EFFICIENT_DST_RULE_CONJUNCTIVE_HPP
