#ifndef EFFICIENT_DST_RULE_CONJUNCTIVE_HPP
#define EFFICIENT_DST_RULE_CONJUNCTIVE_HPP

#include <commonality.hpp>
#include <conjunctive_weight.hpp>
#include <mass.hpp>

namespace efficient_DST{

	template <typename T = double>
	class rule_conjunctive {
		std::string to_string() const {
			return "Conjunctive rule";
		}

	public:

		mass<T> operator()(const mass<T>& m1, const mass<T>& m2) const {
			const powerset_btree<T>& m1_definition = m1.get_definition();
			const powerset_btree<T>& m2_definition = m2.get_definition();

			if (m1_definition.size() * m2_definition.size() < 0.5 * m1_definition.get_FOD_size() * pow(2, m1_definition.get_FOD_size())){
				const std::vector<set_N_value<T>* >& focal_sets_1 = m1_definition.elements();
				const std::vector<set_N_value<T>* >& focal_sets_2 = m2_definition.elements();
				powerset_btree<T> focal_sets_12(m1_definition.get_FOD(), m1_definition.get_block_size());

				for (size_t i1 = 0; i1 < focal_sets_1.size(); ++i1){
					for (size_t i2 = 0; i2 < focal_sets_2.size(); ++i2){
						const boost::dynamic_bitset<>& set = FOD::set_intersection(focal_sets_1[i1]->set, focal_sets_2[i2]->set);
						set_N_value<T>* node = focal_sets_12[set];
						if (node){
							node->value += focal_sets_1[i1]->value * focal_sets_2[i2]->value;
						}else{
							focal_sets_12.insert(set, focal_sets_1[i1]->value * focal_sets_2[i2]->value);
						}
					}
				}
				mass<T> m12(focal_sets_12);
				m12.remove_negligible_values();
				m12.normalize();
				return m12;
			}else{
				mobius_aggregate<T> q1(m1_definition, order_relation_t::superset, mobius_transformation_form_t::additive);
				mobius_aggregate<T> q2(m2_definition, order_relation_t::superset, mobius_transformation_form_t::additive);
				commonality<T> q12 = operator()(*(commonality<T>*) &q1, *(commonality<T>*) &q2);
				return mass<T>(q12);
			}
		}


		commonality<T> operator()(const commonality<T>& q1, const commonality<T>& q2) const {
			const powerset_btree<T>& w1_definition = q1.inversion(mobius_transformation_form_t::multiplicative);
			const powerset_btree<T>& w2_definition = q2.inversion(mobius_transformation_form_t::multiplicative);
			mobius_aggregate<T> q12(
				weight_fusion(w1_definition, w2_definition),
				order_relation_t::superset,
				mobius_transformation_form_t::multiplicative
			);
			return *(commonality<T>*) &q12;
		}


		conjunctive_weight<T> operator()(const conjunctive_weight<T>& w1, const conjunctive_weight<T>& w2) const {
			return conjunctive_weight<T>(weight_fusion(w1.get_definition(), w2.get_definition()));
		}


	protected:

		static T multiply(const T& val1, const T& val2){
			return val1 * val2;
		}

		powerset_btree<T> weight_fusion(const powerset_btree<T>& w1_definition, const powerset_btree<T>& w2_definition) const {
			powerset_btree<T> w12_definition(w1_definition.get_FOD(), w1_definition.get_block_size());
			w12_definition.fill_with_union_of_powersets(w1_definition, w2_definition, multiply, 1);
			decomposition_weight<T>::remove_negligible_values(w12_definition);
			decomposition_weight<T>::normalize(w12_definition);
			return w12_definition;
		}
	};

} // namespace efficient_DST

#endif // EFFICIENT_DST_RULE_CONJUNCTIVE_HPP
