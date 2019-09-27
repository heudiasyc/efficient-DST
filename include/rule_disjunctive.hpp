#ifndef EFFICIENT_DST_RULE_DISJUNCTIVE_HPP
#define EFFICIENT_DST_RULE_DISJUNCTIVE_HPP

#include <implicability.hpp>
#include <disjunctive_weight.hpp>
#include <mass.hpp>

namespace efficient_DST{

	template <typename T = double>
	class rule_disjunctive {
		std::string to_string() const {
			return "Disjunctive rule";
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
						const boost::dynamic_bitset<>& set = FOD::set_union(focal_sets_1[i1]->set, focal_sets_2[i2]->set);
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
				mobius_aggregate<T> b1(m1_definition, order_relation_t::subset, mobius_transformation_form_t::additive);
				mobius_aggregate<T> b2(m2_definition, order_relation_t::subset, mobius_transformation_form_t::additive);
				implicability<T> b12 = operator()(*(implicability<T>*) &b1, *(implicability<T>*) &b2);
				return mass<T>(b12);
			}
		}


		implicability<T> operator()(const implicability<T>& b1, const implicability<T>& b2) const {
			const powerset_btree<T>& v1_definition = b1.inversion(mobius_transformation_form_t::multiplicative);
			const powerset_btree<T>& v2_definition = b2.inversion(mobius_transformation_form_t::multiplicative);
			mobius_aggregate<T> b12(
				weight_fusion(v1_definition, v2_definition),
				order_relation_t::subset,
				mobius_transformation_form_t::multiplicative
			);
			return *(implicability<T>*) &b12;
		}


		disjunctive_weight<T> operator()(const disjunctive_weight<T>& v1, const disjunctive_weight<T>& v2) const {
			return disjunctive_weight<T>(weight_fusion(v1.get_definition(), v2.get_definition()));
		}


	protected:

		static T multiply(const T& val1, const T& val2){
			return val1 * val2;
		}

		powerset_btree<T> weight_fusion(const powerset_btree<T>& v1_definition, const powerset_btree<T>& v2_definition) const {
			powerset_btree<T> v12_definition(v1_definition.get_FOD(), v1_definition.get_block_size());
			v12_definition.fill_with_union_of_powersets(v1_definition, v2_definition, multiply, 1);
			decomposition_weight<T>::remove_negligible_values(v12_definition);
			decomposition_weight<T>::normalize(v12_definition);
			return v12_definition;
		}
	};

} // namespace efficient_DST

#endif // EFFICIENT_DST_RULE_DISJUNCTIVE_HPP
