#ifndef EFFICIENT_DST_DISJUNCTIVE_WEIGHT_HPP
#define EFFICIENT_DST_DISJUNCTIVE_WEIGHT_HPP

#include <decomposition_weight.hpp>
#include <mobius_transform.hpp>
#include <mobius_aggregate.hpp>

namespace efficient_DST{

	template <typename T = double>
	class disjunctive_weight : public decomposition_weight<T> {
	public:

		disjunctive_weight(const disjunctive_weight<T>& v) : mobius_transform<T>(v.get_definition())
		{}

		disjunctive_weight(const powerset_btree<T>& focal_log_sets_values) : mobius_transform<T>(focal_log_sets_values)
		{}

		disjunctive_weight(const FOD& fod) : mobius_transform<T>(fod)
		{}

		disjunctive_weight(const mobius_aggregate<T>& b) : mobius_transform<T>(b.inversion(mobius_transformation_form_t::multiplicative))
		{
			if (b.order_relation != order_relation_t::subset) {
				std::cerr << "The given Möbius aggregate is not the implicability function and thus can only be the commonality one. "
						<< "\nDoing so, you got the conjunctive weight function instead." << std::endl;
			}
			invert_values(this->definition);
		}


		powerset_btree<T> inverted_definition() const {
			powerset_btree<T> inverted_definition(this->definition);
			invert_values(inverted_definition);
			return inverted_definition;
		}

		template <class fusion_rule>
		disjunctive_weight<T> apply(const fusion_rule fusion, const disjunctive_weight<T>& v2) const {
			return fusion(*this, v2);
		}

	protected:

		void invert_values(powerset_btree<T>& values) const {
			std::vector<set_N_value<T>* > elements = values.elements();
			for (size_t i = 0; i < elements.size(); ++i){
				elements[i]->value = 1 / elements[i]->value;
			}
		}
	};

} // namespace efficient_DST

#endif // EFFICIENT_DST_DISJUNCTIVE_WEIGHT_HPP
