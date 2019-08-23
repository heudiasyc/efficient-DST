#ifndef EFFICIENT_DST_CONJUNCTIVE_WEIGHT_HPP
#define EFFICIENT_DST_CONJUNCTIVE_WEIGHT_HPP

#include <decomposition_weight.hpp>
#include <mobius_transform.hpp>
#include <mobius_aggregate.hpp>

namespace efficient_DST{

	template <typename T = double>
	class conjunctive_weight : public decomposition_weight<T> {
	public:

		conjunctive_weight(const conjunctive_weight<T>& w) : mobius_transform<T>(w.get_definition())
		{}

		conjunctive_weight(const powerset_btree<T>& focal_log_sets_values) : mobius_transform<T>(focal_log_sets_values)
		{}

		conjunctive_weight(const FOD& fod) : mobius_transform<T>(fod)
		{}

		conjunctive_weight(const mobius_aggregate<T>& q) : mobius_transform<T>(q.inversion(mobius_transformation_form_t::multiplicative))
		{
			if (q.order_relation != order_relation_t::superset) {
				std::cerr << "The given Möbius aggregate is not the commonality function and thus can only be the implicability one. "
						<< "\nDoing so, you got the disjunctive weight function instead." << std::endl;
			}
			invert_values(this->definition);
		}


		powerset_btree<T> inverted_definition() const {
			powerset_btree<T> inverted_definition(this->definition);
			invert_values(inverted_definition);
			return inverted_definition;
		}

		template <class fusion_rule>
		conjunctive_weight<T> apply(const fusion_rule fusion, const conjunctive_weight<T>& w2) const {
			return fusion(*this, w2);
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

#endif // EFFICIENT_DST_CONJUNCTIVE_WEIGHT_HPP
