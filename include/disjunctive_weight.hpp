#ifndef EFFICIENT_DST_DISJUNCTIVE_WEIGHT_HPP
#define EFFICIENT_DST_DISJUNCTIVE_WEIGHT_HPP

#include <decomposition_weight.hpp>
#include <mobius_transform.hpp>
#include <mobius_aggregate.hpp>

namespace efficient_DST{

	template <typename T = double>
	class disjunctive_weight : public decomposition_weight<T> {
	public:

		disjunctive_weight(const disjunctive_weight<T>& v) : decomposition_weight<T>(v.get_definition())
		{}

		disjunctive_weight(const powerset_btree<T>& focal_log_sets_values) : decomposition_weight<T>(focal_log_sets_values)
		{}

		disjunctive_weight(const FOD& fod) : decomposition_weight<T>(fod)
		{}

		disjunctive_weight(const mobius_aggregate<T>& b) : decomposition_weight<T>(b)
		{
			if (b.order_relation != order_relation_t::subset) {
				std::cerr << "The given Möbius aggregate is not the implicability function and thus can only be the commonality one. "
						<< "\nDoing so, you got the conjunctive weight function instead." << std::endl;
			}
			invert_values(this->definition);
			boost::dynamic_bitset<> emptyset(this->definition.fod->size());
			// fod is not supposed to be defined in the conjunctive weight function
			this->definition.nullify(this->definition[emptyset]);
		}


		/*
		 * The disjunctive weight function is the inverse of the multiplicative Möbius transform of the implicability function.
		 * So, invert its values before computing it.
		 */
		powerset_btree<T> inverted_definition() const {
			powerset_btree<T> inverted_definition(this->definition);
			invert_values(inverted_definition);
			// The following part ensures that v(emptyset) is defined (as it is not when directly building v).
			// Indeed, v(emptyset) is required for the computation of the implicability function.
			boost::dynamic_bitset<> emptyset(inverted_definition.fod->size());
			// inv_v_fod represents v(emptyset)^{-1}, which is also equal to b(fod) and m(fod).
			T inv_v_fod = 1;
			const std::vector<set_N_value<T>* >& supersets = inverted_definition.strict_supersets_of(emptyset);
			for (size_t i = 0; i < supersets.size(); ++i){
				// divide by the images of v^{-1} to obtain the product on images of v which corresponds to b(fod)
				inv_v_fod /= supersets[i]->value;
			}
			inverted_definition.insert(emptyset, inv_v_fod);
			return inverted_definition;
		}


		template <class fusion_rule>
		disjunctive_weight<T> apply(const disjunctive_weight<T>& v2) const {
			const fusion_rule fusion;
			return fusion(*this, v2);
		}

	protected:

		static void invert_values(powerset_btree<T>& values) {
			std::vector<set_N_value<T>* > elements = values.elements();
			for (size_t i = 0; i < elements.size(); ++i){
				elements[i]->value = 1 / elements[i]->value;
			}
		}
	};

} // namespace efficient_DST

#endif // EFFICIENT_DST_DISJUNCTIVE_WEIGHT_HPP
