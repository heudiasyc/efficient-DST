#ifndef EFFICIENT_DST_CONJUNCTIVE_WEIGHT_HPP
#define EFFICIENT_DST_CONJUNCTIVE_WEIGHT_HPP

#include <decomposition_weight.hpp>
#include <mobius_transform.hpp>
#include <mobius_aggregate.hpp>

namespace efficient_DST{

	template <typename T = double>
	class conjunctive_weight : public decomposition_weight<T> {
	public:

		conjunctive_weight(const conjunctive_weight<T>& w) : decomposition_weight<T>(w.get_definition())
		{}

		conjunctive_weight(const powerset_btree<T>& focal_log_sets_values) : decomposition_weight<T>(focal_log_sets_values)
		{}

		conjunctive_weight(const FOD& fod) : decomposition_weight<T>(fod)
		{}

		conjunctive_weight(const mobius_aggregate<T>& q) : decomposition_weight<T>(q)
		{
			if (q.order_relation != order_relation_t::superset) {
				std::cerr << "The given Möbius aggregate is not the commonality function and thus can only be the implicability one. "
						<< "\nDoing so, you got the disjunctive weight function instead." << std::endl;
			}
			invert_values(this->definition);
			boost::dynamic_bitset<> fod(this->definition.fod->size());
			fod.set();
			// fod is not supposed to be defined in the conjunctive weight function
			this->definition.nullify(this->definition[fod]);
		}


		/*
		 * The conjunctive weight function is the inverse of the multiplicative Möbius transform of the commonality function.
		 * So, invert its values before computing it.
		 */
		powerset_btree<T> inverted_definition() const {
			powerset_btree<T> inverted_definition(this->definition);
			invert_values(inverted_definition);
			// The following part ensures that w(fod) is defined (as it is not when directly building w).
			// Indeed, w(fod) is required for the computation of the commonality function.
			boost::dynamic_bitset<> fod(inverted_definition.fod->size());
			fod.set();
			// inv_w_fod represents w(fod)^{-1}, which is also equal to q(fod) and m(fod).
			T inv_w_fod = 1;
			const std::vector<set_N_value<T>* >& subsets = inverted_definition.strict_subsets_of(fod);
			for (size_t i = 0; i < subsets.size(); ++i){
				// divide by the images of w^{-1} to obtain the product on images of w which corresponds to q(fod)
				inv_w_fod /= subsets[i]->value;
			}
			inverted_definition.insert(fod, inv_w_fod);
			return inverted_definition;
		}


		template <class fusion_rule>
		conjunctive_weight<T> apply(const fusion_rule fusion, const conjunctive_weight<T>& w2) const {
			return fusion(*this, w2);
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

#endif // EFFICIENT_DST_CONJUNCTIVE_WEIGHT_HPP
