#ifndef OW_BFT_PLAUSIBILITY_HPP
#define OW_BFT_PLAUSIBILITY_HPP

#include <mass_aggregate.hpp>
#include <implicability.hpp>

namespace ow_bft{

	template <typename T = double>
	class plausibility : public mass_aggregate<T> {
	protected:

		T compute_aggregation_at_emptyset() const {
			return 0;
		}

		T compute_aggregation_at_fod() const {
			return 1;
		}

		T compute_aggregation(const boost::dynamic_bitset<>& key) const {
			/*
			T sum = 0;
			std::vector<set_N_value<T>* > elements = this->focal_elements.elements();
			for (size_t i = 0; i < elements.size(); ++i) {
					boost::dynamic_bitset new_key = this->fod->to_set(elements[i]->fod_elements);
					if(!this->fod->are_disjoint(key, new_key))
						sum += elements[i]->value;
				}
			}
			return sum;
			*/

			/*
			 *  plausibility<FOD, T> pl(bel);
				std::reverse(pl.values().begin(), pl.values().end());
				T conflict_mass = 1 - pl[0];
				BOOST_FOREACH (T& v, pl.values()) {
					v = 1 - v - conflict_mass;
				}
			 */
			return 1 - implicability<T>::compute_aggregation(this->fod->set_negate(key));
		}

		T compute_aggregation(const std::vector<fod_element*>& fod_elements) const {
			return compute_aggregation(this->fod->to_set(fod_elements));
		}

	public:

		plausibility(const mass<T>& m) : mass_aggregate<T>(m)
		{
			this->set_values_for_special_elements();
		}

		plausibility(const mass_aggregate<T>& ma) : plausibility(ma.get_mass_equivalent())
		{}

		plausibility(const plausibility<T>& p) : mass_aggregate<T>(p.get_mass_equivalent())
		{
			this->special_elements.copy(p.get_special_elements());
		}

		plausibility(const FOD& fod) : mass_aggregate<T>(fod)
		{
			this->set_values_for_special_elements();
		}

		plausibility(const FOD& fod, const Special_case s_case) : mass_aggregate<T>(fod, s_case)
		{
			this->set_values_for_special_elements();
		}

		template <class fusion_rule>
		plausibility<T> apply(const fusion_rule fusion, const plausibility<T>& p2) const {
			return fusion(*this, p2);
		}

		mass<T> to_mass(const powerset_btree<T>& f_elements) const {

			std::cerr << "\nUnimplemented method to_mass of plausibility. Impossible to infer masses from plausibility on focal elements only.\n";

			// initialization
			mass<T> m(f_elements.get_FOD());

			return m;
		}
	};

} // namespace ow_bft

#endif // OW_BFT_PLAUSIBILITY_HPP

