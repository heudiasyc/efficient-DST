#ifndef EFFICIENT_DST_PIGNISTIC_PROBABILITY_HPP
#define EFFICIENT_DST_PIGNISTIC_PROBABILITY_HPP

#include <mobius_aggregate.hpp>

namespace efficient_DST{

	template <typename T = double>
	class pignistic_probability : public mobius_aggregate<T> {
	protected:

		T compute_aggregation_at_emptyset() const {
			return 0;
		}

		T compute_aggregation_at_fod() const {
			return 1;
		}

		T compute_aggregation(const boost::dynamic_bitset<>& setA) const {
			T sum = 0;
			const std::vector<set_N_value<T>* >& f_elements = this->mass_equivalent.focal_elements.elements();

			for (size_t i = 0; i < f_elements.size(); ++i) {

				const boost::dynamic_bitset<>& setB = f_elements[i]->set;
				const boost::dynamic_bitset<>& intersection = this->fod.set_intersection(setA, setB);

				sum += f_elements[i]->value * intersection.count() / setB.count();
			}
			return sum;
		}

	public:

		pignistic_probability(const mass<T>& m) : mass_aggregate<T>(m)
		{
			this->set_values_for_special_elements();
		}

		pignistic_probability(const mobius_aggregate<T>& ma) : pignistic_probability(ma.get_mass_equivalent())
		{}

		pignistic_probability(const pignistic_probability<T>& bet_p) : mass_aggregate<T>(bet_p.get_mass_equivalent())
		{
			this->special_elements.copy(bet_p.get_special_elements());
		}

		pignistic_probability(const FOD& fod) : mass_aggregate<T>(fod)
		{
			this->set_values_for_special_elements();
		}

		pignistic_probability(const FOD& fod, const special_case_t s_case) : mass_aggregate<T>(fod, s_case)
		{
			this->set_values_for_special_elements();
		}

		template <class fusion_rule>
		pignistic_probability<T> apply(const fusion_rule fusion, const pignistic_probability<T>& p2) const {
			return fusion(*this, p2);
		}

		mass<T> to_mass(const powerset_btree<T>& f_elements) const {

			std::cerr << "\nUnimplemented method to_mass of pignistic_probability. Impossible to infer masses from pignistic probability on focal elements only.\n";

			// initialization
			mass<T> m(f_elements.get_FOD());

			return m;
		}
	};

} // namespace efficient_DST

#endif // EFFICIENT_DST_PIGNISTIC_PROBABILITY_HPP
