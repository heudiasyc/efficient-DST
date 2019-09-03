#ifndef EFFICIENT_DST_BELIEF_HPP
#define EFFICIENT_DST_BELIEF_HPP

#include <implicability.hpp>

namespace efficient_DST{

	template <typename T = double>
	class belief : public implicability<T> {
	protected:

		static void display_invalid_belief_message(){
			std::cerr << "The Belief function is not defined for a mass function with non-null image at emptyset."
									<< "\nYou got the implicability function corresponding to the given information instead." << std::endl;
		}

	public:

		belief(const mass<T>& m) : implicability<T>(m)
		{
			const boost::dynamic_bitset<> emptyset(m.get_FOD().size());
			if(m.find(emptyset) != 0){
				display_invalid_belief_message();
			}
		}

		belief(const implicability<T>& b) : implicability<T>(b)
		{
			const boost::dynamic_bitset<> emptyset(b.get_FOD().size());
			if(b[emptyset] != 0){
				display_invalid_belief_message();
			}
		}

		belief(const powerset_btree<T>& focal_points_values) : implicability<T>(focal_points_values)
		{
			const boost::dynamic_bitset<> emptyset(focal_points_values.fod->size());
			if(focal_points_values[emptyset] != 0){
				display_invalid_belief_message();
			}
		}

		belief(const std::vector<T>& powerset_values, const FOD& fod) : implicability<T>(powerset_values, fod)
		{
			if(powerset_values[0] != 0){
				display_invalid_belief_message();
			}
		}


		template <class fusion_rule>
		belief<T> apply(const belief<T>& b2) const {
			const fusion_rule fusion;
			return fusion(*this, b2);
		}
	};

} // namespace efficient_DST

#endif // EFFICIENT_DST_BELIEF_HPP
