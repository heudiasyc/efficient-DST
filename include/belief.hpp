#ifndef EFFICIENT_DST_BELIEF_HPP
#define EFFICIENT_DST_BELIEF_HPP

#include <implicability.hpp>

namespace efficient_DST{

	template <typename T, size_t N>
	class belief : public implicability<T, N> {
	protected:

		static void display_invalid_belief_message(){
			std::cerr << "The Belief function is not defined for a mass function with non-null image at emptyset."
									<< "\nYou got the implicability function corresponding to the given information instead." << std::endl;
		}

	public:

		belief(const mass<T, N>& m) : implicability<T, N>(m)
		{
			if(m.at_emptyset() != 0){
				display_invalid_belief_message();
			}
		}

		belief(const implicability<T, N>& b) : implicability<T, N>(b)
		{
			if(b.at_emptyset() != 0){
				display_invalid_belief_message();
			}
		}

		belief(const powerset_btree<T, N>& focal_points_values) : implicability<T, N>(focal_points_values)
		{
			if(focal_points_values.sub_fod_of_size(0)){
				display_invalid_belief_message();
			}
		}

		belief(const std::vector<T>& powerset_values, const FOD<N>& fod) : implicability<T, N>(powerset_values, fod)
		{
			if(powerset_values[0] != 0){
				display_invalid_belief_message();
			}
		}


		template <class fusion_rule>
		belief<T, N> apply(const belief<T, N>& b2) const {
			const fusion_rule fusion;
			return fusion(*this, b2);
		}
	};

} // namespace efficient_DST

#endif // EFFICIENT_DST_BELIEF_HPP
