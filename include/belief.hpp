#ifndef EFFICIENT_DST_BELIEF_HPP
#define EFFICIENT_DST_BELIEF_HPP

#include <implicability.hpp>

namespace efficient_DST{

	template <size_t N, typename T = float>
	class belief : public implicability<N, T> {
	protected:

		static void display_invalid_belief_message(){
			std::cerr << "The Belief function is not defined for a mass function with non-null image at emptyset."
									<< "\nYou got the implicability function corresponding to the given information instead." << std::endl;
		}

	public:

		belief(const mass<N, T>& m) : implicability<N, T>(m)
		{
			if(m.at_emptyset() != 0){
				display_invalid_belief_message();
			}
		}

		belief(const implicability<N, T>& b) : implicability<N, T>(b)
		{
			if(b.at_emptyset() != 0){
				display_invalid_belief_message();
			}
		}

//		belief(
//			const powerset_btree<N, T>& focal_points_values,
//			const scheme_type_t& scheme_type,
//			const std::vector<std::bitset<N> >& iota_sequence,
//			const T& neutral_value
//		) : zeta_transform<T, N, down_inclusion<N, T> >(
//				focal_points_values,
//				scheme_type,
//				iota_sequence,
//				neutral_value
//			)
//		{
//			if(focal_points_values.empty_set()){
//				display_invalid_belief_message();
//			}
//		}


		template <class fusion_rule>
		belief<N, T> fuse_with(const belief<N, T>& b2) const {
			const fusion_rule fusion;
			return fusion(*this, b2);
		}
	};

} // namespace efficient_DST

#endif // EFFICIENT_DST_BELIEF_HPP
