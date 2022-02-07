#ifndef EFFICIENT_DST_PLAUSIBILITY_HPP
#define EFFICIENT_DST_PLAUSIBILITY_HPP

#include <implicability.hpp>

namespace efficient_DST{

	template <size_t N, typename T = float>
	class plausibility : public implicability<N, T> {
	public:
		using typename powerset_function<N, T>::subset;
		using zeta_transform<down_inclusion<N, T>, N, T>::operator[];

		plausibility(const mass<N, T>& m) : implicability<N, T>(m)
		{}

		plausibility(const disjunctive_weight<N, T>& v) : implicability<N, T>(v)
		{}

		plausibility(const implicability<N, T>& b) : implicability<N, T>(b)
		{}

		plausibility(const plausibility<N, T>& pl) : implicability<N, T>(pl)
		{}


		template <class fusion_rule>
		plausibility<N, T> fuse_with(const plausibility<N, T>& p2) const {
			const fusion_rule fusion;
			return fusion(*this, p2);
		}

		T* contour() {
			T contour[N] = {0};
			subset singleton = 1;

			for(size_t i = 0; i < N; ++i){
				contour[i] = (*this)[singleton];
				singleton <<= 1;
			}
			return contour;
		}

		T operator[](const subset& set) const {
			const subset& dual_set = ~set;
			set_N_value<N, T>* set_value = this->definition[dual_set];
			if(set_value){
				return 1-set_value->value;
			}
			return 1-this->find_non_focal_point_image(dual_set);
		}
	};

} // namespace efficient_DST

#endif // EFFICIENT_DST_PLAUSIBILITY_HPP

