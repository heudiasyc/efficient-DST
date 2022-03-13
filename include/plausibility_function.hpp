#ifndef EFFICIENT_DST_PLAUSIBILITY_FUNCTION_HPP
#define EFFICIENT_DST_PLAUSIBILITY_FUNCTION_HPP

#include <implicability_function.hpp>

namespace efficient_DST{

	template <size_t N, typename T = float>
	class plausibility_function : public implicability_function<N, T> {
	public:
		using typename powerset_function<N, T>::subset;
		using powerset_function<N, T>::emptyset;
		using powerset_function<N, T>::fullset;

		plausibility_function(const mass_function<N, T>& m) : implicability_function<N, T>(m)
		{}

		plausibility_function(const weight_function<N, T>& v) : implicability_function<N, T>(v)
		{}

		plausibility_function(const disjunctive_decomposition<N, T>& v) : implicability_function<N, T>(v)
		{}

		plausibility_function(const mass_function<N, T>& m, scheme_type_t scheme_type) : implicability_function<N, T>(m, scheme_type)
		{}

		plausibility_function(const weight_function<N, T>& v, scheme_type_t scheme_type) : implicability_function<N, T>(v, scheme_type)
		{}

		plausibility_function(const disjunctive_decomposition<N, T>& v, scheme_type_t scheme_type) : implicability_function<N, T>(v, scheme_type)
		{}

		plausibility_function(const implicability_function<N, T>& b) : implicability_function<N, T>(b)
		{}

		plausibility_function(const plausibility_function<N, T>& pl) : implicability_function<N, T>(pl)
		{}


		template <class fusion_rule>
		plausibility_function<N, T> fuse_with(const plausibility_function<N, T>& p2) const {
			const fusion_rule fusion;
			return fusion(*this, p2);
		}

		std::unordered_map<size_t, T> contour() {
			std::unordered_map<size_t, T> contour;
			subset singleton = 1;
			T val;

			for(size_t i = 0; i < N; ++i){
				val = (*this)[singleton];
				if (val != 0){
					contour.emplace(i, val);
				}
				singleton <<= 1;
			}
			return contour;
		}

		T* contour_array() {
			T contour[N] = {0};
			subset singleton = 1;

			for(size_t i = 0; i < N; ++i){
				contour[i] = (*this)[singleton];
				singleton <<= 1;
			}
			return contour;
		}

		T at_emptyset() const {
			return (*this)[emptyset];
		}

		T at_fullset() const {
			return (*this)[fullset];
		}

		T operator[](const std::vector<std::string>& labels) const {
			return (*this)[this->outcomes.get_subset(labels)];
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

#endif // EFFICIENT_DST_PLAUSIBILITY_FUNCTION_HPP

