#ifndef EFFICIENT_DST_PLAUSIBILITY_VECTOR_HPP
#define EFFICIENT_DST_PLAUSIBILITY_VECTOR_HPP

#include <implicability_vector.hpp>

namespace efficient_DST{

	template <size_t N, typename T = float>
	class plausibility_vector : public implicability_vector<N, T> {
	public:

		plausibility_vector(const mass_vector<N, T>& m, const bool& core_reduced = true) : implicability_vector<N, T>(m, core_reduced)
		{}

		plausibility_vector(const weight_vector<N, T>& v, const bool& core_reduced = true) : implicability_vector<N, T>(v, core_reduced)
		{}

		plausibility_vector(const disjunctive_decomposition_vector<N, T>& v, const bool& core_reduced = true) : implicability_vector<N, T>(v, core_reduced)
		{}

		plausibility_vector(const implicability_vector<N, T>& b) : implicability_vector<N, T>(b)
		{}

		plausibility_vector(const plausibility_vector<N, T>& pl) : zeta_transform_vector<down_inclusion<N, T>, N, T>(pl)
		{}


		template <class fusion_rule>
		plausibility_vector<N, T> fuse_with(const plausibility_vector<N, T>& p2) const {
			const fusion_rule fusion;
			return fusion(*this, p2);
		}

		std::unordered_map<size_t, T> contour() {
			std::unordered_map<size_t, T> contour;
			size_t singleton = 1;
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
			size_t singleton = 1;

			for(size_t i = 0; i < N; ++i){
				contour[i] = (*this)[singleton];
				singleton <<= 1;
			}
			return contour;
		}

		T at_emptyset() const {
			return (*this)[0];
		}

		T at_fullset() const {
			return (*this)[(1 << N)-1];
		}

		T operator[](const std::vector<std::string>& labels) const {
			return (*this)[this->outcomes.get_subset_index(labels)];
		}

		T operator[](const size_t& set) const {
			return 1-this->definition[set ^ ((1 << N)-1)];
		}
	};

} // namespace efficient_DST

#endif // EFFICIENT_DST_PLAUSIBILITY_VECTOR_HPP

