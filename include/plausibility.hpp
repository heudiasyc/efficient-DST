#ifndef EFFICIENT_DST_PLAUSIBILITY_HPP
#define EFFICIENT_DST_PLAUSIBILITY_HPP

#include <implicability.hpp>

namespace efficient_DST{

	template <typename T, size_t N>
	class plausibility : public implicability<T, N> {
	public:

		plausibility(const mass<T, N>& m) : implicability<T, N>(m)
		{}

		plausibility(const implicability<T, N>& b) : implicability<T, N>(b)
		{}


		template <class fusion_rule>
		plausibility<T, N> apply(const plausibility<T, N>& p2) const {
			const fusion_rule fusion;
			return fusion(*this, p2);
		}

		std::vector<T> get_contour() {
			std::vector<T> contour;
			contour.reserve(N);
			std::bitset<N> singleton = 0;

			for(size_t i = 0; i < N; ++i){
				singleton[i] = true;
				contour.emplace_back(find(singleton));
				singleton[i] = false;
			}
			return contour;
		}

		T at_emptyset() const {
			set_N_value<T, N>* set_value = this->definition.sub_fod_of_size(N);
			if(set_value){
				return 1-set_value->value;
			}
			return 1-this->find_non_focal_point_image(~std::bitset<N>(0));
		}

		T at_fod() const {
			set_N_value<T, N>* set_value = this->definition.sub_fod_of_size(0);
			if(set_value){
				return 1-set_value->value;
			}
			return 1-this->find_non_focal_point_image(std::bitset<N>(0));
		}

		T operator[](const std::vector<std::string>& labels) const {
			return find(this->definition.get_FOD()->to_set(labels));
		}

		T find(const std::bitset<N>& set) const {
			const std::bitset<N>& dual_set = ~set;
			set_N_value<T, N>* set_value = this->definition[dual_set];
			if(set_value){
				return 1-set_value->value;
			}
			return 1-this->find_non_focal_point_image(dual_set);
		}
	};

} // namespace efficient_DST

#endif // EFFICIENT_DST_PLAUSIBILITY_HPP

