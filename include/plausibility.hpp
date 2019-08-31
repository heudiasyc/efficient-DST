#ifndef EFFICIENT_DST_PLAUSIBILITY_HPP
#define EFFICIENT_DST_PLAUSIBILITY_HPP

#include <implicability.hpp>

namespace efficient_DST{

	template <typename T = double>
	class plausibility : public implicability<T> {
	public:

		plausibility(const mass<T>& m) : implicability<T>(m)
		{}

		plausibility(const implicability<T>& b) : implicability<T>(b)
		{}


		template <class fusion_rule>
		plausibility<T> apply(const fusion_rule fusion, const plausibility<T>& p2) const {
			return fusion(*this, p2);
		}

		std::vector<T> get_contour() {
			std::vector<T> contour;
			contour.reserve(this->definition.fod->size());
			boost::dynamic_bitset<> singleton(this->definition.fod->size());

			for(size_t i = 0; i < this->definition.fod->size(); ++i){
				singleton[i] = true;
				contour.emplace_back(find(singleton));
				singleton[i] = false;
			}
			return contour;
		}

		T at_emptyset() const {
			set_N_value<T>* set_value = this->definition.sub_fod_of_size(this->definition.fod->size());
			if(set_value){
				return 1-set_value->value;
			}
			boost::dynamic_bitset<> fod(this->definition.fod->size());
			fod.set();
			return 1-this->compute_aggregation(fod);
		}

		T at_fod() const {
			set_N_value<T>* set_value = this->definition.sub_fod_of_size(0);
			if(set_value){
				return 1-set_value->value;
			}
			boost::dynamic_bitset<> emptyset(this->definition.fod->size());
			return 1-this->compute_aggregation(emptyset);
		}

		T operator[](const std::vector<std::string>& labels) const {
			return find(this->definition.fod->to_set(labels));
		}

		T find(const boost::dynamic_bitset<>& set) const {
			const boost::dynamic_bitset<>& dual_set = FOD::set_negate(set);
			set_N_value<T>* set_value = this->definition[dual_set];
			if(set_value){
				return 1-set_value->value;
			}
			return 1-this->compute_aggregation(dual_set);
		}
	};

} // namespace efficient_DST

#endif // EFFICIENT_DST_PLAUSIBILITY_HPP

