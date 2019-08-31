#ifndef EFFICIENT_DST_PIGNISTIC_PROBABILITY_HPP
#define EFFICIENT_DST_PIGNISTIC_PROBABILITY_HPP

#include <mass.hpp>

namespace efficient_DST{

	template <typename T = double>
	class pignistic_probability {
	protected:
		const mass<T> source_mass;

		T compute_aggregation(const boost::dynamic_bitset<>& setA, const bool& singleton) const {
			T sum = 0;

			if (singleton){
				const std::unordered_map<size_t, std::vector<set_N_value<T>* > >& cardinality_map = this->source_mass.get_definition().supersets_of_by_cardinality(setA);

				for (auto kv : cardinality_map) {
					for (size_t i = 0; i < kv.second.size(); ++i) {
						sum += kv.second[i]->value / kv.first;
					}
				}
			}else{
				const std::unordered_map<size_t, std::vector<set_N_value<T>* > >& cardinality_map = this->source_mass.get_definition().elements_by_set_cardinality();

				for (auto kv : cardinality_map) {
					for (size_t i = 0; i < kv.second.size(); ++i) {
						const set_N_value<T>* B = kv.second[i];
						const boost::dynamic_bitset<>& intersection = FOD::set_intersection(setA, B->set);

						sum += B->value * intersection.count() / kv.first;
					}
				}
			}
			return sum;
		}

	public:

		pignistic_probability(const mass<T>& m) : source_mass(m)
		{}

		pignistic_probability(const pignistic_probability<T>& bet_p) : source_mass(bet_p.source_mass)
		{}

		std::vector<T> get_contour() {
			std::vector<T> contour;
			contour.reserve(this->source_mass.get_FOD().size());
			boost::dynamic_bitset<> singleton(this->source_mass.get_FOD().size());

			for(size_t i = 0; i < this->source_mass.get_FOD().size(); ++i){
				singleton[i] = true;
				contour.emplace_back(compute_aggregation(singleton, true));
				singleton[i] = false;
			}
			return contour;
		}

		T at_emptyset() const {
			return 0;
		}

		T at_fod() const {
			boost::dynamic_bitset<> fod(this->source_mass.get_FOD().size());
			fod.set();
			return compute_aggregation(fod, false);
		}

		T operator[](const std::vector<std::string>& labels) const {
			return find(this->source_mass.get_FOD().to_set(labels));
		}

		T find(const boost::dynamic_bitset<>& set) const {
			bool is_singleton = false;
			if (set.count() == 1){
				is_singleton = true;
			}
			return compute_aggregation(set, is_singleton);
		}
	};

} // namespace efficient_DST

#endif // EFFICIENT_DST_PIGNISTIC_PROBABILITY_HPP
