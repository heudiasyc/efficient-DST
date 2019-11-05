#ifndef EFFICIENT_DST_FOD_HPP
#define EFFICIENT_DST_FOD_HPP

#include <boost/dynamic_bitset.hpp>
#include <unordered_map>
#include <vector>
#include <iostream>

#include <memory_pool.hpp>


namespace efficient_DST{

	enum class order_t: bool {ascending, descending};

	struct fod_element{
		size_t position_in_fod;
		std::string label;

		fod_element(size_t _position_in_fod, std::string _label) :
			position_in_fod(_position_in_fod),
			label(_label)
		{}

		fod_element() :
			position_in_fod(0),
			label("")
		{}
	};

	static std::string to_string(const std::vector<efficient_DST::fod_element*>& fod_elements){
		std::string labels = "{";
		size_t size = fod_elements.size();
		if(size > 0){
			size_t i = 0;
			for (; i < fod_elements.size()-1; ++i) {
				labels += fod_elements[i]->label + ", ";
			}
			labels += fod_elements[i]->label;
		}
		labels += "}";
		return labels;
	}

	std::ostream& operator<<(std::ostream& os, const std::vector<efficient_DST::fod_element*>& fod_elements) {
		os << to_string(fod_elements);
		return os;
	}

	template<size_t N>
	class FOD {
	protected:
		fod_element fod_elements[N];
		std::unordered_map<std::string, fod_element*> elements_by_labels;

	public:

		typedef std::bitset<N> subset;

		FOD(std::string element_labels[N]) {
			for (size_t i = 0; i < N; ++i) {
				if (this->elements_by_labels.find(element_labels[i]) == this->elements_by_labels.end()){
					this->fod_elements[i].position_in_fod = i;
					this->fod_elements[i].label = element_labels[i];
					this->elements_by_labels.emplace(element_labels[i], &this->fod_elements[i]);
				}else{
					std::cerr << "The provided array of FOD elements does not contain " << N << " *distinct* elements.\n";
					exit(1);
				}
			}
		}

		FOD(const FOD<N>& fod) {
			for (size_t i = 0; i < N; ++i) {
				this->fod_elements[i].position_in_fod = i;
				this->fod_elements[i].label = fod.fod_elements[i].label;
				this->elements_by_labels.emplace(fod.fod_elements[i].label, &this->fod_elements[i]);
			}
		}

		const fod_element& elements() const {
			return this->fod_elements;
		}

		fod_element* to_element(const std::string& element_label) const {
			try {
				return this->elements_by_labels.at(element_label);
			}catch(const std::out_of_range& oor){
				std::cerr << "\nUnknwown label " << element_label << " for FOD element:\n";
				std::cerr << "Out of Range error: " << oor.what() << '\n';
				exit(1);
			}
		}

		std::vector<fod_element*> to_elements(const subset& set) const {
			std::vector<fod_element*> fes;
			fes.reserve(N);
			for (size_t i = 0; i < N; ++i){
				if(set[i])
					fes.emplace_back((fod_element*) &this->fod_elements[i]);
			}
			return fes;
		}

		std::vector<fod_element*> to_elements(const std::vector<std::string>& element_labels) const {
			std::vector<fod_element*> fes;
			fes.reserve(element_labels.size());
			for (size_t i = 0; i < element_labels.size(); ++i){
				fes.emplace_back(to_element(element_labels[i]));
			}
			return fes;
		}

		std::vector<std::string> to_labels(const subset& set) const {
			std::vector<std::string> labels;
			labels.reserve(N);
			for (size_t i = 0; i < N; ++i){
				if(set[i])
					labels.emplace_back(this->fod_elements[i].label);
			}
			return labels;
		}

		subset to_set(const std::vector<fod_element*>& elements) const {

			subset set = 0;
			for (size_t i = 0; i < elements.size(); ++i){
				if(elements[i]->position_in_fod >= N){
					std::cerr << "\nElement position " << elements[i]->position_in_fod << " in FOD is outside of this FOD.\n";
					exit(1);
				}else{
					set.set(elements[i]->position_in_fod);
				}
			}
			return set;
		}

		subset to_set(const std::vector<std::string>& element_labels) const {
			subset set = 0;
			for (size_t i = 0; i < element_labels.size(); ++i){
				set.set(to_element(element_labels[i])->position_in_fod);
			}
			return set;
		}

		static inline subset empty_set() {
			return std::bitset<N>(0);
		}

		static inline subset fod_set() {
			return ~std::bitset<N>(0);
		}

		static const size_t size() {
			return N;
		}

		static std::vector<std::string> to_labels(const std::vector<fod_element*>& elements) {
			std::vector<std::string> labels;
			labels.reserve(elements.size());
			for (size_t i = 0; i < elements.size(); ++i){
				labels.emplace_back(elements[i]->label);
			}
			return labels;
		}

		const std::string to_string(const subset& set) const{
			return efficient_DST::to_string(to_elements(set));
		}

		static const subset set_intersection(const subset& set1, const subset& set2) {
			return set1 & set2;
		}

		static const subset set_union(const subset& set1, const subset& set2) {
			return set1 | set2;
		}

		static const subset set_minus(const subset& set1, const subset& set2) {
			return set1 - set2;		// set1 & (set1 ^ set2);
		}

		std::vector<fod_element*> set_minus(const std::vector<fod_element*>& elements1, const std::vector<fod_element*>& elements2) const {
			return to_elements(
				set_minus(to_set(elements1), to_set(elements2))
			);
		}

		static const subset set_xor(const subset& set1, const subset& set2) {
			return set1 ^ set2;
		}

		std::vector<fod_element*> set_xor(const std::vector<fod_element*>& elements1, const std::vector<fod_element*>& elements2) const {
			return to_elements(
				set_xor(to_set(elements1), to_set(elements2))
			);
		}

		static const subset set_negate(const subset& set) {
			return ~set;
		}

		std::vector<fod_element*> set_negate(const std::vector<fod_element*>& elements) const {
			return to_elements(
				set_negate(to_set(elements))
			);
		}

		template<typename value_type>
		void sort_cardinalities(std::vector<size_t>& ordered_cardinalities, const std::unordered_map<size_t, value_type>& map, order_t order) const {
			const size_t& F_card = map.size();
			ordered_cardinalities.reserve(F_card);

			if(F_card > 0 && N <= F_card*log2(F_card)){
				if (order == order_t::ascending){
					for(size_t c = 0; c <= N; ++c) {
						if(map.find(c) != map.end()){
							ordered_cardinalities.push_back(c);
						}
					}
				}else{
					for(size_t c = N; c > 0; --c) {
						if(map.find(c) != map.end()){
							ordered_cardinalities.push_back(c);
						}
					}
					if(map.find(0) != map.end()){
						ordered_cardinalities.push_back(0);
					}
				}
			}else{
				for(auto kv : map) {
					ordered_cardinalities.push_back(kv.first);
				}
				if (order == order_t::ascending){
					// sort cardinalities in ascending order
					std::sort(ordered_cardinalities.begin(), ordered_cardinalities.end());
				}else{
					// sort cardinalities in descending order
					std::sort(ordered_cardinalities.begin(), ordered_cardinalities.end(), std::greater<size_t>());
				}
			}
		}

		static bool is_emptyset(const subset& set) {
			return !set.any();
		}

		static bool is_singleton(const subset& set) {
			size_t i = 0;
			size_t nb_of_elements = 0;
			while(i < set.size()){
				if(set[i]){
					if(++nb_of_elements > 1)
						return false;
				}
				++i;
			}
			return nb_of_elements != 0;
		}

		static bool are_disjoint(const subset& set1, const subset& set2) {
			size_t size_lim = std::min(set1.size(), set2.size());
			for (size_t i = 0; i < size_lim; ++i) {
				if(set1[i] == set2[i]){
					return false;
				}
			}
			return true;
		}

		static bool is_subset_of(const subset& set1, const subset& set2) {
			//return set1.is_subset_of(set2);
			return (set1 & set2) == set1;
		}

		bool is_subset_of(const std::vector<fod_element*>& elements1, const std::vector<fod_element*>& elements2) const {
			return is_subset_of(to_set(elements1), to_set(elements2));
		}

		static bool is_strict_subset_of(const subset& set1, const subset& set2) {
			//return set1.is_proper_subset_of(set2);
			return set1 != set2 && (set1 & set2) == set1;
		}

		bool is_strict_subset_of(const std::vector<fod_element*>& elements1, const std::vector<fod_element*>& elements2) const {
			return is_strict_subset_of(to_set(elements1), to_set(elements2));
		}

		static bool is_superset_of(const subset& set1, const subset& set2) {
			return is_subset_of(set2, set1);
		}

		bool is_superset_of(const std::vector<fod_element*>& elements1, const std::vector<fod_element*>& elements2) const {
			return is_superset_of(to_set(elements1), to_set(elements2));
		}

		static bool is_strict_superset_of(const subset& set1, const subset& set2) {
			return is_strict_subset_of(set2, set1);
		}

		bool is_strict_superset_of(const std::vector<fod_element*>& elements1, const std::vector<fod_element*>& elements2) const {
			return is_strict_superset_of(to_set(elements1), to_set(elements2));
		}
	};


}	// namespace efficient_DST

#endif // EFFICIENT_DST_FOD_HPP
