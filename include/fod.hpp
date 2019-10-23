#ifndef EFFICIENT_DST_FOD_HPP
#define EFFICIENT_DST_FOD_HPP

#include <boost/dynamic_bitset.hpp>
#include <unordered_map>
#include <vector>
#include <iostream>

#include <memory_pool.hpp>


namespace efficient_DST{

	class fod_element{
	public:
		size_t position_in_fod;
		std::string label;

		fod_element(size_t _position_in_fod, std::string _label) :
			position_in_fod(_position_in_fod),
			label(_label)
		{}
	};

	template < typename T = double>
	std::string to_string(const T& n){
		std::ostringstream stm ;
		stm << std::fixed;
		stm << std::setprecision(5);
		stm << n ;
		return stm.str() ;
	}

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

/*	class cardinality_distribution {
	private:
		long double unit_proportion;

		static size_t comb(size_t n, const size_t& k){
			if (k > n) return 0;
			if (k * 2 > n) k = n-k;
			if (k == 0) return 1;

			size_t result = n;
			for (size_t i = 2; i <= k; ++i){
				result *= --n;
				result /= i;
			}
			return result;
		}

		double comb_proportion(size_t n, const size_t& k) const {
			if (k > n) return 0;
			if (k * 2 > n) k = n-k;
			if (k == 0) return this->unit_proportion;

			double result = n/2;
			for (size_t i = 2; i <= k; ++i){
				result *= --n;
				result /= 2*i;
			}
			result *= this->unit_proportion / pow(0.5, k);
			return result;
		}
	public:
		cardinality_distribution(size_t n){
			this->unit_proportion = pow(0.5, n);
		}
	};*/

	class FOD {
	protected:
		static const size_t fod_elements_block_size = 50;
		memory_pool<fod_element> fod_element_pool;
		std::vector<fod_element*> fod_elements;
		std::unordered_map<std::string, fod_element*> elements_by_labels;

		void set_elements(const std::vector<std::string>& element_labels){
			this->fod_elements.clear();
			push_back_elements(element_labels);
		}

		void push_back(const std::string& element_label){
			try {
				this->elements_by_labels.at(element_label);
			}catch(const std::out_of_range& oor){
				fod_element *fe = this->fod_element_pool.emplace(this->size(), element_label);
				this->fod_elements.emplace_back(fe);
				this->elements_by_labels.emplace(element_label, this->fod_elements.back());
			}
		}

		void push_back_elements(const std::vector<std::string>& element_labels){
			for (size_t i = 0; i < element_labels.size(); ++i) {
				push_back(element_labels[i]);
			}
		}

	public:

		FOD(const std::vector<std::string>& element_labels) : fod_element_pool(fod_elements_block_size) {
			set_elements(element_labels);
		}

		FOD(const FOD& fod) : FOD(to_labels(fod.elements()))
		{}

		const std::vector<fod_element*>& elements() const {
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

		std::vector<fod_element*> to_elements(const boost::dynamic_bitset<>& set) const {
			if(set.size() > this->size()){
				std::cerr << "\nError : More elements in bitset of size " << set.size() << " than in FOD of size " << this->size() << ".\n";
				exit(1);
			}
			std::vector<fod_element*> fes;
			for (size_t i = 0; i < set.size(); ++i){
				if(set[i])
					fes.emplace_back(this->fod_elements[i]);
			}
			return fes;
		}

		std::vector<fod_element*> to_elements(const std::vector<std::string>& element_labels) const {
			std::vector<fod_element*> fes;
			for (size_t i = 0; i < element_labels.size(); ++i){
				fes.emplace_back(to_element(element_labels[i]));
			}
			return fes;
		}

		std::vector<std::string> to_labels(const boost::dynamic_bitset<>& set) const {
			if(set.size() > this->size()){
				std::cerr << "\nError : More elements in bitset of size " << set.size() << " than in FOD of size " << this->size() << ".\n";
				exit(1);
			}
			std::vector<std::string> labels;
			for (size_t i = 0; i < set.size(); ++i){
				if(set[i])
					labels.emplace_back(this->fod_elements[i]->label);
			}
			return labels;
		}

		boost::dynamic_bitset<> to_set(const std::vector<fod_element*>& elements) const {

			boost::dynamic_bitset<> set(this->size());
			for (size_t i = 0; i < elements.size(); ++i){
				if(elements[i]->position_in_fod >= this->size()){
					std::cerr << "\nElement position " << elements[i]->position_in_fod << " in FOD is outside of this FOD.\n";
					exit(1);
				}else{
					set.set(elements[i]->position_in_fod);
				}
			}
			return set;
		}

		boost::dynamic_bitset<> to_set(const std::vector<std::string>& element_labels) const {
			boost::dynamic_bitset<> set(this->size());
			for (size_t i = 0; i < element_labels.size(); ++i){
				set.set(to_element(element_labels[i])->position_in_fod);
			}
			return set;
		}

		size_t size() const {
			return this->fod_elements.size();
		}

		static std::vector<std::string> to_labels(const std::vector<fod_element*>& elements) {
			std::vector<std::string> labels;
			for (size_t i = 0; i < elements.size(); ++i){
				labels.emplace_back(elements[i]->label);
			}
			return labels;
		}

		const std::string to_string(const boost::dynamic_bitset<>& set) const{
			return efficient_DST::to_string(to_elements(set));
		}

		static const boost::dynamic_bitset<> set_intersection(const boost::dynamic_bitset<>& set1, const boost::dynamic_bitset<>& set2) {
			return set1 & set2;
		}

		static const boost::dynamic_bitset<> set_union(const boost::dynamic_bitset<>& set1, const boost::dynamic_bitset<>& set2) {
			return set1 | set2;
		}

		static const boost::dynamic_bitset<> set_minus(const boost::dynamic_bitset<>& set1, const boost::dynamic_bitset<>& set2) {
			return set1 - set2;		// set1 & (set1 ^ set2);
		}

		std::vector<fod_element*> set_minus(const std::vector<fod_element*>& elements1, const std::vector<fod_element*>& elements2) const {
			return to_elements(
				set_minus(to_set(elements1), to_set(elements2))
			);
		}

		static const boost::dynamic_bitset<> set_xor(const boost::dynamic_bitset<>& set1, const boost::dynamic_bitset<>& set2) {
			return set1 ^ set2;
		}

		std::vector<fod_element*> set_xor(const std::vector<fod_element*>& elements1, const std::vector<fod_element*>& elements2) const {
			return to_elements(
				set_xor(to_set(elements1), to_set(elements2))
			);
		}

		static const boost::dynamic_bitset<> set_negate(const boost::dynamic_bitset<>& set) {
			return ~set;
		}

		std::vector<fod_element*> set_negate(const std::vector<fod_element*>& elements) const {
			return to_elements(
				set_negate(to_set(elements))
			);
		}

		static bool is_emptyset(const boost::dynamic_bitset<>& set) {
			return !set.any();
		}

		static bool is_singleton(const boost::dynamic_bitset<>& set) {
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

		static bool are_disjoint(const boost::dynamic_bitset<>& set1, const boost::dynamic_bitset<>& set2) {
			size_t size_lim = std::min(set1.size(), set2.size());
			for (size_t i = 0; i < size_lim; ++i) {
				if(set1[i] == set2[i]){
					return false;
				}
			}
			return true;
		}

		static bool is_or_is_subset_of(const boost::dynamic_bitset<>& set1, const boost::dynamic_bitset<>& set2) {
			return set1.is_subset_of(set2);
		}

		bool is_or_is_subset_of(const std::vector<fod_element*>& elements1, const std::vector<fod_element*>& elements2) const {
			return is_or_is_subset_of(to_set(elements1), to_set(elements2));
		}

		static bool is_subset_of(const boost::dynamic_bitset<>& set1, const boost::dynamic_bitset<>& set2) {
			return set1.is_proper_subset_of(set2);
		}

		bool is_subset_of(const std::vector<fod_element*>& elements1, const std::vector<fod_element*>& elements2) const {
			return is_subset_of(to_set(elements1), to_set(elements2));
		}

		static bool is_or_is_superset_of(const boost::dynamic_bitset<>& set1, const boost::dynamic_bitset<>& set2) {
			return is_or_is_subset_of(set2, set1);
		}

		bool is_or_is_superset_of(const std::vector<fod_element*>& elements1, const std::vector<fod_element*>& elements2) const {
			return is_or_is_superset_of(to_set(elements1), to_set(elements2));
		}

		static bool is_superset_of(const boost::dynamic_bitset<>& set1, const boost::dynamic_bitset<>& set2) {
			return is_subset_of(set2, set1);
		}

		bool is_superset_of(const std::vector<fod_element*>& elements1, const std::vector<fod_element*>& elements2) const {
			return is_superset_of(to_set(elements1), to_set(elements2));
		}
	};


}	// namespace efficient_DST

#endif // EFFICIENT_DST_FOD_HPP
