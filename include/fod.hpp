#ifndef OW_BFT_FOD_HPP
#define OW_BFT_FOD_HPP

#include <boost/dynamic_bitset.hpp>
#include <unordered_map>
#include <vector>
#include <iostream>
#include <fod_element.hpp>
#include <memory_pool.hpp>
#include <powerset_function.hpp>

namespace ow_bft{

	template < typename T = double>
	std::string to_string(const T& n){
		std::ostringstream stm ;
		stm << std::fixed;
		stm << std::setprecision(5);
		stm << n ;
		return stm.str() ;
	}

	static std::string to_string(const std::vector<ow_bft::fod_element*>& fod_elements){
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

	std::ostream& operator<<(std::ostream& os, const std::vector<ow_bft::fod_element*>& fod_elements) {
		os << to_string(fod_elements);
		return os;
	}

	class FOD {
	protected:
		static const size_t fod_elements_block_size = 50;
		//static const size_t powerset_functions_block_size = 5;
		memory_pool<fod_element> fod_element_pool;
		std::vector<fod_element*> fod_elements;
		std::unordered_map<std::string, fod_element*> elements_by_labels;
		//memory_pool<powerset_t> powersets;
		//std::vector<powerset_function* > powerset_functions;

	public:

		FOD() :
			fod_element_pool(fod_elements_block_size)
			//powersets(powersets_block_size)
		{
			//powerset_functions.reserve(powerset_functions_block_size);
		}

		FOD(const std::vector<std::string>& element_labels) : FOD() {
			set_elements(element_labels);
		}

		FOD(const FOD& fod) : FOD(to_labels(fod.elements()))
		{}

		/*
		 * powerset functions must have been constructed with this fod as parameter

		void set_powerset_functions(std::vector<powerset_function* > powerset_functions){
			this->powerset_functions = powerset_functions;
		}
*/
		/*
		 * p_function must have been constructed with this fod as parameter

		void push_back_powerset_function(powerset_function& p_function){
			this->powerset_functions.emplace_back(&p_function);
		}
*/
		/*
		template<typename... Ts>
		powerset_t* create_new_powerset(Ts... args){
			powerset_t* pp = this->powersets.emplace(args...);
			this->powerset_pointers.emplace_back(pp);
			return pp;
		}

		void set_powerset(powerset& managed_powerset){
			this->managed_powerset = &managed_powerset;
		}
		*/
/*
		void erase_powerset_function(const powerset_function& p_function){
			auto it = std::find(this->powerset_functions.begin(), this->powerset_functions.end(), &p_function);
			this->powerset_functions.erase(it);
		}
*/
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
				propagate_push_back(element_label);
			}
		}

		void push_back_elements(const std::vector<std::string>& element_labels){
			for (size_t i = 0; i < element_labels.size(); ++i) {
				push_back(element_labels[i]);
			}
		}
/*
		void push_back_elements_if_absent(const std::vector<std::string> element_labels){
			for (size_t i = 0; i < element_labels.size(); ++i){
				try {
					this->elements_by_labels.at(element_labels[i]);
				}catch(const std::out_of_range& oor){
					push_back(element_labels[i]);
				}
			}
		}

		std::vector<fod_element*> push_back_elements_if_absent_and_return_them(const std::vector<std::string> element_labels){
			std::vector<fod_element*> fes;
			fod_element* fe;
			for (size_t i = 0; i < element_labels.size(); ++i){
				try {
					fe = this->elements_by_labels.at(element_labels[i]);
				}catch(const std::out_of_range& oor){
					push_back(element_labels[i]);
					fe = this->fod_elements.back();
				}
				fes.push_back(fe);
			}
			return fes;
		}*/

		void erase(const std::string& element_label){
			try {
				fod_element* fe = this->elements_by_labels.at(element_label);
				erase(fe);
			}catch(const std::out_of_range& oor){
				std::cerr << "\nCan't erase FOD element of label " << element_label << ":\n";
				std::cerr << "Already deleted.\n";
			}
		}

		void erase(const fod_element* element){
			if(element)
				erase(element->position_in_fod);
		}

		void erase(const size_t position){
			if(position >= 0 && position < this->fod_elements.size()){
				/*
				if(position+1 < this->elements.size())
					this->elements[position+1].set_predecessor(this->elements[position].get_predecessor());
				else
					this->tail_element = this->elements[position].get_predecessor();
				*/
				this->elements_by_labels.erase(this->fod_elements[position]->label);
				/*
				for (size_t i = 0; i < this->powerset_functions.size(); ++i) {
					this->powerset_functions[i]->erase_elements_containing_fod_element(position);
				}*/
				this->fod_element_pool.erase(this->fod_elements[position]);
				this->fod_elements.erase(this->fod_elements.begin() + position);
				update_positions_after(position);
			}else{
				std::cerr << "\nCan't erase fod element of index " << position << " : there's only " << this->fod_elements.size() << " in FOD !\n";
			}
		}

		const std::vector<fod_element*>& elements() const {
			return this->fod_elements;
		}

		const std::vector<fod_element*> to_elements(const boost::dynamic_bitset<>& set) const {
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

		fod_element* to_element(const std::string& element_label) const {
			try {
				return this->elements_by_labels.at(element_label);
			}catch(const std::out_of_range& oor){
				std::cerr << "\nUnknwown label " << element_label << " for FOD element:\n";
				std::cerr << "Out of Range error: " << oor.what() << '\n';
				exit(1);
			}
		}

		const std::vector<fod_element*> to_elements(const std::vector<std::string>& element_labels) const {
			std::vector<fod_element*> fes;
			for (size_t i = 0; i < element_labels.size(); ++i){
				fes.emplace_back(to_element(element_labels[i]));
			}
			return fes;
		}

		static std::vector<std::string> to_labels(const std::vector<fod_element*>& elements) {
			std::vector<std::string> labels;
			for (size_t i = 0; i < elements.size(); ++i){
				labels.emplace_back(elements[i]->label);
			}
			return labels;
		}

		const boost::dynamic_bitset<> to_set(const std::vector<fod_element*>& elements) const {

			boost::dynamic_bitset<> key(this->size());
			for (size_t i = 0; i < elements.size(); ++i){
				if(elements[i]->position_in_fod >= this->size()){
					std::cerr << "\nElement position " << elements[i]->position_in_fod << " in FOD is outside of this FOD.\n";
					exit(1);
				}else{
					key.set(elements[i]->position_in_fod);
				}
			}
			return key;
		}

		const boost::dynamic_bitset<> to_set(const std::vector<std::string>& element_labels) const {
			boost::dynamic_bitset<> key(this->size());
			for (size_t i = 0; i < element_labels.size(); ++i){
				key.set(to_element(element_labels[i])->position_in_fod);
			}
			return key;
		}

		size_t size() const {
			return this->fod_elements.size();
		}

		size_t capacity() const {
			return this->size();
		}

	private:

		void propagate_push_back(const std::string& element_label){
			//this->tail_element = &this->elements.back();
			// defines the hashmap of labels
			this->elements_by_labels.emplace(element_label, this->fod_elements.back());
			/*
			for (int p = 0; p < this->powersets.size(); ++p) {
				this->powersets[p].emplace_back();
			}
			*/
		}

		// updates all fod element index (position in FOD)
		// meant to be called after any modification of this->elements
		void update_positions_after(const size_t index){
			for (size_t i = index; i < this->fod_elements.size(); ++i) {
				this->fod_elements[i]->position_in_fod = i;
			}
			//this->tail_element->chained_update();
		}

		// =============================================================================
		// SUBSET OPERATIONS ===========================================================

	public:

		static const boost::dynamic_bitset<> set_intersection(const boost::dynamic_bitset<>& set1, const boost::dynamic_bitset<>& set2) {
			return set1 & set2;
		}

		const std::vector<fod_element*> set_intersection(const std::vector<fod_element*>& elements1, const std::vector<fod_element*>& elements2) const {
			return to_elements(
				set_intersection(to_set(elements1), to_set(elements2))
			);
		}

		static const boost::dynamic_bitset<> set_union(const boost::dynamic_bitset<>& set1, const boost::dynamic_bitset<>& set2) {
			return set1 | set2;
		}

		const std::vector<fod_element*> set_union(const std::vector<fod_element*>& elements1, const std::vector<fod_element*>& elements2) const {
			return to_elements(
				set_union(to_set(elements1), to_set(elements2))
			);
		}

		static const boost::dynamic_bitset<> set_minus(const boost::dynamic_bitset<>& set1, const boost::dynamic_bitset<>& set2) {
			return set1 - set2;		// set1 & (set1 ^ set2);
		}

		const std::vector<fod_element*> set_minus(const std::vector<fod_element*>& elements1, const std::vector<fod_element*>& elements2) const {
			return to_elements(
				set_minus(to_set(elements1), to_set(elements2))
			);
		}

		static const boost::dynamic_bitset<> set_xor(const boost::dynamic_bitset<>& set1, const boost::dynamic_bitset<>& set2) {
			return set1 ^ set2;
		}

		const std::vector<fod_element*> set_xor(const std::vector<fod_element*>& elements1, const std::vector<fod_element*>& elements2) const {
			return to_elements(
				set_xor(to_set(elements1), to_set(elements2))
			);
		}

		static const boost::dynamic_bitset<> set_negate(const boost::dynamic_bitset<>& set) {
			return ~set;
		}

		const std::vector<fod_element*> set_negate(const std::vector<fod_element*>& elements) const {
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


}	// namespace ow_bft

#endif // OW_BFT_FOD_HPP
