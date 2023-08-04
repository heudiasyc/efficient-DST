/*
 * Copyright (C) 2019-2023  Maxime Chaveroche (maxime.chaveroche@gmail.com)
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the CeCILL License, either version 2.1 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * CeCILL License for more details.
 * 
 * You should have received a copy of the CeCILL License
 * along with this program. If not, see <http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html>.
 */
 
#ifndef EFFICIENT_DST_INTERNAL_FOD_HPP
#define EFFICIENT_DST_INTERNAL_FOD_HPP

#include <boost/dynamic_bitset.hpp>
#include <unordered_map>
#include <vector>
#include <iostream>
#include <fod_element.hpp>
#include <fod.hpp>
#include <memory_pool.hpp>
#include <powerset_function.hpp>

class belief_function;

namespace efficient_DST{

	class internal_sub_FOD {
	protected:
		FOD* mother;
		std::vector<fod_element*> fod_elements;
		std::unordered_map<std::string, size_t> indices_by_labels;
		belief_function* managed_bft_function;

	public:

		internal_sub_FOD(const FOD& mother) :
			mother(&mother),
			managed_bft_function(nullptr)
		{}

		/*
		 * powersets must have been constructed with this fod as parameter

		void set_powersets(std::vector<powerset_t* > powersets){
			this->powersets = powersets;
		}
		*/
		/*
		 * powerset must have been constructed with this fod as parameter

		void push_back_powerset(powerset_t& powerset){
			this->powersets.emplace_back(&powerset);
		}*/

		/*
		template<typename... Ts>
		powerset_t* create_new_powerset(Ts... args){
			powerset_t* pp = this->powersets.emplace(args...);
			this->powerset_pointers.emplace_back(pp);
			return pp;
		}*/
		void set_managed_bft_function(belief_function& managed_bft_function){
			this->managed_bft_function = &managed_bft_function;
		}
/*
		void erase_powerset(const powerset_t& powerset){
			auto it = std::find(this->powersets.begin(), this->powersets.end(), &powerset);
			this->powerset_pointers.erase(it);
			this->powersets.erase(&powerset);
		}
*/
		void set_elements(const std::vector<std::string>& element_labels){
			this->fod_elements.clear();
			push_back_elements(element_labels);
		}

		void push_back(const std::string& element_label){
			try {
				this->mother->elements_by_labels.at(element_label);
			}catch(const std::out_of_range& oor){
				this->mother->push_back(element_label);
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
				size_t i = this->indices_by_labels.at(element_label);
				erase(i);
			}catch(const std::out_of_range& oor){
				std::cerr << "Out of Range error: " << oor.what() << '\n';
			}
		}

		std::vector<fod_element*> elements() const {
			return this->fod_elements;
		}

		const std::vector<fod_element*> to_elements(const boost::dynamic_bitset<>& set) const {
			std::vector<fod_element*> fes;
			for (size_t i = 0; i < set.size(); ++i){
				if(set[i])
					fes.emplace_back(this->fod_elements[i]);
			}
			return fes;
		}

		const std::vector<fod_element*> to_elements(const std::vector<std::string>& element_labels) const {
			std::vector<fod_element*> fes;
			for (size_t i = 0; i < element_labels.size(); ++i){
				try {
					fes.emplace_back(this->fod_elements[this->indices_by_labels.at(element_labels[i])]);
				}catch(const std::out_of_range& oor){
					std::cerr << "Out of Range error: " << oor.what() << '\n';
					return fes;
				}
			}
			return fes;
		}
		/*
		std::vector<fod_element*> to_elements(){
			std::vector<fod_element*> fod_elements;
			for (int i = 0; i < this->fod_elements.size(); ++i){
				fod_elements.push_back( &(this->fod_elements[i]) );
			}
			return fod_elements;
		}
		*/
		const boost::dynamic_bitset<> to_set(const std::vector<fod_element*>& elements) const {
			size_t limit = std::min(elements.size(), this->size());
			boost::dynamic_bitset<> key(this->size());
			for (size_t i = 0; i < limit; ++i){
				key.set(elements[i]->position_in_fod);
			}
			return key;
		}

		const boost::dynamic_bitset<> to_set(const std::vector<std::string>& element_labels) const {
			size_t limit = std::min(element_labels.size(), this->size());
			boost::dynamic_bitset<> key(this->size());
			fod_element *fe;
			for (size_t i = 0; i < limit; ++i){
				try {
					fe = this->elements_by_labels.at(element_labels[i]);
					key.set(fe->position_in_fod);
				}catch(const std::out_of_range& oor){
					std::cerr << "Out of Range error: " << oor.what() << '\n';
				}
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
				for (size_t i = 0; i < this->powersets.size(); ++i) {
					this->powersets[i]->erase_elements_containing_fod_element(position);
				}*/
				if(this->managed_powerset)
					this->managed_powerset->erase_elements_containing_fod_element(position);
				this->fod_element_pool.erase(this->fod_elements[position]);
				this->fod_elements.erase(this->fod_elements.begin() + position);
				update_positions_after(position);
			}else{
				std::cerr << "\nCan't erase fod element of index " << position << " : there's only " << this->fod_elements.size() << " in FOD !\n";
			}
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

		bool is_or_is_subset_of(const std::vector<fod_element*>& elements1, const std::vector<fod_element*>& elements2) {
			return is_or_is_subset_of(to_set(elements1), to_set(elements2));
		}

		static bool is_subset_of(const boost::dynamic_bitset<>& set1, const boost::dynamic_bitset<>& set2) {
			return set1.is_proper_subset_of(set2);
		}

		bool is_subset_of(const std::vector<fod_element*>& elements1, const std::vector<fod_element*>& elements2) {
			return is_subset_of(to_set(elements1), to_set(elements2));
		}

		static bool is_or_is_superset_of(const boost::dynamic_bitset<>& set1, const boost::dynamic_bitset<>& set2) {
			return is_or_is_subset_of(set2, set1);
		}

		bool is_or_is_superset_of(const std::vector<fod_element*>& elements1, const std::vector<fod_element*>& elements2) {
			return is_or_is_superset_of(to_set(elements1), to_set(elements2));
		}

		static bool is_superset_of(const boost::dynamic_bitset<>& set1, const boost::dynamic_bitset<>& set2) {
			return is_subset_of(set2, set1);
		}

		bool is_superset_of(const std::vector<fod_element*>& elements1, const std::vector<fod_element*>& elements2) {
			return is_superset_of(to_set(elements1), to_set(elements2));
		}
	};


}	// namespace efficient_DST

#endif // EFFICIENT_DST_INTERNAL_FOD_HPP
