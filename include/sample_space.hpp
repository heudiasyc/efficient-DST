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

#ifndef EFFICIENT_DST_SAMPLE_SPACE_HPP
#define EFFICIENT_DST_SAMPLE_SPACE_HPP

#include <unordered_map>
#include <vector>
#include <iostream>
#include <bitset>


namespace efficient_DST{

//	std::ostream& operator<<(std::ostream& os, const std::vector<std::string>& labels) {
//		os << to_string(labels);
//		return os;
//	}

	template<size_t N>
	class sample_space {
	protected:
		std::string labels[N];
		std::unordered_map<std::string, size_t> indices_by_labels;

	public:

		typedef std::bitset<N> subset;

		sample_space(std::string labels[N]) {
			for (size_t i = 0; i < N; ++i) {
				if (this->indices_by_labels.find(labels[i]) == this->indices_by_labels.end()){
					this->labels[i] = labels[i];
					this->indices_by_labels.emplace(labels[i], i);
				}else{
					std::cerr << "The provided array of outcome labels does not contain " << N << " *distinct* elements.\n";
					exit(1);
				}
			}
		}

		sample_space(std::vector<std::string> labels) {
			for (size_t i = 0; i < N; ++i) {
				if (this->indices_by_labels.find(labels[i]) == this->indices_by_labels.end()){
					this->labels[i] = labels[i];
					this->indices_by_labels.emplace(labels[i], i);
				}else{
					std::cerr << "The provided array of outcome labels does not contain " << N << " *distinct* elements.\n";
					exit(1);
				}
			}
		}

		sample_space(const sample_space<N>& outcomes) {
			for (size_t i = 0; i < N; ++i) {
				this->labels[i] = outcomes.labels[i];
				this->indices_by_labels.emplace(outcomes.labels[i], i);
			}
		}

		const std::string* get_labels() const {
			return this->labels;
		}

		std::vector<std::string> get_labels(const subset& set) const {
			std::vector<std::string> labels;
			//labels.reserve(N);
			subset singleton = 1;
			for (size_t i = 0; i < N; ++i){
				if((set & singleton) != 0)
					labels.emplace_back(this->labels[i]);
				singleton <<= 1;
			}
			return labels;
		}

		std::vector<std::string> get_labels(const size_t& set) const {
			std::vector<std::string> labels;
			//labels.reserve(N);
			size_t singleton = 1;
			for (size_t i = 0; i < N; ++i){
				if((set & singleton) != 0)
					labels.emplace_back(this->labels[i]);
				singleton <<= 1;
			}
			return labels;
		}

//		std::vector<std::string> get_labels(const size_t& index) const {
//			return get_labels((subset) index);
//		}

		size_t get_index(const std::string& label) const {
			try {
				return this->indices_by_labels.at(label);
			}catch(const std::out_of_range& oor){
				std::cerr << "\nUnknwown label " << label << " for outcome:\n";
				std::cerr << "Out of Range error: " << oor.what() << '\n';
				exit(1);
			}
		}

		subset get_subset(const std::vector<std::string>& labels) const {
			subset set = 0;
			subset singleton = 1;
			for (size_t i = 0; i < labels.size(); ++i){
//				set.set(this->get_index(labels[i]));
				set |= singleton << this->get_index(labels[i]);
			}
//			bool set_array[N] = {false};
//			size_t card = labels.size();
//			for (size_t i = 0; i < card; ++i){
//				set_array[this->get_index(labels[i])] = true;
//			}
//			size_t i = 0;
//			size_t j = 0;
//			while (card > 0 && i < N){
//				if (set_array[i]){
//					singleton <<= i - j;
//					set |= singleton;
//					--card;
//					j = i;
//				}
//				++i;
//			}
			return set;
		}

		size_t get_subset_index(const subset& set) const {
			return set.to_ullong();
//			size_t index = 0;
////			for (size_t i = 0; i < labels.size(); ++i){
////				index += pow(2, this->get_index(labels[i]));
////			}
//			subset singleton = 1;
//			size_t power = 1;
//			size_t j = 0;
//			for (size_t i = 0; i < N; ++i){
//				if ((singleton & set) != 0){
//					power *= pow(2, i - j);
//					index += power;
//					j = i;
//				}
//				singleton <<= 1;
//			}
//			return index;
		}

		size_t get_subset_index(const std::vector<std::string>& labels) const {
			size_t index = 0;
			size_t singleton = 1;
			for (size_t i = 0; i < labels.size(); ++i){
//				index += pow(2, this->get_index(labels[i]));
				index |= singleton << this->get_index(labels[i]);
			}
//			bool set_array[N] = {false};
//			size_t card = labels.size();
//			for (size_t i = 0; i < card; ++i){
//				set_array[this->get_index(labels[i])] = true;
//			}
//			size_t i = 0;
//			size_t j = 0;
//			size_t power = 1;
//			while (card > 0 && i < N){
//				if (set_array[i]){
//					power *= pow(2, i - j);
//					index += power;
//					--card;
//					j = i;
//				}
//				++i;
//			}
			return index;
		}

		static inline std::string to_string(const std::vector<std::string>& labels){
			std::string line = "{";
			if(labels.size() > 0){
				size_t i = 0;
				for (; i < labels.size()-1; ++i) {
					line += labels[i] + ", ";
				}
				line += labels[i];
			}
			line += "}";
			return line;
		}

		const std::string to_string(const subset& set) const{
			return to_string(this->get_labels(set));
		}

		const std::string to_string(const size_t& set) const{
			return to_string(this->get_labels((subset) set));
		}
	};


}	// namespace efficient_DST

#endif // EFFICIENT_DST_SAMPLE_SPACE_HPP
