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
			for (size_t i = 0; i < labels.size(); ++i){
				set.set(this->get_index(labels[i]));
			}
//			subset singleton = 1;
//			bool set_array[N];
//			for (size_t i = 0; i < labels.size(); ++i){
//				set_array[this->get_index(labels[i])] = true;
//			}
//			for (size_t i = 0; i < N; ++i){
//				if (set_array[i])
//					set |= singleton;
//				singleton <<= 1;
//			}
			return set;
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
	};


}	// namespace efficient_DST

#endif // EFFICIENT_DST_SAMPLE_SPACE_HPP
