#ifndef EFFICIENT_DST_IOTA_ELEMENTS_HPP
#define EFFICIENT_DST_IOTA_ELEMENTS_HPP

#include <vector>
#include <unordered_set>
#include <map>
#include <stdlib.h>

#include <algorithm>
#include <math.h>
#include <string>
#include <time.h>
#include <boost/functional/hash.hpp>
#include <iomanip>

#include <powerset_btree.hpp>


namespace efficient_DST{

	enum class irreducible_t: bool { join, meet };

	template <class T, size_t N>
	class iota_elements {
	protected:
		std::vector<std::bitset<N> > sequence;
		std::unordered_set<std::bitset<N> > manifest;
		//std::unordered_map<std::bitset<N>, set_N_value<T, N>* > manifest;

	public:
		const irreducible_t irreducible_type;

		iota_elements(
				const irreducible_t& irreducible_type,
				//const std::vector<std::vector<std::bitset<N>>>& support_by_card) :
				const powerset_btree<float, N>& support) :
			irreducible_type(irreducible_type)
		{
			std::binary_function<size_t, size_t, bool> comp;
			if(irreducible_type == irreducible_t::join){
				comp = std::less<size_t>();
			}else{
				comp = std::greater<size_t>();
			}
			std::map<size_t, std::vector<std::bitset<N> >, comp > iota_elements_card_map;
			//std::unordered_map<size_t, std::vector<std::bitset<N>>> iota_elements_card_map;

			if (this->irreducible_type == irreducible_t::join){
				// iota elements (join-irreducible elements)
				std::bitset<N> singleton = 1;
				for (size_t i = 0; i < N; ++i) {
					const std::vector<set_N_value<T, N>* >& support_supersets = support.supersets_of(singleton);

					if (support_supersets.size() > 0) {
						std::bitset<N> iota_element((const std::bitset<N>&) support_supersets[0]->set);

						for (size_t ii = 1; ii < support_supersets.size(); ++ii) {
							iota_element &= support_supersets[ii]->set;
							if (iota_element == singleton) {
								break;
							}
						}
						bool insertion = this->manifest.emplace(iota_element).second;
						if (insertion){
							size_t cardinality = iota_element.count();
							if (iota_elements_card_map.find(cardinality) == iota_elements_card_map.end()){
								iota_elements_card_map.emplace(cardinality, std::vector<std::bitset<N> >());
								//iota_elements_card_map[cardinality].reserve(N);
							}
							iota_elements_card_map[cardinality].emplace_back(iota_element);
							DEBUG({
								std::clog << "\nNEW INSERTED IOTA ELEMENT :\n";
								std::clog << iota_element << std::endl;
							});
						}
					}
					singleton <<= 1;
				}
			} else {
				// dual iota elements (meet-irreducible elements)
				std::bitset<N> singleton = 1;
				std::bitset<N> singleton_dual = ~singleton;
				for (size_t i = 0; i < N; ++i) {
					const std::vector<set_N_value<T, N>* >& support_subsets = support.subsets_of(singleton_dual);

					if (support_subsets.size() > 0) {
						std::bitset<N> iota_element_dual((const std::bitset<N>&) support_subsets[0]->set);
						//std::clog << "Computing iota element dual associated to :" << singleton_dual << "\n";
						//std::clog << iota_element_dual << std::endl;
						for (size_t ii = 1; ii < support_subsets.size(); ++ii) {
							//std::clog << support_subsets[ii]->set << std::endl;
							iota_element_dual |= support_subsets[ii]->set;
							if (iota_element_dual == singleton_dual) {
								break;
							}
						}
						bool insertion = this->manifest.emplace(iota_element_dual).second;
						if (insertion){
							size_t cardinality = iota_element_dual.count();
							if (iota_elements_card_map.find(cardinality) == iota_elements_card_map.end()){
								iota_elements_card_map.emplace(cardinality, std::vector<std::bitset<N> >());
								//iota_elements_card_map[cardinality].reserve(N);
							}
							iota_elements_card_map[cardinality].emplace_back(iota_element_dual);
							DEBUG({
								std::clog << "\nNEW INSERTED IOTA ELEMENT DUAL :\n";
								std::clog << iota_element_dual << std::endl;
							});
						}
					}
					singleton <<= 1;
					singleton_dual = ~singleton;
				}
			}

			//sort_cardinalities(ordered_cardinalities, iota_elements_card_map, order, N);
			this->sequence.reserve(this->manifest.size());

			for (const auto& c_iota_elements : iota_elements_card_map) {
				const std::vector<std::bitset<N> >& elements = c_iota_elements.second;
				for (size_t i = 0; i < elements.size(); ++i) {
					DEBUG(std::clog << elements[i] << std::endl;);
					this->sequence.emplace_back(elements[i]);
				}
			}
		}


		const std::vector<std::bitset<N> >& get_sequence() const {
			return this->sequence;
		}


		bool contains(const std::bitset<N>& set){
			if (this->manifest.find(set) == this->manifest.end()){
				return false;
			}else{
				return true;
			}
		}
	};
}	// namespace efficient_DST

#endif // EFFICIENT_DST_IOTA_ELEMENTS_HPP
