#ifndef EFFICIENT_DST_POWERSET_FUNCTION_HPP
#define EFFICIENT_DST_POWERSET_FUNCTION_HPP

#include <boost/dynamic_bitset.hpp>
#include <vector>
#include <unordered_map>
#include <boost/functional/hash.hpp>
#include <iostream>
#include <iomanip>

#include <fod.hpp>
#include <powerset_btree.hpp>


namespace std {

	template <size_t N>
	struct hash< array<string, N> >{
		size_t operator()(const array<string, N>& key) const{
			size_t seed = 0;
			for (size_t i = 0; i < key.size(); ++i) {
				boost::hash_combine<string>(seed, key[i]);
			}
			return seed;
		}
	};

	template<>
	struct hash< vector<string> >{
		size_t operator()(vector<string>& key) const{
			size_t seed = 0;
			for (size_t i = 0; i < key.size(); ++i) {
				boost::hash_combine<string>(seed, key[i]);
			}
			return seed;
		}
	};
}

namespace efficient_DST{

	template <typename T, size_t N>
	class powerset_function {
	private:
		static constexpr T zero = 0;

	protected:
		/*
		 * Only sets necessary to the definition of this Möbius transform
		 * and their respective images are stored. Their data structure is a binary tree.
		 */
		powerset_btree<T, N> definition;


		T at_emptyset(const T& neutral_value) const {
			set_N_value<T, N>* set_value = this->definition.sub_fod_of_size(0);
			if(set_value)
				return set_value->value;
			else
				return neutral_value;
		}

		T at_fod(const T& neutral_value) const {
			set_N_value<T, N>* set_value = this->definition.sub_fod_of_size(N);
			if(set_value)
				return set_value->value;
			else
				return neutral_value;
		}

		T find(const std::bitset<N>& set, const T& neutral_value) const {
			set_N_value<T, N>* set_value = this->definition[set];
			if(set_value)
				return set_value->value;
			else
				return neutral_value;
		}

	public:
		// allow user to configure the floating-point tolerance
		static const size_t block_size = 1000;
		static constexpr T precision = 1e-10;

		static inline bool is_equivalent_to_zero(const T& value) {
			return (value < zero ? -value : value) <= precision;
		}

		powerset_function (const powerset_btree<T, N>& definition) :
			definition(definition)
		{}

		powerset_function (FOD<N>* fod, const size_t& block_size) :
			definition(fod, block_size)
		{}

		powerset_function (FOD<N>& fod) :
			definition(&fod, block_size)
		{}

		virtual ~powerset_function()
		{}

		// =============================================================================

		const powerset_btree<T, N>& get_definition() const {
			return this->definition;
		}

		const FOD<N>& get_FOD() const {
			return *(this->definition.get_FOD());
		}

		const size_t& get_block_size() const {
			return this->definition.get_block_size();
		}
	};
}		// namespace efficient_DST

#endif // EFFICIENT_DST_POWERSET_FUNCTION_HPP
