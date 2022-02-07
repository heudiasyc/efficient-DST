#ifndef EFFICIENT_DST_POWERSET_FUNCTION_HPP
#define EFFICIENT_DST_POWERSET_FUNCTION_HPP

#include <vector>
#include <unordered_map>
//#include <boost/functional/hash.hpp>
#include <iostream>
#include <iomanip>

#include <sample_space.hpp>
#include <powerset_btree.hpp>


//namespace std {
//
//	template <size_t N>
//	struct hash< array<string, N> >{
//		size_t operator()(const array<string, N>& key) const{
//			size_t seed = 0;
//			for (size_t i = 0; i < key.size(); ++i) {
//				boost::hash_combine<string>(seed, key[i]);
//			}
//			return seed;
//		}
//	};
//
//	template<>
//	struct hash< vector<string> >{
//		size_t operator()(vector<string>& key) const{
//			size_t seed = 0;
//			for (size_t i = 0; i < key.size(); ++i) {
//				boost::hash_combine<string>(seed, key[i]);
//			}
//			return seed;
//		}
//	};
//}

namespace efficient_DST{

	template <size_t N, typename T = float>
	class powerset_function {
	private:
		static constexpr T zero = 0;

	protected:
		typedef typename sample_space<N>::subset subset;
		sample_space<N> outcomes;
		/*
		 * Only sets necessary to the definition of this Möbius transform
		 * and their respective images are stored. Their data structure is a binary tree.
		 */
		powerset_btree<N, T> definition;
		const T default_value;
		const subset emptyset = 0;
		const subset fullset = ~emptyset;

	public:
		// allow user to configure the floating-point tolerance
		static constexpr T precision = 1e-10;

		static inline bool is_equivalent_to_zero(const T& value) {
			return (value < zero ? -value : value) <= precision;
		}

		powerset_function (const sample_space<N>& outcomes, const powerset_btree<N, T>& definition, const T& default_value) :
			outcomes(outcomes),
			definition(definition),
			default_value(default_value)
		{}

		powerset_function (const sample_space<N>& outcomes, const size_t& block_size, const T& default_value) :
			outcomes(outcomes),
			definition(block_size),
			default_value(default_value)
		{}

		powerset_function (const sample_space<N>& outcomes, const T& default_value) :
			outcomes(outcomes),
			definition(),
			default_value(default_value)
		{}

		virtual ~powerset_function()
		{}

		// =============================================================================

		const powerset_btree<N, T>& get_definition() const {
			return this->definition;
		}

		const sample_space<N>& get_sample_space() const {
			return this->outcomes;
		}

		const T& get_default_value() const {
			return this->default_value;
		}

		T at_emptyset() const {
			return (*this)[emptyset];
		}

		T at_fullset() const {
			return (*this)[fullset];
		}

		T operator[](const subset& set) const {
			set_N_value<N, T>* set_value = this->definition[set];
			if(set_value)
				return set_value->value;
			else
				return this->default_value;
		}

		T operator[](const std::vector<std::string>& labels) const {
			return (*this)[this->outcomes.get_subset(labels)];
		}

		std::ostream& print() const {
			std::vector<set_N_value<N, T>* > values = this->definition.elements();
			std::cout << std::endl;
			for (size_t i = 0; i < values.size(); ++i) {
				std::cout << values[i]->to_string(this->outcomes) << std::endl;
			}

			return std::cout;
		}
	};
}		// namespace efficient_DST

#endif // EFFICIENT_DST_POWERSET_FUNCTION_HPP
