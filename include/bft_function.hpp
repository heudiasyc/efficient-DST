#ifndef OW_BFT_BFT_FUNCTION_HPP
#define OW_BFT_BFT_FUNCTION_HPP

#include <boost/dynamic_bitset.hpp>
#include <vector>
#include <unordered_map>
#include <boost/functional/hash.hpp>
#include <iostream>
#include <iomanip>

#include <detail/is_small.hpp>
#include <fod.hpp>
#include <infinite_fod.hpp>
#include <fod_element.hpp>
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

namespace ow_bft{

	template <class T = double>
	static std::string to_string(const set_N_value<T>& s) {
		return to_string<T>(s.value) + "\t <- " + to_string(s.fod_elements);
	}

	template <class T = double>
	std::ostream& print(std::ostream& os, const std::vector<set_N_value<T>* >& values) {
		for (size_t i = 0; i < values.size(); ++i) {
			os << to_string<T>(*(values[i])) << std::endl;
		}
		return os;
	}

	template <class T = double>
	std::ostream& print(std::ostream& os, const powerset_btree<T>& p) {
		std::vector<set_N_value<T>* > values = p.elements();
		std::cerr << std::endl;

		return print<T>(os, values);
	}

	enum Special_case {degenerate, vacuous};

	template <typename T = double>
	class bft_function {
	protected:

		bool is_equivalent_to_zero(const T& value) const {
			return ow_bft::detail::is_small(value, precision);
		}

	public:
		// allow user to configure the floating-point tolerance
		static const size_t block_size = 100;
		const T precision = 1e-10;
/*
		bft_function(FOD& _fod)
		{
			//this->fod->push_back_powerset(this->focal_elements);
		}
*/
		virtual ~bft_function()
		{}

		// =============================================================================

		virtual T at_emptyset() const = 0;

		virtual T at_fod() const = 0;

		virtual T operator[](const std::vector<std::string>& labels) const = 0;

		virtual T find(const boost::dynamic_bitset<>& key) const = 0;
	};
}		// namespace ow_bft

#endif // OW_BFT_BFT_FUNCTION_HPP
