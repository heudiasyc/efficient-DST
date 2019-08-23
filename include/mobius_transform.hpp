#ifndef EFFICIENT_DST_MOBIUS_TRANSFORM_HPP
#define EFFICIENT_DST_MOBIUS_TRANSFORM_HPP

#include <boost/dynamic_bitset.hpp>
#include <vector>
#include <unordered_map>
#include <boost/functional/hash.hpp>
#include <iostream>
#include <iomanip>

#include <detail/is_small.hpp>
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

	template <class T = double>
	static std::string to_string(const set_N_value<T>& s, const FOD& fod) {
		return to_string<T>(s.value) + "\t <- " + fod.to_string(s.set);
	}

	template <class T = double>
	std::ostream& print(std::ostream& os, const std::vector<set_N_value<T>* >& values, const FOD& fod) {
		for (size_t i = 0; i < values.size(); ++i) {
			os << to_string<T>(*(values[i]), fod) << std::endl;
		}
		return os;
	}

	template <class T = double>
	std::ostream& print(std::ostream& os, const powerset_btree<T>& p) {
		std::vector<set_N_value<T>* > values = p.elements();
		std::cerr << std::endl;

		return print<T>(os, values, *p.fod);
	}

	enum special_case_t {degenerate, vacuous};
	enum class order_relation_t: bool { subset, superset };
	enum class mobius_transformation_form_t: bool { additive, multiplicative };

	template <typename T = double>
	class mobius_transform {
	protected:

		bool is_equivalent_to_zero(const T& value) const {
			return efficient_DST::detail::is_small(value, precision);
		}

	public:
		// allow user to configure the floating-point tolerance
		static const size_t block_size = 100;
		const T precision = 1e-10;
		/*
		 * Only sets necessary to the definition of this Möbius transform
		 * and their respective images are stored. Their data structure is a binary tree.
		 */
		powerset_btree<T> definition;

		mobius_transform (const powerset_btree<T>& definition) :
			definition(definition)
		{}

		mobius_transform (const FOD& fod) :
			definition(fod, block_size)
		{}

		virtual ~mobius_transform()
		{}

		// =============================================================================

		virtual T at_emptyset() const = 0;

		virtual T at_fod() const = 0;

		virtual T operator[](const std::vector<std::string>& labels) const = 0;

		virtual T find(const boost::dynamic_bitset<>& set) const = 0;

		const powerset_btree<T>& get_definition() const {
			return this->definition;
		}

		const FOD& get_FOD() const {
			return *(this->definition.fod);
		}
	};
}		// namespace efficient_DST

#endif // EFFICIENT_DST_MOBIUS_TRANSFORM_HPP
