#ifndef EFFICIENT_DST_POWERSET_FUNCTION_HPP
#define EFFICIENT_DST_POWERSET_FUNCTION_HPP

#include <vector>
#include <unordered_map>
#include <iostream>
#include <iomanip>

#include <sample_space.hpp>
#include <powerset_btree.hpp>


namespace efficient_DST{

	template <size_t N, typename T = float>
	class powerset_function {
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
		static constexpr T precision = (T) efficient_DST::precision;

		static inline bool is_equivalent_to_zero(const T& value) {
			return (value < 0 ? -value : value) <= precision;
		}

		powerset_function (const sample_space<N>& outcomes, const powerset_btree<N, T>& definition, const T& default_value) :
			outcomes(outcomes),
			definition(definition),
			default_value(default_value)
		{}

		powerset_function (const sample_space<N>& outcomes, const size_t& init_size, const T& default_value) :
			outcomes(outcomes),
			definition(init_size),
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

		T operator[](const std::vector<std::string>& labels) const {
			return (*this)[this->outcomes.get_subset(labels)];
		}

		T operator[](const subset& set) const {
			size_t index = this->definition[set];
			if(index < definition.number_of_nodes())
				return definition.get_node(index).value;
			else
				return this->default_value;
		}

		std::ostream& print(const bool& including_null = false) const {
			return this->definition.print(this->outcomes, including_null);
		}

//		std::ostream& print() const {
//			std::vector<set_N_value<N, T>* > values = this->definition.elements();
//			std::cout << std::endl;
//			for (size_t i = 0; i < values.size(); ++i) {
//				std::cout << values[i]->to_string(this->outcomes) << std::endl;
//			}
//
//			return std::cout;
//		}
	};
}		// namespace efficient_DST

#endif // EFFICIENT_DST_POWERSET_FUNCTION_HPP
