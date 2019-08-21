#ifndef EFFICIENT_DST_CONVERTER_HPP
#define EFFICIENT_DST_CONVERTER_HPP

#include <powerset_btree.hpp>
#include <detail/is_small.hpp>

namespace efficient_DST{

	template <class T = double>
	class converter {
	protected:

		bool is_equivalent_to_zero(const T& value) const {
			return efficient_DST::detail::is_small(value, precision);
		}

		const T precision;

	public:

		converter(const T& _precision) : precision(_precision)
		{}

		virtual ~converter()
		{}

		virtual powerset_btree<T> convert(const powerset_btree<T>& b_tree) const = 0;
	};

} // namespace efficient_DST

#endif // EFFICIENT_DST_CONVERTER_HPP
