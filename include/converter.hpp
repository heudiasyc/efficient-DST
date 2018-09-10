#ifndef OW_BFT_CONVERTER_HPP
#define OW_BFT_CONVERTER_HPP

#include <powerset_btree.hpp>
#include <detail/is_small.hpp>

namespace ow_bft{

	template <class T = double>
	class converter {
	protected:

		bool is_equivalent_to_zero(const T& value) const {
			return ow_bft::detail::is_small(value, precision);
		}

		const T precision;

	public:

		converter(const T& _precision) : precision(_precision)
		{}

		virtual ~converter()
		{}

		virtual powerset_btree<T> convert(const powerset_btree<T>& b_tree) const = 0;
	};

} // namespace ow_bft

#endif // OW_BFT_CONVERTER_HPP
