#ifndef OW_BFT_CONVERTER_TO_COMMONALITY_HPP
#define OW_BFT_CONVERTER_TO_COMMONALITY_HPP

#include <converter.hpp>

namespace ow_bft{

	template <class bft_function, class T = double>
	class converter_to_commonality_special_elements : public converter<T> {
	public:

		converter_to_commonality_special_elements(const T& _precision) : converter<T>(_precision)
		{}

	};

} // namespace ow_bft

#endif // OW_BFT_CONVERTER_TO_COMMONALITY_HPP
