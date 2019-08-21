#ifndef EFFICIENT_DST_CONVERTER_TO_COMMONALITY_HPP
#define EFFICIENT_DST_CONVERTER_TO_COMMONALITY_HPP

#include <converter.hpp>

namespace efficient_DST{

	template <class bft_function, class T = double>
	class converter_to_commonality_special_elements : public converter<T> {
	public:

		converter_to_commonality_special_elements(const T& _precision) : converter<T>(_precision)
		{}

	};

} // namespace efficient_DST

#endif // EFFICIENT_DST_CONVERTER_TO_COMMONALITY_HPP
