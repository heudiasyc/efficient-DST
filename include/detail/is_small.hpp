#ifndef EFFICIENT_DST_DETAIL_IS_SMALL_HPP
#define EFFICIENT_DST_DETAIL_IS_SMALL_HPP

#include <boost/version.hpp>
#include <boost/test/floating_point_comparison.hpp>

namespace efficient_DST{
	namespace detail{

		template <typename FPT>
		inline bool is_small(FPT fpv, FPT tolerance){
			//if(fpv != NULL){
				#if BOOST_VERSION >= 105700
					return boost::math::fpc::is_small(fpv, tolerance);
				#else
					return test_tools::check_is_small(fpv, tolerance);
				#endif
			//}else{
			//	return true;
			//}
		}

	} // namespace detail

} // namespace efficient_DST

#endif // EFFICIENT_DST_DETAIL_IS_SMALL_HPP
