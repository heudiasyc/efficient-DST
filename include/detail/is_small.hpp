#ifndef OW_BFT_DETAIL_IS_SMALL_HPP
#define OW_BFT_DETAIL_IS_SMALL_HPP

#include <boost/version.hpp>
#include <boost/test/floating_point_comparison.hpp>

namespace ow_bft{
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

} // namespace ow_bft

#endif // OW_BFT_DETAIL_IS_SMALL_HPP
