#ifndef BOOST_BFT_PSEUDO_BFT_FUNCTION_HPP
#define BOOST_BFT_PSEUDO_BFT_FUNCTION_HPP

#include <bft_function.hpp>


namespace ow_bft{

	template <typename T = double>
	class pseudo_bft_function : public bft_function<T> {
	public:
		FOD existence_fod;	// existence fod conditioning the universe fod, i.e. the fod from bft_function
		powerset_btree<T> existence_focal_elements;

		pseudo_bft_function(FOD& _fod) :
			bft_function<T>(_fod),
			existence_fod({
				"Xi",
				"not", 		// nothing, i.e. {emptyset}
				"Lambda"	// something, i.e. {Lambda}
			})
		{
			this->existence_focal_elements = powerset_btree<T>(this->existence_fod, this->block_size);
			this->existence_fod.push_back_powerset(this->existence_focal_elements);
		}
	};
}	// namespace ow_bft
#endif // BOOST_BFT_PSEUDO_BFT_FUNCTION_HPP
