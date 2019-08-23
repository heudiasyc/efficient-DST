#ifndef BOOST_BFT_PSEUDO_BFT_FUNCTION_HPP
#define BOOST_BFT_PSEUDO_BFT_FUNCTION_HPP

#include <belief_function.hpp>


namespace efficient_DST{

	template <typename T = double>
	class pseudo_bft_function : public mobius_transform<T> {
	public:
		FOD existence_fod;	// existence fod conditioning the universe fod, i.e. the fod from bft_function
		powerset_btree<T> existence_focal_elements;

		pseudo_bft_function(FOD& _fod) :
			mobius_transform<T>(_fod),
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
}	// namespace efficient_DST
#endif // BOOST_BFT_PSEUDO_BFT_FUNCTION_HPP
