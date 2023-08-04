/*
 * Copyright (C) 2019-2023  Maxime Chaveroche (maxime.chaveroche@gmail.com)
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the CeCILL License, either version 2.1 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * CeCILL License for more details.
 * 
 * You should have received a copy of the CeCILL License
 * along with this program. If not, see <http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html>.
 */
 
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
