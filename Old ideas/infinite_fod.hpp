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
 
#ifndef INFINITE_FOD_HPP
#define INFINITE_FOD_HPP

#include <unordered_map>
#include <vector>
#include <fod.hpp>

namespace efficient_DST {

class infinite_FOD : public FOD {
protected:
	FOD* mirror_fod;		// FOD with a finite cardinality that virtually limits the one of this FOD
public:

	infinite_FOD() : FOD(), mirror_fod(nullptr){
		FOD::push_back("infinity");
	}

	infinite_FOD(FOD& _mirror_fod) : infinite_FOD()
	{
		this->mirror_fod = &_mirror_fod;
	}

	infinite_FOD(std::vector<std::string> element_labels) : infinite_FOD(){ // @suppress("Class members should be properly initialized")
		FOD::push_back_elements(element_labels);
	}

	void erase(std::string element_label){
		if(element_label != "infinity")
			erase(element_label);
	}

	void erase(fod_element element){
		if(element.position_in_fod != 0)
			erase(element);
	}

	size_t capacity(){
		return 0;
	}
};

}	// namespace efficient_DST

#endif // INFINITE_FOD_HPP
