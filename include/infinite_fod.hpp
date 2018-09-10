#ifndef INFINITE_FOD_HPP
#define INFINITE_FOD_HPP

#include <unordered_map>
#include <vector>
#include <fod.hpp>

namespace ow_bft {

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

}	// namespace ow_bft

#endif // INFINITE_FOD_HPP
