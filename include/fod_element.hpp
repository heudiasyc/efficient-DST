#ifndef OW_BFT_FOD_ELEMENT_HPP
#define OW_BFT_FOD_ELEMENT_HPP

namespace ow_bft{

	class fod_element{
	public:
		size_t position_in_fod;
		std::string label;

		fod_element(size_t _position_in_fod, std::string _label) :
			position_in_fod(_position_in_fod),
			label(_label)
		{}
	};

}	// namespace ow_bft

#endif // OW_BFT_FOD_ELEMENT_HPP
