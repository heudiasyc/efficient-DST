#ifndef EFFICIENT_DST_FOD_ELEMENT_HPP
#define EFFICIENT_DST_FOD_ELEMENT_HPP

namespace efficient_DST{

	class fod_element{
	public:
		size_t position_in_fod;
		std::string label;

		fod_element(size_t _position_in_fod, std::string _label) :
			position_in_fod(_position_in_fod),
			label(_label)
		{}
	};

}	// namespace efficient_DST

#endif // EFFICIENT_DST_FOD_ELEMENT_HPP
