#ifndef OW_BFT_IMPLICABILITY_HPP
#define OW_BFT_IMPLICABILITY_HPP

#include <belief.hpp>

namespace ow_bft{

	template <typename T = double>
	class implicability : public belief<T> {
	protected:

		T compute_aggregation_at_emptyset() const {
			return this->mass_equivalent.at_emptyset();
		}

	public:

		implicability(const mass<T>& m) : mass_aggregate<T>(m)
		{
			this->compute_values_for_special_elements(this->mass_equivalent.get_focal_elements(), this->special_elements);
		}

		implicability(const mass_aggregate<T>& ma) : implicability(ma.get_mass_equivalent())
		{}

		implicability(const implicability<T>& b) : mass_aggregate<T>(b.get_mass_equivalent())
		{
			this->special_elements.copy(b.get_special_elements());
		}

		implicability(const FOD& fod) : mass_aggregate<T>(fod)
		{
			this->compute_values_for_special_elements(this->mass_equivalent.get_focal_elements(), this->special_elements);
		}

		implicability(const FOD& fod, const Special_case s_case) : mass_aggregate<T>(fod, s_case)
		{
			this->compute_values_for_special_elements(this->mass_equivalent.get_focal_elements(), this->special_elements);
		}

		template <class fusion_rule>
		implicability<T> apply(const fusion_rule fusion, const implicability<T>& i2) const {
			return fusion(*this, i2);
		}
	};

} // namespace ow_bft

#endif // OW_BFT_IMPLICABILITY_HPP
