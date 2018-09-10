#ifndef OW_BFT_CONVERTER_TO_MASS_HPP
#define OW_BFT_CONVERTER_TO_MASS_HPP

#include <converter.hpp>

namespace ow_bft{

	template <class bft_function, class T = double>
	class converter_to_mass_focal_elements : public converter<T> {
	public:

		converter_to_mass_focal_elements(const T& _precision) : converter<T>(_precision)
		{}

		powerset_btree<T> convert(const powerset_btree<T>& b_tree) const {

			// initialization
			powerset_btree<T> m_tree(*(b_tree.fod), b_tree.block_size);

			// conversion
			bft_function::to_mass_focal_elements(b_tree, m_tree);

			adjust_values_to_precision(m_tree);

			return m_tree;
		}

		void convert_without_initialization(powerset_btree<T>& m_tree) const {

			// conversion
			bft_function::to_mass_focal_elements_without_initialization(m_tree);

			adjust_values_to_precision(m_tree);
		}

		void adjust_values_to_precision(powerset_btree<T>& m_tree) const {

			const std::vector<set_N_value<T>* >& elements = m_tree.elements();
			T sum = 0;
			T val;

			for (size_t i = 0; i < elements.size(); ++i) {
				val = elements[i]->value;
				if(this->is_equivalent_to_zero(val)){
					m_tree.nullify(elements[i]);
				}else{
					sum += val;
				}
			}

			if(sum != 1){
				// normalize
				const std::vector<set_N_value<T>* >& elements = m_tree.elements();

				for (size_t i = 0; i < elements.size(); ++i) {
					elements[i]->value /= sum;
				}
			}
		}
	};

} // namespace ow_bft

#endif // OW_BFT_CONVERTER_TO_MASS_HPP
