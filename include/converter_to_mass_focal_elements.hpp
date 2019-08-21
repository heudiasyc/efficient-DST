#ifndef EFFICIENT_DST_CONVERTER_TO_MASS_HPP
#define EFFICIENT_DST_CONVERTER_TO_MASS_HPP

#include <converter.hpp>

namespace efficient_DST{

	template <class bft_function, class T = double>
	class converter_to_mass_focal_elements {
	protected:
		static bool is_equivalent_to_zero(const T& value, T precision) {
			return efficient_DST::detail::is_small(value, precision);
		}
	public:

		static powerset_btree<T> convert(const powerset_btree<T>& b_tree, T precision) {

			// initialization
			powerset_btree<T> m_tree(*(b_tree.fod), b_tree.block_size);

			// conversion
			bft_function::to_mass_focal_elements(b_tree, m_tree);

			adjust_values_to_precision(m_tree, precision);

			return m_tree;
		}

		static void convert_without_initialization(powerset_btree<T>& m_tree, T precision) {

			// conversion
			bft_function::to_mass_focal_elements_without_initialization(m_tree);

			adjust_values_to_precision(m_tree, precision);
		}

		static void adjust_values_to_precision(powerset_btree<T>& m_tree, T precision) {

			const std::vector<set_N_value<T>* >& elements = m_tree.elements();
			T sum = 0;
			T val;

			for (size_t i = 0; i < elements.size(); ++i) {
				val = elements[i]->value;
				if(is_equivalent_to_zero(val, precision)){
					m_tree.nullify(elements[i]);
				}else{
					sum += val;
				}
			}

			normalize(m_tree, sum);
		}

		static void normalize(powerset_btree<T>& m_tree, T sum = 0) {

			const std::vector<set_N_value<T>* >& elements = m_tree.elements();

			if(sum == 0){
				for (size_t i = 0; i < elements.size(); ++i) {
					sum += elements[i]->value;
				}
				if(sum == 0){
					std::cerr << "\nSum of mass values equal to 0."
							<< "\nThis means that this mass function is either empty or contains as much positive values as negative values."
							<< "\nEither way, this mass function cannot be normalized into a valid mass function.";
					exit(1);
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

} // namespace efficient_DST

#endif // EFFICIENT_DST_CONVERTER_TO_MASS_HPP
