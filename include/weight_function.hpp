#ifndef EFFICIENT_DST_WEIGHT_FUNCTION_HPP
#define EFFICIENT_DST_WEIGHT_FUNCTION_HPP

#include <mobius_transform.hpp>
#include <zeta_transform.hpp>
#include <weight_vector.hpp>

namespace efficient_DST{

	template <size_t N, typename T = float>
	class weight_function : public mobius_transform<N, T>{
	public:
		using typename powerset_function<N, T>::subset;
		using powerset_function<N, T>::emptyset;
		using powerset_function<N, T>::fullset;

		weight_function(const weight_function<N, T>& w) : mobius_transform<N, T>(w.outcomes, w.definition, 1)
		{}

		weight_function(const weight_vector<N, T>& w) : mobius_transform<N, T>(w.outcomes, w.definition, 1)
		{}

		weight_function(
			const sample_space<N>& outcomes,
			const powerset_btree<N, T>& log_focal_sets
		) : mobius_transform<N, T>(outcomes, log_focal_sets, 1)
		{
			this->remove_negligible_values();
			this->normalize();
		}

		weight_function(
			const sample_space<N>& outcomes
		) : mobius_transform<N, T>(outcomes, 1)
		{}

		weight_function(const zeta_transform<up_inclusion<N, T>, N, T >& q) : weight_function<N, T>(q.get_sample_space(), q.inversion(operation_type_t::multiplication))
		{}

		weight_function(const zeta_transform<down_inclusion<N, T>, N, T >& b) : weight_function<N, T>(b.get_sample_space(), b.inversion(operation_type_t::multiplication))
		{}


		void regularize() {
			this->definition.nullify(this->definition[emptyset]);
			this->normalize();
		}

		void normalize() {
			normalize(this->definition);
		}

		static void normalize(powerset_btree<N, T>& definition) {
			T prod = 1;
			const std::vector<set_N_value<N, T>* >& elements = definition._elements();
			for (size_t i = 0; i < elements.size(); ++i) {
				prod *= elements[i]->value;
			}
			if(prod == 0){
				definition.print(true);
				std::cerr << "\nProduct of weight values equal to 0."
						<< "\nThis means that at least one of these weights is null, which is not allowed in a weight function." << std::endl;
				exit(1);
			}
			if(prod != 1){
				// normalize
				T factor = pow(prod, 1/elements.size());
				for (size_t i = 0; i < elements.size(); ++i) {
					elements[i]->value /= factor;
				}
			}
		}

//		std::ostream& print() const {
//			std::vector<set_N_value<N, T>* > values = this->definition.elements();
//			std::cout << std::endl;
//			for (size_t i = 0; i < values.size(); ++i) {
//				std::cout << set_N_value<N, T>::to_string(values[i]->value) + "\t <- " + this->outcomes.to_string(values[i]->set) << std::endl;
//			}
//
//			return std::cout;
//		}

//		void remove_negligible_values() {
//			remove_negligible_values(this->definition);
//		}
//
//		static void remove_negligible_values(powerset_btree<N, T>& definition) {
//			mobius_transform<N, T>::remove_negligible_values(definition, 1);
//		}
	};

} // namespace efficient_DST

#endif // EFFICIENT_DST_WEIGHT_FUNCTION_HPP
