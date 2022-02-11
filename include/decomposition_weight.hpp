#ifndef EFFICIENT_DST_DECOMPOSITION_WEIGHT_HPP
#define EFFICIENT_DST_DECOMPOSITION_WEIGHT_HPP

#include <mobius_transform.hpp>
#include <zeta_transform.hpp>

namespace efficient_DST{

	template <size_t N, typename T = float>
	class decomposition_weight : public mobius_transform<N, T>{
	public:
		using typename powerset_function<N, T>::subset;
		using powerset_function<N, T>::emptyset;
		using powerset_function<N, T>::fullset;

		decomposition_weight(
			const sample_space<N>& outcomes,
			const powerset_btree<N, T>& log_focal_sets
		) : mobius_transform<N, T>(outcomes, log_focal_sets, 1)
		{
			this->remove_negligible_values();
			this->normalize();
		}

		decomposition_weight(
			const sample_space<N>& outcomes
		) : mobius_transform<N, T>(outcomes, 1)
		{}

		virtual ~decomposition_weight(){}


		T at_emptyset() const {
			return (*this)[emptyset];
		}

		T at_fullset() const {
			return (*this)[fullset];
		}

		T operator[](const std::vector<std::string>& labels) const {
			return (*this)[this->outcomes.get_subset(labels)];
		}

		T operator[](const subset& set) const {
			set_N_value<N, T>* set_value = this->definition[set];
			if(set_value)
				return 1/set_value->value;
			else
				return this->default_value;
		}

		std::ostream& print() const {
			std::vector<set_N_value<N, T>* > values = this->definition.elements();
			std::cout << std::endl;
			for (size_t i = 0; i < values.size(); ++i) {
				std::cout << set_N_value<N, T>::to_string(1/values[i]->value) + "\t <- " + this->outcomes.to_string(values[i]->set) << std::endl;
			}

			return std::cout;
		}

		void regularize() {
			this->definition.nullify(this->definition[emptyset]);
			this->normalize();
		}

		void normalize() {
			normalize(this->definition);
		}

		static void normalize(powerset_btree<N, T>& definition) {
			T prod = 1;
			const std::vector<set_N_value<N, T>* >& elements = definition.elements();
			for (size_t i = 0; i < elements.size(); ++i) {
				prod *= elements[i]->value;
			}
			if(prod == 0){
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

//		void remove_negligible_values() {
//			remove_negligible_values(this->definition);
//		}
//
//		static void remove_negligible_values(powerset_btree<N, T>& definition) {
//			mobius_transform<N, T>::remove_negligible_values(definition, 1);
//		}
	};

} // namespace efficient_DST

#endif // EFFICIENT_DST_DECOMPOSITION_WEIGHT_HPP
