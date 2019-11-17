#ifndef EFFICIENT_DST_DECOMPOSITION_WEIGHT_HPP
#define EFFICIENT_DST_DECOMPOSITION_WEIGHT_HPP

#include <mobius_transform.hpp>
#include <zeta_transform.hpp>

namespace efficient_DST{

	template <typename T, size_t N>
	class decomposition_weight : public mobius_transform<T, N>{

	public:

		decomposition_weight(const powerset_btree<T, N>& focal_log_sets_values) : mobius_transform<T, N>(focal_log_sets_values)
		{}

		decomposition_weight(FOD<N>& fod) : mobius_transform<T, N>(fod)
		{}

		decomposition_weight(const zeta_transform<T, N>& z) : mobius_transform<T, N>(z, operation_t::multiplication)
		{
			this->remove_negligible_values();
			this->normalize();
		}

		virtual ~decomposition_weight(){}


		T at_emptyset() const {
			return powerset_function<T, N>::at_emptyset(1);
		}

		T at_fod() const {
			return powerset_function<T, N>::at_fod(1);
		}

		T operator[](const std::vector<std::string>& labels) const {
			return find(this->definition.get_FOD()->to_set(labels));
		}

		T find(const std::bitset<N>& set) const {
			return powerset_function<T, N>::find(set, 1);
		}

		void regularize() {
			this->definition.nullify(this->definition[std::bitset<N>(0)]);
			this->normalize();
		}

		void normalize() {
			normalize(this->definition);
		}

		static void normalize(powerset_btree<T, N>& definition) {
			T prod = 1;
			const std::vector<set_N_value<T, N>* >& elements = definition.elements();
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

		void remove_negligible_values() {
			remove_negligible_values(this->definition);
		}

		static void remove_negligible_values(powerset_btree<T, N>& definition) {
			mobius_transform<T, N>::remove_negligible_values(definition, 1);
		}
	};

} // namespace efficient_DST

#endif // EFFICIENT_DST_DECOMPOSITION_WEIGHT_HPP
