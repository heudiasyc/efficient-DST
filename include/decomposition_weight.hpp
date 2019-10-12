#ifndef EFFICIENT_DST_DECOMPOSITION_WEIGHT_HPP
#define EFFICIENT_DST_DECOMPOSITION_WEIGHT_HPP

#include <mobius_transform.hpp>
#include <zeta_transform.hpp>

namespace efficient_DST{

	template <typename T = double>
	class decomposition_weight : public mobius_transform<T>{

	public:

		decomposition_weight(const powerset_btree<T>& focal_log_sets_values) : mobius_transform<T>(focal_log_sets_values)
		{}

		decomposition_weight(FOD& fod) : mobius_transform<T>(fod)
		{}

		decomposition_weight(const zeta_transform<T>& z) : mobius_transform<T>(z.inversion(operation_t::multiplication))
		{
			this->remove_negligible_values();
			this->normalize();
		}


		virtual ~decomposition_weight(){}


		T at_emptyset() const {
			set_N_value<T>* set_value = this->definition.sub_fod_of_size(0);
			if(set_value)
				return set_value->value;
			else
				return 1;
		}

		T at_fod() const {
			set_N_value<T>* set_value = this->definition.sub_fod_of_size(this->definition.get_FOD_size());
			if(set_value)
				return set_value->value;
			else
				return 1;
		}

		T operator[](const std::vector<std::string>& labels) const {
			return find(this->definition.get_FOD()->to_set(labels));
		}

		T find(const boost::dynamic_bitset<>& key) const {
			set_N_value<T>* set_value = this->definition[key];
			if(set_value)
				return set_value->value;
			else
				return 1;
		}

		void regularize() {
			this->definition.nullify(this->definition[std::vector<std::string>{}]);
			this->normalize();
		}

		void normalize() {
			normalize(this->definition);
		}

		static void normalize(powerset_btree<T>& definition) {
			T prod = 1;
			const std::vector<set_N_value<T>* >& elements = definition.elements();
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

		static void remove_negligible_values(powerset_btree<T>& definition) {
			const std::vector<set_N_value<T>* >& elements = definition.elements();
			for (size_t i = 0; i < elements.size(); ++i) {
				if(mobius_transform<T>::is_equivalent_to_zero(1-elements[i]->value)){
					definition.nullify(elements[i]);
				}
			}
		}
	};

} // namespace efficient_DST

#endif // EFFICIENT_DST_DECOMPOSITION_WEIGHT_HPP
