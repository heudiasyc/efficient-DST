#ifndef EFFICIENT_DST_DECOMPOSITION_WEIGHT_HPP
#define EFFICIENT_DST_DECOMPOSITION_WEIGHT_HPP

#include <mobius_transform.hpp>
#include <mobius_aggregate.hpp>

namespace efficient_DST{

	template <typename T = double>
	class decomposition_weight : public mobius_transform<T>{

	public:

		decomposition_weight(const powerset_btree<T>& focal_log_sets_values) : mobius_transform<T>(focal_log_sets_values)
		{}

		decomposition_weight(FOD& fod) : mobius_transform<T>(fod)
		{}

		decomposition_weight(const mobius_aggregate<T>& ma) : mobius_transform<T>(ma.inversion(mobius_transformation_form_t::multiplicative))
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

		void normalize() {
			T prod = 1;
			const std::vector<set_N_value<T>* >& elements = this->definition.elements();
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
			const std::vector<set_N_value<T>* >& elements = this->definition.elements();
			for (size_t i = 0; i < elements.size(); ++i) {
				if(this->is_equivalent_to_zero(1-elements[i]->value)){
					this->definition.nullify(elements[i]);
				}
			}
		}
	};

} // namespace efficient_DST

#endif // EFFICIENT_DST_DECOMPOSITION_WEIGHT_HPP
