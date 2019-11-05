#ifndef EFFICIENT_DST_MOBIUS_TRANSFORM_HPP
#define EFFICIENT_DST_MOBIUS_TRANSFORM_HPP

#include <powerset_function.hpp>
#include <zeta_transform.hpp>

namespace efficient_DST{

	template <typename T, size_t N>
	class mobius_transform : public powerset_function<T, N>{
	public:

		mobius_transform(const powerset_btree<T, N>& support_values) : powerset_function<T, N>(support_values)
		{}

		mobius_transform(FOD<N>& fod) : powerset_function<T, N>(fod)
		{}

		mobius_transform(const zeta_transform<T, N>& z, operation_t transform_operation) : powerset_function<T, N>(z.inversion(transform_operation))
		{}

		void clear() {
			this->definition.nullify();
		}

		void nullify(const std::vector<std::string>& labels) {
			this->definition.nullify(this->definition[labels]);
		}

		static void remove_negligible_values(powerset_btree<T, N>& definition, const T& neutral_value) {
			const std::vector<set_N_value<T, N>* >& elements = definition.elements();
			for (size_t i = 0; i < elements.size(); ++i) {
				if(powerset_function<T, N>::is_equivalent_to_zero(elements[i]->value - neutral_value)){
					definition.nullify(elements[i]);
				}
			}
		}

		void set_values(const std::unordered_map<std::vector<std::string>, T>& values) {
			for (std::pair<std::vector<std::string>, T> labels_U_value : values){
				set_value(labels_U_value.first, labels_U_value.second);
			}
		}

		void set_value(const std::vector<std::string>& labels, const T& value) {
			set_value_directly(this->definition.get_FOD()->to_set(labels), value);
		}

		void set_value_directly(const std::bitset<N>& set, const T& value) {
			this->definition.insert(set, value);
		}

		void set_emptyset_value(const T& value) {
			this->definition.insert(std::bitset<N>(0), value);
		}

		void set_fod_value(const T& value) {
			this->definition.insert(~std::bitset<N>(0), value);
		}
	};

} // namespace efficient_DST

#endif // EFFICIENT_DST_MOBIUS_TRANSFORM_HPP
