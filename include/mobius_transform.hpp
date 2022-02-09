#ifndef EFFICIENT_DST_MOBIUS_TRANSFORM_HPP
#define EFFICIENT_DST_MOBIUS_TRANSFORM_HPP

#include <powerset_function.hpp>
#include <zeta_transform.hpp>
#include <mobius_inversion.hpp>


namespace efficient_DST{

	template <size_t N, typename T = float>
	class mobius_transform : public powerset_function<N, T>{
	public:
		using typename powerset_function<N, T>::subset;
		using powerset_function<N, T>::emptyset;
		using powerset_function<N, T>::fullset;
//		using powerset_function<N, T>::operator[];

		mobius_transform(
			const sample_space<N>& outcomes,
			const powerset_btree<N, T>& support_values,
			const T& default_value
		) : powerset_function<N, T>(outcomes, support_values, default_value)
		{}

		mobius_transform(const sample_space<N>& outcomes, const T& default_value) : powerset_function<N, T>(outcomes, default_value)
		{}

//		mobius_transform(
//			const zeta_transform<T, N, up_inclusion<N, T>, mobius_additive_operation<T> >& z
//		) : powerset_function<N, T>(z.additive_inversion())
//		{}

		void clear() {
			this->definition.nullify();
		}

		void nullify(const std::vector<std::string>& labels) {
			this->nullify(this->definition[this->outcomes.get_subset(labels)]);
		}

		void nullify(const subset& set) {
			this->definition.nullify(this->definition[set]);
		}

		static void remove_negligible_values(powerset_btree<N, T>& definition, const T& neutral_value) {
			const std::vector<set_N_value<N, T>* >& elements = definition.elements();
			for (size_t i = 0; i < elements.size(); ++i) {
				if(powerset_function<N, T>::is_equivalent_to_zero(elements[i]->value - neutral_value)){
					definition.nullify(elements[i]);
				}
			}
		}

		void remove_negligible_values() {
			remove_negligible_values(this->definition, this->default_value);
		}

		void set_values(const std::unordered_map<subset, T>& values) {
			for (const auto& labels_U_value : values) {
				this->set_value(labels_U_value.first, labels_U_value.second);
			}
		}

		void set_values(const std::unordered_map<std::vector<std::string>, T>& values) {
			for (const auto& labels_U_value : values) {
				this->set_value(labels_U_value.first, labels_U_value.second);
			}
		}

		void set_value(const subset& set, const T& value) {
			this->definition.update_or_insert(set, value);
		}

		void set_value(const std::vector<std::string>& labels, const T& value) {
			this->set_value(this->outcomes.get_subset(labels), value);
		}

		void set_emptyset_value(const T& value) {
			this->set_value(emptyset, value);
		}

		void set_fullset_value(const T& value) {
			this->set_value(fullset, value);
		}
	};

} // namespace efficient_DST

#endif // EFFICIENT_DST_MOBIUS_TRANSFORM_HPP
