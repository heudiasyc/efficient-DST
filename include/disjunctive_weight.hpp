#ifndef EFFICIENT_DST_DISJUNCTIVE_WEIGHT_HPP
#define EFFICIENT_DST_DISJUNCTIVE_WEIGHT_HPP

#include <weight_function.hpp>
#include <powerset_function.hpp>
#include <zeta_transform.hpp>

namespace efficient_DST{

	template <size_t N, typename T = float>
	class disjunctive_weight : public weight_function<N, T> {
	public:
		using typename powerset_function<N, T>::subset;
		using powerset_function<N, T>::emptyset;
		using powerset_function<N, T>::fullset;

		disjunctive_weight(const disjunctive_weight<N, T>& v) : weight_function<N, T>(v.outcomes, v.definition)
		{}

		disjunctive_weight(
			const sample_space<N>& outcomes,
			const powerset_btree<N, T>& support
		) : weight_function<N, T>(outcomes, support)
		{}

		disjunctive_weight(const sample_space<N>& outcomes) : weight_function<N, T>(outcomes)
		{}

		disjunctive_weight(
			const zeta_transform<down_inclusion<N, T>, N, T >& b
		) : weight_function<N, T>(b.get_sample_space(), b.inversion(operation_type_t::multiplication))
		{
//			invert_values(this->definition);
		}

//		/*
//		 * The disjunctive weight function is the inverse of the multiplicative Möbius transform of the implicability function.
//		 * So, invert its values before computing it.
//		 */
//		powerset_btree<N, T> inverted_definition() const {
//			powerset_btree<N, T> inverted_definition(this->definition);
//			invert_values(inverted_definition);
//			return inverted_definition;
//		}

		void assign(const subset& set, const T& value) {
			if (set != emptyset){
				this->definition.update_or_insert(set, 1/value);
				// The following part ensures that v(emptyset) is defined (as it is not when directly building v).
				// Indeed, v(emptyset) may be required for the computation of the implicability function.
				// v(emptyset)^{-1} is equal to the product of all other weights.
				set_N_value<N, T>* emptyset = this->definition.empty_set();
				if(emptyset){
					emptyset->value *= value;
				}else{
					this->definition.insert(emptyset, value);
				}
			}else{
				std::cout << "Cannot directly assign the emptyset value of a disjunctive decomposition. Ignoring command.\n";
			}
		}

		void assign(const std::unordered_map<subset, T>& values) {
			for (const auto& labels_U_value : values) {
				this->assign(labels_U_value.first, labels_U_value.second);
			}
		}

		void assign(const std::unordered_map<std::vector<std::string>, T>& values) {
			for (const auto& labels_U_value : values) {
				this->assign(labels_U_value.first, labels_U_value.second);
			}
		}

		void assign(const std::vector<std::string>& labels, const T& value) {
			this->assign(this->outcomes.get_subset(labels), value);
		}

		void assign_emptyset(const T& value) {
			this->assign(emptyset, value);
		}

		void assign_fullset(const T& value) {
			this->assign(fullset, value);
		}

		void nullify(const std::vector<std::string>& labels) {
			this->nullify(this->outcomes.get_subset(labels));
		}

		void nullify(const subset& set) {
			if (set != emptyset){
				set_N_value<N, T>* A = this->definition[set];
				if(A){
					set_N_value<N, T>* emptyset = this->definition.empty_set();
					if(emptyset){
						emptyset->value *= A->value;
					}else{
						// The following part ensures that v(emptyset) is defined (as it is not when directly building v).
						// Indeed, v(emptyset) may be required for the computation of the implicability function.
						// v(emptyset)^{-1} is equal to the product of all other weights.
						this->definition.insert(emptyset, A->value);
					}
					this->definition.nullify(A);
				}
			}else{
				std::cout << "Cannot directly assign the emptyset value of a disjunctive decomposition. Ignoring command.\n";
			}
		}

		void compute_emptyset_value_from_definition(){
			compute_emptyset_value_from_definition(this->definition);
		}

		static void compute_emptyset_value_from_definition(powerset_btree<N, T>& definition){
			const std::vector<set_N_value<N, T>* >& focal_log_elements = definition.elements();
			T val = 1;
			for (size_t i = 0; i < focal_log_elements.size(); ++i){
				if (focal_log_elements[i]->set != emptyset)
					val /= focal_log_elements[i]->value;
			}
			definition.update_or_insert(emptyset, val);
		}

		template <class fusion_rule>
		disjunctive_weight<N, T> fuse_with(const disjunctive_weight<N, T>& v2) const {
			const fusion_rule fusion;
			return fusion(*this, v2);
		}

//	protected:
//
//		static void invert_values(powerset_btree<N, T>& values) {
//			std::vector<set_N_value<N, T>* > elements = values.elements();
//			for (size_t i = 0; i < elements.size(); ++i){
//				elements[i]->value = 1 / elements[i]->value;
//			}
//		}
	};

} // namespace efficient_DST

#endif // EFFICIENT_DST_DISJUNCTIVE_WEIGHT_HPP
