#ifndef EFFICIENT_DST_CONJUNCTIVE_WEIGHT_HPP
#define EFFICIENT_DST_CONJUNCTIVE_WEIGHT_HPP

#include <decomposition_weight.hpp>
#include <powerset_function.hpp>
#include <zeta_transform.hpp>

namespace efficient_DST{

	template <size_t N, typename T = float>
	class conjunctive_weight : public decomposition_weight<N, T> {
	public:
		using typename powerset_function<N, T>::subset;
		using powerset_function<N, T>::emptyset;
		using powerset_function<N, T>::fullset;

		conjunctive_weight(const conjunctive_weight<N, T>& w) : decomposition_weight<N, T>(w.outcomes, w.definition)
		{}

		conjunctive_weight(
			const sample_space<N>& outcomes,
			const powerset_btree<N, T>& support
		) : decomposition_weight<N, T>(outcomes, support)
		{}

		conjunctive_weight(const sample_space<N>& outcomes) : decomposition_weight<N, T>(outcomes)
		{}

		conjunctive_weight(
			const zeta_transform<up_inclusion<N, T>, N, T>& q
		) : decomposition_weight<N, T>(q.get_sample_space(), q.inversion(operation_type_t::multiplication))
		{
//			invert_values(this->definition);
		}

//		/*
//		 * The conjunctive weight function is the inverse of the multiplicative Möbius transform of the commonality function.
//		 * So, invert its values before computing it.
//		 */
//		powerset_btree<N, T> inverted_definition() const {
//			powerset_btree<N, T> inverted_definition(this->definition);
//			invert_values(inverted_definition);
//			return inverted_definition;
//		}

		void assign(const subset& set, const T& value) {
			if (set != fullset){
				this->definition.update_or_insert(set, 1/value);
				// The following part ensures that w(fod) is defined (as it is not when directly building w).
				// Indeed, w(fod) may be required for the computation of the commonality function.
				// w(fod)^{-1} is equal to the product of all other weights.
				set_N_value<N, T>* full_set = this->definition[fullset];
				if(full_set){
					full_set->value *= value;
				}else{
					this->definition.insert(fullset, value);
				}
			}else{
				std::cout << "Cannot directly assign the fullset value of a conjunctive decomposition. Ignoring command.\n";
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
			if (set != fullset){
				set_N_value<N, T>* A = this->definition[set];
				if(A){
					set_N_value<N, T>* full_set = this->definition[fullset];
					if(full_set){
						full_set->value *= A->value;
					}else{
						// if there is no fod in definition, this means that its associated value is 1.
						// Thus, as above, nullifying A implies multiplying 1 with A->value,
						// which translates here to inserting A->value at fod.
						this->definition.insert(fullset, A->value);
					}
					this->definition.nullify(A);
				}
			}else{
				std::cout << "Cannot directly assign the fullset value of a conjunctive decomposition. Ignoring command.\n";
			}
		}

		void compute_fullset_value_from_definition(){
			compute_fullset_value_from_definition(this->definition);
		}

		static void compute_fullset_value_from_definition(powerset_btree<N, T>& definition){
			const std::vector<set_N_value<N, T>* >& focal_log_elements = definition.elements();
			T val = 1;
			for (size_t i = 0; i < focal_log_elements.size(); ++i){
				if (focal_log_elements[i]->set != fullset)
					val /= focal_log_elements[i]->value;
			}
			definition.update_or_insert(fullset, val);
		}

		template <class fusion_rule>
		conjunctive_weight<N, T> fuse_with(const conjunctive_weight<N, T>& w2) const {
			const fusion_rule fusion;
			return fusion(*this, w2);
		}

//	protected:
//
//		static inline void invert_values(powerset_btree<N, T>& values) {
//			std::vector<set_N_value<N, T>* > elements = values.elements();
//			for (size_t i = 0; i < elements.size(); ++i){
//				elements[i]->value = 1 / elements[i]->value;
//			}
//		}
	};

} // namespace efficient_DST

#endif // EFFICIENT_DST_CONJUNCTIVE_WEIGHT_HPP
