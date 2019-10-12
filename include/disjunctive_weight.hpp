#ifndef EFFICIENT_DST_DISJUNCTIVE_WEIGHT_HPP
#define EFFICIENT_DST_DISJUNCTIVE_WEIGHT_HPP

#include <decomposition_weight.hpp>
#include <mobius_transform.hpp>
#include <zeta_transform.hpp>

namespace efficient_DST{

	template <typename T = double>
	class disjunctive_weight : public decomposition_weight<T> {
	public:

		disjunctive_weight(const disjunctive_weight<T>& v) : decomposition_weight<T>(v.get_definition())
		{}

		disjunctive_weight(const powerset_btree<T>& focal_log_sets_values) : decomposition_weight<T>(focal_log_sets_values)
		{}

		disjunctive_weight(FOD& fod) : decomposition_weight<T>(fod)
		{}

		disjunctive_weight(const zeta_transform<T>& b) : decomposition_weight<T>(b)
		{
			if (b.order_relation != order_relation_t::subset) {
				std::cerr << "The given Möbius aggregate is not the implicability function and thus can only be the commonality one. "
						<< "\nDoing so, you got the conjunctive weight function instead." << std::endl;
			}
			invert_values(this->definition);
		}


		/*
		 * The disjunctive weight function is the inverse of the multiplicative Möbius transform of the implicability function.
		 * So, invert its values before computing it.
		 */
		powerset_btree<T> inverted_definition() const {
			powerset_btree<T> inverted_definition(this->definition);
			invert_values(inverted_definition);
			return inverted_definition;
		}

		void set_values(const std::unordered_map<std::vector<std::string>, T>& values) {
			for (std::pair<std::vector<std::string>, T> labels_U_value : values){
				set_value(labels_U_value.first, labels_U_value.second);
			}
		}

		void set_value(const std::vector<std::string>& labels, const T& value) {
			this->definition.insert(labels, value);
			// The following part ensures that v(emptyset) is defined (as it is not when directly building v).
			// Indeed, v(emptyset) may be required for the computation of the implicability function.
			// v(emptyset)^{-1} is equal to the product of all other weights.
			set_N_value<T>* emptyset = this->definition.sub_fod_of_size(0);
			if(emptyset){
				emptyset->value /= value;
			}else{
				this->definition.set_value_of_sub_fod_of_size(0, 1/value);
			}
		}

		void set_fod_value(const T& value) {
			set_value(this->definition.get_FOD()->to_set(this->definition.get_FOD()->elements()), value);
		}

		void nullify(const std::vector<std::string>& labels) {
			if(this->definition.get_FOD()->to_set(labels).count() == 0)
				return;
			set_N_value<T>* A = this->definition[labels];
			if(A){
				set_N_value<T>* emptyset = this->definition.sub_fod_of_size(0);
				if(A != emptyset){
					if(emptyset){
						emptyset->value *= A->value;
					}else{
						this->definition.set_value_of_sub_fod_of_size(0, A->value);
					}
				}
				this->definition.nullify(A);
			}
		}

		void compute_emptyset_value_from_definition(){
			compute_emptyset_value_from_definition(this->definition);
		}

		static void compute_emptyset_value_from_definition(powerset_btree<T>& definition){
			boost::dynamic_bitset<> emptyset(definition.get_FOD_size());
			const std::vector<set_N_value<T>* >& focal_log_elements = definition.strict_supersets_of(emptyset);
			T val = 1;
			for (size_t i = 0; i < focal_log_elements.size(); ++i){
				val /= focal_log_elements[i]->value;
			}
			set_N_value<T>* emptyset_set_N_value = definition.sub_fod_of_size(0);
			if(emptyset_set_N_value){
				emptyset_set_N_value->value = val;
			}else{
				definition.set_value_of_sub_fod_of_size(0, val);
			}
		}

		template <class fusion_rule>
		disjunctive_weight<T> apply(const disjunctive_weight<T>& v2) const {
			const fusion_rule fusion;
			return fusion(*this, v2);
		}

	protected:

		static void invert_values(powerset_btree<T>& values) {
			std::vector<set_N_value<T>* > elements = values.elements();
			for (size_t i = 0; i < elements.size(); ++i){
				elements[i]->value = 1 / elements[i]->value;
			}
		}
	};

} // namespace efficient_DST

#endif // EFFICIENT_DST_DISJUNCTIVE_WEIGHT_HPP
