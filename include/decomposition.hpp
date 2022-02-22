#ifndef EFFICIENT_DST_DECOMPOSITION_HPP
#define EFFICIENT_DST_DECOMPOSITION_HPP

#include <weight_function.hpp>
#include <zeta_transform.hpp>

namespace efficient_DST{

	template <class inclusion, size_t N, typename T = float>
	class decomposition : public weight_function<N, T> {
	public:
		using typename powerset_function<N, T>::subset;
		using powerset_function<N, T>::emptyset;
		using powerset_function<N, T>::fullset;

		set_N_value<N, T> normalizing_set_assignment;


		decomposition(const weight_function<N, T>& w) : weight_function<N, T>(w.outcomes, w.definition)
		{
			compute_normalizing_set_assignment();
		}

		decomposition(const sample_space<N>& outcomes) : weight_function<N, T>(outcomes)
		{
			compute_normalizing_set_assignment();
		}

		decomposition(
			const zeta_transform<inclusion, N, T>& q
		) : weight_function<N, T>(q)
		{
			compute_normalizing_set_assignment();
		}


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
			return 1 - weight_function<N, T>::operator [](set);
		}

		std::ostream& print(const bool& including_null = false) const {
			std::vector<set_N_value<N, T>* > values = this->definition.elements(including_null);
			std::cout << std::endl;
			for (size_t i = 0; i < values.size(); ++i) {
				if(values[i]->set != this->normalizing_set_assignment.set){
					if (values[i]->is_null)
						std::cout << "NULL\t <- " + this->outcomes.to_string(values[i]->set) << std::endl;
					else
						std::cout << set_N_value<N, T>::to_string(1-1/values[i]->value) + "\t <- " + this->outcomes.to_string(values[i]->set) << std::endl;
				}
			}
//			this->definition.print_layout();
			return std::cout;
		}

		void assign(const subset& set, const T& mass) {
//			if (set != fullset){
			T weight = 1-mass;
			this->definition.update_or_insert(set, 1/weight);
			this->normalizing_set_assignment.set = inclusion::set_dual_operation(this->normalizing_set_assignment.set, set);
			this->normalizing_set_assignment.value *= weight;
//			}else{
//				std::cout << "Cannot directly assign the fullset value of a conjunctive decomposition. Ignoring command.\n";
//			}
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

		void nullify(const std::vector<std::string>& labels) {
			this->nullify(this->outcomes.get_subset(labels));
		}

		void nullify(const subset& set) {
//			if (set != fullset){
			set_N_value<N, T>* A = this->definition[set];
			if(A){
				this->definition.nullify(A);
				this->compute_normalizing_set_assignment();
			}
//			}else{
//				std::cout << "Cannot directly assign the fullset value of a conjunctive decomposition. Ignoring command.\n";
//			}
		}

//		set_N_value<N, T> compute_encompassing_set_assignment(){
//			compute_encompassing_set_assignment(this->definition);
//		}

		void compute_normalizing_set_assignment(){
			this->normalizing_set_assignment = set_N_value<N, T>(emptyset, 1);
			const std::vector<set_N_value<N, T>* >& focal_log_elements = this->definition.elements();
			for (size_t i = 0; i < focal_log_elements.size(); ++i){
				this->normalizing_set_assignment.set = inclusion::set_dual_operation(this->normalizing_set_assignment.set, focal_log_elements[i]->set);
			}
			for (size_t i = 0; i < focal_log_elements.size(); ++i){
				if (focal_log_elements[i]->set != this->normalizing_set_assignment.set)
					this->normalizing_set_assignment.value /= focal_log_elements[i]->value;
			}
		}
	};

} // namespace efficient_DST

#endif // EFFICIENT_DST_DECOMPOSITION_HPP
