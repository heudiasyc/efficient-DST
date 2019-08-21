#ifndef EFFICIENT_DST_CONJUNCTIVE_WEIGHT_HPP
#define EFFICIENT_DST_CONJUNCTIVE_WEIGHT_HPP

#include <commonality.hpp>
#include <decomposition_weight.hpp>
#include <mass.hpp>

namespace efficient_DST{

	template <typename T = double>
	class conjunctive_weight : public decomposition_weight<T>{
	protected:

		const commonality<T> commonality_equivalent;
		const FOD *fod;
		/*
		 * special_elements are focal points.
		 * All sets that are not assigned in this->special_elements have a conjunctive weight equal to 1.
		 */
		powerset_btree<T> special_elements;
		bool is_consonant;
		bool is_quasi_bayesian;

		void display_message_no_commonality_function_provided(){
			std::clog << "No commonality function provided for this commonality function aggregate." << std::endl
					<< "Initializing with vacuous case." << std::endl;
		}

	public:

		conjunctive_weight(const mass<T>& m) :
			commonality_equivalent(m),
			fod(&(commonality_equivalent.get_FOD())),
			special_elements(*fod, this->block_size),
			is_consonant(false),
			is_quasi_bayesian(false)
		{
			if(!this->mass_equivalent.is_valid()){
				std::cerr << "\nError: The sum of all images of the mass function provided to this mass aggregate isn't equal to 1.";
				exit(1);
			}
			compute_values_for_special_elements();
		}

		conjunctive_weight(const powerset_btree<T>& m_focal_elements) :
			commonality_equivalent(m_focal_elements),
			fod(&(commonality_equivalent.get_FOD())),
			special_elements(*fod, this->block_size),
			is_consonant(false),
			is_quasi_bayesian(false)
		{
			compute_values_for_special_elements();
		}

		conjunctive_weight(const commonality<T>& q) :
			commonality_equivalent(q),
			fod(&(commonality_equivalent.get_FOD())),
			special_elements(*fod, this->block_size),
			is_consonant(false),
			is_quasi_bayesian(false)
		{
			compute_values_for_special_elements();
		}

		conjunctive_weight(const powerset_btree<T>& m_focal_elements, const powerset_btree<T>& q_special_elements) :
			commonality_equivalent(m_focal_elements, q_special_elements),
			fod(&(commonality_equivalent.get_FOD())),
			special_elements(*fod, this->block_size),
			is_consonant(false),
			is_quasi_bayesian(false)
		{
			compute_values_for_special_elements();
		}

		conjunctive_weight(const conjunctive_weight<T>& w) :
			commonality_equivalent(w.get_commonality_equivalent()),
			fod(&(commonality_equivalent.get_FOD())),
			special_elements(w.get_special_elements()),
			is_consonant(false),
			is_quasi_bayesian(false)
		{}

		conjunctive_weight(const powerset_btree<T>& m_focal_elements, const powerset_btree<T>& q_special_elements, const powerset_btree<T>& _special_elements) :
			commonality_equivalent(m_focal_elements, q_special_elements),
			fod(&(commonality_equivalent.get_FOD())),
			special_elements(_special_elements),
			is_consonant(false),
			is_quasi_bayesian(false)
		{}

		conjunctive_weight(const mobius_aggregate<T>& ma) : conjunctive_weight(ma.get_mass_equivalent()) // @suppress("Class members should be properly initialized")
		{}

		conjunctive_weight(const FOD& fod) : conjunctive_weight(fod, vacuous)
		{
			display_message_no_commonality_function_provided();
		}

		conjunctive_weight(const FOD& _fod, const Special_case s_case) :
			commonality_equivalent(_fod, s_case),
			fod(&(commonality_equivalent.get_FOD())),
			special_elements(*fod, this->block_size),
			is_consonant(false),
			is_quasi_bayesian(false)
		{
			compute_values_for_special_elements();
		}


		template <class fusion_rule>
		conjunctive_weight<T> apply(const fusion_rule fusion, const conjunctive_weight<T>& w2) const {
			return fusion(*this, w2);
		}

		const FOD& get_FOD() const {
			return *(this->fod);
		}

		const powerset_btree<T>& get_special_elements() const {
			return this->special_elements;
		}

		void set_special_elements(const powerset_btree<T>& f_elements){
			/*if(&(f_elements.get_FOD()) != this->fod){
				std::cerr << "\nCan't replace special elements with elements using another FOD.\n";
				return;
			}*/
			this->special_elements.nullify();
			this->special_elements.copy(f_elements);
		}

		const commonality<T>& get_commonality_equivalent() const {
			return this->commonality_equivalent;
		}
/*
		void set_commonality_equivalent(const commonality<T>& q, bool update_focal_elements=true){
			this->commonality_equivalent = commonality<T>(q);

			if(update_focal_elements){
				this->fod->erase_powerset(this->focal_elements);
				this->fod = &(q.get_FOD());
				this->fod->push_back_powerset(this->focal_elements);
				this->focal_elements.nullify();
				set_values_for_special_elements();
			}
		}
*/

		void compute_values_for_special_elements() {
			if(this->commonality_equivalent.get_mass_equivalent().is_dogmatic()){
				std::cerr << "\nThe conjunctive decomposition is not defined for a dogmatic BBA.";
				exit(EXIT_FAILURE);
			}

			if(&this->commonality_equivalent.get_FOD() != this->special_elements.fod){
				std::cerr << "\nTrying to compute values of a conjunctive decomposition from a commonality function\n"
						  << "with given special elements that refer to FOD elements different from the ones of this commonality function.\n";
				exit(1);
			}

			// Initialize commonalities on focal points with the ones on all focal sets
			powerset_btree<T> q_on_focal_points(this->commonality_equivalent.get_special_elements());

			// Initialize implicabilities on negative focal points with the ones on all negative focal sets
			powerset_btree<T> neg_b_on_neg_focal_points(*q_on_focal_points.fod, q_on_focal_points.block_size);
			powerset_btree<T>::reverse_powerset_from_to(q_on_focal_points, neg_b_on_neg_focal_points);

			powerset_btree<T> neg_b_on_neg_focal_sets(neg_b_on_neg_focal_points);

			std::unordered_map<size_t, std::vector<set_N_value<T>* > > q_on_focal_points_map = q_on_focal_points.elements_by_set_cardinality();
			std::unordered_map<size_t, std::vector<set_N_value<T>* > > neg_b_on_neg_focal_points_map = neg_b_on_neg_focal_points.elements_by_set_cardinality();

			std::clog << "\nFocal points initialized with focal sets\n";

			std::vector<set_N_value<T>* > q_on_pure_focal_points;

			const bool& all_focal_points_found = decomposition_weight<T>::linear_analysis_of_focal_points(
				this->commonality_equivalent,
				q_on_focal_points, neg_b_on_neg_focal_points,
				q_on_focal_points_map, neg_b_on_neg_focal_points_map,
				q_on_pure_focal_points,
				this->is_quasi_bayesian
			);

			if(!all_focal_points_found){
				std::unordered_map<size_t, std::vector<set_N_value<T>* > > q_on_focal_sets_map = commonality_equivalent.get_special_elements()
																										.elements_by_set_cardinality();
				const std::vector<std::vector<set_N_value<T>* >* >& q_on_focal_sets_ordered_vector = commonality_equivalent.get_special_elements()
																				.get_vector_of_vectors_ordered_by_cardinality(q_on_focal_sets_map);

				decomposition_weight<T>::compute_focal_points(
					this->commonality_equivalent, neg_b_on_neg_focal_sets,
					q_on_focal_points, neg_b_on_neg_focal_points,
					q_on_focal_points_map, neg_b_on_neg_focal_points_map,
					q_on_focal_sets_ordered_vector,
					q_on_pure_focal_points
				);
				// Consonance check
				decomposition_weight<T>::consonance_check(
					q_on_pure_focal_points,
					q_on_focal_sets_ordered_vector,
					q_on_focal_points.fod,
					this->is_consonant
				);
			}

			std::clog << "\nFocal points found: \n";
			print<T>(std::clog, q_on_focal_points);

			powerset_btree<T> neg_v_special_elements(*q_on_focal_points.fod, q_on_focal_points.block_size);

			decomposition_weight<T>::compute_disjunctive_decomposition(
				neg_b_on_neg_focal_sets,
				neg_b_on_neg_focal_points,
				neg_b_on_neg_focal_points_map,
				neg_v_special_elements,
				this->is_quasi_bayesian,
				this->is_consonant
			);

			powerset_btree<T>::reverse_powerset_from_to(neg_v_special_elements, this->special_elements);

			std::clog << "\nWeights on focal points computed\n";
		}

		static void to_commonality_special_elements(const powerset_btree<T>& w_tree, powerset_btree<T>& q_tree) {
			if(w_tree.fod != q_tree.fod){
				std::cerr << "\nTrying to transform a conjunctive decomposition into a commonality function\n"
						  << "with given special elements that refer to FOD elements different from the ones of this conjunctive decomposition.\n";
				exit(1);
			}

			std::unordered_map<size_t, std::vector<set_N_value<T>* > > elements_by_set_cardinality = w_tree.elements_by_set_cardinality();

			to_commonality_special_elements(w_tree, elements_by_set_cardinality, q_tree);
		}

		static void to_commonality_special_elements(const powerset_btree<T>& w_tree, const powerset_btree<T>& m_tree, powerset_btree<T>& q_tree) {

			/*
			 * w_tree contains focal elements as well as other focal points. Therefore, a reference to focal elements only helps : m_tree.
			 */

			if(w_tree.fod != m_tree.fod || w_tree.fod != q_tree.fod){
				std::cerr << "\nTrying to transform a conjunctive decomposition into a commonality function\n"
						  << "with given special elements that refer to FOD elements different from the ones of this conjunctive decomposition.\n";
				exit(1);
			}

			std::unordered_map<size_t, std::vector<set_N_value<T>* > > elements_by_set_cardinality = m_tree.elements_by_set_cardinality();

			to_commonality_special_elements(w_tree, elements_by_set_cardinality, q_tree);
		}

		static void to_mass_focal_elements(const powerset_btree<T>& w_tree, powerset_btree<T>& m_tree) {
			if(w_tree.fod != m_tree.fod){
				std::cerr << "\nTrying to transform a conjunctive decomposition into a mass function\n"
						  << "with given focal elements that refer to FOD elements different from the ones of this conjunctive decomposition.\n";
				exit(1);
			}

			powerset_btree<T> q_tree(*(w_tree.fod), w_tree.block_size);
			to_commonality_special_elements(w_tree, q_tree);
			commonality<T>::to_mass_focal_elements(q_tree, m_tree);
		}

		T at_emptyset() const {
			set_N_value<T>* set_value = this->special_elements.sub_fod_of_size(0);
			if(set_value){
				return set_value->value;
			}
			return compute_aggregation_at_emptyset();
		}

		T at_fod() const {
			set_N_value<T>* set_value = this->special_elements.sub_fod_of_size(this->fod->size());
			if(set_value){
				return set_value->value;
			}
			return compute_aggregation_at_fod();
		}

		T operator[](const std::vector<std::string>& labels) const {
			return find(this->fod->to_set(labels));
		}

		T find(const boost::dynamic_bitset<>& set) const {
			set_N_value<T>* set_value = this->special_elements[set];
			if(set_value){
				return set_value->value;
			}
			return compute_aggregation(set);
		}

	protected:

		static void to_commonality_special_elements(const powerset_btree<T>& w_tree, std::unordered_map<size_t, std::vector<set_N_value<T>* > >& f_elements, powerset_btree<T>& q_tree) {

			const std::vector<set_N_value<T>* >& w_elements = w_tree.elements();

			std::clog << "\nFocal points : \n";
			print<T>(std::clog, w_tree);

			T w_product = 1;
			for (size_t i = 0; i < w_elements.size(); ++i) {
				w_product *= w_elements[i]->value;
			}

			// set FOD to w_product
			q_tree.set_value_of_sub_fod_of_size(q_tree.fod->size(), w_product);

			const std::vector<std::vector<set_N_value<T>* >* >& ordered_vector = w_tree
								.get_vector_of_vectors_ordered_by_cardinality(f_elements);

			for (size_t c = 0; c < ordered_vector.size(); ++c) {
				for (size_t i = 0; i < ordered_vector[c]->size(); ++i) {

					const std::vector<set_N_value<T>* >& w_i_supersets = w_tree.supersets_of((*ordered_vector[c])[i]->set);
					T q_i_val = w_product;

					for (size_t ii = 0; ii < w_i_supersets.size(); ++ii) {
						q_i_val /= w_i_supersets[ii]->value;
					}

					q_tree.insert((*ordered_vector[c])[i]->set, q_i_val);
				}
			}
		}

		T compute_aggregation_at_emptyset() const {
			return 1;
		}

		T compute_aggregation_at_fod() const {
			std::cerr << "The conjunctive weight function is not defined at FOD";
			return -1;
		}

		T compute_aggregation(const boost::dynamic_bitset<>& set) const {
			if(set.count() == this->fod->size())
				// conjunctive weights are only defined for strict subsets of FOD
				return compute_aggregation_at_fod();
			return 1;
		}
	};

} // namespace efficient_DST

#endif // EFFICIENT_DST_CONJUNCTIVE_WEIGHT_HPP
