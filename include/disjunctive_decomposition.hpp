#ifndef OW_BFT_DISJUNCTIVE_DECOMPOSITION_HPP
#define OW_BFT_DISJUNCTIVE_DECOMPOSITION_HPP

#include <canonical_decomposition.hpp>
#include <implicability.hpp>
#include <mass.hpp>

namespace ow_bft{

	template <typename T = double>
	class disjunctive_decomposition : public canonical_decomposition<T>{
	protected:

		const implicability<T> implicability_equivalent;
		const FOD *fod;
		/*
		 * special_elements are focal planes.
		 * All sets that are not assigned in this->special_elements have a disjunctive weight equal to 1.
		 */
		powerset_btree<T> special_elements;
		bool is_consonant;
		bool is_quasi_bayesian;

		void display_message_no_implicability_function_provided(){
			std::clog << "No implicability function provided for this implicability function aggregate." << std::endl
					<< "Initializing with degenerate case." << std::endl;
		}

	public:

		disjunctive_decomposition(const mass<T>& m) :
			implicability_equivalent(m),
			fod(&(implicability_equivalent.get_FOD())),
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

		disjunctive_decomposition(const powerset_btree<T>& m_focal_elements) :
			implicability_equivalent(m_focal_elements),
			fod(&(implicability_equivalent.get_FOD())),
			special_elements(*fod, this->block_size),
			is_consonant(false),
			is_quasi_bayesian(false)
		{
			compute_values_for_special_elements();
		}

		disjunctive_decomposition(const implicability<T>& q) :
			implicability_equivalent(q),
			fod(&(implicability_equivalent.get_FOD())),
			special_elements(*fod, this->block_size),
			is_consonant(false),
			is_quasi_bayesian(false)
		{
			compute_values_for_special_elements();
		}

		disjunctive_decomposition(const powerset_btree<T>& m_focal_elements, const powerset_btree<T>& b_special_elements) :
			implicability_equivalent(m_focal_elements, b_special_elements),
			fod(&(implicability_equivalent.get_FOD())),
			special_elements(*fod, this->block_size),
			is_consonant(false),
			is_quasi_bayesian(false)
		{
			compute_values_for_special_elements();
		}

		disjunctive_decomposition(const disjunctive_decomposition<T>& v) :
			implicability_equivalent(v.get_implicability_equivalent()),
			fod(&(implicability_equivalent.get_FOD())),
			special_elements(v.get_special_elements()),
			is_consonant(false),
			is_quasi_bayesian(false)
		{}

		disjunctive_decomposition(const powerset_btree<T>& m_focal_elements, const powerset_btree<T>& b_special_elements, const powerset_btree<T>& _special_elements) :
			implicability_equivalent(m_focal_elements, b_special_elements),
			fod(&(implicability_equivalent.get_FOD())),
			special_elements(_special_elements),
			is_consonant(false),
			is_quasi_bayesian(false)
		{}

		disjunctive_decomposition(const mass_aggregate<T>& ma) : disjunctive_decomposition(ma.get_mass_equivalent()) // @suppress("Class members should be properly initialized")
		{}

		disjunctive_decomposition(const FOD& fod) : disjunctive_decomposition(fod, degenerate)
		{
			display_message_no_implicability_function_provided();
		}

		disjunctive_decomposition(const FOD& _fod, const Special_case s_case) :
			implicability_equivalent(_fod, s_case),
			fod(&(implicability_equivalent.get_FOD())),
			special_elements(*fod, this->block_size),
			is_consonant(false),
			is_quasi_bayesian(false)
		{
			compute_values_for_special_elements();
		}


		template <class fusion_rule>
		disjunctive_decomposition<T> apply(const fusion_rule fusion, const disjunctive_decomposition<T>& v2) const {
			return fusion(*this, v2);
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

		const implicability<T>& get_implicability_equivalent() const {
			return this->implicability_equivalent;
		}

		void compute_values_for_special_elements() {
			if(this->implicability_equivalent.get_mass_equivalent().is_normal()){
				std::cerr << "\nThe disjunctive decomposition is not defined for a normal BBA.";
				exit(EXIT_FAILURE);
			}

			if(&this->implicability_equivalent.get_FOD() != this->special_elements.fod){
				std::cerr << "\nTrying to compute values of a disjunctive decomposition from an implicability function\n"
						  << "with given special elements that refer to FOD elements different from the ones of this implicability function.\n";
				exit(1);
			}

			// Initialize implicabilities on focal planes with the ones on all focal sets
			powerset_btree<T> b_on_focal_planes(this->implicability_equivalent.get_special_elements());

			// Initialize commonalities on negative focal planes with the ones on all negative focal sets
			powerset_btree<T> neg_q_on_neg_focal_planes(*b_on_focal_planes.fod, b_on_focal_planes.block_size);
			powerset_btree<T>::reverse_powerset_from_to(b_on_focal_planes, neg_q_on_neg_focal_planes);

			//powerset_btree<T> neg_q_on_neg_focal_sets(neg_q_on_neg_focal_planes);

			std::unordered_map<size_t, std::vector<set_N_value<T>* > > b_on_focal_planes_map = b_on_focal_planes.elements_by_set_cardinality();
			std::unordered_map<size_t, std::vector<set_N_value<T>* > > neg_q_on_neg_focal_planes_map = neg_q_on_neg_focal_planes.elements_by_set_cardinality();

			std::clog << "\nFocal planes initialized with focal sets\n";

			powerset_btree<T> neg_m(*b_on_focal_planes.fod, b_on_focal_planes.block_size);
			powerset_btree<T> neg_q(*b_on_focal_planes.fod, b_on_focal_planes.block_size);
			powerset_btree<T>::reverse_powerset_from_to(this->implicability_equivalent.get_mass_equivalent().get_focal_elements(), neg_m);
			powerset_btree<T>::reverse_powerset_from_to(this->implicability_equivalent.get_special_elements(), neg_q);

			commonality<T> neg_commonality_equivalent(
					neg_m,
					neg_q
			);

			std::vector<set_N_value<T>* > neg_q_on_pure_neg_focal_planes;

			const bool& all_focal_points_found = canonical_decomposition<T>::linear_analysis_of_focal_points(
				neg_commonality_equivalent,
				neg_q_on_neg_focal_planes, b_on_focal_planes,
				neg_q_on_neg_focal_planes_map, b_on_focal_planes_map,
				neg_q_on_pure_neg_focal_planes,
				this->is_quasi_bayesian
			);

			if(!all_focal_points_found){
				std::unordered_map<size_t, std::vector<set_N_value<T>* > > neg_q_on_neg_focal_sets_map = neg_commonality_equivalent.get_special_elements()
																										.elements_by_set_cardinality();
				const std::vector<std::vector<set_N_value<T>* >* >& neg_q_on_neg_focal_sets_ordered_vector = neg_commonality_equivalent.get_special_elements()
																				.get_vector_of_vectors_ordered_by_cardinality(neg_q_on_neg_focal_sets_map);

				canonical_decomposition<T>::compute_focal_points(
					neg_commonality_equivalent, this->implicability_equivalent.get_special_elements(),
					neg_q_on_neg_focal_planes, b_on_focal_planes,
					neg_q_on_neg_focal_planes_map, b_on_focal_planes_map,
					neg_q_on_neg_focal_sets_ordered_vector,
					neg_q_on_pure_neg_focal_planes
				);
				// Consonance check
				canonical_decomposition<T>::consonance_check(
					neg_q_on_pure_neg_focal_planes,
					neg_q_on_neg_focal_sets_ordered_vector,
					neg_q_on_neg_focal_planes.fod,
					this->is_consonant
				);
			}

			std::clog << "\nFocal planes found: \n";
			print<T>(std::clog, b_on_focal_planes);

			canonical_decomposition<T>::compute_disjunctive_decomposition(
				this->implicability_equivalent.get_special_elements(),
				b_on_focal_planes,
				b_on_focal_planes_map,
				this->special_elements,
				this->is_quasi_bayesian,
				this->is_consonant
			);

			std::clog << "\nWeights on focal planes computed\n";
		}

		static void to_implicability_special_elements(const powerset_btree<T>& v_tree, powerset_btree<T>& b_tree) {
			if(v_tree.fod != b_tree.fod){
				std::cerr << "\nTrying to transform a disjunctive decomposition into an implicability function\n"
						  << "with given special elements that refer to FOD elements different from the ones of this disjunctive decomposition.\n";
				exit(1);
			}

			std::unordered_map<size_t, std::vector<set_N_value<T>* > > elements_by_set_cardinality = v_tree.elements_by_set_cardinality();

			to_implicability_special_elements(v_tree, elements_by_set_cardinality, b_tree);
		}

		static void to_implicability_special_elements(const powerset_btree<T>& v_tree, const powerset_btree<T>& m_tree, powerset_btree<T>& b_tree) {

			/*
			 * v_tree contains focal elements as well as other focal planes. Therefore, a reference to focal elements only helps : m_tree.
			 */

			if(v_tree.fod != m_tree.fod || v_tree.fod != b_tree.fod){
				std::cerr << "\nTrying to transform a disjunctive decomposition into an implicability function\n"
						  << "with given special elements that refer to FOD elements different from the ones of this disjunctive decomposition.\n";
				exit(1);
			}

			std::unordered_map<size_t, std::vector<set_N_value<T>* > > elements_by_set_cardinality = m_tree.elements_by_set_cardinality();

			to_implicability_special_elements(v_tree, elements_by_set_cardinality, b_tree);
		}

		static void to_mass_focal_elements(const powerset_btree<T>& v_tree, powerset_btree<T>& m_tree) {
			if(v_tree.fod != m_tree.fod){
				std::cerr << "\nTrying to transform a conjunctive decomposition into a mass function\n"
						  << "with given focal elements that refer to FOD elements different from the ones of this conjunctive decomposition.\n";
				exit(1);
			}

			powerset_btree<T> b_tree(*(v_tree.fod), v_tree.block_size);
			to_implicability_special_elements(v_tree, b_tree);
			implicability<T>::to_mass_focal_elements(b_tree, m_tree);
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
			const std::vector<fod_element*>& fod_elements = this->fod->to_elements(labels);
			set_N_value<T>* set_value = this->special_elements[fod_elements];
			if(set_value){
				return set_value->value;
			}
			return compute_aggregation(fod_elements);
		}

		T find(const boost::dynamic_bitset<>& key) const {
			set_N_value<T>* set_value = this->special_elements[key];
			if(set_value){
				return set_value->value;
			}
			return compute_aggregation(key);
		}

	protected:

		static void to_implicability_special_elements(const powerset_btree<T>& v_tree, std::unordered_map<size_t, std::vector<set_N_value<T>* > >& f_elements, powerset_btree<T>& b_tree) {

			const std::vector<set_N_value<T>* >& v_elements = v_tree.elements();

			std::clog << "\nFocal planes : \n";
			print<T>(std::clog, v_tree);

			T v_product = 1;
			for (size_t i = 0; i < v_elements.size(); ++i) {
				v_product *= v_elements[i]->value;
			}

			// set emptyset to v_product
			b_tree.set_value_of_sub_fod_of_size(0, v_product);

			const std::vector<std::vector<set_N_value<T>* >* >& ordered_vector = v_tree
								.get_vector_of_vectors_ordered_by_cardinality(f_elements);

			for (size_t c = 0; c < ordered_vector.size(); ++c) {
				for (size_t i = 0; i < ordered_vector[c]->size(); ++i) {

					const std::vector<set_N_value<T>* >& v_i_subsets = v_tree.subsets_of((*ordered_vector[c])[i]->fod_elements);
					T b_i_val = v_product;

					for (size_t ii = 0; ii < v_i_subsets.size(); ++ii) {
						b_i_val /= v_i_subsets[ii]->value;
					}

					b_tree.insert((*ordered_vector[c])[i]->fod_elements, b_i_val);
				}
			}
		}

		T compute_aggregation_at_emptyset() const {
			std::cerr << "The disjunctive weight function is not defined at emptyset";
			return -1;
		}

		T compute_aggregation_at_fod() const {
			return 1;
		}

		T compute_aggregation(const boost::dynamic_bitset<>& key) const {
			if(key.count() == this->fod->size())
				// conjunctive weights are only defined for strict subsets of FOD
				return compute_aggregation_at_fod();
			return 1;
		}

		T compute_aggregation(const std::vector<fod_element*>& fod_elements) const {
			if(fod_elements.size() == 0)
				// disjunctive weights are only defined for strict supersets of emptyset
				return compute_aggregation_at_emptyset();
			return 1;
		}
/*
		static bool linear_analysis_of_focal_planes(
				const implicability<T>& implicability_equivalent,
				powerset_btree<T>& b_on_focal_planes,
				powerset_btree<T>& neg_q_on_neg_focal_planes,
				std::unordered_map<size_t, std::vector<set_N_value<T>* > >& b_on_focal_planes_map,
				std::unordered_map<size_t, std::vector<set_N_value<T>* > >& neg_q_on_neg_focal_planes_map,
				std::vector<set_N_value<T>* >& b_on_pure_focal_planes
				){

			bool generated_quasi_bayesian = false;
			std::vector<set_N_value<T>* > focal_sets_except_emptyset;
			focal_sets_except_emptyset.reserve(b_on_focal_planes.fod->size()-1);

			for(auto kv : b_on_focal_planes_map) {
			    if(kv.first > 0){
					for (size_t i = 0; i < kv.second.size(); ++i) {
						focal_sets_except_emptyset.emplace_back(kv.second[i]);
					}
			    }
			}

			boost::dynamic_bitset<> neg_U = b_on_focal_planes.fod->to_set(focal_sets_except_emptyset[0]->fod_elements);
			const boost::dynamic_bitset<> emptyset(b_on_focal_planes.fod->size());
			std::clog << "\nneg_U = "<< focal_sets_except_emptyset[0]->fod_elements;

			for (size_t i = 1; i < focal_sets_except_emptyset.size(); ++i) {
				std::clog << "\nneg_I = union(neg_U, "<< focal_sets_except_emptyset[i]->fod_elements << ")\n";
				const boost::dynamic_bitset<>& A = b_on_focal_planes.fod->to_set(focal_sets_except_emptyset[i]->fod_elements);
				const boost::dynamic_bitset<>& neg_I = b_on_focal_planes.fod->set_union(neg_U, A);
				const size_t& I_card = b_on_focal_planes.fod->size() - neg_I.count();
				std::clog << "|I| = " << I_card << std::endl;
				if(I_card <= 1){
					// add it to b_on_focal_planes if it wasn't already there
					if(!b_on_focal_planes[neg_I]){
						set_N_value<T>* inserted_focal_point = insert_focal_point(
							neg_I,
							commonality_equivalent,
							b_on_focal_planes,
							b_on_focal_planes_map
						);
						insert_neg_focal_point(
							b_on_focal_planes.fod->set_negate(neg_I),
							inserted_focal_point->value,
							neg_b_on_neg_focal_points,
							neg_b_on_neg_focal_points_map
						);
						b_on_pure_focal_planes.push_back(inserted_focal_point);
					}

					if(I_card == 1){
						generated_quasi_bayesian = true;
					}
				}else{
					std::clog << "=> Linear analysis aborted.\n";
					return false;
				}
				neg_U = b_on_focal_planes.fod->set_union(
					neg_U,
					A
				);
			}

			if(generated_quasi_bayesian){
				// add also emptyset to q_on_focal_points if it wasn't already there
				if(!b_on_focal_planes[emptyset]){
					set_N_value<T>* inserted_focal_point = insert_focal_point(
						emptyset,
						commonality_equivalent,
						q_on_focal_points,
						q_on_focal_points_map
					);
					insert_neg_focal_point(
						q_on_focal_points.fod->set_negate(emptyset),
						inserted_focal_point->value,
						neg_b_on_neg_focal_points,
						neg_b_on_neg_focal_points_map
					);
					q_on_pure_focal_points.push_back(inserted_focal_point);
				}
			}
			std::clog << "=> Successful linear analysis of focal points.\n";
			return true;
		}
		*/
	};

} // namespace ow_bft

#endif // OW_BFT_DISJUNCTIVE_DECOMPOSITION_HPP
