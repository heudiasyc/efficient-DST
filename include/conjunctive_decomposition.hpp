#ifndef OW_BFT_CONJUNCTIVE_DECOMPOSITION_HPP
#define OW_BFT_CONJUNCTIVE_DECOMPOSITION_HPP

#include <aggregate.hpp>
#include <commonality.hpp>
#include <mass.hpp>

namespace ow_bft{

	/*
	 * Here, this->special_elements won't only designate focal elements but all focal points (see compute_focal_points)
	 * that have a weight different from 1.
	 * All set that is not assigned in this->special_elements has a conjunctive weight equal to 1.
	 */
	template <typename T = double>
	class conjunctive_decomposition : public aggregate<T>{
	protected:

		const commonality<T> commonality_equivalent;

		const FOD *fod;
		/*
		 * special_elements are the only elements (and images through this aggregate) needed to allow for the reconstruction
		 * of the mass images of focal elements (could be only focal elements).
		 * Here, this->special_elements won't only designate focal elements but all focal points (see compute_focal_points)
		 * that have a weight different from 1.
		 * All sets that are not assigned in this->special_elements have a conjunctive weight equal to 1.
		 */
		powerset_btree<T> special_elements;

		void display_message_no_commonality_function_provided(){
			std::clog << "No commonality function provided for this commonality function aggregate." << std::endl
					<< "Initializing with vacuous case." << std::endl;
		}

	public:

		conjunctive_decomposition(const mass<T>& m) :
			commonality_equivalent(m),
			fod(&(commonality_equivalent.get_FOD())),
			special_elements(*fod, this->block_size)
		{
			if(!this->mass_equivalent.is_valid()){
				std::cerr << "\nError: The sum of all images of the mass function provided to this mass aggregate isn't equal to 1.";
				exit(1);
			}
			compute_values_for_special_elements(this->commonality_equivalent, this->special_elements);
		}

		conjunctive_decomposition(const powerset_btree<T>& m_focal_elements) :
			commonality_equivalent(m_focal_elements),
			fod(&(commonality_equivalent.get_FOD())),
			special_elements(*fod, this->block_size)
		{
			compute_values_for_special_elements(this->commonality_equivalent, this->special_elements);
		}

		conjunctive_decomposition(const commonality<T>& q) :
			commonality_equivalent(q),
			fod(&(commonality_equivalent.get_FOD())),
			special_elements(*fod, this->block_size)
		{
			compute_values_for_special_elements(this->commonality_equivalent, this->special_elements);
		}

		conjunctive_decomposition(const powerset_btree<T>& m_focal_elements, const powerset_btree<T>& q_special_elements) :
			commonality_equivalent(m_focal_elements, q_special_elements),
			fod(&(commonality_equivalent.get_FOD())),
			special_elements(*fod, this->block_size)
		{
			compute_values_for_special_elements(this->commonality_equivalent, this->special_elements);
		}

		conjunctive_decomposition(const conjunctive_decomposition<T>& w) :
			commonality_equivalent(w.get_commonality_equivalent()),
			fod(&(commonality_equivalent.get_FOD())),
			special_elements(w.get_special_elements())
		{}

		conjunctive_decomposition(const powerset_btree<T>& m_focal_elements, const powerset_btree<T>& q_special_elements, const powerset_btree<T>& _special_elements) :
			commonality_equivalent(m_focal_elements, q_special_elements),
			fod(&(commonality_equivalent.get_FOD())),
			special_elements(_special_elements)
		{}

		conjunctive_decomposition(const mass_aggregate<T>& ma) : conjunctive_decomposition(ma.get_mass_equivalent()) // @suppress("Class members should be properly initialized")
		{}

		conjunctive_decomposition(const FOD& fod) : conjunctive_decomposition(fod, vacuous)
		{
			display_message_no_commonality_function_provided();
		}

		conjunctive_decomposition(const FOD& _fod, const Special_case s_case) :
			commonality_equivalent(_fod, s_case),
			fod(&(commonality_equivalent.get_FOD())),
			special_elements(*fod, this->block_size)
		{
			compute_values_for_special_elements(this->commonality_equivalent, this->special_elements);
		}


		template <class fusion_rule>
		conjunctive_decomposition<T> apply(const fusion_rule fusion, const conjunctive_decomposition<T>& w2) const {
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

		static void compute_values_for_special_elements(const commonality<T>& commonality_equivalent, powerset_btree<T>& w_special_elements) {
			if(commonality_equivalent.get_mass_equivalent().is_dogmatic()){
				std::cerr << "\nThe conjunctive decomposition is not defined for a dogmatic BBA.";
				exit(EXIT_FAILURE);
			}

			if(&commonality_equivalent.get_FOD() != w_special_elements.fod){
				std::cerr << "\nTrying to compute values of a conjunctive decomposition from a commonality function\n"
						  << "with given special elements that refer to FOD elements different from the ones of this commonality function.\n";
				exit(1);
			}

			// Initialize commonalities on focal points with the ones on all focal sets
			powerset_btree<T> q_on_focal_points(commonality_equivalent.get_special_elements());

			// Initialize implicabilities on negative focal points with the ones on all negative focal sets
			powerset_btree<T> neg_b_on_neg_focal_points(*q_on_focal_points.fod, q_on_focal_points.block_size);
			powerset_btree<T>::reverse_powerset_from_to(q_on_focal_points, neg_b_on_neg_focal_points);

			powerset_btree<T> neg_b_on_neg_focal_sets(neg_b_on_neg_focal_points);

			std::unordered_map<size_t, std::vector<set_N_value<T>* > > q_on_focal_points_map = q_on_focal_points.elements_by_set_cardinality();
			std::unordered_map<size_t, std::vector<set_N_value<T>* > > neg_b_on_neg_focal_points_map = neg_b_on_neg_focal_points.elements_by_set_cardinality();

			std::clog << "\nFocal points initialized with focal sets\n";

			compute_focal_points(
				commonality_equivalent, neg_b_on_neg_focal_sets,
				q_on_focal_points, neg_b_on_neg_focal_points,
				q_on_focal_points_map, neg_b_on_neg_focal_points_map
			);

			std::clog << "\nFocal points found: \n";
			print<T>(std::clog, q_on_focal_points);

			powerset_btree<T> neg_v_special_elements(*q_on_focal_points.fod, q_on_focal_points.block_size);

			compute_neg_v_special_elements(
				neg_b_on_neg_focal_points,
				neg_b_on_neg_focal_points_map,
				neg_v_special_elements,
				false
			);

			powerset_btree<T>::reverse_powerset_from_to(neg_v_special_elements, w_special_elements);

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

					const std::vector<set_N_value<T>* >& w_i_supersets = w_tree.supersets_of((*ordered_vector[c])[i]->fod_elements);
					T q_i_val = w_product;

					for (size_t ii = 0; ii < w_i_supersets.size(); ++ii) {
						q_i_val /= w_i_supersets[ii]->value;
					}

					q_tree.insert((*ordered_vector[c])[i]->fod_elements, q_i_val);
				}
			}

			/*const std::vector<set_N_value<T>* >& elements = q_tree.elements();
			T val;
			for (size_t i = 0; i < elements.size(); ++i) {
				val = elements[i]->value;
				if(this->is_equivalent_to_zero(val)){
					q_tree.nullify(elements[i]);
				}
			}*/

			//commonality<T> q(w_tree.get_FOD());
			//q.set_special_elements(q_tree);
		}

		T compute_aggregation_at_emptyset() const {
			return 1;
		}

		T compute_aggregation_at_fod() const {
			std::cerr << "The conjunctive weight function is not defined at FOD";
			return -1;
		}

		T compute_aggregation(const boost::dynamic_bitset<>& key) const {
			if(key.count() == this->fod->size())
				// conjunctive weights are only defined for strict subsets of FOD
				return compute_aggregation_at_fod();
			return 1;
		}

		T compute_aggregation(const std::vector<fod_element*>& fod_elements) const {
			if(fod_elements.size() == this->fod->size())
				// conjunctive weights are only defined for strict subsets of FOD
				return compute_aggregation_at_fod();
			return 1;
		}

		static set_N_value<T>* insert_neg_focal_point(
				const boost::dynamic_bitset<>& neg_focal_point,
				T value,
				powerset_btree<T>& neg_b_on_neg_focal_points,
				std::unordered_map<size_t, std::vector<set_N_value<T>* > >& neg_b_on_neg_focal_points_map
				) {

			std::clog << "\nNEW INSERTED NEGATIVE FOCAL POINT :\n";
			std::clog << "\nnegative focal point : " << neg_focal_point << std::endl;

			set_N_value<T>* inserted_neg_focal_point = neg_b_on_neg_focal_points.insert(
					neg_focal_point,
					value
			);
			std::clog << inserted_neg_focal_point->fod_elements << std::endl;

			size_t c = inserted_neg_focal_point->fod_elements.size();

			neg_b_on_neg_focal_points_map[c].push_back(
					inserted_neg_focal_point
			);
			return inserted_neg_focal_point;
		}

		static set_N_value<T>* insert_focal_point(
				const boost::dynamic_bitset<>& focal_point,
				const commonality<T>& commonality_equivalent,
				powerset_btree<T>& q_on_focal_points,
				std::unordered_map<size_t, std::vector<set_N_value<T>* > >& q_on_focal_points_map
				) {

			std::clog << "\nNEW INSERTED FOCAL POINT :\n";
			std::clog << "\nfocal point : " << focal_point << std::endl;

			set_N_value<T>* inserted_focal_point = q_on_focal_points.insert(
					focal_point,
					commonality_equivalent.find(focal_point)
			);
			std::clog << inserted_focal_point->fod_elements << std::endl;

			q_on_focal_points_map[inserted_focal_point->fod_elements.size()].push_back(
					inserted_focal_point
			);
			return inserted_focal_point;
		}

		/*
		 * I use here a notion which I called "focal point" which is the intersection of two focal sets.
		 * That focal point is the biggest set influenced by both these two focal sets in the commonality space.
		 * A focal point is unique among its supersets.
		 * A focal set is also a focal point as it results from the intersection of itself with the FOD element
		 * (which is always a focal element in the context of conjunctive weights)
		 *
		 * Compute all focal points.
		 * If F is the number of focal elements, there are at least F*(F-1)/2 operations here => O(F^2)
		 * In the worst case, given the set of focal elements, it's O(2^(2*N)), where N is the FOD size.
		 * The worst case is obtained if all sets of cardinality this->fod->size()-1 are focal sets
		 */
		static void compute_focal_points(
				const commonality<T>& commonality_equivalent,
				powerset_btree<T>& neg_b_on_neg_focal_sets,
				powerset_btree<T>& q_on_focal_points,
				powerset_btree<T>& neg_b_on_neg_focal_points,
				std::unordered_map<size_t, std::vector<set_N_value<T>* > >& q_on_focal_points_map,
				std::unordered_map<size_t, std::vector<set_N_value<T>* > >& neg_b_on_neg_focal_points_map
				) {

			std::unordered_map<size_t, std::vector<set_N_value<T>* > > q_on_focal_sets_map = q_on_focal_points_map;
			const std::vector<std::vector<set_N_value<T>* >* >& q_on_focal_sets_ordered_vector = q_on_focal_points
																			.get_vector_of_vectors_ordered_by_cardinality(q_on_focal_sets_map);
			std::vector<set_N_value<T>* > q_on_pure_focal_points;

			// for each cardinality of focal set, from the smallest to the biggest, except FOD
			for (size_t c = 0; c < q_on_focal_sets_ordered_vector.size()-1; ++c) {
				for (size_t i = 0; i < q_on_focal_sets_ordered_vector[c]->size(); ++i) {

					const boost::dynamic_bitset<>& setA = q_on_focal_points.fod->to_set((*q_on_focal_sets_ordered_vector[c])[i]->fod_elements);

					// search for sets of same cardinality
					for (size_t ii = i+1; ii < q_on_focal_sets_ordered_vector[c]->size(); ++ii) {

						const boost::dynamic_bitset<>& setB = q_on_focal_points.fod->to_set((*q_on_focal_sets_ordered_vector[c])[ii]->fod_elements);

						// compute their focal point
						const boost::dynamic_bitset<>& focal_point = q_on_focal_points.fod->set_intersection(setA, setB);

						// add it to q_on_focal_points if it wasn't already there
						if(!q_on_focal_points[focal_point]){
							set_N_value<T>* inserted_focal_point = insert_focal_point(
								focal_point,
								commonality_equivalent,
								q_on_focal_points,
								q_on_focal_points_map
							);
							insert_neg_focal_point(
								q_on_focal_points.fod->set_negate(focal_point),
								inserted_focal_point->value,
								neg_b_on_neg_focal_points,
								neg_b_on_neg_focal_points_map
							);

							q_on_pure_focal_points.push_back(inserted_focal_point);
						}
					}

					const std::vector<boost::dynamic_bitset<> >& unions_of_not_subsets_of_smaller_than = neg_b_on_neg_focal_sets
													.unions_of_not_subsets_of_smaller_than(
														q_on_focal_points.fod->set_negate(setA),
														q_on_focal_points.fod->size() - (*q_on_focal_sets_ordered_vector[c])[i]->fod_elements.size()-1
													);
					// search for sets of higher cardinality that are not supersets of the current set
					for (size_t ii = 0; ii < unions_of_not_subsets_of_smaller_than.size(); ++ii) {

						const boost::dynamic_bitset<>& neg_focal_point = unions_of_not_subsets_of_smaller_than[ii];

						// add it to q_on_focal_points if it wasn't already there
						if(!neg_b_on_neg_focal_points[neg_focal_point]){
							const boost::dynamic_bitset<>& focal_point = q_on_focal_points.fod->set_negate(unions_of_not_subsets_of_smaller_than[ii]);

							set_N_value<T>* inserted_focal_point = insert_focal_point(
								focal_point,
								commonality_equivalent,
								q_on_focal_points,
								q_on_focal_points_map
							);
							insert_neg_focal_point(
								q_on_focal_points.fod->set_negate(focal_point),
								inserted_focal_point->value,
								neg_b_on_neg_focal_points,
								neg_b_on_neg_focal_points_map
							);

							q_on_pure_focal_points.push_back(inserted_focal_point);
						}
					}
				}
			}

			// for each pure focal point (focal point that is not a focal set)
			for (size_t i = 0; i < q_on_pure_focal_points.size(); ++i) {

				const boost::dynamic_bitset<>& setA = q_on_focal_points.fod->to_set(q_on_pure_focal_points[i]->fod_elements);

				const std::vector<boost::dynamic_bitset<> >& intersections_of_not_subsets_of_smaller_than = commonality_equivalent
												.get_special_elements()
												.intersections_of_not_subsets_of_smaller_than(
													setA,
													q_on_pure_focal_points[i]->fod_elements.size()-1
												);
				// search for sets of lower cardinality that are not subsets of the current set
				for (size_t ii = 0; ii < intersections_of_not_subsets_of_smaller_than.size(); ++ii) {

					const boost::dynamic_bitset<>& focal_point = intersections_of_not_subsets_of_smaller_than[ii];

					// add it to q_on_focal_points if it wasn't already there
					if(!q_on_focal_points[focal_point]){

						set_N_value<T>* inserted_focal_point = insert_focal_point(
							focal_point,
							commonality_equivalent,
							q_on_focal_points,
							q_on_focal_points_map
						);
						insert_neg_focal_point(
							q_on_focal_points.fod->set_negate(focal_point),
							inserted_focal_point->value,
							neg_b_on_neg_focal_points,
							neg_b_on_neg_focal_points_map
						);

						q_on_pure_focal_points.push_back(inserted_focal_point);
					}
				}

				const std::vector<boost::dynamic_bitset<> >& unions_of_not_subsets_of_smaller_than = neg_b_on_neg_focal_sets
												.unions_of_not_subsets_of_smaller_than(
													q_on_focal_points.fod->set_negate(setA),
													q_on_focal_points.fod->size() - q_on_pure_focal_points[i]->fod_elements.size()-1
												);
				// search for sets of higher cardinality that are not supersets of the current set
				for (size_t ii = 0; ii < unions_of_not_subsets_of_smaller_than.size(); ++ii) {

					const boost::dynamic_bitset<>& neg_focal_point = unions_of_not_subsets_of_smaller_than[ii];

					// add it to q_on_focal_points if it wasn't already there
					if(!neg_b_on_neg_focal_points[neg_focal_point]){
						const boost::dynamic_bitset<>& focal_point = q_on_focal_points.fod->set_negate(unions_of_not_subsets_of_smaller_than[ii]);

						set_N_value<T>* inserted_focal_point = insert_focal_point(
							focal_point,
							commonality_equivalent,
							q_on_focal_points,
							q_on_focal_points_map
						);
						insert_neg_focal_point(
							q_on_focal_points.fod->set_negate(focal_point),
							inserted_focal_point->value,
							neg_b_on_neg_focal_points,
							neg_b_on_neg_focal_points_map
						);

						q_on_pure_focal_points.push_back(inserted_focal_point);
					}
				}
			}
		}

		/*
		 * Principle: Recursively compute all weights on focal points
		 */
		static void compute_neg_v_special_elements(
				const powerset_btree<T>& neg_b_on_neg_focal_points,
				std::unordered_map<size_t, std::vector<set_N_value<T>* > >& neg_b_on_neg_focal_points_map,
				powerset_btree<T>& neg_v_special_elements,
				bool is_consonant){

			std::clog << "\nNegative focal points : \n";
			print<T>(std::clog, neg_b_on_neg_focal_points);

			std::clog << "\nBefore computing aggregations\n";

			const std::vector<std::vector<set_N_value<T>* >* >& ordered_vector = neg_b_on_neg_focal_points
																.get_vector_of_vectors_ordered_by_cardinality(neg_b_on_neg_focal_points_map);

			if(is_consonant){
				// when the structure is consonant, there is at most one focal set per cardinality
				for (size_t c = 1; c < ordered_vector.size(); ++c) {
					//std::clog << "\n======================= NEW WEIGHT\n";
					neg_v_special_elements.insert(
						(*ordered_vector[c])[0]->fod_elements,
						(*ordered_vector[c-1])[0]->value / (*ordered_vector[c])[0]->value
					);
					//std::clog << "\n======================= WEIGHT INSERTED\n";
				}
			}else{
				T val;
				for (size_t c = 1; c < ordered_vector.size(); ++c) {
					for (size_t i = 0; i < ordered_vector[c]->size(); ++i) {

						std::clog << "\n======================= NEW WEIGHT\n";
						std::clog << "\nWeight computation for " << to_string((*ordered_vector[c])[i]->fod_elements) << ":\n";
						std::clog << "Initialization with b_emptyset/b_A = " << (*ordered_vector[0])[0]->value << "/" << (*ordered_vector[c])[i]->value
								<< "\nThen, division by weights:\n";

						val = (*ordered_vector[0])[0]->value / (*ordered_vector[c])[i]->value;	// (*ordered_vector[0])[0] is the emptyset node

						const std::vector<set_N_value<T>* >& neg_v_subsets = neg_v_special_elements
																		.strict_subsets_of(
																			(*ordered_vector[c])[i]->fod_elements
																		);
						for (size_t ii = 0; ii < neg_v_subsets.size(); ++ii) {
							std::clog << to_string<T>(*neg_v_subsets[ii]) << std::endl;

							val /= neg_v_subsets[ii]->value;	// (*ordered_vector[0])[0] is the emptyset node
						}
						// insert weight among neg_v_special elements
						neg_v_special_elements.insert(
							(*ordered_vector[c])[i]->fod_elements,
							val
						);
						std::clog << "\n======================= WEIGHT INSERTED\n";
					}

				}
			}
		}

		/*
		 * Principle: check if there is a focal point supserset S of A contained in all the others.
		 * 			- If so, just apply: w(A) = q(A)^(-1) . q(S), and move on to the next weight computation
		 * 			- If not, proceed compute focal points exponents recursively
		 */
		static T compute_focal_point_aggregation(
				const commonality<T>& commonality_equivalent,
				const powerset_btree<T>& used_commonalities,
				//const std::vector<std::vector<set_N_value<T>* > >& used_q_elements_by_cardinality,
				const set_N_value<T>* A,
				powerset_btree<T>& special_elements){

			std::clog << "\n======================= NEW WEIGHT\n";

			const std::vector<fod_element* >& current_set = A->fod_elements;

			// initialize its weight with the inverse of its own commonality
			std::clog << "\nWeight computation for " << to_string(current_set) << ":\n-1\t <- " << to_string(current_set) << std::endl;
			T weight = 1/A->value;

			// check if there is a focal set in supersets contained in every focal set in supersets...

			// get supersets of current_set
			const std::vector<set_N_value<T>* >& supersets = used_commonalities
															.strict_supersets_of(
																	current_set
															);
			boost::dynamic_bitset<> I = supersets[0]->fod_elements;
			const boost::dynamic_bitset<> emptyset(used_commonalities.fod->size(), 0);
			size_t i = 1;
			for (; i < supersets.size(); ++i) {
				if (I != emptyset){
					I = I & used_commonalities.fod->to_set(supersets[i]->fod_elements);
				}else{
					break;
				}
			}
			for (i = 0; i < supersets.size(); ++i) {
				if (used_commonalities.fod->is_or_is_subset_of(used_commonalities.fod->to_set(supersets[i]->fod_elements), I)){
					return weight*supersets[i]->value;
				}
			}

			// ...if there is none, execute the general procedure

			//std::unordered_map<size_t, std::vector<set_N_value<T>* > > supersets_by_cardinality = used_commonalities
			//																					.strict_supersets_of_by_cardinality(
			//																						current_set
			//																					);
			//const std::vector<std::vector<set_N_value<T>* >* >& ordered_vector = used_commonalities
					//.get_vector_of_vectors_ordered_by_cardinality(supersets_by_cardinality);

			powerset_btree<size_t> polarities(*used_commonalities.fod, used_commonalities.block_size);
			size_t polarity;

			for (i = 0; i < supersets.size(); ++i) {
				// By definition of the conjunctive decomposition, the polarity for this
				// superset in the computation of this current weight is :
				// ->  1 	if (sc - c) is odd
				// -> -1 	if (sc - c) is even
				if(((supersets[i]->fod_elements.size() - current_set.size()) & 1) > 0){	// check the first bit for parity check
					polarity = 1;
				}else{
					polarity = -1;
				}
				polarities.insert(supersets[i]->fod_elements, polarity);
			}

			std::clog << "\nFocal points : \n";
			print<T>(std::clog, polarities);

			std::unordered_map<size_t, std::vector<set_N_value<T>* > > polarities_by_cardinality = polarities
																								.strict_supersets_of_by_cardinality(
																									current_set
																								);
			const std::vector<std::vector<set_N_value<T>* >* >& ordered_vector = polarities
								.get_vector_of_vectors_ordered_by_cardinality(polarities_by_cardinality);

			size_t s;
			for (size_t c = 0; c < ordered_vector.size(); ++c) {
				for (size_t i = 0; i < ordered_vector[c]->size(); ++i) {
					s = -1;
					const std::vector<set_N_value<T>* >& polarities_supersets = polarities
																	.strict_subsets_of(
																			(*ordered_vector[c])[i]->fod_elements
																	);
					for (size_t ii = 0; ii < polarities_supersets.size(); ++ii) {
						s += polarities_supersets[ii];
					}
					set_N_value<size_t>* inserted_polarity = polarities.insert((*ordered_vector[c])[i]->fod_elements, -s);
					std::clog << to_string<size_t>(*inserted_polarity) << std::endl;

					weight *= pow(
								used_commonalities[(*ordered_vector[c])[i]->fod_elements],
								-s
							);
				}
			}

			// insert weight among special elements
			special_elements.insert(current_set, weight);
			std::clog << "\n======================= WEIGHT INSERTED\n";
			return weight;
		}
	};

} // namespace ow_bft

#endif // OW_BFT_CONJUNCTIVE_DECOMPOSITION_HPP
