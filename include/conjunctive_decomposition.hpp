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

		static void compute_values_for_special_elements(const commonality<T>& commonality_equivalent, powerset_btree<T>& special_elements) {
			if(commonality_equivalent.get_mass_equivalent().is_dogmatic()){
				std::cerr << "\nThe conjunctive decomposition is not defined for a dogmatic BBA.";
				exit(EXIT_FAILURE);
			}

			if(&commonality_equivalent.get_FOD() != special_elements.fod){
				std::cerr << "\nTrying to compute values of a conjunctive decomposition from a commonality function\n"
						  << "with given special elements that refer to FOD elements different from the ones of this commonality function.\n";
				exit(1);
			}

			// Initialize used_commonalities
			// (commonality value for elements with value of conjunctive weight different from 1)
			// with all focal elements as defined on masses.
			powerset_btree<T> used_commonalities(commonality_equivalent.get_special_elements());

			std::vector<std::vector<set_N_value<T>* > > used_q_elements_by_cardinality = used_commonalities.elements_by_set_cardinality();

			std::clog << "\nused commonalities copied\n";

			compute_focal_points(commonality_equivalent, used_commonalities, used_q_elements_by_cardinality);
			std::clog << "\nfocal points computed\n";

			compute_aggregations(commonality_equivalent, used_commonalities, used_q_elements_by_cardinality, special_elements);
			std::clog << "\naggregations computed\n";
		}

		static void to_commonality_special_elements(const powerset_btree<T>& w_tree, powerset_btree<T>& q_tree) {
			if(w_tree.fod != q_tree.fod){
				std::cerr << "\nTrying to transform a conjunctive decomposition into a commonality function\n"
						  << "with given special elements that refer to FOD elements different from the ones of this conjunctive decomposition.\n";
				exit(1);
			}

			to_commonality_special_elements(w_tree, w_tree.elements_by_set_cardinality(), q_tree);
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

			to_commonality_special_elements(w_tree, m_tree.elements_by_set_cardinality(), q_tree);
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

		static void to_commonality_special_elements(const powerset_btree<T>& w_tree, const std::vector<std::vector<set_N_value<T>* > >& f_elements, powerset_btree<T>& q_tree) {

			const std::vector<set_N_value<T>* >& w_elements = w_tree.elements();

			T w_product = 1;
			for (size_t i = 0; i < w_elements.size(); ++i) {
				w_product *= w_elements[i]->value;
			}

			// set FOD to w_product
			q_tree.set_value_of_sub_fod_of_size(q_tree.fod->size(), w_product);
			// set emptyset to 1
			q_tree.set_value_of_sub_fod_of_size(0, 1);

			for (size_t c = 1; c < f_elements.size()-1; ++c) {
				for (size_t i = 0; i < f_elements[c].size(); ++i) {

					const std::vector<set_N_value<T>* >& w_i_supersets = w_tree.supersets_of(f_elements[c][i]->fod_elements);
					T q_i_val = w_product;

					for (size_t ii = 0; ii < w_i_supersets.size(); ++ii) {
						q_i_val /= w_i_supersets[ii]->value;
					}

					q_tree.insert(f_elements[c][i]->fod_elements, q_i_val);
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

		static void compute_focal_point_and_add_it_to_used_commonalities(
				const commonality<T>& commonality_equivalent,
				const boost::dynamic_bitset<>& setA,
				const boost::dynamic_bitset<>& setB,
				powerset_btree<T>& used_commonalities,
				std::vector<std::vector<set_N_value<T>* > >& used_q_elements_by_cardinality
				) {

			// compute their focal point
			const boost::dynamic_bitset<>& focal_point = used_commonalities.fod->set_intersection(setA, setB);

			std::clog << "\nsetA : " << setA << std::endl;
			std::clog << "\nsetB : " << setB << std::endl;
			std::clog << "\nfocal point : " << focal_point << std::endl;

			// and add it to used_commonalities if it wasn't already there
			if(!used_commonalities[focal_point]){
				std::clog << "\nNEW INSERTED FOCAL POINT :\n";
				set_N_value<T>* inserted_focal_point = used_commonalities.insert(
						focal_point,
						commonality_equivalent.find(focal_point)
				);
				std::clog << inserted_focal_point->fod_elements << std::endl;
				used_q_elements_by_cardinality[inserted_focal_point->fod_elements.size()].push_back(
						inserted_focal_point
				);
			}
		}

		static void compute_focal_points_from_focal_points(
				const commonality<T>& commonality_equivalent,
				powerset_btree<T>& used_commonalities,
				std::vector<std::vector<set_N_value<T>* > >& used_q_elements_by_cardinality,
				size_t c
				//const std::vector<std::vector<set_N_value<T>* > >& focal_q_elements_by_cardinality
				) {

			// for each set of cardinality c
			for (size_t i = 0; i < used_q_elements_by_cardinality[c].size(); ++i) {

				const boost::dynamic_bitset<>& setA = used_commonalities.fod->to_set(used_q_elements_by_cardinality[c][i]->fod_elements);

				// for each other set of cardinality c (that makes a unique pair with it, no matter the order)
				for (size_t ii = i+1; ii < used_q_elements_by_cardinality[c].size(); ++ii) {

					compute_focal_point_and_add_it_to_used_commonalities(
						commonality_equivalent,
						setA,
						used_commonalities.fod->to_set(used_q_elements_by_cardinality[c][ii]->fod_elements),
						used_commonalities,
						used_q_elements_by_cardinality
					);
				}

				// for each cardinality cc > c
				// (Don't evaluate intersections with FOD as it will simply result in the set other than FOD
				// (these results were already considered with the copy of this->commonality_equivalent.focal_elements))
				for (size_t cc = c+1; cc < used_commonalities.fod->size(); ++cc) {

					// for each set of cardinality cc
					for (size_t ii = 0; ii < used_q_elements_by_cardinality[cc].size(); ++ii) {

						compute_focal_point_and_add_it_to_used_commonalities(
							commonality_equivalent,
							setA,
							used_commonalities.fod->to_set(used_q_elements_by_cardinality[cc][ii]->fod_elements),
							used_commonalities,
							used_q_elements_by_cardinality
						);
					}
				}
			}
		}
/*
		static void compute_focal_points_from_focal_points(
				const commonality<T>& commonality_equivalent,
				powerset_btree<T>& used_commonalities,
				std::vector<std::vector<set_N_value<T>* > >& used_q_elements_by_cardinality,
				size_t& c) {

			// for each set of cardinality c
			for (size_t i = 0; i < used_q_elements_by_cardinality[c].size(); ++i) {

				const boost::dynamic_bitset<>& setA = used_commonalities.fod->to_set(used_q_elements_by_cardinality[c][i]->fod_elements);

				// for each other set of cardinality c (that makes a unique pair with it, no matter the order)
				for (size_t ii = i+1; ii < used_q_elements_by_cardinality[c].size(); ++ii) {

					compute_focal_point_and_add_it_to_used_commonalities(
						commonality_equivalent,
						setA,
						used_commonalities.fod->to_set(used_q_elements_by_cardinality[c][ii]->fod_elements),
						used_commonalities,
						used_q_elements_by_cardinality
					);
				}
			}
		}*/

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
				powerset_btree<T>& used_commonalities,
				std::vector<std::vector<set_N_value<T>* > >& used_q_elements_by_cardinality
			) {

			// if there is only 1 or 0 element in FOD, then there can be no focal point other than the FOD
			if(used_commonalities.fod->size() < 2){
				return;
			}

			// for each cardinality
			for (size_t c = used_commonalities.fod->size()-1; c > 0 ; --c) {
				std::clog << "\n-> c = " << c << std::endl;

				compute_focal_points_from_focal_points(
						commonality_equivalent, used_commonalities, used_q_elements_by_cardinality, c);
			}
		}

		static void compute_aggregations(
				const commonality<T>& commonality_equivalent,
				const powerset_btree<T>& used_commonalities,
				const std::vector<std::vector<set_N_value<T>* > >& used_q_elements_by_cardinality,
				powerset_btree<T>& special_elements){

			T weight;


			std::clog << "\nFocal points : \n";
			print<T>(std::clog, used_commonalities);

			std::clog << "\nBefore computing aggregations\n";

			T q_FOD = used_q_elements_by_cardinality[used_commonalities.fod->size()][0]->value;

			// for each cardinality except the one of FOD
			for (size_t c = used_commonalities.fod->size()-1; c > 0; --c) {

				// for each element of used_q_elements_by_cardinality of cardinality c
				for (size_t i = 0; i < used_q_elements_by_cardinality[c].size(); ++i) {

					std::clog << "\n======================= NEW WEIGHT\n";

					const std::vector<fod_element* >& current_set = used_q_elements_by_cardinality[c][i]->fod_elements;

					// get supersets of used_q_elements_by_cardinality[c][i]
					const std::vector<set_N_value<T>* >& supersets = special_elements.strict_supersets_of(current_set);

					// initialize its weight with the inverse of its own commonality

					std::clog << "\nWeight computation for " << to_string(current_set) << ":\n";
					std::clog << "Initialization with q_FOD/q_A, then division by:\n";
					weight = q_FOD/used_q_elements_by_cardinality[c][i]->value;

					for (size_t si = 0; si < supersets.size(); ++si) {

						std::clog << to_string<T>(*supersets[si]) << std::endl;

						weight /= supersets[si]->value;
					}

					// insert weight among special elements
					special_elements.insert(current_set, weight);

					std::clog << "\n======================= WEIGHT INSERTED\n";
				}
			}

			//emptyset weight computation to avoid issues with size_t which can't hold negative values (i.e. avoid --c when c=0 which causes infinite loop)
			if (used_q_elements_by_cardinality[0].size() > 0) {
				std::clog << "\n======================= NEW WEIGHT\n";

				const std::vector<fod_element* >& current_set = used_q_elements_by_cardinality[0][0]->fod_elements;

				// get supersets of used_q_elements_by_cardinality[c][i]
				const std::vector<set_N_value<T>* >& supersets = special_elements.strict_supersets_of(current_set);

				// initialize its weight with the inverse of its own commonality

				std::clog << "\nWeight computation for " << to_string(current_set) << ":\n";
				std::clog << "Initialization with q_FOD/q_A, then division by:\n";
				weight = q_FOD/used_q_elements_by_cardinality[0][0]->value;

				for (size_t si = 0; si < supersets.size(); ++si) {

					std::clog << to_string<T>(*supersets[si]) << std::endl;

					weight /= supersets[si]->value;
				}

				// insert weight among special elements
				special_elements.insert(current_set, weight);

				std::clog << "\n======================= WEIGHT INSERTED\n";
			}
		}

		/*
		 * TODO: add a small code to check if the smallest focal point supserset S of A is contained in all the others:
		 * 			- If so, just apply: w(A) = q(A)^(-1) . q(S), and move on to the next weight computation
		 * 			- If not, proceed as it already is
		 */
		static void compute_focal_point_aggregation(
				const commonality<T>& commonality_equivalent,
				const powerset_btree<T>& used_commonalities,
				//const std::vector<std::vector<set_N_value<T>* > >& used_q_elements_by_cardinality,
				const set_N_value<T>* A,
				powerset_btree<T>& special_elements){

			powerset_btree<long> polarities(*used_commonalities.fod, used_commonalities.block_size);
			T weight;
			long polarity;

			// for each cardinality except the one of FOD
			//for (size_t c = 0; c < used_commonalities.fod->size(); ++c) {

				// for each element of used_q_elements_by_cardinality of cardinality c
				//for (size_t i = 0; i < used_q_elements_by_cardinality[c].size(); ++i) {

			std::clog << "\n======================= NEW WEIGHT\n";

			const std::vector<fod_element* >& current_set = A->fod_elements;

			// get supersets of used_q_elements_by_cardinality[c][i]
			const std::vector<std::vector<set_N_value<T>* > >& supersets = used_commonalities
																		.supersets_of_by_cardinality(
																				current_set
																		);

			if(supersets)
			set_N_value<T>* S = nullptr;
			size_t sc = current_set.size()+1;
			bool found = false;
			while (!found && sc < supersets.size()) {
				for (size_t si = 0; si < supersets[sc].size(); ++si) {
					S = supersets[sc][si];
					found = true;
					break;
				}
				++sc;
			}
			bool is_simple = true;
			while (!found && sc < supersets.size()) {
				for (size_t si = 0; si < supersets[sc].size(); ++si) {
					S = supersets[sc][si];
					found = true;
					break;
				}
				++sc;
			}

			// initialize polarities for the computation of this current weight
			for (size_t sc = current_set.size()+1; sc < supersets.size(); ++sc) {
				for (size_t si = 0; si < supersets[sc].size(); ++si) {
					// By definition of the conjunctive decomposition, the polarity for this
					// superset in the computation of this current weight is :
					// ->  1 	if (sc - c) is odd
					// -> -1 	if (sc - c) is even
					if(((sc - current_set.size()) & 1) > 0){	// check the first bit for parity check
						polarities.insert(supersets[sc][si]->fod_elements, 1);
					}else{
						polarities.insert(supersets[sc][si]->fod_elements, -1);
					}
				}
			}

			std::clog << "\nFocal points : \n";
			print<T>(std::clog, polarities);

			// initialize its weight with the inverse of its own commonality

			std::clog << "\nWeight computation for " << to_string(current_set) << ":\n-1\t <- " << to_string(current_set) << std::endl;
			weight = 1/A->value;

			for (size_t sc = current_set.size() + 1; sc < supersets.size(); ++sc) {
				for (size_t si = 0; si < supersets[sc].size(); ++si) {

					// all supersets have at least used_q_elements_by_cardinality[c][i]->fod_elements as subset,
					// hence the 1 which stands for the polarity -1 of the weight initialization
					polarity = 1;

					const std::vector<set_N_value<long>* >& subsets = polarities.strict_subsets_of(supersets[sc][si]->fod_elements);

					// compensate for the plurality of all polarities corresponding to sets that are influenced by supersets[sc][si]
					for (size_t i = 0; i < subsets.size(); ++i) {
						polarity -= subsets[i]->value;
					}
					const set_N_value<long>* inserted_polarity = polarities.insert(supersets[sc][si]->fod_elements, polarity);
					std::clog << to_string<long>(*inserted_polarity) << std::endl;

					weight *= pow(
							supersets[sc][si]->value,
							polarity
						);
				}
			}

			// insert weight among special elements
			special_elements.insert(current_set, weight);

			std::clog << "\n======================= WEIGHT INSERTED\n";

			polarities.nullify();
				//}
			//}
		}
	};

} // namespace ow_bft

#endif // OW_BFT_CONJUNCTIVE_DECOMPOSITION_HPP
