#ifndef OW_BFT_CANONICAL_DECOMPOSITION_HPP
#define OW_BFT_CANONICAL_DECOMPOSITION_HPP

#include <aggregate.hpp>
#include <commonality.hpp>
#include <mass.hpp>

namespace ow_bft{

	template <typename T = double>
	class canonical_decomposition : public aggregate<T>{

	public:
		virtual ~canonical_decomposition(){}

	protected:

		static set_N_value<T>* insert_neg_focal_p(
				const boost::dynamic_bitset<>& neg_focal_p,
				T value,
				powerset_btree<T>& neg_func_on_neg_focal_p,
				std::unordered_map<size_t, std::vector<set_N_value<T>* > >& neg_func_on_neg_focal_p_map
				) {

			std::clog << "\nNEW INSERTED NEGATIVE FOCAL P :\n";
			std::clog << "\nnegative focal p : " << neg_focal_p << std::endl;

			set_N_value<T>* inserted_neg_focal_p = neg_func_on_neg_focal_p.insert(
					neg_focal_p,
					value
			);
			std::clog << inserted_neg_focal_p->set << std::endl;

			size_t c = inserted_neg_focal_p->set.count();

			neg_func_on_neg_focal_p_map[c].push_back(
					inserted_neg_focal_p
			);
			return inserted_neg_focal_p;
		}

		static set_N_value<T>* insert_focal_p(
				const boost::dynamic_bitset<>& focal_p,
				const mass_aggregate<T>& func_equivalent,
				powerset_btree<T>& func_on_focal_p,
				std::unordered_map<size_t, std::vector<set_N_value<T>* > >& func_on_focal_p_map
				) {

			set_N_value<T>* inserted_focal_p = func_on_focal_p.insert(
					focal_p,
					func_equivalent.find(focal_p)
			);
			std::clog << "\nNEW INSERTED FOCAL POINT :\n";
			std::clog << inserted_focal_p->set << std::endl;

			func_on_focal_p_map[inserted_focal_p->set.count()].push_back(
					inserted_focal_p
			);
			return inserted_focal_p;
		}
/*
		static bool linear_analysis_of_focal_points(
				const commonality<T>& commonality_equivalent,
				powerset_btree<T>& q_on_focal_points,
				powerset_btree<T>& neg_b_on_neg_focal_points,
				std::unordered_map<size_t, std::vector<set_N_value<T>* > >& q_on_focal_points_map,
				std::unordered_map<size_t, std::vector<set_N_value<T>* > >& neg_b_on_neg_focal_points_map,
				std::vector<set_N_value<T>* >& q_on_pure_focal_points
				){

			bool generated_quasi_bayesian = false;
			std::vector<set_N_value<T>* > focal_sets_except_FOD;
			focal_sets_except_FOD.reserve(q_on_focal_points.fod->size()-1);

			for(auto kv : q_on_focal_points_map) {
			    if(kv.first < q_on_focal_points.fod->size()){
					for (size_t i = 0; i < kv.second.size(); ++i) {
						focal_sets_except_FOD.emplace_back(kv.second[i]);
					}
			    }
			}

			boost::dynamic_bitset<> U = q_on_focal_points.fod->to_set(focal_sets_except_FOD[0]->fod_elements);
			const boost::dynamic_bitset<> emptyset(q_on_focal_points.fod->size());
			std::clog << "\nU = "<< focal_sets_except_FOD[0]->fod_elements;

			for (size_t i = 1; i < focal_sets_except_FOD.size(); ++i) {
				std::clog << "\nI = intersection(U, "<< focal_sets_except_FOD[i]->fod_elements << ")\n";
				const boost::dynamic_bitset<>& A = q_on_focal_points.fod->to_set(focal_sets_except_FOD[i]->fod_elements);
				const boost::dynamic_bitset<>& I = q_on_focal_points.fod->set_intersection(U, A);
				const size_t& I_card = I.count();
				std::clog << "|I| = " << I_card << std::endl;
				if(I_card <= 1){
					// add it to q_on_focal_points if it wasn't already there
					if(!q_on_focal_points[I]){
						set_N_value<T>* inserted_focal_point = insert_focal_p(
							I,
							commonality_equivalent,
							q_on_focal_points,
							q_on_focal_points_map
						);
						insert_neg_focal_p(
							q_on_focal_points.fod->set_negate(I),
							inserted_focal_point->value,
							neg_b_on_neg_focal_points,
							neg_b_on_neg_focal_points_map
						);
						q_on_pure_focal_points.push_back(inserted_focal_point);
					}

					if(I_card == 1){
						generated_quasi_bayesian = true;
					}
				}else{
					std::clog << "=> Linear analysis aborted.\n";
					return false;
				}
				U = q_on_focal_points.fod->set_union(
					U,
					A
				);
			}

			if(generated_quasi_bayesian){
				// add also emptyset to q_on_focal_points if it wasn't already there
				if(!q_on_focal_points[emptyset]){
					set_N_value<T>* inserted_focal_point = insert_focal_p(
						emptyset,
						commonality_equivalent,
						q_on_focal_points,
						q_on_focal_points_map
					);
					insert_neg_focal_p(
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

		static void consonance_check(
				std::vector<set_N_value<T>* >& q_on_pure_focal_points,
				const std::vector<std::vector<set_N_value<T>* >* >& q_on_focal_sets_ordered_vector,
				const FOD*& fod,
				bool& is_consonant
				){
			// Consonance check
			if(q_on_pure_focal_points.size() == 0 && q_on_focal_sets_ordered_vector[0]->size() == 1){
				is_consonant = true;
				for (size_t i = 1; i < q_on_focal_sets_ordered_vector.size()-1; ++i) {
					if(q_on_focal_sets_ordered_vector[i]->size() > 1 || !fod->is_or_is_subset_of(
							(*q_on_focal_sets_ordered_vector[i-1])[0]->set,
							(*q_on_focal_sets_ordered_vector[i])[0]->set
						  )
					){
						is_consonant = false;
						break;
					}
				}
			}
		}

		static bool linear_analysis_of_focal_points(
				const commonality<T>& commonality_equivalent,
				powerset_btree<T>& q_on_focal_points,
				powerset_btree<T>& neg_b_on_neg_focal_points,
				std::unordered_map<size_t, std::vector<set_N_value<T>* > >& q_on_focal_points_map,
				std::unordered_map<size_t, std::vector<set_N_value<T>* > >& neg_b_on_neg_focal_points_map,
				std::vector<set_N_value<T>* >& q_on_pure_focal_points,
				bool& is_quasi_bayesian
				){

			is_quasi_bayesian = false;
			bool has_generated_quasi_bayesian = false;
			std::vector<set_N_value<T>* > focal_sets_except_FOD;
			focal_sets_except_FOD.reserve(q_on_focal_points.fod->size()-1);

			for(auto kv : q_on_focal_points_map) {
			    if(kv.first < q_on_focal_points.fod->size()){
					for (size_t i = 0; i < kv.second.size(); ++i) {
						focal_sets_except_FOD.emplace_back(kv.second[i]);
					}
			    }
			}

			boost::dynamic_bitset<> U = focal_sets_except_FOD[0]->set;
			const boost::dynamic_bitset<> emptyset(q_on_focal_points.fod->size());
			std::clog << "\nU = "<< focal_sets_except_FOD[0]->set;

			for (size_t i = 1; i < focal_sets_except_FOD.size(); ++i) {
				std::clog << "\nI = intersection(U, "<< focal_sets_except_FOD[i]->set << ")\n";
				const boost::dynamic_bitset<>& A = focal_sets_except_FOD[i]->set;
				const boost::dynamic_bitset<>& I = q_on_focal_points.fod->set_intersection(U, A);
				const size_t& I_card = I.count();
				std::clog << "|I| = " << I_card << std::endl;
				if(I_card <= 1){
					// add it to q_on_focal_points if it wasn't already there
					if(!q_on_focal_points[I]){
						set_N_value<T>* inserted_focal_point = canonical_decomposition<T>::insert_focal_p(
							I,
							commonality_equivalent,
							q_on_focal_points,
							q_on_focal_points_map
						);
						canonical_decomposition<T>::insert_neg_focal_p(
							q_on_focal_points.fod->set_negate(I),
							inserted_focal_point->value,
							neg_b_on_neg_focal_points,
							neg_b_on_neg_focal_points_map
						);
						q_on_pure_focal_points.push_back(inserted_focal_point);
					}

					if(I_card == 1){
						has_generated_quasi_bayesian = true;
					}
				}else{
					std::clog << "=> Linear analysis aborted.\n";
					return false;
				}
				U = q_on_focal_points.fod->set_union(
					U,
					A
				);
			}

			if(has_generated_quasi_bayesian){
				// add also emptyset to q_on_focal_points if it wasn't already there
				if(!q_on_focal_points[emptyset]){
					set_N_value<T>* inserted_focal_point = canonical_decomposition<T>::insert_focal_p(
						emptyset,
						commonality_equivalent,
						q_on_focal_points,
						q_on_focal_points_map
					);
					canonical_decomposition<T>::insert_neg_focal_p(
						q_on_focal_points.fod->set_negate(emptyset),
						inserted_focal_point->value,
						neg_b_on_neg_focal_points,
						neg_b_on_neg_focal_points_map
					);
					q_on_pure_focal_points.push_back(inserted_focal_point);
				}
			}else{
				is_quasi_bayesian = true;
			}
			std::clog << "=> Successful linear analysis of focal points.\n";
			return true;
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
				const powerset_btree<T>& neg_b_on_neg_focal_sets,
				powerset_btree<T>& q_on_focal_points,
				powerset_btree<T>& neg_b_on_neg_focal_points,
				std::unordered_map<size_t, std::vector<set_N_value<T>* > >& q_on_focal_points_map,
				std::unordered_map<size_t, std::vector<set_N_value<T>* > >& neg_b_on_neg_focal_points_map,
				const std::vector<std::vector<set_N_value<T>* >* >& q_on_focal_sets_ordered_vector,
				std::vector<set_N_value<T>* >& q_on_pure_focal_points
				) {

			// for each cardinality of focal set, from the smallest to the biggest, except FOD
			for (size_t c = 0; c < q_on_focal_sets_ordered_vector.size()-1; ++c) {
				for (size_t i = 0; i < q_on_focal_sets_ordered_vector[c]->size(); ++i) {

					const boost::dynamic_bitset<>& setA = (*q_on_focal_sets_ordered_vector[c])[i]->set;

					// search for sets of same cardinality
					for (size_t ii = 0; ii < i; ++ii) {

						const boost::dynamic_bitset<>& setB = (*q_on_focal_sets_ordered_vector[c])[ii]->set;

						// compute their focal point
						const boost::dynamic_bitset<>& focal_point = q_on_focal_points.fod->set_intersection(setA, setB);

						// add it to q_on_focal_points if it wasn't already there
						if(!q_on_focal_points[focal_point]){
							set_N_value<T>* inserted_focal_point = insert_focal_p(
								focal_point,
								commonality_equivalent,
								q_on_focal_points,
								q_on_focal_points_map
							);
							insert_neg_focal_p(
								q_on_focal_points.fod->set_negate(focal_point),
								inserted_focal_point->value,
								neg_b_on_neg_focal_points,
								neg_b_on_neg_focal_points_map
							);

							q_on_pure_focal_points.push_back(inserted_focal_point);
						}
					}

					const std::vector<boost::dynamic_bitset<> >& intersections_of_not_subsets_of_smaller_than = commonality_equivalent
													.get_special_elements()
													.intersections_of_not_subsets_of_smaller_than(
														setA,
														(*q_on_focal_sets_ordered_vector[c])[i]->set.count()-1
													);
					// search for sets of higher cardinality that are not supersets of the current set
					for (size_t ii = 0; ii < intersections_of_not_subsets_of_smaller_than.size(); ++ii) {

						const boost::dynamic_bitset<>& focal_point = intersections_of_not_subsets_of_smaller_than[ii];

						// add it to q_on_focal_points if it wasn't already there
						if(!q_on_focal_points[focal_point]){
							set_N_value<T>* inserted_focal_point = insert_focal_p(
								focal_point,
								commonality_equivalent,
								q_on_focal_points,
								q_on_focal_points_map
							);
							insert_neg_focal_p(
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

				const boost::dynamic_bitset<>& setA = q_on_pure_focal_points[i]->set;

				const std::vector<boost::dynamic_bitset<> >& intersections_of_not_subsets_of_smaller_than = commonality_equivalent
												.get_special_elements()
												.intersections_of_not_subsets_of_smaller_than(
													setA,
													q_on_pure_focal_points[i]->set.count()-1
												);
				// search for sets of lower cardinality that are not subsets of the current set
				for (size_t ii = 0; ii < intersections_of_not_subsets_of_smaller_than.size(); ++ii) {

					const boost::dynamic_bitset<>& focal_point = intersections_of_not_subsets_of_smaller_than[ii];

					// add it to q_on_focal_points if it wasn't already there
					if(!q_on_focal_points[focal_point]){

						set_N_value<T>* inserted_focal_point = insert_focal_p(
							focal_point,
							commonality_equivalent,
							q_on_focal_points,
							q_on_focal_points_map
						);
						insert_neg_focal_p(
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
													q_on_focal_points.fod->size() - q_on_pure_focal_points[i]->set.count()-1
												);
				// search for sets of higher cardinality that are not supersets of the current set
				for (size_t ii = 0; ii < unions_of_not_subsets_of_smaller_than.size(); ++ii) {

					const boost::dynamic_bitset<>& neg_focal_point = unions_of_not_subsets_of_smaller_than[ii];

					// add it to q_on_focal_points if it wasn't already there
					if(!neg_b_on_neg_focal_points[neg_focal_point]){
						const boost::dynamic_bitset<>& focal_point = q_on_focal_points.fod->set_negate(unions_of_not_subsets_of_smaller_than[ii]);

						set_N_value<T>* inserted_focal_point = insert_focal_p(
							focal_point,
							commonality_equivalent,
							q_on_focal_points,
							q_on_focal_points_map
						);
						insert_neg_focal_p(
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
		 * Principle: Recursively compute all weights on focal planes
		 */
		static void compute_disjunctive_decomposition(
				const powerset_btree<T>& neg_b_on_neg_focal_sets,
				const powerset_btree<T>& neg_b_on_neg_focal_points,
				std::unordered_map<size_t, std::vector<set_N_value<T>* > >& neg_b_on_neg_focal_points_map,
				powerset_btree<T>& neg_v_special_elements,
				const bool& is_quasi_bayesian,
				const bool& is_consonant){

			//std::clog << "\nNegative focal points : \n";
			//print<T>(std::clog, neg_b_on_neg_focal_points);

			T val;
			if(is_quasi_bayesian){
				for(auto kv : neg_b_on_neg_focal_points_map) {
					if(kv.first != 0 && kv.first != neg_b_on_neg_focal_points.fod->size()){
						for (size_t i = 0; i < kv.second.size(); ++i) {

							std::clog << "\n======================= NEW WEIGHT\n";
							std::clog << "\nWeight computation for " << neg_b_on_neg_focal_points.fod->to_string(kv.second[i]->set) << ":\n";
							std::clog << "b(" << neg_b_on_neg_focal_points.fod->to_string(neg_b_on_neg_focal_points_map[0][0]->set) << ") / b(" << neg_b_on_neg_focal_points.fod->to_string(kv.second[i]->set)
									<< ") = " << neg_b_on_neg_focal_points_map[0][0]->value << "/" << kv.second[i]->value << "\n";

							val = neg_b_on_neg_focal_points_map[0][0]->value / kv.second[i]->value;	// neg_b_on_neg_focal_points_map[0][0] is the emptyset node

							std::clog << "---------------------\n";
							std::clog << "= " << val << "\n";

							// insert weight among neg_v_special elements
							neg_v_special_elements.insert(
								kv.second[i]->set,
								val
							);
							std::clog << "\n======================= WEIGHT INSERTED\n";
						}
					}
				}
				if(neg_b_on_neg_focal_points_map.find(neg_b_on_neg_focal_points.fod->size()) != neg_b_on_neg_focal_points_map.end()){
					std::clog << "\n======================= NEW WEIGHT\n";
					std::clog << "\nWeight computation for " << neg_b_on_neg_focal_points.fod->to_string(neg_b_on_neg_focal_points_map[neg_b_on_neg_focal_points.fod->size()][0]->set) << ":\n";
					std::clog << "Initialization with b(" << neg_b_on_neg_focal_points.fod->to_string(neg_b_on_neg_focal_points_map[0][0]->set) << ") / b(" << neg_b_on_neg_focal_points.fod->to_string(neg_b_on_neg_focal_points_map[neg_b_on_neg_focal_points.fod->size()][0]->set)
							<< ") = " << neg_b_on_neg_focal_points_map[0][0]->value << "/" << neg_b_on_neg_focal_points_map[neg_b_on_neg_focal_points.fod->size()][0]->value << "\n";

					val = neg_b_on_neg_focal_points_map[0][0]->value / neg_b_on_neg_focal_points_map[neg_b_on_neg_focal_points.fod->size()][0]->value;	// neg_b_on_neg_focal_points_map[0][0] is the emptyset node

					const std::vector<set_N_value<T>* >& neg_v_subsets = neg_v_special_elements
																	.strict_subsets_of(
																		neg_b_on_neg_focal_points_map[neg_b_on_neg_focal_points.fod->size()][0]->set
																	);
					if(neg_v_subsets.size() > 0){
						std::clog << "Then, division by weights:\n";
					}
					for (size_t ii = 0; ii < neg_v_subsets.size(); ++ii) {
						std::clog << to_string<T>(*neg_v_subsets[ii], *neg_b_on_neg_focal_points.fod) << std::endl;

						val /= neg_v_subsets[ii]->value;	// (*ordered_vector[0])[0] is the emptyset node
					}
					std::clog << "---------------------\n";
					std::clog << "= " << val << "\n";

					// insert weight among neg_v_special elements
					neg_v_special_elements.insert(
						neg_b_on_neg_focal_points_map[neg_b_on_neg_focal_points.fod->size()][0]->set,
						val
					);
					std::clog << "\n======================= WEIGHT INSERTED\n";
				}
			}else{
				const std::vector<std::vector<set_N_value<T>* >* >& ordered_vector = neg_b_on_neg_focal_points
																		.get_vector_of_vectors_ordered_by_cardinality(neg_b_on_neg_focal_points_map);

				if(is_consonant){

					for (size_t c = 1; c < ordered_vector.size(); ++c) {
						std::clog << "\n======================= NEW WEIGHT\n";
						std::clog << "\nWeight computation for " << neg_b_on_neg_focal_points.fod->to_string((*ordered_vector[c])[0]->set) << ":\n";

						val = (*ordered_vector[c-1])[0]->value / (*ordered_vector[c])[0]->value;
						std::clog << "b(" << neg_b_on_neg_focal_points.fod->to_string((*ordered_vector[c-1])[0]->set) << ") / b(" << neg_b_on_neg_focal_points.fod->to_string((*ordered_vector[c])[0]->set)
														<< ") = " << val << "\n";

						neg_v_special_elements.insert(
							(*ordered_vector[c])[0]->set,
							val
						);
						std::clog << "\n======================= WEIGHT INSERTED\n";
					}
				}else{
					for (size_t c = 1; c < ordered_vector.size(); ++c) {
						for (size_t i = 0; i < ordered_vector[c]->size(); ++i) {

							std::clog << "\n======================= NEW WEIGHT\n";
							std::clog << "\nWeight computation for " << neg_b_on_neg_focal_points.fod->to_string((*ordered_vector[c])[i]->set) << ":\n";
							std::clog << "Initialization with b(" << neg_b_on_neg_focal_points.fod->to_string((*ordered_vector[0])[0]->set) << ") / b(" << neg_b_on_neg_focal_points.fod->to_string((*ordered_vector[c])[i]->set)
									<< ") = " << (*ordered_vector[0])[0]->value << "/" << (*ordered_vector[c])[i]->value << "\n";

							val = (*ordered_vector[0])[0]->value / (*ordered_vector[c])[i]->value;	// (*ordered_vector[0])[0] is the emptyset node

							const std::vector<set_N_value<T>* >& neg_v_subsets = neg_v_special_elements
																			.strict_subsets_of(
																				(*ordered_vector[c])[i]->set
																			);
							if(neg_v_subsets.size() > 0){
								std::clog << "Then, division by weights:\n";
							}
							for (size_t ii = 0; ii < neg_v_subsets.size(); ++ii) {
								std::clog << to_string<T>(*neg_v_subsets[ii], *neg_b_on_neg_focal_points.fod) << std::endl;

								val /= neg_v_subsets[ii]->value;	// (*ordered_vector[0])[0] is the emptyset node
							}
							std::clog << "---------------------\n";
							std::clog << "= " << val << "\n";

							// insert weight among neg_v_special elements
							neg_v_special_elements.insert(
								(*ordered_vector[c])[i]->set,
								val
							);
							std::clog << "\n======================= WEIGHT INSERTED\n";
						}
					}
				}
			}
		}

		/*
		 * Principle: check if there is a focal point supserset S of A contained in all the others.
		 * 			- If so, just apply: w(A) = q(A)^(-1) . q(S), and move on to the next weight computation
		 * 			- If not, proceed compute focal points exponents recursively
		 */
		static T compute_conjunctive_weight_at(
				const set_N_value<T>* A,
				const commonality<T>& commonality_equivalent,
				const powerset_btree<T>& q_on_focal_points
				){

			std::clog << "\n======================= NEW WEIGHT\n";

			const boost::dynamic_bitset<>& current_set = A->set;

			// initialize its weight with the inverse of its own commonality
			std::clog << "\nWeight computation for " << q_on_focal_points.fod->to_string(current_set) << ":\n-1\t <- " << q_on_focal_points.fod->to_string(current_set) << std::endl;
			T weight = 1/A->value;

			// check if there is a focal set in supersets contained in every focal set in supersets...

			// get supersets of current_set
			const std::vector<set_N_value<T>* >& supersets = q_on_focal_points
															.strict_supersets_of(
																	current_set
															);
			boost::dynamic_bitset<> I = supersets[0]->set;
			const boost::dynamic_bitset<> emptyset(q_on_focal_points.fod->size(), 0);
			size_t i = 1;
			for (; i < supersets.size(); ++i) {
				if (I != emptyset){
					I = I & supersets[i]->set;
				}else{
					break;
				}
			}
			for (i = 0; i < supersets.size(); ++i) {
				if (supersets[i]->set == I){
					std::clog << "1\t <- " << q_on_focal_points.fod->to_string(supersets[i]->set) << std::endl;
					return weight*supersets[i]->value;
				}
			}

			// ...if there is none, execute the general procedure
			powerset_btree<long int> polarities(*q_on_focal_points.fod, q_on_focal_points.block_size);
			long int polarity;

			for (i = 0; i < supersets.size(); ++i) {
				// By definition of the conjunctive decomposition, the polarity for this
				// superset in the computation of this current weight is :
				// ->  1 	if (sc - c) is odd
				// -> -1 	if (sc - c) is even
				if(((supersets[i]->set.count() - current_set.count()) & 1) > 0){	// check the first bit for parity check
					polarity = 1;
				}else{
					polarity = -1;
				}
				polarities.insert(supersets[i]->set, polarity);
			}

			std::unordered_map<size_t, std::vector<set_N_value<long int>* > > polarities_by_cardinality = polarities
																								.strict_supersets_of_by_cardinality(
																									current_set
																								);
			const std::vector<std::vector<set_N_value<long int>* >* >& ordered_vector = polarities
								.get_vector_of_vectors_ordered_by_cardinality(polarities_by_cardinality);

			long int s;
			for (size_t c = 0; c < ordered_vector.size(); ++c) {
				for (size_t i = 0; i < ordered_vector[c]->size(); ++i) {
					s = -1;
					const std::vector<set_N_value<long int>* >& polarities_supersets = polarities
																	.strict_subsets_of(
																			(*ordered_vector[c])[i]->set
																	);
					for (size_t ii = 0; ii < polarities_supersets.size(); ++ii) {
						s += polarities_supersets[ii]->value;
					}
					set_N_value<long int>* inserted_polarity = polarities.insert((*ordered_vector[c])[i]->set, -s);
					std::clog << to_string<long int>(*inserted_polarity, *q_on_focal_points.fod) << std::endl;

					weight *= pow(
								q_on_focal_points[(*ordered_vector[c])[i]->set]->value,
								-s
							);
				}
			}
			return weight;
		}
	};

} // namespace ow_bft

#endif // OW_BFT_CANONICAL_DECOMPOSITION_HPP
