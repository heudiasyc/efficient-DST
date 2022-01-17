#ifndef EFFICIENT_DST_COMPUTATION_SCHEME_HPP
#define EFFICIENT_DST_COMPUTATION_SCHEME_HPP

#include "macros.hpp"

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <vector>
#include <unordered_set>
#include <iostream>
#include <iomanip>
#include <utility>
#include <functional>

#include <fod.hpp>
#include <powerset_btree.hpp>
#include <powerset_function.hpp>


namespace efficient_DST{

	enum class transform_type_t: bool { zeta, Mobius };
	enum class order_relation_t: bool { subset, superset };
	enum class operation_t: bool { addition, multiplication };
	enum class irreducible_t: bool { join, meet };
	enum class scheme_type_t: int8_t { direct, consonant, semilattice, lattice };

	template <typename T, size_t N>
	class computation_scheme {
	public:

		static void autoset_and_build(
				const powerset_btree<T, N>& support,
				powerset_btree<T, N>& focal_points_tree,
				const order_relation_t& order_relation,
				const operation_t& transform_operation,
				std::vector<std::bitset<N> >& iota_sequence,
				scheme_type_t& scheme_type
		){
			// check if original_structure is almost Bayesian
			DEBUG(std::clog << "Linear analysis:" << std::endl;);
			const bool& is_almost_bayesian = linear_analysis_of_support(
					support,
					focal_points_tree,
					order_relation,
					neutral_value
			);

			if(is_almost_bayesian){
				DEBUG(std::clog << "almost Bayesian." << std::endl;);
				scheme_type = scheme_type_t::direct;
			}else{
				DEBUG(std::clog << "Consonance check:" << std::endl;);
				const bool& is_consonant = consonance_check(support);

				if(is_consonant){
					DEBUG(std::clog << "consonant." << std::endl;);
					focal_points_tree.copy(support);
					scheme_type = scheme_type_t::consonant;
				}else{
					DEBUG(std::clog << "not consonant." << std::endl;);

					if(support.size() < 3 * N){
						DEBUG({
							std::clog << "Number of focal sets equivalent to |FOD|." << std::endl;
							std::clog << "Transform to semilattice:" << std::endl;
						});

						build_semilattice_support(
								support,
								focal_points_tree,
								order_relation,
								neutral_value
						);

						if(focal_points_tree.size() < 3 * N){
							DEBUG(std::clog << "Number of focal points also equivalent to |FOD|." << std::endl;);
							scheme_type = scheme_type_t::direct;
						}else{
							DEBUG(std::clog << "Number of focal points superior to |FOD|." << std::endl;);

							set_semilattice_computation_scheme(
									support,
									order_relation,
									iota_sequence
									//iota_fiber_sequence
							);
							scheme_type = scheme_type_t::semilattice;
						}
					}else{
						DEBUG({
							std::clog << "Number of focal sets superior to |FOD|." << std::endl;
							std::clog << "Transform to lattice:" << std::endl;
						});

						if (order_relation == order_relation_t::subset){
							build_lattice_support_upper_closure(
									support,
									focal_points_tree,
									iota_sequence,
									neutral_value
							);
						}else{
							build_lattice_support_lower_closure(
									support,
									focal_points_tree,
									iota_sequence,
									neutral_value
							);
						}
						scheme_type = scheme_type_t::lattice;
					}
				}
			}
		}


		static void set_semilattice_computation_scheme(
				const powerset_btree<T, N>& support,
				const order_relation_t& order_relation,
				std::vector<std::bitset<N> >& iota_sequence
				//std::vector<size_t>& iota_fiber_sequence
		){
			if(order_relation == order_relation_t::subset){
				compute_iota_elements(irreducible_t::meet, support, iota_sequence);//, iota_fiber_sequence);
				//for (size_t o = 0; o < iota_sequence.size(); ++o){
				//	iota_sequence[o].flip();
				//}
			}else{
				compute_iota_elements(irreducible_t::join, support, iota_sequence);//, iota_fiber_sequence);
			}
		}


		static void execute(
				powerset_btree<T, N>& focal_points_tree,
				const transform_type_t& transform_type,
				const order_relation_t& order_relation,
				const operation_t& transform_operation,
				const std::vector<std::bitset<N> >& iota_sequence,
				//const std::vector<size_t>& iota_fiber_sequence,
				const scheme_type_t& scheme_type
		) {
			T neutral_value;
			std::function<T(const T&, const T&)> range_binary_operator;

			if (transform_operation == operation_t::addition){
				if (transform_type == transform_type_t::zeta)
					range_binary_operator = addition;
				else
					range_binary_operator = subtraction;
				neutral_value = 0;
			} else{
				if (transform_type == transform_type_t::zeta)
					range_binary_operator = multiplication;
				else
					range_binary_operator = division;
				neutral_value = 1;
			}

			switch (scheme_type) {
				case scheme_type_t::direct:
					execute_direct_transformation(
							focal_points_tree,
							transform_type,
							order_relation,
							range_binary_operator
					);
					break;

				case scheme_type_t::consonant:
					execute_consonant_transformation(
							focal_points_tree,
							transform_type,
							order_relation,
							range_binary_operator,
							neutral_value);
					break;

				case scheme_type_t::semilattice:
					if (order_relation == order_relation_t::subset){
						execute_EMT_with_upper_semilattice(
								focal_points_tree,
								transform_type,
								range_binary_operator,
								iota_sequence);
					}else{
						execute_EMT_with_lower_semilattice(
								focal_points_tree,
								transform_type,
								range_binary_operator,
								iota_sequence);
					}
					break;

				case scheme_type_t::lattice:
					if(order_relation == order_relation_t::subset){
						execute_EMT_with_lattice_upper_closure(
									focal_points_tree,
									transform_type,
									range_binary_operator,
									iota_sequence);
					}else{
						execute_EMT_with_lattice_lower_closure(
									focal_points_tree,
									transform_type,
									range_binary_operator,
									iota_sequence);
					}

					break;

				default:
					break;
			}
		}


		static void build_and_execute(
				const powerset_btree<T, N>& support,
				powerset_btree<T, N>& focal_points_tree,
				const transform_type_t& transform_type,
				const order_relation_t& order_relation,
				const operation_t& transform_operation,
				std::vector<std::bitset<N> >& iota_sequence,
				//std::vector<size_t>& iota_fiber_sequence,
				scheme_type_t& scheme_type
		) {
			switch(scheme_type){
				case scheme_type_t::semilattice:
					build_and_execute_EMT_with_semilattice(
							support,
							focal_points_tree,
							transform_type_t::zeta,
							order_relation,
							transform_operation,
							iota_sequence
							//iota_fiber_sequence
					);
					break;
				case scheme_type_t::lattice:
					build_and_execute_EMT_with_lattice(
							support,
							focal_points_tree,
							transform_type_t::zeta,
							order_relation,
							transform_operation,
							iota_sequence
					);
					break;
				default:
					build_and_execute_direct_transformation(
							support,
							focal_points_tree,
							transform_type_t::zeta,
							order_relation,
							transform_operation,
							scheme_type
					);
					break;
			}
		}


		static powerset_btree<T, N> direct_transformation(
				const powerset_btree<T, N>& support,
				const transform_type_t& transform_type,
				const order_relation_t& order_relation,
				const operation_t& transform_operation
		) {
			powerset_btree<T, N> focal_points_tree(support.get_FOD(), 2 * support.size());
			scheme_type_t scheme_type;
			build_and_execute_direct_transformation(
					support,
					focal_points_tree,
					transform_type,
					order_relation,
					transform_operation,
					scheme_type
			);
			return focal_points_tree;
		}


		static powerset_btree<T, N> EMT_with_semilattice(
				const powerset_btree<T, N>& support,
				const transform_type_t& transform_type,
				const order_relation_t& order_relation,
				const operation_t& transform_operation
		) {
			powerset_btree<T, N> focal_points_tree(support.get_FOD(), N * support.size());
			std::vector<std::bitset<N> > iota_sequence;
			build_and_execute_EMT_with_semilattice(
					support,
					focal_points_tree,
					transform_type,
					order_relation,
					transform_operation,
					iota_sequence
					//iota_fiber_sequence
			);
			return focal_points_tree;
		}


		/*
		 * Here, if transform_type = Mobius, then support must contain the image of the lattice support.
		 */
		static powerset_btree<T, N> EMT_with_lattice(
				const powerset_btree<T, N>& support,
				const transform_type_t& transform_type,
				const order_relation_t& order_relation,
				const operation_t& transform_operation
		) {
			powerset_btree<T, N> cropped_lattice_support(support.get_FOD(), N * support.size());
			std::vector<std::bitset<N> > iota_sequence;
			build_and_execute_EMT_with_lattice(
					support,
					cropped_lattice_support,
					transform_type,
					order_relation,
					transform_operation,
					iota_sequence
					//iota_fiber_sequence
			);
			return cropped_lattice_support;
		}


		static powerset_btree<T, N> EMT_with_lattice_Mobius_from_zeta_values(
				const std::vector<T>& vec,
				FOD<N>& fod,
				const transform_type_t& transform_type,
				const order_relation_t& order_relation,
				const operation_t& transform_operation
		) {
			powerset_btree<T, N> cropped_lattice_support(&fod, vec.size());
			extract_semilattice_support(vec, cropped_lattice_support, order_relation);
			std::vector<std::bitset<N> > iota_sequence;
			std::function<T(const T&, const T&)> range_binary_operator;

			if (transform_operation == operation_t::addition){
				if (transform_type == transform_type_t::zeta)
					range_binary_operator = addition;
				else
					range_binary_operator = subtraction;
			} else{
				if (transform_type == transform_type_t::zeta)
					range_binary_operator = multiplication;
				else
					range_binary_operator = division;
			}

			if (order_relation == order_relation_t::subset){
				build_lattice_support_upper_closure(
						vec,
						cropped_lattice_support,
						iota_sequence
				);
				execute_EMT_with_lattice_upper_closure(
						cropped_lattice_support,
						transform_type,
						range_binary_operator,
						iota_sequence
				);
			}else{
				build_lattice_support_lower_closure(
						vec,
						cropped_lattice_support,
						iota_sequence
				);
				execute_EMT_with_lattice_lower_closure(
						cropped_lattice_support,
						transform_type,
						range_binary_operator,
						iota_sequence
				);
			}
			return cropped_lattice_support;
		}

	protected:

		static void build_and_execute_direct_transformation(
				const powerset_btree<T, N>& support,
				powerset_btree<T, N>& focal_points_tree,
				const transform_type_t& transform_type,
				const order_relation_t& order_relation,
				const operation_t& transform_operation,
				scheme_type_t& scheme_type
		) {
			T neutral_value;
			std::function<T(const T&, const T&)> range_binary_operator;

			if (transform_operation == operation_t::addition){
				if (transform_type == transform_type_t::zeta)
					range_binary_operator = addition;
				else
					range_binary_operator = subtraction;
				neutral_value = 0;
			} else{
				if (transform_type == transform_type_t::zeta)
					range_binary_operator = multiplication;
				else
					range_binary_operator = division;
				neutral_value = 1;
			}

			// check if original_structure is almost Bayesian
			const bool& is_almost_bayesian = try_linear_focal_points_computation<
				T,
				N,
				up_inclusion<T, N>,
				value_inplace_operation
			>(
				support,
				focal_points_tree
			);

			scheme_type = scheme_type_t::direct;

			if (!is_almost_bayesian){
				DEBUG(std::clog << "Consonance check:" << std::endl;);
				const bool& is_consonant = consonance_check(support);

				if(is_consonant){
					focal_points_tree.copy(support);
					scheme_type = scheme_type_t::consonant;
				}else{
					build_semilattice_support(
							support,
							focal_points_tree,
							order_relation,
							neutral_value
					);
				}
			}

			if(scheme_type == scheme_type_t::direct){
				execute_direct_transformation(
						focal_points_tree,
						transform_type,
						order_relation,
						range_binary_operator);
			}else{
				execute_consonant_transformation(
						focal_points_tree,
						transform_type,
						order_relation,
						range_binary_operator,
						neutral_value);
			}
		}


		static void build_and_execute_EMT_with_semilattice(
				const powerset_btree<T, N>& support,
				powerset_btree<T, N>& focal_points_tree,
				const transform_type_t& transform_type,
				const order_relation_t& order_relation,
				const operation_t& transform_operation,
				std::vector<std::bitset<N> >& iota_sequence
		) {
			T neutral_value;
			std::function<T(const T&, const T&)> range_binary_operator;

			if (transform_operation == operation_t::addition){
				if (transform_type == transform_type_t::zeta)
					range_binary_operator = addition;
				else
					range_binary_operator = subtraction;
				neutral_value = 0;
			} else{
				if (transform_type == transform_type_t::zeta)
					range_binary_operator = multiplication;
				else
					range_binary_operator = division;
				neutral_value = 1;
			}

			if(order_relation == order_relation_t::subset){
				compute_iota_elements(irreducible_t::meet, support, iota_sequence);
			}else{
				compute_iota_elements(irreducible_t::join, support, iota_sequence); //iota_fiber_sequence);
			}

			build_semilattice_support(
					support,
					focal_points_tree,
					order_relation,
					neutral_value
			);

			if (order_relation == order_relation_t::subset){
				execute_EMT_with_upper_semilattice(
						focal_points_tree,
						transform_type,
						range_binary_operator,
						iota_sequence
				);
			}else{
				execute_EMT_with_lower_semilattice(
						focal_points_tree,
						transform_type,
						range_binary_operator,
						iota_sequence
				);
			}
		}


		static void build_and_execute_EMT_with_lattice(
				const powerset_btree<T, N>& support,
				powerset_btree<T, N>& cropped_lattice_support,
				const transform_type_t& transform_type,
				const order_relation_t& order_relation,
				const operation_t& transform_operation,
				std::vector<std::bitset<N> >& iota_sequence
		) {
			T neutral_value;
			std::function<T(const T&, const T&)> range_binary_operator;

			if (transform_operation == operation_t::addition){
				if (transform_type == transform_type_t::zeta)
					range_binary_operator = addition;
				else
					range_binary_operator = subtraction;
				neutral_value = 0;
			} else{
				if (transform_type == transform_type_t::zeta)
					range_binary_operator = multiplication;
				else
					range_binary_operator = division;
				neutral_value = 1;
			}

			if (order_relation == order_relation_t::subset){
				build_lattice_support_upper_closure(
						support,
						cropped_lattice_support,
						iota_sequence,
						neutral_value
				);
				execute_EMT_with_lattice_upper_closure(
						cropped_lattice_support,
						transform_type,
						range_binary_operator,
						iota_sequence
				);
			}else{
				build_lattice_support_lower_closure(
						support,
						cropped_lattice_support,
						iota_sequence,
						neutral_value
				);
				execute_EMT_with_lattice_lower_closure(
						cropped_lattice_support,
						transform_type,
						range_binary_operator,
						iota_sequence
				);
			}
		}


		static bool consonance_check(
				const powerset_btree<T, N>& support
		){
			// Cardinality map of focal sets in original_structure
			std::unordered_map<size_t, std::vector<set_N_value<T, N>* > > support_card_map = support.elements_by_set_cardinality();
			std::vector<size_t> support_ordered_cardinalities;
			support.get_FOD()->sort_cardinalities(support_ordered_cardinalities, support_card_map, order_t::ascending);

			// if at least one focal point has already been discovered
			// OR if there are at least two elements with the same 0-indexed cardinality (i.e. the least one) in structure,
			// then structure cannot be consonant
			if(support_card_map[support_ordered_cardinalities[0]].size() > 1)
				return false;

			for (size_t i = 1; i < support_ordered_cardinalities.size(); ++i) {
				// if there are at least two elements with the same i-indexed cardinality in structure
				// OR if the set of cardinality index i-1 is not a subset of the one of cardinality index i,
				// then structure cannot be consonant
				if (support_card_map[support_ordered_cardinalities[i]].size() > 1 )
					return false;
				if(!FOD<N>::is_subset_of(
						support_card_map[support_ordered_cardinalities[i-1]][0]->set,
						support_card_map[support_ordered_cardinalities[i]][0]->set
					  )
				){
					return false;
				}
			}
			return true;
		}
	};
}		// namespace efficient_DST

#endif // EFFICIENT_DST_COMPUTATION_SCHEME_HPP
