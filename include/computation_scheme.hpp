#ifndef EFFICIENT_DST_COMPUTATION_SCHEME_HPP
#define EFFICIENT_DST_COMPUTATION_SCHEME_HPP

#include "macros.hpp"

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <boost/dynamic_bitset.hpp>
#include <vector>
#include <unordered_map>
#include <unordered_set>
#include <boost/functional/hash.hpp>
#include <iostream>
#include <iomanip>

#include <fod.hpp>
#include <powerset_btree.hpp>
#include <powerset_function.hpp>


namespace efficient_DST{

	enum class transform_type_t: bool { zeta, Mobius };
	enum class order_relation_t: bool { subset, superset };
	enum class operation_t: bool { addition, multiplication };
	enum class version_t: bool { regular, dual };
	enum class scheme_type_t: int8_t { direct, consonant, semilattice, lattice };

	template <typename T = double>
	class computation_scheme {
	public:

		static void autoset_and_build(
				const powerset_btree<T>& support,
				powerset_btree<T>& focal_points_tree,
				powerset_btree<set_N_value<T>* >& focal_points_dual_tree,
				const order_relation_t& order_relation,
				const operation_t& transform_operation,
				std::vector<boost::dynamic_bitset<> >& iota_sequence,
				//std::vector<size_t>& iota_fiber_sequence,
				scheme_type_t& scheme_type
		){
			T neutral_value;
			if (transform_operation == operation_t::addition){
				neutral_value = 0;
			} else{
				neutral_value = 1;
			}
			powerset_btree<T>::flip_powerset_from_to(focal_points_tree, focal_points_dual_tree);

			// Cardinality map of focal sets in original_structure
			std::unordered_map<size_t, std::vector<set_N_value<T>* > > support_card_map = support.elements_by_set_cardinality();
			std::vector<set_N_value<T>* > elements_generated_from_support;
			elements_generated_from_support.reserve(2 * support.size());
			std::unordered_set<size_t> focal_points;
			focal_points.reserve(2 * support.size());

			// check if original_structure is almost Bayesian
			DEBUG(std::clog << "Linear analysis:" << std::endl;);
			const bool& is_almost_bayesian = linear_analysis_of_support(
					focal_points_tree,
					focal_points_dual_tree,
					support_card_map,
					elements_generated_from_support,
					order_relation,
					neutral_value
			);

			if(is_almost_bayesian){
				DEBUG(std::clog << "almost Bayesian." << std::endl;);
				scheme_type = scheme_type_t::direct;
			}else{
				// sort sets in support by increasing cardinality
				const std::vector<size_t>& support_ordered_cardinalities = powerset_btree<T>::get_sorted_cardinalities(support_card_map, *support.get_FOD());

				DEBUG(std::clog << "Consonance check:" << std::endl;);
				const bool& is_consonant = consonance_check(
					elements_generated_from_support,
					support_card_map,
					support_ordered_cardinalities
				);

				if(is_consonant){
					DEBUG(std::clog << "consonant." << std::endl;);
					scheme_type = scheme_type_t::consonant;
				}else{
					DEBUG(std::clog << "not consonant." << std::endl;);

					focal_points.reserve(2 * support.get_FOD_size() * support.size());

					if(support.size() < 3 * support.get_FOD_size()){
						DEBUG({
							std::clog << "Number of focal sets equivalent to |FOD|." << std::endl;
							std::clog << "Transform to semilattice:" << std::endl;
						});

						build_semilattice_support(
								focal_points_tree,
								focal_points_dual_tree,
								support_card_map,
								support_ordered_cardinalities,
								elements_generated_from_support,
								order_relation,
								neutral_value
						);

						if(focal_points_tree.size() < 3 * support.get_FOD_size()){
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
									focal_points_tree,
									iota_sequence,
									neutral_value
							);
						}else{
							build_lattice_support_lower_closure(
									focal_points_tree,
									focal_points_dual_tree,
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
				const powerset_btree<T>& support,
				const order_relation_t& order_relation,
				std::vector<boost::dynamic_bitset<> >& iota_sequence
				//std::vector<size_t>& iota_fiber_sequence
		){
			if(order_relation == order_relation_t::subset){
				compute_iota_elements(version_t::dual, support, iota_sequence);//, iota_fiber_sequence);
				for (size_t o = 0; o < iota_sequence.size(); ++o){
					iota_sequence[o].flip();
				}
			}else{
				compute_iota_elements(version_t::regular, support, iota_sequence);//, iota_fiber_sequence);
			}
		}


		static void extract_semilattice_support(
			const std::vector<T>& powerset_values,
			powerset_btree<T>& focal_points_tree,
			powerset_btree<set_N_value<T>* >& focal_points_dual_tree,
			const order_relation_t& order_relation
		) {
			size_t powerset_size = pow(2, focal_points_tree.get_FOD_size());
			if(powerset_values.size() != powerset_size){
				std::cerr << "\nThe given vector does not feature the same size as the powerset of the given FOD.\n";
				return;
			}
			std::unordered_map<T, std::vector<boost::dynamic_bitset<> > > focal_points_map;

			for (size_t n = 0; n < powerset_size; ++n){
				boost::dynamic_bitset<> set(focal_points_tree.get_FOD_size(), n);

				if (focal_points_map.find(powerset_values[n]) == focal_points_map.end()){
					focal_points_map.emplace(powerset_values[n], (std::vector<boost::dynamic_bitset<> >) {set});
				}else{
					bool new_point = true;
					std::vector<boost::dynamic_bitset<> >& focal_points = focal_points_map[powerset_values[n]];
					for (size_t i = 0; i < focal_points.size(); ++i){
						if (order_relation == order_relation_t::subset){
							if (FOD::is_subset_of(set, focal_points[i])){
								focal_points[i] = set;
								new_point = false;
							}else if (FOD::is_subset_of(focal_points[i], set)){
								new_point = false;
								break;
							}
						}else{
							if (FOD::is_superset_of(set, focal_points[i])){
								focal_points[i] = set;
								new_point = false;
							}else if (FOD::is_superset_of(focal_points[i], set)){
								new_point = false;
								break;
							}
						}
					}
					if (new_point){
						focal_points.emplace_back(set);
					}
				}
			}
			for (auto kv : focal_points_map){
				for (size_t i = 0; i < kv.second.size(); ++i){
					if(!focal_points_tree[kv.second[i]]){
						insert_focal_point(focal_points_tree, focal_points_dual_tree, kv.second[i], kv.first);
					}
				}
			}
		}


		static void execute(
				powerset_btree<T>& focal_points_tree,
				powerset_btree<set_N_value<T>* >& focal_points_dual_tree,
				const transform_type_t& transform_type,
				const order_relation_t& order_relation,
				const operation_t& transform_operation,
				const std::vector<boost::dynamic_bitset<> >& iota_sequence,
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
				const powerset_btree<T>& support,
				powerset_btree<T>& focal_points_tree,
				powerset_btree<set_N_value<T>* >& focal_points_dual_tree,
				const transform_type_t& transform_type,
				const order_relation_t& order_relation,
				const operation_t& transform_operation,
				std::vector<boost::dynamic_bitset<> >& iota_sequence,
				//std::vector<size_t>& iota_fiber_sequence,
				scheme_type_t& scheme_type
		) {
			switch(scheme_type){
				case scheme_type_t::semilattice:
					build_and_execute_EMT_with_semilattice(
							support,
							focal_points_tree,
							focal_points_dual_tree,
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
							focal_points_dual_tree,
							transform_type_t::zeta,
							order_relation,
							transform_operation,
							iota_sequence
							//iota_fiber_sequence
					);
					break;
				default:
					build_and_execute_direct_transformation(
							support,
							focal_points_tree,
							focal_points_dual_tree,
							transform_type_t::zeta,
							order_relation,
							transform_operation,
							scheme_type
					);
					break;
			}
		}


		static powerset_btree<T> direct_transformation(
				const powerset_btree<T>& support,
				const transform_type_t& transform_type,
				const order_relation_t& order_relation,
				const operation_t& transform_operation
		) {
			powerset_btree<T> focal_points_tree(support);
			powerset_btree<set_N_value<T>* > focal_points_dual_tree(support.get_FOD(), support.get_block_size());
			scheme_type_t scheme_type;
			build_and_execute_direct_transformation(
					support,
					focal_points_tree,
					focal_points_dual_tree,
					transform_type,
					order_relation,
					transform_operation,
					scheme_type
			);
			return focal_points_tree;
		}


		static powerset_btree<T> EMT_with_semilattice(
				const powerset_btree<T>& support,
				const transform_type_t& transform_type,
				const order_relation_t& order_relation,
				const operation_t& transform_operation
		) {
			powerset_btree<T> focal_points_tree(support);
			powerset_btree<set_N_value<T>* > focal_points_dual_tree(support.get_FOD(), support.get_block_size());
			std::vector<boost::dynamic_bitset<> > iota_sequence;
			//std::vector<size_t> iota_fiber_sequence;
			build_and_execute_EMT_with_semilattice(
					support,
					focal_points_tree,
					focal_points_dual_tree,
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
		static powerset_btree<T> EMT_with_lattice(
				const powerset_btree<T>& support,
				const transform_type_t& transform_type,
				const order_relation_t& order_relation,
				const operation_t& transform_operation
		) {
			powerset_btree<T> cropped_lattice_support(support);
			powerset_btree<set_N_value<T>* > cropped_lattice_support_dual(support.get_FOD(), support.get_block_size());
			std::vector<boost::dynamic_bitset<> > iota_sequence;
			//std::vector<size_t> iota_fiber_sequence;
			build_and_execute_EMT_with_lattice(
					support,
					cropped_lattice_support,
					cropped_lattice_support_dual,
					transform_type,
					order_relation,
					transform_operation,
					iota_sequence
					//iota_fiber_sequence
			);
			return cropped_lattice_support;
		}


		static powerset_btree<T> EMT_with_lattice_Mobius_from_zeta_values(
				const std::vector<T>& vec,
				FOD& fod,
				const transform_type_t& transform_type,
				const order_relation_t& order_relation,
				const operation_t& transform_operation
		) {
			powerset_btree<T> focal_points(&fod, powerset_function<T>::block_size);
			powerset_btree<set_N_value<T>* > focal_points_dual(&fod, powerset_function<T>::block_size);
			extract_semilattice_support(vec, focal_points, focal_points_dual, order_relation);
			powerset_btree<T> cropped_lattice_support(focal_points);
			powerset_btree<set_N_value<T>* > cropped_lattice_support_dual(&fod, powerset_function<T>::block_size);
			std::vector<boost::dynamic_bitset<> > iota_sequence;
			//std::vector<size_t> iota_fiber_sequence;
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
						cropped_lattice_support,
						iota_sequence,
						vec
				);
				execute_EMT_with_lattice_upper_closure(
						cropped_lattice_support,
						transform_type,
						range_binary_operator,
						iota_sequence
				);
			}else{
				powerset_btree<T>::flip_powerset_from_to(cropped_lattice_support, cropped_lattice_support_dual);
				build_lattice_support_lower_closure(
						cropped_lattice_support,
						cropped_lattice_support_dual,
						iota_sequence,
						vec
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


		static std::vector<T> FMT(
				const std::vector<T>& powerset_values,
				const size_t& N,
				const transform_type_t& transform_type,
				const order_relation_t& order_relation,
				const operation_t& transform_operation
		) {
			if(powerset_values.size() != pow(2, N)){
				std::cerr << "\nThe size of the given vector is not 2^N, where N is the given size corresponding to the considered FOD.\n";
				return powerset_values;
			}
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

			std::vector<T> transform(powerset_values);
			size_t sub_powerset_size, sub_powerset_dual_size, index;
			for (size_t i = 1; i <= N; ++i){

				sub_powerset_size = pow(2, i);
				for (size_t j = 1; j <= sub_powerset_size; j += 2){

					sub_powerset_dual_size = pow(2, N - i);
					for (size_t k = 0; k <= sub_powerset_dual_size-1; ++k){

						if (order_relation == order_relation_t::subset){
							index = j * sub_powerset_dual_size + k;
							transform[index] = range_binary_operator(transform[index], transform[(j-1) * sub_powerset_dual_size + k]);
						}else{
							index = (j-1) * sub_powerset_dual_size + k;
							transform[index] = range_binary_operator(transform[index], transform[j * sub_powerset_dual_size + k]);
						}
					}
				}
			}
			return transform;
		}

	protected:

		static void execute_direct_transformation(
				powerset_btree<T>& focal_points_tree,
				const transform_type_t& transform_type,
				const order_relation_t& order_relation,
				std::function<T(const T&, const T&)> range_binary_operator
		) {
			if (transform_type == transform_type_t::zeta){
				powerset_btree<T> focal_points_N_initial_values(focal_points_tree);
				T val;
				const std::vector<set_N_value<T>* >& focal_points = focal_points_N_initial_values.elements();
				std::vector<set_N_value<T>* > elements;

				for (size_t i = 0; i < focal_points.size(); ++i) {
					val = focal_points[i]->value;
					if (order_relation == order_relation_t::subset){
						elements = focal_points_N_initial_values.strict_subsets_of(focal_points[i]->set);
					} else{
						elements = focal_points_N_initial_values.strict_supersets_of(focal_points[i]->set);
					}
					for (size_t ii = 0; ii < elements.size(); ++ii) {
						val = range_binary_operator(val, elements[ii]->value);
					}
					focal_points_tree.insert(focal_points[i]->set, val);
				}
			}else{
				T val;
				std::unordered_map<size_t, std::vector<set_N_value<T>* > > focal_points_card_map = focal_points_tree.elements_by_set_cardinality();
				std::vector<size_t> focal_points_ordered_cardinalities = powerset_btree<T>::get_sorted_cardinalities(focal_points_card_map, *focal_points_tree.get_FOD());
				std::vector<set_N_value<T>* > elements;

				if (order_relation == order_relation_t::superset) {
					std::reverse(focal_points_ordered_cardinalities.begin(), focal_points_ordered_cardinalities.end());
				}
				for (size_t c = 0; c < focal_points_ordered_cardinalities.size(); ++c) {
					const std::vector<set_N_value<T>* >& focal_points = focal_points_card_map[focal_points_ordered_cardinalities[c]];
					for (size_t i = 0; i < focal_points.size(); ++i) {
						val = focal_points[i]->value;
						if (order_relation == order_relation_t::subset){
							elements = focal_points_tree.strict_subsets_of(focal_points[i]->set);
						} else{
							elements = focal_points_tree.strict_supersets_of(focal_points[i]->set);
						}
						for (size_t ii = 0; ii < elements.size(); ++ii) {
							val = range_binary_operator(val, elements[ii]->value);
						}
						focal_points_tree.insert(focal_points[i]->set, val);
					}
				}
			}
		}


		static void execute_consonant_transformation(
				powerset_btree<T>& support,
				const transform_type_t& transform_type,
				const order_relation_t& order_relation,
				std::function<T(const T&, const T&)> range_binary_operator,
				const T& neutral_value
		) {
			std::unordered_map<size_t, std::vector<set_N_value<T>* > > support_card_map = support.elements_by_set_cardinality();
			std::vector<size_t> support_ordered_cardinalities = powerset_btree<T>::get_sorted_cardinalities(support_card_map, *support.get_FOD());
			T value, preceding_value = neutral_value;

			if (order_relation == order_relation_t::superset) {
				std::reverse(support_ordered_cardinalities.begin(), support_ordered_cardinalities.end());
			}
			for (size_t i = 0; i < support_ordered_cardinalities.size(); ++i) {
				set_N_value<T>* set_value = support_card_map[support_ordered_cardinalities[i]][0];
				value = range_binary_operator(set_value->value, preceding_value);
				if (transform_type == transform_type_t::zeta){
					set_value->value = value;
					preceding_value = set_value->value;
				}else{
					preceding_value = set_value->value;
					set_value->value = value;
				}
			}
		}


		static void build_and_execute_direct_transformation(
				const powerset_btree<T>& support,
				powerset_btree<T>& focal_points_tree,
				powerset_btree<set_N_value<T>* >& focal_points_dual_tree,
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

			powerset_btree<T>::flip_powerset_from_to(focal_points_tree, focal_points_dual_tree);

			// Cardinality map of focal sets in original_structure
			std::unordered_map<size_t, std::vector<set_N_value<T>* > > support_card_map = support.elements_by_set_cardinality();
			std::vector<set_N_value<T>* > elements_generated_from_support;

			// check if original_structure is almost Bayesian
			const bool& is_almost_bayesian = linear_analysis_of_support(
					focal_points_tree,
					focal_points_dual_tree,
					support_card_map,
					elements_generated_from_support,
					order_relation,
					neutral_value
			);

			scheme_type = scheme_type_t::direct;

			if(!is_almost_bayesian){
				// sort sets in support by increasing cardinality
				const std::vector<size_t>& support_ordered_cardinalities = powerset_btree<T>::get_sorted_cardinalities(support_card_map, *support.get_FOD());

				DEBUG(std::clog << "Consonance check:" << std::endl;);
				const bool& is_consonant = consonance_check(
					elements_generated_from_support,
					support_card_map,
					support_ordered_cardinalities
				);

				if(is_consonant){
					scheme_type = scheme_type_t::consonant;
				}else{
					build_semilattice_support(
							focal_points_tree,
							focal_points_dual_tree,
							support_card_map,
							support_ordered_cardinalities,
							elements_generated_from_support,
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


		static void execute_EMT_with_lower_semilattice(
				powerset_btree<T>& focal_points_tree,
				const transform_type_t& transform_type,
				std::function<T(const T&, const T&)> range_binary_operator,
				const std::vector<boost::dynamic_bitset<> >& iota_sequence
		) {
			if (iota_sequence.size() == 0)
				return;

			std::vector<boost::dynamic_bitset<> > sync_sequence;
			sync_sequence.reserve(iota_sequence.size());

			sync_sequence.emplace_back(iota_sequence[0]);
			for (size_t o = 1; o < iota_sequence.size(); ++o){
				sync_sequence.emplace_back(FOD::set_union(sync_sequence[o-1], iota_sequence[o]));
			}

			std::vector<set_N_value<T>* > focal_points = focal_points_tree.elements();
			const powerset_btree<set_N_value<T>* >& focal_points_superset_map = EMT_superset_map(focal_points_tree, iota_sequence);

			//boost::dynamic_bitset<> emptyset(focal_points_tree.get_FOD_size());
			//emptyset.set(1);

			size_t iota_index;
			for (size_t o = 0; o < iota_sequence.size(); ++o){
				if (transform_type == transform_type_t::zeta)
					iota_index = o;
				else
					iota_index = iota_sequence.size()-1 - o;

				for (size_t i = 0; i < focal_points.size(); ++i){
					const boost::dynamic_bitset<>& set_B = FOD::set_union(focal_points[i]->set, iota_sequence[iota_index]);

					if (focal_points[i]->set != set_B){
						const set_N_value<set_N_value<T>* >* X = focal_points_superset_map.smallest_superset_of(set_B);

						//if(focal_points[i]->set == emptyset){
							//std::cout << set_B << std::endl;
							//if(X)
							//	std::cout << " => " << X->value->value << "\t <- " << focal_points_tree.get_FOD()->to_string(X->set) << std::endl;
						//}

						if (X && FOD::is_or_is_subset_of(X->set, FOD::set_union(focal_points[i]->set, sync_sequence[iota_index]))){
							focal_points[i]->value = range_binary_operator(focal_points[i]->value, X->value->value);
						}
					}
				}
			}
		}


		static void execute_EMT_with_upper_semilattice(
				powerset_btree<T>& focal_points_tree,
				const transform_type_t& transform_type,
				std::function<T(const T&, const T&)> range_binary_operator,
				const std::vector<boost::dynamic_bitset<> >& iota_sequence
		) {
			if (iota_sequence.size() == 0)
				return;

			std::vector<boost::dynamic_bitset<> > sync_sequence;
			sync_sequence.reserve(iota_sequence.size());

			sync_sequence.emplace_back(iota_sequence[0]);
			for (size_t o = 1; o < iota_sequence.size(); ++o){
				sync_sequence.emplace_back(FOD::set_intersection(sync_sequence[o-1], iota_sequence[o]));
			}

			std::vector<set_N_value<T>* > focal_points = focal_points_tree.elements();
			const powerset_btree<set_N_value<T>* >& focal_points_subset_map = EMT_subset_map(focal_points_tree, iota_sequence);

			boost::dynamic_bitset<> fod(focal_points_tree.get_FOD_size());
			fod.set(1);
			fod.flip();

			size_t iota_index;
			std::clog << "HEY !" << std::endl;
			for (size_t o = 0; o < iota_sequence.size(); ++o){
				if (transform_type == transform_type_t::zeta)
					iota_index = o;
				else
					iota_index = iota_sequence.size()-1 - o;

				std::clog << "IOTA: " << iota_sequence[o] << std::endl;

				for (size_t i = 0; i < focal_points.size(); ++i){
					const boost::dynamic_bitset<>& set_B = FOD::set_intersection(focal_points[i]->set, iota_sequence[iota_index]);

					if (focal_points[i]->set != set_B){
						const set_N_value<set_N_value<T>* >* X = focal_points_subset_map.biggest_subset_of(set_B);

						if(focal_points[i]->set == fod){
							std::clog << focal_points[i]->set << " + " << iota_sequence[iota_index] << " = " << set_B << std::endl;
							if(X)
								std::clog << " => " << X->value->value << "\t <- " << focal_points_tree.get_FOD()->to_string(X->set) << std::endl;
						}

						if (X && FOD::is_or_is_superset_of(X->set, FOD::set_intersection(focal_points[i]->set, sync_sequence[iota_index]))){
							focal_points[i]->value = range_binary_operator(focal_points[i]->value, X->value->value);
						}
					}
				}
			}
		}

/*
		static void execute_EMT_with_semilattice(
				powerset_btree<T>& focal_points_tree,
				powerset_btree<set_N_value<T>* >& focal_points_dual_tree,
				const transform_type_t& transform_type,
				const order_relation_t& order_relation,
				std::function<T(const T&, const T&)> range_binary_operator,
				const std::vector<boost::dynamic_bitset<> >& iota_sequence
				//const std::vector<size_t>& iota_fiber_sequence
		) {
			if (order_relation == order_relation_t::subset){
				execute_EMT_with_upper_semilattice(
						focal_points_dual_tree,
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
		}
*/

		static void build_and_execute_EMT_with_semilattice(
				const powerset_btree<T>& support,
				powerset_btree<T>& focal_points_tree,
				powerset_btree<set_N_value<T>* >& focal_points_dual_tree,
				const transform_type_t& transform_type,
				const order_relation_t& order_relation,
				const operation_t& transform_operation,
				std::vector<boost::dynamic_bitset<> >& iota_sequence
				//std::vector<size_t>& iota_fiber_sequence
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

			powerset_btree<T>::flip_powerset_from_to(focal_points_tree, focal_points_dual_tree);

			if(order_relation == order_relation_t::subset){
				compute_iota_elements(version_t::dual, support, iota_sequence); //iota_fiber_sequence);
				//for (size_t o = 0; o < iota_sequence.size(); ++o){
				//	iota_sequence[o].flip();
				//}
			}else{
				compute_iota_elements(version_t::regular, support, iota_sequence); //iota_fiber_sequence);
			}

			// Cardinality map of focal sets in original_structure
			std::unordered_map<size_t, std::vector<set_N_value<T>* > > support_card_map = support.elements_by_set_cardinality();
			std::vector<set_N_value<T>* > elements_generated_from_support;
			// sort sets in support by increasing cardinality
			const std::vector<size_t>& support_ordered_cardinalities = powerset_btree<T>::get_sorted_cardinalities(support_card_map, *support.get_FOD());

			build_semilattice_support(
					focal_points_tree,
					focal_points_dual_tree,
					support_card_map,
					support_ordered_cardinalities,
					elements_generated_from_support,
					order_relation,
					neutral_value
			);

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
		}


		static powerset_btree<set_N_value<T>* > EMT_superset_map(
				const powerset_btree<set_N_value<T>* >& focal_points_dual_tree,
				const std::vector<boost::dynamic_bitset<> >& iota_sequence_dual
		) {
			powerset_btree<set_N_value<T>* > superset_map(focal_points_dual_tree.get_FOD(), focal_points_dual_tree.get_block_size());
			const std::vector<set_N_value<set_N_value<T>* >* >& elements = focal_points_dual_tree.elements();
			std::unordered_set<size_t> to_be_inserted;
			to_be_inserted.reserve(focal_points_dual_tree.get_FOD_size() * focal_points_dual_tree.size());

			for (size_t e = 0; e < elements.size(); ++e) {
				superset_map.insert(elements[e]->set, elements[e]->value);
				for (size_t i = 0; i < iota_sequence_dual.size(); ++i) {
					//superset_map.insert_null_node(FOD::set_union(elements[e]->set, iota_sequence_dual[i]));
					to_be_inserted.insert(FOD::set_union(elements[e]->set, iota_sequence_dual[i]).to_ulong());
				}
			}
			// between focal_points_dual_tree.size() and |lattice support| elements in to_be_inserted
			for (const auto& new_set: to_be_inserted){
				superset_map.insert_null_node(boost::dynamic_bitset<>(focal_points_dual_tree.get_FOD_size(), new_set));
			}

			build_EMT_superset_map(superset_map);
			return superset_map;
		}


		static powerset_btree<set_N_value<T>* > EMT_superset_map(
				const powerset_btree<T>& focal_points_tree,
				const std::vector<boost::dynamic_bitset<> >& iota_sequence
		) {
			powerset_btree<set_N_value<T>* > superset_map(focal_points_tree.get_FOD(), focal_points_tree.get_block_size());
			const std::vector<set_N_value<T>* >& elements = focal_points_tree.elements();

			for (size_t e = 0; e < elements.size(); ++e) {
				superset_map.insert(elements[e]->set, elements[e]);
				for (size_t i = 0; i < iota_sequence.size(); ++i) {
					superset_map.insert_null_node(FOD::set_union(elements[e]->set, iota_sequence[i]));
				}
			}
			build_EMT_superset_map(superset_map);
			return superset_map;
		}


		static void build_EMT_superset_map(powerset_btree<set_N_value<T>* >& superset_map){
			std::unordered_map<size_t, std::vector<set_N_value<set_N_value<T>* >* > > card_map = superset_map.elements_by_set_cardinality(true);
			std::vector<size_t> ordered_cardinalities = powerset_btree<set_N_value<T>* >::get_sorted_cardinalities(card_map, *superset_map.get_FOD());
			// powerset_btree of nodes without right child in this powerset.
			// Each node in terminal_connection_tree is associated with the address of the corresponding node in this powerset.
			powerset_btree<set_N_value<set_N_value<T>* >* > terminal_connection_tree(superset_map.get_FOD(), superset_map.get_block_size());

			size_t c = 0;
			if (ordered_cardinalities[0] == 0){
				c = 1;
			}

			for (; c < ordered_cardinalities.size(); ++c){
				const std::vector<set_N_value<set_N_value<T>* >* >& elements = card_map[ordered_cardinalities[c]];

				for (size_t i = 0; i < elements.size(); ++i){

					if(!elements[i]->is_null){
						const std::vector<set_N_value<set_N_value<set_N_value<T>* >* >* >& subsets = terminal_connection_tree.subsets_of(elements[i]->set);
						node<set_N_value<T>* >* node_i = (efficient_DST::node<set_N_value<T>* >*) elements[i];
						for (size_t s = 0; s < subsets.size(); ++s){
							node<set_N_value<T>* >* subset = (efficient_DST::node<set_N_value<T>* >*) subsets[s]->value;

							if(subset->is_null || node_i->depth > subset->depth){
								subset->right = node_i;
							}
						}
						for (size_t s = 0; s < subsets.size(); ++s){
							node<set_N_value<T>* >* subset = (efficient_DST::node<set_N_value<T>* >*) subsets[s]->value;
							if(subset->is_null || node_i->depth > subset->depth){
								terminal_connection_tree.nullify(subsets[s]);
							}
						}
					}
				}
				for (size_t i = 0; i < elements.size(); ++i){
					node<set_N_value<T>* >* node = (efficient_DST::node<set_N_value<T>* >*) elements[i];
					if (node->is_null && !node->right){
						terminal_connection_tree.insert(elements[i]->set, elements[i]);
					}
				}
			}
		}


		static powerset_btree<set_N_value<T>* > EMT_subset_map(
				const powerset_btree<T>& focal_points_tree,
				const std::vector<boost::dynamic_bitset<> >& iota_sequence
		) {
			powerset_btree<set_N_value<T>* > subset_map(focal_points_tree.get_FOD(), focal_points_tree.get_block_size());
			const std::vector<set_N_value<T>* >& elements = focal_points_tree.elements();

			for (size_t e = 0; e < elements.size(); ++e) {
				subset_map.insert(elements[e]->set, elements[e]);
				for (size_t i = 0; i < iota_sequence.size(); ++i) {
					subset_map.insert_null_node(FOD::set_intersection(elements[e]->set, iota_sequence[i]));
				}
			}
			build_EMT_subset_map(subset_map);
			return subset_map;
		}


		static void build_EMT_subset_map(powerset_btree<set_N_value<T>* >& subset_map){
			std::unordered_map<size_t, std::vector<set_N_value<set_N_value<T>* >* > > card_map = subset_map.elements_by_set_cardinality(true);
			std::vector<size_t> ordered_cardinalities = powerset_btree<set_N_value<T>* >::get_sorted_cardinalities(card_map, *subset_map.get_FOD());
			std::reverse(ordered_cardinalities.begin(), ordered_cardinalities.end());
			// powerset_btree of nodes without right child in this powerset.
			// Each node in terminal_connection_tree is associated with the address of the corresponding node in this powerset.
			powerset_btree<set_N_value<set_N_value<T>* >* > terminal_connection_tree(subset_map.get_FOD(), subset_map.get_block_size());

			size_t c = 0;
			if (ordered_cardinalities[0] == subset_map.get_FOD_size()){
				c = 1;
			}

			for (; c < ordered_cardinalities.size(); ++c){
				const std::vector<set_N_value<set_N_value<T>* >* >& elements = card_map[ordered_cardinalities[c]];

				for (size_t i = 0; i < elements.size(); ++i){

					if(!elements[i]->is_null){
						const std::vector<set_N_value<set_N_value<set_N_value<T>* >* >* >& supersets = terminal_connection_tree.supersets_of(elements[i]->set);
						node<set_N_value<T>* >* node_i = (efficient_DST::node<set_N_value<T>* >*) elements[i];
						for (size_t s = 0; s < supersets.size(); ++s){
							node<set_N_value<T>* >* superset = (efficient_DST::node<set_N_value<T>* >*) supersets[s]->value;

							if(superset->is_null || node_i->depth < superset->depth){
								superset->left = node_i;
							}
						}
						for (size_t s = 0; s < supersets.size(); ++s){
							node<set_N_value<T>* >* superset = (efficient_DST::node<set_N_value<T>* >*) supersets[s]->value;
							if(superset->is_null || node_i->depth < superset->depth){
								terminal_connection_tree.nullify(supersets[s]);
							}
						}
					}
				}
				for (size_t i = 0; i < elements.size(); ++i){
					node<set_N_value<T>* >* node = (efficient_DST::node<set_N_value<T>* >*) elements[i];
					if (!node->left){
						terminal_connection_tree.insert(elements[i]->set, elements[i]);
					}
				}
			}
		}


		static void execute_EMT_with_lattice_upper_closure(
				powerset_btree<T>& lattice_support,
				const transform_type_t& transform_type,
				std::function<T(const T&, const T&)> range_binary_operator,
				const std::vector<boost::dynamic_bitset<> >& iota_sequence
		) {
			clock_t t = clock();
			const std::vector<set_N_value<T>* >& lattice_support_elements = lattice_support.elements();

			size_t iota_index;
			for (size_t o = 0; o < iota_sequence.size(); ++o){
				if (transform_type == transform_type_t::zeta)
					iota_index = o;
				else
					iota_index = iota_sequence.size()-1 - o;

				for (size_t i = 0; i < lattice_support_elements.size(); ++i){
					const boost::dynamic_bitset<>& set_B = FOD::set_union(lattice_support_elements[i]->set, iota_sequence[iota_index]);

					if (set_B != lattice_support_elements[i]->set){
						set_N_value<T>* B = lattice_support[set_B];

						if (B){
							B->value = range_binary_operator(B->value, lattice_support_elements[i]->value);
						}
					}
				}
			}
			t = clock() - t;
			std::cout << "time spent computing = " << ((float)t)/CLOCKS_PER_SEC << " sec" << std::endl;
		}


		static void execute_EMT_with_lattice_lower_closure(
				powerset_btree<T>& lattice_support,
				const transform_type_t& transform_type,
				std::function<T(const T&, const T&)> range_binary_operator,
				const std::vector<boost::dynamic_bitset<> >& iota_sequence
		) {
			clock_t t = clock();
			const std::vector<set_N_value<T>* >& lattice_support_elements = lattice_support.elements();

			size_t iota_index;
			for (size_t o = 0; o < iota_sequence.size(); ++o){
				if (transform_type == transform_type_t::zeta)
					iota_index = o;
				else
					iota_index = iota_sequence.size()-1 - o;

				for (size_t i = 0; i < lattice_support_elements.size(); ++i){
					const boost::dynamic_bitset<>& set_B = FOD::set_intersection(lattice_support_elements[i]->set, iota_sequence[iota_index]);

					if (set_B != lattice_support_elements[i]->set){
						set_N_value<T>* B = lattice_support[set_B];

						if (B){
							B->value = range_binary_operator(B->value, lattice_support_elements[i]->value);
						}
					}
				}
			}
			t = clock() - t;
			std::cout << "time spent computing = " << ((float)t)/CLOCKS_PER_SEC << " sec" << std::endl;
		}


		static void build_and_execute_EMT_with_lattice(
				const powerset_btree<T>& support,
				powerset_btree<T>& cropped_lattice_support,
				powerset_btree<set_N_value<T>* >& cropped_lattice_support_dual,
				const transform_type_t& transform_type,
				const order_relation_t& order_relation,
				const operation_t& transform_operation,
				std::vector<boost::dynamic_bitset<> >& iota_sequence
				//std::vector<size_t>& iota_fiber_sequence
		) {
			clock_t t = clock();
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
						cropped_lattice_support,
						iota_sequence,
						neutral_value
				);
				t = clock() - t;
				std::cout << "time spent building = " << ((float)t)/CLOCKS_PER_SEC << " sec" << std::endl;
				t = clock();
				execute_EMT_with_lattice_upper_closure(
						cropped_lattice_support,
						transform_type,
						range_binary_operator,
						iota_sequence
				);
				t = clock() - t;
				std::cout << "time spent calling the computation = " << ((float)t)/CLOCKS_PER_SEC << " sec" << std::endl;
			}else{
				powerset_btree<T>::flip_powerset_from_to(cropped_lattice_support, cropped_lattice_support_dual);
				build_lattice_support_lower_closure(
						cropped_lattice_support,
						cropped_lattice_support_dual,
						iota_sequence,
						neutral_value
				);
				t = clock() - t;
				std::cout << "time spent building = " << ((float)t)/CLOCKS_PER_SEC << " sec" << std::endl;
				t = clock();
				execute_EMT_with_lattice_lower_closure(
						cropped_lattice_support,
						transform_type,
						range_binary_operator,
						iota_sequence
				);
				t = clock() - t;
				std::cout << "time spent calling the computation = " << ((float)t)/CLOCKS_PER_SEC << " sec" << std::endl;
			}
		}


		static inline T addition(const T& a, const T& b){
			return a+b;
		}

		static inline T subtraction(const T& a, const T& b){
			return a-b;
		}

		static inline T multiplication(const T& a, const T& b){
			return a*b;
		}

		static inline T division(const T& a, const T& b){
			return a/b;
		}


		static inline set_N_value<T>* insert_focal_point(
				powerset_btree<T>& focal_points_tree,
				powerset_btree<set_N_value<T>* >& focal_points_dual_tree,
				const boost::dynamic_bitset<>& focal_point,
				const T& neutral_value
		) {
			set_N_value<T>* inserted_focal_point = focal_points_tree.insert(focal_point, neutral_value);
			focal_points_dual_tree.insert(~focal_point, inserted_focal_point);
			DEBUG({
				std::clog << "\nNEW INSERTED FOCAL POINT :\n";
				std::clog << inserted_focal_point->set << std::endl;
			});
			return inserted_focal_point;
		}

		static inline set_N_value<T>* insert_dual_focal_point(
				powerset_btree<T>& focal_points_tree,
				powerset_btree<set_N_value<T>* >& focal_points_dual_tree,
				const boost::dynamic_bitset<>& dual_focal_point,
				const T& neutral_value
		) {
			set_N_value<T>* inserted_focal_point = focal_points_tree.insert(~dual_focal_point, neutral_value);
			focal_points_dual_tree.insert(dual_focal_point, inserted_focal_point);
			DEBUG({
				std::clog << "\nNEW INSERTED FOCAL POINT :\n";
				std::clog << inserted_focal_point->set << std::endl;
			});
			return inserted_focal_point;
		}


		static bool linear_analysis_of_support(
				powerset_btree<T>& focal_points_tree,
				powerset_btree<set_N_value<T>* >& focal_points_dual_tree,
				std::unordered_map<size_t, std::vector<set_N_value<T>* > >& support_card_map,
				std::vector<set_N_value<T>* >& elements_generated_from_support,
				const order_relation_t& order_relation,
				const T& neutral_value
		){

			if(order_relation == order_relation_t::subset){
				std::vector<set_N_value<T>* > support_except_emptyset;
				support_except_emptyset.reserve(focal_points_tree.size()-1);

				for(auto kv : support_card_map) {
					if(kv.first > 0){
						for (size_t i = 0; i < kv.second.size(); ++i) {
							support_except_emptyset.emplace_back(kv.second[i]);
						}
					}
				}

				boost::dynamic_bitset<> fod(focal_points_tree.get_FOD_size());
				fod.set();

				boost::dynamic_bitset<> neg_U = support_except_emptyset[0]->set;
				DEBUG(std::clog << "\nneg_U = "<< support_except_emptyset[0]->set;);

				for (size_t i = 1; i < support_except_emptyset.size(); ++i) {
					const boost::dynamic_bitset<>& A = support_except_emptyset[i]->set;

					DEBUG(std::clog << "\nneg_I = union(U, "<< A << ")\n";);
					const boost::dynamic_bitset<>& neg_I = FOD::set_union(neg_U, A);

					const size_t& I_card = fod.size() - neg_I.count();
					DEBUG(std::clog << "|I| = " << I_card << std::endl;);

					if(I_card > 1){
						DEBUG(std::clog << "=> Linear analysis aborted.\n";);
						return false;
					}else if(I_card == 1){
						// add it to this->structure if it wasn't already there
						if(!focal_points_tree[neg_I]){
							set_N_value<T>* inserted_focal_point = insert_focal_point(focal_points_tree, focal_points_dual_tree, neg_I, neutral_value);
							elements_generated_from_support.push_back(inserted_focal_point);
						}
					}
					neg_U = FOD::set_intersection(neg_U, A);
				}
				// add also fod to this->structure if it wasn't already there
				if(!focal_points_tree[fod]){
					set_N_value<T>* inserted_focal_point = insert_focal_point(focal_points_tree, focal_points_dual_tree, fod, neutral_value);
					elements_generated_from_support.push_back(inserted_focal_point);
				}
			}else{
				std::vector<set_N_value<T>* > support_except_FOD;
				support_except_FOD.reserve(focal_points_tree.size()-1);

				for(auto kv : support_card_map) {
					if(kv.first < focal_points_tree.get_FOD_size()){
						for (size_t i = 0; i < kv.second.size(); ++i) {
							support_except_FOD.emplace_back(kv.second[i]);
						}
					}
				}

				const boost::dynamic_bitset<> emptyset(focal_points_tree.get_FOD_size());

				boost::dynamic_bitset<> U = support_except_FOD[0]->set;
				DEBUG(std::clog << "\nU = "<< support_except_FOD[0]->set;);

				for (size_t i = 1; i < support_except_FOD.size(); ++i) {
					const boost::dynamic_bitset<>& A = support_except_FOD[i]->set;

					DEBUG(std::clog << "\nI = intersection(U, "<< A << ")\n";);
					const boost::dynamic_bitset<>& I = FOD::set_intersection(U, A);

					const size_t& I_card = I.count();
					DEBUG(std::clog << "|I| = " << I_card << std::endl;);

					if(I_card > 1){
						DEBUG(std::clog << "=> Linear analysis aborted.\n";);
						return false;
					}else if(I_card == 1){
						// add it to this->structure if it wasn't already there
						if(!focal_points_tree[I]){
							set_N_value<T>* inserted_focal_point = insert_focal_point(focal_points_tree, focal_points_dual_tree, I, neutral_value);
							elements_generated_from_support.push_back(inserted_focal_point);
						}
					}
					U = FOD::set_union(U, A);
				}
				// add also emptyset to this->structure if it wasn't already there
				if(!focal_points_tree[emptyset]){
					set_N_value<T>* inserted_focal_point = insert_focal_point(focal_points_tree, focal_points_dual_tree, emptyset, neutral_value);
					elements_generated_from_support.push_back(inserted_focal_point);
				}
			}
			return true;
		}


		static bool consonance_check(
				const std::vector<set_N_value<T>* >& elements_generated_from_support,
				std::unordered_map<size_t, std::vector<set_N_value<T>* > >& support_card_map,
				const std::vector<size_t>& support_ordered_cardinalities
		){
			// if at least one focal point has already been discovered
			// OR if there are at least two elements with the same 0-indexed cardinality (i.e. the least one) in structure,
			// then structure cannot be consonant
			if(elements_generated_from_support.size() > 0 || support_card_map[support_ordered_cardinalities[0]].size() > 1)
				return false;

			for (size_t i = 1; i < support_ordered_cardinalities.size(); ++i) {
				// if there are at least two elements with the same i-indexed cardinality in structure
				// OR if the set of cardinality index i-1 is not a subset of the one of cardinality index i,
				// then structure cannot be consonant
				if(support_card_map[support_ordered_cardinalities[i]].size() > 1 || !FOD::is_or_is_subset_of(
						support_card_map[support_ordered_cardinalities[i-1]][0]->set,
						support_card_map[support_ordered_cardinalities[i]][0]->set
					  )
				){
					return false;
				}
			}
			return true;
		}


		static void compute_iota_elements(
				const version_t& version,
				const powerset_btree<T>& support,
				std::vector<boost::dynamic_bitset<> >& iota_sequence
				//std::vector<size_t>& iota_fiber_sequence
		) {
			//powerset_btree<size_t> iota_tree(support.get_FOD(), support.get_block_size());
			powerset_btree<bool> iota_tree(support.get_FOD(), support.get_block_size());

			if(version == version_t::regular){
				for (size_t i = 0; i < support.get_FOD_size(); ++i) {
					boost::dynamic_bitset<> singleton(support.get_FOD_size());
					singleton.set(i);
					const std::vector<set_N_value<T>* >& support_supersets = support.supersets_of(singleton);

					if (support_supersets.size() > 0) {
						boost::dynamic_bitset<> iota_element((const boost::dynamic_bitset<>&) support_supersets[0]->set);

						for (size_t ii = 1; ii < support_supersets.size(); ++ii) {
							iota_element = FOD::set_intersection(iota_element, support_supersets[ii]->set);
							if (iota_element == singleton) {
								break;
							}
						}
						if (!iota_tree[iota_element]){
							set_N_value<bool>* inserted_iota_element = iota_tree.insert(iota_element, true);//i);
							DEBUG({
								std::clog << "\nNEW INSERTED IOTA ELEMENT :\n";
								std::clog << inserted_iota_element->set << std::endl;
							});
						}
					}
				}
			}else{
				for (size_t i = 0; i < support.get_FOD_size(); ++i) {
					boost::dynamic_bitset<> singleton_dual(support.get_FOD_size());
					singleton_dual.set(i);
					singleton_dual.flip();
					const std::vector<set_N_value<T>* >& support_subsets = support.subsets_of(singleton_dual);

					if (support_subsets.size() > 0) {
						boost::dynamic_bitset<> iota_element_dual((const boost::dynamic_bitset<>&) support_subsets[0]->set);

						for (size_t ii = 1; ii < support_subsets.size(); ++ii) {
							iota_element_dual = FOD::set_union(iota_element_dual, support_subsets[ii]->set);
							if (iota_element_dual == singleton_dual) {
								break;
							}
						}
						if (!iota_tree[iota_element_dual]){
							set_N_value<bool>* inserted_iota_element_dual = iota_tree.insert(iota_element_dual, true);//i);
							DEBUG({
								std::clog << "\nNEW INSERTED IOTA ELEMENT DUAL :\n";
								std::clog << inserted_iota_element_dual->set << std::endl;
							});
						}
					}
				}
			}
			std::unordered_map<size_t, std::vector<set_N_value<bool>* > > iota_elements_card_map = iota_tree.elements_by_set_cardinality();
			std::vector<size_t> ordered_cardinalities = powerset_btree<bool>::get_sorted_cardinalities(iota_elements_card_map, *iota_tree.get_FOD());
			iota_sequence.reserve(iota_tree.size());
			//iota_fiber_sequence.reserve(iota_tree.size());

			size_t cc;
			for (size_t c = 0; c < ordered_cardinalities.size(); ++c) {
				if(version == version_t::regular){
					cc = c;
				}else{
					cc = ordered_cardinalities.size()-1 - c;
				}
				const std::vector<set_N_value<bool>* >& iota_elements_nodes = iota_elements_card_map[ordered_cardinalities[cc]];
				for (size_t i = 0; i < iota_elements_nodes.size(); ++i) {
					DEBUG(std::clog << iota_elements_nodes[i]->set << std::endl;);
					iota_sequence.emplace_back(iota_elements_nodes[i]->set);
					//iota_fiber_sequence.emplace_back(iota_elements_nodes[i]->value);
				}
			}
		}


		/*
		 * I use here a notion that I called "focal point" which is the intersection (resp. union)
		 * of a set of focal sets for the superset order relation (resp. subset order relation).
		 * The image of the focal points defines both the zeta and Möbius transforms entirely.
		 * A focal set is also a focal point.
		 *
		 * Compute all focal points.
		 * If F is the number of focal sets and dot_F the number of focal points, then the upper bound complexity is O(F.dot_F),
		 * where dot_F is in [F, 2^N], and N is the FOD size.
		 * The worst case is obtained if all sets of cardinality N-1 (resp. all singletons) are focal sets
		 * for the superset order relation (resp. subset order relation).
		 */
		static void build_semilattice_support(
				powerset_btree<T>& focal_points_tree,
				powerset_btree<set_N_value<T>* >& focal_points_dual_tree,
				std::unordered_map<size_t, std::vector<set_N_value<T>* > >& support_card_map,
				const std::vector<size_t>& support_ordered_cardinalities,
				std::vector<set_N_value<T>* >& elements_generated_from_support,
				const order_relation_t& order_relation,
				const T& neutral_value
		) {

			std::function<boost::dynamic_bitset<>(const boost::dynamic_bitset<>&, const boost::dynamic_bitset<>&)> domain_binary_operator;
			std::function<std::vector<boost::dynamic_bitset<> >(const powerset_btree<T>& powerset, const boost::dynamic_bitset<>&, size_t)> tree_operator;
			std::function<std::vector<boost::dynamic_bitset<> >(const powerset_btree<set_N_value<T>* >& powerset, const boost::dynamic_bitset<>&, size_t)> tree_operator_dual;

			if(order_relation == order_relation_t::subset){
				domain_binary_operator = FOD::set_union;
				tree_operator = powerset_btree<T>::unions_with_not_subsets_of_smaller_than;
				tree_operator_dual = powerset_btree<set_N_value<T>* >::intersections_with_not_subsets_of_smaller_than;
			}else{
				domain_binary_operator = FOD::set_intersection;
				tree_operator = powerset_btree<T>::intersections_with_not_subsets_of_smaller_than;
				tree_operator_dual = powerset_btree<set_N_value<T>* >::unions_with_not_subsets_of_smaller_than;
			}

			// avoid neutral/absorbing sets for our binary_operator
			size_t c_init = 0, c_end = support_ordered_cardinalities.size();
			if (support_ordered_cardinalities[0] == 0)
				c_init = 1;
			if (support_ordered_cardinalities.back() == focal_points_tree.get_FOD_size())
				--c_end;

			// for each cardinality of focal set, from the smallest to the biggest, except emptyset and FOD
			for (size_t c = c_init; c < c_end; ++c) {
				for (size_t i = 0; i < support_card_map[support_ordered_cardinalities[c]].size(); ++i) {
					const boost::dynamic_bitset<>& setA = support_card_map[support_ordered_cardinalities[c]][i]->set;

					// search for sets of same cardinality
					for (size_t ii = 0; ii < i; ++ii) {
						const boost::dynamic_bitset<>& setB = support_card_map[support_ordered_cardinalities[c]][ii]->set;

						// compute their focal point
						const boost::dynamic_bitset<>& focal_point = domain_binary_operator(setA, setB);

						// add it to this->structure if it wasn't already there
						if(!focal_points_tree[focal_point]){
							set_N_value<T>* inserted_focal_point = insert_focal_point(focal_points_tree, focal_points_dual_tree, focal_point, neutral_value);
							elements_generated_from_support.push_back(inserted_focal_point);
						}
					}
					const std::vector<boost::dynamic_bitset<> >& operations_with_not_subsets_of_smaller_than = tree_operator(focal_points_tree, setA, support_ordered_cardinalities[c]-1);

					// search for sets of lower cardinality that are not subsets of the current set
					for (size_t ii = 0; ii < operations_with_not_subsets_of_smaller_than.size(); ++ii) {
						const boost::dynamic_bitset<>& focal_point = operations_with_not_subsets_of_smaller_than[ii];

						// add it to this->structure if it wasn't already there
						if(!focal_points_tree[focal_point]){
							set_N_value<T>* inserted_focal_point = insert_focal_point(focal_points_tree, focal_points_dual_tree, focal_point, neutral_value);
							elements_generated_from_support.push_back(inserted_focal_point);
						}
					}
				}
			}
			// for each pure focal point (focal point that is not a focal set)
			size_t i = 0;
			while (true) {
				if (i >= elements_generated_from_support.size())
					break;
				const boost::dynamic_bitset<>& setA = elements_generated_from_support[i]->set;
				const std::vector<boost::dynamic_bitset<> >& operations_with_not_subsets_of_smaller_than = tree_operator(focal_points_tree, setA, elements_generated_from_support[i]->cardinality-1);

				// search for sets of lower cardinality that are not subsets of the current set
				for (size_t ii = 0; ii < operations_with_not_subsets_of_smaller_than.size(); ++ii) {
					const boost::dynamic_bitset<>& focal_point = operations_with_not_subsets_of_smaller_than[ii];

					// add it to this->structure if it wasn't already there
					if(!focal_points_tree[focal_point]){
						set_N_value<T>* inserted_focal_point = insert_focal_point(focal_points_tree, focal_points_dual_tree, focal_point, neutral_value);
						elements_generated_from_support.push_back(inserted_focal_point);
					}
				}
				const std::vector<boost::dynamic_bitset<> >& dual_operations_with_not_subsets_of_smaller_than = tree_operator_dual(
						focal_points_dual_tree, ~setA,
						focal_points_tree.get_FOD_size() - elements_generated_from_support[i]->cardinality-1);

				// search for sets of higher cardinality that are not supersets of the current set
				for (size_t ii = 0; ii < dual_operations_with_not_subsets_of_smaller_than.size(); ++ii) {
					const boost::dynamic_bitset<>& dual_focal_point = dual_operations_with_not_subsets_of_smaller_than[ii];

					// add it to this->structure if it wasn't already there
					if(!focal_points_dual_tree[dual_focal_point]){
						set_N_value<T>* inserted_focal_point = insert_dual_focal_point(focal_points_tree, focal_points_dual_tree, dual_focal_point, neutral_value);
						elements_generated_from_support.push_back(inserted_focal_point);
					}
				}
				++i;
			}
			DEBUG({
				std::clog << "\nFocal points: \n";
				print<T>(std::clog, focal_points_tree);
			});
		}


		/*
		 * In this function, this->iota_sequence is supposed to contain all regular iota elements (i.e. join-irreducible of this lattice).
		 */
		static void build_lattice_support_upper_closure(
				powerset_btree<T>& focal_points_tree,
				std::vector<boost::dynamic_bitset<> >& iota_sequence,
				const T& neutral_value
		){
			//std::vector<boost::dynamic_bitset<> > iota_sequence_dual;
			//std::vector<size_t> iota_fiber_sequence_dual;
			compute_iota_elements(
					version_t::regular,
					focal_points_tree,
					iota_sequence
					//iota_fiber_sequence_dual
			);

			for (size_t o = 0; o < iota_sequence.size(); ++o) {
				const std::vector<set_N_value<T>* >& lattice_elements = focal_points_tree.elements();
				for (size_t i = 0; i < lattice_elements.size(); ++i) {
					boost::dynamic_bitset<> new_set = FOD::set_union(iota_sequence[o], lattice_elements[i]->set);
					if(!focal_points_tree[new_set]){
						set_N_value<T>* inserted_focal_point = focal_points_tree.insert(new_set, neutral_value);
						DEBUG({
							std::clog << "\nNEW INSERTED FOCAL POINT :\n";
							std::clog << inserted_focal_point->set << std::endl;
						});
					}
				}
			}
			DEBUG({
				std::clog << "\nCropped lattice support: \n";
				print<T>(std::clog, focal_points_tree);
			});
		}


		static void build_lattice_support_upper_closure(
				powerset_btree<T>& focal_points_tree,
				std::vector<boost::dynamic_bitset<> >& iota_sequence,
				const std::vector<T>& vec
		){
			//std::vector<boost::dynamic_bitset<> > iota_sequence_dual;
			//std::vector<size_t> iota_fiber_sequence_dual;
			compute_iota_elements(
					version_t::regular,
					focal_points_tree,
					iota_sequence
					//iota_fiber_sequence_dual
			);

			for (size_t o = 0; o < iota_sequence.size(); ++o) {
				const std::vector<set_N_value<T>* >& lattice_elements = focal_points_tree.elements();
				for (size_t i = 0; i < lattice_elements.size(); ++i) {
					boost::dynamic_bitset<> new_set = FOD::set_union(iota_sequence[o], lattice_elements[i]->set);
					if(!focal_points_tree[new_set]){
						set_N_value<T>* inserted_focal_point = focal_points_tree.insert(new_set, vec[new_set.to_ulong()]);
						DEBUG({
							std::clog << "\nNEW INSERTED FOCAL POINT :\n";
							std::clog << inserted_focal_point->set << std::endl;
						});
					}
				}
			}
			DEBUG({
				std::clog << "\nCropped lattice support: \n";
				print<T>(std::clog, focal_points_tree);
			});
		}


		/*
		 * In this function, this->iota_sequence is supposed to contain all regular iota elements (i.e. join-irreducible of this lattice).
		 */
		static void build_lattice_support_lower_closure(
				powerset_btree<T>& focal_points_tree,
				powerset_btree<set_N_value<T>* >& focal_points_dual_tree,
				std::vector<boost::dynamic_bitset<> >& iota_sequence,
				const T& neutral_value
		){
			//std::vector<boost::dynamic_bitset<> > iota_sequence_dual;
			//std::vector<size_t> iota_fiber_sequence_dual;
			compute_iota_elements(
					version_t::dual,
					focal_points_tree,
					iota_sequence
					//iota_fiber_sequence_dual
			);
			clock_t t, count, counte;
			counte = 0;
			count = 0;
			clock_t tl = clock();
			for (size_t o = 0; o < iota_sequence.size(); ++o) {
				t = clock();
				const std::vector<set_N_value<T>* >& lattice_elements = focal_points_tree.elements();
				t = clock() - t;
				counte += t;
				for (size_t i = 0; i < lattice_elements.size(); ++i) {
					boost::dynamic_bitset<> new_set = FOD::set_intersection(iota_sequence[o], lattice_elements[i]->set);
					t = clock();
					if(!focal_points_tree[new_set]){
						insert_focal_point(focal_points_tree, focal_points_dual_tree, new_set, neutral_value);
					}
					t = clock() - t;
					count += t;
				}
			}
			tl = clock() - tl;
			std::cout << "time spent retrieving lattice = " << ((float)counte)/CLOCKS_PER_SEC << " sec" << std::endl;
			std::cout << "time spent building lattice = " << ((float)tl)/CLOCKS_PER_SEC << " sec" << std::endl;
			std::cout << "time spent inserting = " << ((float)count)/CLOCKS_PER_SEC << " sec" << std::endl;
			/*
			 * 			for (size_t o = 0; o < iota_sequence.size(); ++o) {
				const std::vector<set_N_value<set_N_value<T>* >* >& lattice_elements = focal_points_dual_tree.elements();
				for (size_t i = 0; i < lattice_elements.size(); ++i) {
					boost::dynamic_bitset<> new_set = FOD::set_union(iota_sequence[iota_sequence.size()-1-o], lattice_elements[i]->set);
					if(!focal_points_dual_tree[new_set]){
						insert_dual_focal_point(focal_points_tree, focal_points_dual_tree, new_set, neutral_value);
					}
				}
			}
			 */
			DEBUG({
				std::clog << "\nCropped lattice support: \n";
				print<T>(std::clog, focal_points_tree);
			});
		}


		static void build_lattice_support_lower_closure(
				powerset_btree<T>& focal_points_tree,
				powerset_btree<set_N_value<T>* >& focal_points_dual_tree,
				std::vector<boost::dynamic_bitset<> >& iota_sequence,
				const std::vector<T>& vec
		){
			//std::vector<boost::dynamic_bitset<> > iota_sequence_dual;
			//std::vector<size_t> iota_fiber_sequence_dual;
			compute_iota_elements(
					version_t::dual,
					focal_points_tree,
					iota_sequence
					//iota_fiber_sequence_dual
			);
			clock_t t, count, counte;
			counte = 0;
			count = 0;
			clock_t tl = clock();
			for (size_t o = 0; o < iota_sequence.size(); ++o) {
				t = clock();
				const std::vector<set_N_value<T>* >& lattice_elements = focal_points_tree.elements();
				t = clock() - t;
				counte += t;
				for (size_t i = 0; i < lattice_elements.size(); ++i) {
					boost::dynamic_bitset<> new_set = FOD::set_intersection(iota_sequence[o], lattice_elements[i]->set);
					t = clock();
					if(!focal_points_tree[new_set]){
						insert_focal_point(focal_points_tree, focal_points_dual_tree, new_set, vec[new_set.to_ulong()]);
					}
					t = clock() - t;
					count += t;
				}
			}
			tl = clock() - tl;
			std::cout << "time spent retrieving lattice = " << ((float)counte)/CLOCKS_PER_SEC << " sec" << std::endl;
			std::cout << "time spent building lattice = " << ((float)tl)/CLOCKS_PER_SEC << " sec" << std::endl;
			std::cout << "time spent inserting = " << ((float)count)/CLOCKS_PER_SEC << " sec" << std::endl;
			/*
			 * 			for (size_t o = 0; o < iota_sequence.size(); ++o) {
				const std::vector<set_N_value<set_N_value<T>* >* >& lattice_elements = focal_points_dual_tree.elements();
				for (size_t i = 0; i < lattice_elements.size(); ++i) {
					boost::dynamic_bitset<> new_set = FOD::set_union(iota_sequence[iota_sequence.size()-1-o], lattice_elements[i]->set);
					if(!focal_points_dual_tree[new_set]){
						insert_dual_focal_point(focal_points_tree, focal_points_dual_tree, new_set, neutral_value);
					}
				}
			}
			 */
			DEBUG({
				std::clog << "\nCropped lattice support: \n";
				print<T>(std::clog, focal_points_tree);
			});
		}
	};
}		// namespace efficient_DST

#endif // EFFICIENT_DST_COMPUTATION_SCHEME_HPP
