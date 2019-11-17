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
	enum class version_t: bool { regular, dual };
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
				//std::vector<size_t>& iota_fiber_sequence,
				scheme_type_t& scheme_type
		){
			T neutral_value;
			if (transform_operation == operation_t::addition){
				neutral_value = 0;
			} else{
				neutral_value = 1;
			}

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
				compute_iota_elements(version_t::dual, support, iota_sequence);//, iota_fiber_sequence);
				//for (size_t o = 0; o < iota_sequence.size(); ++o){
				//	iota_sequence[o].flip();
				//}
			}else{
				compute_iota_elements(version_t::regular, support, iota_sequence);//, iota_fiber_sequence);
			}
		}


		static void extract_semilattice_support(
			const std::vector<T>& powerset_values,
			powerset_btree<T, N>& focal_points_tree,
			const order_relation_t& order_relation
		) {
			size_t powerset_size = pow(2, N);
			if(powerset_values.size() != powerset_size){
				std::cerr << "\nThe given vector does not feature the same size as the powerset of the given FOD.\n";
				return;
			}
			std::unordered_map<T, std::vector<std::bitset<N> > > focal_points_map;
			focal_points_map.reserve(powerset_size);

			for (size_t n = 0; n < powerset_size; ++n){
				std::bitset<N> set(n);

				if (focal_points_map.find(powerset_values[n]) == focal_points_map.end()){
					focal_points_map.emplace(powerset_values[n], (std::vector<std::bitset<N> >) {set});
					focal_points_map[powerset_values[n]].reserve(N);
				}else{
					bool new_point = true;
					std::vector<std::bitset<N> >& focal_points_with_same_image = focal_points_map[powerset_values[n]];
					for (size_t i = 0; i < focal_points_with_same_image.size(); ++i){
						if (order_relation == order_relation_t::subset){
							if (FOD<N>::is_subset_of(set, focal_points_with_same_image[i])){
								focal_points_with_same_image[i] = set;
								new_point = false;
							}else if (FOD<N>::is_subset_of(focal_points_with_same_image[i], set)){
								new_point = false;
								break;
							}
						}else{
							if (FOD<N>::is_superset_of(set, focal_points_with_same_image[i])){
								focal_points_with_same_image[i] = set;
								new_point = false;
							}else if (FOD<N>::is_superset_of(focal_points_with_same_image[i], set)){
								new_point = false;
								break;
							}
						}
					}
					if (new_point){
						focal_points_with_same_image.emplace_back(set);
					}
				}
			}
			size_t to_reserve = 0;
			for (auto kv : focal_points_map){
				to_reserve += kv.second.size();
			}
			//focal_points_register.reserve(to_reserve);
			for (auto kv : focal_points_map){
				for (size_t i = 0; i < kv.second.size(); ++i){
					focal_points_tree.insert(kv.second[i], kv.first);
				}
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


		static std::vector<T> FMT(
				const std::vector<T>& powerset_values,
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
				powerset_btree<T, N>& focal_points_tree,
				const transform_type_t& transform_type,
				const order_relation_t& order_relation,
				std::function<T(const T&, const T&)> range_binary_operator
		) {
			if (transform_type == transform_type_t::zeta){
				powerset_btree<T, N> focal_points_N_initial_values(focal_points_tree);
				T val;
				const std::vector<set_N_value<T, N>* >& focal_points = focal_points_N_initial_values.elements();
				std::vector<set_N_value<T, N>* > elements;

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
				std::unordered_map<size_t, std::vector<set_N_value<T, N>* > > focal_points_card_map = focal_points_tree.elements_by_set_cardinality();
				std::vector<size_t> focal_points_ordered_cardinalities;
				std::vector<set_N_value<T, N>* > elements;

				if (order_relation == order_relation_t::subset) {
					focal_points_tree.get_FOD()->sort_cardinalities(focal_points_ordered_cardinalities, focal_points_card_map, order_t::ascending);
				}else{
					focal_points_tree.get_FOD()->sort_cardinalities(focal_points_ordered_cardinalities, focal_points_card_map, order_t::descending);
				}

				for (size_t c = 0; c < focal_points_ordered_cardinalities.size(); ++c) {
					const std::vector<set_N_value<T, N>* >& focal_points = focal_points_card_map[focal_points_ordered_cardinalities[c]];
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
				powerset_btree<T, N>& support,
				const transform_type_t& transform_type,
				const order_relation_t& order_relation,
				std::function<T(const T&, const T&)> range_binary_operator,
				const T& neutral_value
		) {
			std::unordered_map<size_t, std::vector<set_N_value<T, N>* > > support_card_map = support.elements_by_set_cardinality();
			std::vector<size_t> support_ordered_cardinalities;
			T value, preceding_value = neutral_value;

			if (order_relation == order_relation_t::subset) {
				support.get_FOD()->sort_cardinalities(support_ordered_cardinalities, support_card_map, order_t::ascending);
			}else{
				support.get_FOD()->sort_cardinalities(support_ordered_cardinalities, support_card_map, order_t::descending);
			}
			for (size_t i = 0; i < support_ordered_cardinalities.size(); ++i) {
				set_N_value<T, N>* set_value = support_card_map[support_ordered_cardinalities[i]][0];
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
			const bool& is_almost_bayesian = linear_analysis_of_support(
					support,
					focal_points_tree,
					order_relation,
					neutral_value
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


		static void execute_EMT_with_lower_semilattice(
				powerset_btree<T, N>& focal_points_tree,
				//std::unordered_map<size_t, T>& focal_points_register,
				const transform_type_t& transform_type,
				std::function<T(const T&, const T&)> range_binary_operator,
				const std::vector<std::bitset<N> >& iota_sequence
		) {
			clock_t t;

			if (iota_sequence.size() == 0)
				return;

			std::vector<std::bitset<N> > sync_sequence;
			sync_sequence.reserve(iota_sequence.size());

			sync_sequence.emplace_back(iota_sequence[0]);
			for (size_t i = 1; i < iota_sequence.size(); ++i){
				sync_sequence.emplace_back(sync_sequence[i-1] | iota_sequence[i]);
			}
			std::unordered_map<std::bitset<N>, set_N_value<T, N>* > proxy_map;
			build_proxy_superset_map(proxy_map, focal_points_tree, iota_sequence);

			//std::bitset<N> emptyset(focal_points_tree.get_FOD_size());
			//emptyset.set(1);

			t = clock();
			size_t iota_index;
			const std::vector<set_N_value<T, N>* >& focal_points = focal_points_tree.elements();
			for (size_t i = 0; i < iota_sequence.size(); ++i){
				if (transform_type == transform_type_t::zeta)
					iota_index = i;
				else
					iota_index = iota_sequence.size()-1 - i;

				for (size_t e = 0; e < focal_points.size(); ++e){
					const std::bitset<N>& set_B = focal_points[e]->set | iota_sequence[iota_index];

					if (focal_points[e]->set != set_B){
						const set_N_value<T, N>* X = proxy_map[set_B];

						/*if(focal_points[i]->set == iota_element){
							std::clog << focal_points[i]->set << " + " << iota_sequence[iota_index] << " = " << set_B << std::endl;
							if(X)
								std::clog << " => " << X->value->value << "\t <- " << focal_points_tree.get_FOD()->to_string(X->set) << std::endl;
						}*/

						if (X && FOD<N>::is_subset_of(X->set, focal_points[e]->set | sync_sequence[iota_index])){
							T& val = focal_points[e]->value;
							val = range_binary_operator(val, X->value);
						}
					}
				}
			}
			t = clock() - t;
			std::cout << (((float) t)/CLOCKS_PER_SEC) << std::endl;
		}


		static void execute_EMT_with_upper_semilattice(
				powerset_btree<T, N>& focal_points_tree,
				//std::unordered_map<size_t, T>& focal_points_register,
				const transform_type_t& transform_type,
				std::function<T(const T&, const T&)> range_binary_operator,
				const std::vector<std::bitset<N> >& iota_sequence
		) {
			if (iota_sequence.size() == 0)
				return;

			std::vector<std::bitset<N> > sync_sequence;
			sync_sequence.reserve(iota_sequence.size());

			sync_sequence.emplace_back(iota_sequence[0]);
			for (size_t i = 1; i < iota_sequence.size(); ++i){
				sync_sequence.emplace_back(sync_sequence[i-1] & iota_sequence[i]);
			}

			std::unordered_map<std::bitset<N>, set_N_value<T, N>* > proxy_map;
			build_proxy_subset_map(proxy_map, focal_points_tree, iota_sequence);

			size_t iota_index;
			const std::vector<set_N_value<T, N>* >& focal_points = focal_points_tree.elements();
			DEBUG(std::clog << "STARTING COMPUTATION SEQUENCE:\n";);
			for (size_t i = 0; i < iota_sequence.size(); ++i){
				if (transform_type == transform_type_t::zeta)
					iota_index = i;
				else
					iota_index = iota_sequence.size()-1 - i;

				DEBUG(std::clog << "iota element " << iota_sequence[iota_index] << ":\n";);
				DEBUG(std::clog << "sync element " << sync_sequence[iota_index] << ":\n";);

				for (size_t e = 0; e < focal_points.size(); ++e){
					const std::bitset<N>& set_B = focal_points[e]->set & iota_sequence[iota_index];

					if (focal_points[e]->set != set_B){
						const set_N_value<T, N>* X = proxy_map[set_B];

						DEBUG(std::clog << "\tConsidering " << focal_points[e]->set << " :\n";);
						DEBUG(
							if(X)
								std::clog << "\t\tproxy = " << X->set << "\n";
							else
								std::clog << "\t\tno proxy\n";
						);

						if (X && FOD<N>::is_superset_of(X->set, focal_points[e]->set & sync_sequence[iota_index])){
							DEBUG(std::clog << "\t\t" << focal_points[e]->set << " <- operation("
									<< focal_points[e]->set << " = " << focal_points[e]->value <<  ", " << X->set << " = " << X->value << ")\n";);
							T& val = focal_points[e]->value;
							val = range_binary_operator(val, X->value);
						}
					}
				}
			}
			DEBUG(std::clog << "DONE\n";);
		}
	/*
		static void execute_EMT_with_upper_semilattice(
				//powerset_btree<T, N>& focal_points_tree,
				std::unordered_map<size_t, T>& focal_points_register,
				const size_t& fod_size,
				const transform_type_t& transform_type,
				std::function<T(const T&, const T&)> range_binary_operator,
				const std::vector<std::bitset<N> >& iota_sequence
		) {
			if (iota_sequence.size() == 0)
				return;

			std::vector<std::bitset<N> > sync_sequence;
			sync_sequence.reserve(iota_sequence.size());

			sync_sequence.emplace_back(iota_sequence[0]);
			for (size_t o = 1; o < iota_sequence.size(); ++o){
				sync_sequence.emplace_back(FOD::set_intersection(sync_sequence[o-1], iota_sequence[o]));
			}

			std::vector<set_N_value<T, N>* > focal_points = focal_points_tree.elements();
			std::unordered_map<size_t, std::pair<const size_t, T>* > proxy_map;
			const powerset_btree<set_N_value<T, N>* >& focal_points_subset_map = EMT_subset_map(focal_points_tree, iota_sequence);
			std::unordered_map<size_t, std::pair<const size_t, T>* >& proxy_map,
			const FOD& fod,
			const std::unordered_map<size_t, T>& focal_points_register,
			const std::vector<std::bitset<N> >& iota_sequence,
			const order_relation_t& order_relation

			std::bitset<N> fod(fod_size);
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
					const std::bitset<N>& set_B = FOD::set_intersection(focal_points[i]->set, iota_sequence[iota_index]);

					if (focal_points[i]->set != set_B){
						const set_N_value<set_N_value<T, N>* >* X = focal_points_subset_map.biggest_subset_of(set_B);

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


		static void execute_EMT_with_semilattice(
				powerset_btree<T, N>& focal_points_tree,
				powerset_btree<set_N_value<T, N>* >& focal_points_dual_tree,
				const transform_type_t& transform_type,
				const order_relation_t& order_relation,
				std::function<T(const T&, const T&)> range_binary_operator,
				const std::vector<std::bitset<N> >& iota_sequence
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
				compute_iota_elements(version_t::dual, support, iota_sequence);
			}else{
				compute_iota_elements(version_t::regular, support, iota_sequence); //iota_fiber_sequence);
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


		/*static powerset_btree<set_N_value<T, N>* > EMT_superset_map(
				const powerset_btree<set_N_value<T, N>* >& focal_points_dual_tree,
				const std::vector<std::bitset<N> >& iota_sequence_dual
		) {
			powerset_btree<set_N_value<T, N>* > superset_map(focal_points_dual_tree.get_FOD(), focal_points_dual_tree.get_block_size());
			const std::vector<set_N_value<set_N_value<T, N>* >* >& elements = focal_points_dual_tree.elements();
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
				superset_map.insert_null_node(std::bitset<N>(focal_points_dual_tree.get_FOD_size(), new_set));
			}

			build_EMT_superset_map(superset_map);
			return superset_map;
		}


		static powerset_btree<set_N_value<T, N>* > EMT_superset_map(
				const powerset_btree<T, N>& focal_points_tree,
				const std::vector<std::bitset<N> >& iota_sequence
		) {
			powerset_btree<set_N_value<T, N>* > superset_map(focal_points_tree.get_FOD(), focal_points_tree.get_block_size());
			const std::vector<set_N_value<T, N>* >& elements = focal_points_tree.elements();

			for (size_t e = 0; e < elements.size(); ++e) {
				superset_map.insert(elements[e]->set, elements[e]);
				for (size_t i = 0; i < iota_sequence.size(); ++i) {
					superset_map.insert_null_node(FOD::set_union(elements[e]->set, iota_sequence[i]));
				}
			}
			build_EMT_superset_map(superset_map);
			return superset_map;
		}


		static void build_EMT_superset_map(powerset_btree<set_N_value<T, N>* >& superset_map){
			std::unordered_map<size_t, std::vector<set_N_value<set_N_value<T, N>* >* > > card_map = superset_map.elements_by_set_cardinality(true);
			std::vector<size_t> ordered_cardinalities = powerset_btree<set_N_value<T, N>* >::get_sorted_cardinalities(card_map, *superset_map.get_FOD());
			// powerset_btree of nodes without right child in this powerset.
			// Each node in terminal_connection_tree is associated with the address of the corresponding node in this powerset.
			powerset_btree<set_N_value<set_N_value<T, N>* >* > terminal_connection_tree(superset_map.get_FOD(), superset_map.get_block_size());

			size_t c = 0;
			if (ordered_cardinalities[0] == 0){
				c = 1;
			}

			for (; c < ordered_cardinalities.size(); ++c){
				const std::vector<set_N_value<set_N_value<T, N>* >* >& elements = card_map[ordered_cardinalities[c]];

				for (size_t i = 0; i < elements.size(); ++i){

					if(!elements[i]->is_null){
						const std::vector<set_N_value<set_N_value<set_N_value<T, N>* >* >* >& subsets = terminal_connection_tree.subsets_of(elements[i]->set);
						node<set_N_value<T, N>* >* node_i = (efficient_DST::node<set_N_value<T, N>* >*) elements[i];
						for (size_t s = 0; s < subsets.size(); ++s){
							node<set_N_value<T, N>* >* subset = (efficient_DST::node<set_N_value<T, N>* >*) subsets[s]->value;

							if(subset->is_null || node_i->depth > subset->depth){
								subset->right = node_i;
							}
						}
						for (size_t s = 0; s < subsets.size(); ++s){
							node<set_N_value<T, N>* >* subset = (efficient_DST::node<set_N_value<T, N>* >*) subsets[s]->value;
							if(subset->is_null || node_i->depth > subset->depth){
								terminal_connection_tree.nullify(subsets[s]);
							}
						}
					}
				}
				for (size_t i = 0; i < elements.size(); ++i){
					node<set_N_value<T, N>* >* node = (efficient_DST::node<set_N_value<T, N>* >*) elements[i];
					if (node->is_null && !node->right){
						terminal_connection_tree.insert(elements[i]->set, elements[i]);
					}
				}
			}
		}*/


		/*static void build_proxy_map(
				std::unordered_map<size_t, set_N_value<T, N>* >& proxy_map,
				FOD* fod,
				const powerset_btree<T, N>& focal_points_tree,
				//std::unordered_map<size_t, T>& focal_points_register,
				const std::vector<std::bitset<N> >& iota_sequence,
				const order_relation_t& order_relation
		) {
			proxy_map.reserve(iota_sequence.size() * focal_points_tree.size());

			//std::unordered_map<size_t, std::vector<size_t> > focal_points_register_card_map;
			//focal_points_register_card_map.reserve(fod->size());

			const std::vector<set_N_value<T, N>* >& focal_points = focal_points_tree.elements();
			for (size_t e = 0; e < focal_points.size(); ++e){
				proxy_map.emplace(focal_points[e]->set.to_ulong(), focal_points[e]);

				//const size_t& card = std::bitset<N>(fod->size(), kv.first).count();
				//if (focal_points_register_card_map.find(card) == focal_points_register_card_map.end()){
				//	focal_points_register_card_map[card] = std::vector<size_t>();
				//	if (card == 1)
				//		focal_points_register_card_map[card].reserve(1);
				//	else
				//		focal_points_register_card_map[card].reserve(focal_points_register.size());
				//}
				//focal_points_register_card_map[card].emplace_back(kv.first);
			}

			if (order_relation == order_relation_t::subset)
				build_proxy_subset_map(proxy_map, fod, focal_points_tree, iota_sequence);
			else
				build_proxy_superset_map(proxy_map, fod, focal_points_tree, iota_sequence);
		}*/


		static void build_proxy_subset_map(
				std::unordered_map<std::bitset<N>, set_N_value<T, N>* >& proxy_map,
				const powerset_btree<T, N>& focal_points_tree,
				const std::vector<std::bitset<N> >& iota_sequence
		) {
			proxy_map.reserve(iota_sequence.size() * focal_points_tree.size());
			const std::vector<set_N_value<T, N>* >& focal_points = focal_points_tree.elements();

			for (size_t e = 0; e < focal_points.size(); ++e){
				proxy_map.emplace(focal_points[e]->set, focal_points[e]);
			}

			powerset_btree<bool, N> sets_missing_proxies(focal_points_tree.get_FOD(), iota_sequence.size() * focal_points_tree.size());

			for (size_t i = 0; i < iota_sequence.size(); ++i) {
				for (size_t e = 0; e < focal_points.size(); ++e){
					const std::bitset<N>& set = focal_points[e]->set & iota_sequence[i];
					bool insertion = proxy_map.emplace(set, nullptr).second;
					//set_sequence.emplace_back(set.to_ulong());
					if (insertion)
						sets_missing_proxies.insert(set, true);
				}
			}

			std::unordered_map<size_t, std::vector<set_N_value<T, N>* > > focal_points_card_map = focal_points_tree.elements_by_set_cardinality();
			std::vector<size_t> ordered_cardinalities;
			sets_missing_proxies.get_FOD()->sort_cardinalities(ordered_cardinalities, focal_points_card_map, order_t::descending);

			size_t c = 0;
			if (ordered_cardinalities[0] == N){
				c = 1;
			}
			for (; c < ordered_cardinalities.size(); ++c){
				const std::vector<set_N_value<T, N>* >& focal_points_with_same_size = focal_points_card_map[ordered_cardinalities[c]];

				for (size_t i = 0; i < focal_points_with_same_size.size(); ++i){
					const std::vector<set_N_value<bool, N>* >& supersets = sets_missing_proxies.supersets_of(focal_points_with_same_size[i]->set);

					for (size_t s = 0; s < supersets.size(); ++s){
						proxy_map[supersets[s]->set] = focal_points_with_same_size[i];
						sets_missing_proxies.nullify(supersets[s]);
					}
				}
			}
		}


		static void build_proxy_superset_map(
				std::unordered_map<std::bitset<N>, set_N_value<T, N>* >& proxy_map,
				const powerset_btree<T, N>& focal_points_tree,
				const std::vector<std::bitset<N> >& iota_sequence
		) {
			proxy_map.reserve(iota_sequence.size() * focal_points_tree.size());
			const std::vector<set_N_value<T, N>* >& focal_points = focal_points_tree.elements();

			for (size_t e = 0; e < focal_points.size(); ++e){
				proxy_map.emplace(focal_points[e]->set, focal_points[e]);
			}

			clock_t t,temp_t;
			t = clock();
			size_t count = 0;
			powerset_btree<bool, N> sets_missing_proxies(focal_points_tree.get_FOD(), iota_sequence.size() * focal_points_tree.size());

			for (size_t i = 0; i < iota_sequence.size(); ++i) {
				for (size_t e = 0; e < focal_points.size(); ++e){
					const std::bitset<N>& set = focal_points[e]->set | iota_sequence[i];
					bool insertion = proxy_map.emplace(set, nullptr).second;
					//set_sequence.emplace_back(set.to_ulong());
					if (insertion){
						temp_t = clock();
						sets_missing_proxies.insert(set, true);
						temp_t = clock() - temp_t;
						count += temp_t;
					}
				}
			}
			//t = clock() - t;
			//std::cout << (((float) t)/CLOCKS_PER_SEC) << " ";

			std::unordered_map<size_t, std::vector<set_N_value<T, N>* > > focal_points_card_map = focal_points_tree.elements_by_set_cardinality();
			std::vector<size_t> ordered_cardinalities;
			sets_missing_proxies.get_FOD()->sort_cardinalities(ordered_cardinalities, focal_points_card_map, order_t::ascending);

			size_t c = 0;
			if (ordered_cardinalities[0] == 0){
				c = 1;
			}

			for (; c < ordered_cardinalities.size(); ++c){
				const std::vector<set_N_value<T, N>* >& focal_points_with_same_size = focal_points_card_map[ordered_cardinalities[c]];

				for (size_t i = 0; i < focal_points_with_same_size.size(); ++i){
					const std::vector<set_N_value<bool, N>* >& subsets = sets_missing_proxies.subsets_of(focal_points_with_same_size[i]->set);

					for (size_t s = 0; s < subsets.size(); ++s){
						proxy_map[subsets[s]->set] = focal_points_with_same_size[i];
						sets_missing_proxies.nullify(subsets[s]);
					}
				}
			}
		}


		/*static powerset_btree<set_N_value<T, N>* > EMT_subset_map(
				const powerset_btree<T, N>& focal_points_tree,
				const std::vector<std::bitset<N> >& iota_sequence
		) {
			powerset_btree<set_N_value<T, N>* > subset_map(focal_points_tree.get_FOD(), focal_points_tree.get_block_size());
			const std::vector<set_N_value<T, N>* >& elements = focal_points_tree.elements();

			for (size_t e = 0; e < elements.size(); ++e) {
				subset_map.insert(elements[e]->set, elements[e]);
				for (size_t i = 0; i < iota_sequence.size(); ++i) {
					subset_map.insert_null_node(FOD::set_intersection(elements[e]->set, iota_sequence[i]));
				}
			}
			build_EMT_subset_map(subset_map);
			return subset_map;
		}


		static void build_EMT_subset_map(powerset_btree<set_N_value<T, N>* >& subset_map){
			std::unordered_map<size_t, std::vector<set_N_value<set_N_value<T, N>* >* > > card_map = subset_map.elements_by_set_cardinality(true);
			std::vector<size_t> ordered_cardinalities = powerset_btree<set_N_value<T, N>* >::get_sorted_cardinalities(card_map, *subset_map.get_FOD());
			std::reverse(ordered_cardinalities.begin(), ordered_cardinalities.end());
			// powerset_btree of nodes without right child in this powerset.
			// Each node in terminal_connection_tree is associated with the address of the corresponding node in this powerset.
			powerset_btree<set_N_value<set_N_value<T, N>* >* > terminal_connection_tree(subset_map.get_FOD(), subset_map.get_block_size());

			size_t c = 0;
			if (ordered_cardinalities[0] == subset_map.get_FOD_size()){
				c = 1;
			}

			for (; c < ordered_cardinalities.size(); ++c){
				const std::vector<set_N_value<set_N_value<T, N>* >* >& elements = card_map[ordered_cardinalities[c]];

				for (size_t i = 0; i < elements.size(); ++i){

					if(!elements[i]->is_null){
						const std::vector<set_N_value<set_N_value<set_N_value<T, N>* >* >* >& supersets = terminal_connection_tree.supersets_of(elements[i]->set);
						node<set_N_value<T, N>* >* node_i = (efficient_DST::node<set_N_value<T, N>* >*) elements[i];
						for (size_t s = 0; s < supersets.size(); ++s){
							node<set_N_value<T, N>* >* superset = (efficient_DST::node<set_N_value<T, N>* >*) supersets[s]->value;

							if(superset->is_null || node_i->depth < superset->depth){
								superset->left = node_i;
							}
						}
						for (size_t s = 0; s < supersets.size(); ++s){
							node<set_N_value<T, N>* >* superset = (efficient_DST::node<set_N_value<T, N>* >*) supersets[s]->value;
							if(superset->is_null || node_i->depth < superset->depth){
								terminal_connection_tree.nullify(supersets[s]);
							}
						}
					}
				}
				for (size_t i = 0; i < elements.size(); ++i){
					node<set_N_value<T, N>* >* node = (efficient_DST::node<set_N_value<T, N>* >*) elements[i];
					if (!node->left){
						terminal_connection_tree.insert(elements[i]->set, elements[i]);
					}
				}
			}
		}*/


		static void execute_EMT_with_lattice_upper_closure(
				powerset_btree<T, N>& lattice_support,
				//std::unordered_map<size_t, T>& focal_points_register,
				//const size_t& fod_size,
				const transform_type_t& transform_type,
				std::function<T(const T&, const T&)> range_binary_operator,
				const std::vector<std::bitset<N> >& iota_sequence
		) {
			std::unordered_map<std::bitset<N>, set_N_value<T, N>* > focal_points_register;
			focal_points_register.reserve(lattice_support.size());

			const std::vector<set_N_value<T, N>* >& lattice_support_elements = lattice_support.elements();
			for (size_t e = 0; e < lattice_support_elements.size(); ++e){
				focal_points_register[lattice_support_elements[e]->set] = lattice_support_elements[e];
			}

			size_t iota_index;
			for (size_t i = 0; i < iota_sequence.size(); ++i){
				if (transform_type == transform_type_t::Mobius)
					iota_index = i;
				else
					iota_index = iota_sequence.size()-1 - i;

				/*for (auto kv : focal_points_register){
					const size_t& set_B = (std::bitset<N>(fod_size, kv.first) | iota_sequence[iota_index]).to_ulong();

					if (set_B != kv.first){
						auto B = focal_points_register.find(set_B);

						if (B != focal_points_register.end()){
							B->second = range_binary_operator(B->second, kv.second);
						}
					}
				}
				*/
				for (size_t e = 0; e < lattice_support_elements.size(); ++e){
					const std::bitset<N>& set_B = lattice_support_elements[e]->set | iota_sequence[iota_index];

					if (set_B != lattice_support_elements[e]->set){
						auto B = focal_points_register.find(set_B);

						if (B != focal_points_register.end()){
							T& val = B->second->value;
							val = range_binary_operator(val, lattice_support_elements[e]->value);
						}
					}
				}
			}
		}


		static void execute_EMT_with_lattice_lower_closure(
				powerset_btree<T, N>& lattice_support,
				//std::unordered_map<size_t, T>& focal_points_register,
				//const size_t& fod_size,
				const transform_type_t& transform_type,
				std::function<T(const T&, const T&)> range_binary_operator,
				const std::vector<std::bitset<N> >& iota_sequence
		) {
			std::unordered_map<std::bitset<N>, set_N_value<T, N>* > focal_points_register;
			focal_points_register.reserve(lattice_support.size());

			const std::vector<set_N_value<T, N>* >& lattice_support_elements = lattice_support.elements();
			for (size_t e = 0; e < lattice_support_elements.size(); ++e){
				focal_points_register[lattice_support_elements[e]->set] = lattice_support_elements[e];
			}

			size_t iota_index;
			for (size_t i = 0; i < iota_sequence.size(); ++i){
				if (transform_type == transform_type_t::Mobius)
					iota_index = i;
				else
					iota_index = iota_sequence.size()-1 - i;

				for (size_t e = 0; e < lattice_support_elements.size(); ++e){
					const std::bitset<N>& set_B = lattice_support_elements[e]->set & iota_sequence[iota_index];

					if (set_B != lattice_support_elements[e]->set){
						auto B = focal_points_register.find(set_B);

						if (B != focal_points_register.end()){
							T& val = B->second->value;
							val = range_binary_operator(val, lattice_support_elements[e]->value);
						}
					}
				}

				/*for (auto kv : focal_points_register){
					const size_t& set_B = (std::bitset<N>(fod_size, kv.first) & iota_sequence[iota_index]).to_ulong();

					if (set_B != kv.first){
						auto B = focal_points_register.find(set_B);

						if (B != focal_points_register.end()){
							B->second = range_binary_operator(B->second, kv.second);
						}
					}
				}*/
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


		static inline set_N_value<T, N>* insert_focal_point(
				powerset_btree<T, N>& focal_points_tree,
				powerset_btree<set_N_value<T, N>*, N >& focal_points_dual_tree,
				const std::bitset<N>& focal_point,
				const T& neutral_value
		) {
			set_N_value<T, N>* inserted_focal_point = focal_points_tree.insert(focal_point, neutral_value);
			focal_points_dual_tree.insert(~focal_point, inserted_focal_point);
			DEBUG({
				std::clog << "\nNEW INSERTED FOCAL POINT :\n";
				std::clog << inserted_focal_point->set << std::endl;
			});
			return inserted_focal_point;
		}

		static inline set_N_value<T, N>* insert_dual_focal_point(
				powerset_btree<T, N>& focal_points_tree,
				powerset_btree<set_N_value<T, N>*, N >& focal_points_dual_tree,
				const std::bitset<N>& dual_focal_point,
				const T& neutral_value
		) {
			set_N_value<T, N>* inserted_focal_point = focal_points_tree.insert(~dual_focal_point, neutral_value);
			focal_points_dual_tree.insert(dual_focal_point, inserted_focal_point);
			DEBUG({
				std::clog << "\nNEW INSERTED FOCAL POINT :\n";
				std::clog << inserted_focal_point->set << std::endl;
			});
			return inserted_focal_point;
		}


		static bool linear_analysis_of_support(
				const powerset_btree<T, N>& support,
				powerset_btree<T, N>& focal_points_tree,
				const order_relation_t& order_relation,
				T& neutral_value
		){
			std::unordered_set<std::bitset<N>> focal_points_register;
			focal_points_register.reserve(support.size() + N + 1);
			const std::vector<set_N_value<T, N>* >& elements = support.elements();
			for (size_t i = 0; i < elements.size(); ++i){
				focal_points_register.emplace(elements[i]->set);
				focal_points_tree.insert(elements[i]->set, elements[i]->value);
			}

			if(order_relation == order_relation_t::subset){
				const std::vector<set_N_value<T, N>* >& support_except_emptyset = support.strict_supersets_of(std::bitset<N>(0));

				std::bitset<N> fod = 0;
				fod.set();

				std::bitset<N> neg_U = support_except_emptyset[0]->set;
				DEBUG(std::clog << "\nneg_U = "<< support_except_emptyset[0]->set;);

				for (size_t i = 1; i < support_except_emptyset.size(); ++i) {
					const std::bitset<N>& A = support_except_emptyset[i]->set;

					DEBUG(std::clog << "\nneg_I = union(U, "<< A << ")\n";);
					const std::bitset<N>& neg_I = neg_U | A;

					const size_t& I_card = fod.size() - neg_I.count();
					DEBUG(std::clog << "|I| = " << I_card << std::endl;);

					if(I_card > 1){
						DEBUG(std::clog << "=> Linear analysis aborted.\n";);
						return false;
					}else if(I_card == 1){
						bool insertion = focal_points_register.emplace(neg_I).second;
						if (insertion)
							focal_points_tree.insert(neg_I, neutral_value);
					}
					neg_U &= A;
				}
				bool insertion = focal_points_register.emplace(fod).second;
				if (insertion)
					focal_points_tree.insert(fod, neutral_value);
			}else{
				std::bitset<N> fod = 0;
				fod.set();
				const std::vector<set_N_value<T, N>* >& support_except_FOD = support.strict_subsets_of(fod);
				const std::bitset<N> emptyset = 0;

				std::bitset<N> U = support_except_FOD[0]->set;
				DEBUG(std::clog << "\nU = "<< support_except_FOD[0]->set;);

				for (size_t i = 1; i < support_except_FOD.size(); ++i) {
					const std::bitset<N>& A = support_except_FOD[i]->set;

					DEBUG(std::clog << "\nI = intersection(U, "<< A << ")\n";);
					const std::bitset<N>& I = U & A;

					const size_t& I_card = I.count();
					DEBUG(std::clog << "|I| = " << I_card << std::endl;);

					if(I_card > 1){
						DEBUG(std::clog << "=> Linear analysis aborted.\n";);
						return false;
					}else if(I_card == 1){
						bool insertion = focal_points_register.emplace(I).second;
						if (insertion)
							focal_points_tree.insert(I, neutral_value);
					}
					U |= A;
				}
				bool insertion = focal_points_register.emplace(emptyset).second;
				if (insertion)
					focal_points_tree.insert(emptyset, neutral_value);
			}
			return true;
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


		static void compute_iota_elements(
				const version_t& version,
				const powerset_btree<T, N>& support,
				std::vector<std::bitset<N> >& iota_sequence
		) {
			std::unordered_set<std::bitset<N>> iota_elements;
			iota_elements.reserve(N);
			std::unordered_map<size_t, std::vector<std::bitset<N> > > iota_elements_card_map;
			iota_elements_card_map.reserve(N);

			if(version == version_t::regular){
				for (size_t i = 0; i < N; ++i) {
					std::bitset<N> singleton = 0;
					singleton.set(i);
					const std::vector<set_N_value<T, N>* >& support_supersets = support.supersets_of(singleton);

					if (support_supersets.size() > 0) {
						std::bitset<N> iota_element((const std::bitset<N>&) support_supersets[0]->set);

						for (size_t ii = 1; ii < support_supersets.size(); ++ii) {
							iota_element &= support_supersets[ii]->set;
							if (iota_element == singleton) {
								break;
							}
						}
						bool insertion = iota_elements.emplace(iota_element).second;
						if (insertion){
							size_t cardinality = iota_element.count();
							if (iota_elements_card_map.find(cardinality) == iota_elements_card_map.end()){
								iota_elements_card_map.emplace(cardinality, std::vector<std::bitset<N> >());
								iota_elements_card_map[cardinality].reserve(N);
							}
							iota_elements_card_map[cardinality].emplace_back(iota_element);
							DEBUG({
								std::clog << "\nNEW INSERTED IOTA ELEMENT :\n";
								std::clog << iota_element << std::endl;
							});
						}
					}
				}
			}else{
				support.elements();
				for (size_t i = 0; i < N; ++i) {
					std::bitset<N> singleton_dual = 0;
					singleton_dual.set(i);
					singleton_dual.flip();
					const std::vector<set_N_value<T, N>* >& support_subsets = support.subsets_of(singleton_dual);

					if (support_subsets.size() > 0) {
						std::bitset<N> iota_element_dual((const std::bitset<N>&) support_subsets[0]->set);
						//std::clog << "Computing iota element dual associated to :" << singleton_dual << "\n";
						//std::clog << iota_element_dual << std::endl;
						for (size_t ii = 1; ii < support_subsets.size(); ++ii) {
							//std::clog << support_subsets[ii]->set << std::endl;
							iota_element_dual |= support_subsets[ii]->set;
							if (iota_element_dual == singleton_dual) {
								break;
							}
						}
						bool insertion = iota_elements.emplace(iota_element_dual).second;
						if (insertion){
							size_t cardinality = iota_element_dual.count();
							if (iota_elements_card_map.find(cardinality) == iota_elements_card_map.end()){
								iota_elements_card_map.emplace(cardinality, std::vector<std::bitset<N> >());
								iota_elements_card_map[cardinality].reserve(N);
							}
							iota_elements_card_map[cardinality].emplace_back(iota_element_dual);
							DEBUG({
								std::clog << "\nNEW INSERTED IOTA ELEMENT DUAL :\n";
								std::clog << iota_element_dual << std::endl;
							});
						}
					}
				}
			}
			std::vector<size_t> ordered_cardinalities;
			order_t order = order_t::ascending;
			if(version == version_t::dual){
				order = order_t::descending;
			}
			support.get_FOD()->sort_cardinalities(ordered_cardinalities, iota_elements_card_map, order);
			iota_sequence.reserve(iota_elements.size());

			for (size_t c = 0; c < ordered_cardinalities.size(); ++c) {
				const std::vector<std::bitset<N> >& iota_elements_c = iota_elements_card_map[ordered_cardinalities[c]];
				for (size_t i = 0; i < iota_elements_c.size(); ++i) {
					DEBUG(std::clog << iota_elements_c[i] << std::endl;);
					iota_sequence.emplace_back(iota_elements_c[i]);
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
		 *
		static void build_semilattice_support(
				powerset_btree<T, N>& focal_points_tree,
				powerset_btree<set_N_value<T, N>* >& focal_points_dual_tree,
				std::unordered_map<size_t, std::vector<set_N_value<T, N>* > >& support_card_map,
				const std::vector<size_t>& support_ordered_cardinalities,
				std::vector<set_N_value<T, N>* >& elements_generated_from_support,
				const order_relation_t& order_relation,
				const T& neutral_value
		) {

			std::function<std::bitset<N>(const std::bitset<N>&, const std::bitset<N>&)> domain_binary_operator;
			std::function<std::vector<std::bitset<N> >(const powerset_btree<T, N>& powerset, const std::bitset<N>&, size_t)> tree_operator;
			std::function<std::vector<std::bitset<N> >(const powerset_btree<set_N_value<T, N>* >& powerset, const std::bitset<N>&, size_t)> tree_operator_dual;

			if(order_relation == order_relation_t::subset){
				domain_binary_operator = FOD::set_union;
				tree_operator = powerset_btree<T, N>::unions_with_not_subsets_of_smaller_than;
				tree_operator_dual = powerset_btree<set_N_value<T, N>* >::intersections_with_not_subsets_of_smaller_than;
			}else{
				domain_binary_operator = FOD::set_intersection;
				tree_operator = powerset_btree<T, N>::intersections_with_not_subsets_of_smaller_than;
				tree_operator_dual = powerset_btree<set_N_value<T, N>* >::unions_with_not_subsets_of_smaller_than;
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
					const std::bitset<N>& setA = support_card_map[support_ordered_cardinalities[c]][i]->set;

					// search for sets of same cardinality
					for (size_t ii = 0; ii < i; ++ii) {
						const std::bitset<N>& setB = support_card_map[support_ordered_cardinalities[c]][ii]->set;

						// compute their focal point
						const std::bitset<N>& focal_point = domain_binary_operator(setA, setB);

						// add it to this->structure if it wasn't already there
						if(!focal_points_tree[focal_point]){
							set_N_value<T, N>* inserted_focal_point = insert_focal_point(focal_points_tree, focal_points_dual_tree, focal_point, neutral_value);
							elements_generated_from_support.push_back(inserted_focal_point);
						}
					}
					const std::vector<std::bitset<N> >& operations_with_not_subsets_of_smaller_than = tree_operator(focal_points_tree, setA, support_ordered_cardinalities[c]-1);

					// search for sets of lower cardinality that are not subsets of the current set
					for (size_t ii = 0; ii < operations_with_not_subsets_of_smaller_than.size(); ++ii) {
						const std::bitset<N>& focal_point = operations_with_not_subsets_of_smaller_than[ii];

						// add it to this->structure if it wasn't already there
						if(!focal_points_tree[focal_point]){
							set_N_value<T, N>* inserted_focal_point = insert_focal_point(focal_points_tree, focal_points_dual_tree, focal_point, neutral_value);
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
				const std::bitset<N>& setA = elements_generated_from_support[i]->set;
				const std::vector<std::bitset<N> >& operations_with_not_subsets_of_smaller_than = tree_operator(focal_points_tree, setA, elements_generated_from_support[i]->cardinality-1);

				// search for sets of lower cardinality that are not subsets of the current set
				for (size_t ii = 0; ii < operations_with_not_subsets_of_smaller_than.size(); ++ii) {
					const std::bitset<N>& focal_point = operations_with_not_subsets_of_smaller_than[ii];

					// add it to this->structure if it wasn't already there
					if(!focal_points_tree[focal_point]){
						set_N_value<T, N>* inserted_focal_point = insert_focal_point(focal_points_tree, focal_points_dual_tree, focal_point, neutral_value);
						elements_generated_from_support.push_back(inserted_focal_point);
					}
				}
				const std::vector<std::bitset<N> >& dual_operations_with_not_subsets_of_smaller_than = tree_operator_dual(
						focal_points_dual_tree, ~setA,
						focal_points_tree.get_FOD_size() - elements_generated_from_support[i]->cardinality-1);

				// search for sets of higher cardinality that are not supersets of the current set
				for (size_t ii = 0; ii < dual_operations_with_not_subsets_of_smaller_than.size(); ++ii) {
					const std::bitset<N>& dual_focal_point = dual_operations_with_not_subsets_of_smaller_than[ii];

					// add it to this->structure if it wasn't already there
					if(!focal_points_dual_tree[dual_focal_point]){
						set_N_value<T, N>* inserted_focal_point = insert_dual_focal_point(focal_points_tree, focal_points_dual_tree, dual_focal_point, neutral_value);
						elements_generated_from_support.push_back(inserted_focal_point);
					}
				}
				++i;
			}
			DEBUG({
				std::clog << "\nFocal points: \n";
				print<T>(std::clog, focal_points_tree);
			});
		}*/


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
				const powerset_btree<T, N>& support,
				powerset_btree<T, N>& focal_points_tree,
				const order_relation_t& order_relation,
				const T& neutral_value
		) {
			std::unordered_set<std::bitset<N>> focal_points_register;
			std::vector<std::bitset<N> > focal_points;
			const std::vector<set_N_value<T, N>* >& support_elements = support.elements();
			focal_points_register.reserve(2 * N * support.size());
			focal_points.reserve(2 * N * support.size());

			for (size_t i = 0; i < support_elements.size(); ++i){
				focal_points_register.emplace(support_elements[i]->set);
				focal_points.emplace_back(support_elements[i]->set);
				focal_points_tree.insert(support_elements[i]->set, support_elements[i]->value);
			}

			if(order_relation == order_relation_t::subset){
				for (size_t i = 0; i < support_elements.size(); ++i){
					for (size_t ii = i+1; ii < focal_points.size(); ++ii){
						const std::bitset<N>& focal_point = support_elements[i]->set | focal_points[ii];
						bool insertion = focal_points_register.emplace(focal_point).second;
						if (insertion){
							focal_points.emplace_back(focal_point);
							focal_points_tree.insert(focal_point, neutral_value);
							DEBUG({
								std::clog << "\nNEW INSERTED FOCAL POINT :\n";
								std::clog << focal_point << std::endl;
							});
						}
					}
				}
			}else{
				for (size_t i = 0; i < support_elements.size(); ++i){
					for (size_t ii = i+1; ii < focal_points.size(); ++ii){
						const std::bitset<N>& focal_point = support_elements[i]->set & focal_points[ii];
						bool insertion = focal_points_register.emplace(focal_point).second;

						if (insertion){
							focal_points.emplace_back(focal_point);
							focal_points_tree.insert(focal_point, neutral_value);
							DEBUG({
								std::clog << "\nNEW INSERTED FOCAL POINT :\n";
								std::clog << focal_point << std::endl;
							});
						}
					}
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
		 *
		static void build_semilattice_support(
				std::vector<std::bitset<N> >& focal_points,
				std::unordered_set<size_t, T>& focal_points_register,
				const size_t& fod_size,
				std::unordered_map<size_t, std::vector<set_N_value<T, N>* > >& support_card_map,
				const order_t& order,
				const std::vector<size_t>& support_ordered_cardinalities,
				const order_relation_t& order_relation,
				const T& neutral_value
		) {
			// avoid neutral/absorbing sets for our binary_operator
			size_t c_init = 0, c_end = support_ordered_cardinalities.size();
			if (order == order_t::ascending){
				if (support_ordered_cardinalities[0] == 0)
					c_init = 1;
				if (support_ordered_cardinalities.back() == fod_size)
					--c_end;
			}else{
				if (support_ordered_cardinalities[0] == fod_size)
					c_init = 1;
				if (support_ordered_cardinalities.back() == 0)
					--c_end;
			}

			for (size_t c = c_init; c < c_end; ++c){
				for (size_t i = 0; i < support_card_map[support_ordered_cardinalities[c]].size(); ++i){
					focal_points_register.emplace(support_card_map[support_ordered_cardinalities[c]][i]->set.to_ulong(), support_card_map[support_ordered_cardinalities[c]][i]->value);
					focal_points.emplace_back(support_card_map[support_ordered_cardinalities[c]][i]->set);
				}
			}
			std::vector<std::bitset<N> > support = focal_points;

			if(order_relation == order_relation_t::subset){
				for (size_t i = 0; i < support.size(); ++i){
					for (size_t ii = i+1; ii < focal_points.size(); ++ii){
						const std::bitset<N>& focal_point = support[i] | focal_points[ii];
						bool insertion = focal_points_register.emplace(focal_point.to_ulong(), neutral_value).second;
						if (insertion)
							focal_points.emplace_back(focal_point);
					}
				}
			}else{
				for (size_t i = 0; i < support.size(); ++i){
					for (size_t ii = i+1; ii < focal_points.size(); ++ii){
						const std::bitset<N>& focal_point = support[i] & focal_points[ii];
						bool insertion = focal_points_register.emplace(focal_point.to_ulong(), neutral_value).second;
						if (insertion)
							focal_points.emplace_back(focal_point);
					}
				}
			}

			DEBUG({
				std::clog << "\nFocal points: \n";
				print<T>(std::clog, focal_points_tree);
			});
		}*/


		/*
		 * In this function, this->iota_sequence is supposed to contain all regular iota elements (i.e. join-irreducible of this lattice).
		 */
		static void build_lattice_support_upper_closure(
				const powerset_btree<T, N>& support,
				powerset_btree<T, N>& cropped_lattice_support,
				std::vector<std::bitset<N> >& iota_sequence,
				const T& neutral_value
		){
			compute_iota_elements(
					version_t::regular,
					support,
					iota_sequence
			);

			std::unordered_set<std::bitset<N>> focal_points_register;
			std::vector<std::bitset<N> > focal_points;
			focal_points.reserve(2 * N * support.size());
			focal_points_register.reserve(2 * N * support.size());

			const std::vector<set_N_value<T, N>* >& elements = support.elements();
			for (size_t i = 0; i < elements.size(); ++i){
				focal_points_register.emplace(elements[i]->set);
				focal_points.emplace_back(elements[i]->set);
				cropped_lattice_support.insert(elements[i]->set, elements[i]->value);
			}

			for (size_t i = 0; i < iota_sequence.size(); ++i) {
				for (size_t e = 0; e < focal_points.size(); ++e) {
					std::bitset<N> new_set = iota_sequence[i] | focal_points[e];
					bool insertion = focal_points_register.emplace(new_set).second;
					if (insertion){
						focal_points.emplace_back(new_set);
						cropped_lattice_support.insert(new_set, neutral_value);
						DEBUG({
							std::clog << "\nNEW INSERTED FOCAL POINT :\n";
							std::clog << new_set << std::endl;
						});
					}
				}
			}
			DEBUG({
				std::clog << "\nCropped lattice support: \n";
				cropped_lattice_support.print(std::clog);
			});
		}


		static void build_lattice_support_upper_closure(
				const std::vector<T>& vec,
				powerset_btree<T, N>& cropped_lattice_support,
				std::vector<std::bitset<N> >& iota_sequence
		){
			extract_semilattice_support(
					vec,
					cropped_lattice_support,
					order_relation_t::subset
			);

			std::vector<std::bitset<N> > focal_points;
			std::unordered_set<std::bitset<N>> focal_points_register;
			focal_points.reserve(2 * N * cropped_lattice_support.size());
			focal_points_register.reserve(2 * N * cropped_lattice_support.size());

			const std::vector<set_N_value<T, N>* >& focal_points_elements = cropped_lattice_support.elements();
			for (size_t e = 0; e < focal_points_elements.size(); ++e){
				focal_points.emplace_back(focal_points_elements[e]->set);
				focal_points_register.emplace(focal_points_elements[e]->set);
			}

			compute_iota_elements(
					version_t::regular,
					cropped_lattice_support,
					iota_sequence
			);

			for (size_t i = 0; i < iota_sequence.size(); ++i) {
				for (size_t e = 0; e < focal_points.size(); ++e) {
					std::bitset<N> new_set = iota_sequence[i] | focal_points[e];
					bool insertion = focal_points_register.emplace(new_set).second;
					if(insertion){
						focal_points.emplace_back(new_set);
						cropped_lattice_support.insert(new_set, vec[new_set.to_ulong()]);
						DEBUG({
							std::clog << "\nNEW INSERTED FOCAL POINT :\n";
							std::clog << new_set << std::endl;
						});
					}
				}
			}
			DEBUG({
				std::clog << "\nCropped lattice support: \n";
				cropped_lattice_support.print(std::clog);
			});
		}


		/*
		 * In this function, this->iota_sequence is supposed to contain all regular iota elements (i.e. join-irreducible of this lattice).
		 */
		static void build_lattice_support_lower_closure(
				const powerset_btree<T, N>& support,
				powerset_btree<T, N>& cropped_lattice_support,
				std::vector<std::bitset<N> >& iota_sequence,
				const T& neutral_value
		){
			compute_iota_elements(
					version_t::dual,
					support,
					iota_sequence
			);

			std::vector<std::bitset<N> > focal_points;
			std::unordered_set<std::bitset<N>> focal_points_register;
			focal_points.reserve(2 * N * support.size());
			focal_points_register.reserve(2 * N * support.size());

			const std::vector<set_N_value<T, N>* >& elements = support.elements();
			for (size_t i = 0; i < elements.size(); ++i){
				focal_points_register.emplace(elements[i]->set);
				focal_points.emplace_back(elements[i]->set);
				cropped_lattice_support.insert(elements[i]->set, elements[i]->value);
			}

			for (size_t i = 0; i < iota_sequence.size(); ++i) {
				for (size_t e = 0; e < focal_points.size(); ++e) {
					std::bitset<N> new_set = iota_sequence[i] & focal_points[e];
					bool insertion = focal_points_register.emplace(new_set).second;
					if(insertion){
						focal_points.emplace_back(new_set);
						cropped_lattice_support.insert(new_set, neutral_value);
						DEBUG({
							std::clog << "\nNEW INSERTED FOCAL POINT :\n";
							std::clog << new_set << std::endl;
						});
					}
				}
			}
			/*
			 * 			for (size_t o = 0; o < iota_sequence.size(); ++o) {
				const std::vector<set_N_value<set_N_value<T, N>* >* >& lattice_elements = focal_points_dual_tree.elements();
				for (size_t i = 0; i < lattice_elements.size(); ++i) {
					std::bitset<N> new_set = FOD::set_union(iota_sequence[iota_sequence.size()-1-o], lattice_elements[i]->set);
					if(!focal_points_dual_tree[new_set]){
						insert_dual_focal_point(focal_points_tree, focal_points_dual_tree, new_set, neutral_value);
					}
				}
			}
			*/
			DEBUG({
				std::clog << "\nCropped lattice support: \n";
				cropped_lattice_support.print(std::clog);
			});
		}


		static void build_lattice_support_lower_closure(
				const std::vector<T>& vec,
				powerset_btree<T, N>& cropped_lattice_support,
				std::vector<std::bitset<N> >& iota_sequence
		){
			extract_semilattice_support(
					vec,
					cropped_lattice_support,
					order_relation_t::superset
			);

			std::vector<std::bitset<N> > focal_points;
			std::unordered_set<std::bitset<N>> focal_points_register;
			focal_points.reserve(2 * N * cropped_lattice_support.size());
			focal_points_register.reserve(2 * N * cropped_lattice_support.size());

			const std::vector<set_N_value<T, N>* >& focal_points_elements = cropped_lattice_support.elements();
			for (size_t e = 0; e < focal_points_elements.size(); ++e){
				focal_points.emplace_back(focal_points_elements[e]->set);
				focal_points_register.emplace(focal_points_elements[e]->set);
			}

			compute_iota_elements(
					version_t::dual,
					cropped_lattice_support,
					iota_sequence
			);
			for (size_t i = 0; i < iota_sequence.size(); ++i) {
				for (size_t e = 0; e < focal_points.size(); ++e) {
					std::bitset<N> new_set = iota_sequence[i] & focal_points[e];
					bool insertion = focal_points_register.emplace(new_set).second;
					if(insertion){
						focal_points.emplace_back(new_set);
						cropped_lattice_support.insert(new_set, vec[new_set.to_ulong()]);
						DEBUG({
							std::clog << "\nNEW INSERTED FOCAL POINT :\n";
							std::clog << new_set << std::endl;
						});
					}
				}
			}
			/*
			 * 			for (size_t o = 0; o < iota_sequence.size(); ++o) {
				const std::vector<set_N_value<set_N_value<T, N>* >* >& lattice_elements = focal_points_dual_tree.elements();
				for (size_t i = 0; i < lattice_elements.size(); ++i) {
					std::bitset<N> new_set = FOD::set_union(iota_sequence[iota_sequence.size()-1-o], lattice_elements[i]->set);
					if(!focal_points_dual_tree[new_set]){
						insert_dual_focal_point(focal_points_tree, focal_points_dual_tree, new_set, neutral_value);
					}
				}
			}*/
			DEBUG({
				std::clog << "\nCropped lattice support: \n";
				cropped_lattice_support.print(std::clog);
			});
		}
	};
}		// namespace efficient_DST

#endif // EFFICIENT_DST_COMPUTATION_SCHEME_HPP
