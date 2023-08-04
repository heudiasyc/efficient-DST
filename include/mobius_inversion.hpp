/*
 * Copyright (C) 2019-2023  Maxime Chaveroche (maxime.chaveroche@gmail.com)
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the CeCILL License, either version 2.1 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * CeCILL License for more details.
 * 
 * You should have received a copy of the CeCILL License
 * along with this program. If not, see <http://www.cecill.info/licences/Licence_CeCILL_V2.1-en.html>.
 */

#ifndef EFFICIENT_DST_MOBIUS_INVERSION_TEMPLATE_HPP
#define EFFICIENT_DST_MOBIUS_INVERSION_TEMPLATE_HPP


#include <unordered_set>

#include "macros.hpp"
#include <powerset_btree.hpp>


namespace efficient_DST{

	enum class scheme_type_t: int8_t { direct, consonant, semilattice, lattice, reduced_FMT };

	template<size_t N, typename T = float>
	struct iota_elements {
		typedef typename sample_space<N>::subset subset;

		static inline void join_irreducible(
			const powerset_btree<N, T>& support,
			std::vector<subset >& sequence
		){
			// join-irreducible elements of the lattice support
			std::unordered_set<subset > manifest;
			manifest.reserve(N);
			std::map<size_t, std::vector<subset >, std::less<size_t> > iota_elements_card_map;
			subset singleton = 1;
			std::vector<set_N_value<N, T> const * > support_supersets;
			support_supersets.reserve(support.size());
			for (size_t i = 0; i < N; ++i) {
				support.supersets_of(singleton, support_supersets);

				if (support_supersets.size() > 0) {
					subset iota_element((const subset&) support_supersets[0]->set);

					for (size_t ii = 1; ii < support_supersets.size(); ++ii) {
						iota_element &= support_supersets[ii]->set;
						if (iota_element == singleton) {
							break;
						}
					}
					bool insertion = manifest.emplace(iota_element).second;
					if (insertion){
						size_t cardinality = iota_element.count();
						if (iota_elements_card_map.find(cardinality) == iota_elements_card_map.end()){
							iota_elements_card_map.emplace(cardinality, std::vector<subset >());
						}
						iota_elements_card_map[cardinality].emplace_back(iota_element);
					}
				}
				singleton <<= 1;
				support_supersets.clear();
			}
			sequence.reserve(manifest.size());

			for (const auto& c_iota_elements : iota_elements_card_map) {
				const std::vector<subset >& elements = c_iota_elements.second;
				for (size_t i = 0; i < elements.size(); ++i) {
					sequence.emplace_back(elements[i]);
				}
			}
//			return sequence;
		}

		static inline void meet_irreducible(
			const powerset_btree<N, T>& support,
			std::vector<subset >& sequence
		){
			// meet-irreducible elements of the lattice support
//			support.print();
			std::unordered_set<subset > manifest;
			manifest.reserve(N);
			std::map<size_t, std::vector<subset >, std::greater<size_t> > iota_elements_card_map;
			subset singleton = 1;
			subset singleton_dual = ~singleton;
			std::vector<set_N_value<N, T> const * > support_subsets;
			support_subsets.reserve(support.size());
			for (size_t i = 0; i < N; ++i) {
				support.subsets_of(singleton_dual, support_subsets);

				if (support_subsets.size() > 0) {
					subset iota_element((const subset&) support_subsets[0]->set);
					for (size_t ii = 1; ii < support_subsets.size(); ++ii) {
						iota_element |= support_subsets[ii]->set;
						if (iota_element == singleton_dual) {
							break;
						}
					}
					bool insertion = manifest.emplace(iota_element).second;
					if (insertion){
						size_t cardinality = iota_element.count();
						if (iota_elements_card_map.find(cardinality) == iota_elements_card_map.end()){
							iota_elements_card_map.emplace(cardinality, std::vector<subset >());
						}
						iota_elements_card_map[cardinality].emplace_back(iota_element);
					}
				}
				singleton <<= 1;
				singleton_dual = ~singleton;
				support_subsets.clear();
			}
			sequence.reserve(manifest.size());

			for (const auto& c_iota_elements : iota_elements_card_map) {
				const std::vector<subset >& elements = c_iota_elements.second;
				for (size_t i = 0; i < elements.size(); ++i) {
					sequence.emplace_back(elements[i]);
				}
			}
//			return sequence;
		}
	};

	template<size_t N, typename T = float>
	struct up_inclusion {
		typedef typename sample_space<N>::subset subset;

		static inline void compute_iota_elements(
			const powerset_btree<N, T>& support,
			std::vector<subset >& sequence
		){
			iota_elements<N, T>::meet_irreducible(support, sequence);
		}

		static inline void compute_iota_elements_dual(
			const powerset_btree<N, T>& support,
			std::vector<subset >& sequence
		){
			iota_elements<N, T>::join_irreducible(support, sequence);
		}

		static inline size_t target_index_offset_FMT(){
			return 1;
		}

		static inline size_t source_index_offset_FMT(){
			return 0;
		}

		static inline subset FMT_target(const subset& set) {
			return ~set;
		}

		static inline subset set_operation(const subset& a, const subset& b){
			return a & b;
		}

		static inline size_t set_operation(const size_t& a, const size_t& b){
			return a & b;
		}

		static inline subset set_dual_operation(const subset& a, const subset& b){
			return a | b;
		}

		static inline size_t set_dual_operation(const size_t& a, const size_t& b){
			return a | b;
		}

		static inline bool naturally_ranked(const size_t a, const size_t& b){
			return a >= b;
		}

		static inline std::map<size_t, std::vector<size_t>, std::greater<size_t> > index_card_mapping(powerset_btree<N, T>& tree){
			return tree.elements_indices_by_descending_cardinality();
		}

		static inline std::map<size_t, std::vector<set_N_value<N, T>* >, std::greater<size_t> > card_mapping(powerset_btree<N, T>& tree){
			return tree._elements_by_descending_cardinality();
		}

		static inline std::map<size_t, std::vector<size_t>, std::less<size_t> > index_card_mapping_dual(powerset_btree<N, T>& tree){
			return tree.elements_indices_by_ascending_cardinality();
		}

		static inline std::map<size_t, std::vector<set_N_value<N, T>* >, std::less<size_t> > card_mapping_dual(powerset_btree<N, T>& tree){
			return tree._elements_by_ascending_cardinality();
		}

		static inline void elements_related_to(const powerset_btree<N, T>& tree, const subset& set, std::vector<set_N_value<N, T> const * >& elements){
			tree.supersets_of(set, elements);
		}

		static inline void addresses_related_to(const powerset_btree<N, set_N_value<N, T>* >& tree, const subset& set, std::vector<set_N_value<N, set_N_value<N, T>* > const * >& addresses){
			tree.supersets_of(set, addresses);
		}

		static inline void elements_dually_related_to(const powerset_btree<N, T>& tree, const subset& set, std::vector<set_N_value<N, T> const * >& elements){
			tree.subsets_of(set, elements);
		}

		static inline void sets_dually_related_to(const powerset_btree<N, bool>& tree, const subset& set, std::vector<size_t>& sets){
			tree.subsets_of(set, sets);
		}

		static inline subset absorbing_set_for_operation(){
			return 0;
		}

		static inline subset absorbing_set_for_dual_operation(){
			return ~absorbing_set_for_operation();
		}
	};

	template<size_t N, typename T = float>
	struct down_inclusion {
		typedef typename sample_space<N>::subset subset;

		static inline void compute_iota_elements(
			const powerset_btree<N, T>& support,
			std::vector<subset >& sequence
		){
			iota_elements<N, T>::join_irreducible(support, sequence);
		}

		static inline void compute_iota_elements_dual(
			const powerset_btree<N, T>& support,
			std::vector<subset >& sequence
		){
			iota_elements<N, T>::meet_irreducible(support, sequence);
		}

		static inline size_t target_index_offset_FMT(){
			return 0;
		}

		static inline size_t source_index_offset_FMT(){
			return 1;
		}

		static inline subset FMT_target(const subset& set) {
			return set;
		}

		static inline subset set_operation(const subset& a, const subset& b){
			return a | b;
		}

		static inline size_t set_operation(const size_t& a, const size_t& b){
			return a | b;
		}

		static inline subset set_dual_operation(const subset& a, const subset& b){
			return a & b;
		}

		static inline size_t set_dual_operation(const size_t& a, const size_t& b){
			return a & b;
		}

		static inline bool naturally_ranked(const size_t a, const size_t& b){
			return a <= b;
		}

		static inline std::map<size_t, std::vector<size_t>, std::less<size_t> > index_card_mapping(powerset_btree<N, T>& tree){
			return tree.elements_indices_by_ascending_cardinality();
		}

		static inline std::map<size_t, std::vector<set_N_value<N, T>* >, std::less<size_t> > card_mapping(powerset_btree<N, T>& tree){
			return tree._elements_by_ascending_cardinality();
		}

		static inline std::map<size_t, std::vector<size_t>, std::greater<size_t> > index_card_mapping_dual(powerset_btree<N, T>& tree){
			return tree.elements_indices_by_descending_cardinality();
		}

		static inline std::map<size_t, std::vector<set_N_value<N, T>* >, std::greater<size_t> > card_mapping_dual(powerset_btree<N, T>& tree){
			return tree._elements_by_descending_cardinality();
		}

		static inline void elements_related_to(const powerset_btree<N, T>& tree, const subset& set, std::vector<set_N_value<N, T> const * >& elements){
			tree.subsets_of(set, elements);
		}

		static inline void addresses_related_to(const powerset_btree<N, set_N_value<N, T>* >& tree, const subset& set, std::vector<set_N_value<N, set_N_value<N, T>* > const * >& addresses){
			tree.subsets_of(set, addresses);
		}

		static inline void elements_dually_related_to(const powerset_btree<N, T>& tree, const subset& set, std::vector<set_N_value<N, T> const * >& elements){
			tree.supersets_of(set, elements);
		}

		static inline void sets_dually_related_to(const powerset_btree<N, bool>& tree, const subset& set, std::vector<size_t>& sets){
			tree.supersets_of(set, sets);
		}

		static inline subset absorbing_set_for_operation(){
			return ~absorbing_set_for_dual_operation();
		}

		static inline subset absorbing_set_for_dual_operation(){
			return 0;
		}
	};

	template<typename T = float>
	struct addition {

//		static constexpr T neutral_value = 0;
		static inline T neutral_value() {
			return 0;
		}

		static inline void zeta(T& a, const T& b) {
			a += b;
		}

		static inline void mobius(T& a, const T& b) {
			a -= b;
		}
	};

	template<typename T = float>
	struct multiplication {

//		static constexpr T neutral_value = 1;
		static inline T neutral_value() {
			return 1;
		}

		static inline void zeta(T& a, const T& b) {
			a *= b;
		}

		static inline void mobius(T& a, const T& b) {
			a /= b;
		}
	};

	template<typename T = float>
	struct zeta_additive_operation {

		static constexpr T neutral_value = 0;

		void operator()(T& a, const T& b) const {
			a += b;
		}
	};

	template<typename T = float>
	struct mobius_additive_operation {

		static constexpr T neutral_value = 0;

		void operator()(T& a, const T& b) const {
			a -= b;
		}
	};

	template<typename T = float>
	struct zeta_multiplicative_operation {

		static constexpr T neutral_value = 1;

		void operator()(T& a, const T& b) const {
			a *= b;
		}
	};

	template<typename T = float>
	struct mobius_multiplicative_operation {

		static constexpr T neutral_value = 1;

		void operator()(T& a, const T& b) const {
			a /= b;
		}
	};

	template<class inclusion, class operation_type, size_t N, typename T = float>
	struct zeta_tranformation {

//		static constexpr T neutral_value = operation_type::neutral_value;
		static inline T neutral_value() {
			return operation_type::neutral_value();
		}

		static inline void value_inplace_operation(T& a, const T& b) {
			operation_type::zeta(a, b);
		}

		static inline size_t subgraph_index(const size_t& n, const size_t& i){
			return n - i;
		}

		static inline size_t subgraph_dual_index(const size_t& n, const size_t& i){
			return i;
		}

		//template<class inclusion, class value_inplace_operation>
		static void execute_direct_transformation(
				const powerset_btree<N, T>& support,
				powerset_btree<N, T>& focal_points_tree
		){
			//clock_t t;
			//t = clock();
//			powerset_btree<N, T> focal_points_N_initial_values(focal_points_tree);
			/*
			std::function<T(const subset&)> elements_related_to;
			if(order_relation == order_relation_t::subset){
				elements_related_to = focal_points_N_initial_values.subsets_of;
			}else{
				elements_related_to = focal_points_N_initial_values.supersets_of;
			}*/
			const std::vector<set_N_value<N, T>* >& focal_points = focal_points_tree._elements();
			T val;
			std::vector<set_N_value<N, T> const * > elements;
			elements.reserve(support.size());
			for (size_t i = 0; i < focal_points.size(); ++i) {
				val = focal_points[i]->value;
				inclusion::elements_related_to(
					support,
					focal_points[i]->set,
					elements
				);
				for (size_t ii = 0; ii < elements.size(); ++ii) {
					if (elements[ii]->set != focal_points[i]->set){
						value_inplace_operation(val, elements[ii]->value);
					}
				}
				focal_points[i]->value = val;
				elements.clear();
			}
			//t = clock() - t;
			//std::cout << (((float) t)/CLOCKS_PER_SEC) << std::endl;
		}

//		static inline void consonant_operation(T& value, T& preceding_value){
//			value_inplace_operation(value, preceding_value);
//			preceding_value = value;
//		}

		//template<class inclusion, class value_inplace_operation>
		static void execute_consonant_transformation(
				powerset_btree<N, T>& support
		) {
			T preceding_value = operation_type::neutral_value();
			/*
			std::binary_function<size_t, size_t, bool> comp;
			if(order_relation == order_relation_t::subset){
				comp = std::less<size_t>();
			}else{
				comp = std::greater<size_t>();
			}*/
//			const std::binary_function<size_t, size_t, bool>& comp = card_map_comparator();
//			std::map<size_t, std::vector<set_N_value<N, T>* >, comp > support_card_map = support.elements_by_set_cardinality(comp);
			const auto& support_card_map = inclusion::card_mapping(support);

			for (const auto& c_support_elements : support_card_map) {
				const std::vector<set_N_value<N, T>* >& support_elements = c_support_elements.second;
//				consonant_operation(support_elements[0]->value, preceding_value);
				value_inplace_operation(support_elements[0]->value, preceding_value);
				preceding_value = support_elements[0]->value;
				/*
				value = range_binary_operator(set_value->value, preceding_value);
				if (transform_type == transform_type_t::zeta){
					set_value->value = value;
					preceding_value = set_value->value;
				}else{
					preceding_value = set_value->value;
					set_value->value = value;
				}*/
			}
		}

	};

	template<class inclusion, class operation_type, size_t N, typename T = float>
	struct mobius_tranformation {

//		static constexpr T neutral_value = operation_type::neutral_value;
		static inline T neutral_value() {
			return operation_type::neutral_value();
		}

		static inline void value_inplace_operation(T& a, const T& b) {
			operation_type::mobius(a, b);
		}

		static inline size_t subgraph_index(const size_t& n, const size_t& i){
			return i;
		}

		static inline size_t subgraph_dual_index(const size_t& n, const size_t& i){
			return n - i;
		}

		//template<class inclusion, class value_inplace_operation>
		static void execute_direct_transformation(
				const powerset_btree<N, T>& initial_focal_points_tree,	// useless here but maintains the isomorphism with execute_direct_transformation from zeta_transformation
				powerset_btree<N, T>& focal_points_tree
		){
			//clock_t t;
			//t = clock();
			/*
			std::binary_function<size_t, size_t, bool> comp;
			std::function<T(const subset&)> elements_related_to;
			if(order_relation == order_relation_t::subset){
				comp = std::less<size_t>();
				elements_related_to = focal_points_tree.subsets_of;
			}else{
				comp = std::greater<size_t>();
				elements_related_to = focal_points_tree.supersets_of;
			}
			*/
//			const std::binary_function<size_t, size_t, bool>& comp = this->card_map_comparator();
//			std::map<size_t, std::vector<set_N_value<N, T>* >, comp > focal_points_card_map = focal_points_tree.elements_by_set_cardinality(comp);
//			std::cout << "Direct Mobius transformation\n";
			const auto& focal_points_card_map = inclusion::card_mapping(focal_points_tree);
			T val;
			std::vector<set_N_value<N, T> const * > elements;
			elements.reserve(focal_points_tree.size());
			for (const auto& c_focal_points : focal_points_card_map) {
				const std::vector<set_N_value<N, T>* >& focal_points = c_focal_points.second;
				for (size_t i = 0; i < focal_points.size(); ++i) {
					val = focal_points[i]->value;
					inclusion::elements_related_to(focal_points_tree, focal_points[i]->set, elements);
					for (size_t ii = 0; ii < elements.size(); ++ii) {
						if (elements[ii]->set != focal_points[i]->set){
							value_inplace_operation(val, elements[ii]->value);
						}
					}
					focal_points[i]->value = val;
					elements.clear();
				}
			}
			//t = clock() - t;
			//std::cout << (((float) t)/CLOCKS_PER_SEC) << std::endl;
		}

//		static inline void consonant_operation(T& value, T& preceding_value){
//			T old_value = value;
//			value_inplace_operation(value, preceding_value);
//			preceding_value = old_value;
//		}

		//template<class inclusion, class value_inplace_operation>
		static void execute_consonant_transformation(
				powerset_btree<N, T>& support
		) {
			T preceding_value = operation_type::neutral_value();
			/*
			std::binary_function<size_t, size_t, bool> comp;
			if(order_relation == order_relation_t::subset){
				comp = std::less<size_t>();
			}else{
				comp = std::greater<size_t>();
			}*/
//			const std::binary_function<size_t, size_t, bool>& comp = card_map_comparator();
//			std::map<size_t, std::vector<set_N_value<N, T>* >, comp > support_card_map = support.elements_by_set_cardinality(comp);
			const auto& support_card_map = inclusion::card_mapping(support);
			T old_value;
			for (const auto& c_support_elements : support_card_map) {
				const std::vector<set_N_value<N, T>* >& support_elements = c_support_elements.second;
//				consonant_operation(support_elements[0]->value, preceding_value);
				old_value = support_elements[0]->value;
				value_inplace_operation(support_elements[0]->value, preceding_value);
				preceding_value = old_value;
				/*
				value = range_binary_operator(set_value->value, preceding_value);
				if (transform_type == transform_type_t::zeta){
					set_value->value = value;
					preceding_value = set_value->value;
				}else{
					preceding_value = set_value->value;
					set_value->value = value;
				}*/
			}
		}
	};

	template<size_t N, typename T = float>
	static bool consonance_check(
		const powerset_btree<N, T>& support
	){
		// Cardinality map of focal sets in original_structure
		std::map<size_t, std::vector<set_N_value<N, T> const * >, std::less<size_t> > support_card_map = support.elements_by_ascending_cardinality();

		typename sample_space<N>::subset previous_set = 0;
		for (const auto& c_support_elements : support_card_map) {
			const std::vector<set_N_value<N, T> const * >& support_elements = c_support_elements.second;
			// if there are at least two elements with the same cardinality
			// OR if a lower rank set is not a subset of the current set,
			// then structure cannot be consonant
			if (support_elements.size() > 1){
				return false;
			}
			if((previous_set & support_elements[0]->set) != previous_set){
				return false;
			}
			previous_set = support_elements[0]->set;
		}
		return true;
	}

	template<class inclusion, class transformation, size_t N, typename T = float>
	struct efficient_mobius_inversion {

		typedef typename sample_space<N>::subset subset;

		static bool try_linear_focal_points_computation (
			const powerset_btree<N, T>& support,
			powerset_btree<N, T>& focal_points_tree
		){
			std::cout << "Starting linear analysis\n";
			const subset& absorbing_set = inclusion::absorbing_set_for_operation();
			subset dual_absorbing_set = absorbing_set;
			const std::vector<set_N_value<N, T> const * >& elements = support.elements();
			for (size_t i = 0; i < elements.size(); ++i){
				dual_absorbing_set = inclusion::set_dual_operation(dual_absorbing_set, elements[i]->set);
				focal_points_tree.insert(elements[i]->set, elements[i]->value);
			}

			subset U = elements[0]->set;

			size_t inserted_node_index;

			for (size_t i = 1; i < elements.size(); ++i) {
				const subset& A = elements[i]->set;
				if (A != dual_absorbing_set){
					const subset& I = inclusion::set_operation(U, A);
//					std::cout << I << " = set_operation(" << U << ", "<< A << ")\n";
					inserted_node_index = focal_points_tree.insert_set_if_smaller_than(
						I,
						transformation::neutral_value(),
						1
					);
					if (inserted_node_index >= focal_points_tree.number_of_nodes()){
						return false;
					}
					U = inclusion::set_dual_operation(U, A);
				}
			}
			focal_points_tree.insert_set(
				absorbing_set,
				transformation::neutral_value()
			);
//			std::cout << "END PRINTING\n";
			return true;
		}

		static void execute_direct_transformation(
			const powerset_btree<N, T>& support,
			powerset_btree<N, T>& focal_points_tree
		) {
			std::cout << "Executing direct\n";
			transformation::execute_direct_transformation(support, focal_points_tree);
			std::cout << "--- Done\n";
		}

		static void execute_consonant_transformation(
			powerset_btree<N, T>& support
		) {
			std::cout << "Executing consonant\n";
			transformation::execute_consonant_transformation(support);
			std::cout << "--- Done\n";
		}

		static void execute_FMT(
				std::vector<T>& transform
		) {
			std::cout << "Executing FMT\n";
//			if(transform.size() != pow(2, N)){
//				std::cerr << "\nThe size of the given vector is not 2^N, where N is the given number of outcomes.\n";
////					return powerset_values;
//			}
//				std::vector<T> transform(powerset_values);
			if (transform.size() == 0)
				return;
			size_t n = (size_t) log2(transform.size());
			size_t sub_powerset_size, sub_powerset_dual_size, index;
			for (size_t i = 1; i <= n; ++i){

				sub_powerset_size = 1 << i;
				for (size_t j = 1; j <= sub_powerset_size; j += 2){

					sub_powerset_dual_size = 1 << (n - i);
					for (size_t k = 0; k <= sub_powerset_dual_size-1; ++k){
						index = (j-inclusion::target_index_offset_FMT()) * sub_powerset_dual_size + k;
//						std::cout << "Pour " << transform[index] << "\n";
//						std::cout << "index " << index << " / " << (j-inclusion::source_index_offset_FMT()) * sub_powerset_dual_size + k << "\n";
						transformation::value_inplace_operation(transform[index], transform[(j-inclusion::source_index_offset_FMT()) * sub_powerset_dual_size + k]);
//						std::cout << "= " << transform[index] << "\n";
					}
				}
			}
			std::cout << "--- Done\n";
//				return transform;
		}

		static void build_core_reduced_powerset(
				const powerset_btree<N, T>& support,
				std::vector<subset >& iota_sequence,
				powerset_btree<N, T>& core_reduced_powerset
		) {
			std::cout << "Building core reduced powerset\n";
			const std::vector<set_N_value<N, T> const * >& support_elements = support.elements();
			subset core = 0;
			for(size_t i = 0; i < support_elements.size(); ++i){
				core |= support_elements[i]->set;
			}
//			std::vector<subset > iota_sequence;
			subset singleton = 1;
			for(size_t i = 0; i < N; ++i){
				if ((singleton & core) != 0){
					iota_sequence.emplace_back(singleton);
				}
				singleton <<= 1;
			}

			std::vector<subset > focal_points;
			focal_points.reserve(pow(2, iota_sequence.size()));

			for (size_t i = 0; i < support_elements.size(); ++i){
				focal_points.emplace_back(support_elements[i]->set);
//				core_reduced_powerset.insert(support_elements[i]->set, support_elements[i]->value);
			}

			singleton = 0;
			size_t newly_inserted_index = support[singleton];
			core_reduced_powerset.insert_set(
				singleton,
				transformation::neutral_value()
			);
			if (newly_inserted_index < support.number_of_nodes()){
				focal_points.emplace_back(singleton);
			}

			for (size_t i = 0; i < iota_sequence.size(); ++i) {
				size_t end = focal_points.size();
				for (size_t ii = 0; ii < end; ++ii) {
					const subset& new_set = iota_sequence[i] | focal_points[ii];
					newly_inserted_index = core_reduced_powerset.insert_set(
						new_set,
						transformation::neutral_value()
					);
					if (newly_inserted_index < core_reduced_powerset.number_of_nodes()){
						focal_points.emplace_back(new_set);
					}
				}
			}
			std::cout << "--- Built\n";
//			return iota_sequence;
		}

		static void execute_FMT_reduced_to_core(
				powerset_btree<N, T>& core_reduced_powerset,
				const std::vector<subset >& iota_sequence
		) {
			std::cout << "Executing FMT reduced to core\n";
			std::vector<set_N_value<N, T> const * > powerset_elements = core_reduced_powerset.elements();
			for (size_t i = 0; i < iota_sequence.size(); ++i){
				for (size_t e = 0; e < powerset_elements.size(); ++e){
					subset set_B = inclusion::set_operation(powerset_elements[e]->set, inclusion::FMT_target(iota_sequence[i]));

					if (set_B != powerset_elements[e]->set){
						set_N_value<N, T>& B = core_reduced_powerset._node(core_reduced_powerset[set_B]);
						transformation::value_inplace_operation(B.value, powerset_elements[e]->value);
					}
				}
			}
			std::cout << "--- Done\n";
		}

		/*
		 * In this function, this->iota_sequence is supposed to contain all regular iota elements (i.e. join-irreducible of this lattice).
		 */
		static void build_truncated_lattice_support(
				const powerset_btree<N, T>& support,
				const std::vector<subset >& iota_sequence,
				powerset_btree<N, T>& truncated_lattice_support
		) {
			std::cout << "Building truncated lattice support\n";
			std::vector<subset > focal_points;
			focal_points.reserve(support.size());

			const std::vector<set_N_value<N, T> const * >& elements = support.elements();
			for (size_t i = 0; i < elements.size(); ++i){
				focal_points.emplace_back(elements[i]->set);
//				truncated_lattice_support.insert(elements[i]->set, elements[i]->value);
			}
			update_truncated_lattice_support(
				iota_sequence,
				truncated_lattice_support,
				focal_points,
				transformation::neutral_value()
			);
			std::cout << "--- Built\n";
		}

		static void update_truncated_lattice_support(
				const std::vector<subset >& iota_sequence,
				powerset_btree<N, T>& truncated_lattice_support,
				std::vector<subset > focal_points,
				T default_value
		) {
			std::cout << "Updating truncated lattice support\n";
			size_t newly_inserted_index;

			for (size_t i = 0; i < iota_sequence.size(); ++i) {
				size_t end = focal_points.size();
				for (size_t ii = 0; ii < end; ++ii) {
					const subset& new_set = inclusion::set_operation(iota_sequence[i], focal_points[ii]);
					newly_inserted_index = truncated_lattice_support.insert_set(
						new_set,
						default_value
					);
					if (newly_inserted_index < truncated_lattice_support.number_of_nodes()){
						focal_points.emplace_back(new_set);
					}
				}
			}
//			truncated_lattice_support.print();
		}

		static void execute_EMT_with_lattice(
				powerset_btree<N, T>& lattice_support,
				const std::vector<subset >& iota_sequence
		) {
			std::cout << "Executing EMT with lattice\n";
			const std::vector<set_N_value<N, T> const * >& lattice_support_elements = lattice_support.elements();
//			lattice_support.print();
//			std::cout << "IOTA ELEMENTS:\n";
//			std::cout << iota_sequence[0] << "\n";
//			for (size_t i = 1; i < iota_sequence.size(); ++i){
//				std::cout << iota_sequence[i] << "\n";
//			}

			size_t nb_iota = iota_sequence.size()-1;
//			size_t iota_index;
			for (size_t i = 0; i < iota_sequence.size(); ++i){
//				if (transform_type == transform_type_t::Mobius)
//					iota_index = i;
//				else
//					iota_index = nb_iota - i;
				const size_t& iota_index = transformation::subgraph_index(nb_iota, i);
//				std::cout << "Iota " << iota_sequence[iota_index] << " and iota_index " << iota_index << " generate:\n";

				for (size_t e = 0; e < lattice_support_elements.size(); ++e){
//					const subset& set_B = lattice_support_elements[e]->set | iota_sequence[iota_index];
					const subset& set_B = inclusion::set_operation(lattice_support_elements[e]->set, iota_sequence[iota_index]);
//					std::cout << "\tProxy " << set_B << " with element " << lattice_support_elements[e]->set << "\n";

					if (set_B != lattice_support_elements[e]->set){
						size_t B = lattice_support[set_B];
						if (B < lattice_support.number_of_nodes()){
//							std::cout << "\tSet = " << lattice_support._node(B).set << ", Target = " << lattice_support_elements[e]->set << "\n";
//							std::cout << "\t\t" << B->set << " - " << lattice_support_elements[e]->set << "\n";
							transformation::value_inplace_operation(lattice_support._node(B).value, lattice_support_elements[e]->value);
						}
					}
				}
			}
			std::cout << "--- Done\n";
		}

		static void build_semilattice_support(
				const powerset_btree<N, T>& support,
				powerset_btree<N, T>& focal_points_tree
		) {
			std::cout << "Building semilattice support\n";
			std::vector<subset > focal_points;
			const std::vector<set_N_value<N, T> const * >& support_elements = support.elements();
			focal_points.reserve(N * support.size());

			for (size_t i = 0; i < support_elements.size(); ++i){
				focal_points.emplace_back(support_elements[i]->set);
//				focal_points_tree.insert(support_elements[i]->set, support_elements[i]->value);
			}

			size_t newly_inserted_index;

			for (size_t i = 0; i < support_elements.size(); ++i){
				size_t end = focal_points.size();
				for (size_t ii = i+1; ii < end; ++ii){
					const subset& focal_point = inclusion::set_operation(support_elements[i]->set, focal_points[ii]);
					newly_inserted_index = focal_points_tree.insert_set(
						focal_point,
						transformation::neutral_value()
					);
					if (newly_inserted_index < focal_points_tree.number_of_nodes()){
						focal_points.emplace_back(focal_point);
					}
				}
			}
			std::cout << "--- Built\n";
		}

		static void update_semilattice_support(
				const powerset_btree<N, T>& support,
				powerset_btree<N, T>& focal_points_tree,
				std::vector<subset > focal_points,
				T default_value
		) {
			std::cout << "Updating semilattice support\n";
			const std::vector<set_N_value<N, T> const * >& support_elements = support.elements();

			size_t newly_inserted_index;

			for (size_t i = 0; i < support_elements.size(); ++i){
				size_t end = focal_points.size();
				for (size_t ii = 0; ii < end; ++ii){
					const subset& focal_point = inclusion::set_operation(support_elements[i]->set, focal_points[ii]);
					newly_inserted_index = focal_points_tree.insert_set(
						focal_point,
						default_value
					);
					if (newly_inserted_index < focal_points_tree.number_of_nodes()){
						focal_points.emplace_back(focal_point);
					}
				}
			}
		}

		static void build_bridge_map(
				std::unordered_map<subset, size_t>& bridge_map,
				powerset_btree<N, T>& focal_points_tree,
				const std::vector<subset >& iota_sequence
		) {
			std::cout << "Building bridge map\n";
			const std::vector<size_t>& focal_points = focal_points_tree.elements_indices();

			powerset_btree<N, bool> proxies_missing_targets(iota_sequence.size());
			size_t i;
			for (size_t e = 0; e < focal_points.size(); ++e){
				const subset& set = focal_points_tree.get_node(focal_points[e]).set;
				bridge_map[set] = focal_points[e];
//				std::cout << e << "\n\t";
				for (i = 0; i < iota_sequence.size(); ++i) {
//					std::cout << i << " ";
//					const subset& set = focal_points[e]->set | iota_sequence[i];
					const subset& proxy = inclusion::set_dual_operation(set, iota_sequence[i]);
//					std::cout << "Proxy = " << proxy << "\n";
					if (focal_points_tree[proxy] >= focal_points_tree.number_of_nodes()){
						if (proxies_missing_targets.insert_set(proxy, true) < proxies_missing_targets.number_of_nodes()){
//							std::cout << "Missing\n";
						}
					}
				}
//				std::cout << "\n";
			}

			if (proxies_missing_targets.size() == 0){
				return;
			}

//			const std::binary_function<size_t, size_t, bool>& comp = std::less<size_t>();
//			const std::binary_function<size_t, size_t, bool>& comp = card_map_dual_comparator();
//			std::map<size_t, std::vector<set_N_value<N, T>* >, comp > focal_points_card_map = focal_points_tree.elements_by_set_cardinality(comp);
			const auto& focal_points_card_map = inclusion::index_card_mapping_dual(focal_points_tree);
			//clock_t t;
			//t = clock();
			std::vector<size_t> proxies;
			proxies.reserve(proxies_missing_targets.size());
			size_t nb_targets_to_find, cc = 0;
			for (const auto& c_focal_points : focal_points_card_map) {
				for (; cc < c_focal_points.first; ++cc){
					nb_targets_to_find += proxies_missing_targets.get_nb_sets_of_cardinality(cc);
				}
				if (nb_targets_to_find > 0){
					const std::vector<size_t>& focal_points_with_same_size = c_focal_points.second;

					for (size_t i = 0; i < focal_points_with_same_size.size(); ++i){
						const subset& set = focal_points_tree.get_node(focal_points_with_same_size[i]).set;
//						const std::vector<set_N_value<bool, N>* >& subsets = proxies_missing_targets.subsets_of(focal_points_with_same_size[i]->set);
						inclusion::sets_dually_related_to(proxies_missing_targets, set, proxies);

						for (size_t s = 0; s < proxies.size(); ++s){
//							std::cout << "Proxy = " << proxies_missing_targets.get_node(proxies[s]).set << "\n";
							bridge_map[proxies_missing_targets.get_node(proxies[s]).set] = focal_points_with_same_size[i];
//							std::cout << "Target access = " << bridge_map[proxies_missing_targets.get_node(proxies[s]).set] << "\n";
							proxies_missing_targets.nullify(proxies[s]);
							--nb_targets_to_find;
						}
						if (nb_targets_to_find == 0){
							break;
						}
						proxies.clear();
					}
				}
			}
			const std::vector<set_N_value<N, bool> const * > elements = proxies_missing_targets.elements();
			for (size_t i = 0; i < elements.size(); ++i){
				bridge_map[elements[i]->set] = focal_points_tree.number_of_nodes();
			}
			//t = clock() - t;
			//std::cout << "missing proxy search : " << (((float) t)/CLOCKS_PER_SEC) << " ";
			std::cout << "--- Built\n";
		}


		static void execute_EMT_with_semilattice(
				powerset_btree<N, T>& focal_points_tree,
				const std::vector<subset >& iota_sequence,
				std::unordered_map<subset, size_t> bridge_map
		) {
			//clock_t t;
			std::cout << "Executing EMT with semilattice\n";

			if (iota_sequence.size() == 0)
				return;

			std::vector<subset > sync_sequence;
			sync_sequence.reserve(iota_sequence.size());
			sync_sequence.emplace_back(iota_sequence[0]);
			for (size_t i = 1; i < iota_sequence.size(); ++i){
//				sync_sequence.emplace_back(sync_sequence[i-1] | iota_sequence[i]);
				sync_sequence.emplace_back(inclusion::set_dual_operation(sync_sequence[i-1], iota_sequence[i]));
			}
//			std::unordered_map<subset, set_N_value<N, T> const * > bridge_map;
//			build_bridge_map(bridge_map, focal_points_tree, iota_sequence);
//			std::cout << "Bridge built\n";

			//t = clock();
			size_t nb_iota = iota_sequence.size()-1;
//			size_t iota_index;
			const std::vector<set_N_value<N, T>* >& focal_points = focal_points_tree._elements();
			for (size_t i = 0; i < iota_sequence.size(); ++i){
//				if (transform_type == transform_type_t::zeta)
//					iota_index = i;
//				else
//					iota_index = nb_iota - i;
				const size_t& iota_index = transformation::subgraph_dual_index(nb_iota, i);
//				std::cout << "IOTA = " << iota_sequence[iota_index] << "\n";

				for (size_t e = 0; e < focal_points.size(); ++e){
//					const subset& proxy = focal_points[e]->set | iota_sequence[iota_index];
					const subset& proxy = inclusion::set_dual_operation(focal_points[e]->set, iota_sequence[iota_index]);

					if (focal_points[e]->set != proxy){
						const size_t coupled_set_index = bridge_map[proxy];

//						if (coupled_set && FOD<N>::is_subset_of(coupled_set->set, focal_points[e]->set | sync_sequence[iota_index])){
						if (coupled_set_index < focal_points_tree.number_of_nodes()){
							const set_N_value<N, T>& coupled_set = focal_points_tree.get_node(coupled_set_index);
							const subset& sync_set = inclusion::set_dual_operation(focal_points[e]->set, sync_sequence[iota_index]);
							if(inclusion::set_operation(coupled_set.set, sync_set) == coupled_set.set){
//								std::cout << "\tSet = " << focal_points[e]->set << ", Proxy = " << proxy << ", Target = " << coupled_set_index << "\n";
								transformation::value_inplace_operation(focal_points[e]->value, coupled_set.value);
							}
						}
					}
				}
			}
			//t = clock() - t;
			//std::cout << (((float) t)/CLOCKS_PER_SEC) << std::endl;
			std::cout << "--- Done\n";
		}

		static scheme_type_t autoset_and_build(
				const powerset_btree<N, T>& support,
				powerset_btree<N, T>& focal_points_tree,
				std::vector<subset >& iota_sequence,
				std::unordered_map<subset, size_t>& bridge_map
		){
			// check if original_structure is almost Bayesian
			const bool& is_almost_bayesian = try_linear_focal_points_computation (
				support,
				focal_points_tree
			);

			if(is_almost_bayesian){
//				std::cout << "Almost Bayesian\n";
				return scheme_type_t::direct;
			}else{
				const bool& is_consonant = consonance_check<N, T>(support);

				if(is_consonant){
					focal_points_tree = support;
					return scheme_type_t::consonant;
				}else{

					if(support.size() < 3 * N){

						build_semilattice_support(
							support,
							focal_points_tree
						);

						if(focal_points_tree.size() < 3 * N){
							return scheme_type_t::direct;
						}else{

							inclusion::compute_iota_elements_dual(
								support,
								iota_sequence
							);
							build_bridge_map(bridge_map, focal_points_tree, iota_sequence);
							return scheme_type_t::semilattice;
						}
					}else{
						inclusion::compute_iota_elements(
							support,
							iota_sequence
						);

						build_truncated_lattice_support(
							support,
							iota_sequence,
							focal_points_tree
						);
						return scheme_type_t::lattice;
					}
				}
			}
		}

		static void execute(
			const powerset_btree<N, T>& support,
			powerset_btree<N, T>& focal_points_tree,
			const std::vector<subset >& iota_sequence,
			const std::unordered_map<subset, size_t>& bridge_map,
			const scheme_type_t& scheme_type
		) {

			switch (scheme_type) {
				case scheme_type_t::direct:
					execute_direct_transformation(
						support,
						focal_points_tree
					);
					break;

				case scheme_type_t::consonant:
					execute_consonant_transformation(
						focal_points_tree
					);
					break;

				case scheme_type_t::semilattice:
					execute_EMT_with_semilattice(
						focal_points_tree,
						iota_sequence,
						bridge_map
					);
					break;

				case scheme_type_t::lattice:
					execute_EMT_with_lattice(
						focal_points_tree,
						iota_sequence
					);

					break;

				case scheme_type_t::reduced_FMT:
					execute_FMT_reduced_to_core(
						focal_points_tree,
						iota_sequence
					);

					break;

				default:
					break;
			}
		}

		static void direct_transformation(
				const powerset_btree<N, T>& support,
				powerset_btree<N, T>& focal_points_tree,
				scheme_type_t& scheme_type
		) {
//			powerset_btree<N, T> focal_points_tree(support.size());
			scheme_type = scheme_type_t::direct;
			// check if original_structure is almost Bayesian
			const bool& is_almost_bayesian = try_linear_focal_points_computation (
				support,
				focal_points_tree
			);

			if (!is_almost_bayesian){
				const bool& is_consonant = consonance_check<N, T>(support);
				if(is_consonant){
					focal_points_tree = support;
					scheme_type = scheme_type_t::consonant;
				}else{
					build_semilattice_support(
						support,
						focal_points_tree
					);
				}
			}
			if(scheme_type == scheme_type_t::direct){
				execute_direct_transformation(
					support,
					focal_points_tree
				);
			}else{
				execute_consonant_transformation(
					focal_points_tree
				);
			}
		}

		static void FMT(
				const std::vector<T>& support,
				std::vector<T>& transform
		) {
			transform = support;

			execute_FMT(
				transform
			);
		}

//		static void FMT_reduced_to_core(
//				const T& support,
//				T& core_reduced_powerset,
//				std::vector<subset >& iota_sequence
//		) {
//			iota_sequence = build_core_reduced_powerset(support, core_reduced_powerset);
//			execute_FMT_reduced_to_core(
//					core_reduced_powerset,
//					iota_sequence
//			);
//		}

		static void FMT_reduced_to_core(
				const powerset_btree<N, T>& support,
				powerset_btree<N, T>& core_reduced_powerset,
				std::vector<subset >& iota_sequence
		) {
			build_core_reduced_powerset(support, iota_sequence, core_reduced_powerset);
			execute_FMT_reduced_to_core(
					core_reduced_powerset,
					iota_sequence
			);
		}

		static void EMT_with_semilattice(
			const powerset_btree<N, T>& support,
			powerset_btree<N, T>& focal_points_tree,
			std::vector<subset >& iota_sequence,
			std::unordered_map<subset, size_t>& bridge_map
		) {
//			powerset_btree<N, T> focal_points_tree(support.size());
			build_semilattice_support(
				support,
				focal_points_tree
			);
//			const std::vector<subset >&
			inclusion::compute_iota_elements_dual(
				support,
				iota_sequence
			);
			build_bridge_map(bridge_map, focal_points_tree, iota_sequence);
			execute_EMT_with_semilattice(
					focal_points_tree,
					iota_sequence,
					bridge_map
			);
		}

		static void EMT_with_lattice(
			const powerset_btree<N, T>& support,
			powerset_btree<N, T>& truncated_lattice_support,
			std::vector<subset >& iota_sequence
		) {
//			powerset_btree<N, T> truncated_lattice_support(support.size());
//			const std::vector<subset >&
			inclusion::compute_iota_elements(
				support,
				iota_sequence
			);
			build_truncated_lattice_support(
				support,
				iota_sequence,
				truncated_lattice_support
			);
			execute_EMT_with_lattice(
					truncated_lattice_support,
					iota_sequence
			);
		}
	};
}		// namespace efficient_DST

#endif // EFFICIENT_DST_MOBIUS_INVERSION_TEMPLATE_HPP
