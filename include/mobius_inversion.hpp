#ifndef EFFICIENT_DST_EMT_HPP
#define EFFICIENT_DST_EMT_HPP

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
	class mobius_inversion {
	public:

		virtual ~mobius_inversion()
		{}

	protected:

		T neutral_value;

		static powerset_btree<T, N> FMT(
				const powerset_btree<T, N>& support
		) {
			const std::vector<set_N_value<T, N>* >& support_elements = support.elements();
			std::bitset<N> core = 0;
			for(size_t i = 0; i < support_elements.size(); ++i){
				core |= support_elements[i]->set;
			}

//			size_t reduced_powerset_size = pow(2, core.count());
			powerset_btree<T, N> transform(support.get_FOD(), pow(2, core.count()));
			transform.copy(support);
//			powerset_btree<T, N> transform(support.get_FOD(), reduced_powerset_size);
			std::vector<set_N_value<T, N>* > focal_points = transform.elements();
//			std::unordered_map<std::bitset<N>, set_N_value<T, N>* > focal_points_map;
//			focal_points.reserve(reduced_powerset_size);
//			focal_points_map.reserve(reduced_powerset_size);
//
//			for (size_t i = 0; i < support_elements.size(); ++i){
//				focal_points.emplace_back(support_elements[i]->set);
//				focal_points_map.emplace(support_elements[i]->set, transform.insert(support_elements[i]->set, support_elements[i]->value));
//			}

			size_t i = core._Find_first();
			while(i < N){
				for (size_t e = 0; e < focal_points.size(); ++e) {
					std::bitset<N> new_set = (const std::bitset<N>&) focal_points[e]->set;
					new_set.set(i, true);
					set_N_value<T, N>* node = transform.find(new_set);
					if (!node){
						transform.insert(new_set, neutral_value);
					}
				}
				i = core._Find_next(i);
			}

			std::vector<set_N_value<T, N>* > powerset_elements = transform.elements();

			i = core._Find_first();
			while(i < N){
				for (size_t e = 0; e < powerset_elements.size(); ++e){
					std::bitset<N> set_B = (const std::bitset<N>&) powerset_elements[e]->set;
					set_B.set(i, true);

					if (set_B != powerset_elements[e]->set){
						set_N_value<T, N>* B = transform.find(set_B);

						if (B){
							value_inplace_operation_FMT(B->second->value, powerset_elements[e]->value);
						}
					}
				}
				i = core._Find_next(i);
			}
			return transform;
		}

		static inline void value_inplace_operation_FMT(T& a, T& b);

		static inline void down_inclusion_value_inplace_operation_FMT(T& a, T& b){
			value_inplace_operation(a, b);
		}

		static inline void up_inclusion_value_inplace_operation_FMT(T& a, T& b){
			value_inplace_operation(b, a);
		}


		static std::vector<T> FMT_vectorized(
				const std::vector<T>& powerset_values,
		) {
			if(powerset_values.size() != pow(2, N)){
				std::cerr << "\nThe size of the given vector is not 2^N, where N is the given size corresponding to the considered FOD.\n";
				return powerset_values;
			}

			std::vector<T> transform(powerset_values);
			size_t sub_powerset_size, sub_powerset_dual_size, index;
			for (size_t i = 1; i <= N; ++i){

				sub_powerset_size = pow(2, i);
				for (size_t j = 1; j <= sub_powerset_size; j += 2){

					sub_powerset_dual_size = pow(2, N - i);
					for (size_t k = 0; k <= sub_powerset_dual_size-1; ++k){
						index = (j-target_index_offset_FMT()) * sub_powerset_dual_size + k;
						value_inplace_operation(transform[index], transform[(j-source_index_offset_FMT()) * sub_powerset_dual_size + k]);
					}
				}
			}
			return transform;
		}

		static inline size_t target_index_offset_FMT();

		static inline size_t down_inclusion_target_index_offset_FMT(){
			return 0;
		}

		static inline size_t up_inclusion_target_index_offset_FMT(){
			return 1;
		}

		static inline size_t source_index_offset_FMT();

		static inline size_t down_inclusion_source_index_offset_FMT(){
			return 1;
		}

		static inline size_t up_inclusion_source_index_offset_FMT(){
			return 0;
		}

		static void execute_direct_transformation(powerset_btree<T, N>& focal_points_tree);

		static void execute_direct_Zeta_transformation(
				powerset_btree<T, N>& focal_points_tree
		){
			//clock_t t;
			//t = clock();
			powerset_btree<T, N> focal_points_N_initial_values(focal_points_tree);
			/*
			std::function<T(const std::bitset<N>&)> elements_related_to;
			if(order_relation == order_relation_t::subset){
				elements_related_to = focal_points_N_initial_values.subsets_of;
			}else{
				elements_related_to = focal_points_N_initial_values.supersets_of;
			}*/
			const std::vector<set_N_value<T, N>* >& focal_points = focal_points_N_initial_values.elements();
			T val;
			for (size_t i = 0; i < focal_points.size(); ++i) {
				val = focal_points[i]->value;
				const std::vector<set_N_value<T, N>* >& elements = elements_related_to(focal_points[i]->set);
				DEBUG(std::clog << "Elements related to " << focal_points[i]->set << " : \n";);
				for (size_t ii = 0; ii < elements.size(); ++ii) {
					DEBUG(std::clog << elements[ii]->set << std::endl;);
					if (elements[ii]->set != focal_points[i]->set){
						value_inplace_operation(val, elements[ii]->value);
					}
				}
				focal_points_tree.insert(focal_points[i]->set, val);
			}
			//t = clock() - t;
			//std::cout << (((float) t)/CLOCKS_PER_SEC) << std::endl;
		}

		static void execute_direct_Mobius_transformation(
				powerset_btree<T, N>& focal_points_tree
		){
			//clock_t t;
			//t = clock();
			/*
			std::binary_function<size_t, size_t, bool> comp;
			std::function<T(const std::bitset<N>&)> elements_related_to;
			if(order_relation == order_relation_t::subset){
				comp = std::less<size_t>();
				elements_related_to = focal_points_tree.subsets_of;
			}else{
				comp = std::greater<size_t>();
				elements_related_to = focal_points_tree.supersets_of;
			}
			*/
			const std::binary_function<size_t, size_t, bool>& comp = card_map_comparator();
			std::map<size_t, std::vector<set_N_value<T, N>* >, comp > focal_points_card_map = focal_points_tree.elements_by_set_cardinality(comp);
			T val;
			for (const auto& c_focal_points : focal_points_card_map) {
				const std::vector<std::bitset<N> >& focal_points = c_focal_points.second;
				for (size_t i = 0; i < focal_points.size(); ++i) {
					val = focal_points[i]->value;
					const std::vector<set_N_value<T, N>* >& elements = elements_related_to(focal_points[i]->set);
					for (size_t ii = 0; ii < elements.size(); ++ii) {
						if (elements[ii]->set != focal_points[i]->set){
							value_inplace_operation(val, elements[ii]->value);
						}
					}
					focal_points_tree.insert(focal_points[i]->set, val);
				}
			}
			//t = clock() - t;
			//std::cout << (((float) t)/CLOCKS_PER_SEC) << std::endl;
		}


		void execute_consonant_transformation(
				powerset_btree<T, N>& support,
		) const {
			T value, preceding_value = this->neutral_value;
			/*
			std::binary_function<size_t, size_t, bool> comp;
			if(order_relation == order_relation_t::subset){
				comp = std::less<size_t>();
			}else{
				comp = std::greater<size_t>();
			}*/
			const std::binary_function<size_t, size_t, bool>& comp = card_map_comparator();
			std::map<size_t, std::vector<set_N_value<T, N>* >, comp > support_card_map = support.elements_by_set_cardinality(comp);

			for (const auto& c_support_elements : support_card_map) {
				const std::vector<set_N_value<T, N>* >& support_elements = c_support_elements.second;
				consonant_operation(support_elements[0]->value, preceding_value);
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


		static void build_bridge_map(
				std::unordered_map<std::bitset<N>, set_N_value<T, N>* >& bridge_map,
				const powerset_btree<T, N>& focal_points_tree,
				const std::vector<std::bitset<N> >& iota_sequence
		) {
			const std::vector<set_N_value<T, N>* >& focal_points = focal_points_tree.elements();

			powerset_btree<bool, N> proxies_missing_targets(focal_points_tree.get_FOD(), iota_sequence.size());

			for (size_t i = 0; i < iota_sequence.size(); ++i) {
				for (size_t e = 0; e < focal_points.size(); ++e){
//					const std::bitset<N>& set = focal_points[e]->set | iota_sequence[i];
					const std::bitset<N>& set = set_dual_operation(focal_points[e]->set, iota_sequence[i]);
					if (!focal_points_tree.find(set)){
						proxies_missing_targets.insert(set, true);
					}
				}
			}
			if (proxies_missing_targets.size() == 0){
				return;
			}

//			const std::binary_function<size_t, size_t, bool>& comp = std::less<size_t>();
			const std::binary_function<size_t, size_t, bool>& comp = card_map_dual_comparator();
			std::map<size_t, std::vector<set_N_value<T, N>* >, comp > focal_points_card_map = focal_points_tree.elements_by_set_cardinality(comp);
			//clock_t t;
			//t = clock();
			size_t nb_targets_to_find, cc = 0;
			for (const auto& c_focal_points : focal_points_card_map) {
				for (; cc < c_focal_points.first; ++cc){
					nb_targets_to_find += proxies_missing_targets.get_nb_sets_of_cardinality(cc);
				}
				if (nb_targets_to_find > 0){
					const std::vector<set_N_value<T, N>* >& focal_points_with_same_size = c_focal_points.second;

					for (size_t i = 0; i < focal_points_with_same_size.size(); ++i){
//						const std::vector<set_N_value<bool, N>* >& subsets = proxies_missing_targets.subsets_of(focal_points_with_same_size[i]->set);
						const std::vector<set_N_value<bool, N>* >& proxies = elements_dually_related_to(proxies_missing_targets, focal_points_with_same_size[i]->set);

						for (size_t s = 0; s < proxies.size(); ++s){
							bridge_map[proxies[s]->set] = focal_points_with_same_size[i];
							proxies_missing_targets.nullify(proxies[s]);
							--nb_targets_to_find;
						}
						if (nb_targets_to_find == 0){
							break;
						}
					}
				}
			}
			//t = clock() - t;
			//std::cout << "missing proxy search : " << (((float) t)/CLOCKS_PER_SEC) << " ";
		}


		static void execute_EMT_with_semilattice(
				powerset_btree<T, N>& focal_points_tree,
				const transform_type_t& transform_type,
				const std::vector<std::bitset<N> >& iota_sequence
		) {
			clock_t t;

			if (iota_sequence.size() == 0)
				return;

			std::vector<std::bitset<N> > sync_sequence;
			sync_sequence.reserve(iota_sequence.size());

			sync_sequence.emplace_back(iota_sequence[0]);
			for (size_t i = 1; i < iota_sequence.size(); ++i){
//				sync_sequence.emplace_back(sync_sequence[i-1] | iota_sequence[i]);
				sync_sequence.emplace_back(set_dual_operation(sync_sequence[i-1], iota_sequence[i]));
			}
			std::unordered_map<std::bitset<N>, set_N_value<T, N>* > bridge_map;
			build_bridge_map(bridge_map, focal_points_tree, iota_sequence);

			//t = clock();
			size_t nb_iota = iota_sequence.size()-1;
//			size_t iota_index;
			const std::vector<set_N_value<T, N>* >& focal_points = focal_points_tree.elements();
			for (size_t i = 0; i < iota_sequence.size(); ++i){
//				if (transform_type == transform_type_t::zeta)
//					iota_index = i;
//				else
//					iota_index = nb_iota - i;
				const size_t& iota_index = subgraph_dual_index(nb_iota, i);

				for (size_t e = 0; e < focal_points.size(); ++e){
//					const std::bitset<N>& proxy = focal_points[e]->set | iota_sequence[iota_index];
					const std::bitset<N>& proxy = set_dual_operation(focal_points[e]->set, iota_sequence[iota_index]);

					if (focal_points[e]->set != proxy){
						const set_N_value<T, N>* coupled_set = bridge_map[proxy];

//						if (coupled_set && FOD<N>::is_subset_of(coupled_set->set, focal_points[e]->set | sync_sequence[iota_index])){
						if (coupled_set && FOD<N>::is_subset_of(coupled_set->set, set_dual_operation(focal_points[e]->set, sync_sequence[iota_index]))){
							value_inplace_operation(focal_points[e]->value, coupled_set->value);
						}
					}
				}
			}
			//t = clock() - t;
			//std::cout << (((float) t)/CLOCKS_PER_SEC) << std::endl;
		}


		static void execute_EMT_with_lattice(
				powerset_btree<T, N>& lattice_support,
				const transform_type_t& transform_type,
				const std::vector<std::bitset<N> >& iota_sequence
		) {
			const std::vector<set_N_value<T, N>* >& lattice_support_elements = lattice_support.elements();

			size_t nb_iota = iota_sequence.size()-1;
//			size_t iota_index;
			for (size_t i = 0; i < iota_sequence.size(); ++i){
//				if (transform_type == transform_type_t::Mobius)
//					iota_index = i;
//				else
//					iota_index = nb_iota - i;
				const size_t& iota_index = subgraph_index(nb_iota, i);

				for (size_t e = 0; e < lattice_support_elements.size(); ++e){
//					const std::bitset<N>& set_B = lattice_support_elements[e]->set | iota_sequence[iota_index];
					const std::bitset<N>& set_B = set_operation(lattice_support_elements[e]->set, iota_sequence[iota_index]);

					if (set_B != lattice_support_elements[e]->set){
						set_N_value<T, N>* B = lattice_support.find(set_B);

						if (B){
							value_inplace_operation(B->value, lattice_support_elements[e]->value);
						}
					}
				}
			}
		}


		static inline size_t subgraph_index(const size_t& n, const size_t& i);

		static inline size_t mobius_subgraph_index(const size_t& n, const size_t& i){
			return i;
		}

		static inline size_t zeta_subgraph_index(const size_t& n, const size_t& i){
			return n - i;
		}

		static inline size_t subgraph_dual_index(const size_t& n, const size_t& i);

		static inline size_t mobius_subgraph_dual_index(const size_t& n, const size_t& i){
			return n - i;
		}

		static inline size_t zeta_subgraph_dual_index(const size_t& n, const size_t& i){
			return i;
		}

		static inline void value_inplace_operation(T& a, const T& b);

		static inline void zeta_additive_value_inplace_operation(T& a, const T& b){
			a += b;
		}

		static inline void mobius_additive_value_inplace_operation(T& a, const T& b){
			a -= b;
		}

		static inline void zeta_multiplicative_value_inplace_operation(T& a, const T& b){
			a *= b;
		}

		static inline void mobius_multiplicative_value_inplace_operation(T& a, const T& b){
			a /= b;
		}

		static inline std::bitset<N> set_operation(const std::bitset<N>& a, const std::bitset<N>& b);

		static inline std::bitset<N> down_inclusion_set_operation(const std::bitset<N>& a, const std::bitset<N>& b){
			return a | b;
		}

		static inline std::bitset<N> up_inclusion_set_operation(const std::bitset<N>& a, const std::bitset<N>& b){
			return a & b;
		}

		static inline std::bitset<N> set_dual_operation(const std::bitset<N>& a, const std::bitset<N>& b);

		static inline std::bitset<N> down_inclusion_set_dual_operation(const std::bitset<N>& a, const std::bitset<N>& b){
			return a & b;
		}

		static inline std::bitset<N> up_inclusion_set_dual_operation(const std::bitset<N>& a, const std::bitset<N>& b){
			return a | b;
		}

		static inline std::binary_function<size_t, size_t, bool> card_map_comparator();

		static inline std::binary_function<size_t, size_t, bool> down_inclusion_card_map_comparator(){
			return std::less<size_t>();
		}

		static inline std::binary_function<size_t, size_t, bool> up_inclusion_card_map_comparator(){
			return std::greater<size_t>();
		}

		static inline std::binary_function<size_t, size_t, bool> card_map_dual_comparator();

		static inline std::binary_function<size_t, size_t, bool> down_inclusion_card_map_dual_comparator(){
			return std::greater<size_t>();
		}

		static inline std::binary_function<size_t, size_t, bool> up_inclusion_card_map_dual_comparator(){
			return std::less<size_t>();
		}

		static inline std::vector<set_N_value<T, N>* > elements_related_to(const powerset_btree<T, N>& tree, const std::bitset<N>& set);

		static inline std::vector<set_N_value<T, N>* > elements_down_related_to(const powerset_btree<T, N>& tree, const std::bitset<N>& set){
			return tree.subsets_of(set);
		}

		static inline std::vector<set_N_value<T, N>* > elements_up_related_to(const powerset_btree<T, N>& tree, const std::bitset<N>& set){
			return tree.supersets_of(set);
		}

		static inline std::vector<set_N_value<T, N>* > elements_dually_related_to(const powerset_btree<T, N>& tree, const std::bitset<N>& set);

		static inline std::vector<set_N_value<T, N>* > elements_dually_down_related_to(const powerset_btree<T, N>& tree, const std::bitset<N>& set){
			return tree.supersets_of(set);
		}

		static inline std::vector<set_N_value<T, N>* > elements_dually_up_related_to(const powerset_btree<T, N>& tree, const std::bitset<N>& set){
			return tree.subsets_of(set);
		}

		static inline void consonant_operation(T& value, const T& preceding_value);

		static inline void zeta_consonant_operation(T& value, T& preceding_value){
			value_inplace_operation(value, preceding_value);
			preceding_value = value;
		}

		static inline void mobius_consonant_operation(T& value, T& preceding_value){
			T old_value = value;
			value_inplace_operation(value, preceding_value);
			preceding_value = old_value;
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
	};
}		// namespace efficient_DST

#endif // EFFICIENT_DST_EMT_HPP
