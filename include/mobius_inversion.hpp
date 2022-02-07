#ifndef EFFICIENT_DST_MOBIUS_INVERSION_TEMPLATE_HPP
#define EFFICIENT_DST_MOBIUS_INVERSION_TEMPLATE_HPP


#include <unordered_set>

#include "macros.hpp"
#include <powerset_btree.hpp>


namespace efficient_DST{

	enum class scheme_type_t: int8_t { direct, consonant, semilattice, lattice };

	template<size_t N, typename T = float>
	struct iota_elements {
		typedef typename sample_space<N>::subset subset;

		static inline std::vector<subset > join_irreducible(
			const powerset_btree<N, T>& support
		){
			// join-irreducible elements of the lattice support
			std::vector<subset > sequence;
			std::unordered_set<subset > manifest;
			std::map<size_t, std::vector<subset >, std::less<size_t> > iota_elements_card_map;
			subset singleton = 1;
			for (size_t i = 0; i < N; ++i) {
				const std::vector<set_N_value<N, T>* >& support_supersets = support.supersets_of(singleton);

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
			}
			sequence.reserve(manifest.size());

			for (const auto& c_iota_elements : iota_elements_card_map) {
				const std::vector<subset >& elements = c_iota_elements.second;
				for (size_t i = 0; i < elements.size(); ++i) {
					sequence.emplace_back(elements[i]);
				}
			}
			return sequence;
		}

		static inline std::vector<subset > meet_irreducible(
			const powerset_btree<N, T>& support
		){
			// meet-irreducible elements of the lattice support
			std::vector<subset > sequence;
			std::unordered_set<subset > manifest;
			std::map<size_t, std::vector<subset >, std::greater<size_t> > iota_elements_card_map;
			subset singleton = 1;
			subset singleton_dual = ~singleton;
			for (size_t i = 0; i < N; ++i) {
				const std::vector<set_N_value<N, T>* >& support_subsets = support.subsets_of(singleton_dual);

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
			}
			sequence.reserve(manifest.size());

			for (const auto& c_iota_elements : iota_elements_card_map) {
				const std::vector<subset >& elements = c_iota_elements.second;
				for (size_t i = 0; i < elements.size(); ++i) {
					sequence.emplace_back(elements[i]);
				}
			}
			return sequence;
		}
	};

	template<size_t N, typename T = float>
	struct up_inclusion {
		typedef typename sample_space<N>::subset subset;

		static inline std::vector<subset > compute_iota_elements(
			const powerset_btree<N, T>& support
		){
			return iota_elements<N, T>::meet_irreducible(support);
		}

		static inline std::vector<subset > compute_iota_elements_dual(
			const powerset_btree<N, T>& support
		){
			return iota_elements<N, T>::join_irreducible(support);
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

		static inline subset set_dual_operation(const subset& a, const subset& b){
			return a | b;
		}

		static inline bool naturally_ranked(const size_t a, const size_t& b){
			return a >= b;
		}

		static inline std::map<size_t, std::vector<set_N_value<N, T>* >, std::greater<size_t> > card_mapping(const powerset_btree<N, T>& tree){
			return tree.elements_by_descending_cardinality();
		}

		static inline std::map<size_t, std::vector<set_N_value<N, T>* >, std::less<size_t> > card_mapping_dual(const powerset_btree<N, T>& tree){
			return tree.elements_by_ascending_cardinality();
		}

		static inline std::vector<set_N_value<N, T>* > elements_related_to(const powerset_btree<N, T>& tree, const subset& set){
			return tree.supersets_of(set);
		}

		static inline std::vector<set_N_value<N, set_N_value<N, T>* >* > addresses_related_to(const powerset_btree<N, set_N_value<N, T>* >& tree, const subset& set){
			return tree.supersets_of(set);
		}

		static inline std::vector<set_N_value<N, T>* > elements_dually_related_to(const powerset_btree<N, T>& tree, const subset& set){
			return tree.subsets_of(set);
		}

		static inline std::vector<set_N_value<N, bool>* > sets_dually_related_to(const powerset_btree<N, bool>& tree, const subset& set){
			return tree.subsets_of(set);
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

		static inline std::vector<subset > compute_iota_elements(
			const powerset_btree<N, T>& support
		){
			return iota_elements<N, T>::join_irreducible(support);
		}

		static inline std::vector<subset > compute_iota_elements_dual(
			const powerset_btree<N, T>& support
		){
			return iota_elements<N, T>::meet_irreducible(support);
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

		static inline subset set_dual_operation(const subset& a, const subset& b){
			return a & b;
		}

		static inline bool naturally_ranked(const size_t a, const size_t& b){
			return a <= b;
		}

		static inline std::map<size_t, std::vector<set_N_value<N, T>* >, std::less<size_t> > card_mapping(const powerset_btree<N, T>& tree){
			return tree.elements_by_ascending_cardinality();
		}

		static inline std::map<size_t, std::vector<set_N_value<N, T>* >, std::greater<size_t> > card_mapping_dual(const powerset_btree<N, T>& tree){
			return tree.elements_by_descending_cardinality();
		}

		static inline std::vector<set_N_value<N, T>* > elements_related_to(const powerset_btree<N, T>& tree, const subset& set){
			return tree.subsets_of(set);
		}

		static inline std::vector<set_N_value<N, set_N_value<N, T>* >* > addresses_related_to(const powerset_btree<N, set_N_value<N, T>* >& tree, const subset& set){
			return tree.subsets_of(set);
		}

		static inline std::vector<set_N_value<N, T>* > elements_dually_related_to(const powerset_btree<N, T>& tree, const subset& set){
			return tree.supersets_of(set);
		}

		static inline std::vector<set_N_value<N, bool>* > sets_dually_related_to(const powerset_btree<N, bool>& tree, const subset& set){
			return tree.supersets_of(set);
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
				powerset_btree<N, T>& focal_points_tree
		){
			//clock_t t;
			//t = clock();
			powerset_btree<N, T> focal_points_N_initial_values(focal_points_tree);
			/*
			std::function<T(const subset&)> elements_related_to;
			if(order_relation == order_relation_t::subset){
				elements_related_to = focal_points_N_initial_values.subsets_of;
			}else{
				elements_related_to = focal_points_N_initial_values.supersets_of;
			}*/
			const std::vector<set_N_value<N, T>* >& focal_points = focal_points_N_initial_values.elements();
			T val;
			for (size_t i = 0; i < focal_points.size(); ++i) {
				val = focal_points[i]->value;
				const std::vector<set_N_value<N, T>* >& elements = inclusion::elements_related_to(
					focal_points_N_initial_values,
					focal_points[i]->set
				);
				DEBUG(std::clog << "Elements related to " << focal_points[i]->set << " : \n";);
				for (size_t ii = 0; ii < elements.size(); ++ii) {
					DEBUG(std::clog << elements[ii]->set << std::endl;);
					if (elements[ii]->set != focal_points[i]->set){
						value_inplace_operation(val, elements[ii]->value);
					}
				}
				focal_points_tree.update_or_insert(focal_points[i]->set, val);
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
			const auto& focal_points_card_map = inclusion::card_mapping(focal_points_tree);
			T val;
			for (const auto& c_focal_points : focal_points_card_map) {
				const std::vector<set_N_value<N, T>* >& focal_points = c_focal_points.second;
				for (size_t i = 0; i < focal_points.size(); ++i) {
					val = focal_points[i]->value;
					const std::vector<set_N_value<N, T>* >& elements = inclusion::elements_related_to(focal_points_tree, focal_points[i]->set);
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
		std::map<size_t, std::vector<set_N_value<N, T>* >, std::less<size_t> > support_card_map = support.elements_by_ascending_cardinality();

		typename sample_space<N>::subset previous_set = 0;
		for (const auto& c_support_elements : support_card_map) {
			const std::vector<set_N_value<N, T>* >& support_elements = c_support_elements.second;
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
			const std::vector<set_N_value<N, T>* >& elements = support.elements();
			for (size_t i = 0; i < elements.size(); ++i){
				focal_points_tree.insert(elements[i]->set, elements[i]->value);
			}

			const subset& absorbing_set = inclusion::absorbing_set_for_operation();
			const subset& dual_absorbing_set = inclusion::absorbing_set_for_dual_operation();

			subset U = elements[0]->set;
			DEBUG(std::clog << "\nU = "<< elements[0]->set;);

			set_N_value<N, T>* newly_inserted_address;

			for (size_t i = 1; i < elements.size(); ++i) {
				const subset& A = elements[i]->set;
				if (A != dual_absorbing_set){
					DEBUG(std::clog << "\nI = set_operation(U, "<< A << ")\n";);
					const subset& I = inclusion::set_operation(U, A);

					newly_inserted_address = focal_points_tree.insert_or_update_if_null(
						I,
						transformation::neutral_value()
					);
					if (newly_inserted_address && newly_inserted_address->cardinality > 1){
						DEBUG(std::clog << newly_inserted_address->cardinality << "=> Linear focal points computation aborted.\n";);
						return false;
					}
					U = inclusion::set_dual_operation(U, A);
				}
			}
			focal_points_tree.insert_or_update_if_null(
				absorbing_set,
				transformation::neutral_value()
			);
			return true;
		}

		static void execute_direct_transformation(
			powerset_btree<N, T>& focal_points_tree
		) {
			transformation::execute_direct_transformation(focal_points_tree);
		}

		static void execute_consonant_transformation(
			powerset_btree<N, T>& support
		) {
			transformation::execute_consonant_transformation(support);
		}

		static powerset_btree<N, T> FMT_reduced_to_core(
				powerset_btree<N, T>& support
		) {
			powerset_btree<N, T> core_reduced_powerset(support.size());
			const std::vector<subset >& iota_sequence = build_core_reduced_powerset(support, core_reduced_powerset);

			std::vector<set_N_value<N, T>* > powerset_elements = core_reduced_powerset.elements();
			for (size_t i = 0; i < iota_sequence.size(); ++i){
				for (size_t e = 0; e < powerset_elements.size(); ++e){
					subset set_B = inclusion::set_operation(powerset_elements[e]->set, inclusion::FMT_target(iota_sequence[i]));

					if (set_B != powerset_elements[e]->set){
						set_N_value<N, T>* B = core_reduced_powerset.find(set_B);
						transformation::value_inplace_operation(B->value, powerset_elements[e]->value);
					}
				}
			}
			return core_reduced_powerset;
		}

		static void execute_FMT(
				std::vector<T>& transform
		) {
			if(transform.size() != pow(2, N)){
				std::cerr << "\nThe size of the given vector is not 2^N, where N is the given number of outcomes.\n";
//					return powerset_values;
			}

//				std::vector<T> transform(powerset_values);
			size_t sub_powerset_size, sub_powerset_dual_size, index;
			for (size_t i = 1; i <= N; ++i){

				sub_powerset_size = pow(2, i);
				for (size_t j = 1; j <= sub_powerset_size; j += 2){

					sub_powerset_dual_size = pow(2, N - i);
					for (size_t k = 0; k <= sub_powerset_dual_size-1; ++k){
						index = (j-inclusion::target_index_offset_FMT()) * sub_powerset_dual_size + k;
						transformation::value_inplace_operation(transform[index], transform[(j-inclusion::source_index_offset_FMT()) * sub_powerset_dual_size + k]);
					}
				}
			}
//				return transform;
		}

		static void build_bridge_map(
				std::unordered_map<subset, set_N_value<N, T>* >& bridge_map,
				const powerset_btree<N, T>& focal_points_tree,
				const std::vector<subset >& iota_sequence
		) {
			const std::vector<set_N_value<N, T>* >& focal_points = focal_points_tree.elements();

			powerset_btree<N, bool> proxies_missing_targets(iota_sequence.size());

			for (size_t i = 0; i < iota_sequence.size(); ++i) {
				for (size_t e = 0; e < focal_points.size(); ++e){
//					const subset& set = focal_points[e]->set | iota_sequence[i];
					const subset& set = inclusion::set_dual_operation(focal_points[e]->set, iota_sequence[i]);
					if (!focal_points_tree.find(set)){
						proxies_missing_targets.insert(set, true);
					}
				}
			}
			if (proxies_missing_targets.size() == 0){
				return;
			}

//			const std::binary_function<size_t, size_t, bool>& comp = std::less<size_t>();
//			const std::binary_function<size_t, size_t, bool>& comp = card_map_dual_comparator();
//			std::map<size_t, std::vector<set_N_value<N, T>* >, comp > focal_points_card_map = focal_points_tree.elements_by_set_cardinality(comp);
			const auto& focal_points_card_map = inclusion::card_mapping_dual(focal_points_tree);
			//clock_t t;
			//t = clock();
			size_t nb_targets_to_find, cc = 0;
			for (const auto& c_focal_points : focal_points_card_map) {
				for (; cc < c_focal_points.first; ++cc){
					nb_targets_to_find += proxies_missing_targets.get_nb_sets_of_cardinality(cc);
				}
				if (nb_targets_to_find > 0){
					const std::vector<set_N_value<N, T>* >& focal_points_with_same_size = c_focal_points.second;

					for (size_t i = 0; i < focal_points_with_same_size.size(); ++i){
//						const std::vector<set_N_value<bool, N>* >& subsets = proxies_missing_targets.subsets_of(focal_points_with_same_size[i]->set);
						const std::vector<set_N_value<N, bool>* >& proxies = inclusion::sets_dually_related_to(proxies_missing_targets, focal_points_with_same_size[i]->set);

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
				powerset_btree<N, T>& focal_points_tree,
				const std::vector<subset >& iota_sequence
		) {
			//clock_t t;

			if (iota_sequence.size() == 0)
				return;

			std::vector<subset > sync_sequence;
			sync_sequence.reserve(iota_sequence.size());

			sync_sequence.emplace_back(iota_sequence[0]);
			for (size_t i = 1; i < iota_sequence.size(); ++i){
//				sync_sequence.emplace_back(sync_sequence[i-1] | iota_sequence[i]);
				sync_sequence.emplace_back(inclusion::set_dual_operation(sync_sequence[i-1], iota_sequence[i]));
			}
			std::unordered_map<subset, set_N_value<N, T>* > bridge_map;
			build_bridge_map(bridge_map, focal_points_tree, iota_sequence);

			//t = clock();
			size_t nb_iota = iota_sequence.size()-1;
//			size_t iota_index;
			const std::vector<set_N_value<N, T>* >& focal_points = focal_points_tree.elements();
			for (size_t i = 0; i < iota_sequence.size(); ++i){
//				if (transform_type == transform_type_t::zeta)
//					iota_index = i;
//				else
//					iota_index = nb_iota - i;
				const size_t& iota_index = transformation::subgraph_dual_index(nb_iota, i);

				for (size_t e = 0; e < focal_points.size(); ++e){
//					const subset& proxy = focal_points[e]->set | iota_sequence[iota_index];
					const subset& proxy = inclusion::set_dual_operation(focal_points[e]->set, iota_sequence[iota_index]);

					if (focal_points[e]->set != proxy){
						const set_N_value<N, T>* coupled_set = bridge_map[proxy];
//						if (coupled_set && FOD<N>::is_subset_of(coupled_set->set, focal_points[e]->set | sync_sequence[iota_index])){
						const subset& sync_set = inclusion::set_dual_operation(focal_points[e]->set, sync_sequence[iota_index]);
						if (coupled_set && (coupled_set->set & sync_set) == coupled_set->set){
							transformation::value_inplace_operation(focal_points[e]->value, coupled_set->value);
						}
					}
				}
			}
			//t = clock() - t;
			//std::cout << (((float) t)/CLOCKS_PER_SEC) << std::endl;
		}


		static void execute_EMT_with_lattice(
				powerset_btree<N, T>& lattice_support,
				const std::vector<subset >& iota_sequence
		) {
			const std::vector<set_N_value<N, T>* >& lattice_support_elements = lattice_support.elements();

			size_t nb_iota = iota_sequence.size()-1;
//			size_t iota_index;
			for (size_t i = 0; i < iota_sequence.size(); ++i){
//				if (transform_type == transform_type_t::Mobius)
//					iota_index = i;
//				else
//					iota_index = nb_iota - i;
				const size_t& iota_index = transformation::subgraph_index(nb_iota, i);

				for (size_t e = 0; e < lattice_support_elements.size(); ++e){
//					const subset& set_B = lattice_support_elements[e]->set | iota_sequence[iota_index];
					const subset& set_B = inclusion::set_operation(lattice_support_elements[e]->set, iota_sequence[iota_index]);

					if (set_B != lattice_support_elements[e]->set){
						set_N_value<N, T>* B = lattice_support.find(set_B);

						if (B){
							transformation::value_inplace_operation(B->value, lattice_support_elements[e]->value);
						}
					}
				}
			}
		}

		static void build_semilattice_support(
				const powerset_btree<N, T>& support,
				powerset_btree<N, T>& focal_points_tree
		) {
			std::vector<subset > focal_points;
			const std::vector<set_N_value<N, T>* >& support_elements = support.elements();
			focal_points.reserve(support.size());

			for (size_t i = 0; i < support_elements.size(); ++i){
				focal_points.emplace_back(support_elements[i]->set);
				focal_points_tree.insert(support_elements[i]->set, support_elements[i]->value);
			}

			set_N_value<N, T>* newly_inserted_address;

			for (size_t i = 0; i < support_elements.size(); ++i){
				size_t end = focal_points.size();
				for (size_t ii = i+1; ii < end; ++ii){
					const subset& focal_point = inclusion::set_operation(support_elements[i]->set, focal_points[ii]);
					newly_inserted_address = focal_points_tree.insert_or_update_if_null(
						focal_point,
						transformation::neutral_value()
					);
					if (newly_inserted_address){
						focal_points.emplace_back(focal_point);
					}
				}
			}
		}

		/*
		 * In this function, this->iota_sequence is supposed to contain all regular iota elements (i.e. join-irreducible of this lattice).
		 */
		static const std::vector<subset >& build_truncated_lattice_support(
				const powerset_btree<N, T>& support,
				powerset_btree<N, T>& truncated_lattice_support
		) {
			const std::vector<subset >& iota_sequence = inclusion::compute_iota_elements(
				support
			);

			std::vector<subset > focal_points;
			focal_points.reserve(support.size());

			const std::vector<set_N_value<N, T>* >& elements = support.elements();
			for (size_t i = 0; i < elements.size(); ++i){
				focal_points.emplace_back(elements[i]->set);
				truncated_lattice_support.insert(elements[i]->set, elements[i]->value);
			}

			set_N_value<N, T>* newly_inserted_address;

			for (size_t i = 0; i < iota_sequence.size(); ++i) {
				size_t end = focal_points.size();
				for (size_t ii = 0; ii < end; ++ii) {
					const subset& new_set = inclusion::set_operation(iota_sequence[i], focal_points[ii]);
					newly_inserted_address = truncated_lattice_support.insert_or_update_if_null(
						new_set,
						transformation::neutral_value()
					);
					if (newly_inserted_address){
						focal_points.emplace_back(new_set);
					}
				}
			}
			DEBUG({
				std::clog << "\nCropped lattice support: \n";
				truncated_lattice_support.print(std::clog);
			});
			return iota_sequence;
		}

		static const std::vector<subset >& build_core_reduced_powerset(
				const powerset_btree<N, T>& support,
				powerset_btree<N, T>& core_reduced_powerset
		) {
			const std::vector<set_N_value<N, T>* >& support_elements = support.elements();
			subset core = 0;
			for(size_t i = 0; i < support_elements.size(); ++i){
				core |= support_elements[i]->set;
			}
			std::vector<subset > iota_sequence;
			subset singleton = 1;
			for(size_t i = 0; i < N; ++i){
				if ((singleton & core) != 0){
					iota_sequence.emplace_back(singleton);
				}
				singleton <<= 1;
			}

			std::vector<subset > focal_points;
			focal_points.reserve(pow(2, iota_sequence.size()));

			singleton = 0;
			set_N_value<N, T>* inserted_node = support.find(singleton);
			if (!inserted_node || inserted_node->is_null){
				core_reduced_powerset.insert(singleton, transformation::neutral_value());
			}else{
				core_reduced_powerset.insert(singleton, inserted_node->value);
			}
			focal_points.emplace_back(singleton);

			for (size_t i = 0; i < iota_sequence.size(); ++i) {
				size_t end = focal_points.size();
				for (size_t ii = 0; ii < end; ++ii) {
					const subset& new_set = iota_sequence[i] | focal_points[ii];
					inserted_node = support.find(new_set);
					if (!inserted_node || inserted_node->is_null){
						inserted_node = core_reduced_powerset.insert_or_update_if_null(new_set, transformation::neutral_value());
					}else{
						inserted_node = core_reduced_powerset.insert_or_update_if_null(new_set, inserted_node->value);
					}
					if (inserted_node){
						focal_points.emplace_back(new_set);
					}
				}
			}
			DEBUG({
				std::clog << "\nCropped lattice support: \n";
				core_reduced_powerset.print(std::clog);
			});
			return iota_sequence;
		}

		static scheme_type_t autoset_and_build(
				const powerset_btree<N, T>& support,
				powerset_btree<N, T>& focal_points_tree,
				std::vector<subset >& iota_sequence
		){
			// check if original_structure is almost Bayesian
			DEBUG(std::clog << "Linear analysis:" << std::endl;);
			const bool& is_almost_bayesian = try_linear_focal_points_computation (
				support,
				focal_points_tree
			);

			if(is_almost_bayesian){
				DEBUG(std::clog << "almost Bayesian." << std::endl;);
				return scheme_type_t::direct;
			}else{
				DEBUG(std::clog << "Consonance check:" << std::endl;);
				const bool& is_consonant = consonance_check<N, T>(support);

				if(is_consonant){
					DEBUG(std::clog << "consonant." << std::endl;);
					focal_points_tree.copy(support);
					return scheme_type_t::consonant;
				}else{
					DEBUG(std::clog << "not consonant." << std::endl;);

					if(support.size() < 3 * N){
						DEBUG({
							std::clog << "Number of focal sets equivalent to number of outcomes." << std::endl;
							std::clog << "Transform to semilattice:" << std::endl;
						});

						build_semilattice_support(
							support,
							focal_points_tree
						);

						if(focal_points_tree.size() < 3 * N){
							DEBUG(std::clog << "Number of focal points also equivalent to number of outcomes." << std::endl;);
							return scheme_type_t::direct;
						}else{
							DEBUG(std::clog << "Number of focal points superior to number of outcomes." << std::endl;);

							iota_sequence = inclusion::compute_iota_elements_dual(
								support
							);
							return scheme_type_t::semilattice;
						}
					}else{
						DEBUG({
							std::clog << "Number of focal sets superior to number of outcomes." << std::endl;
							std::clog << "Transform to lattice:" << std::endl;
						});

						iota_sequence = build_truncated_lattice_support(
							support,
							focal_points_tree
						);
						return scheme_type_t::lattice;
					}
				}
			}
		}

		static void execute(
			powerset_btree<N, T>& focal_points_tree,
			const std::vector<subset >& iota_sequence,
			const scheme_type_t& scheme_type
		) {

			switch (scheme_type) {
				case scheme_type_t::direct:
					execute_direct_transformation(
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
						iota_sequence
					);
					break;

				case scheme_type_t::lattice:
					execute_EMT_with_lattice(
						focal_points_tree,
						iota_sequence
					);

					break;

				default:
					break;
			}
		}

		static powerset_btree<N, T> direct_transformation(
				const powerset_btree<N, T>& support
		) {
			powerset_btree<N, T> focal_points_tree(support.size());
			scheme_type_t scheme_type = scheme_type_t::direct;
			// check if original_structure is almost Bayesian
			const bool& is_almost_bayesian = try_linear_focal_points_computation (
				support,
				focal_points_tree
			);

			if (!is_almost_bayesian){
				DEBUG(std::clog << "Consonance check:" << std::endl;);
				const bool& is_consonant = consonance_check<N, T>(support);

				if(is_consonant){
					focal_points_tree.copy(support);
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
					focal_points_tree
				);
			}else{
				execute_consonant_transformation(
					focal_points_tree
				);
			}
			return focal_points_tree;
		}

		static powerset_btree<N, T> EMT_with_semilattice(
			const powerset_btree<N, T>& support
		) {
			powerset_btree<N, T> focal_points_tree(support.size());
			build_semilattice_support(
				support,
				focal_points_tree
			);
			const std::vector<subset >& iota_sequence = inclusion::compute_iota_elements_dual(
				support
			);
			execute_EMT_with_semilattice(
					focal_points_tree,
					iota_sequence
			);
			return focal_points_tree;
		}

		static powerset_btree<N, T> EMT_with_lattice(
			const powerset_btree<N, T>& support
		) {
			powerset_btree<N, T> truncated_lattice_support(support.size());
			const std::vector<subset >& iota_sequence = build_truncated_lattice_support(
				support,
				truncated_lattice_support
			);
			execute_EMT_with_lattice(
					truncated_lattice_support,
					iota_sequence
			);
			return truncated_lattice_support;
		}
	};
}		// namespace efficient_DST

#endif // EFFICIENT_DST_MOBIUS_INVERSION_TEMPLATE_HPP
