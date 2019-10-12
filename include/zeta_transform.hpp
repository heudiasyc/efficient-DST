#ifndef EFFICIENT_DST_ZETA_TRANSFORM_HPP
#define EFFICIENT_DST_ZETA_TRANSFORM_HPP

#include <mobius_transform.hpp>
#include <powerset_btree.hpp>
#include <math.h>


namespace efficient_DST{

	enum class order_relation_t: bool { subset, superset };
	enum class operation_t: bool { addition, multiplication };

	template <class T=double>
	class zeta_transform : public mobius_transform<T> {
	public:
		enum class scheme_type_t: int8_t { direct, consonant, semilattice, lattice };
		enum class version_t: bool { regular, dual };
		enum class transform_type_t: bool { zeta, Mobius };

	protected:
		scheme_type_t scheme_type;
		std::vector<size_t> iota_fiber_sequence;	// fod element indices
		std::vector<boost::dynamic_bitset<> > iota_sequence;
		powerset_btree<set_N_value<T>* > dual_definition;
		powerset_btree<T> original_mobius_transform;
		operation_t original_operation;
		std::unordered_map<size_t, powerset_btree<set_N_value<T>* > > definition_by_cardinality;
		std::vector<size_t> ordered_cardinalities_in_definition;

	public:
		const order_relation_t order_relation;

		zeta_transform(const zeta_transform<T>& z) :
			mobius_transform<T>(z.definition),
			scheme_type(z.scheme_type),
			dual_definition(z.dual_definition),
			original_mobius_transform(z.original_mobius_transform),
			original_operation(z.original_operation),
			definition_by_cardinality(z.definition_by_cardinality),
			order_relation(z.order_relation)
		{}


		/*
		 * Constructor when you have a zeta transform such as the commonality or implicability ones without knowledge about its focal points.
		 * - powerset_values contains all images for the whole powerset of FOD.
		 * - fod is the corresponding FOD.
		 * - order_relation is the order relation of this zeta transform (e.g. commonality->superset or implicability->subset).
		 */
		zeta_transform(
				const std::vector<T>& powerset_values,
				FOD& fod,
				const order_relation_t& order_relation) :
			mobius_transform<T>(fod),
			scheme_type(scheme_type_t::semilattice),
			dual_definition(&fod, mobius_transform<T>::block_size),
			original_mobius_transform(&fod, mobius_transform<T>::block_size),
			original_operation(operation_t::addition),
			order_relation(order_relation)
		{
			if(powerset_values.size() != pow(2, fod.size())){
				std::cerr << "\nThe given vector does not feature the same size as the powerset of the given FOD.\n";
				return;
			}
			to_semilattice_support(
				powerset_values
			);
			set_semilattice_computation_scheme(this->definition);
			set_definition_by_cardinality();
		}


		/*
		 * Constructor when you have a zeta transform such as the commonality or implicability ones AND know its focal points.
		 * - focal_points_tree is supposed to contain all focal points and their image.
		 * - order_relation is the order relation of this zeta transform (e.g. commonality->superset or implicability->subset).
		 */
		zeta_transform(
				const powerset_btree<T>& focal_points_tree,
				const order_relation_t& order_relation) :
			mobius_transform<T>(focal_points_tree),
			scheme_type(scheme_type_t::semilattice),
			dual_definition(focal_points_tree.get_FOD(), focal_points_tree.get_block_size()),
			original_mobius_transform(focal_points_tree.get_FOD(), focal_points_tree.get_block_size()),
			original_operation(operation_t::addition),
			order_relation(order_relation)
		{
			if (this->order_relation == order_relation_t::subset)
				powerset_btree<T>::flip_powerset_from_to(this->definition, this->dual_definition);
			set_semilattice_computation_scheme(this->definition);
			set_definition_by_cardinality();
		}


		/*
		 * Constructor when you have a Möbius transform such as the mass or conjunctive/disjunctive weight function.
		 * - support is supposed to contain all focal sets and their image.
		 * - order_relation is the order relation of this zeta transform (e.g. commonality->superset or implicability->subset).
		 * - transform_operation is the operation of the zeta transform (e.g. in DST, we usually use the addition on the mass function
		 * and the multiplication on the conjunctive/disjunctive weight function).
		 */
		zeta_transform(
				const powerset_btree<T>& support,
				const order_relation_t& order_relation,
				const operation_t& transform_operation) :
			mobius_transform<T>(support),
			scheme_type(scheme_type_t::direct),
			dual_definition(support.get_FOD(), support.get_block_size()),
			original_mobius_transform(support.get_FOD(), support.get_block_size()),
			original_operation(transform_operation),
			order_relation(order_relation)
		{
			powerset_btree<T>::flip_powerset_from_to(this->definition, this->dual_definition);

			T neutral_value;
			if (transform_operation == operation_t::addition){
				neutral_value = 0;
			} else{
				neutral_value = 1;
			}
			set_computation_scheme(support, neutral_value);
			execute_scheme(
					transform_operation,
					transform_type_t::zeta,
					this->definition,
					this->dual_definition);
			set_definition_by_cardinality();
		}


		powerset_btree<T> inversion(const operation_t& transform_operation) const {
			if (transform_operation == this->original_operation && this->original_mobius_transform.size() > 0){
				return this->original_mobius_transform;
			}
			powerset_btree<T> mobius_transform_definition(this->definition);
			powerset_btree<set_N_value<T>* > mobius_transform_dual_definition(this->definition.get_FOD(), this->definition.get_block_size());
			powerset_btree<T>::flip_powerset_from_to(mobius_transform_definition, mobius_transform_dual_definition);
			execute_scheme(
					transform_operation,
					transform_type_t::Mobius,
					mobius_transform_definition,
					mobius_transform_dual_definition);
			return mobius_transform_definition;
		}

		const powerset_btree<T>& get_original_mobius_transform() const {
			return this->original_mobius_transform;
		}

		T at_emptyset() const {
			set_N_value<T>* set_value = this->definition.sub_fod_of_size(0);
			if(set_value){
				return set_value->value;
			}
			boost::dynamic_bitset<> emptyset(this->definition.get_FOD_size());
			return find_non_focal_point_image(emptyset);
		}

		T at_fod() const {
			set_N_value<T>* set_value = this->definition.sub_fod_of_size(this->definition.get_FOD_size());
			if(set_value){
				return set_value->value;
			}
			boost::dynamic_bitset<> fod(this->definition.get_FOD_size());
			fod.set();
			return find_non_focal_point_image(fod);
		}

		T operator[](const std::vector<std::string>& labels) const {
			return find(this->definition.get_FOD()->to_set(labels));
		}

		T find(const boost::dynamic_bitset<>& set) const {
			set_N_value<T>* set_value = this->definition[set];
			if(set_value){
				return set_value->value;
			}
			return find_non_focal_point_image(set);
		}

	protected:

		T find_non_focal_point_image(const boost::dynamic_bitset<>& set) const {
			set_N_value<set_N_value<T>* >* A = nullptr;
			size_t card = set.count();

			if(this->order_relation == order_relation_t::subset){
				if(card >= this->ordered_cardinalities_in_definition.back()){
					for(size_t c = 0; c < this->ordered_cardinalities_in_definition.size(); ++c){
						if(this->ordered_cardinalities_in_definition[c] <= card){
							card = c;
							break;
						}
					}
					for(size_t c = card; c < this->ordered_cardinalities_in_definition.size(); ++c){
						const powerset_btree<set_N_value<T>* >& p = this->definition_by_cardinality.at(this->ordered_cardinalities_in_definition[c]);
						const std::vector<set_N_value<set_N_value<T>* >* >& subsets = p.subsets_of(set);

						if(subsets.size() > 0){
							A = subsets[0];
							break;
						}
					}
				}
			}else{
				if(card <= this->ordered_cardinalities_in_definition.back()){
					for(size_t c = 0; c < this->ordered_cardinalities_in_definition.size(); ++c){
						if(this->ordered_cardinalities_in_definition[c] >= card){
							card = c;
							break;
						}
					}
					for(size_t c = card; c < this->ordered_cardinalities_in_definition.size(); ++c){
						const powerset_btree<set_N_value<T>* >& p = this->definition_by_cardinality.at(this->ordered_cardinalities_in_definition[c]);
						std::vector<set_N_value<set_N_value<T>* >* > supersets = p.supersets_of(set);

						if(supersets.size() > 0){
							A = supersets[0];
							break;
						}
					}
				}
			}
			if (A){
				return A->value->value;
			}else{
				if(this->original_operation == operation_t::addition){
					return 0;
				}else{
					return 1;
				}
			}
		}


		void set_definition_by_cardinality(){
			std::unordered_map<size_t, std::vector<set_N_value<T>* > > definition_card_map = this->definition.elements_by_set_cardinality();
			this->ordered_cardinalities_in_definition = powerset_btree<T>::get_sorted_cardinalities(definition_card_map, *this->definition.get_FOD());
			if (this->order_relation == order_relation_t::subset) {
				std::reverse(this->ordered_cardinalities_in_definition.begin(), this->ordered_cardinalities_in_definition.end());
			}

			this->definition_by_cardinality.reserve(this->ordered_cardinalities_in_definition.size());

			for(size_t c = 0; c < this->ordered_cardinalities_in_definition.size(); ++c){
				this->definition_by_cardinality.emplace(std::piecewise_construct, std::make_tuple(this->ordered_cardinalities_in_definition[c]),
						std::make_tuple(this->definition.get_FOD(), this->definition.get_block_size()));
				powerset_btree<set_N_value<T>* >& p_c = this->definition_by_cardinality[this->ordered_cardinalities_in_definition[c]];
				const std::vector<set_N_value<T>* >& elements = definition_card_map[this->ordered_cardinalities_in_definition[c]];
				for(size_t i = 0; i < elements.size(); ++i){
					p_c.insert(elements[i]->set, elements[i]);
				}
			}
		}


		void set_computation_scheme(
				const powerset_btree<T>& support,
				const T& neutral_value
				){
			// Cardinality map of focal sets in original_structure
			std::unordered_map<size_t, std::vector<set_N_value<T>* > > support_card_map = support.elements_by_set_cardinality();
			std::vector<set_N_value<T>* > elements_generated_from_support;

			// check if original_structure is almost Bayesian
			std::clog << "Linear analysis:" << std::endl;
			const bool& is_almost_bayesian = linear_analysis_of_support(
				support_card_map,
				elements_generated_from_support,
				neutral_value
			);

			if(is_almost_bayesian){
				std::clog << "almost Bayesian." << std::endl;
				this->scheme_type = scheme_type_t::direct;
			}else{
				// sort sets in support by increasing cardinality
				const std::vector<size_t>& support_ordered_cardinalities = powerset_btree<T>::get_sorted_cardinalities(support_card_map, *support.get_FOD());

				std::clog << "Consonance check:" << std::endl;
				const bool& is_consonant = consonance_check(
					elements_generated_from_support,
					support_card_map,
					support_ordered_cardinalities
				);

				if(is_consonant){
					std::clog << "consonant." << std::endl;
					this->scheme_type = scheme_type_t::consonant;
				}else{
					std::clog << "not consonant." << std::endl;

					if(support.size() < 3 * support.get_FOD_size()){
						std::clog << "Number of focal sets equivalent to |FOD|." << std::endl;
						std::clog << "Transform to semilattice:" << std::endl;

						to_semilattice_support(
							support_card_map,
							support_ordered_cardinalities,
							elements_generated_from_support,
							neutral_value
						);

						if(this->definition.size() < 3 * support.get_FOD_size()){
							std::clog << "Number of focal points also equivalent to |FOD|." << std::endl;
							this->scheme_type = scheme_type_t::direct;
						}else{
							std::clog << "Number of focal points superior to |FOD|." << std::endl;
							set_semilattice_computation_scheme(support);
							this->scheme_type = scheme_type_t::semilattice;
						}
					}else{
						std::clog << "Number of focal sets superior to |FOD|." << std::endl;
						std::clog << "Transform to lattice:" << std::endl;

						compute_iota_elements(version_t::regular, support);
						to_lattice_support(neutral_value);
						this->scheme_type = scheme_type_t::lattice;
					}
				}
			}
		}


		void set_semilattice_computation_scheme(const powerset_btree<T>& support){
			if(this->order_relation == order_relation_t::subset){
				compute_iota_elements(version_t::dual, support);
				for (size_t o = 0; o < this->iota_sequence.size(); ++o){
					this->iota_sequence[o].flip();
				}
			}else{
				compute_iota_elements(version_t::regular, support);
			}
		}


		void execute_scheme(
				const operation_t& transform_operation,
				const transform_type_t& transform_type,
				powerset_btree<T>& focal_points_tree,
				powerset_btree<set_N_value<T>* >& focal_points_dual_tree
				) const {
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

			switch (this->scheme_type) {
				case scheme_type_t::direct:
					if (transform_type == transform_type_t::zeta){
						directly_compute_transform(
								focal_points_tree,
								range_binary_operator);
					}else{
						directly_compute_inversion(
								focal_points_tree,
								range_binary_operator);
					}
					break;

				case scheme_type_t::consonant:
					consonant_computations(
							focal_points_tree,
							transform_type,
							range_binary_operator,
							neutral_value);
					break;

				case scheme_type_t::semilattice:
					EMT_with_semilattice(
							focal_points_tree,
							focal_points_dual_tree,
							transform_type,
							range_binary_operator);
					break;

				case scheme_type_t::lattice:
					EMT_with_lattice(
							focal_points_tree,
							transform_type,
							range_binary_operator);
					break;

				default:
					break;
			}
		}


		static const T addition(const T& a, const T& b){
			return a+b;
		}

		static const T subtraction(const T& a, const T& b){
			return a-b;
		}

		static const T multiplication(const T& a, const T& b){
			return a*b;
		}

		static const T division(const T& a, const T& b){
			return a/b;
		}


		void directly_compute_transform(
				powerset_btree<T>& focal_points_tree,
				std::function<T(const T&, const T&)> range_binary_operator
				) const {
			powerset_btree<T> focal_points_N_initial_values(focal_points_tree);
			T val;
			const std::vector<set_N_value<T>* >& focal_points = focal_points_N_initial_values.elements();
			std::vector<set_N_value<T>* > elements;

			for (size_t i = 0; i < focal_points.size(); ++i) {
				val = focal_points[i]->value;
				if (this->order_relation == order_relation_t::subset){
					elements = focal_points_N_initial_values.strict_subsets_of(focal_points[i]->set);
				} else{
					elements = focal_points_N_initial_values.strict_supersets_of(focal_points[i]->set);
				}
				for (size_t ii = 0; ii < elements.size(); ++ii) {
					val = range_binary_operator(val, elements[ii]->value);
				}
				focal_points_tree.insert(focal_points[i]->set, val);
			}
		}


		void directly_compute_inversion(
				powerset_btree<T>& focal_points_tree,
				std::function<T(const T&, const T&)> range_binary_operator
				) const {
			T val;
			std::unordered_map<size_t, std::vector<set_N_value<T>* > > focal_points_card_map = focal_points_tree.elements_by_set_cardinality();
			std::vector<size_t> focal_points_ordered_cardinalities = powerset_btree<T>::get_sorted_cardinalities(focal_points_card_map, *focal_points_tree.get_FOD());
			std::vector<set_N_value<T>* > elements;

			if (this->order_relation == order_relation_t::superset) {
				std::reverse(focal_points_ordered_cardinalities.begin(), focal_points_ordered_cardinalities.end());
			}
			for (size_t c = 0; c < focal_points_ordered_cardinalities.size(); ++c) {
				const std::vector<set_N_value<T>* >& focal_points = focal_points_card_map[focal_points_ordered_cardinalities[c]];
				for (size_t i = 0; i < focal_points.size(); ++i) {
					val = focal_points[i]->value;
					if (this->order_relation == order_relation_t::subset){
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


		void consonant_computations(
				powerset_btree<T>& support,
				const transform_type_t& transform_type,
				std::function<T(const T&, const T&)> range_binary_operator,
				const T& neutral_value
				) const {
			std::unordered_map<size_t, std::vector<set_N_value<T>* > > support_card_map = support.elements_by_set_cardinality();
			std::vector<size_t> support_ordered_cardinalities = powerset_btree<T>::get_sorted_cardinalities(support_card_map, *support.get_FOD());
			T value, preceding_value = neutral_value;

			if (this->order_relation == order_relation_t::superset) {
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


		void EMT_with_lower_semilattice(
				powerset_btree<T>& focal_points_tree,
				const transform_type_t& transform_type,
				std::function<T(const T&, const T&)> range_binary_operator,
				std::vector<boost::dynamic_bitset<> > sync_sequence
				) const {

			std::vector<set_N_value<T>* > focal_points = focal_points_tree.elements();
			const powerset_btree<set_N_value<T>* >& focal_points_superset_map = powerset_btree<T>::EMT_superset_map(focal_points_tree);

			size_t iota_index;
			for (size_t o = 0; o < this->iota_sequence.size(); ++o){
				if (transform_type == transform_type_t::zeta)
					iota_index = o;
				else
					iota_index = this->iota_sequence.size()-1 - o;

				for (size_t i = 0; i < focal_points.size(); ++i){
					if (!focal_points[i]->set[this->iota_fiber_sequence[iota_index]]){
						boost::dynamic_bitset<> set_B = (const boost::dynamic_bitset<>&) focal_points[i]->set;
						set_B.set(this->iota_fiber_sequence[iota_index]);
						const set_N_value<set_N_value<T>* >* X = focal_points_superset_map.closest_node_containing(set_B);

						if (X && FOD::is_or_is_subset_of(X->set, FOD::set_union(focal_points[i]->set, sync_sequence[iota_index]))){
							focal_points[i]->value = range_binary_operator(focal_points[i]->value, X->value->value);
						}
					}
				}
			}
		}


		void EMT_with_upper_semilattice(
				powerset_btree<set_N_value<T>* >& focal_points_dual_tree,
				const transform_type_t& transform_type,
				std::function<T(const T&, const T&)> range_binary_operator,
				std::vector<boost::dynamic_bitset<> > sync_sequence
				) const {

			std::vector<set_N_value<set_N_value<T>* >* > focal_points_dual = focal_points_dual_tree.elements();
			const powerset_btree<set_N_value<T>* >& focal_points_dual_superset_map = powerset_btree<T>::EMT_superset_map(focal_points_dual_tree);

			size_t iota_index;
			for (size_t o = 0; o < this->iota_sequence.size(); ++o){
				if (transform_type == transform_type_t::zeta)
					iota_index = o;
				else
					iota_index = this->iota_sequence.size()-1 - o;

				for (size_t i = 0; i < focal_points_dual.size(); ++i){
					if (!focal_points_dual[i]->set[this->iota_fiber_sequence[iota_index]]){
						boost::dynamic_bitset<> set_B = (const boost::dynamic_bitset<>&) focal_points_dual[i]->set;
						set_B.set(this->iota_fiber_sequence[iota_index]);
						const set_N_value<set_N_value<T>* >* X = focal_points_dual_superset_map.closest_node_containing(set_B);

						if (X && FOD::is_or_is_subset_of(X->set, FOD::set_union(focal_points_dual[i]->set, sync_sequence[iota_index]))){
							focal_points_dual[i]->value->value = range_binary_operator(focal_points_dual[i]->value->value, X->value->value);
						}
					}
				}
			}
		}


		void EMT_with_semilattice(
				powerset_btree<T>& focal_points_tree,
				powerset_btree<set_N_value<T>* >& focal_points_dual_tree,
				const transform_type_t& transform_type,
				std::function<T(const T&, const T&)> range_binary_operator
				) const {

			if (this->iota_sequence.size() == 0)
				return;

			std::vector<boost::dynamic_bitset<> > sync_sequence;
			sync_sequence.reserve(this->iota_sequence.size());

			sync_sequence.emplace_back(this->iota_sequence[0]);
			for (size_t o = 1; o < this->iota_sequence.size(); ++o){
				sync_sequence.emplace_back(FOD::set_union(sync_sequence[o-1], this->iota_sequence[o]));
			}

			if (this->order_relation == order_relation_t::subset){
				EMT_with_upper_semilattice(
						focal_points_dual_tree,
						transform_type,
						range_binary_operator,
						sync_sequence);
			}else{
				EMT_with_lower_semilattice(
						focal_points_tree,
						transform_type,
						range_binary_operator,
						sync_sequence);
			}
		}


		void EMT_with_lattice(
				powerset_btree<T>& lattice_support,
				const transform_type_t& transform_type,
				std::function<T(const T&, const T&)> range_binary_operator
				) const {
			const std::vector<set_N_value<T>* >& lattice_support_elements = lattice_support.elements();

			bool reverse_sequence = false;
			if(this->order_relation == order_relation_t::subset){
				if(transform_type == transform_type_t::zeta){
					reverse_sequence = true;
				}
			}else{
				if(transform_type == transform_type_t::Mobius){
					reverse_sequence = true;
				}
			}

			size_t iota_index;
			for (size_t o = 0; o < this->iota_sequence.size(); ++o){
				if (reverse_sequence)
					iota_index = this->iota_sequence.size()-1 - o;
				else
					iota_index = o;

				for (size_t i = 0; i < lattice_support_elements.size(); ++i){
					const boost::dynamic_bitset<>& set_B = FOD::set_union(lattice_support_elements[i]->set, this->iota_sequence[iota_index]);

					if (set_B != lattice_support_elements[i]->set){
						set_N_value<T>* B = lattice_support[set_B];

						if (B){
							if (this->order_relation == order_relation_t::subset){
								B->value = range_binary_operator(B->value, lattice_support_elements[i]->value);
							}else{
								lattice_support_elements[i]->value = range_binary_operator(lattice_support_elements[i]->value, B->value);
							}
						}
					}
				}
			}
		}


		set_N_value<T>* insert_focal_point(const boost::dynamic_bitset<>& focal_point, const T& neutral_value) {
			set_N_value<T>* inserted_focal_point = this->definition.insert(focal_point, neutral_value);
			this->dual_definition.insert(~focal_point, inserted_focal_point);
			std::clog << "\nNEW INSERTED FOCAL POINT :\n";
			std::clog << inserted_focal_point->set << std::endl;
			return inserted_focal_point;
		}

		set_N_value<T>* insert_dual_focal_point(const boost::dynamic_bitset<>& dual_focal_point, const T& neutral_value) {
			set_N_value<T>* inserted_focal_point = this->definition.insert(~dual_focal_point, neutral_value);
			this->dual_definition.insert(dual_focal_point, inserted_focal_point);
			std::clog << "\nNEW INSERTED FOCAL POINT :\n";
			std::clog << inserted_focal_point->set << std::endl;
			return inserted_focal_point;
		}


		bool linear_analysis_of_support(
				std::unordered_map<size_t, std::vector<set_N_value<T>* > >& support_card_map,
				std::vector<set_N_value<T>* >& elements_generated_from_support,
				const T& neutral_value
				){

			if(this->order_relation == order_relation_t::subset){
				std::vector<set_N_value<T>* > support_except_emptyset;
				support_except_emptyset.reserve(this->definition.size()-1);

				for(auto kv : support_card_map) {
					if(kv.first > 0){
						for (size_t i = 0; i < kv.second.size(); ++i) {
							support_except_emptyset.emplace_back(kv.second[i]);
						}
					}
				}

				boost::dynamic_bitset<> fod(this->definition.get_FOD_size());
				fod.set();

				boost::dynamic_bitset<> neg_U = support_except_emptyset[0]->set;
				std::clog << "\nneg_U = "<< support_except_emptyset[0]->set;

				for (size_t i = 1; i < support_except_emptyset.size(); ++i) {
					const boost::dynamic_bitset<>& A = support_except_emptyset[i]->set;

					std::clog << "\nneg_I = union(U, "<< A << ")\n";
					const boost::dynamic_bitset<>& neg_I = FOD::set_union(neg_U, A);

					const size_t& I_card = fod.size() - neg_I.count();
					std::clog << "|I| = " << I_card << std::endl;

					if(I_card > 1){
						std::clog << "=> Linear analysis aborted.\n";
						return false;
					}else if(I_card == 1){
						// add it to this->structure if it wasn't already there
						if(!this->definition[neg_I]){
							set_N_value<T>* inserted_focal_point = insert_focal_point(neg_I, neutral_value);
							elements_generated_from_support.push_back(inserted_focal_point);
						}
					}
					neg_U = FOD::set_intersection(neg_U, A);
				}
				// add also fod to this->structure if it wasn't already there
				if(!this->definition[fod]){
					set_N_value<T>* inserted_focal_point = insert_focal_point(fod, neutral_value);
					elements_generated_from_support.push_back(inserted_focal_point);
				}
			}else{
				std::vector<set_N_value<T>* > support_except_FOD;
				support_except_FOD.reserve(this->definition.size()-1);

				for(auto kv : support_card_map) {
					if(kv.first < this->definition.get_FOD_size()){
						for (size_t i = 0; i < kv.second.size(); ++i) {
							support_except_FOD.emplace_back(kv.second[i]);
						}
					}
				}

				const boost::dynamic_bitset<> emptyset(this->definition.get_FOD_size());

				boost::dynamic_bitset<> U = support_except_FOD[0]->set;
				std::clog << "\nU = "<< support_except_FOD[0]->set;

				for (size_t i = 1; i < support_except_FOD.size(); ++i) {
					const boost::dynamic_bitset<>& A = support_except_FOD[i]->set;

					std::clog << "\nI = intersection(U, "<< A << ")\n";
					const boost::dynamic_bitset<>& I = FOD::set_intersection(U, A);

					const size_t& I_card = I.count();
					std::clog << "|I| = " << I_card << std::endl;

					if(I_card > 1){
						std::clog << "=> Linear analysis aborted.\n";
						return false;
					}else if(I_card == 1){
						// add it to this->structure if it wasn't already there
						if(!this->definition[I]){
							set_N_value<T>* inserted_focal_point = insert_focal_point(I, neutral_value);
							elements_generated_from_support.push_back(inserted_focal_point);
						}
					}
					U = FOD::set_union(U, A);
				}
				// add also emptyset to this->structure if it wasn't already there
				if(!this->definition[emptyset]){
					set_N_value<T>* inserted_focal_point = insert_focal_point(emptyset, neutral_value);
					elements_generated_from_support.push_back(inserted_focal_point);
				}
			}
			return true;
		}


		bool consonance_check(
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


		void to_semilattice_support(
			const std::vector<T>& powerset_values
		) {
			std::unordered_map<T, std::vector<boost::dynamic_bitset<> > > focal_points_map;

			for (size_t n = 0; n < pow(2, this->definition.get_FOD_size()); ++n){
				boost::dynamic_bitset<> set(this->definition.get_FOD_size(), n);

				if (focal_points_map.find(powerset_values[n]) == focal_points_map.end()){
					focal_points_map.emplace(powerset_values[n], (std::vector<boost::dynamic_bitset<> >) {set});
				}else{
					bool new_point = true;
					std::vector<boost::dynamic_bitset<> >& focal_points = focal_points_map[powerset_values[n]];
					for (size_t i = 0; i < focal_points.size(); ++i){
						if (this->order_relation == order_relation_t::subset){
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
					if(!this->definition[kv.second[i]]){
						insert_focal_point(kv.second[i], kv.first);
					}
				}
			}
		}


		void compute_iota_elements(const version_t& version, const powerset_btree<T>& support){
			powerset_btree<size_t> iota_tree(support.get_FOD(), support.get_block_size());

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
							set_N_value<size_t>* inserted_iota_element = iota_tree.insert(iota_element, i);
							std::clog << "\nNEW INSERTED IOTA ELEMENT :\n";
							std::clog << inserted_iota_element->set << std::endl;
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
							set_N_value<size_t>* inserted_iota_element_dual = iota_tree.insert(iota_element_dual, i);
							std::clog << "\nNEW INSERTED IOTA ELEMENT DUAL :\n";
							std::clog << inserted_iota_element_dual->set << std::endl;
						}
					}
				}
			}
			std::unordered_map<size_t, std::vector<set_N_value<size_t>* > > iota_elements_card_map = iota_tree.elements_by_set_cardinality();
			std::vector<size_t> ordered_cardinalities = powerset_btree<size_t>::get_sorted_cardinalities(iota_elements_card_map, *iota_tree.get_FOD());
			this->iota_sequence.reserve(iota_tree.size());
			this->iota_fiber_sequence.reserve(iota_tree.size());

			size_t cc;
			for (size_t c = 0; c < ordered_cardinalities.size(); ++c) {
				if(version == version_t::regular){
					cc = c;
				}else{
					cc = ordered_cardinalities.size()-1 - c;
				}
				const std::vector<set_N_value<size_t>* >& iota_elements_nodes = iota_elements_card_map[ordered_cardinalities[cc]];
				for (size_t i = 0; i < iota_elements_nodes.size(); ++i) {
					std::clog << iota_elements_nodes[i]->set << std::endl;
					this->iota_sequence.emplace_back(iota_elements_nodes[i]->set);
					this->iota_fiber_sequence.emplace_back(iota_elements_nodes[i]->value);
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
		void to_semilattice_support(
			std::unordered_map<size_t, std::vector<set_N_value<T>* > >& support_card_map,
			const std::vector<size_t>& support_ordered_cardinalities,
			std::vector<set_N_value<T>* >& elements_generated_from_support,
			const T& neutral_value
		) {

			std::function<boost::dynamic_bitset<>(const boost::dynamic_bitset<>&, const boost::dynamic_bitset<>&)> domain_binary_operator;
			std::function<std::vector<boost::dynamic_bitset<> >(const powerset_btree<T>& powerset, const boost::dynamic_bitset<>&, size_t)> tree_operator;
			std::function<std::vector<boost::dynamic_bitset<> >(const powerset_btree<set_N_value<T>* >& powerset, const boost::dynamic_bitset<>&, size_t)> tree_operator_dual;

			if(this->order_relation == order_relation_t::subset){
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
			if (support_ordered_cardinalities.back() == this->definition.get_FOD_size())
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
						if(!this->definition[focal_point]){
							set_N_value<T>* inserted_focal_point = insert_focal_point(focal_point, neutral_value);
							elements_generated_from_support.push_back(inserted_focal_point);
						}
					}
					const std::vector<boost::dynamic_bitset<> >& operations_with_not_subsets_of_smaller_than = tree_operator(this->definition, setA, support_ordered_cardinalities[c]-1);

					// search for sets of lower cardinality that are not subsets of the current set
					for (size_t ii = 0; ii < operations_with_not_subsets_of_smaller_than.size(); ++ii) {
						const boost::dynamic_bitset<>& focal_point = operations_with_not_subsets_of_smaller_than[ii];

						// add it to this->structure if it wasn't already there
						if(!this->definition[focal_point]){
							set_N_value<T>* inserted_focal_point = insert_focal_point(focal_point, neutral_value);
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
				const std::vector<boost::dynamic_bitset<> >& operations_with_not_subsets_of_smaller_than = tree_operator(this->definition, setA, setA.count()-1);

				// search for sets of lower cardinality that are not subsets of the current set
				for (size_t ii = 0; ii < operations_with_not_subsets_of_smaller_than.size(); ++ii) {
					const boost::dynamic_bitset<>& focal_point = operations_with_not_subsets_of_smaller_than[ii];

					// add it to this->structure if it wasn't already there
					if(!this->definition[focal_point]){
						set_N_value<T>* inserted_focal_point = insert_focal_point(focal_point, neutral_value);
						elements_generated_from_support.push_back(inserted_focal_point);
					}
				}
				const std::vector<boost::dynamic_bitset<> >& dual_operations_with_not_subsets_of_smaller_than = tree_operator_dual(
						this->dual_definition, ~setA,
						this->definition.get_FOD_size() - setA.count()-1);

				// search for sets of higher cardinality that are not supersets of the current set
				for (size_t ii = 0; ii < dual_operations_with_not_subsets_of_smaller_than.size(); ++ii) {
					const boost::dynamic_bitset<>& dual_focal_point = dual_operations_with_not_subsets_of_smaller_than[ii];

					// add it to this->structure if it wasn't already there
					if(!this->dual_definition[dual_focal_point]){
						set_N_value<T>* inserted_focal_point = insert_dual_focal_point(dual_focal_point, neutral_value);
						elements_generated_from_support.push_back(inserted_focal_point);
					}
				}
				++i;
			}
			std::clog << "\nFocal points found: \n";
			print<T>(std::clog, this->definition);
		}


		/*
		 * In this function, this->iota_sequence is supposed to contain all regular iota elements (i.e. join-irreducible of this lattice).
		 */
		void to_lattice_support(const T& neutral_value){
			if (this->iota_sequence.size() == 0)
				return;

			if(this->order_relation == order_relation_t::subset){
				for (size_t o = 0; o < this->iota_sequence.size(); ++o) {
					const std::vector<set_N_value<T>* >& lattice_elements = this->definition.elements();
					for (size_t i = 0; i < lattice_elements.size(); ++i) {
						boost::dynamic_bitset<> new_set = FOD::set_union(this->iota_sequence[o], lattice_elements[i]->set);
						if(!this->definition[new_set]){
							insert_focal_point(new_set, neutral_value);
						}
					}
				}
			}else{
				for (size_t o = 0; o < this->iota_sequence.size(); ++o) {
					const std::vector<set_N_value<set_N_value<T>* >* >& lattice_elements = this->dual_definition.elements();
					for (size_t i = 0; i < lattice_elements.size(); ++i) {
						boost::dynamic_bitset<> new_set = FOD::set_union(~this->iota_sequence[this->iota_sequence.size()-1-o], lattice_elements[i]->set);
						if(!this->dual_definition[new_set]){
							insert_dual_focal_point(new_set, neutral_value);
						}
					}
				}
			}
			std::clog << "\nCropped lattice support: \n";
			print<T>(std::clog, this->definition);
		}
	};
}	// namespace efficient_DST

#endif // EFFICIENT_DST_ZETA_TRANSFORM_HPP
