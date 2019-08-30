#ifndef EFFICIENT_DST_MOBIUS_AGGREGATE_HPP
#define EFFICIENT_DST_MOBIUS_AGGREGATE_HPP

#include <mobius_transform.hpp>
#include <powerset_btree.hpp>
#include <math.h>


namespace efficient_DST{

	enum class scheme_type_t: int8_t { consonant, semilattice, lattice, other };

	template <class T=double>
	class mobius_aggregate : public mobius_transform<T> {
	public:
		scheme_type_t scheme_type;
		const order_relation_t order_relation;
		std::vector<boost::dynamic_bitset<> > sequence;
		powerset_btree<T> dual_definition;
		powerset_btree<T> original_mobius_transform;
		mobius_transformation_form_t original_mobius_transform_form;


		mobius_aggregate(const mobius_aggregate<T>& mobius_aggregation) :
			mobius_transform<T>(mobius_aggregation.definition),
			scheme_type(mobius_aggregation.scheme_type),
			order_relation(mobius_aggregation.order_relation),
			dual_definition(mobius_aggregation.dual_definition),
			original_mobius_transform(mobius_aggregation.original_mobius_transform),
			original_mobius_transform_form(mobius_aggregation.original_mobius_transform_form)
		{}


		/*
		 * Constructor when you have an aggregate function such as the commonality or implicability ones without knowledge about its focal points.
		 * - powerset_values contains all images for the whole powerset of FOD.
		 * - _fod is the corresponding FOD.
		 * - order_relation is the order relation that has been used to aggregate the given function (commonality or implicability).
		 */
		mobius_aggregate(
				const std::vector<T>& powerset_values,
				const FOD& fod,
				const order_relation_t& order_relation) :
			mobius_transform<T>(fod),
			scheme_type(scheme_type_t::other),
			order_relation(order_relation),
			dual_definition(fod, mobius_transform<T>::block_size),
			original_mobius_transform(fod, mobius_transform<T>::block_size),
			original_mobius_transform_form(mobius_transformation_form_t::additive)
		{
			if(powerset_values.size() != pow(2, fod.size())){
				std::cerr << "\nThe given vector does not feature the same size as the powerset of the given FOD.\n";
				return;
			}
			to_lower_semilattice(
				powerset_values,
				order_relation
			);
			this->original_mobius_transform.copy(this->inversion(mobius_transformation_form_t::additive));
			this->original_mobius_transform_form = mobius_transformation_form_t::additive;
		}


		/*
		 * Constructor when you have an aggregate function such as the commonality or implicability ones AND know its focal points.
		 * - structure is supposed to contain all focal points and their image.
		 * - order_relation is the order relation that has been used to aggregate the given function (commonality or implicability).
		 */
		mobius_aggregate(
				const powerset_btree<T>& structure,
				const order_relation_t& order_relation) :
			mobius_transform<T>(structure),
			scheme_type(scheme_type_t::other),
			order_relation(order_relation),
			dual_definition(*structure.fod, structure.block_size),
			original_mobius_transform(*structure.fod, structure.block_size),
			original_mobius_transform_form(mobius_transformation_form_t::additive)
		{
			powerset_btree<T>::reverse_powerset_from_to(this->definition, this->dual_definition);
			this->original_mobius_transform.copy(this->inversion(mobius_transformation_form_t::additive));
			this->original_mobius_transform_form = mobius_transformation_form_t::additive;
		}


		/*
		 * Constructor when you have an original function such as the mass or conjunctive/disjunctive weight ones.
		 * - original_structure is supposed to contain all focal sets and their image.
		 * - order_relation is the order relation that will be used to compute the aggregate function (commonality or implicability) resulting from the call to run().
		 * - mobius_transformation_form is a boolean indicating whether original_structure will be used to compute its
		 * additive Möbius transform (means that original_structure is the mass function) or its
		 * multiplicative Möbius transform (means that original_structure is the conjunctive/disjunctive weight function).
		 */
		mobius_aggregate(
				const powerset_btree<T>& original_structure,
				const order_relation_t& order_relation,
				const mobius_transformation_form_t& mobius_transformation_form) :
			mobius_transform<T>(original_structure),
			scheme_type(scheme_type_t::other),
			order_relation(order_relation),
			dual_definition(*original_structure.fod, original_structure.block_size),
			original_mobius_transform(original_structure),
			original_mobius_transform_form(mobius_transformation_form)
		{
			T default_value;
			if (mobius_transformation_form == mobius_transformation_form_t::additive){
				default_value = 0;
			} else{
				default_value = 1;
			}
			compute_focal_points(original_structure, default_value);
			powerset_btree<T> transformed_structure(this->definition);
			powerset_btree<T> transformed_structure_dual(this->dual_definition);
			execute_scheme(
					mobius_transformation_form,
					false,
					transformed_structure,
					transformed_structure_dual);
			this->definition.copy(transformed_structure);
			this->dual_definition.copy(transformed_structure_dual);
		}


		powerset_btree<T> inversion(const mobius_transformation_form_t& mobius_transformation_form) const {
			if (mobius_transformation_form == this->original_mobius_transform_form && this->original_mobius_transform.size() > 0){
				//return this->original_mobius_transform;
			}
			powerset_btree<T> transformed_structure(this->definition);
			powerset_btree<T> transformed_structure_dual(this->dual_definition);
			execute_scheme(
					mobius_transformation_form,
					true,
					transformed_structure,
					transformed_structure_dual);
			return transformed_structure;
		}

		const powerset_btree<T>& get_original_mobius_transform() const {
			return this->original_mobius_transform;
		}

		T at_emptyset() const {
			set_N_value<T>* set_value = this->definition.sub_fod_of_size(0);
			if(set_value){
				return set_value->value;
			}
			boost::dynamic_bitset<> emptyset(this->definition.fod->size());
			return compute_aggregation(emptyset);
		}

		T at_fod() const {
			set_N_value<T>* set_value = this->definition.sub_fod_of_size(this->definition.fod->size());
			if(set_value){
				return set_value->value;
			}
			boost::dynamic_bitset<> fod(this->definition.fod->size());
			fod.set();
			return compute_aggregation(fod);
		}

		T operator[](const std::vector<std::string>& labels) const {
			return find(this->definition.fod->to_set(labels));
		}

		T find(const boost::dynamic_bitset<>& set) const {
			set_N_value<T>* set_value = this->definition[set];
			if(set_value){
				return set_value->value;
			}
			return compute_aggregation(set);
		}

	protected:

		T compute_aggregation(const boost::dynamic_bitset<>& set) const {
			T sum = 0;
			std::vector<set_N_value<T>* > elements;
			std::function<T(const T&, const T&)> range_binary_operator;

			if (this->order_relation == order_relation_t::subset){
				elements = this->original_mobius_transform.subsets_of(set);
			} else {
				elements = this->original_mobius_transform.supersets_of(set);
			}
			if (this->original_mobius_transform_form == mobius_transformation_form_t::additive){
				range_binary_operator = addition;
			} else {
				range_binary_operator = multiplication;
			}

			for (size_t i = 0; i < elements.size(); ++i) {
				sum = range_binary_operator(sum, elements[i]->value);
			}
			return sum;
		}


		void compute_focal_points(
				const powerset_btree<T>& original_structure,
				const T& default_value
				){
			powerset_btree<T>::reverse_powerset_from_to(original_structure, this->dual_definition);
			powerset_btree<T> original_dual_structure(this->dual_definition);

			// Cardinality map of focal sets in original_structure
			std::unordered_map<size_t, std::vector<set_N_value<T>* > > original_structure_card_map = original_structure.elements_by_set_cardinality();
			std::vector<set_N_value<T>* > new_elements_in_structure;

			// check if original_structure is almost Bayesian
			std::clog << "Linear analysis:" << std::endl;
			const bool& is_almost_bayesian = linear_analysis_of_structure(
				original_structure_card_map,
				new_elements_in_structure,
				default_value
			);

			if(is_almost_bayesian){
				std::clog << "almost Bayesian." << std::endl;
				this->scheme_type = scheme_type_t::other;
			}else{
				// sort focal sets in original_structure by increasing cardinality
				const std::vector<size_t>& original_ordered_cardinalities = powerset_btree<T>::get_sorted_cardinalities(original_structure_card_map, *original_structure.fod);

				std::clog << "Consonance check:" << std::endl;
				const bool& is_consonant = consonance_check(
					new_elements_in_structure,
					original_structure_card_map,
					original_ordered_cardinalities
				);

				if(is_consonant){
					std::clog << "consonant." << std::endl;
					this->scheme_type = scheme_type_t::consonant;
				}else{
					std::clog << "not consonant." << std::endl;
					if(original_structure.size() < 5 * original_structure.fod->size()){
						std::clog << "Number of focal sets equivalent to |FOD|." << std::endl;
						std::clog << "Transform to semilattice:" << std::endl;
						to_semilattice(
							original_structure_card_map,
							original_ordered_cardinalities,
							new_elements_in_structure,
							default_value
						);
						if(this->definition.size() < 5 * original_structure.fod->size()){
							std::clog << "Number of focal points also equivalent to |FOD|." << std::endl;
							this->scheme_type = scheme_type_t::other;
						}else{
							std::clog << "Number of focal points superior to |FOD|." << std::endl;
							if(this->order_relation == order_relation_t::subset){
								compute_core_sequence(original_dual_structure);
							}else{
								compute_core_sequence(original_structure);
							}
							this->scheme_type = scheme_type_t::semilattice;
						}
					}else{
						std::clog << "Number of focal sets superior to |FOD|." << std::endl;
						std::clog << "Transform to lattice:" << std::endl;
						if(this->order_relation == order_relation_t::subset){
							compute_focal_atoms(original_structure);
							to_lattice(false, default_value);
						}else{
							compute_focal_atoms(original_dual_structure);
							to_lattice(true, default_value);
							this->sequence.clear();
							compute_focal_atoms(original_structure);
						}
						this->scheme_type = scheme_type_t::lattice;
					}
				}
			}
		}


		void execute_scheme(
				const mobius_transformation_form_t& mobius_transformation_form,
				const bool& inversion,
				powerset_btree<T>& working_structure,
				powerset_btree<T>& working_structure_dual
				) const {
			T default_value;
			std::function<T(const T&, const T&)> range_binary_operator;

			if (mobius_transformation_form == mobius_transformation_form_t::additive){
				if (inversion)
					range_binary_operator = subtraction;
				else
					range_binary_operator = addition;
				default_value = 0;
			} else{
				if (inversion)
					range_binary_operator = division;
				else
					range_binary_operator = multiplication;
				default_value = 1;
			}

			switch (this->scheme_type) {
				case scheme_type_t::other:
					std::clog << "other" << std::endl;
					if (inversion){
						directly_compute_inversion(
								range_binary_operator,
								working_structure);
					}else{
						directly_compute_aggregations(
								range_binary_operator,
								working_structure);
					}
					break;

				case scheme_type_t::consonant:
					std::clog << "consonant" << std::endl;
					consonant_computations(
							working_structure,
							inversion,
							range_binary_operator,
							default_value);
					break;

				case scheme_type_t::semilattice:
					std::clog << "semilattice" << std::endl;
					if (this->order_relation == order_relation_t::subset){
						EMT_with_lower_semilattice(
								working_structure_dual,
								this->sequence,
								inversion,
								range_binary_operator);
						powerset_btree<T>::reverse_powerset_from_to(working_structure_dual, working_structure);
					}else{
						EMT_with_lower_semilattice(
								working_structure,
								this->sequence,
								inversion,
								range_binary_operator);
					}
					break;

				case scheme_type_t::lattice:
					std::clog << "lattice" << std::endl;
					EMT_with_lattice(
							working_structure,
							this->sequence,
							inversion,
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


		void directly_compute_aggregations(
				std::function<T(const T&, const T&)> range_binary_operator,
				powerset_btree<T>& transformed_structure
				) const {
			powerset_btree<T> initial_structure(transformed_structure);
			T val;
			const std::vector<set_N_value<T>* >& structure_elements = initial_structure.elements();
			std::vector<set_N_value<T>* > elements;

			for (size_t i = 0; i < structure_elements.size(); ++i) {
				val = structure_elements[i]->value;
				if (this->order_relation == order_relation_t::subset){
					elements = initial_structure.strict_subsets_of(structure_elements[i]->set);
				} else{
					elements = initial_structure.strict_supersets_of(structure_elements[i]->set);
				}
				for (size_t ii = 0; ii < elements.size(); ++ii) {
					val = range_binary_operator(val, elements[ii]->value);
				}
				transformed_structure.insert(structure_elements[i]->set, val);
			}
		}


		void directly_compute_inversion(
				std::function<T(const T&, const T&)> range_binary_operator,
				powerset_btree<T>& transformed_structure
				) const {
			T val;
			std::unordered_map<size_t, std::vector<set_N_value<T>* > > structure_card_map = transformed_structure.elements_by_set_cardinality();
			std::vector<size_t> ordered_cardinalities = powerset_btree<T>::get_sorted_cardinalities(structure_card_map, *transformed_structure.fod);
			std::vector<set_N_value<T>* > elements;

			if (this->order_relation == order_relation_t::superset) {
				std::reverse(ordered_cardinalities.begin(), ordered_cardinalities.end());
			}
			for (size_t c = 0; c < ordered_cardinalities.size(); ++c) {
				const std::vector<set_N_value<T>* >& structure_elements = structure_card_map[ordered_cardinalities[c]];
				for (size_t i = 0; i < structure_elements.size(); ++i) {
					val = structure_elements[i]->value;
					if (this->order_relation == order_relation_t::subset){
						elements = transformed_structure.strict_subsets_of(structure_elements[i]->set);
					} else{
						elements = transformed_structure.strict_supersets_of(structure_elements[i]->set);
					}
					for (size_t ii = 0; ii < elements.size(); ++ii) {
						val = range_binary_operator(val, elements[ii]->value);
					}
					transformed_structure.insert(structure_elements[i]->set, val);
				}
			}
		}


		void consonant_computations(
				powerset_btree<T>& focal_points_tree,
				bool inversion,
				std::function<T(const T&, const T&)> range_binary_operator,
				const T& default_value
				) const {
			std::unordered_map<size_t, std::vector<set_N_value<T>* > > structure_card_map = focal_points_tree.elements_by_set_cardinality();
			std::vector<size_t> ordered_cardinalities = powerset_btree<T>::get_sorted_cardinalities(structure_card_map, *focal_points_tree.fod);
			T value, preceding_value = default_value;

			if (this->order_relation == order_relation_t::superset) {
				std::reverse(ordered_cardinalities.begin(), ordered_cardinalities.end());
			}
			for (size_t i = 0; i < ordered_cardinalities.size(); ++i) {
				set_N_value<T>* set_value = structure_card_map[ordered_cardinalities[i]][0];
				value = range_binary_operator(set_value->value, preceding_value);
				if (inversion){
					preceding_value = set_value->value;
					set_value->value = value;
				}else{
					set_value->value = value;
					preceding_value = set_value->value;
				}
			}
		}


		void EMT_with_lower_semilattice(
				powerset_btree<T>& focal_points_tree,
				const std::vector<boost::dynamic_bitset<> >& sequence,
				const bool& inversion,
				std::function<T(const T&, const T&)> range_binary_operator
				) const {
			powerset_btree<set_N_value<T>* > focal_points_superset_map = focal_points_tree.superset_map();
			const std::vector<set_N_value<T>* >& focal_points = focal_points_tree.elements();
			boost::dynamic_bitset<> fod_cum(focal_points_tree.fod->size());

			if (inversion){
				fod_cum.set();
			}

			for (size_t o = 0; o < sequence.size(); ++o){

				if (!inversion)
					fod_cum = FOD::set_union(fod_cum, sequence[o]);

				for (size_t i = 0; i < focal_points.size(); ++i){
					const boost::dynamic_bitset<>& B = FOD::set_union(focal_points[i]->set, sequence[o]);

					if (B != focal_points[i]->set){
						const set_N_value<set_N_value<T>* >* X = focal_points_superset_map.closest_node_containing(B);

						if (X && FOD::is_or_is_subset_of(X->set, FOD::set_union(focal_points[i]->set, fod_cum))){
							focal_points[i]->value = range_binary_operator(focal_points[i]->value, X->value->value);
						}
					}
				}

				if (inversion)
					fod_cum = FOD::set_minus(fod_cum, sequence[o]);
			}
		}


		void EMT_with_lattice(
				powerset_btree<T>& focal_points_tree,
				const std::vector<boost::dynamic_bitset<> >& sequence,
				bool inversion,
				std::function<T(const T&, const T&)> range_binary_operator
				) const {
			const std::vector<set_N_value<T>* >& focal_points = focal_points_tree.elements();
			std::vector<boost::dynamic_bitset<> > used_sequence(sequence);

			if (inversion) {
				std::reverse(used_sequence.begin(), used_sequence.end());
			}

			for (size_t o = 0; o < used_sequence.size(); ++o){
				for (size_t i = 0; i < focal_points.size(); ++i){
					const boost::dynamic_bitset<>& set_B = FOD::set_union(focal_points[i]->set, used_sequence[o]);

					if (set_B != focal_points[i]->set){
						set_N_value<T>* B = focal_points_tree[set_B];

						if (B){
							if (this->order_relation == order_relation_t::subset){
								B->value = range_binary_operator(B->value, focal_points[i]->value);
							}else{
								focal_points[i]->value = range_binary_operator(focal_points[i]->value, B->value);
							}
						}
					}
				}
			}
		}


		set_N_value<T>* insert_focal_point(const boost::dynamic_bitset<>& focal_point, const T& default_value) {
			set_N_value<T>* inserted_focal_point = this->definition.insert(focal_point, default_value);
			this->dual_definition.insert(~focal_point, default_value);
			std::clog << "\nNEW INSERTED FOCAL POINT :\n";
			std::clog << inserted_focal_point->set << std::endl;
			return inserted_focal_point;
		}

		set_N_value<T>* insert_dual_focal_point(const boost::dynamic_bitset<>& dual_focal_point, const T& default_value) {
			this->dual_definition.insert(dual_focal_point, default_value);
			set_N_value<T>* inserted_focal_point = this->definition.insert(~dual_focal_point, default_value);
			std::clog << "\nNEW INSERTED FOCAL POINT :\n";
			std::clog << inserted_focal_point->set << std::endl;
			return inserted_focal_point;
		}


		bool linear_analysis_of_structure(
				std::unordered_map<size_t, std::vector<set_N_value<T>* > >& structure_card_map,
				std::vector<set_N_value<T>* >& new_elements_in_structure,
				const T& default_value
				){

			if(this->order_relation == order_relation_t::subset){
				std::vector<set_N_value<T>* > original_structure_except_emptyset;
				original_structure_except_emptyset.reserve(this->definition.size()-1);

				for(auto kv : structure_card_map) {
					if(kv.first > 0){
						for (size_t i = 0; i < kv.second.size(); ++i) {
							original_structure_except_emptyset.emplace_back(kv.second[i]);
						}
					}
				}

				boost::dynamic_bitset<> fod(this->definition.fod->size());
				fod.set();

				boost::dynamic_bitset<> neg_U = original_structure_except_emptyset[0]->set;
				std::clog << "\nneg_U = "<< original_structure_except_emptyset[0]->set;

				for (size_t i = 1; i < original_structure_except_emptyset.size(); ++i) {
					const boost::dynamic_bitset<>& A = original_structure_except_emptyset[i]->set;

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
							set_N_value<T>* inserted_focal_point = insert_focal_point(neg_I, default_value);
							new_elements_in_structure.push_back(inserted_focal_point);
						}
					}
					neg_U = FOD::set_intersection(neg_U, A);
				}
				// add also fod to this->structure if it wasn't already there
				if(!this->definition[fod]){
					set_N_value<T>* inserted_focal_point = insert_focal_point(fod, default_value);
					new_elements_in_structure.push_back(inserted_focal_point);
				}
			}else{
				std::vector<set_N_value<T>* > original_structure_except_FOD;
				original_structure_except_FOD.reserve(this->definition.size()-1);

				for(auto kv : structure_card_map) {
					if(kv.first < this->definition.fod->size()){
						for (size_t i = 0; i < kv.second.size(); ++i) {
							original_structure_except_FOD.emplace_back(kv.second[i]);
						}
					}
				}

				const boost::dynamic_bitset<> emptyset(this->definition.fod->size());

				boost::dynamic_bitset<> U = original_structure_except_FOD[0]->set;
				std::clog << "\nU = "<< original_structure_except_FOD[0]->set;

				for (size_t i = 1; i < original_structure_except_FOD.size(); ++i) {
					const boost::dynamic_bitset<>& A = original_structure_except_FOD[i]->set;

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
							set_N_value<T>* inserted_focal_point = insert_focal_point(I, default_value);
							new_elements_in_structure.push_back(inserted_focal_point);
						}
					}
					U = FOD::set_union(U, A);
				}
				// add also emptyset to this->structure if it wasn't already there
				if(!this->definition[emptyset]){
					set_N_value<T>* inserted_focal_point = insert_focal_point(emptyset, default_value);
					new_elements_in_structure.push_back(inserted_focal_point);
				}
			}
			return true;
		}


		bool consonance_check(
				const std::vector<set_N_value<T>* >& new_elements_in_structure,
				std::unordered_map<size_t, std::vector<set_N_value<T>* > >& structure_card_map,
				const std::vector<size_t>& ordered_cardinalities
				){

			// if at least one focal point has already been discovered
			// OR if there are at least two elements with the same 0-indexed cardinality (i.e. the least one) in structure,
			// then structure cannot be consonant
			if(new_elements_in_structure.size() > 0 || structure_card_map[ordered_cardinalities[0]].size() > 1)
				return false;

			for (size_t i = 1; i < ordered_cardinalities.size(); ++i) {
				// if there are at least two elements with the same i-indexed cardinality in structure
				// OR if the set of cardinality index i-1 is not a subset of the one of cardinality index i,
				// then structure cannot be consonant
				if(structure_card_map[ordered_cardinalities[i]].size() > 1 || !FOD::is_or_is_subset_of(
						structure_card_map[ordered_cardinalities[i-1]][0]->set,
						structure_card_map[ordered_cardinalities[i]][0]->set
					  )
				){
					return false;
				}
			}

			return true;
		}


		/*
		 * I use here a notion which I called "focal point" which is the intersection (resp. union)
		 * of a set of focal sets for the superset order relation (resp. subset order relation).
		 * The image of a focal point through q (resp. b) is unique among its supersets (resp. subsets).
		 * A focal set is also a focal point.
		 *
		 * Compute all focal points.
		 * If F is the number of focal sets and dot_F the number of focal points, then there are O(F.dot_F) operations,
		 * where dot_F is in [F, 2^N], where N is the FOD size.
		 * The worst case is obtained if all sets of cardinality N-1 (resp. all singletons) are focal sets
		 * for the superset order relation (resp. subset order relation).
		 */
		void to_semilattice(
			std::unordered_map<size_t, std::vector<set_N_value<T>* > >& original_structure_card_map,
			const std::vector<size_t>& original_ordered_cardinalities,
			std::vector<set_N_value<T>* >& new_elements_in_structure,
			const T& default_value
		) {

			std::function<boost::dynamic_bitset<>(const boost::dynamic_bitset<>&, const boost::dynamic_bitset<>&)> domain_binary_operator;
			std::function<std::vector<boost::dynamic_bitset<> >(const powerset_btree<T>& powerset, const boost::dynamic_bitset<>&, size_t)> tree_operator;
			std::function<std::vector<boost::dynamic_bitset<> >(const powerset_btree<T>& powerset, const boost::dynamic_bitset<>&, size_t)> tree_operator_dual;

			if(this->order_relation == order_relation_t::subset){
				domain_binary_operator = FOD::set_union;
				tree_operator = powerset_btree<T>::unions_with_not_subsets_of_smaller_than;
				tree_operator_dual = powerset_btree<T>::intersections_with_not_subsets_of_smaller_than;
			}else{
				domain_binary_operator = FOD::set_intersection;
				tree_operator = powerset_btree<T>::intersections_with_not_subsets_of_smaller_than;
				tree_operator_dual = powerset_btree<T>::unions_with_not_subsets_of_smaller_than;
			}

			// avoid neutral/absorbing sets for our binary_operator
			size_t c_init = 0, c_end = original_ordered_cardinalities.size();
			if (original_ordered_cardinalities[0] == 0)
				c_init = 1;
			if (original_ordered_cardinalities.back() == this->definition.fod->size())
				--c_end;

			// for each cardinality of focal set, from the smallest to the biggest, except emptyset and FOD
			for (size_t c = c_init; c < c_end; ++c) {
				for (size_t i = 0; i < original_structure_card_map[original_ordered_cardinalities[c]].size(); ++i) {
					const boost::dynamic_bitset<>& setA = original_structure_card_map[original_ordered_cardinalities[c]][i]->set;

					// search for sets of same cardinality
					for (size_t ii = 0; ii < i; ++ii) {
						const boost::dynamic_bitset<>& setB = original_structure_card_map[original_ordered_cardinalities[c]][ii]->set;

						// compute their focal point
						const boost::dynamic_bitset<>& focal_point = domain_binary_operator(setA, setB);

						// add it to this->structure if it wasn't already there
						if(!this->definition[focal_point]){
							set_N_value<T>* inserted_focal_point = insert_focal_point(focal_point, default_value);
							new_elements_in_structure.push_back(inserted_focal_point);
						}
					}
					const std::vector<boost::dynamic_bitset<> >& operations_with_not_subsets_of_smaller_than = tree_operator(this->definition, setA, original_ordered_cardinalities[c]-1);

					// search for sets of lower cardinality that are not subsets of the current set
					for (size_t ii = 0; ii < operations_with_not_subsets_of_smaller_than.size(); ++ii) {
						const boost::dynamic_bitset<>& focal_point = operations_with_not_subsets_of_smaller_than[ii];

						// add it to this->structure if it wasn't already there
						if(!this->definition[focal_point]){
							set_N_value<T>* inserted_focal_point = insert_focal_point(focal_point, default_value);
							new_elements_in_structure.push_back(inserted_focal_point);
						}
					}
				}
			}
			// for each pure focal point (focal point that is not a focal set)
			for (size_t i = 0; i < new_elements_in_structure.size(); ++i) {
				const boost::dynamic_bitset<>& setA = new_elements_in_structure[i]->set;
				const std::vector<boost::dynamic_bitset<> >& operations_with_not_subsets_of_smaller_than = tree_operator(this->definition, setA, setA.count()-1);

				// search for sets of lower cardinality that are not subsets of the current set
				for (size_t ii = 0; ii < operations_with_not_subsets_of_smaller_than.size(); ++ii) {
					const boost::dynamic_bitset<>& focal_point = operations_with_not_subsets_of_smaller_than[ii];

					// add it to this->structure if it wasn't already there
					if(!this->definition[focal_point]){
						set_N_value<T>* inserted_focal_point = insert_focal_point(focal_point, default_value);
						new_elements_in_structure.push_back(inserted_focal_point);
					}
				}
				const std::vector<boost::dynamic_bitset<> >& dual_operations_with_not_subsets_of_smaller_than = tree_operator_dual(
						this->dual_definition, ~setA,
						this->definition.fod->size() - setA.count()-1);

				// search for sets of higher cardinality that are not supersets of the current set
				for (size_t ii = 0; ii < dual_operations_with_not_subsets_of_smaller_than.size(); ++ii) {
					const boost::dynamic_bitset<>& dual_focal_point = dual_operations_with_not_subsets_of_smaller_than[ii];

					// add it to this->structure if it wasn't already there
					if(!this->dual_definition[dual_focal_point]){
						set_N_value<T>* inserted_focal_point = insert_dual_focal_point(dual_focal_point, default_value);
						new_elements_in_structure.push_back(inserted_focal_point);
					}
				}
			}
			std::clog << "\nFocal points found: \n";
			print<T>(std::clog, this->definition);
		}


		void to_semilattice(
			const std::vector<T>& powerset_values
		) {

		}


		void compute_core_sequence(const powerset_btree<T>& original_structure){
			const std::vector<set_N_value<T>* >& focal_elements = original_structure.elements();
			boost::dynamic_bitset<> core = focal_elements[0]->set;
			for (size_t i = 1; i < focal_elements.size(); ++i){
				core = FOD::set_union(core, focal_elements[i]->set);
			}
			this->sequence.reserve(original_structure.fod->size());
			for(size_t i = 0; i < core.size(); ++i){
				if(core[i]){
					this->sequence.emplace_back(original_structure.fod->size());
					this->sequence.back().set(i);
				}
			}
		}


		void compute_focal_atoms(const powerset_btree<T>& original_structure){
			powerset_btree<bool> focal_atoms_tree(*original_structure.fod, original_structure.block_size);
			const FOD& fod = *original_structure.fod;

			for (size_t i = 0; i < fod.size(); ++i) {
				boost::dynamic_bitset<> singleton(fod.size());
				singleton.set(i);
				const std::vector<set_N_value<T>* >& focal_supersets = original_structure.supersets_of(singleton);

				if (focal_supersets.size() > 0) {
					boost::dynamic_bitset<> focal_atom((const boost::dynamic_bitset<>&) focal_supersets[0]->set);

					for (size_t ii = 1; ii < focal_supersets.size(); ++ii) {
						focal_atom = FOD::set_intersection(focal_atom, focal_supersets[ii]->set);
						if (focal_atom == singleton) {
							break;
						}
					}
					if (!focal_atoms_tree[focal_atom]){
						set_N_value<bool>* inserted_focal_atom = focal_atoms_tree.insert(focal_atom, true);
						std::clog << "\nNEW INSERTED FOCAL ATOM :\n";
						std::clog << inserted_focal_atom->set << std::endl;
					}
				}
			}
			std::unordered_map<size_t, std::vector<set_N_value<bool>* > > focal_atoms_card_map = focal_atoms_tree.elements_by_set_cardinality();
			const std::vector<size_t>& ordered_cardinalities = powerset_btree<bool>::get_sorted_cardinalities(focal_atoms_card_map, *focal_atoms_tree.fod);
			this->sequence.reserve(focal_atoms_tree.size());

			for (size_t c = 0; c < ordered_cardinalities.size(); ++c) {
				const std::vector<set_N_value<bool>* >& focal_atoms_nodes = focal_atoms_card_map[ordered_cardinalities[c]];
				for (size_t i = 0; i < focal_atoms_nodes.size(); ++i) {
					this->sequence.emplace_back(focal_atoms_nodes[i]->set);
				}
			}

			if (this->order_relation == order_relation_t::subset){
				std::reverse(this->sequence.begin(), this->sequence.end());
			}
		}


		void to_lattice(const bool& dual, const T& default_value){
			if (this->sequence.size() == 0)
				return;

			if(this->order_relation == order_relation_t::subset || dual){
				powerset_btree<T>* working_tree;
				if(dual)
					working_tree = &this->dual_definition;
				else
					working_tree = &this->definition;

				for (size_t o = 0; o < this->sequence.size(); ++o) {
					const std::vector<set_N_value<T>* >& lattice_elements = working_tree->elements();
					for (size_t i = 0; i < lattice_elements.size(); ++i) {
						boost::dynamic_bitset<> new_set = FOD::set_union(this->sequence[o], lattice_elements[i]->set);
						if(dual){
							if(!this->dual_definition[new_set]){
								// this->dual_structure.insert(new_set, this->default_value);
								std::clog << "!!!!!!!!!!!!!!!!!!!!!! = " << new_set << std::endl;
								insert_dual_focal_point(new_set, default_value);
							}
						}else{
							if(!this->definition[new_set]){
								// this->structure.insert(new_set, this->default_value);
								insert_focal_point(new_set, default_value);
							}
						}
					}
				}
			}else{
				// WARNING: less efficient than for the order relation subset

				// compute the intersection of all focal atoms
				boost::dynamic_bitset<> intersection_of_all((const boost::dynamic_bitset<>&) this->sequence[0]);
				for (size_t o = 1; o < this->sequence.size(); ++o) {
					intersection_of_all = FOD::set_intersection(intersection_of_all, this->sequence[o]);
				}

				if(!this->definition[intersection_of_all]){
					insert_focal_point(intersection_of_all, default_value);
				}

				for (size_t o = 0; o < this->sequence.size(); ++o) {
					const std::vector<set_N_value<T>* >& lattice_elements = this->definition.elements();
					for (size_t i = 0; i < lattice_elements.size(); ++i) {
						boost::dynamic_bitset<> new_set = FOD::set_union(this->sequence[o], lattice_elements[i]->set);
						if(!this->definition[new_set]){
							//this->structure.insert(new_set, this->default_value);
							insert_focal_point(new_set, default_value);
						}
					}
				}
			}
		}
	};
}	// namespace efficient_DST

#endif // EFFICIENT_DST_MOBIUS_AGGREGATE_HPP
