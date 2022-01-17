#ifndef EFFICIENT_DST_ZETA_TRANSFORM_HPP
#define EFFICIENT_DST_ZETA_TRANSFORM_HPP

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include <powerset_function.hpp>
#include <computation_scheme.hpp>


namespace efficient_DST{

	template <class T, size_t N>
	class zeta_transform : public powerset_function<T, N> {
	protected:
		scheme_type_t scheme_type;
		//std::vector<size_t> iota_fiber_sequence;	// fod element indices
		std::vector<std::bitset<N> > iota_sequence;
		powerset_btree<T, N> original_mobius_transform;
		operation_t original_operation;
		std::unordered_map<size_t, powerset_btree<set_N_value<T, N>*, N > > definition_by_cardinality;
		std::vector<size_t> ordered_cardinalities_in_definition;

	public:
		const order_relation_t order_relation;

		zeta_transform(const zeta_transform<T, N>& z) :
			powerset_function<T, N>(z.definition),
			scheme_type(z.scheme_type),
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
				FOD<N>& fod,
				const order_relation_t& order_relation) :
			powerset_function<T, N>(&fod, powerset_values.size()),
			scheme_type(scheme_type_t::semilattice),
			original_mobius_transform(&fod, powerset_function<T, N>::block_size),
			original_operation(operation_t::addition),
			order_relation(order_relation)
		{
//			computation_scheme<T, N>::extract_semilattice_support(
//				powerset_values,
//				this->definition,
//				order_relation
//			);
			computation_scheme<T, N>::set_semilattice_computation_scheme(
					this->definition,
					order_relation,
					this->iota_sequence
			);
			set_definition_by_cardinality();
		}


		/*
		 * Constructor when you have a zeta transform such as the commonality or implicability ones AND know its focal points.
		 * - focal_points_tree is supposed to contain all focal points and their image.
		 * - order_relation is the order relation of this zeta transform (e.g. commonality->superset or implicability->subset).
		 */
		zeta_transform(
				const powerset_btree<T, N>& focal_points_tree,
				const order_relation_t& order_relation) :
			powerset_function<T, N>(focal_points_tree.get_FOD(), focal_points_tree.size()),
			scheme_type(scheme_type_t::semilattice),
			original_mobius_transform(focal_points_tree.get_FOD(), focal_points_tree.get_block_size()),
			original_operation(operation_t::addition),
			order_relation(order_relation)
		{
			computation_scheme<T, N>::set_semilattice_computation_scheme(
					this->definition,
					order_relation,
					this->iota_sequence
					//this->iota_fiber_sequence
			);
			set_definition_by_cardinality();
		}


		/*
		 * Constructor when you have a Möbius transform such as the mass or conjunctive/disjunctive weight function AND you want a specific computation scheme.
		 * - support is supposed to contain all focal sets and their image.
		 * - order_relation is the order relation of this zeta transform (e.g. commonality->superset or implicability->subset).
		 * - transform_operation is the operation of the zeta transform (e.g. in DST, we usually use the addition on the mass function
		 * and the multiplication on the conjunctive/disjunctive weight function)
		 * - scheme_type: computation scheme to use.
		 */
		zeta_transform(
				const powerset_btree<T, N>& support,
				const order_relation_t& order_relation,
				const operation_t& transform_operation,
				const scheme_type_t scheme_type) :
			powerset_function<T, N>(support.get_FOD(), N * support.size()),
			scheme_type(scheme_type),
			original_mobius_transform(support),
			original_operation(transform_operation),
			order_relation(order_relation)
		{
			computation_scheme<T, N>::build_and_execute(
					support,
					this->definition,
					transform_type_t::zeta,
					order_relation,
					transform_operation,
					this->iota_sequence,
					//this->iota_fiber_sequence,
					this->scheme_type
			);
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
				const powerset_btree<T, N>& support,
				const order_relation_t& order_relation,
				const operation_t& transform_operation) :
			powerset_function<T, N>(support.get_FOD(), N * support.size()),
			scheme_type(scheme_type_t::direct),
			original_mobius_transform(support),
			original_operation(transform_operation),
			order_relation(order_relation)
		{
			computation_scheme<T, N>::autoset_and_build(
					support,
					this->definition,
					this->order_relation,
					transform_operation,
					this->iota_sequence,
					//this->iota_fiber_sequence,
					this->scheme_type
			);
			computation_scheme<T, N>::execute(
					this->definition,
					transform_type_t::zeta,
					this->order_relation,
					transform_operation,
					this->iota_sequence,
					//this->iota_fiber_sequence,
					this->scheme_type
			);
			set_definition_by_cardinality();
		}


		powerset_btree<T, N> inversion(const operation_t& transform_operation) const {
			if (transform_operation == this->original_operation && this->original_mobius_transform.size() > 0){
				//return this->original_mobius_transform;
			}
			powerset_btree<T, N> mobius_transform_definition(this->definition);
			computation_scheme<T, N>::execute(
					mobius_transform_definition,
					transform_type_t::Mobius,
					this->order_relation,
					transform_operation,
					this->iota_sequence,
					//this->iota_fiber_sequence,
					this->scheme_type
			);
			return mobius_transform_definition;
		}

		const powerset_btree<T, N>& get_original_mobius_transform() const {
			return this->original_mobius_transform;
		}

		T at_emptyset() const {
			set_N_value<T, N>* set_value = this->definition.sub_fod_of_size(0);
			if(set_value){
				return set_value->value;
			}
			std::bitset<N> emptyset = 0;
			return find_non_focal_point_image(emptyset);
		}

		T at_fod() const {
			set_N_value<T, N>* set_value = this->definition.sub_fod_of_size(N);
			if(set_value){
				return set_value->value;
			}
			std::bitset<N> fod = 0;
			fod.set();
			return find_non_focal_point_image(fod);
		}

		T operator[](const std::vector<std::string>& labels) const {
			return find(this->definition.get_FOD()->to_set(labels));
		}

		T find(const std::bitset<N>& set) const {
			set_N_value<T, N>* set_value = this->definition[set];
			if(set_value){
				return set_value->value;
			}
			return find_non_focal_point_image(set);
		}

	protected:

		T find_non_focal_point_image(const std::bitset<N>& set) const {
			set_N_value<set_N_value<T, N>*, N >* A = nullptr;
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
						const powerset_btree<set_N_value<T, N>*, N >& p = this->definition_by_cardinality.at(this->ordered_cardinalities_in_definition[c]);
						const std::vector<set_N_value<set_N_value<T, N>*, N >* >& subsets = p.subsets_of(set);

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
						const powerset_btree<set_N_value<T, N>*, N >& p = this->definition_by_cardinality.at(this->ordered_cardinalities_in_definition[c]);
						std::vector<set_N_value<set_N_value<T, N>*, N >* > supersets = p.supersets_of(set);

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
//			std::unordered_map<size_t, std::vector<set_N_value<T, N>* > > definition_card_map = this->definition.elements_by_set_cardinality();
//			if (this->order_relation == order_relation_t::subset) {
//				this->definition.get_FOD()->sort_cardinalities(this->ordered_cardinalities_in_definition, definition_card_map, order_t::descending);
//			}else{
//				this->definition.get_FOD()->sort_cardinalities(this->ordered_cardinalities_in_definition, definition_card_map, order_t::ascending);
//			}
//
//			this->definition_by_cardinality.reserve(this->ordered_cardinalities_in_definition.size());
//
//			for(size_t c = 0; c < this->ordered_cardinalities_in_definition.size(); ++c){
//				this->definition_by_cardinality.emplace(std::piecewise_construct, std::make_tuple(this->ordered_cardinalities_in_definition[c]),
//						std::make_tuple(this->definition.get_FOD(), this->definition.get_block_size()));
//				powerset_btree<set_N_value<T, N>*, N >& p_c = this->definition_by_cardinality[this->ordered_cardinalities_in_definition[c]];
//				const std::vector<set_N_value<T, N>* >& elements = definition_card_map[this->ordered_cardinalities_in_definition[c]];
//				for(size_t i = 0; i < elements.size(); ++i){
//					p_c.insert(elements[i]->set, elements[i]);
//				}
//			}
		}
	};
}	// namespace efficient_DST

#endif // EFFICIENT_DST_ZETA_TRANSFORM_HPP
