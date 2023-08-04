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
 
#ifndef EFFICIENT_DST_ZETA_TRANSFORM_HPP
#define EFFICIENT_DST_ZETA_TRANSFORM_HPP

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include <powerset_function.hpp>
#include <mobius_inversion.hpp>


namespace efficient_DST{

	enum class operation_type_t: bool { addition, multiplication };

	template <class inclusion, size_t N, typename T = float>
	class zeta_transform : public powerset_function<N, T> {
	protected:
		using typename powerset_function<N, T>::subset;
		using powerset_function<N, T>::emptyset;
		using powerset_function<N, T>::fullset;

		scheme_type_t scheme_type;
		std::vector<subset > iota_sequence;
		std::vector<subset > sync_sequence;
		std::unordered_map<subset, size_t> bridge_map;
		std::map<size_t, powerset_btree<N, set_N_value<N, T>* > > definition_by_cardinality;

	public:

		zeta_transform(const zeta_transform<inclusion, N, T>& z) :
			powerset_function<N, T>(z.outcomes, z.definition, z.default_value),
			scheme_type(z.scheme_type),
			iota_sequence(z.iota_sequence),
			sync_sequence(z.sync_sequence),
			bridge_map(z.bridge_map),
			definition_by_cardinality(z.definition_by_cardinality)
		{}

//		zeta_transform(
//			const sample_space<N>& outcomes,
//			const powerset_btree<N, T>& definition,
//			const scheme_type_t& scheme_type,
//			const std::vector<subset >& iota_sequence,
//			const T& neutral_value
//		) :
//			powerset_function<N, T>(outcomes, definition),
//			scheme_type(scheme_type),
//			iota_sequence(iota_sequence),
//			neutral_value(neutral_value)
//		{
//			this->set_definition_by_cardinality();
//		}

		/*
		 * Constructor when you have a Möbius transform such as the mass or conjunctive/disjunctive weight function.
		 * - support is supposed to contain all focal sets and their image.
		 * - order_relation is the order relation of this zeta transform (e.g. commonality->superset or implicability->subset).
		 * - transform_operation is the operation of the zeta transform (e.g. in DST, we usually use the addition on the mass function
		 * and the multiplication on the conjunctive/disjunctive weight function).
		 */
		// TODO: ajouter la possibilité de modifier les iota elements, la bridhe map, etc
		// ou faire ne méthode qui renvoie une zeta transform dans le subspace demandé.
		// Dans rule_conjunctive, calculer la fusion (n'ajouter que les éléments en plus de powerset1) et
		// renvoyer le core, ainsi qu'un vecteur contenant tous les subsets ajoutés
		// (ceux-ci auront déjà été intégrés dans l'arbre powerset12,
		// mais seront utilisés pour calculer l'update du semilattice et de la bridge_map).
		// => à faire dans la fusion de commonalités (de zeta_transforms plus généralement).
		// Dans le même block, on calculera les nouveaux iota elements par simple intersection des iota de chaque
		// Déplacer la méthode de fusion de poids dans weight_function (comme pour celle de mass_function)
		// Le truc avec la zeta transform in subspace ne doit avoir lieu que dans zeta transform en soi, où l'on souhaite
		// fusionner deux zeta transform (et où donc on a accès au semilattices et co)
		zeta_transform(
			const sample_space<N>& outcomes,
			const powerset_btree<N, T>& support,
			const T& default_value,
			operation_type_t operation_type
		) :
//			powerset_function<N, T>(outcomes, N * support.number_of_nodes(), default_value),
			powerset_function<N, T>(outcomes, support, default_value),
			scheme_type(scheme_type_t::direct)
		{
			if (operation_type == operation_type_t::addition){
				this->compute<addition<T>>(support);
			} else {
				this->compute<multiplication<T>>(support);
			}
			this->set_definition_by_cardinality();
		}

		/*
		 * Constructor when you have a Möbius transform such as the mass or conjunctive/disjunctive weight function AND you want a specific computation scheme
		 * - support is supposed to contain all focal sets and their image.
		 * - order_relation is the order relation of this zeta transform (e.g. commonality->superset or implicability->subset).
		 * - transform_operation is the operation of the zeta transform (e.g. in DST, we usually use the addition on the mass function
		 * and the multiplication on the conjunctive/disjunctive weight function).
		 */
		zeta_transform(
			const sample_space<N>& outcomes,
			const powerset_btree<N, T>& support,
			const T& default_value,
			operation_type_t operation_type,
			scheme_type_t scheme_type
		) :
//			powerset_function<N, T>(outcomes, N * support.number_of_nodes(), default_value),
			powerset_function<N, T>(outcomes, support, default_value),
			scheme_type(scheme_type)
		{
			if (operation_type == operation_type_t::addition){
				this->compute<addition<T>, false>(support);
			} else {
				this->compute<multiplication<T>, false>(support);
			}
			this->set_definition_by_cardinality();
		}

//		template<class operation_type>
//		powerset_btree<N, T> inversion() const {
		powerset_btree<N, T> inversion(const operation_type_t& operation_type) const {
			std::cout << "INVERSION\n";
			powerset_btree<N, T> mobius_transform_definition(this->definition);
			if (operation_type == operation_type_t::addition){
				efficient_mobius_inversion<inclusion, mobius_tranformation<inclusion, addition<T>, N, T>, N, T >::execute(
						this->definition,
						mobius_transform_definition,
						this->iota_sequence,
						this->sync_sequence,
						this->bridge_map,
						this->scheme_type
				);
			} else {
				efficient_mobius_inversion<inclusion, mobius_tranformation<inclusion, multiplication<T>, N, T>, N, T >::execute(
						this->definition,
						mobius_transform_definition,
						this->iota_sequence,
						this->sync_sequence,
						this->bridge_map,
						this->scheme_type
				);
			}
			return mobius_transform_definition;
		}

		powerset_btree<N, T> inversion_in_subspace(const operation_type_t& operation_type, subset subspace) const {
			std::cout << "INVERSION\n";
			powerset_btree<N, T> mobius_transform_definition(this->definition);
			if (operation_type == operation_type_t::addition){
				efficient_mobius_inversion<inclusion, mobius_tranformation<inclusion, addition<T>, N, T>, N, T >::execute(
						this->definition,
						mobius_transform_definition,
						this->iota_sequence,
						this->sync_sequence,
						this->bridge_map,
						this->scheme_type
				);
			} else {
				efficient_mobius_inversion<inclusion, mobius_tranformation<inclusion, multiplication<T>, N, T>, N, T >::execute(
						this->definition,
						mobius_transform_definition,
						this->iota_sequence,
						this->sync_sequence,
						this->bridge_map,
						this->scheme_type
				);
			}
			return mobius_transform_definition;
		}

		T at_emptyset() const {
			return (*this)[emptyset];
		}

		T at_fullset() const {
			return (*this)[fullset];
		}

		T operator[](const std::vector<std::string>& labels) const {
			return (*this)[this->outcomes.get_subset(labels)];
		}

		T operator[](const subset& set) const {
			size_t index = this->definition[set];
			if(index < this->definition.number_of_nodes()){
				return this->definition.get_node(index).value;
			}
			return this->find_non_focal_point_image(set);
		}

		zeta_transform<inclusion, N, T> natural_fusion_with(const zeta_transform<inclusion, N, T>& z2) const {
			size_t max_def_size = std::max(this->definition.size(), z2.definition.size());
			zeta_transform<inclusion, N, T> z12(this->outcomes, this->default_value);
			powerset_btree<N, T> w12_definition;
			subset core;
			std::vector<subset> new_sets;
			if (z2.definition.size() == max_def_size){
				z12.scheme_type = z2.scheme_type;
//				z12.sync_sequence = z2.sync_sequence;
//				z12.bridge_map = z2.bridge_map;
//				z12.definition_by_cardinality = z2.definition_by_cardinality;
				const powerset_btree<N, T>& w2 = z2.inversion(operation_type_t::multiplication);
				w12_definition = w2;
				new_sets = weight_function<N, T>::template natural_fusion<inclusion>(
					this->inversion(operation_type_t::multiplication),
					w2,
					w12_definition,
					core
				);
			} else {
				z12.scheme_type = this->scheme_type;
//				z12.sync_sequence = this->sync_sequence;
//				z12.bridge_map = this->bridge_map;
//				z12.definition_by_cardinality = this->definition_by_cardinality;
				const powerset_btree<N, T>& w1 = this->inversion(operation_type_t::multiplication);
				w12_definition = w1;
				new_sets = weight_function<N, T>::template natural_fusion<inclusion>(
					z2.inversion(operation_type_t::multiplication),
					w1,
					w12_definition,
					core
				);
			}
			z12.definition = w12_definition;
			powerset_btree<N, T> iota_elements(N);
			for (size_t i = 0; i < this->iota_sequence; ++i){
				iota_elements.update_or_insert(this->iota_sequence, 0);
			}
			for (size_t i = 0; i < z2.iota_sequence; ++i){
				iota_elements.update_or_insert(z2.iota_sequence, 0);
			}
			switch (z12.scheme_type){
				case scheme_type_t::semilattice:
					inclusion::compute_iota_elements_dual(
						iota_elements,
						z12.iota_sequence
					);
					efficient_mobius_inversion<
						inclusion, zeta_tranformation<inclusion, multiplication<T>, N, T>, N, T
					>::build_sync_sequence(z12.iota_sequence, z12.sync_sequence);
					efficient_mobius_inversion<
						inclusion, zeta_tranformation<inclusion, multiplication<T>, N, T>, N, T
					>::update_semilattice_support(
						z12.definition.elements(),
						z12.definition,
						new_sets,
						1
					);
					efficient_mobius_inversion<
						inclusion, zeta_tranformation<inclusion, multiplication<T>, N, T>, N, T
					>::build_subspace_bridge_map(
						z12.bridge_map,
						z12.definition,
						z12.iota_sequence,
						core
					);
					efficient_mobius_inversion<
						inclusion, zeta_tranformation<inclusion, multiplication<T>, N, T>, N, T
					>::execute_EMT_with_semilattice(
							z12.definition,
							z12.iota_sequence,
							z12.sync_sequence,
							z12.bridge_map
					);
					break;
				case scheme_type_t::lattice:
					inclusion::compute_iota_elements(
						iota_elements,
						z12.iota_sequence
					);
					efficient_mobius_inversion<
						inclusion, zeta_tranformation<inclusion, multiplication<T>, N, T>, N, T
					>::update_truncated_lattice_support(
						z12.iota_sequence,
						z12.definition,
						new_sets,
						1
					);
					efficient_mobius_inversion<
						inclusion, zeta_tranformation<inclusion, multiplication<T>, N, T>, N, T
					>::execute_EMT_with_lattice(
							z12.definition,
							z12.iota_sequence
					);
					break;
				case scheme_type_t::reduced_FMT:
					efficient_mobius_inversion<
						inclusion, zeta_tranformation<inclusion, multiplication<T>, N, T>, N, T
					>::FMT_reduced_to_core(
							support,
							z12.definition,
							z12.iota_sequence
					);
					break;
				default:
					efficient_mobius_inversion<
						inclusion, zeta_tranformation<inclusion, multiplication<T>, N, T>, N, T
					>::direct_transformation(
							support,
							z12.definition,
							z12.scheme_type
					);
					break;
			}
			z12.set_definition_by_cardinality();
//			weight_function<N, T> w12(w1.get_sample_space(), w12_definition);
//			powerset_btree<N, T>& w12_definition = w12.definition_();
//			set_N_value<N, T>* node;
//			T val;
//			for (const auto& set : manifest) {
////				w12.assign(set, w1[set] * w2[set]);
//				val = 1;
//				node = w1_definition[set];
//				if (node){
//					val = node->value;
//				}
//				node = w2_definition[set];
//				if (node){
//					val *= node->value;
//				}
//				if (!powerset_function<N, T>::is_equivalent_to_zero(val - 1))
//					w12_definition.insert(set, val);
//			}
			return z12;
		}

	protected:

		T find_non_focal_point_image(const subset& set) const {
			const size_t& card = set.count();
			std::vector<set_N_value<N, set_N_value<N, T>* >* >& related_elements;
			related_elements.reserve(1);
			for (const auto& c_focal_points : this->definition_by_cardinality) {
				if(inclusion::naturally_ranked(c_focal_points.first, card)){
					const powerset_btree<N, set_N_value<N, T>* >& p = c_focal_points.second;
					inclusion::addresses_related_to(p, set, related_elements);

					if(related_elements.size() > 0){
						return related_elements[0]->value->value;
					}
				}
			}
			return this->default_value;
		}


		void set_definition_by_cardinality(){
			const auto& definition_card_map = inclusion::card_mapping_dual(this->definition);
			// create a powerset_btree in this->definition_by_cardinality for each vector of elements in definition_card_map
			for (const auto& c_focal_points : definition_card_map) {
				this->definition_by_cardinality.emplace(
					std::piecewise_construct,
					std::make_tuple(c_focal_points.first),
					std::make_tuple(1.5 * c_focal_points.second.size())
				);
				powerset_btree<N, set_N_value<N, T>* >& p_c = this->definition_by_cardinality[c_focal_points.first];
				const std::vector<set_N_value<N, T>* >& elements = c_focal_points.second;
				for(size_t i = 0; i < elements.size(); ++i){
					p_c.insert(elements[i]->set, elements[i]);
				}
			}
		}

		template<class operation_type, bool autoset = true>
		void compute(const powerset_btree<N, T>& support){
			if (this->default_value != operation_type::neutral_value()){
				std::cerr << "Ill-defined support: the default value of your compact definition must match the neutral value for the operator you chose.\n";
				exit(1);
			}
			if (autoset){
				this->scheme_type = efficient_mobius_inversion<
					inclusion, zeta_tranformation<inclusion, operation_type, N, T>, N, T
				>::autoset_and_build(
						support,
						this->definition,
						this->iota_sequence,
						this->sync_sequence,
						this->bridge_map
				);
				efficient_mobius_inversion<
					inclusion, zeta_tranformation<inclusion, operation_type, N, T>, N, T
				>::execute(
						support,
						this->definition,
						this->iota_sequence,
						this->sync_sequence,
						this->bridge_map,
						this->scheme_type
				);
			}else{
				switch (this->scheme_type){
					case scheme_type_t::semilattice:
						efficient_mobius_inversion<
							inclusion, zeta_tranformation<inclusion, operation_type, N, T>, N, T
						>::EMT_with_semilattice(
								support,
								this->definition,
								this->iota_sequence,
								this->sync_sequence,
								this->bridge_map
						);
						break;
					case scheme_type_t::lattice:
						efficient_mobius_inversion<
							inclusion, zeta_tranformation<inclusion, operation_type, N, T>, N, T
						>::EMT_with_lattice(
								support,
								this->definition,
								this->iota_sequence
						);
						break;
					case scheme_type_t::reduced_FMT:
						efficient_mobius_inversion<
							inclusion, zeta_tranformation<inclusion, operation_type, N, T>, N, T
						>::FMT_reduced_to_core(
								support,
								this->definition,
								this->iota_sequence
						);
						break;
					default:
						efficient_mobius_inversion<
							inclusion, zeta_tranformation<inclusion, operation_type, N, T>, N, T
						>::direct_transformation(
								support,
								this->definition,
								this->scheme_type
						);
						break;
				}
			}
//			this->default_value = operation_type::neutral_value();
		}

		void normalize(const subset& normalizing_set, const T& normalizing_value){
			size_t normalizing_assignment_index = this->definition[normalizing_set];
			const std::vector<set_N_value<N, T>* >& focal_log_elements = this->definition._elements();
			if(normalizing_assignment_index < this->definition.number_of_nodes()){
				const set_N_value<N, T>& normalizing_assignment = this->definition.get_node(normalizing_assignment_index);
				if (normalizing_assignment.value != normalizing_value){
					for (size_t i = 0; i < focal_log_elements.size(); ++i){
						if(focal_log_elements[i]->set != normalizing_set){
							focal_log_elements[i]->value /= normalizing_assignment.value;
						}
					}
				}
			}
			if(normalizing_assignment_index >= this->definition.number_of_nodes() || this->definition.get_node(normalizing_assignment_index).value != normalizing_value){
				for (size_t i = 0; i < focal_log_elements.size(); ++i){
					if(focal_log_elements[i]->set != normalizing_set){
						focal_log_elements[i]->value *= normalizing_value;
					}
				}
				bool insertion = (this->definition[normalizing_set] >= this->definition.number_of_nodes());
				this->definition.update_or_insert(normalizing_set, normalizing_value);
				if (insertion){
					switch (this->scheme_type){
						case scheme_type_t::semilattice:
							efficient_mobius_inversion<
								inclusion, zeta_tranformation<inclusion, addition<T>, N, T>, N, T
							>::update_semilattice_support(
									this->definition.elements(),
									this->definition,
									{normalizing_set},
									normalizing_value
							);
							this->iota_sequence.clear();
							inclusion::compute_iota_elements_dual(
									this->definition,
									this->iota_sequence
							);
							this->bridge_map.clear();
							efficient_mobius_inversion<
								inclusion, zeta_tranformation<inclusion, addition<T>, N, T>, N, T
							>::build_bridge_map(this->bridge_map, this->definition, this->iota_sequence);
							break;
						case scheme_type_t::lattice:
							iota_sequence.clear();
							inclusion::compute_iota_elements(
									this->definition,
									this->iota_sequence
							);
							efficient_mobius_inversion<
								inclusion, zeta_tranformation<inclusion, addition<T>, N, T>, N, T
							>::update_truncated_lattice_support(
									this->iota_sequence,
									this->definition,
									{normalizing_set},
									normalizing_value
							);
							break;
						case scheme_type_t::reduced_FMT:
							std::cout << "Not implemented !\n";
							break;
						default:
							break;
					}
				}
			}
		}
	};
}	// namespace efficient_DST

#endif // EFFICIENT_DST_ZETA_TRANSFORM_HPP
