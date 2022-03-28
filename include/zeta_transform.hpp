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
		std::map<size_t, powerset_btree<N, set_N_value<N, T>* > > definition_by_cardinality;

	public:

		zeta_transform(const zeta_transform<inclusion, N, T>& z) :
			powerset_function<N, T>(z.outcomes, z.definition, z.default_value),
			scheme_type(z.scheme_type),
			iota_sequence(z.iota_sequence),
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
		zeta_transform(
			const sample_space<N>& outcomes,
			const powerset_btree<N, T>& support,
			const T& default_value,
			operation_type_t operation_type
		) :
			powerset_function<N, T>(outcomes, N * support.number_of_nodes(), default_value),
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
			powerset_function<N, T>(outcomes, N * support.number_of_nodes(), default_value),
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
			powerset_btree<N, T> mobius_transform_definition(this->definition);
			if (operation_type == operation_type_t::addition){
				efficient_mobius_inversion<inclusion, mobius_tranformation<inclusion, addition<T>, N, T>, N, T >::execute(
						this->definition,
						mobius_transform_definition,
						this->iota_sequence,
						this->scheme_type
				);
			} else {
				efficient_mobius_inversion<inclusion, mobius_tranformation<inclusion, multiplication<T>, N, T>, N, T >::execute(
						this->definition,
						mobius_transform_definition,
						this->iota_sequence,
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
						this->iota_sequence
				);
				efficient_mobius_inversion<
					inclusion, zeta_tranformation<inclusion, operation_type, N, T>, N, T
				>::execute(
						support,
						this->definition,
						this->iota_sequence,
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
								this->iota_sequence
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
				this->definition.update_or_insert(normalizing_set, normalizing_value);
			}
		}
	};
}	// namespace efficient_DST

#endif // EFFICIENT_DST_ZETA_TRANSFORM_HPP
