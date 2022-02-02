#ifndef EFFICIENT_DST_ZETA_TRANSFORM_HPP
#define EFFICIENT_DST_ZETA_TRANSFORM_HPP

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include <powerset_function.hpp>
#include <mobius_inversion.hpp>


namespace efficient_DST{

	enum class operation_type_t: bool { addition, multiplication };

	template <class T, size_t N, class inclusion>
	class zeta_transform : public powerset_function<T, N> {
	protected:
		scheme_type_t scheme_type;
		std::vector<std::bitset<N> > iota_sequence;
		std::map<size_t, powerset_btree<set_N_value<T, N>*, N > > definition_by_cardinality;
		T neutral_value;

	public:

		zeta_transform(const zeta_transform<T, N, inclusion>& z) :
			powerset_function<T, N>(z.definition),
			scheme_type(z.scheme_type),
			iota_sequence(z.iota_sequence),
			definition_by_cardinality(z.definition_by_cardinality),
			neutral_value(z.neutral_value)
		{}

		zeta_transform(
			const powerset_btree<T, N>& definition,
			const scheme_type_t& scheme_type,
			const std::vector<std::bitset<N> >& iota_sequence,
			const T& neutral_value
		) :
			powerset_function<T, N>(definition),
			scheme_type(scheme_type),
			iota_sequence(iota_sequence),
			neutral_value(neutral_value)
		{
			set_definition_by_cardinality();
		}

		/*
		 * Constructor when you have a Möbius transform such as the mass or conjunctive/disjunctive weight function.
		 * - support is supposed to contain all focal sets and their image.
		 * - order_relation is the order relation of this zeta transform (e.g. commonality->superset or implicability->subset).
		 * - transform_operation is the operation of the zeta transform (e.g. in DST, we usually use the addition on the mass function
		 * and the multiplication on the conjunctive/disjunctive weight function).
		 */
		zeta_transform(const powerset_btree<T, N>& support, operation_type_t operation_type) :
			powerset_function<T, N>(support.get_FOD(), N * support.size()),
			scheme_type(scheme_type_t::direct)
		{
			if (operation_type == operation_type_t::addition){
				this->compute<addition<T>>(support);
			} else {
				this->compute<multiplication<T>>(support);
			}
			set_definition_by_cardinality();
		}


		powerset_btree<T, N> inversion(const operation_type_t& operation_type) const {
			powerset_btree<T, N> mobius_transform_definition(this->definition);
			if (operation_type == operation_type_t::addition){
				efficient_mobius_inversion<T, N, inclusion, mobius_tranformation<T, N, inclusion, addition<T>> >::execute(
						mobius_transform_definition,
						this->iota_sequence,
						this->scheme_type
				);
			} else {
				efficient_mobius_inversion<T, N, inclusion, mobius_tranformation<T, N, inclusion, multiplication<T>> >::execute(
						mobius_transform_definition,
						this->iota_sequence,
						this->scheme_type
				);
			}
			return mobius_transform_definition;
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
			return this->find(this->definition.get_FOD()->to_set(labels));
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
			const size_t& card = set.count();

			for (const auto& c_focal_points : this->definition_by_cardinality) {
				if(inclusion::naturally_ranked(c_focal_points.first, card)){
					const powerset_btree<set_N_value<T, N>*, N >& p = c_focal_points.second;
					const std::vector<set_N_value<set_N_value<T, N>*, N >* >& related_elements = inclusion::addresses_related_to(p, set);

					if(related_elements.size() > 0){
						A = related_elements[0];
						break;
					}
				}
			}
			if (A){
				return A->value->value;
			}else{
				return this->neutral_value;
			}
		}


		void set_definition_by_cardinality(){
			const auto& definition_card_map = inclusion::card_mapping_dual(this->definition);
			// create a powerset_btree in this->definition_by_cardinality for each vector of elements in definition_card_map
			for (const auto& c_focal_points : definition_card_map) {
				this->definition_by_cardinality.emplace(
					std::piecewise_construct,
					std::make_tuple(c_focal_points.first),
					std::make_tuple(this->definition.get_FOD(), this->definition.get_block_size())
				);
				powerset_btree<set_N_value<T, N>*, N >& p_c = this->definition_by_cardinality[c_focal_points.first];
				const std::vector<set_N_value<T, N>* >& elements = c_focal_points.second;
				for(size_t i = 0; i < elements.size(); ++i){
					p_c.insert(elements[i]->set, elements[i]);
				}
			}
		}

		template<class operation_type>
		void compute(const powerset_btree<T, N>& support){
			this->scheme_type = efficient_mobius_inversion<
				T, N, inclusion, zeta_tranformation<T, N, inclusion, operation_type>
			>::autoset_and_build(
					support,
					this->definition,
					this->iota_sequence
			);

			efficient_mobius_inversion<
				T, N, inclusion, zeta_tranformation<T, N, inclusion, operation_type>
			>::execute(
					this->definition,
					this->iota_sequence,
					this->scheme_type
			);
			this->neutral_value = operation_type::neutral_value();
		}
	};
}	// namespace efficient_DST

#endif // EFFICIENT_DST_ZETA_TRANSFORM_HPP
