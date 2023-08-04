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

#ifndef EFFICIENT_DST_ZETA_TRANSFORM_VECTOR_HPP
#define EFFICIENT_DST_ZETA_TRANSFORM_VECTOR_HPP

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include <powerset_vector.hpp>
#include <mobius_inversion.hpp>
#include <zeta_transform.hpp>


namespace efficient_DST{

	template <class inclusion, size_t N, typename T = float>
	class zeta_transform_vector : public powerset_vector<N, T> {
	protected:
		size_t zero = 0;
		size_t core = (1 << N)-1;
		bool core_reduced;

	public:

		zeta_transform_vector(const zeta_transform_vector<inclusion, N, T>& z) :
			powerset_vector<N, T>(z.outcomes, z.definition, z.default_value),
			core(z.core),
			core_reduced(z.core_reduced)
		{}

		/*
		 * Constructor when you have a Möbius transform such as the mass or conjunctive/disjunctive weight function.
		 * - support is supposed to contain all focal sets and their image.
		 * - order_relation is the order relation of this zeta transform (e.g. commonality->superset or implicability->subset).
		 * - transform_operation is the operation of the zeta transform (e.g. in DST, we usually use the addition on the mass function
		 * and the multiplication on the conjunctive/disjunctive weight function).
		 */
		zeta_transform_vector(
			const sample_space<N>& outcomes,
			const std::vector<T>& support,
			const T& default_value,
			operation_type_t operation_type,
			const bool& core_reduced = true
		) :
			powerset_vector<N, T>(outcomes, default_value),
			core_reduced(core_reduced)
		{
			if (operation_type == operation_type_t::addition){
				this->compute<addition<T>>(support);
			} else {
				this->compute<multiplication<T>>(support);
			}
		}

//		template<class operation_type>
//		powerset_btree<N, T> inversion() const {
		std::vector<T> inversion(const operation_type_t& operation_type, const bool& core_reduced = false) const {
			std::vector<T> mobius_transform_definition;
			if (core_reduced && !this->core_reduced){
				size_t zero = augmented_zero(this->definition);
				size_t core = reduced_core(this->definition);
				size_t n = std::bitset<N>(core ^ zero).count();
				mobius_transform_definition = std::vector<T>(1 << n, this->default_value);
				size_t j = 0;
				for (size_t i = 0; i < (1 << N); ++i){
					if ((i | zero) == i && (i & core) == i){
						mobius_transform_definition[j] = this->definition[i];
						++j;
					}
				}
			}else{
				mobius_transform_definition = this->definition;
			}
			if (operation_type == operation_type_t::addition){
				efficient_mobius_inversion<inclusion, mobius_tranformation<inclusion, addition<T>, N, T>, N, T >::execute_FMT(
						mobius_transform_definition
				);
				if (this->definition.size() < (1 << N))
					mobius_transform_definition = expansion_from_core(mobius_transform_definition, addition<T>::neutral_value());
			} else {
				efficient_mobius_inversion<inclusion, mobius_tranformation<inclusion, multiplication<T>, N, T>, N, T >::execute_FMT(
						mobius_transform_definition
				);
				if (this->definition.size() < (1 << N))
					mobius_transform_definition = expansion_from_core(mobius_transform_definition, multiplication<T>::neutral_value());
			}
			return mobius_transform_definition;
		}


		std::ostream& print(const bool& including_null = false) const {
			std::cout << std::endl;
			size_t j = 0;
			for (size_t i = 0; i < (1 << N); ++i) {
				if ((i | this->zero) == i && (i & this->core) == i){
//					if (including_null || !powerset_function<N, T>::is_equivalent_to_zero(this->definition[j] - this->default_value)){
					std::cout << set_N_value<N, T>::to_string(this->definition[j]) + "\t <- " + this->outcomes.to_string(i) << std::endl;
//					}
					++j;
				}
			}
			return std::cout;
		}

		size_t reduced_core(const std::vector<T>& support) const {
			size_t core = 0;
			for (size_t i = 0; i < support.size(); ++i){
				if (!powerset_function<N, T>::is_equivalent_to_zero(support[i] - this->default_value)){
					core |= i;
				}
			}
			return core;
		}

		size_t augmented_zero(const std::vector<T>& support) const {
			size_t zero = (1 << N)-1;
			for (size_t i = 0; i < support.size(); ++i){
				if (!powerset_function<N, T>::is_equivalent_to_zero(support[i] - this->default_value)){
					zero &= i;
				}
			}
			return zero;
		}

		static std::vector<T> expansion_from_core(const std::vector<T>& core_definition, const T& default_value, const size_t& core, const size_t& zero) {
			if (core_definition.size() == (1 << N))
				return core_definition;
			std::vector<T> expansion(1 << N, default_value);
			size_t j = 0;
			for (size_t i = 0; i < (1 << N); ++i){
				if ((i | zero) == i && (i & core) == i){
					expansion[i] = core_definition[j];
					++j;
				}
			}
			return expansion;
		}

	protected:

		void core_reduction(const std::vector<T>& support){
			this->zero = augmented_zero(support);
			this->core = reduced_core(support);
//			std::cout << "Core = " << this->core << "\n";
			size_t n = std::bitset<N>(this->core ^ this->zero).count();
			this->definition = std::vector<T>(1 << n, this->default_value);
			size_t j = 0;
			for (size_t i = 0; i < (1 << N); ++i){
				if ((i | this->zero) == i && (i & this->core) == i){
					this->definition[j] = support[i];
					++j;
				}
			}
		}

		std::vector<T> expansion_from_core(const std::vector<T>& core_definition, const T& default_value) const {
			return expansion_from_core(core_definition, default_value, this->core, this->zero);
		}

		template<class operation_type>
		void compute(const std::vector<T>& support){
			if (this->default_value != operation_type::neutral_value()){
				std::cerr << "Ill-defined support: the default value of your compact definition must match the neutral value for the operator you chose.\n";
				exit(1);
			}
			if (this->core_reduced){
				core_reduction(support);
			}else{
				this->definition = support;
//				this->core = (1 << N)-1;
			}
			efficient_mobius_inversion<
				inclusion, zeta_tranformation<inclusion, operation_type, N, T>, N, T
			>::execute_FMT(
				this->definition
			);
		}

		void normalize(const size_t& index, const T& normalizing_value){
			if(this->definition[index] != normalizing_value){
				if(this->definition[index] != this->default_value){
					for (size_t i = 0; i < this->definition.size(); ++i){
						if(i != index){
							this->definition[i] /= this->definition[index];
						}
					}
				}
				for (size_t i = 0; i < this->definition.size(); ++i){
					if(i != index){
						this->definition[i] *= normalizing_value;
					}
				}
				this->definition[index] = normalizing_value;
			}
		}
	};
}	// namespace efficient_DST

#endif // EFFICIENT_DST_ZETA_TRANSFORM_VECTOR_HPP
