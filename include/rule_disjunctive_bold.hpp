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
 
#ifndef EFFICIENT_DST_RULE_DISJUNCTIVE_BOLD_HPP
#define EFFICIENT_DST_RULE_DISJUNCTIVE_BOLD_HPP

#include <implicability.hpp>
#include <disjunctive_weight.hpp>
#include <mass.hpp>

namespace efficient_DST{

	template <typename T, size_t N>
	class rule_disjunctive_bold {
	public:

		std::string to_string() const {
			return "Bold disjunctive rule";
		}

		mass<T, N> operator()(const mass<T, N>& m1, const mass<T, N>& m2) const {
			zeta_transform<T, N, down_inclusion<T, N> > b1(m1.get_definition(), operation_type_t::addition);
			zeta_transform<T, N, down_inclusion<T, N> > b2(m2.get_definition(), operation_type_t::addition);
			implicability<T, N> b12 = operator()(*(implicability<T, N>*) &b1, *(implicability<T, N>*) &b2);
			return mass<T, N>(b12);
		}


		implicability<T, N> operator()(const implicability<T, N>& b1, const implicability<T, N>& b2) const {
			const disjunctive_weight<T, N>& v1(b1);
			const disjunctive_weight<T, N>& v2(b2);
			return implicability<T, N>(operator ()(v1, v2));
		}


		disjunctive_weight<T, N> operator()(const disjunctive_weight<T, N>& v1, const disjunctive_weight<T, N>& v2) const {
			return disjunctive_weight<T, N>(weight_fusion(v1.get_definition(), v2.get_definition()));
		}


	protected:

		static T max(const T& val1, const T& val2){
			return std::max(val1, val2);
		}

		powerset_btree<T, N> weight_fusion(const powerset_btree<T, N>& v1_definition, const powerset_btree<T, N>& v2_definition) const {
			powerset_btree<T, N> v12_definition(v1_definition.get_FOD(), v1_definition.get_block_size());
			v12_definition.fill_with_union_of_powersets(v1_definition, v2_definition, max, 1);
			decomposition_weight<T, N>::remove_negligible_values(v12_definition);
			disjunctive_weight<T, N>::compute_emptyset_value_from_definition(v12_definition);
			return v12_definition;
		}
	};

} // namespace efficient_DST

#endif // EFFICIENT_DST_RULE_DISJUNCTIVE_BOLD_HPP
