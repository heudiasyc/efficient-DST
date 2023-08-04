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
 
#ifndef EFFICIENT_DST_RULE_DISJUNCTIVE_HPP
#define EFFICIENT_DST_RULE_DISJUNCTIVE_HPP

#include <rule_classic_general.hpp>
#include <implicability.hpp>
#include <disjunctive_weight.hpp>
#include <mass.hpp>

namespace efficient_DST{

	template <typename T, size_t N>
	class rule_disjunctive : public rule_classic_general<T, N> {
	public:

		std::string to_string() const {
			return "Disjunctive rule";
		}

		mass<T, N> operator()(const mass<T, N>& m1, const mass<T, N>& m2) const {
			const powerset_btree<T, N>& m1_definition = m1.get_definition();
			const powerset_btree<T, N>& m2_definition = m2.get_definition();

			if (m1_definition.size() * m2_definition.size() < 0.5 * m1_definition.get_FOD_size() * pow(2, m1_definition.get_FOD_size())){
				return rule_classic_general<T, N>::operator ()(m1, m2, FOD<N>::set_union);
			}else{
				zeta_transform<T, N, down_inclusion<T, N> > b1(m1_definition, operation_type_t::addition);
				zeta_transform<T, N, down_inclusion<T, N> > b2(m2_definition, operation_type_t::addition);
				implicability<T, N> b12 = operator()(*(implicability<T, N>*) &b1, *(implicability<T, N>*) &b2);
				return mass<T, N>(b12);
			}
		}


		implicability<T, N> operator()(const implicability<T, N>& b1, const implicability<T, N>& b2) const {
			const powerset_btree<T, N>& v1_inverted_definition = b1.inversion(operation_type_t::multiplication);
			const powerset_btree<T, N>& v2_inverted_definition = b2.inversion(operation_type_t::multiplication);
			zeta_transform<T, N, down_inclusion<T, N>> b12(
				rule_classic_general<T, N>::weight_fusion(v1_inverted_definition, v2_inverted_definition),
				operation_type_t::multiplication
			);
			return *(implicability<T, N>*) &b12;
		}


		disjunctive_weight<T, N> operator()(const disjunctive_weight<T, N>& v1, const disjunctive_weight<T, N>& v2) const {
			return disjunctive_weight<T, N>(rule_classic_general<T, N>::weight_fusion(v1.get_definition(), v2.get_definition()));
		}
	};

} // namespace efficient_DST

#endif // EFFICIENT_DST_RULE_DISJUNCTIVE_HPP
