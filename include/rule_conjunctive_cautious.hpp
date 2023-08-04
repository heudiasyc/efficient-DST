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
 
#ifndef EFFICIENT_DST_RULE_CONJUNCTIVE_CAUTIOUS_HPP
#define EFFICIENT_DST_RULE_CONJUNCTIVE_CAUTIOUS_HPP

#include <commonality.hpp>
#include <conjunctive_weight.hpp>
#include <mass.hpp>

namespace efficient_DST{

	template <typename T, size_t N>
	class rule_conjunctive_cautious {
	public:

		std::string to_string() const {
			return "Cautious conjunctive rule";
		}

		mass<T, N> operator()(const mass<T, N>& m1, const mass<T, N>& m2) const {
			zeta_transform<T, N, up_inclusion<T, N> > q1(m1.get_definition(), operation_type_t::addition);
			zeta_transform<T, N, up_inclusion<T, N> > q2(m2.get_definition(), operation_type_t::addition);
			commonality<T, N> q12 = operator()(*(commonality<T, N>*) &q1, *(commonality<T, N>*) &q2);
			return mass<T, N>(q12);
		}


		commonality<T, N> operator()(const commonality<T, N>& q1, const commonality<T, N>& q2) const {
			const conjunctive_weight<T, N>& w1(q1);
			const conjunctive_weight<T, N>& w2(q2);
			return commonality<T, N>(operator ()(w1, w2));
		}


		conjunctive_weight<T, N> operator()(const conjunctive_weight<T, N>& w1, const conjunctive_weight<T, N>& w2) const {
			return conjunctive_weight<T, N>(weight_fusion(w1.get_definition(), w2.get_definition()));
		}


	protected:

		static T min(const T& val1, const T& val2){
			return std::min(val1, val2);
		}

		powerset_btree<T, N> weight_fusion(const powerset_btree<T, N>& w1_definition, const powerset_btree<T, N>& w2_definition) const {
			powerset_btree<T, N> w12_definition(w1_definition.get_FOD(), w1_definition.get_block_size());
			w12_definition.fill_with_union_of_powersets(w1_definition, w2_definition, min, 1);
			decomposition_weight<T, N>::remove_negligible_values(w12_definition);
			conjunctive_weight<T, N>::compute_fod_value_from_definition(w12_definition);
			return w12_definition;
		}
	};

} // namespace efficient_DST

#endif // EFFICIENT_DST_RULE_CONJUNCTIVE_CAUTIOUS_HPP
