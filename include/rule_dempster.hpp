#ifndef EFFICIENT_DST_RULE_DEMPSTER_HPP
#define EFFICIENT_DST_RULE_DEMPSTER_HPP

#include <commonality.hpp>
#include <conjunctive_weight.hpp>
#include <mass.hpp>
#include <rule_conjunctive.hpp>

namespace efficient_DST{

	template <typename T, size_t N>
	class rule_Dempster : public rule_conjunctive<T, N> {
	public:

		std::string to_string() const {
			return "Dempster's rule";
		}

	public:

		mass<T, N> operator()(const mass<T, N>& m1, const mass<T, N>& m2) const {
			mass<T, N> m12 = rule_conjunctive<T, N>::operator ()(m1, m2);
			m12.regularize();
			return m12;
		}

		commonality<T, N> operator()(const commonality<T, N>& q1, const commonality<T, N>& q2) const {
			const conjunctive_weight<T, N>& w1(q1);
			const conjunctive_weight<T, N>& w2(q2);
			return commonality<T, N>(this->operator()(w1, w2));
		}

		conjunctive_weight<T, N> operator()(const conjunctive_weight<T, N>& w1, const conjunctive_weight<T, N>& w2) const {
			conjunctive_weight<T, N> w12 = rule_conjunctive<T, N>::operator ()(w1, w2);
			w12.regularize();
			return w12;
		}
	};

} // namespace efficient_DST

#endif // EFFICIENT_DST_RULE_DEMPSTER_HPP
