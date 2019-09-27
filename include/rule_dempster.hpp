#ifndef EFFICIENT_DST_RULE_DEMPSTER_HPP
#define EFFICIENT_DST_RULE_DEMPSTER_HPP

#include <commonality.hpp>
#include <conjunctive_weight.hpp>
#include <mass.hpp>
#include <rule_conjunctive.hpp>

namespace efficient_DST{

	template <typename T = double>
	class rule_Dempster : public rule_conjunctive<T> {
	public:

		std::string to_string() const {
			return "Dempster's rule";
		}

	public:

		mass<T> operator()(const mass<T>& m1, const mass<T>& m2) const {
			mass<T> m12 = rule_conjunctive<T>::operator ()(m1, m2);
			m12.regularize();
			return m12;
		}

		commonality<T> operator()(const commonality<T>& q1, const commonality<T>& q2) const {
			const conjunctive_weight<T>& w1(q1);
			const conjunctive_weight<T>& w2(q2);
			return commonality<T>(this->operator()(w1, w2));
		}

		conjunctive_weight<T> operator()(const conjunctive_weight<T>& w1, const conjunctive_weight<T>& w2) const {
			conjunctive_weight<T> w12 = rule_conjunctive<T>::operator ()(w1, w2);
			w12.regularize();
			return w12;
		}
	};

} // namespace efficient_DST

#endif // EFFICIENT_DST_RULE_DEMPSTER_HPP
