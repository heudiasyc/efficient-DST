// Copyright (c) 2017-2020
// Maxime Chaveroche

#ifndef BOOST_BFT_RULE_CONJUNCTIVE_CAUTIOUS_DISCOUNTED_HPP
#define BOOST_BFT_RULE_CONJUNCTIVE_CAUTIOUS_DISCOUNTED_HPP

#include <boost/bft/conjunctive_decomposition.hpp>
#include <boost/bft/mass.hpp>
#include <boost/bft/min.hpp>
#include <boost/bft/rule_base.hpp>
#include <boost/bft/rule_conjunctive.hpp>
#include <boost/bft/detail/is_small.hpp>

namespace boost
{
namespace bft
{

struct rule_conjunctive_cautious_discounted : public rule_base
{
    std::string to_string() const
    {
        return "discounted cautious conjunctive rule";
    }

    const double alpha = 0;
    long n = 0;
    mass<FOD, T>& m01 = 0;
    mass<FOD, T>& m02 = 0;

    template <class FOD, typename T>
    mass<FOD, T>
    operator()(const mass<FOD, T>& m1, const mass<FOD, T>& m2) const
    {
        // compute canonical conjunctive decompositions w_1(A), w_2(A)
        conjunctive_decomposition decomposition;
        const bft_function<FOD, T> w1 = m1.apply(decomposition);
        const bft_function<FOD, T> w2 = m2.apply(decomposition);
        const double eta1;
        const double eta2;

        bft_function<FOD, T> w_min12 = min(w1, w2);

        // compute conjunctive combination of GSBBA A^{w_1(A) \wedge w_2(A)}
        mass<FOD, T> m_result(vacuous);
        rule_conjunctive rule;
        for (std::size_t A = 0; A < FOD::powerset_size - 1; ++A) {
            if (detail::is_small(1 - w_min12[A], detail::tolerance)) {
                continue;
            }
            m_result = m_result.apply(
                rule, mass<FOD, T>::create_mass_from_conjunctive_weight(
                          A, w_min12[A]));
        }
        return m_result;
    }
};

} // namespace bft

} // namespace boost

#endif // BOOST_BFT_RULE_CONJUNCTIVE_CAUTIOUS_DISCOUNTED_HPP
