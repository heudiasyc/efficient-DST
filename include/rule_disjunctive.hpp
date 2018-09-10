// Copyright (c) 2011-2014
// Marek Kurdej
//
// Distributed under the Boost Software License, Version 1.0.
// See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt

#ifndef BOOST_BFT_RULE_DISJUNCTIVE_HPP
#define BOOST_BFT_RULE_DISJUNCTIVE_HPP

#include <boost/bft/implicability.hpp>
#include <boost/bft/mass.hpp>
#include <boost/bft/rule_base.hpp>
#include <boost/bft/to_implicability.hpp>
#include <boost/bft/to_mass.hpp>

namespace boost
{
namespace bft
{

struct rule_disjunctive : public rule_base
{
    std::string to_string() const
    {
        return "disjunctive rule";
    }

    template <class FOD, typename T>
    implicability<FOD, T> operator()(const implicability<FOD, T>& b1,
                                     const implicability<FOD, T>& b2) const
    {
        implicability<FOD, T> b12;
        for (std::size_t i = 0; i < FOD::powerset_size; ++i) {
            b12[i] = b1[i] * b2[i];
        }
        return b12;
    }

    template <class FOD, typename T>
    mass<FOD, T>
    operator()(const mass<FOD, T>& m1, const mass<FOD, T>& m2) const
    {
        implicability<FOD, T> b1 = to_implicability(m1);
        implicability<FOD, T> b2 = to_implicability(m2);
        implicability<FOD, T> b12 = operator()(b1, b2);
        return to_mass(b12);
    }
};

} // namespace bft

} // namespace boost

#endif // BOOST_BFT_RULE_CONJUNCTIVE_HPP
