// Copyright (c) 2011-2014
// Marek Kurdej
//
// Distributed under the Boost Software License, Version 1.0.
// See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt

#ifndef BOOST_BFT_DISCOUNTING_HPP
#define BOOST_BFT_DISCOUNTING_HPP

#include <boost/bft/mass.hpp>
#include <boost/foreach.hpp>

namespace boost
{
namespace bft
{

/// Performs \alpha-discounting, i.e. non-Omega masses are discounted by the
/// factor 1-\alpha.
///
/// \alpha-discounting gives:
/// {}^\alpha m(A) = (1-\alpha) m(A)         , \forall A \subset Omega
/// {}^\alpha m(\Omega) = (1-\alpha) m(\Omega) + \alpha
struct discounting
{
    typedef double value_type;

    discounting(value_type alpha)
        : m_alpha(alpha)
    {
        BOOST_ASSERT(this->m_alpha >= 0);
        BOOST_ASSERT(this->m_alpha <= 1);
    }

    template <class FOD, typename T>
    void operator()(mass<FOD, T>& m) const
    {
        BOOST_FOREACH (T& v, m.values()) {
            v *= (1 - this->m_alpha);
        }
        m.values().back() += this->m_alpha;
    }

    template <class FOD, typename T>
    mass<FOD, T> operator()(const mass<FOD, T>& m) const
    {
        mass<FOD, T> m_out(m);
        operator()(m_out);
        return m_out;
    }

private:
    value_type m_alpha;
};

} // namespace bft

} // namespace boost

#endif // BOOST_BFT_DISCOUNTING_HPP
