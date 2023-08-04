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
 
// Copyright (c) 2017-2020
// Maxime Chaveroche

#ifndef BOOST_BFT_N_DISCOUNTING_HPP
#define BOOST_BFT_N_DISCOUNTING_HPP

#include <cmath>

#include <boost/bft/mass.hpp>
#include <boost/bft/commonality.hpp>
#include <boost/foreach.hpp>

namespace boost
{
namespace bft
{

/// Performs \alpha-discounting n times, i.e. non-Omega masses are discounted by the
/// factor (1-\alpha)^n.
///
/// \alpha-discounting on masses gives:
/// {}^{n,\alpha} m(A) = (1-\alpha)^n m(A)         , \forall A \subset Omega
/// {}^{n,\alpha} m(\Omega) = 1 - (1 - \alpha)^n [1 - {}^{0,\alpha} m(\Omega)]

/// \alpha-discounting on commonalities gives:
/// {}^{n,\alpha} m(A) = 1 - (1 - \alpha)^n [1 - {}^{0,\alpha} m(A)] , \forall A \subseteq Omega
struct n_discounting
{
    typedef double value_type;

    n_discounting(value_type alpha)
        : alpha(alpha)
    {
        BOOST_ASSERT(this->alpha >= 0);
        BOOST_ASSERT(this->alpha <= 1);
    }

    // For masses
    template <class FOD, typename T>
    void operator()(int n, mass<FOD, T>& m) const
    {
    	value_type n_alpha = pow(1-this->alpha, n);
        BOOST_FOREACH (T& v, m.values()) {
            v *= n_alpha;
        }
        m.values().back() += 1-n_alpha;
    }

    template <class FOD, typename T>
    mass<FOD, T> operator()(int n, const mass<FOD, T>& m) const
    {
        mass<FOD, T> m_out(m);
        operator()(n, m_out);
        return m_out;
    }

    // For commonalities
    template <class FOD, typename T>
    void operator()(int n, commonality<FOD, T>& q) const
    {
    	value_type n_alpha = pow(1-this->alpha, n);
        BOOST_FOREACH (T& v, q.values()) {
            v = 1-n_alpha*(1-v);
        }
    }

    template <class FOD, typename T>
    commonality<FOD, T> operator()(int n, const commonality<FOD, T>& q) const
    {
    	commonality<FOD, T> q_out(q);
        operator()(n, q_out);
        return q_out;
    }


private:
    value_type alpha;
};

} // namespace bft

} // namespace boost

#endif // BOOST_BFT_N_DISCOUNTING_HPP
