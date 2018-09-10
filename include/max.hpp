// Copyright (c) 2011-2014
// Marek Kurdej
//
// Distributed under the Boost Software License, Version 1.0.
// See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt

#ifndef BOOST_BFT_MAX_HPP
#define BOOST_BFT_MAX_HPP

#include <boost/bft/bft_function.hpp>
#include <boost/math/tr1.hpp> // BOOST_PREVENT_MACRO_SUBSTITUTION

namespace boost
{
namespace bft
{

template <class FOD, typename T>
bft_function<FOD, T>(max)
    BOOST_PREVENT_MACRO_SUBSTITUTION(const bft_function<FOD, T>& f1,
                                     const bft_function<FOD, T>& f2)
{
    bft_function<FOD, T> f_result;
    for (int A = 0; A < FOD::powerset_size; ++A) {
        f_result[A] = std::max(f1[A], f2[A]);
    }
    return f_result;
}

} // namespace bft

} // namespace boost

#endif // BOOST_BFT_MAX_HPP
