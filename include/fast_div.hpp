///
/// @file  fast_div.hpp
/// @brief Integer division of small types is much faster than
///        integer division of large types on most CPUs. The
///        fast_div(x, y) function tries to take advantage of
///        this by casting x and y to smaller types (if possible)
///        before doing the division.
///
/// Copyright (C) 2017 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#ifndef FAST_DIV_HPP
#define FAST_DIV_HPP

#include <int128_t.hpp>

#include <cassert>
#include <limits>
#include <type_traits>

namespace primesum {

template <typename T>
struct fastdiv
{
  typedef typename std::conditional<sizeof(T) / 2 <= sizeof(uint8_t), uint8_t,
          typename std::conditional<sizeof(T) / 2 == sizeof(uint16_t), uint16_t,
          typename std::conditional<sizeof(T) / 2 == sizeof(uint32_t), uint32_t,
          uint64_t>::type>::type>::type type;
};

template <typename X, typename Y>
typename std::enable_if<(sizeof(X) == sizeof(Y)), X>::type
fast_div(X x, Y y)
{
  static_assert(prt::is_integral<X>::value &&
                prt::is_integral<Y>::value,
                "fast_div(x, y): types must be integral");

  using fastdiv_t = typename fastdiv<X>::type;

  if (x <= std::numeric_limits<fastdiv_t>::max() &&
      y <= std::numeric_limits<fastdiv_t>::max())
  {
    return (fastdiv_t) x / (fastdiv_t) y;
  }

  return x / y;
}

template <typename X, typename Y>
typename std::enable_if<(sizeof(X) > sizeof(Y)), X>::type
fast_div(X x, Y y)
{
  static_assert(prt::is_integral<X>::value &&
                prt::is_integral<Y>::value,
                "fast_div(x, y): types must be integral");

  using fastdiv_t = typename fastdiv<X>::type;

  if (x <= std::numeric_limits<fastdiv_t>::max())
    return (fastdiv_t) x / (fastdiv_t) y;

  return x / y;
}

} // namespace

#endif
