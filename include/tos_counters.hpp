///
/// @file   tos_counters.hpp
/// @brief  This file contains functions to initialize, update and query
///         the special tree data structure used for summing the
///         number of unsieved elements in the Lagarias-Miller-Odlyzko
///         and Deleglise-Rivat prime summing algorithms.
///
///         The implementation is a modified version of the
///         algorithm described in the paper:
///
///         Tomás Oliveira e Silva, Computing pi(x): the combinatorial method,
///         Revista do DETUA, vol. 4, no. 6, March 2006, pp. 767-768.
///         http://sweet.ua.pt/tos/bib/5.4.pdf
///
/// Copyright (C) 2016 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#ifndef TOS_COUNTERS_HPP
#define TOS_COUNTERS_HPP

#include <stdint.h>
#include <vector>

namespace primesum {

/// Initialize the counters from the sieve array.
/// @pre segment_size is a power of 2.
/// @pre sieve[i] = 1 for unsieved elements and sieve[i] = 0
///      for crossed-off elements.
/// Runtime: O(N log N).
///
template <typename T1, typename T2>
inline void cnt_finit(const T1& sieve,
                      std::vector<T2>& cnt,
                      int64_t low,
                      int64_t segment_size)
{
  for (T2 i = 0; i < segment_size; i++)
  {
    cnt[i] = (low + i) * sieve[i];
    for (T2 k = (i + 1) & ~i, j = i; k >>= 1; j &= j - 1)
      cnt[i] += cnt[j - 1];
  }
}

/// Update the counters after that an element has been
/// crossed-off for the first time in the sieve array.
/// @pre segment_size is a power of 2.
/// Runtime: O(log N).
///
template <typename T>
inline void cnt_update(std::vector<T>& cnt,
                       int64_t n,
                       int64_t low,
                       int64_t segment_size)
{
  int64_t i = n - low;
  do
  {
    cnt[i] -= n;
    i |= i + 1;
  }
  while (i < segment_size);
}

/// Get the sum of the unsieved elements <= pos
/// in the current segment (sieve array).
/// Runtime: O(log N).
///
template <typename T>
inline int64_t cnt_query(const std::vector<T>& cnt, int64_t pos)
{
  int64_t sum = cnt[pos++];
  for (; pos &= pos - 1; sum += cnt[pos - 1]);
  return sum;
}

} // namespace

#endif
