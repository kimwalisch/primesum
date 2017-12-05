///
/// @file  SumBits.hpp
/// @brief Class that returns the sum of the 1 bit indexes
///        inside of a 16-bit word.
///
/// Copyright (C) 2017 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#ifndef SUMBITS_HPP
#define SUMBITS_HPP

#include <stdint.h>

namespace primesum {

class SumBits
{
public:
  SumBits()
  {
    for (int i = 0; i < (1 << 16); i++)
    {
      sum_bits_[i] = 0;
      for (int j = i; j != 0; j &= j - 1)
        sum_bits_[i] += (uint8_t) (__builtin_ctzll(j) * 2);
    }
  }

  template <typename T>
  uint8_t operator[](T i) const
  {
    return sum_bits_[i];
  }

private:
  uint8_t sum_bits_[1 << 16];
};

} // namespace

#endif
