///
/// @file  BitSieve.hpp
/// @brief The BitSieve class is a bit array for prime sieving
///        that packs 128 numbers into 8 bytes i.e. each bit
///        corresponds to an odd integer.
///
/// Copyright (C) 2017 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#ifndef BITSIEVE_HPP
#define BITSIEVE_HPP

#include <int128_t.hpp>
#include <stdint.h>

#include <cassert>
#include <cstddef>
#include <array>
#include <vector>

namespace primesum {

class BitSieve
{
public:
  BitSieve(std::size_t size);

  /// Pre-sieve the multiples >= low of the first c primes.
  /// @warning Also removes the first c primes.
  /// @pre c < 10
  ///
  void pre_sieve(uint64_t c,
                 uint64_t low);

  /// Compute the sum of the unsieved numbers inside [start, stop]
  int128_t sum(uint64_t low,
               uint64_t start,
               uint64_t stop) const;

  /// Sum the unsieved numbers inside [start, stop].
  /// As an optimization this method counts either forwards or
  /// backwards depending on what's faster.
  ///
  int128_t sum(uint64_t start,
               uint64_t stop,
               uint64_t low,
               uint64_t high,
               int128_t sum_0_start,
               int128_t sum_low_high) const
  {
    if (start > stop)
      return 0;

    if (stop - start < (high - low) - stop)
      return sum(low, start, stop);
    else
      return sum_low_high - sum_0_start - sum(low, stop + 1, (high - 1) - low);
  }

  void set(uint64_t pos)
  {
    assert(pos < size_);
    assert(pos % 2 != 0);
    sieve_[pos >> 7] |= set_bit_[pos & 127];
  }

  void unset(uint64_t pos)
  {
    assert(pos < size_);
    assert(pos % 2 != 0);
    sieve_[pos >> 7] &= ~set_bit_[pos & 127];
  }

  bool operator[](uint64_t pos) const
  {
    assert(pos < size_);
    assert(pos % 2 != 0);
    return (sieve_[pos >> 7] & set_bit_[pos & 127]) != 0;
  }

  std::size_t size() const
  {
    return size_;
  }

private:
  static const std::array<uint64_t, 128> set_bit_;
  std::vector<uint64_t> sieve_;
  std::size_t size_;
};

} // namespace

#endif
