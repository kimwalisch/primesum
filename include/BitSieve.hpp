///
/// @file  BitSieve.hpp
/// @brief The BitSieve class is a bit array for prime sieving
///        that packs 64 numbers into 8 bytes i.e. each bit
///        corresponds to one integer.
///
/// Copyright (C) 2016 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#ifndef BITSIEVE_HPP
#define BITSIEVE_HPP

#include <cassert>
#include <cstddef>
#include <vector>
#include <int128.hpp>
#include <stdint.h>

namespace primesum {

class BitSieve
{
public:
  BitSieve(std::size_t size);

  /// Pre-sieve the multiples (>= low) of the first c primes.
  /// @warning Also removes the first c primes.
  /// @pre c < 10
  ///
  void pre_sieve(uint64_t c,
                 uint64_t low);

  /// Count the number of 1 bits inside [start, stop]
  uint64_t count(uint64_t start,
                 uint64_t stop) const;

  /// Count the number of 1 bits inside [0, stop]
  uint64_t count(uint64_t stop) const
  {
    return count(0, stop);
  }

  /// Compute the sum of the unsieved numbers inside [start, stop]
  maxint_t sum(uint64_t low,
               uint64_t start,
               uint64_t stop) const;

  /// Sum the unsieved numbers inside [start, stop].
  /// As an optimization this method counts either forwards or
  /// backwards depending on what's faster.
  ///
  maxint_t sum(uint64_t start,
               uint64_t stop,
               uint64_t low,
               uint64_t high,
               uint64_t sum_0_start,
               uint64_t sum_low_high) const
  {
    if (start > stop)
      return 0;

    if (stop - start < high - low - stop)
      return sum(low, start, stop);
    else
      // optimization, same as sum(low, start, stop)
      return sum_low_high - sum_0_start - sum(low, stop + 1, (high - 1) - low);
  }

  void set(uint64_t pos)
  {
    assert(pos < size_);
    sieve_[pos >> 6] |= ((uint64_t) 1) << (pos & 63);
  }

  void unset(uint64_t pos)
  {
    assert(pos < size_);
    sieve_[pos >> 6] &= unset_bit_[pos & 63];
  }

  bool operator[](uint64_t pos) const
  {
    assert(pos < size_);
    return (sieve_[pos >> 6] >> (pos & 63)) & 1;
  }

  std::size_t size() const
  {
    return size_;
  }
private:
  static const uint64_t unset_bit_[64];
  std::vector<uint64_t> sieve_;
  std::size_t size_;
};

} // namespace

#endif
