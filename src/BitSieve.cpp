///
/// @file  BitSieve.cpp
/// @brief The BitSieve class is a bit array for prime sieving
///        that packs 128 numbers into 8 bytes i.e. each bit
///        corresponds to an odd integer.
///
/// Copyright (C) 2017 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <BitSieve.hpp>
#include <popcnt.hpp>
#include <imath.hpp>
#include <int128_t.hpp>
#include <SumBits.hpp>

#include <stdint.h>
#include <algorithm>
#include <cassert>
#include <vector>

using namespace std;
using namespace primesum;

namespace {

const SumBits sumBits;

const uint64_t primes[] = { 0, 2, 3, 5, 7, 11, 13, 17, 19, 23 };

/// bitmasks with multiples of the i-th prime
const uint64_t masks[] =
{
  0x0000000000000000ull,
  0x5555555555555555ull, // 2
  0x9249249249249249ull, // 3
  0x1084210842108421ull, // 5
  0x8102040810204081ull, // 7
  0x0080100200400801ull, // 11
  0x0010008004002001ull, // 13
  0x0008000400020001ull, // 17
  0x0200004000080001ull, // 19
  0x0000400000800001ull  // 23
};

uint64_t fast_modulo(uint64_t x, uint64_t y)
{
  assert(x < y * 2);
  x = (x < y) ? x : x - y;
  return x;
}

uint64_t sum_bits(uint64_t bits, uint64_t& low)
{
  uint64_t sum = 0;
  uint64_t bits0 = bits & 0xffff;
  uint64_t bits1 = (bits >> 16) & 0xffff;
  uint64_t bits2 = (bits >> 32) & 0xffff;
  uint64_t bits3 = (bits >> 48);

  sum += low * popcnt64(bits0) + sumBits[bits0];
  sum += (low + 32) * popcnt64(bits1) + sumBits[bits1];
  sum += (low + 64) * popcnt64(bits2) + sumBits[bits2];
  sum += (low + 96) * popcnt64(bits3) + sumBits[bits3];

  low += 128;

  return sum;
}

int128_t sum_bits(const uint64_t* bits, uint64_t size, uint64_t& low)
{
  int128_t sum = 0;

  for (uint64_t i = 0; i < size; i++)
    sum += sum_bits(bits[i], low);

  return sum;
}

} // namespace

namespace primesum {

const uint64_t BitSieve::set_bit_[128] =
{
  (1ull <<  0), (1ull <<  0), (1ull <<  1), (1ull <<  1),
  (1ull <<  2), (1ull <<  2), (1ull <<  3), (1ull <<  3),
  (1ull <<  4), (1ull <<  4), (1ull <<  5), (1ull <<  5),
  (1ull <<  6), (1ull <<  6), (1ull <<  7), (1ull <<  7),
  (1ull <<  8), (1ull <<  8), (1ull <<  9), (1ull <<  9),
  (1ull << 10), (1ull << 10), (1ull << 11), (1ull << 11),
  (1ull << 12), (1ull << 12), (1ull << 13), (1ull << 13),
  (1ull << 14), (1ull << 14), (1ull << 15), (1ull << 15),
  (1ull << 16), (1ull << 16), (1ull << 17), (1ull << 17),
  (1ull << 18), (1ull << 18), (1ull << 19), (1ull << 19),
  (1ull << 20), (1ull << 20), (1ull << 21), (1ull << 21),
  (1ull << 22), (1ull << 22), (1ull << 23), (1ull << 23),
  (1ull << 24), (1ull << 24), (1ull << 25), (1ull << 25),
  (1ull << 26), (1ull << 26), (1ull << 27), (1ull << 27),
  (1ull << 28), (1ull << 28), (1ull << 29), (1ull << 29),
  (1ull << 30), (1ull << 30), (1ull << 31), (1ull << 31),
  (1ull << 32), (1ull << 32), (1ull << 33), (1ull << 33),
  (1ull << 34), (1ull << 34), (1ull << 35), (1ull << 35),
  (1ull << 36), (1ull << 36), (1ull << 37), (1ull << 37),
  (1ull << 38), (1ull << 38), (1ull << 39), (1ull << 39),
  (1ull << 40), (1ull << 40), (1ull << 41), (1ull << 41),
  (1ull << 42), (1ull << 42), (1ull << 43), (1ull << 43),
  (1ull << 44), (1ull << 44), (1ull << 45), (1ull << 45),
  (1ull << 46), (1ull << 46), (1ull << 47), (1ull << 47),
  (1ull << 48), (1ull << 48), (1ull << 49), (1ull << 49),
  (1ull << 50), (1ull << 50), (1ull << 51), (1ull << 51),
  (1ull << 52), (1ull << 52), (1ull << 53), (1ull << 53),
  (1ull << 54), (1ull << 54), (1ull << 55), (1ull << 55),
  (1ull << 56), (1ull << 56), (1ull << 57), (1ull << 57),
  (1ull << 58), (1ull << 58), (1ull << 59), (1ull << 59),
  (1ull << 60), (1ull << 60), (1ull << 61), (1ull << 61),
  (1ull << 62), (1ull << 62), (1ull << 63), (1ull << 63)
};

BitSieve::BitSieve(std::size_t size) :
  sieve_(ceil_div(size, 128)),
  size_(size)
{ }

/// Pre-sieve the multiples >= low of the first c primes.
/// @warning Also removes the first c primes.
/// @pre c < 10
///
void BitSieve::pre_sieve(uint64_t c, uint64_t low)
{
  assert(c < primes.size());

  if (sieve_.empty())
    return;

  // reset all bits
  sieve_[0] = ~0ull;

  uint64_t last = 1;
  uint64_t sieved = last;
  uint64_t sieve_size = sieve_.size();

  // pre-sieve multiples of first c primes
  for (uint64_t i = 2; i <= c; i++)
  {
    uint64_t prime = primes[i];
    uint64_t end_copy = min(sieved * prime, sieve_size);

    // pre-sieve multiples of primes < i-th prime
    // by copying a small pre-sieved buffer
    while (last < end_copy)
    {
      uint64_t copy_words = min(sieved, sieve_size - last);
      copy(sieve_.begin(),
           sieve_.begin() + copy_words,
           sieve_.begin() + last);
      last += copy_words;
    }

    // first multiple (of prime) >= low
    uint64_t multiple = ceil_div(low, prime) * prime;
    multiple += prime * (~multiple & 1);
    uint64_t p64 = prime - 64 % prime;
    uint64_t mask = masks[i];
    uint64_t shift = (multiple - low) / 2;

    // pre-sieve multiples of the i-th prime
    for (uint64_t j = 0; j < last; j++)
    {
      sieve_[j] &= ~(mask << shift);
      shift = fast_modulo(shift + p64, prime);
    }

    sieved = last;
  }

  // fill up the rest of the sieve
  while (last < sieve_size)
  {
    uint64_t copy_words = min(sieved, sieve_size - last);
    copy(sieve_.begin(),
         sieve_.begin() + copy_words,
         sieve_.begin() + last);
    last += copy_words;
  }
}

/// Compute the sum of the unsieved numbers inside [start, stop]
int128_t BitSieve::sum(uint64_t low,
                       uint64_t start,
                       uint64_t stop) const
{
  // start & stop must be odd
  start += start % 2;
  stop -= stop % 2;

  if (start > stop)
    return 0;

  assert(stop < size_);

  uint64_t start_idx = start / 128;
  uint64_t stop_idx = stop / 128;
  uint64_t m1 = ~(set_bit_[start % 128] - 1);
  uint64_t m2 = ((set_bit_[stop % 128] - 1) << 1) | 1;

  int128_t sum;
  low += start - start % 128;

  if (start_idx == stop_idx)
    sum = sum_bits(sieve_[start_idx] & (m1 & m2), low);
  else
  {
    sum = sum_bits(sieve_[start_idx] & m1, low);
    sum += sum_bits(&sieve_[start_idx + 1], stop_idx - (start_idx + 1), low);
    sum += sum_bits(sieve_[stop_idx] & m2, low);
  }

  return sum;
}

} // namespace
