///
/// @file  fast_div.cpp
/// @brief Test fast_div(x, y) function
///
/// Copyright (C) 2026 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <fast_div.hpp>
#include <int128_t.hpp>

#include <stdint.h>
#include <cstdlib>
#include <iostream>
#include <random>

using namespace primesum;

void check(bool OK)
{
  std::cout << "   " << (OK ? "OK" : "ERROR") << "\n";
  if (!OK)
    std::exit(1);
}

template <typename X, typename Y>
void test_boundaries()
{
  using UX = typename pstd::make_unsigned<X>::type;
  X max = pstd::numeric_limits<X>::max();
  Y divisors[] = { 1, 2, 3, pstd::numeric_limits<Y>::max() };

  for (Y y : divisors)
  {
    check(fast_div(X(0), y) == 0);
    check(fast_div(max, y) == max / y);

    for (UX x = 1; x <= UX(max); x <<= 1)
    {
      check(fast_div(X(x - 1), y) == X(x - 1) / y);
      check(fast_div(X(x), y) == X(x) / y);
      check(fast_div(X(x + 1), y) == X(x + 1) / y);
      if (x > UX(max) / 2)
        break;
    }

    if (sizeof(Y) < sizeof(X))
    {
      // Test quotient boundaries at 2^32 and 2^64.
      int bits = (sizeof(X) == sizeof(uint64_t)) ? 32 : 64;
      UX boundary = UX(y) << bits;
      if (boundary > 0 && boundary < UX(max))
      {
        check(fast_div(X(boundary - 1), y) == X(boundary - 1) / y);
        check(fast_div(X(boundary), y) == X(boundary) / y);
        check(fast_div(X(boundary + 1), y) == X(boundary + 1) / y);
      }
    }
  }
}

int main()
{
  test_boundaries<int64_t, int32_t>();
  test_boundaries<int64_t, uint32_t>();
  test_boundaries<uint64_t, int32_t>();
  test_boundaries<uint64_t, uint32_t>();
  test_boundaries<int64_t, int64_t>();
  test_boundaries<uint64_t, uint64_t>();
  test_boundaries<int128_t, int32_t>();
  test_boundaries<int128_t, uint32_t>();
  test_boundaries<uint128_t, int32_t>();
  test_boundaries<uint128_t, uint32_t>();
  test_boundaries<int128_t, int64_t>();
  test_boundaries<int128_t, uint64_t>();
  test_boundaries<uint128_t, int64_t>();
  test_boundaries<uint128_t, uint64_t>();
  test_boundaries<int128_t, int128_t>();
  test_boundaries<uint128_t, uint128_t>();

  std::random_device rd;
  std::mt19937 gen(rd());

  std::uniform_int_distribution<int32_t> dist_i32(1, pstd::numeric_limits<int32_t>::max());
  std::uniform_int_distribution<uint64_t> dist_u64(0, pstd::numeric_limits<uint64_t>::max());

  // Test unsigned/signed
  for (int i = 0; i < 10000; i++)
  {
    uint64_t x = dist_i32(gen);
     int32_t y = dist_i32(gen);
    uint64_t res = fast_div(x, y);

    std::cout << "fast_div(" << x << ", " << y << ") = " << res;
    check(res == x / y);

    x = dist_u64(gen);
    y = dist_i32(gen);
    res = fast_div(x, y);

    std::cout << "fast_div(" << x << ", " << y << ") = " << res;
    check(res == x / y);
  }

  std::uniform_int_distribution<uint64_t> dist_u62(0, uint64_t(1ull << 62));

  // Test signed/signed
  for (int i = 0; i < 10000; i++)
  {
    // Test x < 2^64
    int128_t x = dist_u64(gen);
     int32_t y = dist_i32(gen);
    int128_t res = fast_div(x, y);

    std::cout << "fast_div(" << x << ", " << y << ") = " << res;
    check(res == x / y);

    // Test x > 2^64
    int128_t low = dist_u64(gen);
    int128_t high = int128_t(dist_u62(gen)) << 64;
    x = high | low;
    y = dist_i32(gen);
    res = fast_div(x, y);

    std::cout << "fast_div(" << x << ", " << y << ") = " << res;
    check(res == x / y);
  }

  std::cout << std::endl;
  std::cout << "All tests passed successfully!" << std::endl;

  return 0;
}
