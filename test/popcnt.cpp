///
/// @file  popcnt.cpp
/// @brief Test scalar population counts.
///
/// Copyright (C) 2026 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <popcnt.hpp>

#include <stdint.h>
#include <cstdlib>
#include <iostream>
#include <random>

void check(uint64_t x)
{
  uint64_t expected = 0;
  for (uint64_t bits = x; bits; bits >>= 1)
    expected += bits & 1;

  if (popcnt64(x) != expected)
  {
    std::cerr << "Incorrect popcnt64(" << x << ")" << std::endl;
    std::exit(1);
  }
}

int main()
{
  check(0);
  check(~0ull);
  check(0x5555555555555555);
  check(0xaaaaaaaaaaaaaaaa);

  for (int i = 0; i < 64; i++)
  {
    check(1ull << i);
    check(~(1ull << i));
  }

  std::mt19937_64 gen(0);
  for (int i = 0; i < 10000; i++)
    check(gen());

  std::cout << "All tests passed successfully!" << std::endl;
  return 0;
}
