///
/// @file   count_primes2.cpp
/// @brief  Count the primes within [10^i, 10^i + 10^8]
///         for i = 12 to 19
///
/// Copyright (C) 2018 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <primesieve.hpp>

#include <stdint.h>
#include <array>
#include <cstdlib>
#include <cmath>
#include <iostream>
#include <iomanip>

using namespace std;
using namespace primesieve;

const array<uint64_t, 6> pix =
{
  3618282, // pi[10^12, 10^12+10^8]
  3342093, // pi[10^13, 10^13+10^8]
  3102679, // pi[10^14, 10^14+10^8]
  2893937, // pi[10^15, 10^15+10^8]
  2714904, // pi[10^16, 10^16+10^8]
  2555873  // pi[10^17, 10^17+10^8]
};

void check(bool OK)
{
  cout << "   " << (OK ? "OK" : "ERROR") << "\n";
  if (!OK)
    exit(1);
}

int main()
{
  cout << left;

  for (size_t i = 0; i < pix.size(); i++)
  {
    size_t j = i + 12;
    cout << "Sieving the primes within [10^" << j << ", 10^" << j << " + 10^8]" << endl;
    uint64_t start = (uint64_t) pow(10.0, j);
    uint64_t stop = start + (uint64_t) 1e8;
    uint64_t count = count_primes(start, stop);
    cout << "\rPrime count: " << setw(7) << count;
    check(count == pix[i]);
  }

  cout << endl;
  cout << "All tests passed successfully!" << endl;

  return 0;
}
