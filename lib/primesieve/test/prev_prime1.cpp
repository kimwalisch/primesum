///
/// @file   prev_prime1.cpp
/// @brief  Test prev_prime() of primesieve::iterator.
///
/// Copyright (C) 2017 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <primesieve.hpp>

#include <stdint.h>
#include <cstdlib>
#include <exception>
#include <iostream>
#include <vector>

using namespace std;

void check(bool OK)
{
  cout << "   " << (OK ? "OK" : "ERROR") << "\n";
  if (!OK)
    exit(1);
}

int main()
{
  vector<uint64_t> primes;
  primesieve::generate_primes(100000, &primes);
  primesieve::iterator it;
  uint64_t back = primes.size() - 1;
  uint64_t prime;

  for (uint64_t i = back; i > 0; i--)
  {
    it.skipto(primes[i] + 1);
    prime = it.prev_prime();
    cout << "prev_prime(" << primes[i] + 1 << ") = " << prime;
    check(prime == primes[i]);

    it.skipto(primes[i]);
    prime = it.prev_prime();
    cout << "prev_prime(" << primes[i] << ") = " << prime;
    check(prime == primes[i - 1]);
  }

  it.skipto(100000000);
  prime = it.prev_prime();
  uint64_t sum = 0;

  // iterate over the primes below 10^8
  for (; prime > 0; prime = it.prev_prime())
    sum += prime;

  cout << "Sum of the primes below 10^8 = " << sum;
  check(sum == 279209790387276ull);

  for (uint64_t i = 0; i < 1000; i++)
  {
    prime = it.prev_prime();
    cout << "prev_prime(0) = " << prime;
    check(prime == 0);
  }

  for (uint64_t i = 0; i < 1000; i++)
  {
    uint64_t old = prime;
    prime = it.next_prime();
    cout << "next_prime(" << old << ") = " << prime;
    check(prime == primes[i]);
  }

  it.skipto(primes.back());

  for (uint64_t i = 0; i < 1000; i++)
  {
    prime = it.prev_prime();
    uint64_t p1 = primes.size() - (i + 1);
    uint64_t p2 = primes.size() - (i + 2);
    cout << "prev_prime(" << primes[p1] << ") = " << prime;
    check(prime == primes[p2]);
  }

  for (uint64_t i = 0; i < 1000; i++)
  {
    uint64_t old = prime;
    uint64_t j = primes.size() - 1000 + i;
    prime = it.next_prime();
    cout << "next_prime(" << old << ") = " << prime;
    check(prime == primes[j]);
  }

  it.skipto(18446744073709551615ull, 18446744073709551557ull);
  prime = it.prev_prime();
  cout << "prev_prime(" << 18446744073709551615ull << ") = " << prime;
  check(prime == 18446744073709551557ull);

  cout << endl;
  cout << "All tests passed successfully!" << endl;

  return 0;
}
