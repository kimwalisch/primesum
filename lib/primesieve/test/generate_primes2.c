///
/// @file   generate_primes2.c
/// @brief  Test prime number generation.
///
/// Copyright (C) 2017 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <primesieve.h>

#include <stdint.h>
#include <inttypes.h>
#include <stdio.h>
#include <stdlib.h>

// primes inside [0, 100]
const uint64_t small_primes[25] =
{
   2,  3,  5,  7, 11,
  13, 17, 19, 23, 29,
  31, 37, 41, 43, 47,
  53, 59, 61, 67, 71,
  73, 79, 83, 89, 97
};

// primes inside [18446744073709550681, 18446744073709551533]
const uint64_t large_primes[19] =
{
  18446744073709550681ull,
  18446744073709550717ull,
  18446744073709550719ull,
  18446744073709550771ull,
  18446744073709550773ull,
  18446744073709550791ull,
  18446744073709550873ull,
  18446744073709551113ull,
  18446744073709551163ull,
  18446744073709551191ull,
  18446744073709551253ull,
  18446744073709551263ull,
  18446744073709551293ull,
  18446744073709551337ull,
  18446744073709551359ull,
  18446744073709551427ull,
  18446744073709551437ull,
  18446744073709551521ull,
  18446744073709551533ull
};

void check(int OK)
{
  if (OK)
    printf("   OK\n");
  else
  {
    printf("   ERROR\n");
    exit(1);
  }
}

int main()
{
  size_t i;
  size_t size = 0;
  uint64_t* primes = (uint64_t*) primesieve_generate_primes(0, 100, &size, UINT64_PRIMES);
  printf("primes.size = %zu", size);
  check(size == 25);

  for (i = 0; i < size; i++)
  {
    printf("primes[%zu] = %" PRIu64, i, primes[i]);
    check(primes[i] == small_primes[i]);
  }

  primesieve_free(primes);
  primes = (uint64_t*) primesieve_generate_primes(18446744073709550672ull, 18446744073709551556ull, &size, UINT64_PRIMES);
  printf("primes.size = %zu", size);
  check(size == 19);

  for (i = 0; i < size; i++)
  {
    printf("primes[%zu] = %" PRIu64, i, primes[i]);
    check(primes[i] == large_primes[i]);
  }

  primesieve_free(primes);
  printf("\n");
  printf("All tests passed successfully!\n");

  return 0;
}
