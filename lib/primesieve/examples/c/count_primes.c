/** @example count_primes.c
 *  C program that shows how to count primes. */

#include <primesieve.h>
#include <inttypes.h>
#include <stdio.h>

int main()
{
  uint64_t count = primesieve_count_primes(0, 1000);
  printf("Primes below 1000 = %" PRIu64 "\n", count);

  return 0;
}
