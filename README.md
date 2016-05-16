primesum
========
[![Build Status](https://travis-ci.org/kimwalisch/primesum.svg)](https://travis-ci.org/kimwalisch/primesum)
[![GitHub license](https://img.shields.io/badge/license-BSD%202-blue.svg)](https://github.com/kimwalisch/primesum/blob/master/COPYING)

primesum is a command-line program (and C++ library) that computes the
sum of the primes below an integer x&nbsp;≤&nbsp;10<sup>31</sup> as quickly
as possible. **primesum** is in fact a modified version of the author's
[primecount](https://github.com/kimwalisch/primecount) program.

Binaries
--------
Will come within a week or so...

Usage examples
--------------
Open a terminal and run primesum using e.g.:
```sh
# Sum the primes below 10^14
$ ./primesum 1e14

# Print progress and status information during computation
$ ./primesum 1e20 --status

# Use 4 threads
$ ./primesum 1e14 --threads=4 --time
```

Command-line options
--------------------
```
Usage: primesum x [OPTION]...
Sum the primes below x <= 10^31 using fast implementations of the
combinatorial prime counting function.

Options:

  -d,    --deleglise_rivat  Sum primes using Deleglise-Rivat algorithm
  -l,    --lmo              Sum primes using Lagarias-Miller-Odlyzko
  -s[N], --status[=N]       Show computation progress 1%, 2%, 3%, ...
                            [N] digits after decimal point e.g. N=1, 99.9%
         --test             Run various correctness tests and exit
         --time             Print the time elapsed in seconds
  -t<N>, --threads=<N>      Set the number of threads, 1 <= N <= CPU cores
  -v,    --version          Print version and license information
  -h,    --help             Print this help menu

Advanced Deleglise-Rivat options:

  -a<N>, --alpha=<N>        Tuning factor, 1 <= alpha <= x^(1/6)
         --P2               Only compute the 2nd partial sieve function
         --S1               Only compute the ordinary leaves
         --S2_trivial       Only compute the trivial special leaves
         --S2_easy          Only compute the easy special leaves
         --S2_hard          Only compute the hard special leaves
```

Benchmarks
----------
TODO.

Build instructions (Unix-like OSes)
-----------------------------------
You need to have installed a C++ compiler, GNU make and GNU Autotools
(automake, autoconf, libtool) to build ```primesum```.

```sh
$ ./build.sh
$ make check
$ sudo make install
```

References
----------
1. M. Deléglise and Jean-Louis Nicolas, Maximal product of primes whose sum is bounded, Bull. Proc. of the Steklov Institute 17 (2013), 82-112.
