primesum 256-bit
================
[![Build Status](https://travis-ci.org/kimwalisch/primesum.svg)](https://travis-ci.org/kimwalisch/primesum)
[![GitHub license](https://img.shields.io/badge/license-BSD%202-blue.svg)](https://github.com/kimwalisch/primesum/blob/master/COPYING)

primesum 256-bit allows to compute prime sums for values of
x&nbsp;>&nbsp;10<sup>20</sup>. No binaries are provided, it must be build
from source. primesum 256-bit has already been used to compute several
world records.

Build instructions (Unix-like OSes)
-----------------------------------
You need to have installed a C++ compiler which supports OpenMP 4.0 or
later (e.g. GCC ≥ 5.0, Clang ≥ 3.8), GNU make and the
<a href="http://www.boost.org/">Boost C++ libraries</a>.

Download 
[primesum-256-bit.tar.gz](https://github.com/kimwalisch/primesum/archive/256-bit.tar.gz)
and build it:

```sh
# Ubuntu/Debian prerequisites
sudo apt-get install g++ make cmake libboost-all-dev

cmake .
make -j8
sudo make install
```

Usage examples
--------------
Open a terminal and run primesum using e.g.:
```sh
# Sum the primes below 10^14
./primesum 1e14

# Print progress and status information during computation
./primesum 1e20 --status

# Use 4 threads
./primesum 1e14 --threads=4 --time
```

Command-line options
--------------------
```
Usage: primesum x [OPTION]...
Sum the primes below x <= 10^31 using fast implementations of the
combinatorial prime summing function.

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

References
----------
1. M. Deléglise and Jean-Louis Nicolas, Maximal product of primes whose sum is bounded, 3 Jul 2012, http://arxiv.org/abs/1207.0603.
