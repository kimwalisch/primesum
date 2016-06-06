primesum 256-bits
=================
[![Build Status](https://travis-ci.org/kimwalisch/primesum.svg)](https://travis-ci.org/kimwalisch/primesum)
[![GitHub license](https://img.shields.io/badge/license-BSD%202-blue.svg)](https://github.com/kimwalisch/primesum/blob/master/COPYING)

With the 256-bits version of primesum you can compute prime sums for
values of x&nbsp;>&nbsp;10<sup>20</sup>. The 256-bits version of primesum
has already been used to compute several world records.

Build instructions (Unix-like OSes)
-----------------------------------
You need to have installed a C++ compiler which supports OpenMP 4.0 or
later (e.g. GCC >= 5.0), GNU make and the
<a href="http://www.boost.org/">Boost C++ libraries</a>.

Download
[primesum-256-bits-1.0.tar.gz](https://dl.bintray.com/kimwalisch/primesum/primesum-256-bits-1.0.tar.gz)
and build it using:

```sh
$ ./build.sh
$ make check
$ sudo make install
```

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

[A046731](https://oeis.org/A046731) world records
-------------------------------------------------

<table>
  <tr align="center">
    <td><b>x</b></td>
    <td><b>Sum of the primes below x</b></td>
  </tr>
  </tr>
  <tr align="right">
    <td>10<sup>21</sup></td>
    <td>10,449,550,362,130,704,786,220,283,253,063,405,651,965</td>
  </tr>
  </tr>
  <tr align="right">
    <td>10<sup>22</sup></td>
    <td>996,973,504,763,259,668,279,213,971,353,794,878,368,213</td>
  </tr>
</table>

Note that the sum of the primes below 10^21 was already correctly
computed by Marc Deléglise in 2009 but when he verified the result
using his program he found a different result (off by 1) so he
withdrew his result in 2011.

[A099824](https://oeis.org/A099824) world records
-------------------------------------------------

<table>
  <tr align="center">
    <td><b>n</b></td>
    <td><b>Sum of the first n primes</b></td>
  </tr>
  </tr>
  <tr align="right">
    <td>10<sup>18</sup></td>
    <td>21,849,887,810,843,912,935,127,758,942,358,047,227</td>
  </tr>
  <tr align="right">
    <td>10<sup>19</sup></td>
    <td>2,302,808,849,326,957,165,657,230,565,155,878,163,277</td>
  </tr>
  <tr align="right">
    <td>10<sup>20</sup></td>
    <td>242,048,824,942,159,504,049,568,772,767,666,927,073,373</td>
  </tr>
</table>

References
----------
1. M. Deléglise and Jean-Louis Nicolas, Maximal product of primes whose sum is bounded, 3 Jul 2012, http://arxiv.org/abs/1207.0603.
