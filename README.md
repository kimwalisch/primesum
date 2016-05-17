primesum
========
[![Build Status](https://travis-ci.org/kimwalisch/primesum.svg)](https://travis-ci.org/kimwalisch/primesum)
[![GitHub license](https://img.shields.io/badge/license-BSD%202-blue.svg)](https://github.com/kimwalisch/primesum/blob/master/COPYING)

**primesum** is a command-line program (and C++ library) that computes the
sum of the primes below an integer x&nbsp;≤&nbsp;10<sup>31</sup> as quickly
as possible. **primesum** is a modified version of the author's
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

Benchmark
---------
<table>
  <tr align="center">
    <td><b>x</b></td>
    <td><b>Prime sum</b></td>
    <td><b>Time elapsed</b></td>
  </tr>
  <tr align="right">
    <td>10<sup>10</sup></td>
    <td>2,220,822,432,581,729,238</td>
    <td>0.02s</td>
  </tr>
  </tr>
  <tr align="right">
    <td>10<sup>11</sup></td>
    <td>201,467,077,743,744,681,014</td>
    <td>0.06s</td>
  </tr>
  </tr>
  <tr align="right">
    <td>10<sup>12</sup></td>
    <td>18,435,588,552,550,705,911,377</td>
    <td>0.07s</td>
  </tr>
  </tr>
  <tr align="right">
    <td>10<sup>13</sup></td>
    <td>1,699,246,443,377,779,418,889,494</td>
    <td>0.26s</td>
  </tr>
  </tr>
  <tr align="right">
    <td>10<sup>14</sup></td>
    <td>157,589,260,710,736,940,541,561,021</td>
    <td>1.03s</td>
  </tr>
  <tr align="right">
    <td>10<sup>15</sup></td>
    <td>14,692,398,516,908,006,398,225,702,366</td>
    <td>4.58s</td>
  </tr>
  <tr align="right">
    <td>10<sup>16</sup></td>
    <td>1,376,110,854,313,351,899,159,632,866,552</td>
    <td>24.93s</td>
  </tr>
  <tr align="right">
    <td>10<sup>17</sup></td>
    <td>129,408,626,276,669,278,966,252,031,311,350</td>
    <td>144.45s</td>
  </tr>
  <tr align="right">
  <td>10<sup>18</sup></td>
  <td>12,212,914,292,949,226,570,880,576,733,896,687</td>
    <td>833.35s</td>
  </tr>
</table>

The benchmarks above were run on an Intel Core i7-6700 CPU (4 x 3.4 GHz) from
2015 using a Linux x64 operating system and primesum was compiled using
GCC 5.2.

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
