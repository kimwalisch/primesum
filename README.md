primesum
========
[![Build Status](https://travis-ci.org/kimwalisch/primesum.svg)](https://travis-ci.org/kimwalisch/primesum)
[![GitHub license](https://img.shields.io/badge/license-BSD%202-blue.svg)](https://github.com/kimwalisch/primesum/blob/master/COPYING)

**primesum** is a command-line program that computes the sum of the
primes below an integer x&nbsp;≤&nbsp;10<sup>20</sup> as quickly as
possible using a modified version of the combinatorial prime counting
function algorithm <a href="#references">[1]</a>.

**primesum** is a modified version of the author's
[primecount](https://github.com/kimwalisch/primecount) program.

Binaries
--------
Below are the latest precompiled binaries for Windows 64-bit and Linux x64.
These binaries are statically linked and require a CPU which supports the
POPCNT instruction (2010 or later).

* [primesum-1.1-win64.zip](https://github.com/kimwalisch/primesum/releases/download/v1.1/primesum-1.1-win64.zip), 382K
* [primesum-1.1-linux-x64.tar.gz](https://github.com/kimwalisch/primesum/releases/download/v1.1/primesum-1.1-linux-x64.tar.gz), 894K

primesum 256-bit
----------------
[primesum 256-bit](https://github.com/kimwalisch/primesum/tree/256-bit)
allows to compute prime sums for values of x&nbsp;>&nbsp;10<sup>20</sup>
but it runs only at about half the speed of the primesum 128-bit
version due to slower 256-bit integer arithmetic. primesum 256-bit has
already been used to compute many new prime sum world records!

* [primesum-1.1-256-win64.zip](https://github.com/kimwalisch/primesum/releases/download/v1.1-256-bit/primesum-1.1-256-win64.zip), 498K
* [primesum-1.1-256-linux-x64.tar.gz](https://github.com/kimwalisch/primesum/releases/download/v1.1-256-bit/primesum-1.1-256-linux-x64.tar.gz), 1005K

Build instructions
------------------
You need to have installed a C++ compiler, cmake and make. Ideally
primesum should be compiled using the GCC compiler as GCC supports both
OpenMP and 128-bit integers.

```sh
cmake .
make -j
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
Sum the primes below x <= 10^20 using fast implementations of the
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

Benchmark
---------
<table>
  <tr align="center">
    <td><b>x</b></td>
    <td><b>Sum of the primes below x</b></td>
    <td><b>Time elapsed</b></td>
  </tr>
  <tr align="right">
    <td>10<sup>10</sup></td>
    <td>2,220,822,432,581,729,238</td>
    <td>0.02s</td>
  </tr>
  <tr align="right">
    <td>10<sup>11</sup></td>
    <td>201,467,077,743,744,681,014</td>
    <td>0.03s</td>
  </tr>
  <tr align="right">
    <td>10<sup>12</sup></td>
    <td>18,435,588,552,550,705,911,377</td>
    <td>0.04s</td>
  </tr>
  <tr align="right">
    <td>10<sup>13</sup></td>
    <td>1,699,246,443,377,779,418,889,494</td>
    <td>0.13s</td>
  </tr>
  <tr align="right">
    <td>10<sup>14</sup></td>
    <td>157,589,260,710,736,940,541,561,021</td>
    <td>0.44s</td>
  </tr>
  <tr align="right">
    <td>10<sup>15</sup></td>
    <td>14,692,398,516,908,006,398,225,702,366</td>
    <td>1.36s</td>
  </tr>
  <tr align="right">
    <td>10<sup>16</sup></td>
    <td>1,376,110,854,313,351,899,159,632,866,552</td>
    <td>5.03s</td>
  </tr>
  <tr align="right">
    <td>10<sup>17</sup></td>
    <td>129,408,626,276,669,278,966,252,031,311,350</td>
    <td>24.05s</td>
  </tr>
  <tr align="right">
    <td>10<sup>18</sup></td>
    <td>12,212,914,292,949,226,570,880,576,733,896,687</td>
    <td>110.96s</td>
  </tr>
  <tr align="right">
    <td>10<sup>19</sup></td>
    <td>1,156,251,260,549,368,082,781,614,413,945,980,126</td>
    <td>438.50s</td>
  </tr>
  <tr align="right">
    <td>10<sup>20</sup></td>
    <td>109,778,913,483,063,648,128,485,839,045,703,833,541</td>
    <td>1909.077s</td>
  </tr>
</table>

The benchmarks above were run on an Intel Core i7-6700 CPU (4 x 3.4 GHz) from
2015 using a Linux x64 operating system and primesum was compiled using
GCC 5.2.

[A046731](https://oeis.org/A046731) world records
-------------------------------------------------

<table>
  <tr align="center">
    <td><b>x</b></td>
    <td><b>Sum of the primes below x</b></td>
    <td><b>Date</b></td>
    <td><b>Computed by</b></td>
  </tr>
  <tr align="right">
    <td>10<sup>21</sup></td>
    <td>10,449,550,362,130,704,786,220,283,253,063,405,651,965</td>
    <td>June 6, 2016</td>
    <td>Kim Walisch</td>
  </tr>
  <tr align="right">
    <td>10<sup>22</sup></td>
    <td>996,973,504,763,259,668,279,213,971,353,794,878,368,213</td>
    <td>June 6, 2016</td>
    <td>Kim Walisch</td>
  </tr>
  <tr align="right">
    <td>10<sup>23</sup></td>
    <td>95,320,530,117,168,404,458,544,684,912,403,185,555,509,650</td>
    <td>June 11, 2016</td>
    <td>Kim Walisch</td>
  </tr>
  <tr align="right">
    <td>10<sup>24</sup></td>
    <td>9,131,187,511,156,941,634,384,410,084,928,380,134,453,142,199</td>
    <td>June 17, 2016</td>
    <td>David Baugh</td>
  </tr>
  <tr align="right">
    <td>10<sup>25</sup></td>
    <td>876,268,031,750,623,105,684,911,815,303,505,535,704,119,354,853</td>
    <td>Oct. 16, 2016</td>
    <td>David Baugh</td>
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
    <td><b>Date</b></td>
    <td><b>Computed by</b></td>
  </tr>
  <tr align="right">
    <td>10<sup>18</sup></td>
    <td>21,849,887,810,843,912,935,127,758,942,358,047,227</td>
    <td>June 5, 2016</td>
    <td>Kim Walisch</td>
  </tr>
  <tr align="right">
    <td>10<sup>19</sup></td>
    <td>2,302,808,849,326,957,165,657,230,565,155,878,163,277</td>
    <td>June 5, 2016</td>
    <td>Kim Walisch</td>
  </tr>
  <tr align="right">
    <td>10<sup>20</sup></td>
    <td>242,048,824,942,159,504,049,568,772,767,666,927,073,373</td>
    <td>June 5, 2016</td>
    <td>Kim Walisch</td>
  </tr>
  <tr align="right">
    <td>10<sup>21</sup></td>
    <td>25,380,411,223,557,757,489,768,734,384,174,904,646,137,001</td>
    <td>June 11, 2016</td>
    <td>Kim Walisch</td>
  </tr>
  <tr align="right">
    <td>10<sup>22</sup></td>
    <td>2,655,479,563,137,417,712,148,525,630,666,125,075,397,977,159</td>
    <td>June 22, 2016</td>
    <td>David Baugh</td>
  </tr>
  <tr align="right">
    <td>10<sup>23</sup></td>
    <td>277,281,399,946,560,013,844,427,926,899,019,949,823,102,890,613</td>
    <td>Sep. 26, 2016</td>
    <td>David Baugh</td>
  </tr>
</table>

References
----------
1. M. Deléglise and Jean-Louis Nicolas, Maximal product of primes whose sum is bounded, 3 Jul 2012, http://arxiv.org/abs/1207.0603.
