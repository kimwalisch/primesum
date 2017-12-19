# primesum

[![Build Status](https://travis-ci.org/kimwalisch/primesum.svg)](https://travis-ci.org/kimwalisch/primesum)
[![Github Releases](https://img.shields.io/github/release/kimwalisch/primesum.svg)](https://github.com/kimwalisch/primesum/releases)

**primesum** is a command-line program that computes the sum of the
primes below an integer x&nbsp;≤&nbsp;10<sup>31</sup> as quickly as
possible using a modified version of the combinatorial prime counting
function algorithm <a href="#references">[1]</a>. primesum has already
been used to compute many new [prime sum world records](#a046731-world-records)!

**primesum** is a modified version of the author's
[primecount](https://github.com/kimwalisch/primecount) program.

## Binaries

Below are the latest precompiled binaries for Windows, Linux and macOS.
These binaries are statically linked and require a CPU which supports the
POPCNT instruction (2010 or later).

* [primesum-1.3-win64.zip](https://github.com/kimwalisch/primesum/releases/download/v1.3/primesum-1.3-win64.zip), 529 KB
* [primesum-1.3-linux-x64.tar.gz](https://github.com/kimwalisch/primesum/releases/download/v1.3/primesum-1.3-linux-x64.tar.gz), 1.1 MB
* [primesum-1.3-macOS-x64.zip](https://github.com/kimwalisch/primesum/releases/download/v1.3/primesum-1.3-macOS-x64.zip), 956 KB

## Build instructions

You need to have installed a C++ compiler, cmake and make. Ideally
primesum should be compiled using the GCC compiler as GCC supports both
OpenMP and 128-bit integers.

```sh
cmake .
make -j
sudo make install
```

## Usage examples

Open a terminal and run primesum using e.g.:
```sh
# Sum the primes below 10^14
./primesum 1e14

# Print progress and status information during computation
./primesum 1e18 --status

# Use 4 threads
./primesum 1e14 --threads=4 --time
```

## Command-line options

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

## Benchmark

<table>
  <tr align="center">
    <td><b>x</b></td>
    <td><b>Sum of the primes below x</b></td>
    <td><b>Time elapsed</b></td>
  </tr>
  <tr align="right">
    <td>10<sup>10</sup></td>
    <td>2,220,822,432,581,729,238</td>
    <td>0.01s</td>
  </tr>
  <tr align="right">
    <td>10<sup>11</sup></td>
    <td>201,467,077,743,744,681,014</td>
    <td>0.02s</td>
  </tr>
  <tr align="right">
    <td>10<sup>12</sup></td>
    <td>18,435,588,552,550,705,911,377</td>
    <td>0.04s</td>
  </tr>
  <tr align="right">
    <td>10<sup>13</sup></td>
    <td>1,699,246,443,377,779,418,889,494</td>
    <td>0.11s</td>
  </tr>
  <tr align="right">
    <td>10<sup>14</sup></td>
    <td>157,589,260,710,736,940,541,561,021</td>
    <td>0.36s</td>
  </tr>
  <tr align="right">
    <td>10<sup>15</sup></td>
    <td>14,692,398,516,908,006,398,225,702,366</td>
    <td>1.16s</td>
  </tr>
  <tr align="right">
    <td>10<sup>16</sup></td>
    <td>1,376,110,854,313,351,899,159,632,866,552</td>
    <td>3.66s</td>
  </tr>
  <tr align="right">
    <td>10<sup>17</sup></td>
    <td>129,408,626,276,669,278,966,252,031,311,350</td>
    <td>14.60s</td>
  </tr>
  <tr align="right">
    <td>10<sup>18</sup></td>
    <td>12,212,914,292,949,226,570,880,576,733,896,687</td>
    <td>66.66s</td>
  </tr>
  <tr align="right">
    <td>10<sup>19</sup></td>
    <td>1,156,251,260,549,368,082,781,614,413,945,980,126</td>
    <td>330.01s</td>
  </tr>
  <tr align="right">
    <td>10<sup>20</sup></td>
    <td>109,778,913,483,063,648,128,485,839,045,703,833,541</td>
    <td>1486.87s</td>
  </tr>
</table>

The benchmarks above were run on an Intel Core i7-6700 CPU (4 x 3.4 GHz) from
2015 using a Linux x64 operating system and primesum was compiled using
GCC 6.4.

## [A046731](https://oeis.org/A046731) world records

<table>
  <tr align="center">
    <td><b>x</b></td>
    <td><b>Sum of the primes below x</b></td>
    <td><b>Date</b></td>
    <td><b>Computed by</b></td>
  </tr>
  <tr align="right">
    <td>10<sup>21</sup></td>
    <td>10449550362130704786220283253063405651965</td>
    <td>June 6, 2016</td>
    <td>Kim Walisch</td>
  </tr>
  <tr align="right">
    <td>10<sup>22</sup></td>
    <td>996973504763259668279213971353794878368213</td>
    <td>June 6, 2016</td>
    <td>Kim Walisch</td>
  </tr>
  <tr align="right">
    <td>10<sup>23</sup></td>
    <td>95320530117168404458544684912403185555509650</td>
    <td>June 11, 2016</td>
    <td>Kim Walisch</td>
  </tr>
  <tr align="right">
    <td>10<sup>24</sup></td>
    <td>9131187511156941634384410084928380134453142199</td>
    <td>June 17, 2016</td>
    <td>David Baugh</td>
  </tr>
  <tr align="right">
    <td>10<sup>25</sup></td>
    <td>876268031750623105684911815303505535704119354853</td>
    <td>Oct. 16, 2016</td>
    <td>David Baugh</td>
  </tr>
</table>

Note that the sum of the primes below 10^21 was already correctly
computed by Marc Deléglise in 2009 but when he verified the result
using his program he found a different result (off by 1) so he
withdrew his result in 2011.

## [A099824](https://oeis.org/A099824) world records

<table>
  <tr align="center">
    <td><b>n</b></td>
    <td><b>Sum of the first n primes</b></td>
    <td><b>Date</b></td>
    <td><b>Computed by</b></td>
  </tr>
  <tr align="right">
    <td>10<sup>18</sup></td>
    <td>21849887810843912935127758942358047227</td>
    <td>June 5, 2016</td>
    <td>Kim Walisch</td>
  </tr>
  <tr align="right">
    <td>10<sup>19</sup></td>
    <td>2302808849326957165657230565155878163277</td>
    <td>June 5, 2016</td>
    <td>Kim Walisch</td>
  </tr>
  <tr align="right">
    <td>10<sup>20</sup></td>
    <td>242048824942159504049568772767666927073373</td>
    <td>June 5, 2016</td>
    <td>Kim Walisch</td>
  </tr>
  <tr align="right">
    <td>10<sup>21</sup></td>
    <td>25380411223557757489768734384174904646137001</td>
    <td>June 11, 2016</td>
    <td>Kim Walisch</td>
  </tr>
  <tr align="right">
    <td>10<sup>22</sup></td>
    <td>2655479563137417712148525630666125075397977159</td>
    <td>June 22, 2016</td>
    <td>David Baugh</td>
  </tr>
  <tr align="right">
    <td>10<sup>23</sup></td>
    <td>277281399946560013844427926899019949823102890613</td>
    <td>Sep. 26, 2016</td>
    <td>David Baugh</td>
  </tr>
</table>

## References

1. M. Deléglise and Jean-Louis Nicolas, Maximal product of primes whose sum is bounded, 3 Jul 2012, http://arxiv.org/abs/1207.0603.
