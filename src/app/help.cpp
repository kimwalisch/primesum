///
/// @file  help.cpp
///
/// Copyright (C) 2017 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include <primesum.hpp>

#include <iostream>
#include <cstdlib>
#include <string>

using namespace std;

namespace {

const string helpMenu(
  "Usage: primesum x [OPTION]...\n"
  "Sum the primes below x <= 10^20 using fast implementations of the\n"
  "combinatorial prime summing function.\n"
  "\n"
  "Options:\n"
  "\n"
  "  -d,    --deleglise_rivat  Sum primes using Deleglise-Rivat algorithm\n"
  "  -l,    --lmo              Sum primes using Lagarias-Miller-Odlyzko\n"
  "  -s[N], --status[=N]       Show computation progress 1%, 2%, 3%, ...\n"
  "                            [N] digits after decimal point e.g. N=1, 99.9%\n"
  "         --test             Run various correctness tests and exit\n"
  "         --time             Print the time elapsed in seconds\n"
  "  -t<N>, --threads=<N>      Set the number of threads, 1 <= N <= CPU cores\n"
  "  -v,    --version          Print version and license information\n"
  "  -h,    --help             Print this help menu\n"
  "\n"
  "Advanced Deleglise-Rivat options:\n"
  "\n"
  "  -a<N>, --alpha=<N>        Tuning factor, 1 <= alpha <= x^(1/6)\n"
  "         --P2               Only compute the 2nd partial sieve function\n"
  "         --S1               Only compute the ordinary leaves\n"
  "         --S2_trivial       Only compute the trivial special leaves\n"
  "         --S2_easy          Only compute the easy special leaves\n"
  "         --S2_hard          Only compute the hard special leaves\n"
  "\n"
  "Examples:\n"
  "\n"
  "  primesum 1e13\n"
  "  primesum 1e13 --status --threads=4"
);

const string versionInfo(
  "primesum " PRIMESUM_VERSION ", <https://github.com/kimwalisch/primesum>\n"
  "Copyright (C) 2016 - 2017 Kim Walisch\n"
  "BSD 2-Clause License <http://opensource.org/licenses/BSD-2-Clause>"
);

} // end namespace

namespace primesum {

void help()
{
  cout << helpMenu << endl;
  exit(1);
}

void version()
{
  cout << versionInfo << endl;
  exit(1);
}

} // namespace
