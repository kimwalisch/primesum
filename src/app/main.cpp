///
/// @file   main.cpp
/// @brief  primesum console application.
///
/// Copyright (C) 2018 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include "cmdoptions.hpp"

#include <primesum-internal.hpp>
#include <primesum.hpp>
#include <imath.hpp>
#include <int128_t.hpp>
#include <int256_t.hpp>
#include <PhiTiny.hpp>
#include <S1.hpp>
#include <S2.hpp>

#include <stdint.h>
#include <exception>
#include <iostream>
#include <limits>
#include <string>

using namespace std;
using namespace primesum;

namespace primesum {

int64_t int64_cast(int128_t x)
{
  if (x > numeric_limits<int64_t>::max())
    throw primesum_error("this is a 63-bit function, x must be < 2^63");
  return (int64_t) x;
}

int256_t P2(int128_t x, int threads)
{
  if (x < 1)
    return 0;

  double alpha = get_alpha_deleglise_rivat(x);
  string limit = get_max_x(alpha);

  if (x > to_int128(limit))
    throw primesum_error("P2(x): x must be <= " + limit);

  if (is_print())
    set_print_variables(true);

  int64_t y = (int64_t) (iroot<3>(x) * alpha);
  return P2(x, y, threads);
}

int256_t S1(int128_t x, int threads)
{
  if (x < 1)
    return 0;

  double alpha = get_alpha_deleglise_rivat(x);
  string limit = get_max_x(alpha);

  if (x > to_int128(limit))
    throw primesum_error("S1(x): x must be <= " + limit);

  if (is_print())
    set_print_variables(true);

  int64_t y = (int64_t) (iroot<3>(x) * alpha);
  int64_t c = PhiTiny::get_c(y);

  return S1(x, y, c, threads);
}

int256_t S2_trivial(int128_t x, int threads)
{
  if (x < 1)
    return 0;

  double alpha = get_alpha_deleglise_rivat(x);
  string limit = get_max_x(alpha);

  if (x > to_int128(limit))
    throw primesum_error("S2_trivial(x): x must be <= " + limit);

  if (is_print())
    set_print_variables(true);

  int64_t y = (int64_t) (iroot<3>(x) * alpha);
  int64_t z = (int64_t) (x / y);
  int64_t c = PhiTiny::get_c(y);

  if (x <= numeric_limits<int64_t>::max())
    return S2_trivial((int64_t) x, y, z, c, threads);
  else
    return S2_trivial(x, y, z, c, threads);
}

int256_t S2_easy(int128_t x, int threads)
{
  if (x < 1)
    return 0;

  double alpha = get_alpha_deleglise_rivat(x);
  string limit = get_max_x(alpha);

  if (x > to_int128(limit))
    throw primesum_error("S2_easy(x): x must be <= " + limit);

  if (is_print())
    set_print_variables(true);

  int64_t y = (int64_t) (iroot<3>(x) * alpha);
  int64_t z = (int64_t) (x / y);
  int64_t c = PhiTiny::get_c(y);

  if (x <= numeric_limits<int64_t>::max())
    return S2_easy((int64_t) x, y, z, c, threads);
  else
    return S2_easy(x, y, z, c, threads);
}

int256_t S2_hard(int128_t x, int threads)
{
  if (x < 1)
    return 0;

  double alpha = get_alpha_deleglise_rivat(x);
  string limit = get_max_x(alpha);

  if (x > to_int128(limit))
    throw primesum_error("S2_hard(x): x must be <= " + limit);

  if (is_print())
    set_print_variables(true);

  int64_t y = (int64_t) (iroot<3>(x) * alpha);
  int64_t z = (int64_t) (x / y);
  int64_t c = PhiTiny::get_c(y);

  return S2_hard(x, y, z, c, threads);
}

} // namespace

int main (int argc, char* argv[])
{
  PrimeSumOptions pco = parseOptions(argc, argv);
  double time = get_time();

  int128_t x = pco.x;
  int256_t res = 0;
  int threads = get_num_threads();

  try
  {
    switch (pco.option)
    {
      case OPTION_DELEGLISE_RIVAT:
        res = pi_deleglise_rivat(x, threads); break;
      case OPTION_DELEGLISE_RIVAT_PARALLEL1:
        res = pi_deleglise_rivat_parallel1(x, threads); break;
      case OPTION_LMO:
        res = pi_lmo(int64_cast(x), threads); break;
      case OPTION_LMO1:
        res = pi_lmo1(int64_cast(x)); break;
      case OPTION_LMO2:
        res = pi_lmo2(int64_cast(x)); break;
      case OPTION_LMO3:
        res = pi_lmo3(int64_cast(x)); break;
      case OPTION_LMO4:
        res = pi_lmo4(int64_cast(x)); break;
      case OPTION_LMO5:
        res = pi_lmo5(int64_cast(x)); break;
      case OPTION_LMO_PARALLEL1:
        res = pi_lmo_parallel1(int64_cast(x), threads); break;
      case OPTION_P2:
        res = P2(x, threads); break;
      case OPTION_PI:
        res = pi(x, threads); break;
      case OPTION_S1:
        res = S1(x, threads); break;
      case OPTION_S2_EASY:
        res = S2_easy(x, threads); break;
      case OPTION_S2_HARD:
        res = S2_hard(x, threads); break;
      case OPTION_S2_TRIVIAL:
        res = S2_trivial(x, threads); break;
    }
  }
  catch (bad_alloc&)
  {
    cerr << "Error: failed to allocate memory, your system most likely does" << endl
         << "       not have enough memory to run this computation." << endl;
    return 1;
  }
  catch (exception& e)
  {
    cerr << "Error: " << e.what() << endl;
    return 1;
  }

  if (print_result())
  {
    if (is_print())
      cout << endl;
    cout << res << endl;
    if (pco.time)
      print_seconds(get_time() - time);
  }

  return 0;
}
