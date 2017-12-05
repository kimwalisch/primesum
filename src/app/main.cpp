///
/// @file   main.cpp
/// @brief  primesum console application.
///
/// Copyright (C) 2016 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include "cmdoptions.hpp"

#include <primesum-internal.hpp>
#include <primesum.hpp>
#include <imath.hpp>
#include <int128_t.hpp>
#include <PhiTiny.hpp>
#include <S1.hpp>
#include <S2.hpp>

#include <stdint.h>
#include <exception>
#include <iostream>
#include <limits>
#include <string>

#ifdef HAVE_MPI
  #include <mpi.h>
#endif

using namespace std;
using namespace primesum;

namespace primesum {

int64_t int64_cast(maxint_t x)
{
  if (x > numeric_limits<int64_t>::max())
    throw primesum_error("this is a 63-bit function, x must be < 2^63");
  return (int64_t) x;
}

maxint_t P2(maxint_t x, int threads)
{
  if (x < 1)
    return 0;

  double alpha = get_alpha_deleglise_rivat(x);
  string limit = get_max_x(alpha);

  if (x > to_maxint(limit))
    throw primesum_error("P2(x): x must be <= " + limit);

  if (print_status())
    set_print_variables(true);

  int64_t y = (int64_t) (iroot<3>(x) * alpha);
  return P2(x, y, threads);
}

maxint_t S1(maxint_t x, int threads)
{
  if (x < 1)
    return 0;

  double alpha = get_alpha_deleglise_rivat(x);
  string limit = get_max_x(alpha);

  if (x > to_maxint(limit))
    throw primesum_error("S1(x): x must be <= " + limit);

  if (print_status())
    set_print_variables(true);

  int64_t y = (int64_t) (iroot<3>(x) * alpha);
  int64_t c = PhiTiny::get_c(y);

  return S1(x, y, c, threads);
}

maxint_t S2_trivial(maxint_t x, int threads)
{
  if (x < 1)
    return 0;

  double alpha = get_alpha_deleglise_rivat(x);
  string limit = get_max_x(alpha);

  if (x > to_maxint(limit))
    throw primesum_error("S2_trivial(x): x must be <= " + limit);

  if (print_status())
    set_print_variables(true);

  int64_t y = (int64_t) (iroot<3>(x) * alpha);
  int64_t z = (int64_t) (x / y);
  int64_t c = PhiTiny::get_c(y);

  if (x <= numeric_limits<int64_t>::max())
    return S2_trivial((int64_t) x, y, z, c, threads);
  else
    return S2_trivial(x, y, z, c, threads);
}

maxint_t S2_easy(maxint_t x, int threads)
{
  if (x < 1)
    return 0;

  double alpha = get_alpha_deleglise_rivat(x);
  string limit = get_max_x(alpha);

  if (x > to_maxint(limit))
    throw primesum_error("S2_easy(x): x must be <= " + limit);

  if (print_status())
    set_print_variables(true);

  int64_t y = (int64_t) (iroot<3>(x) * alpha);
  int64_t z = (int64_t) (x / y);
  int64_t c = PhiTiny::get_c(y);

  if (x <= numeric_limits<int64_t>::max())
    return S2_easy((int64_t) x, y, z, c, threads);
  else
    return S2_easy(x, y, z, c, threads);
}

maxint_t S2_hard(maxint_t x, int threads)
{
  if (x < 1)
    return 0;

  double alpha = get_alpha_deleglise_rivat(x);
  string limit = get_max_x(alpha);

  if (x > to_maxint(limit))
    throw primesum_error("S2_hard(x): x must be <= " + limit);

  if (print_status())
    set_print_variables(true);

  int64_t y = (int64_t) (iroot<3>(x) * alpha);
  int64_t z = (int64_t) (x / y);
  int64_t c = PhiTiny::get_c(y);

  return S2_hard(x, y, z, c, threads);
}

} // namespace

int main (int argc, char* argv[])
{
#ifdef HAVE_MPI
  MPI_Init(&argc, &argv);
#endif

  PrimeSumOptions pco = parseOptions(argc, argv);
  double time = get_wtime();

  maxint_t x = pco.x;
  maxint_t res = 0;
  int threads = pco.threads;

  try
  {
    if (x > numeric_limits<uint64_t>::max())
      throw primesum_error("this primesum version only works up to 2^64-1");

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
#ifdef HAVE_MPI
    MPI_Finalize();
#endif
    cerr << "Error: failed to allocate memory, your system most likely does" << endl
         << "       not have enough memory to run this computation." << endl;
    return 1;
  }
  catch (exception& e)
  {
#ifdef HAVE_MPI
    MPI_Finalize();
#endif
    cerr << "Error: " << e.what() << endl;
    return 1;
  }

  if (print_result())
  {
    if (print_status())
      cout << endl;
    cout << res << endl;
    if (pco.time)
      print_seconds(get_wtime() - time);
  }

#ifdef HAVE_MPI
    MPI_Finalize();
#endif

  return 0;
}
