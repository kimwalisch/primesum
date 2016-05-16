///
/// @file  cmdoptions.hpp
///
/// Copyright (C) 2016 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#ifndef PRIMESUM_CMDOPTIONS_HPP
#define PRIMESUM_CMDOPTIONS_HPP

#include <primesum.hpp>
#include <int128.hpp>
#include <stdint.h>

namespace primesum {

enum OptionValues
{
  OPTION_ALPHA,
  OPTION_DELEGLISE_RIVAT,
  OPTION_DELEGLISE_RIVAT_PARALLEL3,
  OPTION_HELP,
  OPTION_LMO,
  OPTION_LMO1,
  OPTION_LMO2,
  OPTION_LMO4,
  OPTION_LMO5,
  OPTION_LMO_PARALLEL2,
  OPTION_NUMBER,
  OPTION_P2,
  OPTION_PI,
  OPTION_S1,
  OPTION_S2_EASY,
  OPTION_S2_HARD,
  OPTION_S2_TRIVIAL,
  OPTION_STATUS,
  OPTION_TEST,
  OPTION_TIME,
  OPTION_THREADS,
  OPTION_VERSION
};

struct PrimeSumOptions
{
  maxint_t x;
  int64_t option;
  bool time;
  int threads;
  PrimeSumOptions() :
    x(-1),
    option(OPTION_PI),
    time(false),
    threads(get_num_threads())
  { }
};

PrimeSumOptions parseOptions(int, char**);

} // namespace primesum

#endif
