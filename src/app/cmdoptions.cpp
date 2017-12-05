///
/// @file   cmdoptions.cpp
/// @brief  Parse command-line options for the primesum console
///         (terminal) application.
///
/// Copyright (C) 2015 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#include "cmdoptions.hpp"
#include <primesum-internal.hpp>
#include <int128_t.hpp>

#include <stdint.h>
#include <vector>
#include <string>
#include <sstream>
#include <map>
#include <exception>
#include <cstdlib>
#include <cstddef>

using std::string;
using std::istringstream;

namespace primesum {

void help();
void version();
bool test();

double to_double(const string& s)
{
  istringstream i(s);
  double x;
  if (!(i >> x))
    throw primesum_error("failed to convert \"" + s + "\" to double");
  return x;
}

/// e.g. id = "--threads", value = "4"
struct Option
{
  string id;
  string value;
  template <typename T>
  T getValue() const
  {
    return (T) to_int128(value);
  }
};

/// Command-line options
std::map<string, OptionValues> optionMap;

void initOptionMap()
{
  optionMap["-a"]                          = OPTION_ALPHA;
  optionMap["--alpha"]                     = OPTION_ALPHA;
  optionMap["-d"]                          = OPTION_DELEGLISE_RIVAT;
  optionMap["--deleglise_rivat"]           = OPTION_DELEGLISE_RIVAT;
  optionMap["--deleglise_rivat_parallel1"] = OPTION_DELEGLISE_RIVAT_PARALLEL1;
  optionMap["-h"]                          = OPTION_HELP;
  optionMap["--help"]                      = OPTION_HELP;
  optionMap["-l"]                          = OPTION_LMO;
  optionMap["--lmo"]                       = OPTION_LMO;
  optionMap["--lmo1"]                      = OPTION_LMO1;
  optionMap["--lmo2"]                      = OPTION_LMO2;
  optionMap["--lmo3"]                      = OPTION_LMO3;
  optionMap["--lmo4"]                      = OPTION_LMO4;
  optionMap["--lmo5"]                      = OPTION_LMO5;
  optionMap["--lmo_parallel1"]             = OPTION_LMO_PARALLEL1;
  optionMap["--number"]                    = OPTION_NUMBER;
  optionMap["--P2"]                        = OPTION_P2;
  optionMap["--pi"]                        = OPTION_PI;
  optionMap["--S1"]                        = OPTION_S1;
  optionMap["--S2_easy"]                   = OPTION_S2_EASY;
  optionMap["--S2_hard"]                   = OPTION_S2_HARD;
  optionMap["--S2_trivial"]                = OPTION_S2_TRIVIAL;
  optionMap["-s"]                          = OPTION_STATUS;
  optionMap["--status"]                    = OPTION_STATUS;
  optionMap["--test"]                      = OPTION_TEST;
  optionMap["--time"]                      = OPTION_TIME;
  optionMap["-t"]                          = OPTION_THREADS;
  optionMap["--threads"]                   = OPTION_THREADS;
  optionMap["-v"]                          = OPTION_VERSION;
  optionMap["--version"]                   = OPTION_VERSION;
}

/// e.g. "--threads=8" -> { id = "--threads", value = "8" }
Option makeOption(const string& str)
{
  Option option;
  size_t delimiter = string::npos;
  if (optionMap.count(str) == 0)
    delimiter = str.find_first_of("=0123456789");

  if (delimiter == string::npos)
    option.id = str;
  else
  {
    option.id = str.substr(0, delimiter);
    option.value = str.substr(delimiter + (str.at(delimiter) == '=' ? 1 : 0));
  }
  if (option.id.empty() && !option.value.empty())
    option.id = "--number";
  if (optionMap.count(option.id) == 0)
    option.id = "--help";

  return option;
}

PrimeSumOptions parseOptions(int argc, char** argv)
{
  initOptionMap();
  PrimeSumOptions pco;
  std::vector<int128_t> numbers;

  try
  {
    // iterate over the command-line options
    for (int i = 1; i < argc; i++)
    {
      Option option = makeOption(argv[i]);
      switch (optionMap[option.id])
      {
        case OPTION_ALPHA:   set_alpha(to_double(option.value)); break;
        case OPTION_NUMBER:  numbers.push_back(option.getValue<int128_t>()); break;
        case OPTION_THREADS: pco.threads = option.getValue<int>(); break;
        case OPTION_HELP:    help(); break;
        case OPTION_STATUS:  set_print(true);
                             if (!option.value.empty())
                                set_status_precision(option.getValue<int>());
                             pco.time = true;
                             break;
        case OPTION_TIME:    pco.time = true; break;
        case OPTION_TEST:    if (test()) exit(0); exit(1);
        case OPTION_VERSION: version(); break;
        default:             pco.option = optionMap[option.id];
      }
    }
  }
  catch (std::exception&)
  {
    help();
  }

  if (numbers.size() == 1)
    pco.x = numbers[0];
  else
    help();

  return pco;
}

} // namespace
