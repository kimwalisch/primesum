///
/// @file   primesieve_error.hpp
/// @brief  The primesieve_error class is used for all
///         exceptions within primesieve.
///
/// Copyright (C) 2017 Kim Walisch, <kim.walisch@gmail.com>
///
/// This file is distributed under the BSD License. See the COPYING
/// file in the top level directory.
///

#ifndef PRIMESIEVE_ERROR_HPP
#define PRIMESIEVE_ERROR_HPP

#include <stdexcept>
#include <string>

namespace primesieve {

/// primesieve throws a primesieve_error exception
/// if an error occurs e.g. prime > 2^64.
///
class primesieve_error : public std::runtime_error
{
public:
  primesieve_error(const std::string& msg)
    : std::runtime_error(msg)
  { }
};

} // namespace

#endif
