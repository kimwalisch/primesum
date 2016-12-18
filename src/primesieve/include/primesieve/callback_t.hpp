/**
 *  @file   callback_t.hpp
 *  @brief  Callback types with C++ linkage and extern "C" linkage.
 * 
 *  Copyright (C) 2015 Kim Walisch, <kim.walisch@gmail.com>
 * 
 *  This file is distributed under the BSD License. See the COPYING
 *  file in the top level directory.
 */

#ifndef CALLBACK_t_HPP
#define CALLBACK_t_HPP

#include <stdint.h>

// C++ linkage
typedef void (*callback_t)(uint64_t);

extern "C"
{
  typedef void (*callback_c_t)(uint64_t);
}

#endif
