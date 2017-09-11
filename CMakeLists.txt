cmake_minimum_required(VERSION 3.1)
project(primesum)
set(PRIMESUM_VERSION_MAJOR 1)
set(PRIMESUM_VERSION_MINOR 2)
set(PRIMESUM_VERSION "${PRIMESUM_VERSION_MAJOR}.${PRIMESUM_VERSION_MINOR}")
include_directories(include src/primesieve/include)
set(CMAKE_BUILD_TYPE Release)

# Build options ######################################################

option(WITH_POPCNT       "Enable x86 POPCNT instruction"  ON)
option(WITH_LIBDIVIDE    "Use libdivide.h"                ON)
option(WITH_OPENMP       "Enable OpenMP support"          ON)
option(BUILD_STATIC_LIBS "Build static libprimesieve"     ON)
option(BUILD_SHARED_LIBS "Build shared libprimesieve"     OFF)

# Compiler must support C++11 or later ###############################

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Include Check* #####################################################

include(CheckCXXCompilerFlag)
include(CheckTypeSize)

# primesum binary source files #######################################

set(SRC src/app/cmdoptions.cpp
        src/app/help.cpp
        src/app/main.cpp
        src/BitSieve.cpp
        src/FactorTable.cpp
        src/generate.cpp
        src/nth_prime.cpp
        src/P2.cpp
        src/PhiTiny.cpp
        src/PiTable.cpp
        src/phi.cpp
        src/phi_sum.cpp
        src/pi_primesieve.cpp
        src/pi_legendre.cpp
        src/primesum.cpp
        src/print.cpp
        src/S1.cpp
        src/S2LoadBalancer.cpp
        src/S2Status.cpp
        src/test.cpp
        src/Wheel.cpp
        src/deleglise-rivat/S2_trivial.cpp
        src/deleglise-rivat/S2_hard.cpp
        src/deleglise-rivat/pi_deleglise_rivat_parallel1.cpp
        src/lmo/pi_lmo1.cpp
        src/lmo/pi_lmo2.cpp
        src/lmo/pi_lmo3.cpp
        src/lmo/pi_lmo4.cpp
        src/lmo/pi_lmo5.cpp
        src/lmo/pi_lmo_parallel1.cpp)

# libprimesieve ######################################################

option(BUILD_PRIMESIEVE "Build primesieve binary" OFF)
set(BUILD_TESTS OFF CACHE BOOL "Build primesieve tests" FORCE)

add_subdirectory(lib/primesieve)

# Check for OpenMP ###################################################

if(WITH_OPENMP)
    find_package(OpenMP QUIET)
    if(OPENMP_FOUND)
        message(STATUS "Found OpenMP: ${OpenMP_CXX_FLAGS}")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    endif()
endif()

# Check -mpopcnt compiler flag #######################################

if(WITH_POPCNT)
    set(COPY_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
    set(CMAKE_CXX_FLAGS -Werror)
    check_cxx_compiler_flag(-mpopcnt mpopcnt)
    set(CMAKE_CXX_FLAGS "${COPY_CXX_FLAGS}")
    if(mpopcnt)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mpopcnt")
    endif()
endif()

# Check if int128_t type exists ######################################

check_type_size("int128_t" INT128_T)
if(HAVE_INT128_T)
    add_definitions(-DHAVE_INT128_T)
else()
    check_type_size("__int128_t" __INT128_T)
    if(HAVE___INT128_T)
        add_definitions(-DHAVE___INT128_T)
    else()
        message(FATAL_ERROR 
                "Your compiler does not support __int128_t!")
    endif()
endif()

# Check if libdivide.h works #########################################

if(WITH_LIBDIVIDE)
    message(STATUS "Performing Test libdivide")

    try_run(RUN_LIBDIVIDE COMPILE_LIBDIVIDE
            ${CMAKE_BINARY_DIR}
            ${CMAKE_CURRENT_SOURCE_DIR}/test/libdivide.cpp
            COMPILE_DEFINITIONS -I${CMAKE_CURRENT_SOURCE_DIR}/include)

    if(RUN_LIBDIVIDE EQUAL 0)
        message(STATUS "Performing Test libdivide - Success")
    else()
        message(STATUS "Performing Test libdivide - Failed")
    endif()
endif()

# Include S2_easy* source files ######################################

if(RUN_LIBDIVIDE EQUAL 0)
    set(SRC ${SRC} src/deleglise-rivat/S2_easy_libdivide.cpp)
else()
    set(SRC ${SRC} src/deleglise-rivat/S2_easy.cpp)
endif()

# Build binary #######################################################

add_executable(primesum ${SRC})
target_link_libraries(primesum libprimesieve)
install(TARGETS primesum DESTINATION bin)
