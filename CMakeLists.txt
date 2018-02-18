cmake_minimum_required(VERSION 3.4)
project(primesum CXX)
set(PRIMESUM_VERSION_MAJOR 1)
set(PRIMESUM_VERSION_MINOR 4)
set(PRIMESUM_VERSION "${PRIMESUM_VERSION_MAJOR}.${PRIMESUM_VERSION_MINOR}")

# Build options ######################################################

option(WITH_POPCNT       "Enable x86 POPCNT instruction"  ON)
option(WITH_LIBDIVIDE    "Use libdivide.h"                ON)
option(WITH_OPENMP       "Enable OpenMP support"          ON)
option(BUILD_STATIC_LIBS "Build static libprimesieve"     ON)
option(BUILD_SHARED_LIBS "Build shared libprimesieve"     OFF)
option(BUILD_PRIMESIEVE  "Build primesieve binary"        OFF)

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

# Include Check* #####################################################

include(CheckTypeSize)
include(CheckCXXCompilerFlag)
include(CheckCXXSourceCompiles)
include(CMakePushCheckState)
include(GNUInstallDirs)

# Set default build type to Release ##################################

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING
        "Choose the type of build, options are: None Debug Release RelWithDebInfo MinSizeRel."
        FORCE)
endif()

# Check if int128_t type exists ######################################

check_type_size("__int128_t" __INT128_T LANGUAGE CXX)
if(NOT HAVE___INT128_T)
    check_type_size("int128_t" INT128_T LANGUAGE CXX)
    if(NOT HAVE_INT128_T)
        message(FATAL_ERROR "Compiler does not support int128_t!")
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
target_link_libraries(primesum PRIVATE libprimesieve)
target_include_directories(primesum PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
install(TARGETS primesum DESTINATION ${CMAKE_INSTALL_BINDIR})

target_compile_features(primesum
PRIVATE
    cxx_auto_type
    cxx_lambdas
    cxx_nullptr
    cxx_range_for
    cxx_uniform_initialization)

# Check for OpenMP ###################################################

if(WITH_OPENMP)
    find_package(OpenMP QUIET)

    if(OPENMP_FOUND)
        message(STATUS "Found OpenMP: ${OpenMP_CXX_FLAGS}")

        target_compile_options(primesum PRIVATE "${OpenMP_CXX_FLAGS}")
        target_link_libraries(primesum PRIVATE "${OpenMP_CXX_FLAGS}")

        cmake_push_check_state()
        set(CMAKE_REQUIRED_INCLUDES "${CMAKE_CURRENT_SOURCE_DIR}/include")
        set(CMAKE_REQUIRED_FLAGS "${CMAKE_CXX11_STANDARD_COMPILE_OPTION}")

        check_cxx_source_compiles("
            #include <omp.h>
            #include <int256_t.hpp>
            int main() {
                primesum::int256_t sum;
                #pragma omp parallel for reduction(+: sum)
                for (int i = 0; i < 100; i++)
                    sum += i;
                return 0;
            }"
            omp_declare_reduction)

        if(NOT omp_declare_reduction)
            message(FATAL_ERROR "Compiler does not support OpenMP 4.0 or later!")
        endif()

        cmake_pop_check_state()
    endif()
endif()

# Check -mpopcnt compiler flag #######################################

if(WITH_POPCNT)
    cmake_push_check_state()
    set(CMAKE_REQUIRED_FLAGS -Werror)
    check_cxx_compiler_flag(-mpopcnt mpopcnt)
    cmake_pop_check_state()

    if(mpopcnt)
        target_compile_options(primesum PRIVATE -mpopcnt)
    endif()
endif()

# libprimesieve ######################################################

add_subdirectory(lib/primesieve)
