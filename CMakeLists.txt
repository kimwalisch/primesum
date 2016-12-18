cmake_minimum_required(VERSION 2.8.4)
project(primesum)
set(PRIMESUM_VERSION_MAJOR 1)
set(PRIMESUM_VERSION_MINOR 1)
set(PRIMESUM_VERSION "${PRIMESUM_VERSION_MAJOR}.${PRIMESUM_VERSION_MINOR}")
include_directories(include src/primesieve/include)
set(CMAKE_BUILD_TYPE Release)

# Build options ######################################################

option(ENABLE_POPCNT "Enable POPCNT instruction" ON)

# Include Check* #####################################################

include(CheckCXXCompilerFlag)
include(CheckCXXSourceRuns)
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

# primesieve library source files ####################################

set(PRIMESIEVE_PATH src/primesieve/src/primesieve)

set(SRC ${SRC}
        ${PRIMESIEVE_PATH}/EratBig.cpp
        ${PRIMESIEVE_PATH}/EratMedium.cpp
        ${PRIMESIEVE_PATH}/EratSmall.cpp
        ${PRIMESIEVE_PATH}/ParallelPrimeSieve.cpp
        ${PRIMESIEVE_PATH}/PreSieve.cpp
        ${PRIMESIEVE_PATH}/PrimeFinder.cpp
        ${PRIMESIEVE_PATH}/PrimeGenerator.cpp
        ${PRIMESIEVE_PATH}/PrimeSieve-nthPrime.cpp
        ${PRIMESIEVE_PATH}/PrimeSieve.cpp
        ${PRIMESIEVE_PATH}/SieveOfEratosthenes.cpp
        ${PRIMESIEVE_PATH}/WheelFactorization.cpp
        ${PRIMESIEVE_PATH}/iterator.cpp
        ${PRIMESIEVE_PATH}/popcount.cpp
        ${PRIMESIEVE_PATH}/primesieve-api.cpp
        ${PRIMESIEVE_PATH}/primesieve_test.cpp)

# Check for OpenMP ###################################################

find_package(OpenMP QUIET)
if (OPENMP_FOUND)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

# Check for -mpopcnt compiler option #################################

if(ENABLE_POPCNT)
    CHECK_CXX_COMPILER_FLAG(-mpopcnt MPOPCNT_FLAG)
    if(MPOPCNT_FLAG)
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

# Check if libdivide works ###########################################

set(CMAKE_REQUIRED_INCLUDES ${CMAKE_SOURCE_DIR}/include)

check_cxx_source_runs("
    #include <libdivide.h>
    #include <stdint.h>
    #include <limits>

    int main()
    {
        uint64_t tests[20] = {
            0, 1, 2, 3, 4, 5, 6, 7, 8, 9,
            63, 101, 511,
            1 << 5, 1 << 9, 1 << 20,
            std::numeric_limits<int32_t>::max(),
            std::numeric_limits<uint32_t>::max(),
            std::numeric_limits<int64_t>::max(),
            std::numeric_limits<uint64_t>::max() };

        for (int i = 0; i < 20; i++)
        {
            uint64_t x = tests[i];

            for (int j = 2; j < 1000; j++)
            {
                libdivide::divider<uint64_t, libdivide::BRANCHFREE> divider(j);
                if (x / j != x / divider)
                    return 1;
            }
        }

        return 0;
    }"
    HAVE_LIBDIVIDE)

# Include S2_easy* source files ######################################

if(HAVE_LIBDIVIDE)
    set(SRC ${SRC} src/deleglise-rivat/S2_easy_libdivide.cpp)
else()
    set(SRC ${SRC} src/deleglise-rivat/S2_easy.cpp)
endif()

# Build binary #######################################################

add_executable(primesum ${SRC})
install(TARGETS primesum DESTINATION bin)

