include(CheckCXXSourceCompiles)
include(CMakePushCheckState)

cmake_push_check_state()

# Our <int128_t.hpp> requires C++11 or later
if(NOT compiler_supports_cpp11)
    if(CMAKE_CXX11_EXTENSION_COMPILE_OPTION)
        set(CMAKE_REQUIRED_FLAGS "${CMAKE_CXX11_EXTENSION_COMPILE_OPTION}")
    elseif(CMAKE_CXX11_STANDARD_COMPILE_OPTION)
        set(CMAKE_REQUIRED_FLAGS "${CMAKE_CXX11_STANDARD_COMPILE_OPTION}")
    endif()
endif()

set(CMAKE_REQUIRED_INCLUDES "${PROJECT_SOURCE_DIR}/include")

check_cxx_source_compiles("
    #include <int128_t.hpp>
    #include <algorithm>

    using namespace primesum;

    int main() {
        int128_t x = int128_t(1) << 100;
        int128_t y = 1000;
        unsigned divider = 123;
        x = x / divider;
        x /= y;
        x = x + y;

        if (std::min(x, y) != y)
            return 1;

        long double z = ((long double) y) + 1.59867;
        int128_t iz = (int128_t) z;

        if (std::max(y, iz) != iz)
            return 1;

        return 0;
    }" int128)

cmake_pop_check_state()

if(NOT int128)
    message(FATAL_ERROR "Compiler does not support int128_t!")
endif()
