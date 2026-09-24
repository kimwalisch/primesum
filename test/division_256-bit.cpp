#include <iostream>
#include <cstdlib>
#include <random>
#include <limits>

#include <int256_t.hpp>
#include <boost/multiprecision/cpp_int.hpp>

using namespace primesum;

typedef boost::multiprecision::int256_t boost_int256_t;

int main(int argc, char** argv)
{
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_int_distribution<int64_t> dist(1, std::numeric_limits<int64_t>::max());

    std::random_device rd2;
    std::mt19937 gen2(rd2());
    std::uniform_int_distribution<int64_t> dist2(1, std::numeric_limits<int16_t>::max());

    int iters = 10000;

    if (argc > 1)
        iters = std::atoi(argv[1]);

    // Test random dividends and small quotients    
    for (int i = 0; i < iters; i++)
    {
        int256_t x1 = 1;
        boost_int256_t y1 = 1;

        for (int i = 0; i < 4; i++)
        {
            int64_t a = dist(gen);
            int256_t x2 = 1;
            boost_int256_t y2 = 1;

            x1 *= a; x1 += 1;
            y1 *= a; y1 += 1;

            for (int j = 1; j < 17; j++)
            {
                x2 = j;
                y2 = j;

                auto res1 = x1 / x2;
                auto res2 = y1 / y2;

                // std::cout << y1 << " / " << y2 << "\n";

                while (res2 != 0)
                {
                    if ((uint64_t) res1 != (uint64_t) res2)
                    {
                        std::cerr << x1 << std::endl;
                        std::cerr << x2 << std::endl;
                        std::cerr << y1 << std::endl;
                        std::cerr << y2 << std::endl;
                        std::cerr << (x1 / x2) << " != " << (y1 / y2) << std::endl;
                        return 1;
                    }

                    res1 >>= 64;
                    res2 >>= 64;
                }
            }
        }
    }

    // Test random dividends and quotients
    for (int i = 0; i < iters; i++)
    {
        int256_t x1 = 1;
        boost_int256_t y1 = 1;

        for (int i = 0; i < 4; i++)
        {
            int64_t a = dist(gen);
            int256_t x2 = 1;
            boost_int256_t y2 = 1;

            x1 *= a; x1 += 1;
            y1 *= a; y1 += 1;

            for (int j = 0; j < 17; j++)
            {
                int64_t b = dist2(gen);

                x2 *= b; x2 += 1;
                y2 *= b; y2 += 1;

                auto res1 = x1 / x2;
                auto res2 = y1 / y2;

                // std::cout << y1 << " / " << y2 << "\n";

                while (res2 != 0)
                {
                    if ((uint64_t) res1 != (uint64_t) res2)
                    {
                        std::cerr << x1 << std::endl;
                        std::cerr << x2 << std::endl;
                        std::cerr << y1 << std::endl;
                        std::cerr << y2 << std::endl;
                        std::cerr << (x1 / x2) << " != " << (y1 / y2) << std::endl;
                        return 1;
                    }

                    res1 >>= 64;
                    res2 >>= 64;
                }
            }
        }
    }

    return 0;
}
