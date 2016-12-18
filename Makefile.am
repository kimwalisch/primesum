AM_CPPFLAGS = -I$(top_srcdir)/include

lib_LTLIBRARIES = libprimesum.la
libprimesum_la_CXXFLAGS = $(OPENMP_CXXFLAGS) $(POPCNT_FLAGS)
libprimesum_la_LDFLAGS = -version-info @primesum_lib_version@

if IS_WINDOWS
libprimesum_la_LDFLAGS += -no-undefined 
endif

bin_PROGRAMS = primesum
primesum_LDADD = libprimesum.la

# If installation directory is /usr/local/lib
# execute ldconfig after installation
install-exec-hook:
	test "$(prefix)" = "/usr/local" && command -v ldconfig 2>/dev/null && ldconfig /usr/local/lib 2>/dev/null || true
	test -f primesieve-*/primesieve && cd primesieve-* && make install && cd .. || true

uninstall-hook:
	test -f primesieve-*/primesieve && cd primesieve-* && make uninstall && cd .. || true

check: all
	./primesum --test

include_HEADERS = \
	include/primesum.hpp

EXTRA_DIST = \
	autogen.sh \
	build.sh \
	doc/alpha-tuning-factor.nb \
	doc/alpha-tuning-factor.pdf \
	scripts/find_fastest_alpha.sh \
	scripts/worktodo.sh \
	README.md

primesum_SOURCES = \
	src/app/main.cpp \
	src/app/help.cpp \
	src/app/cmdoptions.cpp \
	src/app/cmdoptions.hpp

libprimesum_la_SOURCES = \
	src/BitSieve.cpp \
	src/FactorTable.cpp \
	src/generate.cpp \
	src/nth_prime.cpp \
	src/P2.cpp \
	src/PhiTiny.cpp \
	src/PiTable.cpp \
	src/phi.cpp \
	src/phi_sum.cpp \
	src/pi_primesieve.cpp \
	src/pi_legendre.cpp \
	src/primesum.cpp \
	src/print.cpp \
	src/S1.cpp \
	src/S2LoadBalancer.cpp \
	src/S2Status.cpp \
	src/test.cpp \
	src/Wheel.cpp \
	src/deleglise-rivat/S2_trivial.cpp \
	src/deleglise-rivat/S2_hard.cpp \
	src/deleglise-rivat/pi_deleglise_rivat_parallel1.cpp \
	src/lmo/pi_lmo1.cpp \
	src/lmo/pi_lmo2.cpp \
	src/lmo/pi_lmo3.cpp \
	src/lmo/pi_lmo4.cpp \
	src/lmo/pi_lmo5.cpp \
	src/lmo/pi_lmo_parallel1.cpp \
	include/aligned_vector.hpp \
	include/bit_scan_forward.hpp \
	include/BitSieve.hpp \
	include/calculator.hpp \
	include/FactorTable.hpp \
	include/fast_div.hpp \
	include/generate.hpp \
	include/int128_t.hpp \
	include/isqrt.hpp \
	include/min_max.hpp \
	include/popcount.hpp \
	include/imath.hpp \
	include/print.hpp \
	include/primesum-internal.hpp \
	include/PhiTiny.hpp \
	include/PiTable.hpp \
	include/S1.hpp \
	include/S2.hpp \
	include/S2LoadBalancer.hpp \
	include/S2Status.hpp \
	include/SumBits.hpp \
	include/tos_sums.hpp \
	include/Wheel.hpp

if HAVE_LIBDIVIDE

libprimesum_la_SOURCES += \
	src/deleglise-rivat/S2_easy_libdivide.cpp \
	include/libdivide.h

else

libprimesum_la_SOURCES += \
	src/deleglise-rivat/S2_easy.cpp

endif
