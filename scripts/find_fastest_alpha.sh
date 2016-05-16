#!/bin/bash

# Find the fastest alpha tuning factors for primesum.
# Usage:
#   ./find_fastest_alpha.sh [--start=n] [--stop=n] [-t=threads]
# Description:
#   This script calculates pi(10^n) for start <= n <= stop using different
#   alpha tuning factors and prints out the fastest alphas. The fastest
#   alpha factors will vary slightly for different CPUs. Once a list of
#   fast alpha factors has been generated you can follow the
#   instructions in doc/alpha-factor-tuning.pdf and generate a new
#   alpha tuning function for use in primesum's source code.

# Execute in base directory
if [ "$(basename $(pwd))" = "scripts" ]
then
    cd ..
fi

command -v ./primesum >/dev/null 2>/dev/null
if [ $? -ne 0 ]
then
    echo "Error: no primesum binary in current directory."
    exit 1
fi

command -v bc >/dev/null 2>/dev/null
if [ $? -ne 0 ]
then
    echo "Error: GNU bc is not installed."
    exit 1
fi

start=1
stop=21
seconds=0
repeat=3
threads=$(./primesum 1e18 --S2_trivial -s | grep threads | cut -d'=' -f2 | cut -d' ' -f2)

for i in "$@"
do
    case $i in
        -t=*)
        threads="${i#*=}"
        shift
        ;;
        --start=*)
        start="${i#*=}"
        shift
        ;;
        --stop=*)
        stop="${i#*=}"
        shift
        ;;
        *)
        echo "Find the fastest alpha tuning factors for primesum."
        echo "Usage:"
        echo "  ./find_fastest_alpha.sh [--start=n] [--stop=n] [-t=threads]"
        echo "Description:"
        echo "  This script calculates pi(10^n) for start <= n <= stop using different"
        echo "  alpha tuning factors and prints out the fastest alphas. The fastest"
        echo "  alpha factors will vary slightly for different CPUs. Once a list of"
        echo "  fast alpha factors has been generated you can follow the"
        echo "  instructions in doc/alpha-factor-tuning.pdf and generate a new"
        echo "  alpha tuning function for use in primesum's source code."
        exit 1
        ;;
    esac
done

# Returns 1 if $1 == $2, else 0
function is_equal
{
    echo $1'=='$2 | bc -l
}

# Returns 1 if $1 < $2, else 0
function is_smaller
{
    echo $1'<'$2 | bc -l
}

# Returns 1 if $1 <= $2, else 0
function is_smaller_equal
{
    echo $1'<='$2 | bc -l
}

# Returns 1 if $1 > $2, else 0
function is_greater
{
    echo $1'>'$2 | bc -l
}

# Floating point calculator
# $1: String containing an arithmetic expression
function calc
{
    echo "scale=3; $1" | bc -l
}

# Given 2 numbers, returns the greatest number 
function maximum
{
    if [ $(is_greater $1 $2) -eq 1 ]
    then
        echo $1
    else
        echo $2
    fi
}

# Given 2 numbers, returns the greatest number 
function minimum
{
    if [ $(is_smaller $1 $2) -eq 1 ]
    then
        echo $1
    else
        echo $2
    fi
}

# $1: primesum args
function get_primesum_alpha
{
    alpha=$(./primesum $1 --S2_trivial -s | grep alpha | cut -d'=' -f2 | cut -d' ' -f2)
    echo $alpha
}

# $1: primesum args
function get_primesum_seconds
{
    sleep 0.2
    seconds=$(./primesum $1 --time | grep Seconds | cut -d':' -f2 | cut -d' ' -f2)
    echo $seconds
}

# Calculate pi(10^i) for start <= i <= stop
for ((i = start; i <= stop; i++))
do
    alpha=$(get_primesum_alpha "1e$i")
    fastest_seconds=$(calc "10^9")
    fastest_alpha="1.000"

    # This loops tries to get close (< 25%) to the fastest alpha
    for ((j = 0; j < repeat; j++))
    do
        pivot=$(calc "$alpha / 2")
        max_alpha=$(calc "$alpha + $pivot")
        new_alpha=$(calc "$alpha - $pivot")
        new_alpha=$(maximum 1.000 $new_alpha)
        increment=$(calc "($max_alpha - $new_alpha) / 4")
        increment=$(maximum 0.1 $increment)

        while [ $(is_smaller_equal $new_alpha $max_alpha) -eq 1 ]
        do
            seconds=$(get_primesum_seconds "1e$i -t$threads -a$new_alpha")

            if [ $(is_smaller $seconds $fastest_seconds) -eq 1 ]
            then
                fastest_alpha=$new_alpha
                fastest_seconds=$seconds
            fi

            # Reduce the number of long running primesum benchmarks
            if [ $(is_greater $seconds 20) -eq 1 ]
            then
                repeat=1
            elif [ $(is_greater $seconds 3) -eq 1 ]
            then
                repeat=2
            fi

            # Benchmark runs too quickly for this small input
            if [ $(is_equal $seconds 0) -eq 1 ]
            then
                break
            fi

            new_alpha=$(calc "$new_alpha + $increment")
        done
    done

    if [ $(is_greater $fastest_seconds 0) -eq 1 ]
    then
        # This loops tries to get very close (< 1.6%) to the fastest alpha
        for ((j = 0; j < repeat; j++))
        do
            pivot=$(calc "$fastest_alpha / 8")
            max_alpha=$(calc "$fastest_alpha + $pivot")
            new_alpha=$(calc "$fastest_alpha - $pivot")
            new_alpha=$(maximum 1.000 $new_alpha)
            increment=$(calc "($max_alpha - $new_alpha) / 8")
            increment=$(maximum 0.1 $increment)

            while [ $(is_smaller_equal $new_alpha $max_alpha) -eq 1 ]
            do
                seconds=$(get_primesum_seconds "1e$i -t$threads -a$new_alpha")

                if [ $(is_smaller $seconds $fastest_seconds) -eq 1 ]
                then
                    fastest_alpha=$new_alpha
                    fastest_seconds=$seconds
                fi

                new_alpha=$(calc "$new_alpha + $increment")
            done
        done
    fi

    # When fastest_alpha=1.000 the benchmark ran too quickly to get
    # any meaningful result, so we report undef
    if [ $(is_equal $fastest_alpha 1) -eq 1 ] || \
       [ $(is_equal $fastest_seconds 0) -eq 1 ]
    then
        fastest_alpha="undef"
    fi

    # Print fastest alpha found for pi(10^$i)
    printf '%-11s %-12s %-18s %-22s %-20s\n' "pi(10^$i)" "threads=$threads" "old_alpha=$alpha" "fastest_alpha=$fastest_alpha" "seconds=$fastest_seconds"
done
