#!/bin/bash

# Usage: ./worktodo.sh
# Iterates over the numbers in worktodo.txt (one number per line with
# optional flags) and processes them using primesum.

command -v ./primesum >/dev/null 2>/dev/null
if [ $? -ne 0 ]
then
    echo "Error: no primesum binary in current directory."
    exit 1
fi

# check if primesum supports --log option
log_flag=$(./primesum --help | grep Log > /dev/null && echo --log)

while read first_line < worktodo.txt
do
    # Skip empty lines and comments
    if [ "$first_line" != "" ] && [ ${first_line:0:1} != "#" ]
    then
        ./primesum $first_line $log_flag

        # Check if primesum exited successfully
        if [ $? -ne 0 ]
        then
            echo ""
            echo "Error in primesum_from_file.sh:"
            echo "The following command failed: ./primesum $first_line"
            exit 1
        fi
    fi

    # delete first line from worktodo.txt
    tail -n +2 worktodo.txt > .tmp_worktodo.txt
    mv -f .tmp_worktodo.txt worktodo.txt
done
