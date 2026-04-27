#!/bin/bash

cd "$(dirname "$0")"

# If $1 is given, use it as $chosen, else prompt the user
if [ -n "$1" ]; then
    chosen="$1"
else
    echo "Directories under bin/:"
    for dir in ../bin/*/; do
        echo "$(basename "$dir")"
    done
    echo -n "Enter a directory name from the above list: "
    read chosen
fi

printf "%s" "$chosen" > bin_working
echo "bin_working updated."