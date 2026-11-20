#!/bin/bash

if [ "$#" -ne 2 ]; then
    echo "Usage: $0 <suffix> <number_of_tables>"
    echo "Example: $0 KroupaIMF_ParticleMass32000_RandomSeed00123_Z0.02 443"
    exit 1
fi

SUFFIX=$1
NUM_TABLES=$2
DIR="Tables__${SUFFIX}"

if [ ! -d "$DIR" ]; then
    echo "Error: Source directory '$DIR' not found."
    exit 1
fi

mkdir -p "FeedbackTables"
for (( i=1; i<=NUM_TABLES; i++ )); do
    src_num=$(printf "%04d" $i)
    link_num=$(printf "%06d" $i)

    src_file="../${DIR}/Table_CCSNe_with_deterministic_fixed_count_${src_num}__${SUFFIX}.txt"
    link_name="FeedbackTables/FeedbackTable_Resolved_SNeII_N${link_num}"

    ln -s "$src_file" "$link_name"
done

echo "Created $NUM_TABLES symlinks for suffix: $SUFFIX"
