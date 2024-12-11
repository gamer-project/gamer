#!/bin/bash
# Create softlinks for interpolation tables

# GAMER expects interpolation tables to have filenames in the format "N=i.bin"
# where i is the size of the input vector to be interpolated
# Specify the directory where your original files are located
original_directory="."

# Specify the directory where you want to create the soft links
soft_links_directory="."

# Create soft links for the files with the specified naming convention
for file in N=*_m=*_nDelta=*_nd=*_Gamma=*_g=*.bin; do
    # Extract the value of N from the filename
    n_value=$(echo "$file" | sed -n 's/.*N=\([0-9]*\)_m.*/\1/p')

    # Create a soft link with the desired naming convention
    ln -f -s "$original_directory/$file" "$soft_links_directory/N=$n_value.bin"
done
