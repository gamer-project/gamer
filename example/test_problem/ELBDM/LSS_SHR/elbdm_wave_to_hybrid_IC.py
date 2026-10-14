#!/usr/bin/env python

import os
import h5py
import subprocess
import re
import sys
import numpy as np
import argparse

# Load the command-line parameters
parser = argparse.ArgumentParser(description='Convert GAMER ELBDM wave (RE/IM) IC to hybrid (DENS/PHASE) IC. \n \
                                 Example usage: python elbdm_wave_to_hybrid_IC.py -resolution 256 -input UM_IC_wave -output UM_IC_hybrid.\n\
                                 Conversion is well-defined only if initial conditions do not contain vortices, i.e. at high redshift for cosmological IC.')

# Define command-line arguments
parser.add_argument('-float8', action='store', required=False, type=bool, dest='float8',
                    help='assume double precision for input and output file', default=False)
parser.add_argument('-resolution', action='store', required=True, type=int, dest='resolution',
                    help='input file resolution')
parser.add_argument('-input', action='store', required=True, type=str, dest='input',
                    help='input file')
parser.add_argument('-output', action='store', required=True, type=str, dest='output',
                    help='output file')

args = parser.parse_args()

# Print the command-line arguments for reference
print('\nCommand-line arguments:')
print('-------------------------------------------------------------------')
for t in range(len(sys.argv)):
    print(str(sys.argv[t]), end=' ')
print('')
print('-------------------------------------------------------------------\n')

# Set up interpolation parameters
input_file  = args.input
output_file = args.output
resolution  = args.resolution
float8      = args.float8  # Enable double precision

if float8:
    dtype = np.double
else:
    dtype = np.single

print("Reading wave IC from input file...")
binary = np.fromfile(input_file, dtype=dtype)
binary = binary.reshape((2, resolution, resolution, resolution))

re = binary[0, :, :, :]
im = binary[1, :, :, :]
ph = np.arctan2(im, re)
de = re**2 + im**2

# Unwrap phase
ph = np.unwrap(ph, axis=0)
ph = np.unwrap(ph, axis=1)
ph = np.unwrap(ph, axis=2)
binary[0, :, :, :] = de
binary[1, :, :, :] = ph

print("Writing hybrid IC to output file...")
with open(output_file, "wb") as f:
    binary.tofile(f)
    f.close()
