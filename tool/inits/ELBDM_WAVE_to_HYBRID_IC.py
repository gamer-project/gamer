#!/usr/bin/env python

import os, h5py, subprocess, re, sys
import numpy as np
import argparse

# load the command-line parameters
parser = argparse.ArgumentParser( description='Convert GAMER ELBDM wave (RE/IM) IC to hybrid (DENS/PHASE) IC. Example usage: python ELBDM_WAVE_to_HYBRID_IC.py -resolution 256 -input UM_IC_wave -output UM_IC_hybrid. Conversion only well-defined if IC do not contain vortices, i.e. at high redshift for cosmological IC.')

parser.add_argument( '-float8', action='store', required=False,  type=bool, dest='float8',
                            help='assume double precision for input and output file', default=False )
parser.add_argument( '-resolution', action='store', required=True, type=int, dest='N',
                            help='input file resolution')
parser.add_argument( '-input', action='store', required=True,  type=str, dest='input',
                      help='input file' )
parser.add_argument( '-output', action='store', required=True,  type=str, dest='output',
                      help='output file' )

args=parser.parse_args()

# take note
print( '\nCommand-line arguments:' )
print( '-------------------------------------------------------------------' )
for t in range( len(sys.argv) ):
    print( str(sys.argv[t]) ),
print( '' )
print( '-------------------------------------------------------------------\n' )

# set up interpolation parameters
input_file  = args.input
output_file = args.output
N           = args.N
FLOAT_8     = args.float8           # enable double precision
dtype		= np.single
if FLOAT_8:
   dtype = np.double

print("Reading wave IC from input file...")
binary = np.fromfile(input_file, dtype=dtype)
binary = binary.reshape((2, N, N, N))

re = binary[0, :, :, : ]
im = binary[1, :, :, : ]
ph = np.arctan2(im, re)
de = re**2 + im**2

# unwrap phase
ph = np.unwrap(ph, axis=0)
ph = np.unwrap(ph, axis=1)
ph = np.unwrap(ph, axis=2)
binary[0, :, :, :] = de
binary[1, :, :, :] = ph

print("Writing hybrid IC to output file...")
with open(output_file,"wb") as f:
     binary.tofile(f)
     f.close()

