import argparse
import os, h5py, subprocess, re, sys
import numpy as np

# load the command-line parameters
parser = argparse.ArgumentParser( description='Convert GAMER ELBDM wave (RE/IM) IC to hybrid (DENS/PHASE) IC' )

parser.add_argument( '-float8', action='store', required=False,  type=bool, dest='float8',
                            help='assume double precision for input and output file', default=False )
parser.add_argument( '-N', action='store', required=True, type=int, dest='N',
                            help='input file resolution')
parser.add_argument( '-input', action='store', required=True,  type=str, dest='input',
                      help='input file' )
parser.add_argument( '-output', action='store', required=True,  type=str, dest='output',
                      help='output file' )

print("Convert ELBDM wave IC to hybrid IC!")
print("Note that this conversion is only well-defined as long as the IC do not contain vortices")
print("If IC contain vortices, the phase field is not globally smooth and the conversion might not work as expected")

###
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

