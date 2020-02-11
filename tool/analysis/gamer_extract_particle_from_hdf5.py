from __future__ import print_function
import argparse
import sys
import h5py
import numpy as np
import os.path


# load the command-line parameters
parser = argparse.ArgumentParser( description='Extract particle data from a GAMER HDF5 output' )

parser.add_argument( '-i', action='store', required=True, type=str, dest='filename_in',
                     help='input filename' )
parser.add_argument( '-o', action='store', required=True, type=str, dest='filename_out',
                     help='output filename' )
parser.add_argument( '-t', '--text', action='store_true', dest='output_text',
                     help='output text instead of binary file [False]' )
parser.add_argument( '-d', '--double', action='store_true', dest='float64',
                     help='use double precision [False]' )

args=parser.parse_args()

# take note
print( '\nCommand-line arguments:' )
print( '-------------------------------------------------------------------' )
for t in range( len(sys.argv) ):
   print( str(sys.argv[t]), end=' ' )
print( '' )
print( '-------------------------------------------------------------------\n' )

# short names
filename_in  = args.filename_in
filename_out = args.filename_out
output_text  = args.output_text
float_type   = 'float64' if args.float64 else 'float32'


# open and check files
handle_in = h5py.File( filename_in, 'r' )

# check whether the input file has particle data
assert 'Particle' in handle_in, 'No particle data can be found'

# check whether the outfile file already exists
assert not os.path.isfile( filename_out ), 'output file \"'+filename_out+'\" already exists!'

handle_out = open( filename_out, "a" )


# get the particle attribute list
par_group = handle_in['Particle']
att_list = [ v for v in par_group.keys() ]


# output simulation info
time = handle_in['Info']['KeyInfo']['Time'][0]
step = handle_in['Info']['KeyInfo']['Step']
npar = handle_in['Info']['KeyInfo']['Par_NPar']
print ( "%-19s : %20.14e"%("Time", time) )
print ( "%-19s : %ld"    %("Step", step) )
print ( "%-19s : %ld"    %("# of particles", npar) )
print ( "%-19s :"        %("Particle attributes"), end=' ' ),
for v in att_list: print( "%s"%v, end=' ' ),
print( '' )


# output text file
if output_text:
#  load all particles
   par_data = {}
   for v in att_list:
      par_data[v] = np.asarray( par_group[v], dtype=float_type )

#  output header
   handle_out.write( "#Time %20.14e   Step %13ld   Active Particles %13ld\n\n" % (time, step, npar) )
   handle_out.write( "#" )
   for i, v in enumerate( att_list ):
      handle_out.write( "  %*s"%(20 if i==0 else 21, v) );
   handle_out.write( "\n" )

#  output all attributes of a single particle at a time
   for p in range( npar ):
      for v in att_list:
         handle_out.write( "  %21.14e"%par_data[v][p] );
      handle_out.write( "\n" )


# output binary file
else:
#  output one attribute at a time
   for v in att_list:
      par_data = np.asarray( par_group[v], dtype=float_type )
      par_data.tofile( handle_out )


# close files
handle_out.close()
handle_in.close()
