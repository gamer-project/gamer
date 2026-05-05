#!/usr/bin/env python3.9
################################################################################

# This code outputs the density, cumulative mass, circular velocity and velocity dispersion profiles
# of a target halo. The code uses an initial estimate (center_first_guess) to search the vicinity within
# a specified radius (vicinity) for the coordinates of the maximum density,
# corresponding to the center of the target halo.

# carefully adjust vicinity for performance.

# Example usage after running simulation with low resolution IC:
#    python3 Compute_profiles.py -s 36 -e 36

# Outputs:
#   Halo_Parameter
#   prof_dens/Data_xxxxxx_0_profile_data
#   prof_mass/Data_xxxxxx_0_mass_accumule
#   prof_dens/Data_xxxxxx_0_profile_data
#   prof_circular_vel/Data_xxxxxx_0_circular_velocity
#   prof_veldisp/Data_xxxxxx_0_veldisp_haloRestFrame

################################################################################
import argparse
import numpy as np
import yt
import sys
from Profile_Functions import *


# load the command-line parameters
parser = argparse.ArgumentParser( description='Plot profile and output Halo_parameter' )

parser.add_argument( '-s', action='store', required=True,  type=int, dest='idx_start',
                      help='first data index' )
parser.add_argument( '-e', action='store', required=True,  type=int, dest='idx_end',
                      help='last data index' )
parser.add_argument( '-d', action='store', required=False, type=int, dest='didx',
                      help='delta data index [%(default)d]', default=1 )

# take note
print( '\nCommand-line arguments:' )
print( '-------------------------------------------------------------------' )
print( ' '.join(map(str, sys.argv)) )
print( '-------------------------------------------------------------------\n' )

args        = parser.parse_args()
idx_start   = args.idx_start
idx_end     = args.idx_end
didx        = args.didx

halo_id      = 0   # pick an id for your target halo

for file_id in range( idx_start, idx_end+1, didx ):

    ds                 = yt.load( "../Data_%06d"%file_id )
    center_first_guess = np.array( [0.295, 9.522, 8.27] ) # in cMpc/h. First guess for target halo of low resolution IC at z=0
    vicinity           = 0.3                              # radius in cMpc/h
    compute_profile( ds, center_first_guess, vicinity, halo_id, '.' )

print( "Done !" )
