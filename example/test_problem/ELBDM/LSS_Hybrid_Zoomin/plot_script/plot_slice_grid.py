################################################################################

# This code output the density slice of the simulation with a pre-determined center.

################################################################################

import argparse
import sys
import yt
import numpy as np

# global parameter settings
field         = ('gas', 'density')
axes          = ["x", "y", "z"]                            # directions
lv            = 10                                         # maximum level for sampling AMR grid
dpi           = 300
colormap_dens = 'algae'
center_mode   = 'c'
zooms         = [1,7]

# load the command-line parameters
parser = argparse.ArgumentParser(description='Slice of mass density')

parser.add_argument( '-i', action='store', required=False, type=str, dest='prefix',
                     help='path prefix [%(default)s]', default='../' )
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

args      = parser.parse_args()
idx_start = args.idx_start
idx_end   = args.idx_end
didx      = args.didx
prefix    = args.prefix

center = [0.25817871, 9.3939209,  8.20983887] # center of the target halo at z=0 for low resolution IC

yt.enable_parallelism()
ts = yt.DatasetSeries( [prefix+'/Data_%06d'%idx for idx in range( idx_start, idx_end+1, didx )] )

for ds in ts.piter():
   num = '%s'%ds
   num = int( num[5:11] )
   for ax in axes:
      for zoom in zooms:
         pz_dens = yt.SlicePlot( ds, ax, field, center=center )

         pz_dens.set_cmap( field, colormap_dens )
         pz_dens.annotate_timestamp( time_unit='Gyr', redshift=True, corner='upper_right' )
         pz_dens.set_axes_unit( 'Mpc/h' )
         pz_dens.zoom( zoom )
         pz_dens.save( 'Data_%06d_Slice_%s_%s_x%d.png'%(num, ax, field[1],zoom), mpl_kwargs={"dpi":dpi} )

         pz_dens.annotate_grids()
         pz_dens.save( 'Data_%06d_Slice_%s_%s_x%d_grid.png'%(num, ax, field[1],zoom), mpl_kwargs={"dpi":dpi} )
         pz_dens.zoom( 1./zoom )
