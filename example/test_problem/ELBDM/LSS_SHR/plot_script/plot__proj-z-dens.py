#!/usr/bin/env python3

import argparse
import sys

import pandas as pd
import yt

# load the command-line parameters
parser = argparse.ArgumentParser( description='Projection of mass density' )

parser.add_argument( '-i', action='store', required=False, type=str, dest='prefix',
                     help='path prefix [%(default)s]', default='../' )
parser.add_argument( '-s', action='store', required=True,  type=int, dest='idx_start',
                     help='first data index' )
parser.add_argument( '-e', action='store', required=True,  type=int, dest='idx_end',
                     help='last data index' )
parser.add_argument( '-d', action='store', required=False, type=int, dest='didx',
                     help='delta data index [%(default)d]', default=1 )
parser.add_argument( '-halo', action='store', required=False, type=int, dest='halo',
                     help='which halo [%(default)d]', default=1 )
parser.add_argument( '-axis', action='store', required=False, type=str, dest='axis',
                     help='projection axis [%(default)s]', default='z' )

args=parser.parse_args()

# take note
print( '\nCommand-line arguments:' )
print( '-------------------------------------------------------------------' )
for t in range( len(sys.argv) ):
   print (str(sys.argv[t]),end = ' ')
print( '' )
print( '-------------------------------------------------------------------\n' )


idx_start = args.idx_start
idx_end   = args.idx_end
didx      = args.didx
prefix    = args.prefix
halo      = args.halo
axis      = args.axis

field         = 'density'
colormap_dens = 'algae'
center_mode   = 'c'
dpi           = 150

yt.enable_parallelism()
ts = yt.DatasetSeries( [ prefix+'/Data_%06d'%idx for idx in range(idx_start, idx_end+1, didx)] )

for ds in ts.piter():

   num = '%s'%ds
   num = int(num[5:11])

   pz_dens = yt.ProjectionPlot( ds, axis, field, center=center_mode )
   pz_dens.set_axes_unit( 'code_length' )
   pz_dens.set_zlim( field, 1.0e-5, 4.0e-2 )
   pz_dens.set_cmap( field, colormap_dens )
   pz_dens.annotate_timestamp( time_unit='Gyr', corner='upper_right' )
   pz_dens.save( 'proj/Data_%06d_Projection_%s_density.png'%(num, axis),mpl_kwargs={"dpi":dpi} )

   # zoom in projection
   # df = pd.read_csv( 'Halo_Parameter_%d'%halo , sep = r'\s+' , header = 0 , index_col='#')
   # coordinate_x = df['x'][num]
   # coordinate_y = df['y'][num]
   # coordinate_z = df['z'][num]
   # # coordinate_x = 0.7956176757812499777955
   # # coordinate_y = 1.5774780273437496447286
   # # coordinate_z = 0.1320190429687499722444
   # pz_dens = yt.ProjectionPlot( ds, axis, field, center=[coordinate_x, coordinate_y, coordinate_z] )
   # pz_dens.zoom(8)
   # pz_dens.set_axes_unit( 'kpc' )
   # pz_dens.set_cmap( field, colormap_dens )
   # pz_dens.set_zlim( field, 1.0e-5, 4.0e-2 )
   # # pz_dens.annotate_grids()
   # pz_dens.save( 'proj/Data_%06d_%d_Projection_%s_density_zoom.png'%(num, halo, axis),mpl_kwargs={"dpi":dpi} )
