import argparse
import sys

import pandas as pd
import yt

# load the command-line parameters
parser = argparse.ArgumentParser( description='Slice of mass density' )

parser.add_argument( '-p', action='store', required=False, type=str, dest='prefix',
                     help='prefix [%(default)s]', default='../' )
parser.add_argument( '-s', action='store', required=True,  type=int, dest='idx_start',
                     help='first data index' )
parser.add_argument( '-e', action='store', required=True,  type=int, dest='idx_end',
                     help='last data index' )
parser.add_argument( '-d', action='store', required=False, type=int, dest='didx',
                     help='delta data index [%(default)d]', default=1 )
parser.add_argument( '-halo', action='store', required=False, type=int, dest='halo',
                     help='which halo [%(default)d]', default=1 )


args=parser.parse_args()

# take note
print( '\nCommand-line arguments:' )
print( '-------------------------------------------------------------------' )
for t in range( len(sys.argv) ):
   print(sys.argv[t], end = ' '),
print( '' )
print( '-------------------------------------------------------------------' )


idx_start   = args.idx_start
idx_end     = args.idx_end
didx        = args.didx
prefix      = args.prefix
halo        = args.halo

field       = 'density'
colormap    = 'algae'
center_mode = 'c'
dpi         = 150


yt.enable_parallelism()


df = pd.read_csv( 'Halo_Parameter_%d'%halo , sep = '\s+' , header = 0 , index_col='#')

ts = yt.DatasetSeries( [ prefix+'/Data_%06d'%idx for idx in range(idx_start, idx_end+1, didx) ] )

for ds in ts.piter():

   num = '%s'%ds
   num = int(num[9:11])

   coordinate_x = df['x'][num]
   coordinate_y = df['y'][num]
   coordinate_z = df['z'][num]
   time_z = df['time'][num]
   time_a = 1/(1+time_z)
   halo_radius_max = df['halo_radius'].iloc[-1]/time_a

   sz_dens = yt.SlicePlot( ds, 'z', field, center=[coordinate_x, coordinate_y, coordinate_z] )
   sz_dens.zoom(ds.domain_width.in_units("kpc")[0].d/halo_radius_max/3)
   sz_dens.set_axes_unit( 'kpc' )
   sz_dens.annotate_sphere([coordinate_x, coordinate_y, coordinate_z], radius=(df['halo_radius'][num]/time_a, "kpc"))
   sz_dens.set_zlim( field, 1e-28 ,5.0e-22 )
   sz_dens.set_cmap( field, colormap )
   sz_dens.annotate_timestamp( time_unit='Gyr', corner='upper_right' )
   sz_dens.annotate_grids()
   sz_dens.save('slice/Data_%06d_%d_Slice_z_density_halo.png'%(num, halo), mpl_kwargs={"dpi":dpi} )

   sz_dens.clear_annotations()
   sz_dens.annotate_grids()
   sz_dens.annotate_timestamp( time_unit='Gyr', corner='upper_right' )
   sz_dens.zoom(10)
   sz_dens.save('slice/Data_%06d_%d_Slice_z_density_soliton.png'%(num, halo), mpl_kwargs={"dpi":dpi} )
