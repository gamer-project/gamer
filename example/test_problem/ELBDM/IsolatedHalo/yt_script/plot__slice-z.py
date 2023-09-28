import argparse
import sys
import yt

# load the command-line parameters
parser = argparse.ArgumentParser( description='Slice of mass density' )

parser.add_argument( '-p', action='store', required=False, type=str, dest='prefix',
                     help='path prefix [%(default)s]', default='../' )
parser.add_argument( '-s', action='store', required=True,  type=int, dest='idx_start',
                     help='first data index' )
parser.add_argument( '-e', action='store', required=True,  type=int, dest='idx_end',
                     help='last data index' )
parser.add_argument( '-d', action='store', required=False, type=int, dest='didx',
                     help='delta data index [%(default)d]', default=1 )

args=parser.parse_args()

# take note
print( '\nCommand-line arguments:' )
print( '-------------------------------------------------------------------' )
print( ' '.join(map(str, sys.argv)) )
print( '-------------------------------------------------------------------\n' )


idx_start = args.idx_start
idx_end   = args.idx_end
didx      = args.didx
prefix    = args.prefix

field         = 'density'
colormap_dens = 'algae'
center_mode   = 'c'
dpi           = 150

yt.enable_parallelism()
ts = yt.DatasetSeries( [ prefix+'/Data_%06d'%idx for idx in range(idx_start, idx_end+1, didx) ] )

for ds in ts.piter():

   sz_dens = yt.SlicePlot( ds, 'z', field, center=center_mode )

   sz_dens.set_zlim( field, 1.0e-31, 1.0e-24 )
   sz_dens.set_cmap( field, colormap_dens )
   sz_dens.annotate_timestamp( time_unit='Gyr', corner='upper_right' )
   sz_dens.annotate_grids()

   sz_dens.save( mpl_kwargs={"dpi":dpi} )
