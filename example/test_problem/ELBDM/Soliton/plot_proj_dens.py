import argparse
import sys
import yt


# load the command-line parameters
parser = argparse.ArgumentParser( description='Plot the gas slices and projections' )

parser.add_argument( '-p', action='store', required=False, type=str, dest='prefix',
                     help='prefix [%(default)s]', default='./' )
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


idx_start   = args.idx_start
idx_end     = args.idx_end
didx        = args.didx
prefix      = args.prefix

field       = 'density'
colormap    = 'algae'
center_mode = 'c'
dpi         = 150


yt.enable_parallelism()

ts = yt.DatasetSeries( [ prefix+'/Data_%06d'%idx for idx in range(idx_start, idx_end+1, didx) ] )

for ds in ts.piter():

   pz = yt.ProjectionPlot( ds, 'z', field, center=center_mode, weight_field='density' )
   pz.set_zlim( field, 1.0e-12, 1.0e-4 )
   pz.set_cmap( field, colormap )
   pz.set_unit( field, 'code_mass/code_length**3' )
   pz.set_axes_unit( 'code_length' )
   pz.annotate_timestamp( time_unit='code_time', corner='upper_right', time_format='t = {time:5.1f} {units}' )
   pz.save( mpl_kwargs={"dpi":dpi} )

