import argparse
import sys
import yt

# load the command-line parameters
parser = argparse.ArgumentParser( description='Plot gas density slices for the blast wave test' )

parser.add_argument( '-s', action='store', required=True,  type=int, dest='idx_start',
                     help='first data index' )
parser.add_argument( '-e', action='store', required=True,  type=int, dest='idx_end',
                     help='last data index' )
parser.add_argument( '-d', action='store', required=False, type=int, dest='didx',
                     help='delta data index [%(default)d]', default=1 )
parser.add_argument( '-i', action='store', required=False,  type=str, dest='prefix',
                     help='data path prefix [%(default)s]', default='./' )

args=parser.parse_args()

# take note
print( '\nCommand-line arguments:' )
print( '-------------------------------------------------------------------' )
for t in range( len(sys.argv) ):
   print( str(sys.argv[t]) ),
print( '' )
print( '-------------------------------------------------------------------\n' )


idx_start   = args.idx_start
idx_end     = args.idx_end
didx        = args.didx
prefix      = args.prefix

colormap    = 'arbre'
field       = 'density'    # to change the target field, one must modify set_unit() accordingly
center_mode = 'c'
dpi         = 150


yt.enable_parallelism()

ts = yt.load( [ prefix+'Data_%06d'%idx for idx in range(idx_start, idx_end+1, didx) ] )

for ds in ts.piter():

   sz = yt.SlicePlot( ds, 'z', field, center=center_mode  )
   sz.set_zlim( field, 0.0, 3.5 )
   sz.set_log( field, False )
   sz.set_cmap( field, colormap )
   sz.set_unit( field, 'code_mass/code_length**3' )
   sz.set_axes_unit( 'code_length' )
   sz.annotate_timestamp( time_unit='code_time', corner='upper_right', time_format='t = {time:.4f} {units}' )
   sz.annotate_grids( periodic=False )
   sz.save( mpl_kwargs={"dpi":dpi} )
