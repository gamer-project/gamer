import argparse
import sys
import yt


# load the command-line parameters
parser = argparse.ArgumentParser( description='Plot the gas slices and projections' )

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
for t in range( len(sys.argv) ):
   print str(sys.argv[t]),
print( '' )
print( '-------------------------------------------------------------------\n' )


idx_start   = args.idx_start
idx_end     = args.idx_end
didx        = args.didx
prefix      = args.prefix

colormap    = 'afmhot'
center_mode = 'c'
dpi         = 150


yt.enable_parallelism()

ts = yt.load( [ prefix+'/Data_%06d'%idx for idx in range(idx_start, idx_end+1, didx) ] )

for ds in ts.piter():

   sz = yt.SlicePlot( ds, 'z', 'temperature', center=center_mode, width=[80,80] )
   sz.set_unit( 'temperature', 'keV', equivalency='thermal' )
   sz.set_zlim( 'temperature', 1.0e0, 2.0e2 )
   sz.set_cmap( 'temperature', colormap )
   sz.annotate_timestamp( time_unit='Myr', corner='upper_right' )
   sz.save( mpl_kwargs={"dpi":dpi} )

