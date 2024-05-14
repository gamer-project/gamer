import argparse
import sys
import yt

# load the command-line parameters
parser = argparse.ArgumentParser( description='Plot magnetic field slices for the blast wave test' )

parser.add_argument( '-s', action='store', required=True,  type=int, dest='idx_start',
                     help='first data index' )
parser.add_argument( '-e', action='store', required=True,  type=int, dest='idx_end',
                     help='last data index' )
parser.add_argument( '-d', action='store', required=False, type=int, dest='didx',
                     help='delta data index [%(default)d]', default=1 )
parser.add_argument( '-i', action='store', required=False, type=str, dest='prefix',
                     help='data path prefix [%(default)s]', default='../' )

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

colormap    = 'arbre'
# to change the target field, one must modify set_unit() accordingly
field       = ('gas', 'magnetic_energy_density') if yt.__version__.split(".")[0] == "4" else ('gas', 'magnetic_energy')
center_mode = 'c'
dpi         = 150


yt.enable_parallelism()

ts = yt.DatasetSeries( [ prefix+'/Data_%06d'%idx for idx in range(idx_start, idx_end+1, didx) ] )

for ds in ts.piter():

   sz = yt.SlicePlot( ds, 'z', field, center=center_mode )
   sz.set_zlim( field, 1.0e-4, 1.0e0 )
#  sz.set_log( field, False )
   sz.set_cmap( field, colormap )
#  sz.set_unit( field, 'code_mass/code_length**3' )
   sz.set_axes_unit( 'code_length' )
   sz.annotate_timestamp( time_unit='code_time', corner='upper_right', time_format='t = {time:.4f} {units}' )
   sz.annotate_grids( periodic=False )
   sz.annotate_quiver( ("gas", "magnetic_field_x"), ("gas", "magnetic_field_y"), factor=16, normalize=True )
   sz.save( mpl_kwargs={"dpi":dpi} )
