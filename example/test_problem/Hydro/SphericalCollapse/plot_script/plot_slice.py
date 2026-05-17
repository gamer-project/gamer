import argparse
import sys
import yt

# load the command-line parameters
parser = argparse.ArgumentParser( description='Plot gas slices' )

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

field       = ['temperature', 'DualStatus']
center_mode = 'c'
dpi         = 150


yt.enable_parallelism()

ts = yt.DatasetSeries( [ prefix+'/Data_%06d'%idx for idx in range(idx_start, idx_end+1, didx) ] )

for ds in ts.piter():

   sz = yt.SlicePlot( ds, 'z', field, center=center_mode  )

   sz.set_log(  field[0], True )
   sz.set_cmap( field[0], 'arbre' )
   sz.set_zlim( field[0], 1.0e-3, 1.0e4 )
#  sz.set_unit( field[0], 'code_mass/code_length**3' )

   sz.set_log(  field[1], False )
   sz.set_cmap( field[1], 'algae' )
   sz.set_zlim( field[1], 0.0, 5.0 )
#  sz.set_unit( field[1], 'code_mass/code_length**3' )

#  sz.set_axes_unit( 'code_length' )
   sz.annotate_timestamp( corner='upper_right', redshift=True, time=False, text_args={'color':'k'})
   sz.annotate_grids( periodic=False )
   sz.save( mpl_kwargs={"dpi":dpi} )
