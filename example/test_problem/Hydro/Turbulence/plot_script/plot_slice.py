import argparse
import sys
import yt
import numpy as np


# -------------------------------------------------------------------------------------------------------------------------
# user-specified parameters
colormap    = 'algae'
dpi         = 150


#-------------------------------------------------------------------------------------------------------------------------
# load the command-line parameters
parser = argparse.ArgumentParser( description='Plot slices' )

parser.add_argument( '-s', action='store', required=True,  type=int, dest='idx_start',
                     help='first data index' )
parser.add_argument( '-e', action='store', required=True,  type=int, dest='idx_end',
                     help='last data index' )
parser.add_argument( '-d', action='store', required=False, type=int, dest='didx',
                     help='delta data index [%(default)d]', default=1 )

args=parser.parse_args()

idx_start   = args.idx_start
idx_end     = args.idx_end
didx        = args.didx

# print command-line parameters
print( '\nCommand-line arguments:' )
print( '-------------------------------------------------------------------' )
for t in range( len(sys.argv) ):
   print( str(sys.argv[t]))
print( '' )
print( '-------------------------------------------------------------------\n' )


# -------------------------------------------------------------------------------------------------------------------------
# output figures
yt.enable_parallelism()
ts = yt.DatasetSeries( [ '../Data_%06d'%idx for idx in range(idx_start, idx_end+1, didx) ] )

field = ('gas', 'vorticity_magnitude')

for ds in ts.piter():
   slc = yt.SlicePlot( ds, 0, fields = field, center = 'c' )
   slc.set_background_color( field )
   slc.set_zlim( field, 1.0e-1, 1.0e+2, dynamic_range=None)
   slc.set_cmap( field, colormap )
   slc.set_font( {'size':16} )
   slc.set_axes_unit( 'code_length' )
   slc.annotate_grids( periodic=False )
   slc.annotate_timestamp( time_unit='code_time', corner='upper_right', text_args={'color':'k'} )
   slc.save( mpl_kwargs={"dpi":dpi} )

