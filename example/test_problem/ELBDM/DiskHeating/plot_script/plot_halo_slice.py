import argparse
import sys
import yt
import numpy as np

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
   print( str(sys.argv[t]))
print( '' )
print( '-------------------------------------------------------------------\n' )


idx_start   = args.idx_start
idx_end     = args.idx_end
didx        = args.didx
prefix      = args.prefix

colormap    = 'algae'
dpi         = 150
center_mode   = 'max'

yt.enable_parallelism()

field = ('Dens')

for idx in range(idx_start, idx_end+1, didx):
   ds = yt.load( prefix+'/Data_%06d'%idx )
   ds.periodicity = (True, True, True)
   dens = yt.SlicePlot( ds, 0, fields = field, center = center_mode, width = (60 ,'kpc'))
   dens.set_background_color( field )
   dens.set_zlim( field, 1.0e+1, 1.0e+9, dynamic_range=None)
   dens.set_cmap( field, colormap )
   dens.set_font( {'size':16} )
   dens.set_axes_unit( 'kpc' )
   dens.annotate_grids()
   dens.annotate_sphere( center = Center[Id[idx],3:6], radius=(0.05, 'kpc'),circle_args={'color':'red'})
   dens.annotate_timestamp( time_unit='Myr', corner='upper_right', text_args={'color':'k'} )
   dens.save( mpl_kwargs={"dpi":dpi} )

