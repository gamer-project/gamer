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
parser = argparse.ArgumentParser( description='Plot the halo slices' )

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

field       = ('Dens')
center_mode = 'max'
Center      = np.loadtxt('../../Record__Center', skiprows=1, dtype=float)
if Center.ndim == 1:
   Center = Center.reshape(1,len(Center)) # reshape the array if there is only one row

for idx in range(idx_start, idx_end+1, didx):
   ds             = yt.load( '../../Data_%06d'%idx )
   if sys.version_info[0] == 2:
      ds.periodicity = (True, True, True)
   current_step   = ds.parameters["Step"]
   print("Current Simulation Time = %.5e [code units]"%ds.parameters["Time"][0])
   print("Current Simulation Step = %i"%current_step)

   dens = yt.SlicePlot( ds, 0, fields = field, center = Center[current_step,3:6], width = (60 ,'kpc'))
   dens.set_background_color( field )
   dens.set_zlim( field, 1.0e+1, 1.0e+9, dynamic_range=None)
   dens.set_cmap( field, colormap )
   dens.set_font( {'size':16} )
   dens.set_axes_unit( 'kpc' )
   dens.annotate_grids()
   dens.annotate_sphere( center = Center[current_step,3:6], radius=(0.05, 'kpc'),circle_args={'color':'red'})
   dens.annotate_timestamp( time_unit='Myr', corner='upper_right', text_args={'color':'k'} )
   dens.save( mpl_kwargs={"dpi":dpi} )
