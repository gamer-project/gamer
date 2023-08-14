import argparse
import sys
import yt
import numpy as np


# -------------------------------------------------------------------------------------------------------------------------
# user-specified parameters
dpi         = 150
nbin        = 256*2


#-------------------------------------------------------------------------------------------------------------------------
# load the command-line parameters
parser = argparse.ArgumentParser( description='Get the halo density profile' )

parser.add_argument( '-s', action='store', required=True,  type=int, dest='idx_start',
                     help='first data index' )
parser.add_argument( '-e', action='store', required=True,  type=int, dest='idx_end',
                     help='last data index' )
parser.add_argument( '-d', action='store', required=False, type=int, dest='didx',
                     help='delta data index [%(default)d]', default=1 )

args=parser.parse_args()

idx_start = args.idx_start
idx_end   = args.idx_end
didx      = args.didx

# print command-line parameters
print( '\nCommand-line arguments:' )
print( '-------------------------------------------------------------------' )
for t in range( len(sys.argv) ):
   print( str(sys.argv[t])),
print( '' )
print( '-------------------------------------------------------------------\n' )

# -------------------------------------------------------------------------------------------------------------------------
# output figures
field       = ('gamer', 'Dens')
center_mode = 'max'
Center      = np.loadtxt('../../Record__Center', skiprows=1, dtype=float)

for idx in range(idx_start, idx_end+1, didx):
   ds             = yt.load( '../../Data_%06d'%idx )
   current_step   = ds.parameters["Step"]
   print("Current Simulation Time = %.5e [code units]"%ds.parameters["Time"][0])
   print("Current Simulation Step = %i "%current_step)

   sp = ds.sphere( Center[current_step,3:6], 0.5*ds.domain_width.to_value().max() )
   prof = yt.ProfilePlot( sp, 'radius', field, weight_field='cell_volume', n_bins=nbin, x_log=True, y_log={field:True} )
   prof.set_unit( 'radius', 'kpc' )
   prof.set_unit( field, 'g/cm**3' )

   A = prof.profiles[0].x.d
   B = prof.profiles[0][field].d
   C = [A,B]
   Data = np.asarray(C)
   Data = np.delete(Data,  np.nonzero(Data[1]==0)[0] , axis = 1)
   np.save("Halo_Dens_"+str(ds),Data)

   prof.set_xlim( 5.0e-2, 1.0e+2 )
   prof.save( mpl_kwargs={"dpi":dpi} )
