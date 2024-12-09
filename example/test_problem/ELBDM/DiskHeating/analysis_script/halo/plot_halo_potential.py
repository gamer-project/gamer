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
parser = argparse.ArgumentParser( description='Get the halo potential profile' )

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
field       = ('gamer', 'Pote')
center_mode = 'max'
Center      = np.loadtxt('../../Record__Center', skiprows=1, dtype=float)
if Center.ndim == 1:
   Center = Center.reshape(1,len(Center)) # reshape the array if there is only one row

for idx in range(idx_start, idx_end+1, didx):
   ds = yt.load( '../../Data_%06d'%idx )
   current_step   = ds.parameters["Step"]
   CellSize       = ds.parameters["CellSize"]
   MaxLevel       = ds.parameters["MaxLevel"]
   Resolution     = CellSize[MaxLevel]*ds.parameters["Unit_L"]/3.08568e+21 # spatial resolution in kpc
   print("Current Simulation Time = %.5e [code units]"%ds.parameters["Time"][0])
   print("Current Simulation Step = %i"%current_step)

   sp = ds.sphere( Center[current_step,3:6], 0.5*ds.domain_width.to_value().max() )
   prof = yt.create_profile( sp, 'radius', field,
                             weight_field='cell_volume',
                             n_bins=nbin,
                             units={"radius":"kpc", field:"cm**2/s**2"},
                             logs={"radius":True, field:False},
                             extrema={"radius":((Resolution/2., "kpc"), (0.5*ds.domain_width.to_value().max(), "code_length"))})
   plot = yt.ProfilePlot.from_profiles(prof)
   plot.set_log(field, False)
   plot.save( mpl_kwargs={"dpi":dpi} )
   A = prof.x
   B = prof[field]
   C = [A,B]
   Data = np.asarray(C)
   Data = np.delete(Data,  np.nonzero(Data[1]==0)[0] , axis = 1)
   np.save("Halo_Pote_"+str(ds),Data)
