import argparse
import sys
import yt
import numpy as np

# load the command-line parameters
parser = argparse.ArgumentParser( description='Density profile' )

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
   print( str(sys.argv[t])),
print( '' )
print( '-------------------------------------------------------------------\n' )

idx_start = args.idx_start
idx_end   = args.idx_end
didx      = args.didx
prefix    = args.prefix

field = ('gamer', 'Pote')

center_mode = 'max'
dpi         = 150
nbin        = 256*2

Center = np.loadtxt('../Record__Center', skiprows=1, dtype=float)

#==========================================================================
# input the steps of the output snapshots in Record__Center
step = [0, 117, 234, 351, 468, 585, 702, 819, 936, 1053, 1170, 
      1287, 1404, 1521, 1638, 1755, 1872, 1989, 2106, 2223, 2340,
      2457, 2574, 2691, 2808, 2925]
#==========================================================================
for idx in range(idx_start, idx_end+1, didx):
   ds = yt.load( prefix+'/Data_%06d'%idx )
   sp = ds.sphere( Center[step[idx],3:6], 0.5*ds.domain_width.to_value().max() )
   prof = yt.ProfilePlot( sp, 'radius', field, weight_field='cell_volume', n_bins=nbin, x_log=True, y_log={field:True} )
   prof.set_unit( 'radius', 'kpc' )
#   prof.set_xlim( 6.0e-2, 1.0e+2 )
#   prof.set_ylim( field, 1.0e-27, 1.0e-21 )
   prof.set_unit( field, 'cm**2/s**2' )
   prof.save( mpl_kwargs={"dpi":dpi} )
   A = prof.profiles[0].x.d
   B = prof.profiles[0][field].d
   C = [A,B]
   Data = np.asarray(C)
   np.save("Halo_Pote_"+str(ds),Data)

