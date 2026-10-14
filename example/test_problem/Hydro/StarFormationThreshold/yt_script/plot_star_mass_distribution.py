import argparse
import sys
import yt
import numpy as np

# load the command-line parameters
parser = argparse.ArgumentParser( description='Plot the star mass distribution' )

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
print( ' '.join(map(str, sys.argv)) )
print( '-------------------------------------------------------------------\n' )


idx_start   = args.idx_start
idx_end     = args.idx_end
didx        = args.didx
prefix      = '../'

colormap    = 'algae'
dpi         = 150


yt.enable_parallelism()

ts = yt.DatasetSeries( [ prefix+'/Data_%06d'%idx for idx in range(idx_start, idx_end+1, didx) ] )

x_lim_min = 0.0
x_lim_max = 2.5e2
y_lim_min = 0.0
y_lim_max = 4.0e3

for ds in ts.piter():

#  check particle
   if len( ds.all_data()['all', 'particle_ones'] ) == 0:
      print( 'WARNING: There is no particle in %s !!'%ds )
      continue

#  plot
   p = yt.PhasePlot( ds, ('all', 'ParCreTime'), ('all', 'particle_mass'), ('all', 'particle_ones'),
                     weight_field=None, x_bins=256, y_bins=256 )

   p.set_unit( 'ParCreTime',    'Myr'  )
   p.set_unit( 'particle_mass', 'Msun' )
   p.set_log(  'ParCreTime',     False )
   p.set_log(  'particle_mass',  False )
   p.set_xlim( x_lim_min, x_lim_max )
   p.set_ylim( y_lim_min, y_lim_max )
   p.set_cmap( ('all', 'particle_ones'), colormap )
   p.save( mpl_kwargs={'dpi':dpi} )
