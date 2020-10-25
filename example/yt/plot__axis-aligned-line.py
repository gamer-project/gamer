
#ref: https://yt-project.org/docs/dev/visualizing/manual_plotting.html#line-plots

import argparse
import sys
import yt
import numpy as np
from matplotlib import pyplot as plt


# hard-coded parameters (in code units)
fields = [ 'Dens', 'MomX' ]
axis   = 'x'
center = [ 0.0, 0.0 ]
dpi    = 150


# load the command-line parameters
parser = argparse.ArgumentParser( description='Print and plot the data along an axis-aligned line' )

parser.add_argument( '-i', action='store', required=False, type=str, dest='prefix_in',
                     help='input path prefix [%(default)s]', default='./' )
parser.add_argument( '-o', action='store', required=False, type=str, dest='prefix_out_fig',
                     help='filename prefix of output figures [%(default)s]', default='fig__line' )
parser.add_argument( '-t', action='store', required=False, type=str, dest='prefix_out_tab',
                     help='filename prefix of output tables [%(default)s]', default='tab__line' )
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

idx_start      = args.idx_start
idx_end        = args.idx_end
didx           = args.didx
prefix_in      = args.prefix_in
prefix_out_fig = args.prefix_out_fig
prefix_out_tab = args.prefix_out_tab


yt.enable_parallelism()
ts = yt.load( [ prefix_in+'/Data_%06d'%idx for idx in range(idx_start, idx_end+1, didx) ] )

for ds in ts.piter():

#  1. get the data along the target line
   ray    = ds.ortho_ray( axis, center )  # unsorted data
   r_sort = np.argsort( ray[axis] )       # sorted coordinates


#  2. plot the first field as an example
   plt.semilogy(  np.array( ray[axis][r_sort] ), np.array( ray[ fields[0] ][r_sort] )  )
   plt.xlabel( axis )
   plt.ylabel( fields[0] )

   plt.savefig( prefix_out_fig+'_'+ds.basename+'.png', bbox_inches='tight', pad_inches=0.20, dpi=dpi )
   plt.clf()


#  3. save the table
   tab = open( prefix_out_tab+'_'+ds.basename, 'w' )

#  header
   tab.write( '#%19s' % axis )
   for v in range( len(fields) ):  tab.write( ' %13s' % fields[v] )
   tab.write( '\n' )

#  data
   for t in range( r_sort.size ):
      tab.write( '%20.14e' % ray[axis][r_sort][t] )

      for v in range( len(fields) ):
         tab.write( ' %13.6e' % ray[ fields[v] ][r_sort][t] )

      tab.write( '\n' )

   tab.close()
