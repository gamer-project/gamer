
#ref: https://yt-project.org/doc/visualizing/plts.html#off-axis-projection-plts

import argparse
import sys
import yt


# user-specified parameters
field    = ( 'gas', 'density' )  # target field(s)
normal   = [+1,+1,+1]            # projection direction
north    = [-1,-1,+2]            # 'up' direction in the image
colormap = 'arbre'
dpi      = 150


# load the command-line parameters
parser = argparse.ArgumentParser( description='Plot the gas slices and projections' )

parser.add_argument( '-i', action='store', required=False, type=str, dest='prefix_in',
                     help='input path prefix [%(default)s]', default='./' )
parser.add_argument( '-o', action='store', required=False, type=str, dest='prefix_out',
                     help='output filename prefix [%(default)s]', default='fig__off-axis-projection' )
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

idx_start  = args.idx_start
idx_end    = args.idx_end
didx       = args.didx
prefix_in  = args.prefix_in
prefix_out = args.prefix_out


yt.enable_parallelism()
ts = yt.load( [ prefix_in+'/Data_%06d'%idx for idx in range(idx_start, idx_end+1, didx) ] )

for ds in ts.piter():

#  project data
   plt = yt.OffAxisProjectionPlot( ds, normal=normal, fields=field, north_vector=north )

#  save image
   plt.set_cmap( field, colormap )
   plt.annotate_timestamp( time_unit='code_time', corner='upper_right' )
   plt.save( prefix_out+'_'+ds.basename+'.png', mpl_kwargs={'dpi':dpi} )

#  access the image buffer
   data = plt._frb.data[field]
   print( "units = %s" % data.units )
   print( "max   = %s" % data.d.max() )
