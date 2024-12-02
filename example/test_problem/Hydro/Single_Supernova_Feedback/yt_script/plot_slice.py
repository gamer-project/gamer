import argparse
import sys
import yt
import matplotlib.pyplot as plt

# load the command-line parameters
# parser = argparse.ArgumentParser( description='Plot gas density slices for the blast wave test' )
parser = argparse.ArgumentParser( description='Plot gas density projection for the feedback test' )

parser.add_argument( '-s', action='store', required=True,  type=int, dest='idx_start',
                     help='first data index' )
parser.add_argument( '-e', action='store', required=True,  type=int, dest='idx_end',
                     help='last data index' )
parser.add_argument( '-d', action='store', required=False, type=int, dest='didx',
                     help='delta data index [%(default)d]', default=1 )
parser.add_argument( '-i', action='store', required=False,  type=str, dest='prefix',
                     help='data path prefix [%(default)s]', default='./' )

args=parser.parse_args()

# take note
print( '\nCommand-line arguments:' )
print( '-------------------------------------------------------------------' )
for t in range( len(sys.argv) ):
   print( str(sys.argv[t]) ),
print( '' )
print( '-------------------------------------------------------------------\n' )


idx_start   = args.idx_start
idx_end     = args.idx_end
didx        = args.didx
prefix      = args.prefix

colormap    = 'viridis' #density
#colormap    = 'magma'  #temperature
field       = 'density'    # to change the target field, one must modify set_unit() accordingly
#field       = 'temperature'    # to change the target field, one must modify set_unit() accordingly
#field       = 'pressure'    # to change the target field, one must modify set_unit() accordingly
center_mode = 'c'
dpi         = 300


yt.enable_parallelism()

ts = yt.load( [ prefix+'/Data_%06d'%idx for idx in range(idx_start, idx_end+1, didx) ] )

for ds in ts.piter():

   fig = plt.figure()
   sz = yt.ProjectionPlot( ds, 'z', field, center=center_mode  )
#   sz = yt.ProjectionPlot( ds, 'z', field, weight_field='density', center=center_mode  )
#   sz = yt.SlicePlot( ds, 'z', field, center=center_mode  )
   sz.set_zlim( field, 3.0e0, 8.0e0 ) #density
#   sz.set_zlim( field, 1.0e2, 1.0e5 ) #temperature
   sz.set_log( field, True )
   sz.set_cmap( field, colormap )
#   sz.set_unit( field, 'code_mass/code_length**2' )
   sz.set_unit( field, 'Msun/pc**2' )
#   sz.set_unit( field, 'dyn/cm**2' )
#   sz.set_unit( field, 'K' )
   sz.set_axes_unit( 'pc' )
   sz.annotate_timestamp( time_unit='Myr', corner='upper_right', time_format='t = {time:.2f} {units}' )
#   sz.annotate_grids( periodic=False )
   plt.tight_layout()
   sz.save( mpl_kwargs={"dpi":dpi} )
