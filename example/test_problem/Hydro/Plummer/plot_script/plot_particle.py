import argparse
import sys
import yt

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
   print str(sys.argv[t]),
print( '' )
print( '-------------------------------------------------------------------\n' )


idx_start   = args.idx_start
idx_end     = args.idx_end
didx        = args.didx
prefix      = args.prefix

colormap    = 'arbre'
dpi         = 150
field       = 'particle_mass'


yt.enable_parallelism()

ts = yt.load( [ prefix+'/Data_%06d'%idx for idx in range(idx_start, idx_end+1, didx) ] )

for ds in ts.piter():

   par = yt.ParticlePlot( ds, 'particle_position_x', 'particle_position_y', field, center='c' )

   par.set_background_color( field )
#  par.set_zlim( field, 1.0e9, 1.0e12 )
   par.set_cmap( field, colormap )
   par.set_font( {'size':16} )
   par.set_axes_unit( 'code_length' )
   par.set_unit( field, 'code_mass' )
   par.annotate_timestamp( time_unit='code_time', corner='upper_right', text_args={'color':'k'} )
   par.save( mpl_kwargs={"dpi":dpi} )

