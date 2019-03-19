import argparse
import sys
import yt


# hard-coded parameters (in code units)
field     = 'particle_mass'
proj_axis = 'z'
center    = [0.5, 0.5, 0.5]
width_x   = 1.0e-2
width_y   = 1.0e-2
width_z   = 1.0e-2
colormap  = 'arbre'
dpi       = 150


# load the command-line parameters
parser = argparse.ArgumentParser( description='Plot the projected particle mass' )

parser.add_argument( '-i', action='store', required=False, type=str, dest='prefix_in',
                     help='input path prefix [%(default)s]', default='./' )
parser.add_argument( '-o', action='store', required=False, type=str, dest='prefix_out',
                     help='output filename prefix [%(default)s]', default='fig__thin-slice-projection' )
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
#ts = yt.load( 'Data_??????' )


for ds in ts.piter():

#  plot
   p = yt.ParticleProjectionPlot( ds, proj_axis, fields=field, center=center, width=(width_x,width_y), depth=width_z )
#  p.set_unit( field, 'Msun' )
#  p.set_zlim( field, 1.0e6, 1.0e8 )
   p.set_cmap( field, colormap )
   p.annotate_timestamp( time_unit='code_time', corner='upper_right', text_args={'color':'k'} )

#  save the image
   p.save( prefix_out+'_'+ds.basename+'.png', mpl_kwargs={'dpi':dpi} )
