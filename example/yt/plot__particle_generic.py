
# ref: https://yt-project.org/docs/dev/examining/generic_particle_data.html

import argparse
import sys
import numpy as np
import yt
from yt.units import kpc, Msun, Gyr


# hard-coded parameters (in code units)
field       = 'particle_mass'
time        = 1.23
proj_axis   = 'z'
plot_center = [ 0.0, 0.0, 0.0 ]
plot_width  = [ 0.4, 0.4, 0.4 ]
unit_l      = kpc
unit_m      = Msun
unit_t      = Gyr
nref        = 256
bext        = 1.0
colormap    = 'arbre'
dpi         = 150


# load the command-line parameters
parser = argparse.ArgumentParser( description='Plot the projected particle mass' )

parser.add_argument( '-i', action='store', required=False, type=str, dest='prefix_in',
                     help='input path prefix [%(default)s]', default='./' )
parser.add_argument( '-o', action='store', required=False, type=str, dest='prefix_out',
                     help='output filename prefix [%(default)s]', default='fig__particle' )
parser.add_argument( '-s', action='store', required=True,  type=int, dest='idx_start',
                     help='first data index' )
parser.add_argument( '-e', action='store', required=True,  type=int, dest='idx_end',
                     help='last data index' )
parser.add_argument( '-d', action='store', required=False, type=int, dest='didx',
                     help='delta data index [%(default)d]', default=1 )
parser.add_argument( '-f', action='store', required=False, type=bool, dest='float',
                     help='single precision [%(default)r]', default=False )

args=parser.parse_args()

# take note
print( '\nCommand-line arguments:' )
print( '-------------------------------------------------------------------' )
for t in range( len(sys.argv) ):
   print str(sys.argv[t]),
print( '' )
print( '-------------------------------------------------------------------\n' )

prefix_in  = args.prefix_in
prefix_out = args.prefix_out
idx_start  = args.idx_start
idx_end    = args.idx_end
didx       = args.didx
ftype      = np.float if args.float else np.double


for idx in range(idx_start, idx_end+1, didx):

#  load data
   file_in = prefix_in + '/Data_%06d'%idx
   data    = np.fromfile( file=file_in, dtype=ftype )
   data    = data.reshape( 7, -1 )


#  send data to yt
   m = data[0]
   x = data[1]
   y = data[2]
   z = data[3]

   data_dic = { 'particle_mass'      : m,
                'particle_position_x': x,
                'particle_position_y': y,
                'particle_position_z': z }

   bbox_width = [ max(x)-min(x), max(y)-min(y), max(z)-min(z) ]
   bbox       = np.array(  [  [ min(x)-bext*bbox_width[0], max(x)+bext*bbox_width[0] ],
                              [ min(y)-bext*bbox_width[1], max(y)+bext*bbox_width[1] ],
                              [ min(z)-bext*bbox_width[2], max(z)+bext*bbox_width[2] ]  ]  )

   ds = yt.load_particles( data_dic, length_unit=unit_l, mass_unit=unit_m, time_unit=unit_t,
                           n_ref=nref, bbox=bbox, sim_time=time, periodicity=(False,False,False) )


#  plot
   p = yt.ParticleProjectionPlot( ds, proj_axis, fields=field, center=plot_center,
                                  width=(plot_width[0],plot_width[1]), depth=plot_width[2] )
#  p = yt.ParticlePlot( ds, 'particle_position_x', 'particle_position_y', 'particle_mass' )

   p.set_unit( field, 'Msun' )
   p.set_unit( 'particle_position_x', 'Mpc' )
   p.set_unit( 'particle_position_y', 'Mpc' )
   p.set_zlim( field, 1.0e0, 6.0e0 )
   p.set_cmap( field, colormap )
   p.annotate_timestamp( time_unit='Myr', corner='upper_right', text_args={'color':'k'} )

   p.save( prefix_out+'_%06d'%idx+'.png', mpl_kwargs={'dpi':dpi} )
