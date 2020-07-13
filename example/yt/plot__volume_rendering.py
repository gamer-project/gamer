
# ref: https://yt-project.org/doc/visualizing/volume_rendering.html#

import argparse
import sys
import yt


# user-specified parameters
field      = 'Dens'           # target field
zoom       = 2.0              # zoom-in factor
resolution = 512              # image resolution
clip       = 8.0              # floor bright pixels for better contrast
bounds     = (1.0e0, 1.0e6)   # bounds of the transfer function


# load the command-line parameters
parser = argparse.ArgumentParser( description='Volume rendering' )

parser.add_argument( '-i', action='store', required=False, type=str, dest='prefix_in',
                     help='input path prefix [%(default)s]', default='./' )
parser.add_argument( '-o', action='store', required=False, type=str, dest='prefix_out',
                     help='output filename prefix [%(default)s]', default='fig__volume-rendering' )
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

#  create a scene
   sc = yt.create_scene( ds, field=field, lens_type='perspective' )

#  get a reference to the VolumeSource associated with this scene
   src = sc[0]

#  set the target field
   src.set_field( field )

#  set that the transfer function should be evaluated in log space
   src.tfh.set_log( True )

#  make underdense regions appear opaque
#  src.tfh.grey_opacity = True

#  zoom-in
   sc.camera.zoom( zoom )

#  set the image resolution
   sc.camera.resolution = resolution

#  annotate grids
#  sc.annotate_grids( ds, alpha=0.01 )

#  set the bounds of the transfer function
   src.tfh.set_bounds( bounds )

#  plot the transfer function, along with the CDF of the target field to
#  see how the transfer function corresponds to structure in the CDF
   src.tfh.plot( 'fig__transfer-function.png', profile_field=field )

#  begin rendering
   sc.render()

#  save the image, flooring especially bright pixels for better contrast
   sc.save( prefix_out+'_'+ds.basename+'.png', sigma_clip=clip )

