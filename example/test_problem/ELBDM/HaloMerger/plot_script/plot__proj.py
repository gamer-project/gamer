import argparse
import sys
import yt

# load the command-line parameters
parser = argparse.ArgumentParser( description='Projection of mass density' )

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
print( ' '.join(map(str, sys.argv)) )
print( '-------------------------------------------------------------------\n' )


idx_start = args.idx_start
idx_end   = args.idx_end
didx      = args.didx
prefix    = args.prefix

colormap_dens = 'algae'
dpi           = 150

yt.enable_parallelism()

ts = yt.DatasetSeries( [ prefix+'/Data_%06d'%idx for idx in range(idx_start, idx_end+1, didx) ] )
for ds in ts.piter():
   for center_mode in ['c','m']:
      for direction in ['y', 'z']:
         for field in ['density']:
            pz_dens = yt.ProjectionPlot( ds, direction, field, center=center_mode )

            pz_dens.set_axes_unit( 'kpc' )
            pz_dens.set_unit( field, "Msun/kpc**2" )
            pz_dens.set_zlim( field, 1.0e+3, 1.0e8 )
            pz_dens.set_cmap( field, colormap_dens )
            pz_dens.set_background_color( field )
            pz_dens.annotate_timestamp( time_unit='Gyr', corner='upper_right' )
            if center_mode == 'm':
               pz_dens.zoom(4)

            pz_dens.save("%s_%s"%(ds, center_mode), mpl_kwargs={"dpi":dpi} )
            pz_dens.annotate_grids()
            pz_dens.save( "%s_%s_grids"%(ds, center_mode), mpl_kwargs={"dpi":dpi} )
