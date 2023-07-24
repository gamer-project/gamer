import argparse
import sys
import yt

# load the command-line parameters
parser = argparse.ArgumentParser( description='Plot the gas slices and projections' )

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
prefix      = '../'

colormap    = 'algae'
width_kpc   = 30
center_mode = 'c'
dpi         = 150


yt.enable_parallelism()

ts = yt.DatasetSeries( [ prefix+'/Data_%06d'%idx for idx in range(idx_start, idx_end+1, didx) ] )
#ts = yt.load( 'Data_??????' )

for ds in ts.piter():

#  define density^2 for calculating the weighted temperature
   def _density_square( field, data ):
      return data["density"]**2

   ds.add_field( ("gas", "density_square"), function=_density_square, sampling_type="cell", units="g**2/cm**6" )


#  density slice -- face-on
   sz_dens = yt.SlicePlot( ds, 'z', 'density', center=center_mode, width=(width_kpc,'kpc') )
   sz_dens.set_zlim( 'density', 1.0e-26, 3.0e-23 )
   sz_dens.set_cmap( 'density', colormap )
   sz_dens.annotate_timestamp( time_unit='Myr', corner='upper_right' )
   sz_dens.annotate_grids( periodic=False )
   sz_dens.save( mpl_kwargs={"dpi":dpi} )


#  density slice -- edge-on
   sx_dens = yt.SlicePlot( ds, 'x', 'density', center=center_mode, width=(width_kpc,'kpc') )
   sx_dens.set_zlim( 'density', 1.0e-26, 3.0e-23 )
   sx_dens.set_cmap( 'density', colormap )
   sx_dens.annotate_timestamp( time_unit='Myr', corner='upper_right' )
   sx_dens.annotate_grids( periodic=False )
   sx_dens.save( mpl_kwargs={"dpi":dpi} )


#  density projection -- face-on
   pz_dens = yt.ProjectionPlot( ds, 'z', 'density', center=center_mode, width=(width_kpc,'kpc'), method='integrate', weight_field=None )
   pz_dens.set_zlim( 'density', 1.0e-5, 1.0e-1 )
   pz_dens.set_cmap( 'density', colormap )
   pz_dens.annotate_timestamp( time_unit='Myr', corner='upper_right' )
   pz_dens.save( mpl_kwargs={"dpi":dpi} )


#  density projection -- edge-on
   pz_dens = yt.ProjectionPlot( ds, 'x', 'density', center=center_mode, width=(width_kpc,'kpc'), method='integrate', weight_field=None )
   pz_dens.set_zlim( 'density', 1.0e-5, 1.0e-1 )
   pz_dens.set_cmap( 'density', colormap )
   pz_dens.annotate_timestamp( time_unit='Myr', corner='upper_right' )
   pz_dens.save( mpl_kwargs={"dpi":dpi} )


#  temperature slice -- face-on
   sz_temp = yt.SlicePlot( ds, 'z', 'temperature', center=center_mode, width=(width_kpc,'kpc') )
   sz_temp.set_zlim( 'temperature', 1.0e2, 1.0e4 )
   sz_temp.set_cmap( 'temperature', colormap )
   sz_temp.annotate_timestamp( time_unit='Myr', corner='upper_right' )
   sz_temp.annotate_grids( periodic=False )
   sz_temp.save( mpl_kwargs={"dpi":dpi} )


#  temperature slice -- edge-on
   sx_temp = yt.SlicePlot( ds, 'x', 'temperature', center=center_mode, width=(width_kpc,'kpc') )
   sx_temp.set_zlim( 'temperature', 1.0e2, 1.0e4 )
   sx_temp.set_cmap( 'temperature', colormap )
   sx_temp.annotate_timestamp( time_unit='Myr', corner='upper_right' )
   sx_temp.annotate_grids( periodic=False )
   sx_temp.save( mpl_kwargs={"dpi":dpi} )


#  temperature projection -- face-on
   pz_temp = yt.ProjectionPlot( ds, 'z', 'temperature', center=center_mode, width=(width_kpc,'kpc'),
                                method='integrate', weight_field='density_square' )
   pz_temp.set_zlim( 'temperature', 1.0e2, 1.0e4 )
   pz_temp.set_cmap( 'temperature', colormap )
   pz_temp.annotate_timestamp( time_unit='Myr', corner='upper_right' )
   pz_temp.save( mpl_kwargs={"dpi":dpi} )


#  temperature projection -- edge-on
   px_temp = yt.ProjectionPlot( ds, 'x', 'temperature', center=center_mode, width=(width_kpc,'kpc'),
                                method='integrate', weight_field='density_square' )
   px_temp.set_zlim( 'temperature', 1.0e2, 1.0e4 )
   px_temp.set_cmap( 'temperature', colormap )
   px_temp.annotate_timestamp( time_unit='Myr', corner='upper_right' )
   px_temp.save( mpl_kwargs={"dpi":dpi} )


#  velocity magnitude slice -- face-on
   sz_vabs = yt.SlicePlot( ds, 'z', 'velocity_magnitude', center=center_mode, width=(width_kpc,'kpc') )
   sz_vabs.set_unit( 'velocity_magnitude', 'km/s' )
   sz_vabs.set_log ( 'velocity_magnitude', False )
   sz_vabs.set_zlim( 'velocity_magnitude', 5.0e1, 2.3e2 )
   sz_vabs.set_cmap( 'velocity_magnitude', colormap )
   sz_vabs.annotate_quiver('velocity_x', 'velocity_y', 16)
   sz_vabs.annotate_timestamp( time_unit='Myr', corner='upper_right' )
   sz_vabs.save( mpl_kwargs={"dpi":dpi} )


#  velocity magnitude slice -- edge-on
   sz_vabs = yt.SlicePlot( ds, 'x', 'velocity_magnitude', center=center_mode, width=(width_kpc,'kpc') )
   sz_vabs.set_unit( 'velocity_magnitude', 'km/s' )
   sz_vabs.set_log ( 'velocity_magnitude', False )
   sz_vabs.set_zlim( 'velocity_magnitude', 5.0e1, 2.3e2 )
   sz_vabs.set_cmap( 'velocity_magnitude', colormap )
   sz_vabs.annotate_timestamp( time_unit='Myr', corner='upper_right' )
   sz_vabs.save( mpl_kwargs={"dpi":dpi} )

