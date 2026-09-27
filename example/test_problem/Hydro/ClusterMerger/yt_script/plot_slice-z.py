import argparse
import sys
import yt
from yt.units import Mpc

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
print( ' '.join(map(str, sys.argv)) )
print( '-------------------------------------------------------------------\n' )


idx_start = args.idx_start
idx_end   = args.idx_end
didx      = args.didx
prefix    = args.prefix

colormap_dens = 'algae'
colormap_temp = 'afmhot'
colormap_cdm  = 'arbre'
center_mode   = 'c'
dpi           = 150


# define the particle filter for SMBH particles
def smbh( pfilter, data ):
   filter = data[ "all", "ParType" ] == 4
   return filter

yt.add_particle_filter( "smbh", function=smbh, filtered_type="all", requires=["ParType"] )


yt.enable_parallelism()
ts = yt.DatasetSeries( [ prefix+'/Data_%06d'%idx for idx in range(idx_start, idx_end+1, didx) ] )

for ds in ts.piter():

#  apply the particle filters to the dataset
   ds.add_particle_filter( 'smbh' )


#  density
   field   = 'density'
   sz_dens = yt.SlicePlot( ds, 'z', field, center=center_mode )

   sz_dens.set_zlim( field, 1.0e-30, 1.0e-25 )
   sz_dens.set_cmap( field, colormap_dens )
   sz_dens.annotate_timestamp( time_unit='Myr', corner='upper_right' )
   sz_dens.annotate_grids( periodic=False )

   sz_dens.save( mpl_kwargs={"dpi":dpi} )


#  temperature
   field   = 'temperature'
   sz_temp = yt.SlicePlot( ds, 'z', field, center=center_mode )

   sz_temp.set_zlim( field, 2.0e7, 1.5e8 )
   sz_temp.set_cmap( field, colormap_temp )
   sz_temp.annotate_timestamp( time_unit='Myr', corner='upper_right' )

   sz_temp.save( mpl_kwargs={"dpi":dpi} )


#  CDM particles
   field    = 'particle_mass'
   cdm_mass = yt.ParticlePlot( ds, 'particle_position_x', 'particle_position_y', field,
                               center=center_mode, depth=(8.0*Mpc) )

   cdm_mass.set_background_color( field )
   cdm_mass.set_unit( field, 'Msun' )
   cdm_mass.set_zlim( field, 1.0e9, 1.0e12 )
   cdm_mass.set_cmap( field, colormap_cdm )
   cdm_mass.set_colorbar_label( field, r'Dark matter mass [$M_{\odot}$]' )
   cdm_mass.annotate_timestamp( time_unit='Myr', corner='upper_right' )
   cdm_mass.annotate_particles( width=(10.0*Mpc), p_size=10.0, col='b', ptype='smbh' )

   cdm_mass.save( mpl_kwargs={"dpi":dpi} )

