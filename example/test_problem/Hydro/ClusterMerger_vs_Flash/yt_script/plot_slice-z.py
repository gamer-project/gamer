import yt
from yt.units import Mpc

yt.enable_parallelism()

idx_start = 0
idx_end   = 100
didx      = 1
prefix    = '../Data_'

ts = yt.load( [ prefix+'%06d'%idx for idx in range(idx_start, idx_end+1, didx) ] )

for ds in ts.piter():

#  density
   sz_dens = yt.SlicePlot( ds, 'z', 'density', center='c' )

   sz_dens.set_zlim( 'density', 1.0e-30, 1.0e-25 )
   sz_dens.set_cmap('density', 'algae')
   sz_dens.annotate_timestamp( time_unit='Myr', corner='upper_right' )
   sz_dens.annotate_grids( periodic=False )

   sz_dens.save()


#  temperature
   sz_temp = yt.SlicePlot( ds, 'z', 'temperature', center='c' )

   sz_temp.set_zlim( 'temperature', 2.0e7, 1.5e8 )
   sz_temp.set_cmap('temperature', 'afmhot')
   sz_temp.annotate_timestamp( time_unit='Myr', corner='upper_right' )

   sz_temp.save()


   cdm_mass = yt.ParticlePlot( ds, 'particle_position_x', 'particle_position_y', 'particle_mass', center='c', depth=(8.0*Mpc) )

   cdm_mass.set_background_color( 'particle_mass' )
   cdm_mass.set_zlim( 'particle_mass', 1.0e9, 1.0e12 )
   cdm_mass.set_unit( 'particle_mass', 'Msun' )
   cdm_mass.set_cmap( 'particle_mass', 'arbre')
   cdm_mass.set_colorbar_label( 'particle_mass', 'Dark matter mass [$M_{\odot}$]' )
   cdm_mass.annotate_timestamp( time_unit='Myr', corner='upper_right' )

   cdm_mass.save()

