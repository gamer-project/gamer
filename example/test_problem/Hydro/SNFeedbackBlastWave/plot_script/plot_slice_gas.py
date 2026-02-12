import argparse
import sys
import yt
import numpy as np

# load the command-line parameters
parser = argparse.ArgumentParser( description='Plot gas density slices for the SN blast wave test' )

parser.add_argument( '-s', action='store', required=True,  type=int, dest='idx_start',
                     help='first data index' )
parser.add_argument( '-e', action='store', required=True,  type=int, dest='idx_end',
                     help='last data index' )
parser.add_argument( '-d', action='store', required=False, type=int, dest='didx',
                     help='delta data index [%(default)d]', default=1 )
parser.add_argument( '-i', action='store', required=False, type=str, dest='prefix',
                     help='data path prefix [%(default)s]', default='../' )

args=parser.parse_args()

# take note
print( '\nCommand-line arguments:' )
print( '-------------------------------------------------------------------' )
print( ' '.join(map(str, sys.argv)) )
print( '-------------------------------------------------------------------\n' )


idx_start   = args.idx_start
idx_end     = args.idx_end
didx        = args.didx
prefix      = args.prefix

colormap    = 'arbre'
center_mode = 'c'
dpi         = 150


yt.enable_parallelism()

ts = yt.DatasetSeries( [ prefix+'/Data_%06d'%idx for idx in range(idx_start, idx_end+1, didx) ] )

for ds in ts.piter():

   ad            = ds.all_data()


   # define the derived field
   def _thermal_energy(field, data):
       return data['thermal_energy_density']*data['cell_volume']

   ds.add_field( ('gamer','thermal_energy'), function=_thermal_energy, units='code_mass*code_velocity**2', sampling_type='cell' )

   def _kinetic_energy(field, data):
       return data['kinetic_energy_density']*data['cell_volume']

   ds.add_field( ('gamer','kinetic_energy'), function=_kinetic_energy, units='code_mass*code_velocity**2', sampling_type='cell' )

   def _metal_density(field, data):
       return data['Metal']*data.ds.units.code_density

   ds.add_field( ('gamer','metal_density'), function=_metal_density, units='code_density', sampling_type='cell' )

   def _cell_metal_mass(field, data):
       return data['metal_density']*data['cell_volume']

   ds.add_field( ('gamer','cell_metal_mass'), function=_cell_metal_mass, units='code_mass', sampling_type='cell' )


   # Density
   sz = yt.SlicePlot( ds, 'z', 'density', center=center_mode )
   sz.set_cmap( 'density', colormap )
   sz.set_unit( 'density', 'g/cm**3' )
   sz.set_zlim( 'density', 0.0, 2.0e-22 )
   sz.set_log(  'density', False )
   sz.set_axes_unit( 'pc' )
   sz.annotate_timestamp( time_unit='Myr', corner='upper_right', time_format='t = {time:.1f} {units}' )
   sz.annotate_grids( periodic=False )
   sz.annotate_particles( width=ds.domain_width[2], p_size=80.0, col='y', marker='*' )
   sz.annotate_text( (-400, -370), 'Total internal energy = {: >10.2e}'.format( ad.quantities.total_quantity("thermal_energy").in_units("erg") ), coord_system='plot' )
   sz.annotate_text( (-400, -400), 'Total kinetic energy  = {: >10.2e}'.format( ad.quantities.total_quantity("kinetic_energy").in_units("erg") ), coord_system='plot' )
   sz.save( mpl_kwargs={"dpi":dpi} )
   sz.annotate_cell_edges()
   sz.zoom( 16 )
   sz.save( '%s_Slice_z_density_zoom-in.png'%ds, mpl_kwargs={"dpi":dpi} )


   # Temparature
   sz = yt.SlicePlot( ds, 'z', 'temperature', center=center_mode  )
   sz.set_cmap( 'temperature', colormap )
   sz.set_unit( 'temperature', 'K' )
   sz.set_zlim( 'temperature', 1.0e0, 1.0e6 )
   sz.set_log(  'temperature', True )
   sz.set_axes_unit( 'pc' )
   sz.annotate_timestamp( time_unit='Myr', corner='upper_right', time_format='t = {time:.1f} {units}' )
   sz.annotate_grids( periodic=False )
   sz.annotate_particles( width=ds.domain_width[2], p_size=80.0, col='y', marker='*' )
   sz.annotate_text( (-400, -370), 'Total internal energy = {: >10.2e}'.format( ad.quantities.total_quantity("thermal_energy").in_units("erg") ), coord_system='plot' )
   sz.annotate_text( (-400, -400), 'Total kinetic energy  = {: >10.2e}'.format( ad.quantities.total_quantity("kinetic_energy").in_units("erg") ), coord_system='plot' )
   sz.save( mpl_kwargs={"dpi":dpi} )
   sz.annotate_cell_edges()
   sz.zoom( 16 )
   sz.save( '%s_Slice_z_temperature_zoom-in.png'%ds, mpl_kwargs={"dpi":dpi} )


   # Metal density
   sz = yt.SlicePlot( ds, 'z', 'metal_density', center=center_mode  )
   sz.set_cmap( 'metal_density', colormap )
   sz.set_unit( 'metal_density', 'g/cm**3' )
   sz.set_zlim( 'metal_density', 1.0e-30, 1.0e-24 )
   sz.set_log(  'metal_density', True )
   sz.set_axes_unit( 'pc' )
   sz.annotate_timestamp( time_unit='Myr', corner='upper_right', time_format='t = {time:.1f} {units}' )
   sz.annotate_grids( periodic=False )
   sz.annotate_particles( width=ds.domain_width[2], p_size=80.0, col='y', marker='*' )
   sz.annotate_text( (-400, -370), 'Total gas  metal mass = {: >10.2e}'.format( ad.quantities.total_quantity("cell_metal_mass").in_units("Msun") ), coord_system='plot' )
   sz.annotate_text( (-400, -400), 'Total star metal mass = {: >10.2e}'.format( (ad[("ParMass")]*ad[("ParMetalFrac")]).sum().in_units("Msun") ), coord_system='plot' )
   sz.annotate_text( (-400, -430), 'Total star mass           = {: >10.2e}'.format( ad[("ParMass")].sum().in_units("Msun") ), coord_system='plot' )
   sz.save( mpl_kwargs={"dpi":dpi} )
