import argparse
import sys
import yt
import numpy as np

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
print( ' '.join(map(str, sys.argv)) )
print( '-------------------------------------------------------------------\n' )


idx_start   = args.idx_start
idx_end     = args.idx_end
didx        = args.didx
prefix      = '../'

colormap    = 'algae'
center_mode = 'c'
dpi         = 150


yt.enable_parallelism()

ts = yt.DatasetSeries( [ prefix+'/Data_%06d'%idx for idx in range(idx_start, idx_end+1, didx) ] )


fields   = [ 'density', 'temperature', 'thermal_energy_density', 'number_density', 'sound_speed', 'pressure', 'jeans_length', 'jeans_mass', 'particle_density_on_grid' ]

units    = { 'density'                 : 'g/cm**3',
             'temperature'             : 'K',
             'thermal_energy_density'  : 'erg/cm**3',
             'number_density'          : 'cm**-3',
             'sound_speed'             : 'cm/s',
             'pressure'                : 'erg/cm**3',
             'jeans_length'            : 'kpc',
             'jeans_mass'              : 'Msun',
             'particle_density_on_grid': 'Msun/pc**3',
           }

zlim_min = { 'density'                 : 1.0e-29,
             'temperature'             : 1.0e0,
             'thermal_energy_density'  : 1.0e-21,
             'number_density'          : 1.0e-5,
             'sound_speed'             : 1.0e4,
             'pressure'                : 1.0e-21,
             'jeans_length'            : 1.0e-3,
             'jeans_mass'              : 1.0e2,
             'particle_density_on_grid': 1.0e-2,
           }

zlim_max = { 'density'                 : 1.0e-21,
             'temperature'             : 1.0e8,
             'thermal_energy_density'  : 1.0e-5,
             'number_density'          : 1.0e3,
             'sound_speed'             : 1.0e8,
             'pressure'                : 1.0e-5,
             'jeans_length'            : 1.0e5,
             'jeans_mass'              : 1.0e17,
             'particle_density_on_grid': 1.0e2,
           }

for ds in ts.piter():

   def _jeans_length( field, data ):
      return np.sqrt( (np.pi * data["sound_speed"]**2 )/( data.ds.units.newtons_constant * data["density"] ) )

   ds.add_field( ("gas", "jeans_length"), function=_jeans_length, sampling_type="cell", units="cm" )

   for field in fields:

      sz = yt.SlicePlot( ds, 'z', field, center=center_mode )
      sz.set_axes_unit( 'kpc' )
      sz.set_unit( field, units[field] )
      sz.set_zlim( field, zlim_min[field], zlim_max[field] )
      sz.set_cmap( field, colormap )
      sz.annotate_timestamp( time_unit='Myr', corner='upper_right' )
      sz.annotate_grids( periodic=False )
      sz.save( mpl_kwargs={'dpi':dpi} )


   lambda_J_max = ds.parameters['SF_CreateStar_MaxGasJeansL'] * ds.domain_width[0]/ds.domain_dimensions[0]
   rho_cutoff   = ds.quan( 1.0e-21, 'g/cm**3' )
   T_cutoff     = rho_cutoff * lambda_J_max**2 * ( ds.units.newtons_constant*ds.mu*ds.units.hydrogen_mass )/( ds.units.boltzmann_constant * np.pi * ds.gamma )

   print('lambda_J_max = ', lambda_J_max )
   print('rho_cutoff   = ', rho_cutoff   )
   print('T_cutoff     = ', T_cutoff     )
