import argparse
import sys
import yt
import numpy as np

# load the command-line parameters
parser = argparse.ArgumentParser( description='Plot the gas and star projections' )

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


yt.enable_parallelism()

ts = yt.DatasetSeries( [ prefix+'/Data_%06d'%idx for idx in range(idx_start, idx_end+1, didx) ] )


# target fields and their units and limits
fields   = [ 'density',
             'temperature',
             'thermal_energy_density',
             'number_density',
             'sound_speed',
             'pressure',
             'jeans_length',
             'jeans_mass',
             'particle_density_on_grid',
           ]

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

   # define the derived fields

   # lambda_J = sqrt( (pi*Cs^2) / (G*rho) )
   def _jeans_length( field, data ):
      return np.sqrt( ( np.pi * data["sound_speed"]**2 )/( data.ds.units.newtons_constant * data["density"] ) )

   ds.add_field( ("gas", "jeans_length"), function=_jeans_length, sampling_type="cell", units="cm" )


   # get the corresponding temerature from the density for the given Jeans length
   def get_T_for_MaxJeansLength( rho ):
      lambda_J_max = ds.parameters['SF_CreateStar_MaxGasJeansL'] * ds.domain_width[0] / ds.domain_dimensions[0]
      T            = rho * lambda_J_max**2 * ( ds.units.newtons_constant*ds.mu*ds.units.hydrogen_mass )/( ds.units.boltzmann_constant * np.pi * ds.gamma )
      return T


   # plot the projections
   for field in fields:

      ad = ds.all_data()
      ad.min_level = 0
      ad.max_level = 0

      # Projection weighted by 1 -> get the average; assume no refined level
      pz = yt.ProjectionPlot( ds, 'z', field, weight_field=('index','ones'), center='c', data_source=ad )
      pz.set_axes_unit( 'kpc' )
      pz.set_unit( field, units[field] )
      pz.set_zlim( field, zlim_min[field], zlim_max[field] )
      pz.set_cmap( field, 'algae' )
      pz.annotate_timestamp( time_unit='Myr', corner='upper_right' )
      pz.annotate_grids( periodic=False )

      # annotate the formed star particles
      pz.annotate_particles( width=ds.domain_width[2], p_size=10.0, col='m', marker='*' )

      # annotate the line indicated the star formation threshold for the projected star particle density
      if field == 'particle_density_on_grid':

         if ds.parameters['SF_CreateStar_Scheme'] == 1:   # with minimum density threshold

            x_0 = np.log10( ds.arr( ds.parameters['SF_CreateStar_MinGasDens'], 'code_density' ).in_units('g/cm**3').d / zlim_min['density'] ) / np.log10( zlim_max['density'] / zlim_min['density'] ) * ds.domain_width[0]
            y_0 = ds.domain_left_edge[1]
            x_1 = x_0
            y_1 = ds.domain_right_edge[1]

         elif ds.parameters['SF_CreateStar_Scheme'] == 2: # with maximum Jeans length thredshold

            x_0 = ds.domain_left_edge[0]
            y_0 = np.log10( get_T_for_MaxJeansLength( ds.quan( zlim_min['density'], 'g/cm**3' ) ) )
            x_1 = ds.domain_right_edge[0]
            y_1 = np.log10( get_T_for_MaxJeansLength( ds.quan( zlim_max['density'], 'g/cm**3' ) ) )

         else:
            raise RuntimeError('Unsupported SF_CreateStar_Scheme = %d !!'%(ds.parameters['SF_CreateStar_Scheme']))

         pz.annotate_line( [ x_0, y_0,  ds.domain_right_edge[2] ],
                           [ x_1, y_1,  ds.domain_right_edge[2] ],
                           coord_system="data", color='red', linewidth=3, linestyle='--' )

      # save the image
      pz.save( mpl_kwargs={'dpi':150} )

