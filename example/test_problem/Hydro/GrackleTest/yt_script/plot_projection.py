import argparse
import sys
import yt
import numpy as np

# load the command-line parameters
parser = argparse.ArgumentParser( description='Plot the gas projections' )

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


yt.enable_parallelism()

ts = yt.DatasetSeries( [ prefix+'/Data_%06d'%idx for idx in range(idx_start, idx_end+1, didx) ] )


# target fields and their colorbar limits
zlim = {
         'density'                        : (  1.0e-29,  1.0e-21 ),
         'temperature'                    : (  1.0e+00,  1.0e+08 ),
         'thermal_energy_density'         : (  1.0e-21,  1.0e-05 ),
         'sound_speed'                    : (  1.0e+04,  1.0e+08 ),
         'electron_density'               : (  1.0e-29,  1.0e-21 ),
         'DI_density'                     : (  1.0e-29,  1.0e-21 ),
         'DII_density'                    : (  1.0e-29,  1.0e-21 ),
         'H2I_density'                    : (  1.0e-29,  1.0e-21 ),
         'H2II_density'                   : (  1.0e-29,  1.0e-21 ),
         'HDI_density'                    : (  1.0e-29,  1.0e-21 ),
         'HI_density'                     : (  1.0e-29,  1.0e-21 ),
         'HII_density'                    : (  1.0e-29,  1.0e-21 ),
         'HM_density'                     : (  1.0e-29,  1.0e-21 ),
         'HeI_density'                    : (  1.0e-29,  1.0e-21 ),
         'HeII_density'                   : (  1.0e-29,  1.0e-21 ),
         'HeIII_density'                  : (  1.0e-29,  1.0e-21 ),
         'metal_density'                  : (  1.0e-29,  1.0e-21 ),
         'grackle_T_over_mu'              : (  1.0e+00,  1.0e+08 ),
         'grackle_temperature'            : (  1.0e+00,  1.0e+08 ),
         'grackle_mu'                     : (  5.0e-01,  1.5e+00 ),
         'grackle_cooling_rate'           : ( -1.0e-24,  1.0e-24 ),
         'grackle_cooling_strength'       : ( -1.0e+03,  1.0e+03 ),
         'grackle_cooling_length'         : (  1.0e-04,  1.0e+00 ),
         'grackle_cooling_length_over_dh' : (  5.0e-02,  2.0e+01 ),
       }


def set_derived_fields( ds ):

    def _grackle_T_over_mu( field, data ):
        return data['grackle_temperature'] / data['grackle_mu']
    ds.add_field( ('gas', 'grackle_T_over_mu'), function=_grackle_T_over_mu, sampling_type='cell', units='K', display_name='Grackle $T/\mu$'  )

    def _grackle_cooling_rate( field, data ):
        return data['thermal_energy_density'] / data['grackle_cooling_time']
    ds.add_field( ('gas', 'grackle_cooling_rate'), function=_grackle_cooling_rate, sampling_type='cell', units='erg/cm**3/s' )

    def _grackle_cooling_strength( field, data ):
        return 1.0 / data['grackle_cooling_time']
    ds.add_field( ('gas', 'grackle_cooling_strength'), function=_grackle_cooling_strength, sampling_type='cell', units='Myr**-1' )

    def _grackle_cooling_length( field, data ):
        factor = np.where( data['grackle_cooling_time'] < 0.0, 1.0, np.inf )
        return factor * np.abs( data['grackle_cooling_time'] ) * data['sound_speed']
    ds.add_field( ('gas', 'grackle_cooling_length'), function=_grackle_cooling_length, sampling_type='cell', units='kpc' )

    def _grackle_cooling_length_over_dh( field, data ):
        return data['grackle_cooling_length'] / data['dx']
    ds.add_field( ('gas', 'grackle_cooling_length_over_dh'), function=_grackle_cooling_length_over_dh, sampling_type='cell', units='dimensionless', display_name='Grackle Cooling Length / dh' )


for ds in ts.piter():

    set_derived_fields( ds )

    ad = ds.all_data()

    # plot the projections
    for field in zlim.keys():

        if ds.parameters["Grackle_Primordial"] < 3 and (field == 'DI_density' or field == 'DII_density' or field == 'HDI_density'):
            continue
        if ds.parameters["Grackle_Primordial"] < 2 and (field == 'H2I_density' or field == 'H2II_density' or field == 'HM_density'):
            continue
        if ds.parameters["Grackle_Primordial"] < 1 and (field == 'HI_density' or field == 'HII_density' or field == 'HeI_density' or field == 'HeII_density' or field == 'HeIII_density' or field == 'electron_density'):
            continue
        if ds.parameters["Grackle_Metal"] == 0 and field == 'metal_density':
            continue

        # Projection weighted by 1 -> get the average
        pz = yt.ProjectionPlot( ds, 'z', field, weight_field=('index','ones'), center='c', data_source=ad )
        pz.set_axes_unit( 'kpc' )
        pz.set_zlim( field, zlim[field][0], zlim[field][1] )

        pz.annotate_text( (0.1,0.9),
                          "max = {: >11.3e}\nmin  = {: >11.3e}".format(ad[field].max(),ad[field].min()),
                          coord_system='axis', text_args={"color": "grey"} )

        pz.annotate_timestamp( time_unit='Myr', corner='upper_right' )

        pz.save( 'fig_%s_Projection_z_%s.png'%(ds, field), mpl_kwargs={'dpi':150} )

        if field == 'grackle_cooling_length' or field == 'grackle_cooling_length_over_dh':
            pz.annotate_grids()
            pz.save( 'fig_%s_Projection_z_%s_withgrids.png'%(ds, field), mpl_kwargs={'dpi':150} )
