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

        # projection weighted by 1 -> get the average
        pz = yt.ProjectionPlot( ds, 'z', field, weight_field=('index','ones'), center='c', data_source=ad )
        pz.set_axes_unit( 'kpc' )
        pz.set_zlim( field, zlim[field][0], zlim[field][1] )

        pz.annotate_text( (0.1,0.9),
                          "max = {: >11.3e}\nmin  = {: >11.3e}".format(ad[field].max(),ad[field].min()),
                          coord_system='axis', text_args={"color": "grey"} )

        pz.annotate_timestamp( time_unit='Myr', corner='upper_right' )
        if field == 'grackle_cooling_length' or field == 'grackle_cooling_length_over_dh':
            pz.annotate_grids()

        # gas IC is uniform if either density and temperature vary by less than a factor of 5
        dens_ratio = ad["density"].max() / ad["density"].min()
        temp_ratio = ad["grackle_T_over_mu"].max() / ad["grackle_T_over_mu"].min()
        dataset_is_uniform = (dens_ratio <= 5.0) or (temp_ratio <= 5.0)
        
        pz._setup_plots()
        if not dataset_is_uniform: # secondary axis labels and ticks for gas properties
            # 1. grab the actual data limits and plot limits
            dens_min   = ad["density"].min().in_units("g/cm**3").v # .v gets the raw float value
            dens_max   = ad["density"].max().in_units("g/cm**3").v
            t_mu_min   = ad["grackle_T_over_mu"].min().in_units("K").v
            t_mu_max   = ad["grackle_T_over_mu"].max().in_units("K").v
            plot_axes  = pz.plots[field].axes
            xmin, xmax = plot_axes.get_xlim()
            ymin, ymax = plot_axes.get_ylim()

            # 2. axis mapping functions
            def kpc_to_logdens(x):
                return np.log10(dens_min) + (x - xmin) * (np.log10(dens_max) - np.log10(dens_min)) / (xmax - xmin)
            def logdens_to_kpc(x):
                return xmin + (x - np.log10(dens_min)) * (xmax - xmin) / (np.log10(dens_max) - np.log10(dens_min))
            def kpc_to_logT(y):
                return np.log10(t_mu_min) + (y - ymin) * (np.log10(t_mu_max) - np.log10(t_mu_min)) / (ymax - ymin)
            def logT_to_kpc(y):
                return ymin + (y - np.log10(t_mu_min)) * (ymax - ymin) / (np.log10(t_mu_max) - np.log10(t_mu_min))

            # 3. apply changes to secondary axes
            ax_top     = plot_axes.secondary_xaxis('top', functions=(kpc_to_logdens, logdens_to_kpc))
            ax_right   = plot_axes.secondary_yaxis('right', functions=(kpc_to_logT, logT_to_kpc))

            # 3-1. grab the full font objects from the primary axis
            primary_label_font = plot_axes.xaxis.get_label().get_fontproperties()
            primary_tick_font  = plot_axes.xaxis.get_ticklabels()[0].get_fontproperties()

            # 3-2. style top x-axis (gas density)
            ax_top.set_xlabel(r'$\log_{10} (\rho_{\mathrm{gas}} \ [\mathrm{g/cm^3}])$', labelpad=12)
            ax_top.xaxis.get_label().set_fontproperties(primary_label_font) # Force exact font match
            ax_top.tick_params(axis='x', direction='out', which='both')
            for label in ax_top.get_xticklabels():
                label.set_fontproperties(primary_tick_font)

            # 3-3. style right y-axis (T/mu)
            ax_right.set_ylabel(r'$\log_{10} (T/\mu \ [\mathrm{K}])$', labelpad=12)
            ax_right.yaxis.get_label().set_fontproperties(primary_label_font) # Force exact font match
            ax_right.tick_params(axis='y', direction='out', which='both')
            for label in ax_right.get_yticklabels():
                label.set_fontproperties(primary_tick_font)

            # 3-4. sync the math font
            ax_top.xaxis.get_label().set_math_fontfamily(plot_axes.xaxis.label.get_math_fontfamily())
            ax_right.yaxis.get_label().set_math_fontfamily(plot_axes.xaxis.label.get_math_fontfamily())

            # 4. shift colorbar rightwards for secondary y-axis and save images
            colorbar_axes = pz.plots[field].cb.ax
            pos           = colorbar_axes.get_position()
            new_pos       = [pos.x0 + 0.08, pos.y0, pos.width, pos.height]
            colorbar_axes.set_position(new_pos)
            fig           = pz.plots[field].figure
            #fig           = plot_axes.get_figure()
            fig.savefig('fig_%s_Projection_z_%s.png' % (ds, field), dpi=150, bbox_inches='tight')
        else:
            pz.save( 'fig_%s_Projection_z_%s.png'%(ds, field), mpl_kwargs={'dpi':150} )
