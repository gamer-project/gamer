"""
Please edit the section of `Parameters` and `Plot`.
"""
#====================================================================================================
# Imports
#====================================================================================================
import yt
import matplotlib as mpl
import matplotlib.pylab as plt
from matplotlib.colors import LogNorm
import numpy as np
import argparse
import sys


#====================================================================================================
# Constants
#====================================================================================================
WIDTH_RES     = 800
PARTICLE_MASS = 1.6726231e-24 # proton mass in gram
AX_LABEL      = ["x", "y", "z"]
AX_IDX        = {"x":0, "y":1, "z":2}
CM2INCH       = 1/2.54



#====================================================================================================
# Derived fields
#====================================================================================================
def number_density( field, data ):
    return data["Dens"] / PARTICLE_MASS / yt.YTQuantity(1.0, "g")

def pressure_sr( field, data ):
    from yt.units import speed_of_light_cgs, boltzmann_constant_cgs
    unit_v = speed_of_light_cgs

    pres  = data["frame_density"] * data["Temp"] * yt.YTQuantity(1.0, "K")
    pres /= PARTICLE_MASS * yt.YTQuantity(1.0, "g")
    pres /= speed_of_light_cgs**2
    pres *= unit_v**2
    pres *= boltzmann_constant_cgs
    return pres

def pressure_cr( field, data ):
    from yt.units import speed_of_light_cgs
    unit_v = speed_of_light_cgs
    unit_d = data.ds.mass_unit / data.ds.length_unit**3

    NormalizedConst_CRPres = 7

    p_cr  = data["CRay"] / 3 / data["lorentz_factor"]
    p_cr *= unit_d
    p_cr *= unit_v**2
    return p_cr * NormalizedConst_CRPres

def temperature_keV( field, data ):
    from yt.units import boltzmann_constant_cgs
    erg2eV = 6.242e+11

    temp  = data["Temp"] * boltzmann_constant_cgs * yt.YTQuantity(1.0, "K")
    temp *= erg2eV
    temp /= 1e3
    return temp



#====================================================================================================
# Functions
#====================================================================================================
def get_slice_frb_data( ds, slice_axis, plot_var, center, width, height ):
    field      = plot_var["field"]
    field_unit = plot_var["var_unit"]
    field_min  = plot_var["var_min"]

    ax_0_idx = (AX_IDX[slice_axis]+0) % 3
    ax_1_idx = (AX_IDX[slice_axis]+1) % 3
    ax_2_idx = (AX_IDX[slice_axis]+2) % 3

    sli      = yt.SlicePlot( ds, slice_axis, field, center=center )
    frb_arr  = sli.data_source.to_frb( width, (WIDTH_RES, height/width*WIDTH_RES), height=height )
    plot_arr = np.array( frb_arr[field].to( field_unit ) )
    if field_min != None: plot_arr = np.clip( plot_arr, field_min, None ) # apply the minimum allowed value

    plot_center = center - ds.domain_width/2
    plot_extend = np.array( [plot_center[ax_1_idx] - plot_width /2,
                             plot_center[ax_1_idx] + plot_width /2,
                             plot_center[ax_2_idx] - plot_height/2,
                             plot_center[ax_2_idx] + plot_height/2] )
    return plot_arr, plot_center, plot_extend



#====================================================================================================
# Main
#====================================================================================================
parser = argparse.ArgumentParser( description='Plot gas pressure, number density, temperature, and cosmic rays pressure slice for Fermi bubble' )

parser.add_argument( '-s', action='store', required=True,  type=int, dest='idx_start',
                     help='first data index' )
parser.add_argument( '-e', action='store', required=True,  type=int, dest='idx_end',
                     help='last data index' )
parser.add_argument( '-d', action='store', required=False, type=int, dest='didx',
                     help='delta data index [%(default)d]', default=1 )
parser.add_argument( '-i', action='store', required=False, type=str, dest='prefix',
                     help='data path prefix [%(default)s]', default='./' )

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

for i in range(idx_start, idx_end+1, didx):
    ds = yt.load( prefix + "Data_%06d"%i )
    ds.add_field( name=("gamer", "Pres_SR"),  function=pressure_sr,     sampling_type="local", units="dyne/cm**2" )
    ds.add_field( name=("gamer", "P_CR"),     function=pressure_cr,     sampling_type="local", units="dyne/cm**2" )
    ds.add_field( name=("gamer", "num_dens"), function=number_density,  sampling_type="local", units="1/cm**3"    )
    ds.add_field( name=("gamer", "Temp_keV"), function=temperature_keV, sampling_type="local", units="keV"        )

    slice_axis = "x"
    fig_width  = 28.0 * CM2INCH
    fig_height = 11.0 * CM2INCH
    fig_dpi    = 400

    ax_0_idx = (AX_IDX[slice_axis]+0) % 3
    ax_1_idx = (AX_IDX[slice_axis]+1) % 3
    ax_2_idx = (AX_IDX[slice_axis]+2) % 3

    plot_vars = [ { "field":("gamer", "Pres_SR"),  "cmap":"plasma",        "log_scale":True, "var_unit":"erg/cm**3", "title":"P"     , "var_min":None      },
                  { "field":("gamer", "Temp_keV"), "cmap":"afmhot",        "log_scale":True, "var_unit":"keV",       "title":"k_BT"  , "var_min":None      },
                  { "field":("gamer", "num_dens"), "cmap":"gist_earth",    "log_scale":True, "var_unit":"1/cm**3",   "title":"n"     , "var_min":None      },
                  { "field":("gamer", "P_CR"),     "cmap":"nipy_spectral", "log_scale":True, "var_unit":"erg/cm**3", "title":"P_{CR}", "var_min":1.111e-18 } ]
    plot_savename  = "Fermi_Bubble_Slice_%06d.png"%(i)
    plot_axis_unit = "kpc"
    plot_time_unit = "Myr"
    plot_width     = ds.domain_width[ax_1_idx] # yt.YTQuantity(2, 'kpc')
    plot_height    = ds.domain_width[ax_2_idx] # yt.YTQuantity(4, 'kpc')

    N_vars = len( plot_vars )
    fig, ax = plt.subplots( 1, N_vars, figsize=(fig_width, fig_height), sharey=True )

    for v, var in enumerate(plot_vars):
        plot_cmap      = plot_vars[v]["cmap"]
        plot_log_scale = plot_vars[v]["log_scale"]
        plot_var_unit  = plot_vars[v]["var_unit"]
        plot_title     = plot_vars[v]["title"]

        plot_arr, plot_center, plot_extend = get_slice_frb_data( ds, slice_axis, var, ds.domain_center, plot_width, plot_height )

        norm_func   = mpl.colors.LogNorm if plot_log_scale else mpl.colors.Normalize
        plot_extend = plot_extend * ds.length_unit.to( plot_axis_unit )
        plot_norm   = norm_func( vmin=plot_arr.min(), vmax=plot_arr.max() )
        plot_time   = ds.current_time * ds.time_unit.to( plot_time_unit )

        im1 = ax[v].imshow( plot_arr, extent=plot_extend, norm=plot_norm, cmap=mpl.colormaps[plot_cmap] )
        cbar = fig.colorbar( im1, ax=ax[v], pad=0.0 )

        ax[v].tick_params( which="both", direction="in", color="w" )
        cbar.ax.tick_params( which="both", direction="in" )

        ax[v].set( title="$%s$ (%s)"%(plot_title, plot_var_unit) )
        ax[v].set( xlabel="$%s$ (%s)"%(AX_LABEL[ax_1_idx], plot_axis_unit) )

    ax[0].set( ylabel="$%s$ (%s)"%(AX_LABEL[ax_2_idx], plot_axis_unit) )
    plt.suptitle( "%.2f (%s)"%(plot_time, plot_time_unit) )
    plt.tight_layout()
    plt.savefig( plot_savename, dpi=fig_dpi )
    plt.close()
