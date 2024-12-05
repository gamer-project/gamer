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
LINE_RES      = 8192
PARTICLE_MASS = 1.6726231e-24 # proton mass in gram
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

def temperature_erg( field, data ):
    from yt.units import boltzmann_constant_cgs
    erg2eV = 6.242e+11

    temp  = data["Temp"] * yt.YTQuantity(1.0, "K")
    temp *= boltzmann_constant_cgs
    return temp



#====================================================================================================
# Functions
#====================================================================================================
def get_line_data( ds, start, end, plot_var ):
    field      = plot_var["field"]
    field_unit = plot_var["var_unit"]

    line      = yt.LineBuffer( ds, start, end, LINE_RES )
    plot_arr  = np.array( line[field].to( field_unit ) )
    plot_r    = np.sqrt( (line[("index", "x")]-ds.domain_center[0])**2 +
                         (line[("index", "y")]-ds.domain_center[1])**2 +
                         (line[("index", "z")]-ds.domain_center[2])**2 )
    return plot_arr, plot_r



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

for i in range(idx_start, idx_end+1, didx):
    ds = yt.load( prefix + "Data_%06d"%i )
    ds.add_field( name=("gamer", "Pres_SR"),  function=pressure_sr,     sampling_type="local", units="dyne/cm**2" )
    ds.add_field( name=("gamer", "P_CR"),     function=pressure_cr,     sampling_type="local", units="dyne/cm**2" )
    ds.add_field( name=("gamer", "num_dens"), function=number_density,  sampling_type="local", units="1/cm**3"    )
    ds.add_field( name=("gamer", "Temp_erg"), function=temperature_erg, sampling_type="local", units="erg"        )

    fig_width  = 22.0 * CM2INCH
    fig_height = 28.0 * CM2INCH
    fig_dpi    = 400

    start      = [7, 7, 0]
    end        = [7, 7, 14]
    start_zoom = [7, 7, 9.833333]
    end_zoom   = [7, 7, 10.66666]

    plot_vars = [ { "field":("gamer", "Pres_SR"),  "cmap":"plasma",        "log_scale":True, "var_unit":"erg/cm**3", "title":"P"     , "var_min":None      },
                  { "field":("gamer", "Temp_erg"), "cmap":"afmhot",        "log_scale":True, "var_unit":"keV",       "title":"k_BT"  , "var_min":None      },
                  { "field":("gamer", "num_dens"), "cmap":"gist_earth",    "log_scale":True, "var_unit":"1/cm**3",   "title":"n"     , "var_min":None      } ]
    plot_savename  = "Fermi_Bubble_Profile_%06d.png"%(i)
    plot_axis_unit = "kpc"
    plot_time_unit = "Myr"

    N_vars = len( plot_vars )
    fig, ax = plt.subplots( N_vars, 2, figsize=(fig_width, fig_height), sharex="col" )

    for v, var in enumerate(plot_vars):
        plot_log_scale   = plot_vars[v]["log_scale"]
        plot_var_unit    = plot_vars[v]["var_unit"]
        plot_title       = plot_vars[v]["title"]
        plot_arr, plot_r = get_line_data( ds, start, end, var )
        plot_time        = ds.current_time * ds.time_unit.to( plot_time_unit )

        ax[v, 0].scatter( plot_r, plot_arr, s=5 )

        ax[v, 0].tick_params( which="both", direction="in", top=True, right=True )

        ax[v, 0].set( ylabel=r"$%s$ (%s)"%(plot_title, plot_var_unit) )
        if plot_log_scale: ax[v, 0].set( yscale="log" )

    # zoom-in
    for v, var in enumerate(plot_vars):
        plot_log_scale   = plot_vars[v]["log_scale"]
        plot_var_unit    = plot_vars[v]["var_unit"]
        plot_title       = plot_vars[v]["title"]
        plot_arr, plot_r = get_line_data( ds, start_zoom, end_zoom, var )
        plot_time        = ds.current_time * ds.time_unit.to( plot_time_unit )

        ax[v, 1].scatter( plot_r, plot_arr, s=5 )

        ax[v, 1].tick_params( which="both", direction="in", top=True, right=True, labelleft=False, labelright=True )

        if plot_log_scale: ax[v, 1].set( yscale="log" )

    ax[N_vars-1, 0].set( xlabel=r"$%s$ (%s)"%("r", plot_axis_unit) )
    ax[N_vars-1, 1].set( xlabel=r"$%s$ (%s)"%("r", plot_axis_unit) )
    plt.suptitle( "%.2f (%s)"%(plot_time, plot_time_unit) )
    plt.tight_layout()
    plt.savefig( plot_savename, dpi=fig_dpi )
    plt.close()
