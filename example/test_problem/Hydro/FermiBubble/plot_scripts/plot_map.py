#====================================================================================================
# Imports
#====================================================================================================
import h5py
import matplotlib as mpl
import matplotlib.pylab as plt
from matplotlib.colors import LogNorm, ListedColormap
import numpy as np
import argparse
import sys



#====================================================================================================
# Constants
#====================================================================================================
CM2INCH         = 1/2.54
PROJECT_METHOD  = "hammer"



#====================================================================================================
# Parser
#====================================================================================================
parser = argparse.ArgumentParser( description='Plot perspective projection map' )

parser.add_argument( '-t', action='store', required=True, type=str, dest='plot_type',
                     choices=["x_ray", "gamma_ray"], help='Plot map type.' )

args=parser.parse_args()

# take note
print( '\nCommand-line arguments:' )
print( '-------------------------------------------------------------------' )
print( ' '.join(map(str, sys.argv)) )
print( '-------------------------------------------------------------------\n' )


#====================================================================================================
# Parameters
#====================================================================================================
plot_type = args.plot_type
plot_args = {}
if plot_type == "x_ray":
    plot_args["normal_const"]           = 26e26
    plot_args["observed_energy"]        = 1
    plot_args["scale_cr_engy"]          = 1
    plot_args["offset"]                 = 4.5
    plot_args["lower_limit"]            = 0.0
    plot_args["filename"]               = "./Projected_FRB_Data_000035_XRay.h5"
    plot_args["b_max_idx"]              = 3
    plot_args["b_min_idx"]              = 4
    plot_args["l_max_idx"]              = 5
    plot_args["l_min_idx"]              = 6
    plot_args["field"]                  = "ProjectedXray_08_keV"
    img = plt.imread("cubehelix.png")
    img = img[0,:,:]
    plot_args["cmap"]                   = ListedColormap( img )
    plot_args["log_scale"]              = True
    plot_args["savename"]               = "XRay_map_0.8keV.png"
    plot_args["title"]                  = r"Count rate (photons s$^{-1}$deg$^{-2}$)"
    plot_args["cbar_max"]               = None
    plot_args["cbar_min"]               = None
    plot_args["xtick_label_dark_range"] = [3, 11]
elif plot_type == "gamma_ray":
    plot_args["normal_const"]           = 1.687e29
    plot_args["observed_energy"]        = 108.61e9
    plot_args["scale_cr_engy"]          = 250
    plot_args["offset"]                 = 0.0
    plot_args["lower_limit"]            = 1.e-13
    plot_args["filename"]               = "./Projected_FRB_Data_000035_100e9_1e6_2.4.h5"
    plot_args["b_max_idx"]              = 7
    plot_args["b_min_idx"]              = 8
    plot_args["l_max_idx"]              = 9
    plot_args["l_min_idx"]              = 10
    plot_args["field"]                  = "ProjectedLeptonicGammaRay"
    plot_args["cmap"]                   = mpl.colormaps["nipy_spectral"]
    plot_args["log_scale"]              = True
    plot_args["savename"]               = "GammaRay_100e9_1e6.png"
    plot_args["title"]                  = r"Photon flux (GeV$^{-1}$cm$^{-2}$s$^{-1}$sr$^{-1}$)"
    plot_args["cbar_max"]               = 1.e-11
    plot_args["cbar_min"]               = 1.e-9
    plot_args["xtick_label_dark_range"] = [5, 9]

fig_width  = 20.0 * CM2INCH
fig_height = 11.0 * CM2INCH
fig_dpi    = 400



#====================================================================================================
# Load data
#====================================================================================================
hf = h5py.File( plot_args["filename"], 'r')
group1 = hf.get("Map")
group2 = hf.get("Info")

b_max = np.array(group2.get('Keys')).tolist()[plot_args["b_max_idx"]]
b_min = np.array(group2.get('Keys')).tolist()[plot_args["b_min_idx"]]
l_max = np.array(group2.get('Keys')).tolist()[plot_args["l_max_idx"]]
l_min = np.array(group2.get('Keys')).tolist()[plot_args["l_min_idx"]]

image = group1.get( plot_args["field"] )
image = np.array(image[0]) * plot_args["normal_const"] / plot_args["observed_energy"] * plot_args["scale_cr_engy"] + plot_args["offset"]
plot_arr = image
plot_arr = np.where( plot_arr > plot_args["lower_limit"], plot_arr, plot_args["lower_limit"]  )
plot_arr = np.roll( plot_arr, int(plot_arr.shape[1]/2), axis=1 )

# meshgrid
x    = np.linspace( l_min, l_max, plot_arr.shape[1] ) * np.pi/180.0
y    = np.linspace( b_min, b_max, plot_arr.shape[0] ) * np.pi/180.0
X, Y = np.meshgrid( x, y )



#====================================================================================================
# Plot
#====================================================================================================
norm_func         = mpl.colors.LogNorm if plot_args["log_scale"] else mpl.colors.Normalize
if plot_args["cbar_min"] is None: plot_args["cbar_min"] = plot_arr.min()
if plot_args["cbar_max"] is None: plot_args["cbar_max"] = plot_arr.max()
plot_norm         = norm_func( vmin=plot_args["cbar_min"], vmax=plot_args["cbar_max"] )
plot_x_tick       = np.array( [ -150, -120, -90, -60, -45, -30, -15, 0, 15, 30, 45, 60, 90, 120, 150 ] )
plot_x_tick_label = np.array( [ "10h", "8h", "6h", "4h", "3h", "2h", "1h", "0h", "23h", "22h", "21h", "20h", "18h", "16h", "14h" ] )
plot_y_tick       = np.array( [ -70+10*i for i in range(15) ] )
plot_y_tick_label = np.array( [ r"%d$^{\circ}$"%i for i in plot_y_tick ] )

fig = plt.figure( figsize=(fig_width, fig_height) )

gs = fig.add_gridspec( 1, 2, width_ratios=(20, 1), wspace=0.1 )
ax = fig.add_subplot( gs[0, 0], projection=PROJECT_METHOD )
im = ax.pcolormesh( X, Y, plot_arr, norm=plot_norm, cmap=plot_args["cmap"], shading="auto" )
ax.grid(color='silver', linestyle='-', linewidth=0.5, axis='x')
ax.grid(color='silver', linestyle='-', linewidth=0.5, axis='y')
ax.get_xaxis().set_ticks( plot_x_tick*np.pi/180 )
ax.get_yaxis().set_ticks( plot_y_tick*np.pi/180 )
ax.set_xticklabels( plot_x_tick_label, rotation=90 )
ax.set_yticklabels( plot_y_tick_label )
for i in range(len(plot_x_tick)):
    color = "black" if i > plot_args["xtick_label_dark_range"][0] and i < plot_args["xtick_label_dark_range"][1] else "white"
    ax.get_xticklabels()[i].set_color( color )

ax   = fig.add_subplot( gs[0, 1] )
cbar = fig.colorbar( im, cax=ax, use_gridspec=True )

plt.suptitle( plot_args["title"] )

plt.savefig( plot_args["savename"], dpi=fig_dpi )
plt.show()
