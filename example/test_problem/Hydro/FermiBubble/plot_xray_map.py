#====================================================================================================
# Imports
#====================================================================================================
import yt
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
NORMAL_CONST    = 26e26
OBSERVED_ENERGY = 1
SCALE_CR_ENGY   = 1
OFFSET          = 4.5



#====================================================================================================
# Load data
#====================================================================================================
hf = h5py.File("./Projected_FRB_Data_000035_XRay.h5", 'r')
group1 = hf.get("Map")
group2 = hf.get("Info")

b_max = np.array(group2.get('Keys')).tolist()[3]
b_min = np.array(group2.get('Keys')).tolist()[4]
l_max = np.array(group2.get('Keys')).tolist()[5]
l_min = np.array(group2.get('Keys')).tolist()[6]

image = group1.get("ProjectedXray_08_keV")
image = np.array(image[0]) * NORMAL_CONST / OBSERVED_ENERGY * SCALE_CR_ENGY + OFFSET
plot_arr = image
plot_arr = np.roll( plot_arr, int(plot_arr.shape[1]/2), axis=1 )

# meshgrid
x    = np.linspace( l_min, l_max, plot_arr.shape[1] ) * np.pi/180.0
y    = np.linspace( b_min, b_max, plot_arr.shape[0] ) * np.pi/180.0
X, Y = np.meshgrid( x, y )



#====================================================================================================
# Plot
#====================================================================================================
img = plt.imread("cubehelix.png")
img = img[0,:,:]
plot_cmap         = ListedColormap( img )
plot_log_scale    = True
norm_func         = mpl.colors.LogNorm if plot_log_scale else mpl.colors.Normalize
plot_norm         = norm_func( vmin=plot_arr.min(), vmax=plot_arr.max() )
plot_savename     = "XRay_map_0.8keV.png"
plot_x_tick       = np.array( [ -150, -120, -90, -60, -45, -30, -15, 0, 15, 30, 45, 60, 90, 120, 150 ] )
plot_x_tick_label = np.array( [ "10h", "8h", "6h", "4h", "3h", "2h", "1h", "0h", "23h", "22h", "21h", "20h", "18h", "16h", "14h" ] )
plot_y_tick       = np.array( [ -70+10*i for i in range(15) ] )
plot_y_tick_label = np.array( [ "%d$^{\circ}$"%i for i in plot_y_tick ] )

fig_width  = 16.0 * CM2INCH
fig_height = 11.0 * CM2INCH
fig_dpi    = 400

fig = plt.figure( figsize=(fig_width, fig_height) )

gs = fig.add_gridspec( 1, 2, width_ratios=(20, 1), wspace=0.1 )
ax = fig.add_subplot( gs[0, 0], projection=PROJECT_METHOD )
im = ax.pcolormesh( X, Y, plot_arr, norm=plot_norm, cmap=plot_cmap, shading="auto" )
ax.grid(color='silver', linestyle='-', linewidth=0.5, axis='x')
ax.grid(color='silver', linestyle='-', linewidth=0.5, axis='y')
ax.get_xaxis().set_ticks( plot_x_tick*np.pi/180 )
ax.get_yaxis().set_ticks( plot_y_tick*np.pi/180 )
ax.set_xticklabels( plot_x_tick_label )
ax.set_yticklabels( plot_y_tick_label )
for i in range(len(plot_x_tick)):
    color = "blue" if i > 3 and i < 11 else "white"
    ax.get_xticklabels()[i].set_color( color )

ax   = fig.add_subplot( gs[0, 1] )
cbar = fig.colorbar( im, cax=ax, use_gridspec=True )

plt.suptitle( "Count rate (photons s$^{-1}$deg$^{-2}$)" )

plt.savefig( plot_savename, dpi=fig_dpi )
plt.show()
