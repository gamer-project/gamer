import argparse
import h5py
import numpy as np
import matplotlib.pyplot as plt
import os


# load the command-line parameters
parser = argparse.ArgumentParser(description="Plot dust temperature")
parser.add_argument("-s", type=int, required=True, help="Starting index")
parser.add_argument("-e", type=int, required=True, help="Ending index")
parser.add_argument("-d", type=int, required=True, help="Index step")
args = parser.parse_args()

# Constants
Const_cm = 1.0
Const_amu = 1.660539040e-24  
UNIT_L = 3.08567758149e21    
UNIT_M = 1.9885e42           
UNIT_D = UNIT_M / UNIT_L**3  
UNIT_T = 3.15569252e13


# Configuration: output filenames
fileout  = "fig__GasDensity_plot"
fig_name = "Gas Density v.s Time"

prefix      = '../'
markersize = 5.0
line_width = 1
dpi = 150




gas_density_all = []
time_all = []
for idx in range(args.s, args.e+1, args.d):
    f = h5py.File(os.path.join(prefix, 'Data_%06d' % idx))
    gas_density = f["GridData"]["Dens"][0][0][0][0] * UNIT_D / (Const_amu / Const_cm**3)
    time = f["Info"]["KeyInfo"]["Time"][0]
    gas_density_all.append(gas_density)
    time_all.append(time)


#  plot
f, ax = plt.subplots( 1, 1 )
f.subplots_adjust( wspace=0.4 )
[ f.axes[t].set_xlim( 0.0, 250.0 ) for t in range(0,1,1) ]


ax.set_xlabel( "$\mathrm{t\ [Myr]}$", fontsize="large" )
ax.set_title(fig_name)
ax.plot( time_all, gas_density_all, 'r-o', lw=line_width, mec='none', ms=markersize )
# ax.set_yscale('log')
ax.set_xscale('log')
# ax.set_ylim( 1.0e1, 1.0e4 )
ax.set_xlim( 1.0e-2, 20 )
ax.yaxis.set_minor_locator( plt.LogLocator(base=10.0, subs=[2.0,5.0,8.0]) )
ax.set_ylabel('$\mathrm{Gas Density\ [amu/cm^3]}$', fontsize='large' )

#  show/save figure
plt.savefig( fileout+".png", bbox_inches='tight', pad_inches=0.05, dpi=dpi )
#  plt.show()