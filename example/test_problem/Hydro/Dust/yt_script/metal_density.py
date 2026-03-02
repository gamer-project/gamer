import argparse
import h5py
import numpy as np
import matplotlib.pyplot as plt
import os



# load the command-line parameters
parser = argparse.ArgumentParser(description="Plot metal density evolution from GAMER/Grackle HDF5 outputs.")
parser.add_argument("-s", type=int, required=True, help="Starting index")
parser.add_argument("-e", type=int, required=True, help="Ending index")
parser.add_argument("-d", type=int, required=True, help="Index step")
args = parser.parse_args()


# Constants
Const_cm   = 1.0
Const_amu  = 1.660539040e-24  
UNIT_L     = 3.08567758149e21 
UNIT_M     = 1.9885e42        
UNIT_D     = UNIT_M / UNIT_L**3 
UNIT_T     = 3.15569252e13

# Configuration: output filenames
fileout  = "fig__MetalDensity_plot"
fig_name = "Metal Density v.s Time"

prefix      = '../'
line_width = 1
markersize = 5.0
dpi        = 150


# Load data
density_all = []
time_all = []

for idx in range(args.s, args.e + 1, args.d):
    file_path = os.path.join(prefix, f'Data_{idx:06d}')
    if not os.path.isfile(file_path):
        continue
    with h5py.File(file_path, "r") as f:
        density = f["GridData"]["Metal"][0][0][0][0] * UNIT_D / (Const_amu / Const_cm**3)
        time = f["Info"]["KeyInfo"]["Time"][0]
        density_all.append(density)
        time_all.append(time)


# Plot
fig, ax = plt.subplots(1, 1)
fig.subplots_adjust(wspace=0.4)

ax.set_title(fig_name)
ax.set_xlabel("$\\mathrm{t\\ [Myr]}$", fontsize="large")
ax.set_ylabel("$\\mathrm{Metal\\ density\\ [amu/cm^3]}$", fontsize="large")
ax.plot(time_all, density_all, 'r-o', lw=line_width, mec='none', ms=markersize, label="Numerical")
ax.set_xscale('log')
ax.set_xlim(1.0e-2, 20)
ax.yaxis.set_minor_locator(plt.LogLocator(base=10.0, subs=[2.0, 5.0, 8.0]))

# Save outputs
plt.savefig(fileout + ".png", bbox_inches='tight', pad_inches=0.05, dpi=dpi)
# plt.show()
