import argparse
import h5py
import numpy as np
import matplotlib.pyplot as plt
import os

# load the command-line parameters
parser = argparse.ArgumentParser(description="Plot temperature evolution from GAMER/Grackle HDF5 outputs.")
parser.add_argument("-s", type=int, required=True, help="Starting index")
parser.add_argument("-e", type=int, required=True, help="Ending index")
parser.add_argument("-d", type=int, required=True, help="Index step")
parser.add_argument("-option", type=str, required=True, choices=["edot_0", "edot_const"],
                    help="Analytic cooling mode: edot_0 or edot_const")

args = parser.parse_args()


# Configuration: output filenames
fileout  = "fig__GasTemp_plot"
fig_name = "Gas Temperature v.s Time"

prefix      = '../'
markersize  = 5.0
line_width  = 1
dpi         = 150
k_myr_inv   = 0.4


# Load data
temp_all = []
time_all = []

for idx in range(args.s, args.e + 1, args.d):
    file_path = os.path.join(prefix, f'Data_{idx:06d}')
    if not os.path.isfile(file_path):
        continue
    with h5py.File(file_path, "r") as f:
        temp = f["GridData"]["Temp"][0][0][0][0]
        time = f["Info"]["KeyInfo"]["Time"][0]
        temp_all.append(temp)
        time_all.append(time)


# Plot
# --------------------------------------------
fig, ax = plt.subplots(1, 1)
fig.subplots_adjust(wspace=0.4)

ax.set_xlabel("$\\mathrm{t\\ [Myr]}$", fontsize="large")
ax.set_title(fig_name)
ax.plot(time_all, temp_all, 'ro', lw=line_width, mec='none', ms=markersize, label="Numerical")

# Analytic solution
if args.option == "edot_0":
    T0 = float(temp_all[0])
    T_analytic = np.full_like(time_all, T0, dtype=float)
    ax.plot(time_all, T_analytic, 'b--', label="Analytic")
elif args.option == "edot_const":
    T0 = float(temp_all[0])
    time_all = np.array(time_all)
    T_analytic = T0 * np.exp(-k_myr_inv * time_all)
    ax.plot(time_all, T_analytic, 'b--', label="Analytic")

# Axis settings
ax.set_yscale('log')
ax.set_ylim(1.0e2, 2.0e6)
ax.set_xlim(1.0e-1, 20)
ax.yaxis.set_minor_locator(plt.LogLocator(base=10.0, subs=[2.0, 5.0, 8.0]))
ax.set_ylabel("$\\mathrm{T\\ [K]}$", fontsize='large')
ax.legend()

# Save outputs
plt.savefig(fileout + ".png", bbox_inches='tight', pad_inches=0.05, dpi=dpi)
