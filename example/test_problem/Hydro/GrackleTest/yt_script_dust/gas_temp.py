import argparse
import h5py
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import LogLocator
from matplotlib.ticker import FixedLocator, FixedFormatter
import os

# load the command-line parameters
parser = argparse.ArgumentParser(description="Plot temperature evolution from GAMER/Grackle HDF5 outputs.")
parser.add_argument("-s", type=int, required=True, help="Starting index")
parser.add_argument("-e", type=int, required=True, help="Ending index")
parser.add_argument("-d", type=int, required=True, help="Index step")
args = parser.parse_args()


# Configuration: output filenames
FILEOUT    = "fig__GasTemp_plot"
FIG_NAME   = "Gas Temperature v.s Time"

PREFIX     = '../'
MARKERSIZE = 5.0
LINE_WIDTH = 1
DPI        = 150
MYR_IN_SEC = 3.15576e13   # 1 Myr in seconds

# Load K_MYR_INV from the simulation's HDF5 output to ensure consistency
first_file = os.path.join(PREFIX, 'Data_%06d' % args.s)
with h5py.File(first_file, "r") as f:
    K_MYR_INV = f["Info"]["InputTest"]["GrackleTest_ExpCoolCoeff"][()]

print("k [Myr^-1] =", K_MYR_INV)
if K_MYR_INV != 0:
    T_COOL_MYR = 1.0 / K_MYR_INV
    print("t_cool [Myr] =", T_COOL_MYR)


# Load data
temp_all = []
time_all = []
unit_t   = None
for idx in range(args.s, args.e + 1, args.d):
    file_path = os.path.join(PREFIX, 'Data_%06d'%idx)
    if not os.path.isfile(file_path):
        continue
    with h5py.File(file_path, "r") as f:
        temp = f["GridData"]["Temp"][0][0][0][0]
        time = f["Info"]["KeyInfo"]["Time"][0]
        if unit_t is None:
            unit_t = float(f["Info"]["InputPara"]["Unit_T"])  # code unit -> seconds
        temp_all.append(temp)
        time_all.append(time)

# Convert time from code units to Myr, since K_MYR_INV is defined in Myr^-1
time_all = np.array(time_all) * unit_t / MYR_IN_SEC
if K_MYR_INV == 0:
    time_cool = time_all
else:
    time_cool = time_all / T_COOL_MYR

# Plot
# --------------------------------------------
fig, ax = plt.subplots(1, 1)
fig.subplots_adjust(wspace=0.4)

ax.set_title(FIG_NAME)
ax.set_xlabel(r"$t/t_{\rm cool}$", fontsize="large")
ax.plot(time_cool, temp_all, 'ro', lw=LINE_WIDTH, mec='none', ms=MARKERSIZE, label="Numerical")

T0 = float(temp_all[0])
if K_MYR_INV == 0:
    T_ref = np.full_like(time_all, T0, dtype=float)
    ax.plot(time_cool, T_ref, 'b--', label="Reference")
else:
    T_ref = T0 * np.exp(-K_MYR_INV * time_all)
    ax.plot(time_cool, T_ref, 'b--', label="Reference")

# Axis settings
ax.set_xlim(0, time_cool[-1]*1.05)
ax.set_yscale('log')
ax.set_ylim(1.0e4, 2.0e6)

yticks = [1e5, 2e5, 3e5, 5e5, 1e6, 2e6]
ytick_labels = [
    r"$10^5$", r"$2\times10^5$", r"$3\times10^5$",
    r"$5\times10^5$", r"$10^6$", r"$2\times10^6$"
]

ax.yaxis.set_major_locator(FixedLocator(yticks))
ax.yaxis.set_major_formatter(FixedFormatter(ytick_labels))
ax.set_ylabel("$\\mathrm{T\\ [K]}$", fontsize='large')
ax.legend()

# Save outputs
plt.savefig(FILEOUT + ".png", bbox_inches='tight', pad_inches=0.05, dpi=DPI)
# plt.show()
