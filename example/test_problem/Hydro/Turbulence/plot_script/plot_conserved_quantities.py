import argparse
import subprocess
import numpy as np
import matplotlib.pyplot as plt

# Load data
MHD = subprocess.check_output( ["awk", '$1 == "MHD" {print $2}', "../Record__Note"], text=True ).strip() == "ON"

table = np.loadtxt("../Record__Conservation")
time = table[:, 0]
ekin = table[:,32]
emag = table[:,41] if MHD else np.zeros_like(time)
momx = table[:,14]
momy = table[:,17]
momz = table[:,20]

# Plot
f, ax = plt.subplots(3, 1, figsize=(6,10))
f.subplots_adjust(wspace=0.4)
ax[0].plot(time, emag, label = r"$E_{\rm mag}$")
ax[0].set_yscale('log')
ax[0].set_xlim(0, 100)
ax[0].set_ylabel(r'$E_{\rm mag}$', fontsize='large')

ax[1].plot(time[1:], ekin[1:], label = r"$E_{\rm kin}$")
ax[1].set_yscale('log')
ax[1].set_xlim(0, 100)
ax[1].set_ylim(1e-4, 1e-1)
ax[1].set_ylabel(r'$E_{\rm kin}$', fontsize='large')

ax[2].plot(time, momx, label = "MomX")
ax[2].plot(time, momy, label = "MomY")
ax[2].plot(time, momz, label = "MomZ")
ax[2].set_xlim(0, 100)
ax[2].set_xlabel(r"$t$", fontsize="large")
ax[2].set_ylabel(r'Mom', fontsize='large')
ax[2].legend(loc='upper left')
plt.savefig("fig_conserved_quantities.png", bbox_inches='tight', pad_inches=0.05, dpi=150)
plt.close()


