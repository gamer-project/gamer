import numpy as np
import matplotlib.pylab as plt

# target files
DATA_ID   = 100
data      = np.loadtxt("Xline_y0.000_z0.000_%06d"%DATA_ID)
anal      = np.loadtxt("ReferenceSolution/TM_EoS/SRHD_UR")

# data variables
X         = data[:, 3 ]
Y         = data[:, 4 ]
Z         = data[:, 5 ]
DENS      = data[:, 6 ]
MOMX      = data[:, 7 ]
MOMY      = data[:, 8 ]
MOMZ      = data[:, 9 ]
ENGY      = data[:, 10]
PRES      = data[:, 11]
L_F       = data[:, 12]
VX        = data[:, 13]
VY        = data[:, 14]
VZ        = data[:, 15]

# derived variables
RHO       = DENS/L_F
UX        = VX*L_F
UY        = VY*L_F
UZ        = VZ*L_F

# analytical variables
X_ANAL    = anal[:, 0 ]
RHO_ANAL  = anal[:, 1 ]
UX_ANAL   = anal[:, 2 ]
UY_ANAL   = anal[:, 3 ]
UZ_ANAL   = anal[:, 4 ]
PRES_ANAL = anal[:, 5 ]

# plot
fig, ax = plt.subplots(3, 1, figsize=(15, 10))
plt.subplots_adjust(hspace=0.1, wspace=0)

ax[0].plot(    X_ANAL, PRES_ANAL,          color='C3',                         label="Analytical" )
ax[0].scatter( X,      PRES,         s=15, facecolors='none', edgecolors='C0', label="Simulation" )
ax[0].set(yscale='log')
ax[0].set(ylabel=r"$P$")
ax[0].legend()

ax[1].plot(    X_ANAL, RHO_ANAL,           color='C3',                         label="Analytical" )
ax[1].scatter( X,      RHO,          s=15, facecolors='none', edgecolors='C0', label="Simulation" )
ax[1].set(yscale='log')
ax[1].set(ylabel=r"$\rho$")
ax[1].legend()

ax[2].plot(    X_ANAL, UX_ANAL/1.0e6,      color='C3',                         label="Analytical" )
ax[2].scatter( X,      UX/1.0e6,     s=15, facecolors='none', edgecolors='C0', label="Simulation" )
ax[2].set(ylabel=r"$U_x/c\ [10^6]$")
ax[2].set(xlabel=r"$x$")
ax[2].legend()

plt.tight_layout()

plt.savefig("Fig__Riemann_SRHD_%06d.png"%DATA_ID)
plt.close()
