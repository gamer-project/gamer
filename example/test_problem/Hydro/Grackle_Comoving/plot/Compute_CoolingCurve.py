#!/bin/python

# reconstruct the cooling rate from the a_scale evolution of the temperature
# and compare with the cooling rate table
import numpy as np
import matplotlib.pyplot as plt
import glob


# Load the data
data = np.loadtxt('../Record__GrackleComoving', skiprows=9)

# data columns
a_scale       = data[:,0] # scale factor
dt            = data[:,2] # sec
nden          = data[:,3] # cm^-3
mu            = data[:,4] # mean molecular weight
temp          = data[:,5] # K
eden          = data[:,6] # erg/cm^3
lcool_grackle = data[:,7] # erg/cm^3/s (cooling rate computed by Grackle)


# reconstruct the cooling rate from Record__GrackleComoving
dedt = np.zeros_like(a_scale)
dadt = np.zeros_like(a_scale)
for i in range(1, len(a_scale)-2):
    dedt[i] = (eden   [i+1] - eden   [i-1]) / (dt[i] + dt[i+1])
    dadt[i] = (a_scale[i+1] - a_scale[i-1]) / (dt[i] + dt[i+1])


dedt_com = eden * (-5 / a_scale * dadt) # loss due to cosmic expansion
dedt_rad = dedt - dedt_com              # loss due to radiation

L_cool = (-dedt_rad) / nden**2

fig, ax = plt.subplots(1, 2, figsize=(10, 5))


# cooling rate table
tables = glob.glob('../coolingrate_z*.dat')
redshifts = [float(table.split('_')[-1].split('.')[0][1:]) for table in tables]
redshifts.sort()
for z in redshifts:
    cr = np.loadtxt('../coolingrate_z{}.dat'.format(z))
    ax[0].plot(cr[:,0], cr[:,1], '--', lw=0.5, label='z={}'.format(z))


# cooling rate computed by Grackle during the simulation
ax[0].plot(temp, lcool_grackle, label='Grackle', ms=0, ls='dashed', lw=1, color='gray')


# result from GAMER
# --> skip the first and last two data points to avoid the boundary effect
for z in redshifts:
    a        = 1 / (1 + z)
    idx0     = np.where(a_scale <= a)[0][-1]
    idx1     = idx0 + 1
    if idx0 == 0:
        ax[0].scatter(temp[ 1], L_cool[ 1], s=3, label='GAMER (reconstructed) (z={})'.format(z), zorder=10)
    elif idx0 > len(a_scale)-3:
        ax[0].scatter(temp[-3], L_cool[-3], s=3, label='GAMER (reconstructed) (z={})'.format(z), zorder=10)
    else:
        a0       = a_scale[idx0]
        a1       = a_scale[idx1]
        temp_a   = (temp  [idx0] * (a1 - a) + temp  [idx1] * (a - a0)) / (a1 - a0)
        L_cool_a = (L_cool[idx0] * (a1 - a) + L_cool[idx1] * (a - a0)) / (a1 - a0)
        ax[0].scatter(temp_a, L_cool_a, s=3, label='GAMER (reconstructed) (z={})'.format(z), zorder=10)

# ax[0].plot(temp, L_cool, label='GAMER (reconstructed)', ms=0, ls='solid', lw=0.5, color='black')
ax[0].plot(temp[1:-2], L_cool[1:-2], label='GAMER (reconstructed)', ms=0, ls='solid', lw=0.5, color='black')

ax[0].set_xlabel('Temperature [K]')
ax[0].set_ylabel(r'$\Lambda$ [erg cm$^3$ s$^{-1}$]')
ax[0].set_xscale('log')
ax[0].set_yscale('log')
ax[0].set_xlim(1e4, 2e8)
ax[0].set_ylim(3e-25, 2e-21)
ax[0].legend(ncol=2, fontsize=6)

ax[1].plot(a_scale, temp, 'r-')
ax[1].set_xlabel('a')
ax[1].set_ylabel('Temperature [K]')
ax[1].set_yscale('log')


plt.savefig('CoolingCurve.png', bbox_inches='tight', dpi=300)
