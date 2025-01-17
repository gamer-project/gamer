#!/bin/python

# reconstruct the cooling rate from the time evolution of the temperature
# and compare with the cooling rate table
import numpy as np
import matplotlib.pyplot as plt

# Load the data
data = np.loadtxt('../Record__GrackleComoving', skiprows=2)

time          = data[:,0] # scale factor
dt            = data[:,2] # sec
nden          = data[:,3] # cm^-3
mu            = data[:,4] # mean molecular weight
temp          = data[:,5] # K
eden          = data[:,6] # erg/cm^3
lcool_grackle = data[:,7] # erg/cm^3/s (cooling rate computed by Grackle)

# reconstruct the cooling rate from Record__GrackleComoving
dedt = np.gradient(eden) / dt
dadt = np.gradient(time) / dt

dedt_com = eden * (-5 / time * dadt) # loss due to cosmic expansion
dedt_rad = dedt - dedt_com           # loss due to radiation

L_cool = (-dedt_rad) / nden**2

fig, ax = plt.subplots(1, 2, figsize=(10, 5))
# ax[0].plot(temp, lcool, 'r-')

# cooling rate table
for z in range(0, 10):
    cr = np.loadtxt('../coolingrate_z{:.1f}.dat'.format(z))
    ax[0].plot(cr[:,0], cr[:,1], '--', lw=0.5, label='z={}'.format(z))

# cooling rate computed by Grackle during the simulation
ax[0].plot(temp, lcool_grackle, label='Grackle', ms=0, ls='dashed', lw=1, color='gray')

# result from GAMER
ax[0].scatter(temp[-1], L_cool[-1], s=3, label='GAMER (reconstructed) (z=0)')
for z in range(1, 9):
    a = 1 / (1 + z)
    idx = np.where(time <= a)[0][-1]
    a0 = time[idx]
    a1 = time[idx+1]
    temp_a = (temp[idx] * (a1 - a) + temp[idx+1] * (a - a0)) / (a1 - a0)
    L_cool_a = (L_cool[idx] * (a1 - a) + L_cool[idx+1] * (a - a0)) / (a1 - a0)
    ax[0].scatter(temp_a, L_cool_a, s=3, label='GAMER (reconstructed) (z={})'.format(z), zorder=10)
ax[0].scatter(temp[0], L_cool[0], s=3, label='GAMER (reconstructed) (z=9)', zorder=10)

ax[0].plot(temp, L_cool, label='GAMER (reconstructed)', ms=0, ls='solid', lw=0.5, color='black')

ax[0].set_xlabel('Temperature [K]')
ax[0].set_ylabel(r'$\Lambda$ [erg cm$^3$ s$^{-1}$]')
ax[0].set_xscale('log')
ax[0].set_yscale('log')
ax[0].set_xlim(1e4, 2e8)
ax[0].set_ylim(3e-25, 2e-21)
ax[0].legend(ncol=2, fontsize=6)

ax[1].plot(time, temp, 'r-')
ax[1].set_xlabel('a')
ax[1].set_ylabel('Temperature [K]')
ax[1].set_yscale('log')


plt.savefig('CoolingCurve.png', bbox_inches='tight', dpi=300)
