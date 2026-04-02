#!/usr/bin/env python3.7

import sys
import numpy as np
import math
import yt
import matplotlib.pyplot as plt
from matplotlib.ticker import LogFormatterSciNotation, LogLocator, MultipleLocator, NullFormatter

particle_lv      = 7
a_dump           = np.genfromtxt("../../Record__Dump", skip_header=1)
a_dump           = a_dump[:,1]

data_music       = np.loadtxt("../../input_powerspec.txt")
k_music          = data_music[:,0]
PS_music         = data_music[:,-2]

data_camb        = np.loadtxt("powerspec_z0_Planck18.dat")
k_camb           = data_camb[:,0]
PS_camb          = data_camb[:,1]

data_hires       = np.loadtxt("powerspec_highres_Data_000072.dat")
k_hires          = data_hires[:,0]
PS_hires         = data_hires[:,1] / (8*math.pi**3)

ts0                  = yt.load('../../Data_000000')
base_level_cell_num  = int(ts0.domain_dimensions[1])
max_AMR_level        = int(ts0.parameters["MaxLevel"])
box_length           = float(ts0.domain_width[1].in_units("Mpc/h"))  # [Mpc/h] = [code_length]
NPar_base            = int(round(np.power(int(ts0["Par_NPar"]), 0.3333333)))

User_Defined_List = np.arange(0, 73, 4)
colors = plt.cm.jet(np.linspace(0,1,len(User_Defined_List)))

k_gadget_list       = []
PS_gadget_list      = []

plt.figure(figsize=(7,3), dpi=320)
plt.subplots_adjust(wspace=0.25)
ax1 = plt.subplot(121)
ax1.loglog(k_music, PS_music, c='k', ls="--", label="MUSIC Input PS", lw=0.5)
ax1.loglog(k_camb, PS_camb, c='r', ls="--", label="CAMB PS (z=0; Planck18)", lw=0.5)
ax1.loglog(k_hires, PS_hires, c='b', ls="--", label="$z$ = 0.0 (NMesh=512)", lw=0.5)
k_min = 2.*np.pi/box_length
k_max = k_min*2**(particle_lv)/2.
ax1.loglog([k_min,k_min],[10**-11,10**2.], c="orange", lw=0.5, ls=":")
ax1.loglog([k_max,k_max],[10**-11,10**2.], c="orange", lw=0.5, ls=":")
ax1.text(0.9e-1, 1.17e-8, "[NX0_TOT:%i$^3$][AMR_MAX:%i]"%(base_level_cell_num,max_AMR_level), fontsize=5.5)

for counter, snapshot_idx in enumerate(User_Defined_List):
    a_snapshot = a_dump[snapshot_idx]
    data_temp  = np.loadtxt("../../PowerSpec_%06d"%snapshot_idx, skiprows=1, dtype=float)

    k_gadget   = data_temp[:,0]
    PS_gadget  = data_temp[:,1]/(8*math.pi**3) # PowerSpectrum_MUSIC * 8.0*\pi^3 = PowerSpectrum_GAMER

    ax1.loglog(k_gadget, PS_gadget, ls='-', c=colors[counter], label="$z$ = %.3f"%(1./a_snapshot-1.), lw=0.5)

minor_ticks_x = LogLocator(base = 10.0, subs = np.arange(1.0, 10.0) * 0.1, numticks = 20)
minor_ticks_y = LogLocator(base = 10.0, subs = np.arange(1.0, 10.0) * 0.1, numticks = 20)

ax1.set_xlabel("$k$ (Mpc/h)$^{-1}$", fontsize = 8.)
ax1.set_ylabel("$P(k)$ (Mpc/h)$^3$", fontsize = 8.)
ax1.set_xlim([10**-1.1,10**3.1])
ax1.set_xticks(10.**np.arange(-1,3.1,1.))
ax1.set_xticklabels(10.**np.arange(-1,3.1,1.), fontsize=7.)
ax1.xaxis.set_minor_locator(minor_ticks_x)
ax1.xaxis.set_major_formatter(LogFormatterSciNotation(base=10.0))
ax1.xaxis.set_minor_formatter(NullFormatter())
ax1.set_ylim([10**-8.1,10**0.6])
ax1.set_yticks(10.**np.arange(-8,1.1,2.))
ax1.set_yticklabels(10.**np.arange(-8,1.1,2.), fontsize=7.)
ax1.yaxis.set_minor_locator(minor_ticks_y)
ax1.yaxis.set_major_formatter(LogFormatterSciNotation(base=10.0))
ax1.yaxis.set_minor_formatter(NullFormatter())
ax1.set_title("GAMER-%i$^3$: Power Spectrum"%NPar_base, fontsize=10.)
ax1.legend(loc="best", prop={'size':2.5})

plt.savefig("./Power_Spectrum_Output.png", bbox_inches='tight', dpi=320, pad_inches=0.05)
