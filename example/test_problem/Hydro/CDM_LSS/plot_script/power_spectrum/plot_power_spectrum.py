#!/usr/bin/env python3.7

import sys
import numpy as np
import math
import yt
import matplotlib.pyplot as plt
from matplotlib.ticker import LogFormatterSciNotation, LogLocator, MultipleLocator, NullFormatter
from scipy.optimize import curve_fit

def linear_fit(x, a, b):
   return a*x + b

box_length       = 30. # Mpc/h
particle_lv      = 7
a_dump           = np.genfromtxt("../../Input__DumpTable", skip_header=1, skip_footer=1)
a_dump           = a_dump[:,1]

data_music       = np.loadtxt("/work1/koarakawaii/music/bin/GENERATE_GADGET_for_CDM/input_powerspec.txt") 
k_music          = data_music[:,0]
PS_music         = data_music[:,-2]

ts0                  = yt.load('../../Data_000000')
base_level_cell_num  = int(ts0.domain_dimensions[1])
max_AMR_level        = int(ts0.parameters["MaxLevel"])

User_Defined_List = np.array(range(20))
colors = plt.cm.jet(np.linspace(0,1,len(User_Defined_List)))

k_gadget_list       = []
PS_gadget_list      = []
scaling_factor_list = []   
PS_lowest_freq_list = []

plt.figure(figsize=(7,3), dpi=320)
plt.subplots_adjust(wspace=0.25)
ax1 = plt.subplot(121)
ax1.loglog(k_music, PS_music, c='k', ls="--", label="MUSIC Input PS", lw=0.5)
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
    scaling_factor_list.append(a_snapshot)
    PS_lowest_freq_list.append(PS_gadget[0])

popt, pcov      = curve_fit(linear_fit, np.log10(np.array(scaling_factor_list)), np.log10(np.array(PS_lowest_freq_list)))
slope,intercept = popt[0], popt[1]

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

ax1.set_title("GAMER-128$^3$: Power Spectrum", fontsize=10.)
ax1.legend(loc="best", prop={'size':2.5})

a_range            = np.linspace(a_dump[0], a_dump[-1], 64)
PS_lowest_freq_fit = 10.**(slope*np.log10(a_range) + intercept)

minor_ticks_x = LogLocator(base = 10.0, subs = np.arange(1.0, 10.0) * 0.1, numticks = 20)
minor_ticks_y = LogLocator(base = 10.0, subs = np.arange(1.0, 10.0) * 0.1, numticks = 20)

plt.savefig("./Power_Spectrum_Output.png", bbox_inches='tight', dpi=320, pad_inches=0.05)


