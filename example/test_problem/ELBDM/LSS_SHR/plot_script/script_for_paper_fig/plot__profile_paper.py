#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import yt

# set plot parameters
figure_width = 3.375
font_size_regular = 8
font_size_small = 6
font_size_large = 8
line_width = 0.6

path = '../'
compare_path = ''
Nbody_path = ''
halo = 1
idx = 68
softening_length = 0.3125



def plot_profile(path, name, core_is_true, NFW_is_true, order):

    def soliton(x, core_radius, particle_mass):   
        return ((1.9*(particle_mass/10**-23)**-2*(core_radius**-4))/((1 + 9.1*10**-2*(x/core_radius)**2)**8))*10**9

    # read data
    df = pd.read_csv( path+'/prof_dens/Data_%06d_%d_profile_data'%(idx, halo) , sep = '\t' , header = 0 )
    df_halo_parameter = pd.read_csv( path+'/Halo_Parameter_%d'%1 , sep = r'\s+' , header = 0 , index_col='#')

    current_time_z = df_halo_parameter['time'][idx]
    current_time_a = 1/(current_time_z+1)
    radius = df['radius(kpccm)'][:]*current_time_a
    density = df['density(Msun/kpccm**3)'][:]/current_time_a**3
    dens = density
    halo_radius = df_halo_parameter['halo_radius'][idx]/current_time_a
    
    # plot

    if order == 1: # high-res
        ax.plot( radius, dens, '-', color = 'saddlebrown',lw = 2, label=name)
    elif order == 2: # low-res
        ax.plot( radius, dens,ls = '--', color = 'orchid', lw = 1.5, label=name)
    elif order == 3: # N-body
        ax.plot( radius[radius>softening_length], dens[radius>softening_length],ls =(0,(1,1)), color = 'forestgreen', lw = 2.5, label=name)

    # plt.plot( [halo_radius, halo_radius], [1e-1, 1e9], '--', label= name+' halo radius')
    
    # if FDM
    if (core_is_true):
        particle_mass = df_halo_parameter['mass'][idx]

        core_radius = df_halo_parameter['core_radius_1'][idx]
        x = np.logspace(-1, 3, num=50)
        if order == 1:
            ax.plot( x, soliton(x, core_radius, particle_mass), ':',lw = 1, color = 'peru', label='FDM high-res soliton fit')
        elif order == 2:
            ax.plot( x, soliton(x, core_radius, particle_mass), ':',lw = 1, color = 'pink', label='FDM low-res soliton fit')

    # if NFW
    if (NFW_is_true):
        rho0 = df_halo_parameter['rho0'][idx]
        Rs = df_halo_parameter['Rs'][idx]
        x = np.logspace(-1, 3, num=50)
        ax.plot(x, NFW_dens(x, [rho0, Rs]), ls = '-.',color = 'orange',lw =1.3, label='NFW')

def NFW_dens(r, dens_parameter):
    rho0 = dens_parameter[0]
    Rs = dens_parameter[1]
    return rho0/(r/Rs*(1+r/Rs)**2)


fig, ax = plt.subplots(figsize=(figure_width, 2.2))

### plot density profile
# plot_profile(comparepath, 'FDM low-res', True, False, 2)
# plot_profile(Nbody_path, 'N-body', False, False, 3)
plot_profile(path, 'FDM high-res', True, True, 1)

### set plot
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlim(2e-1,6e2)
ax.set_ylim(5e2,1e10)

# set labels
ax.set_ylabel('density $(\mathrm{M}_{\odot}\mathrm{kpc}^{-3})$', fontsize = font_size_regular,labelpad= 1)
ax.set_xlabel('radius $r$ (kpc)', fontsize = font_size_regular,labelpad= 1)
plt.setp(ax.get_xticklabels(), fontsize = font_size_regular)
plt.setp(ax.get_yticklabels(), fontsize = font_size_regular)

# set ticks
ax.tick_params(bottom=True, top=True, left=True, right=True,  which='both',direction ='in')
ax.tick_params('both', length=2, width=0.5, which='major')
ax.tick_params('both', length=1, width=0.2, which='minor')

# Set y-axis ticks and labels with scientific notation
from matplotlib.ticker import LogLocator

y_major = LogLocator(base=10,numticks=9)
y_minor = LogLocator(base=10, subs=[2.0, 5.0, 8.0], numticks=10)
ax.yaxis.set_major_locator(y_major)
ax.yaxis.set_minor_locator(y_minor)


# spines
for axis in ['top','bottom','left','right']:
    ax.spines[axis].set_linewidth(0.5)

### set legend
from matplotlib.legend_handler import HandlerTuple
from matplotlib.lines import Line2D

custom_legend = [
     Line2D([0], [0], color='saddlebrown', lw = 2)
    # ,Line2D([0], [0], ls = '--', color = 'orchid', lw = 1.5)
    # ,Line2D([0], [0], ls =(0,(1,1)), color = 'forestgreen', lw = 2.5)
    ,Line2D([0], [0], color='gray', linestyle=':', lw = 1)
    ,Line2D([0], [0], color='orange', linestyle='-.', lw =1.3)
    ]

legend_labels = [
    'FDM (high-res)'
    # ,'FDM (low-res)'
    # ,'N-body'
    ,'Soliton fit'
    ,'NFW fit'
    ]

leg = ax.legend(custom_legend, legend_labels, loc='upper right', fontsize = font_size_small,
                handler_map={tuple: HandlerTuple(ndivide=None)}, handlelength=2.2)
leg.get_frame().set_linewidth(0.3)

plt.tight_layout(pad=0.05, w_pad=0, h_pad=0)

FileOut = 'fig_profile_density'+'.pdf'
plt.savefig( FileOut, dpi = 150)
plt.close()
