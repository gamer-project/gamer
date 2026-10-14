#!/usr/bin/env python3

import sys

import load_simulation
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.gridspec import GridSpec

sys.path.append('../shr/')
import SHR

figure_width = 3.375
font_size_regular = 8
font_size_small = 5.6
font_size_large = 8
line_width = 0.6

frames, factor_name = load_simulation.load_sim()
red = plt.cm.Reds(np.linspace(0, 1, 40))
blue = plt.cm.Blues(np.linspace(0, 1, 40))
green = plt.cm.Greens(np.linspace(0, 1, 45))

# Combine the dataframes into one
combined_df = pd.concat(frames)

combined_df = combined_df[combined_df.index %3 == 2]
shr_calculator = SHR.SHR_calculator('planck18')
Nbody = False


fig = plt.figure(figsize=(figure_width, 3.5))
grid = GridSpec(2, 1, height_ratios=[0.9,0.7]) 
ax1 = fig.add_subplot(grid[0, 0])
ax2 = fig.add_subplot(grid[1, 0])

### ax1

### plot c
m = np.logspace(7, 14, 100)

c = shr_calculator.concentration_para_CDM(m, 0)
ax1.plot(m, c, 'k--', label = 'CDM', lw = line_width)
c_1 = shr_calculator.concentration_para_FDM(m, 0, 0.1)
ax1.plot(m, c_1, 'r-', label = 'FDM $m_{22} = 0.1$', lw = line_width+0.1)
c_2 = shr_calculator.concentration_para_FDM(m, 0, 0.2)
ax1.plot(m, c_2, 'b-', label = 'FDM $m_{22} = 0.2$', lw = line_width+0.1)
c_8 = shr_calculator.concentration_para_FDM(m, 0, 0.8)
ax1.plot(m, c_8, 'g-', label = 'FDM $m_{22} = 0.8$', lw = line_width+0.1)

### half mode mass
halfmode_mass_1 = SHR.half_mode_mass(0.1)
ax1.plot([halfmode_mass_1,halfmode_mass_1],[20,15],color = 'r', ls = '--', lw = 0.7,alpha=0.8)
ax1.arrow(halfmode_mass_1,15,0,-1,width = 0,length_includes_head = True,\
    head_width=halfmode_mass_1/6,head_length = 1, color = 'r',alpha=0.8 ,lw = 0)
halfmode_mass_2 = SHR.half_mode_mass(0.2)
ax1.plot([halfmode_mass_2,halfmode_mass_2],[20,15],color = 'b', ls = '--', lw = 0.7,alpha=0.8)
ax1.arrow(halfmode_mass_2,15,0,-1,width = 0,length_includes_head = True,\
    head_width=halfmode_mass_2/6,head_length = 1, color = 'b',alpha=0.8 ,lw = 0)
halfmode_mass_8 = SHR.half_mode_mass(0.8)
ax1.plot([halfmode_mass_8,halfmode_mass_8],[20,15],color = 'g', ls = '--', lw = 0.7,alpha=0.8)
ax1.arrow(halfmode_mass_8,15,0,-1,width = 0,length_includes_head = True,\
    head_width=halfmode_mass_8/6,head_length = 1, color = 'g',alpha=0.8 ,lw = 0)

### plot c_fit
df_range = combined_df[combined_df.index == 68]

for i in range(len(df_range)):
    df = df_range.iloc[i]
    if df['particle_mass']== 8.0e-23:
        ax1.plot(df['halo_mass'], df['c_fit'], marker='^', ls = 'None', mfc = green[20], mec = 'None',\
                        mew = 0.2, markersize = 4)
        if Nbody:
            ax1.plot(df['halo_mass'], df['c_CDM'], marker=(6, 2, 0), ls = 'None', mfc = green[15], mec = green[30],\
                        mew = 1, markersize = 4)
    elif df['particle_mass']== 1.0e-23:
        ax1.plot(df['halo_mass'], df['c_fit'], marker='o', ls = 'None', mfc = red[20], mec = 'None',\
                        mew = 0.2, markersize = 4)
        if Nbody:
            ax1.plot(df['halo_mass'], df['c_CDM'], marker='x', ls = 'None', mfc = red[15], mec = red[30],\
                        mew = 1, markersize = 4)
    else:
        ax1.plot(df['halo_mass'], df['c_fit'], marker='s', ls = 'None', mfc = blue[20], mec = 'None',\
                        mew = 0.2, markersize = 4)
        if Nbody:
            ax1.plot(df['halo_mass'], df['c_CDM'], marker='+', ls = 'None', mfc = blue[15], mec = blue[30],\
                        mew = 1, markersize = 4)

from matplotlib.legend_handler import HandlerTuple
### legend ax1
from matplotlib.lines import Line2D

custom_legend = [
    (Line2D([0], [0], marker='o', mfc = red[20], mec = 'None', markersize=4, linestyle='None'),
    Line2D([0], [0], marker='s', mfc = blue[20], mec = 'None', markersize=4, linestyle='None'),
    Line2D([0], [0], marker='^', mfc = green[20], mec = 'None', markersize=4, linestyle='None')),
    (Line2D([0], [0],marker='x', ls = 'None', mec = red[30], mew = 1, markersize=4,linestyle='None'),
    Line2D([0], [0], marker='+', ls = 'None', mec=blue[30],mew = 1, markersize=4, linestyle='None'),
    Line2D([0], [0], marker=(6, 2, 0), ls = 'None', mec=green[30],mew = 1, markersize=4, linestyle='None')),
    (Line2D([0], [0], color='r', lw = line_width+0.1),
    Line2D([0], [0], color='b', lw = line_width+0.1),
    Line2D([0], [0], color='g', lw = line_width+0.1)),
    Line2D([0], [0], color='k', ls = '--', lw = line_width)  
]

legend_labels = ['FDM simulations', 'N-body simulations','Laroche 2022, ${m}_{22} = 0.1, 0.2, 0.8$','Ishiyama 2021, CDM']

leg1 = ax1.legend(custom_legend, legend_labels, loc='lower right', fontsize=font_size_small,handler_map={tuple: HandlerTuple(ndivide=None)}, handlelength=2.6, labelspacing=0.2)
leg1.get_frame().set_linewidth(0.3)

### annotation z = 0
ax1.annotate('$z = 0$', xy=(0.85, 0.9), xycoords='axes fraction', fontsize=font_size_regular)

### axis
ax1.set_xscale('log')
ax1.set_yscale('log')
ax1.set_xlim(1e10, 2e13)
ax1.set_ylim(1, 24)
# labels
ax1.set_xlabel('${M}_\mathrm{h}$ ($\mathrm{M}_{\odot}$)', fontsize = font_size_regular,labelpad=0.5)
ax1.set_ylabel('concentration $c$', fontsize = font_size_regular,labelpad=0.5)
plt.yticks([1.0,1.1,1.2,1.3,1.4,1.5],[1.0,1.1,1.2,1.3,1.4,1.5])
plt.setp(ax1.get_xticklabels(), fontsize=font_size_regular)
plt.setp(ax1.get_yticklabels(), fontsize=font_size_regular)
# ticks
ax1.tick_params(bottom=True, top=True, left=True, right=True,  which='both',direction ='in')
ax1.tick_params('both', length=2, width=0.5, which='major')
ax1.tick_params('both', length=1, width=0.2, which='minor')
# spines
for axis in ['top','bottom','left','right']:
    ax1.spines[axis].set_linewidth(0.5)

### ax2

load_simulation.plot_m22_redshift( combined_df, 'c_fit', 'v_inner_QP_ave_QP', ax2)
res = load_simulation.fit_regres(combined_df['c_fit'], combined_df['v_inner_QP_ave_QP'], (True, False))
x = np.linspace(1, 20, 100)
# y = res.intercept + res.slope*np.log10(x)
print(res.slope, res.intercept)
y = 0.285*np.log10(x)+1.05  # from the fitting of 42 halos
ax2.plot(x, y, '-.',lw = line_width+0.4, color = 'orange')

from matplotlib.legend_handler import HandlerTuple
### legend ax2
from matplotlib.lines import Line2D

custom_legend = [
    (Line2D([0], [0], marker='o', mfc = red[20], mec = 'None', markersize=4, linestyle='None'),
    Line2D([0], [0], marker='s', mfc = blue[20], mec = 'None', markersize=4, linestyle='None'),
    Line2D([0], [0], marker='^', mfc = green[20], mec = 'None', markersize=4, linestyle='None')),

    Line2D([0], [0], ls = '-.', lw = line_width+0.4,color='orange')
]
legend_labels = ['$m_{22} = 0.1,\ 0.2,\ 0.8$',
                 '$y = 0.285\ \log{x}+1.05$']

leg2 = ax2.legend(custom_legend, legend_labels, loc='lower right', fontsize=font_size_small,handler_map={tuple: HandlerTuple(ndivide=None)}, handlelength=2.9, labelspacing=0.2)
leg2.get_frame().set_linewidth(0.3)


### colorbar
cmap = matplotlib.cm.Greys_r
norm = matplotlib.colors.Normalize(vmin=0, vmax=2)

# Define the range for the custom colormap (half the range)
custom_cmap_range = np.linspace(0, 0.85, 256)

# Create a custom colormap representing only half of the original colormap
custom_cmap = matplotlib.colors.LinearSegmentedColormap.from_list("custom_Greys_r", cmap(custom_cmap_range))

from mpl_toolkits.axes_grid1 import make_axes_locatable

# Create a divider for the lower subplot
divider = make_axes_locatable(ax2)

# Define size and position for the colorbar
cax = divider.append_axes('right', size='10%', pad=0.1, aspect=10)

cbar = fig.colorbar(matplotlib.cm.ScalarMappable(norm=norm, cmap=custom_cmap),
                orientation='vertical', ticks = [0,1,2], alpha = 0.9, cax = cax)
cbar.ax.invert_yaxis()
cbar.set_label('redshift $z$', fontsize = font_size_regular)
cbar.ax.tick_params(labelsize=font_size_regular, width=0.5, length = 2)
cbar.outline.set_linewidth(0.5)

### axis
ax2.set_xscale('log')
ax2.set_ylim(1,1.5)
ax2.set_xlim(1.5,18)
# labels
ax2.set_xlabel('concentration $c$', fontsize = font_size_regular,labelpad=0.5)
ax2.set_ylabel(r'${w_\mathrm{h, in}}\ /\ {\langle w \rangle}_\mathrm{h}$', fontsize = font_size_regular, labelpad=0.5)
plt.setp(ax2.get_xticklabels(), fontsize=font_size_regular)
plt.setp(ax2.get_yticklabels(), fontsize=font_size_regular)
# ticks
ax2.tick_params(bottom=True, top=True, left=True, right=True,  which='both',direction ='in')
ax2.tick_params('both', length=2, width=0.5, which='major')
ax2.tick_params('both', length=1, width=0.2, which='minor')
# spines
for axis in ['top','bottom','left','right']:
    ax2.spines[axis].set_linewidth(0.5)

ax1.annotate('(a)', xy=(0.03, 0.92), xycoords='axes fraction', fontsize=font_size_regular)
ax2.annotate('(b)', xy=(0.03, 0.92), xycoords='axes fraction', fontsize=font_size_regular)


plt.tight_layout(pad=0.05, w_pad=0, h_pad=0)

plt.savefig('fig_c_nonisothermality.pdf', dpi=300)