#!/usr/bin/env python3

import sys

import load_simulation
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.append('../shr/')
import SHR
from matplotlib.gridspec import GridSpec

figure_width = 3.375
font_size_regular = 8
font_size_small = 5.6
font_size_large = 8
line_width = 0.6

### load data
frames, factor_name = load_simulation.load_sim()

# Combine the dataframes into one
combined_df = pd.concat(frames)
combined_df = combined_df[combined_df.index %3 == 2]

shr_calculator = SHR.SHR_calculator('planck18')


fig = plt.figure(figsize=(figure_width, 3.5))
grid = GridSpec(2, 1, height_ratios=[0.55,1]) 
ax1 = fig.add_subplot(grid[1, 0])
ax2 = fig.add_subplot(grid[0, 0])

red = plt.cm.Reds(np.linspace(0, 1, 40))
blue = plt.cm.Blues(np.linspace(0, 1, 40))
green = plt.cm.Greens(np.linspace(0, 1, 45))

# combined_df = combined_df[combined_df.index %4==0]
combined_df['ratio_revised'] = combined_df['core_mass']*0
combined_df['ratio_schive'] = combined_df['core_mass']*0
ratio_revised = []
ratio_shcive = []
for i in range(len(combined_df)):
    ratio_revised.append(combined_df['core_mass'].iloc[i]/shr_calculator.revised_theo_c_FDM_Ms(combined_df['time_z'].iloc[i], combined_df['halo_mass'].iloc[i], combined_df['m22'].iloc[i]))
    ratio_shcive.append(combined_df['core_mass'].iloc[i]/shr_calculator.theo_TH_Ms(combined_df['time_z'].iloc[i], combined_df['halo_mass'].iloc[i], combined_df['m22'].iloc[i]))

ratio_revised = np.array(ratio_revised)
ratio_shcive = np.array(ratio_shcive)
combined_df['ratio_revised'] = ratio_revised
combined_df['ratio_schive'] = ratio_shcive

### plot data
load_simulation.plot_m22_redshift(combined_df, 'halo_mass', 'ratio_revised', ax1)

def plot_m22_schive(combined_df, x_axis, fac, ax):

    m22_8 = combined_df[combined_df['particle_mass'] == 8.0e-23]
    m22_2 = combined_df[combined_df['particle_mass'] == 2.0e-23]
    m22_1 = combined_df[combined_df['particle_mass'] == 1.0e-23]


    for j in range(len(m22_8)):
        # ax.plot(m22_8[x_axis].iloc[j], m22_8[fac].iloc[j], marker='^', ls = 'None', mfc = 'None', mec = green[30], mew = 0.4, markersize = 4)
        ax.plot(m22_8[x_axis].iloc[j], m22_8[fac].iloc[j], marker='^', ls = 'None', mfc = 'None', mec = green[m22_8.index[j]-25], mew = 0.4, markersize = 4)
    for j in range(len(m22_2)):
        # ax.plot(m22_2[x_axis].iloc[j], m22_2[fac].iloc[j], marker='s', ls = 'None', mfc = 'None', mec = blue[30], mew = 0.4, markersize = 4)
        ax.plot(m22_2[x_axis].iloc[j], m22_2[fac].iloc[j], marker='s', ls = 'None', mfc = 'None', mec = blue[m22_2.index[j]-30], mew = 0.4, markersize = 4)
    for j in range(len(m22_1)):
        # ax.plot(m22_1[x_axis].iloc[j], m22_1[fac].iloc[j], marker='o', ls = 'None', mfc = 'None', mec = red[30], mew = 0.4, markersize = 4)
        ax.plot(m22_1[x_axis].iloc[j], m22_1[fac].iloc[j], marker='o', ls = 'None', mfc = 'None', mec = red[m22_1.index[j]-30], mew = 0.4, markersize = 4)

plot_m22_schive(combined_df, 'halo_mass', 'ratio_schive', ax1)


interval, x_ticks = load_simulation.bin(8, 'halo_mass', 'log', combined_df)
combined_df['bins'] = pd.cut(combined_df['halo_mass'], bins=interval)
print(combined_df.groupby('bins')['ratio_revised'].size())
median = combined_df.groupby('bins')['ratio_revised'].quantile(q = 0.5).values
up = combined_df.groupby('bins')['ratio_revised'].quantile(q = 0.84).values
low = combined_df.groupby('bins')['ratio_revised'].quantile(q = 0.16).values

ax1.plot(x_ticks, median,'--', color = 'black', lw = line_width)
ax1.fill_between(x_ticks, low, up, color = 'gold', alpha = 0.4)


### axis
ax1.set_xscale('log')
ax1.set_yscale('log')
ax1.set_xlim(6e9, 4e12)
ax1.set_ylim( 0.5, 2)

# ax1is labels
ax1.set_xlabel( r'${M}_{\mathrm{h}}\:(\mathrm{M}_{\odot})$', fontsize = font_size_regular, labelpad=1)
ax1.set_ylabel( r'${{M}_{\mathrm{s, sim}}}\ /\ {{M}_{\mathrm{s, theory}}}$', fontsize = font_size_regular , rotation=90, labelpad=0)
# ax1.setp(ax1.get_xticklabels(), fontsize=font_size_regular)
# ax1.setp(ax1.get_yticklabels(), fontsize=font_size_regular)
# ticks
ax1.tick_params(bottom=True, top=True, left=True, right=True,  which='both',direction ='in',labelsize=font_size_regular)
ax1.tick_params('both', length=2, width=0.5, which='major')
ax1.tick_params('both', length=1, width=0.2, which='minor')
ax1.set_yticks([0.5,1,2], [0.5,1,2])
ax1.minorticks_on()
from matplotlib.ticker import NullFormatter

ax1.yaxis.set_minor_formatter(NullFormatter())
# spines
for axis in ['top','bottom','left','right']:
    ax1.spines[axis].set_linewidth(0.5)



### colorbar
cmap = matplotlib.cm.Greys_r
norm = matplotlib.colors.Normalize(vmin=0, vmax=2)
# Define the range for the custom colormap (half the range)
custom_cmap_range = np.linspace(0, 0.85, 256)
# Create a custom colormap representing only half of the original colormap
custom_cmap = matplotlib.colors.LinearSegmentedColormap.from_list("custom_Greys_r", cmap(custom_cmap_range))

from mpl_toolkits.axes_grid1 import make_axes_locatable

# Create a divider for the lower subplot
divider = make_axes_locatable(ax1)
cax = divider.append_axes('right', size='10%', pad=0.1, aspect=10)

cbar = fig.colorbar(matplotlib.cm.ScalarMappable(norm=norm, cmap=custom_cmap),
                orientation='vertical', ticks = [0,1,2], alpha = 0.9, cax = cax)
cbar.ax.invert_yaxis()
cbar.set_label('redshift $z$', fontsize = font_size_regular)
cbar.ax.tick_params(labelsize=font_size_regular, width=0.5, length = 2)
cbar.outline.set_linewidth(0.5)

import matplotlib.legend as mlegend
from matplotlib.legend_handler import HandlerTuple
### legend
from matplotlib.lines import Line2D

custom_legend_1 = [
    (Line2D([0], [0], marker='o', color=red[20], mew = 0, markersize=4, linestyle='None'),
    Line2D([0], [0], marker='s', color=blue[20], mew = 0, markersize=4, linestyle='None'),
    Line2D([0], [0], marker='^', color=green[20], mew = 0, markersize=4, linestyle='None')),
    (Line2D([0], [0], marker='o', color='w', mec=red[20],mew = 0.4, markersize=4, linestyle='None'),
    Line2D([0], [0], marker='s', color='w', mec=blue[20],mew = 0.4, markersize=4, linestyle='None'),
    Line2D([0], [0], marker='^', color='w', mec=green[20],mew = 0.4, markersize=4, linestyle='None'))
]
custom_legend_2 = [
    (Line2D([0], [0], color=red[20], lw = 1),
    Line2D([0], [0], color=blue[20], lw = 1),
    Line2D([0], [0], color=green[20], lw = 1)),
    (Line2D([0], [0], color='white', lw = 1),
    Line2D([0], [0], color='white', lw = 1),
    Line2D([0], [0], color='white', lw = 1))
]
legend_labels_1 = [ 'This work','Schive 2014', '           ' ]

legend_labels_2 = [ '$m_{22}=0.1, 0.2, 0.8$', '      ']

leg_1 = ax1.legend(custom_legend_1, legend_labels_1, loc='lower right', fontsize=font_size_small,handler_map={tuple: HandlerTuple(ndivide=None)}, handlelength=2.2,labelspacing=0.2, ncol=2, columnspacing=1.5)
leg_2 = mlegend.Legend(ax1,custom_legend_2, legend_labels_2, loc='upper center', fontsize=font_size_small,handler_map={tuple: HandlerTuple(ndivide=None)}, handlelength=2.3,labelspacing=0.2, ncol=2, columnspacing=1.5)
leg_1._legend_box._children.append(leg_2._legend_box._children[1])
leg_1.get_frame().set_linewidth(0.3)


halo_mass_x = np.logspace(9.7,12.7,50)

m23_2_revised = shr_calculator.revised_theo_c_FDM_Ms(0, halo_mass_x, 0.2)
m23_2_TH = shr_calculator.theo_TH_Ms(0, halo_mass_x, 0.2)
ratio = m23_2_revised/m23_2_TH
ax2.plot(halo_mass_x, m23_2_revised, lw = line_width,color = blue[30])
ax2.plot(halo_mass_x, m23_2_TH,'--',lw = line_width,color = blue[30])

m23_1_revised = shr_calculator.revised_theo_c_FDM_Ms(0, halo_mass_x, 0.1)
m23_1_TH = shr_calculator.theo_TH_Ms(0, halo_mass_x, 0.1)
ratio = m23_1_revised/m23_1_TH
ax2.plot(halo_mass_x, m23_1_revised,lw = line_width, color = red[30])
ax2.plot(halo_mass_x, m23_1_TH, '--',lw = line_width,color = red[30])

m23_8_revised = shr_calculator.revised_theo_c_FDM_Ms(0, halo_mass_x, 0.8)
m23_8_TH = shr_calculator.theo_TH_Ms(0, halo_mass_x, 0.8)
ratio = m23_8_revised/m23_8_TH
ax2.plot(halo_mass_x, m23_8_revised,lw = line_width, color = green[30])
ax2.plot(halo_mass_x, m23_8_TH,'--',lw = line_width, color = green[30])

load_simulation.plot_single(frames, 'halo_mass', 'core_mass',ax2, group_m22 = True, time_range=0.98)

### axis
ax2.set_xscale('log')
ax2.set_yscale('log')
ax2.set_xlim(1e10, 4e12)
ax2.set_ylim( 3e7, 8e9)
# labels
ax2.set_xlabel( r'${M}_{\mathrm{h}}\:(\mathrm{M}_{\odot})$', fontsize = font_size_regular, labelpad=1)
ax2.set_ylabel( r'${M}_{\mathrm{s}}\:(\mathrm{M}_{\odot})$', fontsize = font_size_regular , rotation=90, labelpad=0)


# ticks
ax2.tick_params(bottom=True, top=True, left=True, right=True,  which='both',direction ='in',labelsize=font_size_regular)
ax2.tick_params('both', length=2, width=0.5, which='major')
ax2.tick_params('both', length=1, width=0.2, which='minor')
# plt.yticks([0.5,1,2], [0.5,1,2])
plt.minorticks_on()
from matplotlib.ticker import NullFormatter

ax2.yaxis.set_minor_formatter(NullFormatter())
# spines
for axis in ['top','bottom','left','right']:
    ax2.spines[axis].set_linewidth(0.5)

custom_legend1 = [
    (Line2D([0], [0], marker='o', color=red[20], mew = 0, markersize=4, linestyle='None'),
    Line2D([0], [0], marker='s', color=blue[20], mew = 0, markersize=4, linestyle='None'),
    Line2D([0], [0], marker='^', color=green[20], mew = 0, markersize=4, linestyle='None')),
    (Line2D([0], [0], marker='o', color='w', mec=red[20],mew = 0, markersize=0, linestyle='None'))
]
custom_legend2 = [
    Line2D([0], [0],ls = '-', color='black', lw = line_width),
    Line2D([0], [0],ls= '--', color='black', lw = line_width)
    ]

legend_labels1 = ['$M_{s, sim} , m_{22}=0.1, 0.2, 0.8$','']
legend_labels2 = [ 'This work', 'Schive 2014']

leg = ax2.legend(custom_legend1, legend_labels1, loc='lower right', fontsize=font_size_small,handler_map={tuple: HandlerTuple(ndivide=None)}, handlelength=2.2, labelspacing=0.1,ncol=2, columnspacing=-0.4)
leg2 = mlegend.Legend(ax2,custom_legend2, legend_labels2, loc='upper center', fontsize=font_size_small,handler_map={tuple: HandlerTuple(ndivide=None)}, handlelength=2.2,labelspacing=0.2, ncol=2, columnspacing=1)

leg._legend_box._children.append(leg2._legend_box._children[1])
leg.get_frame().set_linewidth(0.3)

ax2.annotate('z = 0', xy=(0.13, 0.8), xycoords='axes fraction', fontsize=font_size_regular)

ax1.annotate('(b)', xy=(0.02, 0.92), xycoords='axes fraction', fontsize=font_size_regular)
ax2.annotate('(a)', xy=(0.02, 0.9), xycoords='axes fraction', fontsize=font_size_regular)

plt.tight_layout(pad=0.03, w_pad=0, h_pad=0)
plt.savefig('fig_SHR.pdf', dpi = 250)
plt.close()

