#!/usr/bin/env python3

import sys

import load_simulation
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
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

### load data
frames, factor_name = load_simulation.load_sim()
red = plt.cm.Reds(np.linspace(0, 1, 40))
blue = plt.cm.Blues(np.linspace(0, 1, 40))
green = plt.cm.Greens(np.linspace(0, 1, 45))

# Combine the dataframes into one
combined_df = pd.concat(frames)

combined_df = combined_df[combined_df.index %3 == 2]

fig = plt.figure(figsize=(figure_width, 2.5))
grid = GridSpec(5, 5, height_ratios=[0.05,1,0.3,1,0.3], width_ratios=[1,0.35, 1,0.25,0.1], wspace=0, hspace=0)

# # Create subplots using GridSpec indexing
ax1 = fig.add_subplot(grid[1, 0])
ax2 = fig.add_subplot(grid[1, 2])
ax3 = fig.add_subplot(grid[3, 0])
ax4 = fig.add_subplot(grid[3, 2])
ax5 = fig.add_subplot(grid[:, 4])  # Span all rows in third column


x_axis = 'halo_mass'

### ax1 virial ratio

fac = 'Ek2_Ep'
load_simulation.plot_m22_redshift(combined_df, x_axis, fac, ax1)
ax1.plot([7e9, 5e12], [1, 1], 'k--', lw=line_width)
ax1.plot([7e9, 5e12], [1.35, 1.35], ':',color = 'orange', lw=line_width)

ax1.text(0.98,0.15,'virial condition',ha='right', va='top', fontsize=font_size_small+0.5,transform=ax1.transAxes)
ax1.text(0.98,0.85, 'virialized halos', color = 'darkorange',ha='right', va='top', fontsize=font_size_small+0.5,transform=ax1.transAxes)

ax1.set_xscale('log')
ax1.set_xlim(7e9, 5e12)
ax1.set_ylim(0.65, 1.65)

ax1.set_xlabel(r'${M}_{\mathrm{h}}\:(\mathrm{M}_{\odot})$', fontsize = font_size_regular, labelpad = 0)
ax1.set_ylabel(r'$2E_{\mathrm{k}}\ /\ |E_{\mathrm{p}}|$', fontsize = font_size_regular, labelpad = 0)
plt.setp(ax1.get_xticklabels(), fontsize=font_size_small)
plt.setp(ax1.get_yticklabels(), fontsize=font_size_small)

ax1.tick_params(bottom=True, top=True, left=True, right=True,  which='both',direction ='in')
ax1.tick_params('both', length=2, width=0.5, which='major')
ax1.tick_params('both', length=1, width=0.2, which='minor')
ax1.yaxis.set_minor_locator(ticker.AutoMinorLocator())

for axis in ['top','bottom','left','right']:
    ax1.spines[axis].set_linewidth(0.5)

### ax2 equipartition

combined_df['v_ave_Bulk'] = (combined_df['v_ave_Both']**2 - combined_df['v_ave_QP']**2)**0.5
combined_df['v_ave_QP_Bulk'] = combined_df['v_ave_QP']/combined_df['v_ave_Bulk']

fac = 'v_ave_QP_Bulk'
load_simulation.plot_m22_redshift(combined_df, x_axis, fac, ax2)
ax2.plot([7e9, 5e12], [1, 1], 'k--', lw=line_width)

ax2.text(0.98,0.15,'equipartition',ha='right', va='top', fontsize=font_size_small+0.5,transform=ax2.transAxes)

ax2.set_xscale('log')
ax2.set_xlim(7e9, 5e12)
ax2.set_ylim(0.5, 1.5)

ax2.set_xlabel(r'${M}_{\mathrm{h}}\:(\mathrm{M}_{\odot})$', fontsize = font_size_regular, labelpad = 0)
ax2.set_ylabel(r'$\langle w \rangle_{\mathrm{h}}\ /\ \langle v \rangle_{\mathrm{h}}$', fontsize = font_size_regular, labelpad = 1)
plt.setp(ax2.get_xticklabels(), fontsize=font_size_small)
plt.setp(ax2.get_yticklabels(), fontsize=font_size_small)

ax2.tick_params(bottom=True, top=True, left=True, right=True,  which='both',direction ='in')
ax2.tick_params('both', length=2, width=0.5, which='major')
ax2.tick_params('both', length=1, width=0.2, which='minor')
ax2.yaxis.set_minor_locator(ticker.AutoMinorLocator())

for axis in ['top','bottom','left','right']:
    ax2.spines[axis].set_linewidth(0.5)

### ax3 Thermal equlibrium

fac = 'v_soliton_QP_inner_QP'
load_simulation.plot_m22_redshift(combined_df, x_axis, fac, ax3)
print(np.mean(combined_df[fac]))
# 0.89 is from the average of 42 halos
ax3.plot([7e9, 5e12], [0.89, 0.89], 'k--', lw=line_width)

ax3.text(0.98,0.15,'thermal equilibrium',ha='right', va='top', fontsize=font_size_small+0.5,transform=ax3.transAxes)

ax3.set_xscale('log')
ax3.set_xlim(7e9, 5e12)
ax3.set_ylim(0.5, 1.5)

ax3.set_xlabel(r'${M}_{\mathrm{h}}\:(\mathrm{M}_{\odot})$', fontsize = font_size_regular, labelpad = 0)
ax3.set_ylabel(r'$\langle w \rangle_{\mathrm{s}}\ /\ w_{\mathrm{h}, \mathrm{in}}$', fontsize = font_size_regular, labelpad = 0)
plt.setp(ax3.get_xticklabels(), fontsize=font_size_small)
plt.setp(ax3.get_yticklabels(), fontsize=font_size_small)

ax3.tick_params(bottom=True, top=True, left=True, right=True,  which='both',direction ='in')
ax3.tick_params('both', length=2, width=0.5, which='major')
ax3.tick_params('both', length=1, width=0.2, which='minor')
ax3.yaxis.set_minor_locator(ticker.AutoMinorLocator())

for axis in ['top','bottom','left','right']:
    ax3.spines[axis].set_linewidth(0.5)

### ax4 Soliton integrity

core_mass_infer = []
for i in range(len(combined_df)):
    core_mass_infer.append(SHR.soliton_m_div_v(combined_df['m22'].iloc[i])* combined_df['v_soliton_QP'].iloc[i]/SHR.kpc2km)
combined_df['core_mass_infer'] = core_mass_infer
combined_df['core_mass_core_mass_infer'] = combined_df['core_mass']/combined_df['core_mass_infer']
fac = 'core_mass_core_mass_infer'
load_simulation.plot_m22_redshift(combined_df, x_axis, fac, ax4)
ax4.plot([7e9, 5e12], [1, 1], 'k--', lw=line_width)

ax4.text(0.98,0.15,'soliton fidelity',ha='right', va='top', fontsize=font_size_small+0.5,transform=ax4.transAxes)

ax4.set_xscale('log')
ax4.set_xlim(7e9, 5e12)
ax4.set_ylim(0.5, 1.5)

ax4.set_xlabel(r'${M}_{\mathrm{h}}\:(\mathrm{M}_{\odot})$', fontsize = font_size_regular, labelpad = 0)
ax4.set_ylabel(r'$M_{\mathrm{s, sim}}\ /\ M_{\mathrm{s}, \langle w \rangle_{\mathrm{s}}}$', fontsize = font_size_regular, labelpad = 1)
plt.setp(ax4.get_xticklabels(), fontsize=font_size_small)
plt.setp(ax4.get_yticklabels(), fontsize=font_size_small)

ax4.tick_params(bottom=True, top=True, left=True, right=True,  which='both',direction ='in')
ax4.tick_params('both', length=2, width=0.5, which='major')
ax4.tick_params('both', length=1, width=0.2, which='minor')
ax4.yaxis.set_minor_locator(ticker.AutoMinorLocator())

for axis in ['top','bottom','left','right']:
    ax4.spines[axis].set_linewidth(0.5)

## colorbar

cmap = matplotlib.cm.Greys_r
norm = matplotlib.colors.Normalize(vmin=0, vmax=2)

# Define the range for the custom colormap (half the range)
custom_cmap_range = np.linspace(0, 0.85, 256)

# Create a custom colormap representing only half of the original colormap
custom_cmap = matplotlib.colors.LinearSegmentedColormap.from_list("custom_Greys_r", cmap(custom_cmap_range))

# Plot the custom colormap
ax5.axis('off') 
cbar = fig.colorbar(matplotlib.cm.ScalarMappable(norm=norm, cmap=custom_cmap),ax = ax5,\
                orientation='vertical', ticks = [0,1,2], alpha = 0.9, \
                fraction=1.1, anchor = (0,0.9), aspect = 15, pad = 0)
cbar.ax.invert_yaxis()
cbar.set_label('redshift $z$', fontsize = font_size_regular, labelpad = 1)
cbar.ax.tick_params(labelsize=font_size_small, width=0.3, length = 1.5)
cbar.outline.set_linewidth(0.5)


### legend
from matplotlib.legend_handler import HandlerTuple
from matplotlib.lines import Line2D

custom_legend = [
    Line2D([0], [0], marker='o', mfc = red[20], mec = 'None', markersize=4, linestyle='None'),
    Line2D([0], [0], marker='s', mfc = blue[20], mec = 'None', markersize=4, linestyle='None'),
    Line2D([0], [0], marker='^', mfc = green[20], mec = 'None', markersize=4, linestyle='None'),
    Line2D([0], [0],ls = '--', color = 'k', lw = 0.6)]
legend_labels = ['$m_{22} = 0.1$', '$m_{22} = 0.2$','$m_{22} = 0.8$', 'Fiducial\nvalues']

leg = fig.legend(custom_legend, legend_labels, fontsize=font_size_small-0.5,handler_map={tuple: HandlerTuple(ndivide=None)},\
     handlelength=1.1, bbox_to_anchor=(0.96, 0.1), loc="lower right", borderaxespad=-1.5)

leg.get_frame().set_linewidth(0.3)

ax1.annotate('(a)', xy=(0.03, 0.87), xycoords='axes fraction', fontsize=font_size_small+1)
ax2.annotate('(b)', xy=(0.03, 0.87), xycoords='axes fraction', fontsize=font_size_small+1)
ax3.annotate('(c)', xy=(0.03, 0.87), xycoords='axes fraction', fontsize=font_size_small+1)
ax4.annotate('(d)', xy=(0.03, 0.87), xycoords='axes fraction', fontsize=font_size_small+1)

plt.tight_layout(pad=0, w_pad=0, h_pad=0)

plt.savefig('fig_four_ratio.pdf')



