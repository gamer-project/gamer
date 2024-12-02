#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np
import yt
from mpl_toolkits.axes_grid1 import AxesGrid



###################################################################################################
# Genereal Setting
###################################################################################################
dpi        = 600
FONT_SIZE  = 8.0
LINE_WIDTH = 1.0

plt.rcParams['font.size']         = FONT_SIZE
plt.rcParams['figure.titlesize']  = FONT_SIZE
plt.rcParams['figure.dpi']        = dpi

plt.rcParams['axes.titlesize']    = FONT_SIZE
plt.rcParams['axes.labelsize']    = FONT_SIZE
plt.rcParams['axes.labelpad']     = FONT_SIZE/8.0*0.05
plt.rcParams['axes.linewidth']    = LINE_WIDTH

plt.rcParams['legend.fontsize']   = FONT_SIZE
plt.rcParams['lines.linewidth']   = LINE_WIDTH

plt.rcParams['xtick.major.size']  = 4.0*LINE_WIDTH
plt.rcParams['xtick.major.width'] = LINE_WIDTH
plt.rcParams['xtick.minor.size']  = 2.0*LINE_WIDTH
plt.rcParams['xtick.minor.width'] = 0.5*LINE_WIDTH
plt.rcParams['xtick.labelsize']   = FONT_SIZE

plt.rcParams['ytick.major.size']  = 4.0*LINE_WIDTH
plt.rcParams['ytick.major.width'] = LINE_WIDTH
plt.rcParams['ytick.minor.size']  = 2.0*LINE_WIDTH
plt.rcParams['ytick.minor.width'] = 0.5*LINE_WIDTH
plt.rcParams['ytick.labelsize']   = FONT_SIZE

plt.rcParams['font.family']       = 'STIXGeneral'
plt.rcParams['mathtext.fontset']  = 'custom'
plt.rcParams['mathtext.rm']       = 'STIXGeneral:regular'
plt.rcParams['mathtext.it']       = 'STIXGeneral:italic'
plt.rcParams['mathtext.bf']       = 'STIXGeneral:italic:bold'
###################################################################################################


###################################################################################################
# Setting for Plotting
###################################################################################################
n_rows             = 4                 # number of rows in the figure panels
n_cols             = 4                 # number of columns in the figure panels

# fields to plot the projections
#field_toplot       = ('gas','density') # target field
field_toplot       = ('gas','temperature') # target field


Plotting_UNIT_L    = 'pc'              # length unit in the plot
Plotting_UNIT_D    = 'Msun/pc**3'      # density unit in the plot
Plotting_UNIT_T    = 'Myr'             # time unit in the plot

#UNIT_D_todisplay   = r'$\mathit{\rho}\ ({\rm M}_{\odot}{\rm pc}^{\rm -3})$'
UNIT_D_todisplay   = 'T (K)'

axis               = 'z'               # direction for projection
zoom               = 1                 # region for main plot

slab_thickness     = 400.0             # slab thickness for the projection (in plotting unit)
annotated_scale    = 100               # length scale to be annotated in the figure (in plotting unit)
#annotated_radius   = 10.0             # radius of sphere to be annotated in the figure (in plotting unit)

#colormap           = 'viridis'         # ['viridis', 'plasma', 'inferno', 'magma', 'cividis']
#color_lim_min      = 8.e-3             # color bar upper limit for projection (in plotting unit)
#color_lim_max      = 2.2e-2            # color bar lower limit for projection (in plotting unit)

colormap           = 'magma'
color_lim_min      = 5.e2              # color bar upper limit for projection (in plotting unit)
color_lim_max      = 1.1e5             # color bar lower limit for projection (in plotting unit)

buff_size          = 1024              # buffer size for the projection figure

grid_rect_left     =  0.01             # left boundary of the figure panels
grid_rect_bottom   = -0.07             # bottom boundary of the figure panels
grid_rect_width    =  0.90             # width of the figure panels
grid_rect_height   =  1.14             # height of the figure panels

cm                 = 1/2.54            # centimeters in inches
fig_size_x         = 16*cm             # output figure size in the horizontal direction
fig_size_y         = (fig_size_x-1.8*cm)/n_cols*n_rows
fig                = plt.figure( 1,(fig_size_x, fig_size_y), dpi=dpi )
###################################################################################################


###################################################################################################
# Read the data and parameters
###################################################################################################
# path
Data_path    = np.empty( n_rows*n_cols, dtype=object )
Data_path[0] = './Data_000354'
Data_path[1] = './Data_000400'
Data_path[2] = './Data_000500'
Data_path[3] = './Data_000600'
Data_path[4] = '../SSN_momentum/Data_000354'
Data_path[5] = '../SSN_momentum/Data_000400'
Data_path[6] = '../SSN_momentum/Data_000500'
Data_path[7] = '../SSN_momentum/Data_000600'
Data_path[8] = '../../../enzo-dev/run/StarParticle/StarParticleSingleTest/DD0354/data0354'
Data_path[9] = '../../../enzo-dev/run/StarParticle/StarParticleSingleTest/DD0400/data0400'
Data_path[10] = '../../../enzo-dev/run/StarParticle/StarParticleSingleTest/DD0500/data0500'
Data_path[11] = '../../../enzo-dev/run/StarParticle/StarParticleSingleTest/DD0600/data0600'
Data_path[12] = '../../../enzo-dev/run/StarParticle/StarParticleSingleTest/momentum_feedback/DD0354/data0354'
Data_path[13] = '../../../enzo-dev/run/StarParticle/StarParticleSingleTest/momentum_feedback/DD0400/data0400'
Data_path[14] = '../../../enzo-dev/run/StarParticle/StarParticleSingleTest/momentum_feedback/DD0500/data0500'
Data_path[15] = '../../../enzo-dev/run/StarParticle/StarParticleSingleTest/momentum_feedback/DD0600/data0600'

# data ID
Data_idx     = np.empty( n_rows*n_cols, dtype=int )
Data_idx[0]  = 0
Data_idx[1]  = 1
Data_idx[2]  = 2
Data_idx[3]  = 3
Data_idx[4]  = 4
Data_idx[5]  = 5
Data_idx[6]  = 6
Data_idx[7]  = 7
###################################################################################################



def main() -> None:

    # Create the AxesGrid
    grid = AxesGrid( fig, rect=[grid_rect_left, grid_rect_bottom, grid_rect_width, grid_rect_height],
                     nrows_ncols=(n_rows, n_cols),
                     axes_pad=(0.05,0.05), label_mode="L", share_all=True,
                     cbar_location="right", cbar_mode="single", cbar_size="2%", cbar_pad="1%" )

    # Loop for each panel
    for panel_idx in range(0, n_rows*n_cols, 1):

        ###################################################################################################
        # Load the dataset
        #ds = yt.load( Data_path[panel_idx]+'Data_%06d'%Data_idx[panel_idx] )
        ds = yt.load( Data_path[panel_idx] )

        ###################################################################################################
        # Annotated information
        #Data_Annotated_radius = ds.arr( annotated_radius, Plotting_UNIT_L )
        Data_Annotated_scale  = ds.arr( annotated_scale,  Plotting_UNIT_L )

        ###################################################################################################
        # Region to plot
        Data_Center_x         = ds.domain_center[0].in_units('code_length').d
        Data_Center_y         = ds.domain_center[1].in_units('code_length').d
        Data_Center_z         = ds.domain_center[2].in_units('code_length').d
        Data_Center           = ds.arr( [Data_Center_x, Data_Center_y, Data_Center_z], 'code_length' )
        Data_Width            = ds.domain_width[0]/zoom
        Data_slab_thickness   = ds.quan( slab_thickness , Plotting_UNIT_L )
        Data_loc              = Data_Center
        Data_corner_L         = Data_loc - ds.arr( [0.5*Data_Width, 0.5*Data_Width, 0.5*Data_slab_thickness] )
        Data_corner_R         = Data_loc + ds.arr( [0.5*Data_Width, 0.5*Data_Width, 0.5*Data_slab_thickness] )
        Data_region           = ds.box( Data_corner_L, Data_corner_R )

        ###################################################################################################
        # Projection plot
        #proj = yt.ProjectionPlot( ds, axis, field_toplot,
                                  #weight_field=('index','ones'),
                                  #data_source=Data_region,
                                  #center=Data_loc, width=(Data_Width, Data_Width),
                                  #buff_size=(buff_size, buff_size) )

        proj = yt.ProjectionPlot( ds, axis, field_toplot,
                                  weight_field=('density'),
                                  data_source=Data_region,
                                  center=Data_loc, width=(Data_Width, Data_Width),
                                  buff_size=(buff_size, buff_size) )

        # Set the font
        proj.set_font( { "size":FONT_SIZE, "math_fontfamily":'custom' } )

        # Set the axes unit
        proj.set_axes_unit( Plotting_UNIT_L )

        # Set the field unit
        #proj.set_unit( field_toplot, Plotting_UNIT_D )
        proj.set_unit( field_toplot, 'K' )

        # Set the colorbar limits
        proj.set_zlim( field_toplot, color_lim_min, color_lim_max )

        # Set the colormap
        proj.set_cmap( field_toplot, colormap )

        # Set the x label
        proj.set_xlabel( '' )

        # Set the y label
        proj.set_ylabel( '' )

        # Set the colorbar label
        proj.set_colorbar_label( field_toplot, UNIT_D_todisplay )

        # Annotate the time
        proj.annotate_timestamp( time=True, redshift=False, corner='upper_right',
                                 time_format='$t$ = {time:.1f} {units}', time_unit=Plotting_UNIT_T,
                                 redshift_format='$z$ = {redshift:.1f}',
                                 text_args={ "size":FONT_SIZE, "color":"white" } )

        # Annotate the text for title
        if panel_idx == 0:
            proj.annotate_text( [0.025, 0.80], '\nGAMER\nenergy',
               	                coord_system="axis",
                       	        text_args={ "size":1.25*FONT_SIZE, "color":"white", "weight":"normal", "bbox":dict(boxstyle="round",ec='white',fc='white',alpha=0.0) } )

        if panel_idx == 4:
            proj.annotate_text( [0.025, 0.80], '\nGAMER\nmomentum',
                                coord_system="axis",
                                text_args={ "size":1.25*FONT_SIZE, "color":"white", "weight":"normal", "bbox":dict(boxstyle="round",ec='white',fc='white',alpha=0.0) } )

        if panel_idx == 8:
            proj.annotate_text( [0.025, 0.80], '\nEnzo\nenergy',
                                coord_system="axis",
                                text_args={ "size":1.25*FONT_SIZE, "color":"white", "weight":"normal", "bbox":dict(boxstyle="round",ec='white',fc='white',alpha=0.0) } )

        if panel_idx == 12:
            proj.annotate_text( [0.025, 0.80], '\nEnzo\nmomentum',
                                coord_system="axis",
                                text_args={ "size":1.25*FONT_SIZE, "color":"white", "weight":"normal", "bbox":dict(boxstyle="round",ec='white',fc='white',alpha=0.0) } )

        # Annotate the scale
        if panel_idx%4 == 3:
            proj.annotate_scale( corner='upper_right',
                                 coeff=Data_Annotated_scale.in_units(Plotting_UNIT_L).d, unit=Plotting_UNIT_L,
                                 text_args={'size':0.8*FONT_SIZE}, pos=[0.2, 0.1], coord_system='axis' )

        # Annotate the sphere
        #proj.annotate_sphere( center=Data_Center, radius=Data_Annotated_radius,
        #                      circle_args={ "color":"dimgrey", "linestyle":"--", "linewidth":0.5*LINE_WIDTH } )

        # Annotate the text for the sphere
        #proj.annotate_text( (1.2*Data_Annotated_radius, -1.2*Data_Annotated_radius), r'$r$',
        #                    coord_system="plot",
        #                    text_args={ "size":FONT_SIZE, "color":"dimgrey", "weight":"normal" } )

        ###################################################################################################
        # Set the plots to the grid
        plot        = proj.plots[field_toplot]
        plot.figure = fig
        plot.axes   = grid[panel_idx].axes
        plot.cax    = grid.cbar_axes[panel_idx]

        # Set the ticks
        grid[panel_idx].axes.tick_params( which='both', color='black', size=0.0 )
        grid[panel_idx].axes.tick_params( which='both', left=False, right=False, bottom=False, top=False )
        grid[panel_idx].axes.tick_params( labelleft=False, labelright=False, labelbottom=False, labeltop=False )

        proj._setup_plots()

        ###################################################################################################

    # Set y ticks
    grid.axes_llc.set_yticks( [] )

    # Output the figure to files
    fig.set_dpi( dpi )
    fig.set_size_inches( fig_size_x, fig_size_y )

    fig.savefig( "fig_standard_paper_figure_2D_t.png", dpi=dpi )
    fig.savefig( "fig_standard_paper_figure_2D_t.pdf", dpi=dpi )


if __name__ == '__main__':
    main()
