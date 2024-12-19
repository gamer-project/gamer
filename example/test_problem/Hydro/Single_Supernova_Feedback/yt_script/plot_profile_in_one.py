import argparse
import sys
import yt
import matplotlib.pyplot as plt
import numpy as np

# load the command-line parameters
parser = argparse.ArgumentParser( description='Plot the gas slices and projections' )

parser.add_argument( '-i1', action='store', required=False, type=str, dest='prefix_1',
                     help='path prefix [%(default)s]', default='../' )
parser.add_argument( '-i2', action='store', required=False, type=str, dest='prefix_2',
                     help='path prefix [%(default)s]', default='../' )
#parser.add_argument( '-s', action='store', required=True,  type=int, dest='idx_start',
#                     help='first data index' )
#parser.add_argument( '-e', action='store', required=True,  type=int, dest='idx_end',
#                     help='last data index' )
#parser.add_argument( '-d', action='store', required=False, type=int, dest='didx',
#                     help='delta data index [%(default)d]', default=1 )

args=parser.parse_args()

# take note
print( '\nCommand-line arguments:' )
print( '-------------------------------------------------------------------' )
for t in range( len(sys.argv) ):
   print str(sys.argv[t]),
print( '' )
print( '-------------------------------------------------------------------\n' )


#idx_start   = args.idx_start
#idx_end     = args.idx_end
#didx        = args.didx
prefix_1    = args.prefix_1
prefix_2    = args.prefix_2

center_mode = [0.50390625, 0.50390625, 0.50390625]
dpi         = 300


yt.enable_parallelism()


#load data
ts_1 = yt.load( prefix_1 + "/Data_000352" )
es_1 = yt.load( prefix_2 + "/DD0352/data0352" )
ts_2 = yt.load( prefix_1 + "/Data_000354" )
es_2 = yt.load( prefix_2 + "/DD0354/data0354" )
ts_3 = yt.load( prefix_1 + "/Data_000394" )
es_3 = yt.load( prefix_2 + "/DD0394/data0394" )
ts_4 = yt.load( prefix_1 + "/Data_000490" )
es_4 = yt.load( prefix_2 + "/DD0490/data0490" )
ts_5 = yt.load( prefix_1 + "/Data_000600" )
es_5 = yt.load( prefix_2 + "/DD0600/data0600" )

#add new field for weighting
def _density_square( field, data ):
   return data["density"]**2

es_1.add_field( ("gas", "density_square"), function=_density_square, sampling_type="cell", units="Msun**2/pc**6" )
ts_1.add_field( ("gas", "density_square"), function=_density_square, sampling_type="cell", units="Msun**2/pc**6" )
es_2.add_field( ("gas", "density_square"), function=_density_square, sampling_type="cell", units="Msun**2/pc**6" )
ts_2.add_field( ("gas", "density_square"), function=_density_square, sampling_type="cell", units="Msun**2/pc**6" )
es_3.add_field( ("gas", "density_square"), function=_density_square, sampling_type="cell", units="Msun**2/pc**6" )
ts_3.add_field( ("gas", "density_square"), function=_density_square, sampling_type="cell", units="Msun**2/pc**6" )
es_4.add_field( ("gas", "density_square"), function=_density_square, sampling_type="cell", units="Msun**2/pc**6" )
ts_4.add_field( ("gas", "density_square"), function=_density_square, sampling_type="cell", units="Msun**2/pc**6" )
es_5.add_field( ("gas", "density_square"), function=_density_square, sampling_type="cell", units="Msun**2/pc**6" )
ts_5.add_field( ("gas", "density_square"), function=_density_square, sampling_type="cell", units="Msun**2/pc**6" )

#create sphere
GAMER_sphere_1 = ts_1.sphere( center_mode, 0.5*ts_1.domain_width.to_value().max() )
Enzo_sphere_1  = es_1.sphere( center_mode, 0.5*es_1.domain_width.to_value().max() )
GAMER_sphere_2 = ts_2.sphere( center_mode, 0.5*ts_2.domain_width.to_value().max() )
Enzo_sphere_2  = es_2.sphere( center_mode, 0.5*es_2.domain_width.to_value().max() )
GAMER_sphere_3 = ts_3.sphere( center_mode, 0.5*ts_3.domain_width.to_value().max() )
Enzo_sphere_3  = es_3.sphere( center_mode, 0.5*es_3.domain_width.to_value().max() )
GAMER_sphere_4 = ts_4.sphere( center_mode, 0.5*ts_4.domain_width.to_value().max() )
Enzo_sphere_4  = es_4.sphere( center_mode, 0.5*es_4.domain_width.to_value().max() )
GAMER_sphere_5 = ts_5.sphere( center_mode, 0.5*ts_5.domain_width.to_value().max() )
Enzo_sphere_5  = es_5.sphere( center_mode, 0.5*es_5.domain_width.to_value().max() )

   #create profiles

#density
#GAMER_profile_1 = yt.create_profile( GAMER_sphere_1, 'radius', fields='density', weight_field='cell_volume', n_bins=256 )
#Enzo_profile_1  = yt.create_profile( Enzo_sphere_1,  'radius', fields='density', weight_field='cell_volume', n_bins=256 )
#GAMER_profile_2 = yt.create_profile( GAMER_sphere_2, 'radius', fields='density', weight_field='cell_volume', n_bins=256 )
#Enzo_profile_2  = yt.create_profile( Enzo_sphere_2,  'radius', fields='density', weight_field='cell_volume', n_bins=256 )
#GAMER_profile_3 = yt.create_profile( GAMER_sphere_3, 'radius', fields='density', weight_field='cell_volume', n_bins=256 )
#Enzo_profile_3  = yt.create_profile( Enzo_sphere_3,  'radius', fields='density', weight_field='cell_volume', n_bins=256 )
#GAMER_profile_4 = yt.create_profile( GAMER_sphere_4, 'radius', fields='density', weight_field='cell_volume', n_bins=256 )
#Enzo_profile_4  = yt.create_profile( Enzo_sphere_4,  'radius', fields='density', weight_field='cell_volume', n_bins=256 )
#GAMER_profile_5 = yt.create_profile( GAMER_sphere_5, 'radius', fields='density', weight_field='cell_volume', n_bins=256 )
#Enzo_profile_5  = yt.create_profile( Enzo_sphere_5,  'radius', fields='density', weight_field='cell_volume', n_bins=256 )

#remove zero element from profile
#array_GAMER_r_1     = GAMER_profile_1.x.in_units('pc').d
#array_GAMER_dens_1  = GAMER_profile_1['density'].in_units('Msun/pc**3').d
#array_GAMER_r_2     = GAMER_profile_2.x.in_units('pc').d
#array_GAMER_dens_2  = GAMER_profile_2['density'].in_units('Msun/pc**3').d
#array_GAMER_r_3     = GAMER_profile_3.x.in_units('pc').d
#array_GAMER_dens_3  = GAMER_profile_3['density'].in_units('Msun/pc**3').d
#array_GAMER_r_4     = GAMER_profile_4.x.in_units('pc').d
#array_GAMER_dens_4  = GAMER_profile_4['density'].in_units('Msun/pc**3').d
#array_GAMER_r_5     = GAMER_profile_5.x.in_units('pc').d
#array_GAMER_dens_5  = GAMER_profile_5['density'].in_units('Msun/pc**3').d

#array_Enzo_r_1     = Enzo_profile_1.x.in_units('pc').d
#array_Enzo_dens_1  = Enzo_profile_1['density'].in_units('Msun/pc**3').d
#array_Enzo_r_2     = Enzo_profile_2.x.in_units('pc').d
#array_Enzo_dens_2  = Enzo_profile_2['density'].in_units('Msun/pc**3').d
#array_Enzo_r_3     = Enzo_profile_3.x.in_units('pc').d
#array_Enzo_dens_3  = Enzo_profile_3['density'].in_units('Msun/pc**3').d
#array_Enzo_r_4     = Enzo_profile_4.x.in_units('pc').d
#array_Enzo_dens_4  = Enzo_profile_4['density'].in_units('Msun/pc**3').d
#array_Enzo_r_5     = Enzo_profile_5.x.in_units('pc').d
#array_Enzo_dens_5  = Enzo_profile_5['density'].in_units('Msun/pc**3').d

#print( array_GAMER_radius_1 )
#print( array_GAMER_density_1 )

#array_GAMER_r_1    = array_GAMER_r_1[ array_GAMER_dens_1>0.0 ]
#array_GAMER_dens_1 = array_GAMER_dens_1[ array_GAMER_dens_1>0.0 ]
#array_GAMER_r_2    = array_GAMER_r_2[ array_GAMER_dens_2>0.0 ]
#array_GAMER_dens_2 = array_GAMER_dens_2[ array_GAMER_dens_2>0.0 ]
#array_GAMER_r_3    = array_GAMER_r_3[ array_GAMER_dens_3>0.0 ]
#array_GAMER_dens_3 = array_GAMER_dens_3[ array_GAMER_dens_3>0.0 ]
#array_GAMER_r_4    = array_GAMER_r_4[ array_GAMER_dens_4>0.0 ]
#array_GAMER_dens_4 = array_GAMER_dens_4[ array_GAMER_dens_4>0.0 ]
#array_GAMER_r_5    = array_GAMER_r_5[ array_GAMER_dens_5>0.0 ]
#array_GAMER_dens_5 = array_GAMER_dens_5[ array_GAMER_dens_5>0.0 ]

#array_Enzo_r_1    = array_Enzo_r_1[ array_Enzo_dens_1>0.0 ]
#array_Enzo_dens_1 = array_Enzo_dens_1[ array_Enzo_dens_1>0.0 ]
#array_Enzo_r_2    = array_Enzo_r_2[ array_Enzo_dens_2>0.0 ]
#array_Enzo_dens_2 = array_Enzo_dens_2[ array_Enzo_dens_2>0.0 ]
#array_Enzo_r_3    = array_Enzo_r_3[ array_Enzo_dens_3>0.0 ]
#array_Enzo_dens_3 = array_Enzo_dens_3[ array_Enzo_dens_3>0.0 ]
#array_Enzo_r_4    = array_Enzo_r_4[ array_Enzo_dens_4>0.0 ]
#array_Enzo_dens_4 = array_Enzo_dens_4[ array_Enzo_dens_4>0.0 ]
#array_Enzo_r_5    = array_Enzo_r_5[ array_Enzo_dens_5>0.0 ]
#array_Enzo_dens_5 = array_Enzo_dens_5[ array_Enzo_dens_5>0.0 ]

#print( array_GAMER_radius_1_remove_zero )
#print( array_GAMER_density_1_remove_zero )

#temperature
GAMER_profile_1 = yt.create_profile( GAMER_sphere_1, 'radius', fields='temperature', weight_field='density_square', n_bins=256 )
Enzo_profile_1  = yt.create_profile( Enzo_sphere_1,  'radius', fields='temperature', weight_field='density_square', n_bins=256 )
GAMER_profile_2 = yt.create_profile( GAMER_sphere_2, 'radius', fields='temperature', weight_field='density_square', n_bins=256 )
Enzo_profile_2  = yt.create_profile( Enzo_sphere_2,  'radius', fields='temperature', weight_field='density_square', n_bins=256 )
GAMER_profile_3 = yt.create_profile( GAMER_sphere_3, 'radius', fields='temperature', weight_field='density_square', n_bins=256 )
Enzo_profile_3  = yt.create_profile( Enzo_sphere_3,  'radius', fields='temperature', weight_field='density_square', n_bins=256 )
GAMER_profile_4 = yt.create_profile( GAMER_sphere_4, 'radius', fields='temperature', weight_field='density_square', n_bins=256 )
Enzo_profile_4  = yt.create_profile( Enzo_sphere_4,  'radius', fields='temperature', weight_field='density_square', n_bins=256 )
GAMER_profile_5 = yt.create_profile( GAMER_sphere_5, 'radius', fields='temperature', weight_field='density_square', n_bins=256 )
Enzo_profile_5  = yt.create_profile( Enzo_sphere_5,  'radius', fields='temperature', weight_field='density_square', n_bins=256 )

array_GAMER_r_1     = GAMER_profile_1.x.in_units('pc').d
array_GAMER_temp_1  = GAMER_profile_1['temperature'].in_units('K').d
array_GAMER_r_2     = GAMER_profile_2.x.in_units('pc').d
array_GAMER_temp_2  = GAMER_profile_2['temperature'].in_units('K').d
array_GAMER_r_3     = GAMER_profile_3.x.in_units('pc').d
array_GAMER_temp_3  = GAMER_profile_3['temperature'].in_units('K').d
array_GAMER_r_4     = GAMER_profile_4.x.in_units('pc').d
array_GAMER_temp_4  = GAMER_profile_4['temperature'].in_units('K').d
array_GAMER_r_5     = GAMER_profile_5.x.in_units('pc').d
array_GAMER_temp_5  = GAMER_profile_5['temperature'].in_units('K').d

array_Enzo_r_1     = Enzo_profile_1.x.in_units('pc').d
array_Enzo_temp_1  = Enzo_profile_1['temperature'].in_units('K').d
array_Enzo_r_2     = Enzo_profile_2.x.in_units('pc').d
array_Enzo_temp_2  = Enzo_profile_2['temperature'].in_units('K').d
array_Enzo_r_3     = Enzo_profile_3.x.in_units('pc').d
array_Enzo_temp_3  = Enzo_profile_3['temperature'].in_units('K').d
array_Enzo_r_4     = Enzo_profile_4.x.in_units('pc').d
array_Enzo_temp_4  = Enzo_profile_4['temperature'].in_units('K').d
array_Enzo_r_5     = Enzo_profile_5.x.in_units('pc').d
array_Enzo_temp_5  = Enzo_profile_5['temperature'].in_units('K').d

array_GAMER_r_1    = array_GAMER_r_1[ array_GAMER_temp_1>0.0 ]
array_GAMER_temp_1 = array_GAMER_temp_1[ array_GAMER_temp_1>0.0 ]
array_GAMER_r_2    = array_GAMER_r_2[ array_GAMER_temp_2>0.0 ]
array_GAMER_temp_2 = array_GAMER_temp_2[ array_GAMER_temp_2>0.0 ]
array_GAMER_r_3    = array_GAMER_r_3[ array_GAMER_temp_3>0.0 ]
array_GAMER_temp_3 = array_GAMER_temp_3[ array_GAMER_temp_3>0.0 ]
array_GAMER_r_4    = array_GAMER_r_4[ array_GAMER_temp_4>0.0 ]
array_GAMER_temp_4 = array_GAMER_temp_4[ array_GAMER_temp_4>0.0 ]
array_GAMER_r_5    = array_GAMER_r_5[ array_GAMER_temp_5>0.0 ]
array_GAMER_temp_5 = array_GAMER_temp_5[ array_GAMER_temp_5>0.0 ]

array_Enzo_r_1    = array_Enzo_r_1[ array_Enzo_temp_1>0.0 ]
array_Enzo_temp_1 = array_Enzo_temp_1[ array_Enzo_temp_1>0.0 ]
array_Enzo_r_2    = array_Enzo_r_2[ array_Enzo_temp_2>0.0 ]
array_Enzo_temp_2 = array_Enzo_temp_2[ array_Enzo_temp_2>0.0 ]
array_Enzo_r_3    = array_Enzo_r_3[ array_Enzo_temp_3>0.0 ]
array_Enzo_temp_3 = array_Enzo_temp_3[ array_Enzo_temp_3>0.0 ]
array_Enzo_r_4    = array_Enzo_r_4[ array_Enzo_temp_4>0.0 ]
array_Enzo_temp_4 = array_Enzo_temp_4[ array_Enzo_temp_4>0.0 ]
array_Enzo_r_5    = array_Enzo_r_5[ array_Enzo_temp_5>0.0 ]
array_Enzo_temp_5 = array_Enzo_temp_5[ array_Enzo_temp_5>0.0 ]

# create the figure
fig = plt.figure()
ax  = fig.add_subplot(111)
ax.xaxis.get_ticklocs(minor=True)
ax.yaxis.get_ticklocs(minor=True)
ax.minorticks_on()
ax.xaxis.set_ticks_position('both')
ax.yaxis.set_ticks_position('both')    
ax.tick_params(which='both',direction='in')

# Add a title for each column in the legend
def add_col_titles( ax, num_rows, num_cols, col_title ):
    # original handles and labels
    handles, labels = ax.get_legend_handles_labels()

    # empty plot
    handles_title      = [ ax.plot( [], marker="", ls="" )[0] ]*(num_cols)

    # Add the title
    handles_with_title = []
    labels_with_title  = []
    for j in range(num_cols):
        handles_with_title += [handles_title[j]] + handles[j*num_rows:(j+1)*num_rows]
        labels_with_title  +=     [col_title[j]] +  labels[j*num_rows:(j+1)*num_rows]

    # plot the legend
    legend = ax.legend(handles_with_title, labels_with_title, ncol=num_cols, loc='upper right', columnspacing=1.0, handlelength=5.5)
    for vpack in legend._legend_handle_box.get_children():
        for hpack in vpack.get_children()[:1]:
            hpack.get_children()[0].set_width(0)

# plot the profiles

# line styles
LStyle_Dot        = [1, 2]
LStyle_Dash       = [4, 2]
LStyle_DashDot    = [4, 2, 1, 2]
LStyle_DotDashDot = [1, 2, 4, 2, 1, 2]

#density
#ax.plot( array_GAMER_r_1, array_GAMER_dens_1, c=plt.cm.autumn_r( 1.00 ), linewidth=2, label=' ', zorder=2.5 )
#ax.plot( array_GAMER_r_2, array_GAMER_dens_2, c=plt.cm.autumn_r( 0.82 ), linewidth=2, label=' ', zorder=2.5 )
#ax.plot( array_GAMER_r_3, array_GAMER_dens_3, c=plt.cm.autumn_r( 0.64 ), linewidth=2, label=' ', zorder=2.5 )
#ax.plot( array_GAMER_r_4, array_GAMER_dens_4, c=plt.cm.autumn_r( 0.46 ), linewidth=2, label=' ', zorder=2.5 )
#ax.plot( array_GAMER_r_5, array_GAMER_dens_5, c=plt.cm.autumn_r( 0.28 ), linewidth=2, label=' ', zorder=2.5 )
#ax.plot( array_Enzo_r_1, array_Enzo_dens_1, c=plt.cm.winter_r( 1.00 ), dashes=LStyle_Dash, linewidth=2.5, label='t = 0.0 Myr', zorder=2 )
#ax.plot( array_Enzo_r_2, array_Enzo_dens_2, c=plt.cm.winter_r( 0.82 ), dashes=LStyle_Dash, linewidth=2.5, label='t = 3.5 Myr', zorder=2 )
#ax.plot( array_Enzo_r_3, array_Enzo_dens_3, c=plt.cm.winter_r( 0.64 ), dashes=LStyle_Dash, linewidth=2.5, label='t = 4.0 Myr', zorder=2 )
#ax.plot( array_Enzo_r_4, array_Enzo_dens_4, c=plt.cm.winter_r( 0.46 ), dashes=LStyle_Dash, linewidth=2.5, label='t = 5.0 Myr', zorder=2 )
#ax.plot( array_Enzo_r_5, array_Enzo_dens_5, c=plt.cm.winter_r( 0.28 ), dashes=LStyle_Dash, linewidth=2.5, label='t = 6.0 Myr', zorder=2 )

#add_col_titles( ax, 5, 2, ["GAMER","Enzo"] )

#ax.set_yscale( 'log' )
#ax.set_xscale( 'log' )
#ax.set_xlim( 3.0e0, 2.0e2 )
#ax.set_ylim( 2.0e-5, 1.0e0 )
#ax.set_xlabel( 'r (pc)' )
#ax.set_ylabel( r'$\rho\ (M_\odot/{\rm pc}^3)$' )

#temperature
ax.plot( array_GAMER_r_1, array_GAMER_temp_1, c=plt.cm.autumn_r( 1.00 ), linewidth=2, label=' ', zorder=2.5 ) 
ax.plot( array_GAMER_r_2, array_GAMER_temp_2, c=plt.cm.autumn_r( 0.82 ), linewidth=2, label=' ', zorder=2.5 )
ax.plot( array_GAMER_r_3, array_GAMER_temp_3, c=plt.cm.autumn_r( 0.64 ), linewidth=2, label=' ', zorder=2.5 )
ax.plot( array_GAMER_r_4, array_GAMER_temp_4, c=plt.cm.autumn_r( 0.46 ), linewidth=2, label=' ', zorder=2.5 )
ax.plot( array_GAMER_r_5, array_GAMER_temp_5, c=plt.cm.autumn_r( 0.28 ), linewidth=2, label=' ', zorder=2.5 )
ax.plot( array_Enzo_r_1, array_Enzo_temp_1, c=plt.cm.winter_r( 1.00 ), dashes=LStyle_Dash, linewidth=2.5, label='t = 0.0 Myr', zorder=2 )
ax.plot( array_Enzo_r_2, array_Enzo_temp_2, c=plt.cm.winter_r( 0.82 ), dashes=LStyle_Dash, linewidth=2.5, label='t = 3.5 Myr', zorder=2 )
ax.plot( array_Enzo_r_3, array_Enzo_temp_3, c=plt.cm.winter_r( 0.64 ), dashes=LStyle_Dash, linewidth=2.5, label='t = 4.0 Myr', zorder=2 )
ax.plot( array_Enzo_r_4, array_Enzo_temp_4, c=plt.cm.winter_r( 0.46 ), dashes=LStyle_Dash, linewidth=2.5, label='t = 5.0 Myr', zorder=2 )
ax.plot( array_Enzo_r_5, array_Enzo_temp_5, c=plt.cm.winter_r( 0.28 ), dashes=LStyle_Dash, linewidth=2.5, label='t = 6.0 Myr', zorder=2 )

add_col_titles( ax, 5, 2, ["GAMER","Enzo"] )

ax.set_yscale( 'log' )
ax.set_xscale( 'log' )
ax.set_xlim( 3.0e0, 2.0e2 )
ax.set_ylim( 4.0e1, 2.0e9 )
ax.set_xlabel( 'r (pc)' )
ax.set_ylabel( 'T (K)' )

# save the figure
plt.tight_layout()
fig.savefig( 'temperature_profiles_all_time.pdf' )

