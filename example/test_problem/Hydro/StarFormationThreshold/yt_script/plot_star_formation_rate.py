# ref: https://yt-project.org/docs/dev/cookbook/calculating_information.html#using-particle-filters-to-calculate-star-formation-rates

import yt
import numpy as np
from yt.data_objects.particle_filters import add_particle_filter
from matplotlib import pyplot as plt


filein  = "../Data_000050"
fileout = "fig__star_formation_rate"
nbin    = 50
dpi     = 150


# load data
ds = yt.load( filein )


# define the particle filter for the newly formed stars
def new_star( pfilter, data ):
   filter = data[ "all", "ParCreTime" ] > 0
   return filter

add_particle_filter( "new_star", function=new_star, filtered_type="all", requires=["ParCreTime"] )
ds.add_particle_filter( "new_star" )


# get the mass and creation time of the new stars
ad            = ds.all_data()
mass          = ad[ "new_star", "ParMass"    ].in_units( "Msun" )
creation_time = ad[ "new_star", "ParCreTime" ].in_units( "Myr" )


# bin the data
t_start   = 0.0
t_end     = ds.current_time.in_units( "Myr" )
t_bin     = np.linspace( start=t_start, stop=t_end, num=nbin+1 )
upper_idx = np.digitize( creation_time, bins=t_bin, right=True )
time      = 0.5*( t_bin[:-1] + t_bin[1:] )

assert np.all( upper_idx > 0 ) and np.all( upper_idx < len(t_bin) ), "incorrect upper_idx !!"


# calculate the star formation rate
Myr2yr = 1.0e6
sfr    = np.array(  [ mass[upper_idx == j+1].sum() / ( (t_bin[j+1] - t_bin[j])*Myr2yr )
                    for j in range(len(time)) ]  )
sfr[sfr == 0] = np.nan


# plot
plt.plot( time, sfr )
plt.ylim( 0.0, 1.0e1 )
plt.xlabel( "$\mathrm{t\ [Myr]}$",               fontsize="large" )
plt.ylabel( "$\mathrm{SFR\ [M_\odot yr^{-1}]}$", fontsize="large" )


# show/save figure
plt.savefig( fileout+".png", bbox_inches="tight", pad_inches=0.05, dpi=dpi )
#plt.show()
