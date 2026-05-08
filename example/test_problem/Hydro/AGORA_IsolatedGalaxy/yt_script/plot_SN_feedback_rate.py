# ref: https://yt-project.org/docs/dev/cookbook/calculating_information.html#using-particle-filters-to-calculate-star-formation-rates

import yt
import numpy as np
from yt.data_objects.particle_filters import add_particle_filter
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt


filein  = "../Data_000050"
fileout = "fig__SN_feedback_rate"
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


# define the particle filter for the exploded SNe
def exploded_SNe( pfilter, data ):
   filter = data[ "all", "ParSNIITime" ] <= 0 # ParSNIITime is set as negative value of explosion time for exploded particle
   return filter

add_particle_filter( "exploded_SNe", function=exploded_SNe, filtered_type="all", requires=["ParSNIITime"] )
ds.add_particle_filter( "exploded_SNe" )


# get the creation time of the new stars and explosion time of teh SNe
ad            = ds.all_data()

star_mass     = ad[ "new_star", "ParMass" ].in_units( "Msun" )
creation_time = ad[ "new_star", "ParCreTime" ].in_units( "Myr" )

SNe_mass      = ad[ "exploded_SNe", "ParMass" ].in_units( "Msun" )
SNe_expl_time = -1.0*( ad[ "exploded_SNe", "ParSNIITime" ]*ds.units.code_time ).in_units( "Myr" ) # ParSNIITime is set as negative value of explosion time for exploded particle

print( 'Total number of stars = %d'%len(creation_time) )
print( 'Total number of SNe   = %d'%len(SNe_expl_time) )


# bin the data
t_start        = 0.0
t_end          = ds.current_time.in_units( "Myr" )
t_bin          = np.linspace( start=t_start, stop=t_end, num=nbin+1 )
time           = 0.5*( t_bin[:-1] + t_bin[1:] )
star_upper_idx = np.digitize( creation_time, bins=t_bin, right=True )
SNe_upper_idx  = np.digitize( SNe_expl_time, bins=t_bin, right=True )


assert np.all( star_upper_idx > 0 ) and np.all( star_upper_idx < len(t_bin) ), "incorrect star_upper_idx !!"
assert np.all( SNe_upper_idx > 0 )  and np.all( SNe_upper_idx < len(t_bin) ),  "incorrect SNe_upper_idx !!"


# calculate the star formation and SNe mass rate
Myr2yr = 1.0e6
sfr    = np.array(  [ star_mass[star_upper_idx == j+1].sum() / ( (t_bin[j+1] - t_bin[j])*Myr2yr ) for j in range(len(time)) ]  )
sfr[sfr == 0] = np.nan

SNr     = np.array(  [ SNe_mass[SNe_upper_idx == j+1].sum()  / ( (t_bin[j+1] - t_bin[j])*Myr2yr ) for j in range(len(time)) ]  )
SNr[SNr == 0] = np.nan


# calulate the conversion factor
StarsPerSN  = 1.0/(ds.parameters['FB_ResolvedSNeII_NPerMass']*ds.parameters['SF_CreateStar_MinStarMass'])
SNDelayTime = ds.quan( ds.parameters['FB_ResolvedSNeII_DelayTime'], 'code_time' ).in_units('Myr').d


# plot
plt.plot( time,               sfr,                  label='Stars' )
plt.plot( time,               SNr,                  label='SNe'   )
plt.plot( time,               SNr*StarsPerSN, '--', label=r'SNe, $\times$ %.2f'%(StarsPerSN) )
plt.plot( time.d-SNDelayTime, SNr*StarsPerSN, '--', label=r'SNe, $\times$ %.2f, shifted %.1f Myr'%(StarsPerSN, SNDelayTime) )
plt.ylim( 0.0, 1.0e1 )
plt.legend()
plt.xlabel( "$\mathrm{t\ [Myr]}$",        fontsize="large" )
plt.ylabel( "$\mathrm{[M_\odot yr^{-1}]}$", fontsize="large" )


# show/save figure
plt.savefig( fileout+".png", bbox_inches="tight", pad_inches=0.05, dpi=dpi )
#plt.show()
