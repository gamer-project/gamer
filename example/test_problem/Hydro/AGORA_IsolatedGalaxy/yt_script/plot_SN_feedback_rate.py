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


# get the creation time of the new stars and explosion time of teh SNe
ad            = ds.all_data()

assert ad[ "all", "ParSNIINxtE" ].d.astype('int').max() == 10001
feedback_table  = np.genfromtxt( '../FeedbackYieldTable_Resolved_SNeII_N000001', delimiter=None, comments='#', names=['explosion_time', 'progenitor_mass', 'ejected_mass', 'ejected_metals', 'ejected_energy'], dtype=None, encoding=None )
feedback_record = np.genfromtxt( '../Record__FB_Resolved_SNeII', delimiter=None, comments='#', names=['Rank', 'TID', 'lv', 'TimeOld', 'TimeNew', 'SNII_Time_first', 'SNII_Time_last', 'SNII_AccumNum', 'SNII_Energy', 'SNII_Mass', 'SNII_Metal', 'FB_Diameter', 'FB_Flu_Mass', 'Par_PUID', 'Par_NxtE', 'Par_Mass', 'Par_PosX', 'Par_PosY', 'Par_PosZ', 'Par_VelX', 'Par_VelY', 'Par_VelZ', 'Par_MetalMass', 'Flu_Dens', 'Flu_MomX', 'Flu_MomY', 'Flu_MomZ', 'Flu_Engy', 'Flu_Metal'], dtype=None, encoding=None )

star_mass     = ad[ "new_star", "ParMass" ].in_units( "Msun" ) + ds.quan( feedback_table['ejected_mass'], 'Msun' )
creation_time = ad[ "new_star", "ParCreTime" ].in_units( "Myr" )
SNe_mass      = ds.arr( feedback_record['Par_Mass'], 'code_mass' ).in_units( "Msun" )
SNe_expl_time = ds.arr( feedback_record['SNII_Time_last'], 'code_time' ).in_units( "Myr" )

print( 'Total number of stars = %d'%len(creation_time) )
print( 'Total number of SNe   = %d'%len(SNe_expl_time) )


# bin the data
t_start        = 0.0
t_end          = ds.current_time.in_units( "Myr" ).d
t_bin          = np.linspace( start=t_start, stop=t_end, num=nbin+1 )
time           = 0.5*( t_bin[:-1] + t_bin[1:] )
star_upper_idx = np.digitize( creation_time.in_units( "Myr" ).d, bins=t_bin, right=True )
SNe_upper_idx  = np.digitize( SNe_expl_time.in_units( "Myr" ).d, bins=t_bin, right=True )


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
SNDelayTime = ds.quan( feedback_table['explosion_time'], 'Myr' ).d


# plot
plt.plot( time,             sfr,                  label='Stars' )
plt.plot( time,             SNr,                  label='SNe'   )
plt.plot( time,             SNr*StarsPerSN, '--', label=r'SNe, $\times$ %.2f'%(StarsPerSN) )
plt.plot( time-SNDelayTime, SNr*StarsPerSN, '--', label=r'SNe, $\times$ %.2f, shifted %.1f Myr'%(StarsPerSN, SNDelayTime) )
plt.ylim( 0.0, 1.0e1 )
plt.legend()
plt.xlabel( "$\mathrm{t\ [Myr]}$",        fontsize="large" )
plt.ylabel( "$\mathrm{[M_\odot yr^{-1}]}$", fontsize="large" )


# show/save figure
plt.savefig( fileout+".png", bbox_inches="tight", pad_inches=0.05, dpi=dpi )
#plt.show()
