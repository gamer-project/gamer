import yt
import numpy as np
from yt.data_objects.particle_filters import add_particle_filter
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt


filein  = '../Data_000050'
fileout = 'fig__SN_feedback_rate'
nbin    = 50
dpi     = 150
feedback_table_file_prefix = '../FeedbackTables/FeedbackTable_Resolved_SNeII'
feedback_record_file_path  = '../Record__FB_Resolved_SNeII'


# define functions to get the feedback quantities
def get_particle_feedback_amount(ds, obj, ptype, table_file_prefix, quantity_name):

    SNIINxtE = obj[ptype, 'ParSNIINxtE'].d.astype('int')
    number_of_total_SNe    = SNIINxtE // 10000
    number_of_exploded_SNe = SNIINxtE %  10000

    returned_quantity_in_table_units = np.zeros_like(SNIINxtE)
    for nSNe in range(1, np.max(number_of_total_SNe) + 1):
        feedback_table = np.genfromtxt(table_file_prefix+'_N%06d'%nSNe,
                                       delimiter=None, comments='#',
                                       names=['explosion_time_in_Myr',
                                              'progenitor_mass_in_Msun',
                                              'ejected_gas_mass_in_Msun',
                                              'ejected_metal_mass_in_Msun',
                                              'ejected_energy_in_erg'],
                                       dtype=None, encoding=None)

        mask = (number_of_total_SNe == nSNe)
        if quantity_name == 'explosion_time_in_Myr':
            # mean of the explosion times of all the events,
            # including those have not happened yet
            returned_quantity_in_table_units[mask] = np.mean( feedback_table[quantity_name] )
        else:
            # sum of the ejected quantities of the happened events
            returned_quantity_in_table_units[mask] = np.insert(
                                                     np.cumsum( feedback_table[quantity_name] ),
                                                     0, 0 )[number_of_exploded_SNe[mask]]

    units = quantity_name.split('_')[-1]

    return ds.arr( returned_quantity_in_table_units, units )


def get_reconstructed_stellar_mass(ds, obj, ptype, table_file_prefix):
    ejected_mass = get_particle_feedback_amount(ds, obj, ptype,
                                                table_file_prefix,
                                                'ejected_gas_mass_in_Msun')

    return obj[ ptype, 'ParMass' ].in_units('Msun') + ejected_mass


def get_mean_SN_time(ds, obj, ptype, table_file_prefix):
    mean_SN_delay_time = get_particle_feedback_amount(ds, obj, ptype,
                                                      table_file_prefix,
                                                      'explosion_time_in_Myr')

    return obj[ ptype, 'ParCreTime' ].in_units('Myr') + mean_SN_delay_time


# load data
ds = yt.load( filein )


# define the particle filter for the newly formed stars
def new_star( pfilter, data ):
   filter = data[ 'all', 'ParCreTime' ] > 0
   return filter

add_particle_filter( 'new_star', function=new_star, filtered_type='all', requires=['ParCreTime'] )
ds.add_particle_filter( 'new_star' )


# get the explosion time
ad = ds.all_data()
star_expected_SN_time = get_mean_SN_time(ds, ad, 'new_star', feedback_table_file_prefix).in_units('Myr').d
star_masses = get_reconstructed_stellar_mass(ds, ad, 'new_star', feedback_table_file_prefix).in_units('Msun').d
star_count = len(star_expected_SN_time)

feedback_record = np.genfromtxt(feedback_record_file_path,
                                delimiter=None, comments='#',
                                names=['Rank', 'TID', 'lv', 'TimeOld', 'TimeNew', 'SNII_Time_first', 'SNII_Time_last', 'SNII_AccumNum', 'SNII_Energy', 'SNII_Mass', 'SNII_Metal', 'FB_Diameter', 'FB_Flu_Mass', 'Par_PUID', 'Par_NxtE', 'Par_Mass', 'Par_PosX', 'Par_PosY', 'Par_PosZ', 'Par_VelX', 'Par_VelY', 'Par_VelZ', 'Par_MetalMass', 'Flu_Dens', 'Flu_MomX', 'Flu_MomY', 'Flu_MomZ', 'Flu_Engy', 'Flu_Metal'],
                                dtype=None, encoding=None)
SN_explosion_time = ds.arr( feedback_record['SNII_Time_last'], 'code_time' ).in_units('Myr').d
SN_ones = np.ones_like( SN_explosion_time )
SN_count = len(SN_explosion_time)

print('Total number of stars = %d'%star_count)
print('Total number of SNe   = %d'%SN_count)

# bin the data
t_buffer = 40.0  # in Myr
t_start  = 0.0
t_end    = ds.current_time.in_units('Myr').d + t_buffer
t_bin    = np.linspace( start=t_start, stop=t_end, num=nbin+1 )
time     = 0.5*( t_bin[:-1] + t_bin[1:] )
star_upper_idx = np.digitize( star_expected_SN_time, bins=t_bin, right=True )
SN_upper_idx   = np.digitize( SN_explosion_time, bins=t_bin, right=True )
assert np.all( star_upper_idx > 0 ) and np.all( star_upper_idx < len(t_bin) ), 'incorrect star_upper_idx !!'
assert np.all( SN_upper_idx > 0 ) and np.all( SN_upper_idx < len(t_bin) ), 'incorrect SN_upper_idx !!'


# calculate the star formation rate and SNe number rate
Myr2yr = 1.0e6
sfr = np.array( [ star_masses[star_upper_idx == j+1].sum() /
                  ( (t_bin[j+1] - t_bin[j])*Myr2yr ) for j in range(len(time)) ] )
sfr[sfr == 0] = np.nan

SNrate = np.array( [ SN_ones[SN_upper_idx == j+1].sum() /
                     ( (t_bin[j+1] - t_bin[j])*Myr2yr ) for j in range(len(time)) ] )
SNrate[SNrate == 0] = np.nan


# plot
NumSNPerStarMass = ds.quan( ds.parameters['FB_ResolvedSNeII_NPerMass'], 'code_mass**-1' ).in_units('Msun**-1').v
plt.plot(time, sfr * NumSNPerStarMass, label='Mean-time expectation from SFR')
plt.plot(time, SNrate, label='Exploded SNe')
plt.xlim( 0.0, t_end-2.0*t_buffer )
plt.legend()
plt.xlabel(r'$\mathrm{t\ [Myr]}$', fontsize='large')
plt.ylabel(r'Event number rate $\mathrm{[yr^{-1}]}$', fontsize='large')


plt.savefig( fileout+'.png', bbox_inches='tight', pad_inches=0.05, dpi=dpi )
