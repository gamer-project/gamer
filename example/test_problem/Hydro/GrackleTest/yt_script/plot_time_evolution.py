import argparse
import sys
import yt
import matplotlib.pyplot as plt
import numpy as np

# load the command-line parameters
parser = argparse.ArgumentParser( description='Plot the gas time evolution' )

parser.add_argument( '-p', action='store', required=False, type=str, dest='prefix',
                     help='path prefix [%(default)s]', default='../' )
parser.add_argument( '-s', action='store', required=True,  type=int, dest='idx_start',
                     help='first data index' )
parser.add_argument( '-e', action='store', required=True,  type=int, dest='idx_end',
                     help='last data index' )
parser.add_argument( '-d', action='store', required=False, type=int, dest='didx',
                     help='delta data index [%(default)d]', default=1 )

args=parser.parse_args()

# take note
print( '\nCommand-line arguments:' )
print( '-------------------------------------------------------------------' )
print( ' '.join(map(str, sys.argv)) )
print( '-------------------------------------------------------------------\n' )


idx_start = args.idx_start
idx_end   = args.idx_end
didx      = args.didx
prefix    = args.prefix


yt.enable_parallelism()

ts = yt.DatasetSeries( [ prefix+'/Data_%06d'%idx for idx in range(idx_start, idx_end+1, didx) ] )


time               = np.array([])
density            = np.array([])
metal_density      = np.array([])
electron_density   = np.array([])
HI_density         = np.array([])
HII_density        = np.array([])
HeI_density        = np.array([])
HeII_density       = np.array([])
HeIII_density      = np.array([])
HM_density         = np.array([])
H2I_density        = np.array([])
H2II_density       = np.array([])
DI_density         = np.array([])
DII_density        = np.array([])
HDI_density        = np.array([])
temperature        = np.array([])
gr_temperature     = np.array([])
gr_mu              = np.array([])
temperature_max    = np.array([])
temperature_min    = np.array([])
gr_temperature_max = np.array([])
gr_temperature_min = np.array([])


# collect the results
for ds in ts.piter():

    time                 = np.append( time              , ds.current_time.in_units('Myr').d )
    density              = np.append( density           , np.average(ds.all_data()['density'            ].in_units('g/cm**3').d) )

    if ds.parameters['Grackle_Metal'] == 1:
        metal_density    = np.append( metal_density     , np.average(ds.all_data()['metal_density'      ].in_units('g/cm**3').d) )

    if ds.parameters['Grackle_Primordial'] >= 1:
        electron_density = np.append( electron_density  , np.average(ds.all_data()['electron_density'   ].in_units('g/cm**3').d) )
        HI_density       = np.append( HI_density        , np.average(ds.all_data()['HI_density'         ].in_units('g/cm**3').d) )
        HII_density      = np.append( HII_density       , np.average(ds.all_data()['HII_density'        ].in_units('g/cm**3').d) )
        HeI_density      = np.append( HeI_density       , np.average(ds.all_data()['HeI_density'        ].in_units('g/cm**3').d) )
        HeII_density     = np.append( HeII_density      , np.average(ds.all_data()['HeII_density'       ].in_units('g/cm**3').d) )
        HeIII_density    = np.append( HeIII_density     , np.average(ds.all_data()['HeIII_density'      ].in_units('g/cm**3').d) )

    if ds.parameters['Grackle_Primordial'] >= 2:
        HM_density       = np.append( HM_density        , np.average(ds.all_data()['HM_density'         ].in_units('g/cm**3').d) )
        H2I_density      = np.append( H2I_density       , np.average(ds.all_data()['H2I_density'        ].in_units('g/cm**3').d) )
        H2II_density     = np.append( H2II_density      , np.average(ds.all_data()['H2II_density'       ].in_units('g/cm**3').d) )

    if ds.parameters['Grackle_Primordial'] >= 3:
        DI_density       = np.append( DI_density        , np.average(ds.all_data()['DI_density'         ].in_units('g/cm**3').d) )
        DII_density      = np.append( DII_density       , np.average(ds.all_data()['DII_density'        ].in_units('g/cm**3').d) )
        HDI_density      = np.append( HDI_density       , np.average(ds.all_data()['HDI_density'        ].in_units('g/cm**3').d) )

    temperature          = np.append( temperature       , np.average(ds.all_data()['temperature'        ].in_units('K').d)       )
    gr_temperature       = np.append( gr_temperature    , np.average(ds.all_data()['grackle_temperature'].in_units('K').d)       )
    gr_mu                = np.append( gr_mu             , np.average(ds.all_data()['grackle_mu'         ].d)                     )

    temperature_max      = np.append( temperature_max   , np.max(    ds.all_data()['temperature'        ].in_units('K').d)       )
    temperature_min      = np.append( temperature_min   , np.min(    ds.all_data()['temperature'        ].in_units('K').d)       )
    gr_temperature_max   = np.append( gr_temperature_max, np.max(    ds.all_data()['grackle_temperature'].in_units('K').d)       )
    gr_temperature_min   = np.append( gr_temperature_min, np.min(    ds.all_data()['grackle_temperature'].in_units('K').d)       )


# plot the results
f, ax = plt.subplots( 1, 2, figsize=(12.8, 4.8) )


# plot the temperature evolution
ax[0].plot( time,  gr_temperature,       'k-s', label='Grackle T'    )
ax[0].plot( time,  temperature,          'k--', label='GAMER T'      )
ax[0].plot( time,  gr_temperature/gr_mu, 'k:',  label='Grackle T/mu' )

if ds.parameters['GrackleTest_DefaultTestMode'] == 3:
    # plot the expected temperature evolution for the given heating and cooling rates
    r_heat =  ds.parameters['GrackleTest_HeatingRate'] * ds.parameters['Grackle_HydrogenMFrac'] * (ds.parameters['Gamma']-1) * ds.parameters['MolecularWeight'] / ds.units.boltzmann_constant.in_units('erg/K').v *(1.*ds.units.Myr).in_units('s').v
    r_cool = (density[0]*ds.parameters['Grackle_HydrogenMFrac']/ds.units.hydrogen_mass.in_units('g').v*ds.parameters['GrackleTest_CoolingRate']) * ds.parameters['Grackle_HydrogenMFrac'] * (ds.parameters['Gamma']-1) * ds.parameters['MolecularWeight'] / ds.units.boltzmann_constant.in_units('erg/K').v *(1.*ds.units.Myr).in_units('s').v
    ax[0].plot( time,  temperature_max[0]+time*r_heat,  'g-.', label='Analytical max(T)' )
    ax[0].plot( time,  temperature_min[0]-time*r_cool,  'y-.', label='Analytical min(T)' )

    # plot the maximum and minimum temperature evolution as the results of the given heating and cooling rates
    ax[0].plot( time,  temperature_max, 'r--', label='max(T)' )
    ax[0].plot( time,  temperature_min, 'b--', label='min(T)' )

ax[0].set_yscale( 'log' )
ax[0].set_xlim( 0.0, 1.33*np.max(time) )
ax[0].set_ylim( 0.1*np.min(gr_temperature_min), 10.0*np.max(gr_temperature_max) )
ax[0].set_xlabel( 'Time (Myr)',      fontsize='large' )
ax[0].set_ylabel( 'Temperature (K)', fontsize='large' )
ax[0].legend()


# plot the density evolution
ax[1].plot(     time,  density,           'k-s', label=''      )

if ds.parameters['Grackle_Metal'] == 1:
    ax[1].plot( time,  metal_density,     'g-',  label='Metal' )

if ds.parameters['Grackle_Primordial'] >= 1:
    ax[1].plot( time,  electron_density,  'g--', label='e'     )
    ax[1].plot( time,  HI_density,        'r-',  label='HI'    )
    ax[1].plot( time,  HII_density,       'r--', label='HII'   )
    ax[1].plot( time,  HeI_density,       'b-',  label='HeI'   )
    ax[1].plot( time,  HeII_density,      'b--', label='HeII'  )
    ax[1].plot( time,  HeIII_density,     'b-.', label='HEIII' )

if ds.parameters['Grackle_Primordial'] >= 2:
    ax[1].plot( time,  HM_density,        'r:',  label='HM'    )
    ax[1].plot( time,  H2I_density,       'm--', label='H2I'   )
    ax[1].plot( time,  H2II_density,      'm-',  label='H2II'  )

if ds.parameters['Grackle_Primordial'] >= 3:
    ax[1].plot( time,  DI_density,        'c-',  label='DI'    )
    ax[1].plot( time,  DII_density,       'c--', label='DII'   )
    ax[1].plot( time,  HDI_density,       'y-',  label='HDI'   )

ax[1].set_yscale( 'log' )
ax[1].set_xlim( 0.0, 1.33*np.max(time) )
ax[1].set_ylim( 3.0e-6*np.min(density), 3.0*np.max(density) )
ax[1].set_xlabel( 'Time (Myr)',         fontsize='large' )
ax[1].set_ylabel( 'Density (g/cm$^3$)', fontsize='large' )
ax[1].legend()


f.subplots_adjust( wspace=0.2 )
f.savefig( 'fig_time_evolution.png', bbox_inches='tight', pad_inches=0.05, dpi=150 )
