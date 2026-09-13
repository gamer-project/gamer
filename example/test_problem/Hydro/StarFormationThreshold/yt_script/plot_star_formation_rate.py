# ref: https://yt-project.org/docs/dev/cookbook/calculating_information.html#using-particle-filters-to-calculate-star-formation-rates

import yt
import numpy as np
from yt.data_objects.particle_filters import add_particle_filter
from matplotlib import pyplot as plt
import warnings


filein  = "../Data_%06d"%np.genfromtxt("../Record__Dump")[-1][0]
nbin    = int(np.genfromtxt("../Record__Dump")[-1][0]) -2
dpi     = 150


# load data
ds = yt.load( filein )

# reload data for comoving simulation
if ds.parameters['Comoving']:
    import ComovingGAMER_RevisedDataset
    ds = yt.DatasetSeries(ComovingGAMER_RevisedDataset.RevisedDatasets([ filein ]))[0]
    ComovingGAMER_RevisedDataset.check_comoving_units( ds )


# define the derived fields for the star formation rate

# gas density divided by code_density
def _density_code_density( field, data ):
   return data["density"] / data.ds.units.code_density

ds.add_field( ("gas", "density_code_density"), function=_density_code_density, sampling_type="cell", units="dimensionless" )

# jeans length divided by cell size
def _jeans_length_dh( field, data ):
   return np.sqrt( (np.pi * data["sound_speed"]**2 )/( data.ds.units.newtons_constant * data["density"] ) )/data["dx"]

ds.add_field( ("gas", "jeans_length_dh"), function=_jeans_length_dh, sampling_type="cell", units="dimensionless" )

# star mass formed per unit time (in code_mass/code_time)
def _SchmidtLaw_star_formation_rate( field, data ):
   t_ff       = np.sqrt( ( 3.0 * np.pi )/( 32.0 * data.ds.units.newtons_constant * data["density"] ) )
   efficiency = data.ds.parameters['SF_CreateStar_MassEff']
   return efficiency * data["cell_mass"] / t_ff

ds.add_field( ("gas", "SchmidtLaw_star_formation_rate"), function=_SchmidtLaw_star_formation_rate, sampling_type="cell", units="g/s" )


# define the particle filter for the newly formed stars
def new_star( pfilter, data ):
   filter = data[ "all", "ParCreTime" ] > 0
   return filter

add_particle_filter( "new_star", function=new_star, filtered_type="all", requires=["ParCreTime"] )
ds.add_particle_filter( "new_star" )


# get the mass and creation time of the new stars
ad            = ds.all_data()
mass          = ad[ "new_star", "ParMass" ].in_units( "Msun" )
if ds.parameters['Comoving']:
    from yt.utilities.cosmology import Cosmology
    co = Cosmology(hubble_constant=ds.parameters['Hubble0'], omega_matter=ds.parameters['OmegaM0'], omega_lambda=1.0-ds.parameters['OmegaM0'])
    # convert from scale factor to physical time
    creation_time = co.t_from_a( ad[ "new_star", "ParCreTime" ].in_units( "code_time" ).value ).in_units( "Myr" )
else:
    creation_time = ad[ "new_star", "ParCreTime" ].in_units( "Myr" )


# bin the data
if ds.parameters['Comoving']:
    z_start   = 1.0/ds.parameters['A_Init'] - 1.0
    z_end     = ds.current_redshift
    z_bin     = np.linspace( start=z_start, stop=z_end, num=nbin+1 )
    redshift  = 0.5*( z_bin[:-1] + z_bin[1:] )
    t_bin     = np.array( [ co.t_from_z(z).in_units( "Myr" ).value for z in z_bin ] )
else:
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


# calculate the accumulated star mass
starmass = np.array(  [ mass[upper_idx <= j+1].sum() for j in range(len(time)) ]  )
starmass[starmass == 0] = np.nan


# calculate the analytical star formation rate
if   ds.parameters['SF_CreateStar_Scheme'] == 1:   # with minimum density threshold
   a_now = ds.scale_factor if ds.parameters['Comoving'] else 1.0
   star_formation_region = ds.cut_region( ds.all_data(), ["obj['gas', 'density_code_density'] >= %.2e"%(ds.parameters['SF_CreateStar_MinGasDens']*a_now**3)] )

elif ds.parameters['SF_CreateStar_Scheme'] == 2:   # with maximum Jeans length thredshold
   star_formation_region = ds.cut_region( ds.all_data(), ["obj['gas', 'jeans_length_dh'] <= %.2e"%ds.parameters['SF_CreateStar_MaxGasJeansL']] )

else:
   raise RuntimeError('Unsupported SF_CreateStar_Scheme = %d !!'%(ds.parameters['SF_CreateStar_Scheme']))

sfr_analytical = star_formation_region.quantities.total_quantity('SchmidtLaw_star_formation_rate').in_units('Msun/yr').d


# The analytical star formation rate is calculated only by the last data and assume the density is a constant
if ds.parameters['Opt__ResetFluid'] == 0:
   warnings.warn('WARNING: The analytical solution assume a constant-density environment !!')


# plot the star formation rate
if ds.parameters['Comoving']:
    # assume the SFR scales as a^3 due to the increasing box volume
    plt.plot( redshift, sfr_analytical*((1+ds.current_redshift)/(1+redshift))**3, color='r', linestyle='--', label='Analytical' )
    plt.plot( redshift, sfr, label='Simulation' )
    plt.xlabel( "$z$", fontsize="large" )
    plt.xlim( z_end, z_start )
else:
    plt.axhline( sfr_analytical, color='r', linestyle='--', label='Analytical' )
    plt.plot( time, sfr, label='Simulation' )
    plt.ylim( 0.0, 500 )
    plt.xlabel( "$\mathrm{t\ [Myr]}$", fontsize="large" )

plt.legend()
plt.ylabel( "$\mathrm{SFR\ [M_\odot yr^{-1}]}$", fontsize="large" )
plt.savefig( "fig__star_formation_rate.png", bbox_inches="tight", pad_inches=0.05, dpi=dpi )
plt.close()


# plot the accumulated stellar mass
if ds.parameters['Comoving']:
    def IndefIntegral(a):
        """
        \int a^3 dt = \int a^3 dt/da da = \int a^3 / (a*H(a)) da = \int a^2 / H(a) da,
        assuming H(a) = H0 * sqrt( Omega_M0*a^-3 + Omega_Lambda0 )
        """
        H0            = ds.quan( ds.parameters['Hubble0']*100.0, 'km/s/Mpc' ).in_units('yr**-1').v
        Omega_M0      = ds.parameters['OmegaM0']
        Omega_Lambda0 = 1.0 - Omega_M0
        return 1.0/(3.0 * H0 * Omega_Lambda0**1.5) * ( np.sqrt( Omega_Lambda0 * a**3 * (Omega_M0 + Omega_Lambda0 * a**3) ) - Omega_M0 * np.arcsinh( np.sqrt( Omega_Lambda0 * a**3 / Omega_M0 ) ) )

    # in the default setup, the cell size exceeds the Jeans length after a=0.224
    a_start  = 0.224 if ds.parameters['SF_CreateStar_Scheme'] == 2 else ds.parameters['A_Init']
    a_end    = ds.scale_factor
    Integral = IndefIntegral(a_end) - IndefIntegral(a_start)

    # assume the SFR scales as a^3; calculate \int_{a_start}^{a_end} SFR(a) dt
    plt.axhline( sfr_analytical*(1+ds.current_redshift)**3*Integral, color='r', linestyle='--', label='Expected Final' )
    plt.plot( redshift, starmass, label='Simulation' )
    plt.xlabel( "$z$", fontsize="large" )
    plt.xlim( z_end, z_start )
else:
    # assume the SFR is a constant
    plt.axhline( sfr_analytical*ds.current_time.in_units("yr").v, color='r', linestyle='--', label='Expected Final' )
    plt.plot( time, starmass, label='Simulation' )
    plt.xlabel( "$\mathrm{t\ [Myr]}$", fontsize="large" )

plt.legend()
plt.ylabel( "$M_*$ [$M_\odot$]", fontsize="large" )
plt.yscale("log")
plt.savefig( "fig__stellar_mass.png", bbox_inches="tight", pad_inches=0.05, dpi=dpi )
plt.close()
