import cluster_generator as cg
import unyt as u
from numpy.random import RandomState
import numpy as np

# Note that cluster_generator does not use unyt units for speed and simplicity,
# so mass units are Msun, length units are kpc, and time units are Myr

# Put the two clusters at a redshift z = 0.1
z = 0.1

# M200 for both clusters
M200_1 = 6.0e14 # in Msun
M200_2 = 2.0e14 # in Msun

conc = 4.0 # A good approximation to the concentration parameter for both clusters

# Find r200 for both clusters
r200_1 = cg.find_overdensity_radius(M200_1, 200.0, z=z)
r200_2 = cg.find_overdensity_radius(M200_2, 200.0, z=z)

# Scale radii to be used for the sNFW profiles
a1 = r200_1/conc
a2 = r200_2/conc

# For the total mass density profile, we will use a "super-NFW" profile, which
# is very similar to the NFW profile but falls off slightly faster (Lilley, E. J.,
# Wyn Evans, N., & Sanders, J.L. 2018, MNRAS)

# Determine the total mass for each sNFW profile
M1 = cg.snfw_total_mass(M200_1, r200_1, a1)
M2 = cg.snfw_total_mass(M200_2, r200_2, a2)

# Use this total mass to construct total mass profiles for each cluster
Mt1 = cg.snfw_mass_profile(M1, a1)
Mt2 = cg.snfw_mass_profile(M2, a2)

# Use the total mass profiles to determine r500/M500 and r2500/M2500 for
# each cluster
r500_1, M500_1 = cg.find_radius_mass(Mt1, z=z, delta=500.0)
r2500_1, M2500_1 = cg.find_radius_mass(Mt1, z=z, delta=2500.0)
r500_2, M500_2 = cg.find_radius_mass(Mt2, z=z, delta=500.0)
r2500_2, M2500_2 = cg.find_radius_mass(Mt2, z=z, delta=2500.0)

# Total mass density profiles for each cluster
rhot1 = cg.snfw_density_profile(M1, a1)
rhot2 = cg.snfw_density_profile(M2, a2)

# Sprinkle some stars in--2% of the total mass for each cluster
rhos1 = 0.02*rhot1
rhos2 = 0.02*rhot2

# Find the gas mass fraction within R500 (using the relationship between
# M500 and fgas from Vikhlinin, A., et al. 2009, ApJ, 692, 1033
f_g1 = cg.f_gas(M500_1)
f_g2 = cg.f_gas(M500_2)

# This sets the gas density profile using the functional form from Vikhlinin, A.,
# Kravtsov, A., Forman, W., et al. 2006, ApJ, 640, 691 for the first cluster. We
# set the scale density to 1.0 first and will rescale it in the next line by the
# gas mass within r500
rhog1 = cg.vikhlinin_density_profile(1.0, 0.2*r2500_1, 0.67*r200_1, 1.0, 0.67, 3.0)
rhog1 = cg.rescale_profile_by_mass(rhog1, f_g1*M500_1, r500_1)

# Same as above for the second cluster
rhog2 = cg.vikhlinin_density_profile(1.0, 0.2*r2500_2, 0.67*r200_2, 1.0, 0.67, 3.0)
rhog2 = cg.rescale_profile_by_mass(rhog2, f_g2*M500_2, r500_2)

# This is the plasma beta parameter for the ratio of the thermal pressure to the
# magnetic pressure
beta = 100.0

# This sets up the profiles for the first cluster assuming hydrostatic equilibrium,
# taking the gas density, total mass density, and stellar density as input
hse1 = cg.ClusterModel.from_dens_and_tden(0.1, 20000.0, rhog1, rhot1,
                                          stellar_density=rhos1)

# This sets a radial magnetic field strength profile using the beta parameter and
# the pressure in the profile, assuming p_B = B^2/s (thus gaussian=False)
hse1.set_magnetic_field_from_beta(beta, gaussian=False)

# These lines are the same as above for the second cluster
hse2 = cg.ClusterModel.from_dens_and_tden(0.1, 20000.0, rhog2, rhot2,
                                          stellar_density=rhos2)
hse2.set_magnetic_field_from_beta(beta, gaussian=False)

# Write the profiles for each cluster to files
hse1.write_model_to_h5("profile1.h5", overwrite=True)
hse2.write_model_to_h5("profile2.h5", overwrite=True)

# Set a random number generator for the generation of the magnetic field
# vector potential in 3D
prng = RandomState(24)

# This is the width of the GAMER simulation box and its center
w = 15000.0 # in kpc
center = np.array([0.5*w]*3)

# This determines the centers of the clusters, assuming a distance of
# 3 Mpc and zero impact parameter, centered on the box center
d = 3000.0 # in kpc
b = 0.0 # in kpc
center1, center2 = cg.compute_centers_for_binary(center, d, b)

# This sets up a 3D magnetic vector potential which GAMER will take the curl
# of on the AMR grid to get the initial B-field. It is a tangled field which
# uses a Kolmogorov spectrum with a large-scale cutoff of 500 kpc, a
# small-scale cutoff of 10 kpc, and is proportional on average to the pressure
# everywhere (given by the magnetic field profile of the clusters from above).
# Outside of r_max = 5000.0 kpc from each cluster center the average B-field
# is constant
left_edge = center-0.5*w
right_edge = center+0.5*w
dims = (256,)*3
bfield = cg.RadialRandomMagneticVectorPotential(left_edge, right_edge, dims,
                                                10.0, 500.0, center1,
                                                "profile1.h5", ctr2=center2,
                                                profile2="profile2.h5", r_max=5000.0)

# Write the 3D vector potential to the B_IC file
bfield.write_to_h5("B_IC", overwrite=True, length_unit="Mpc",
                   field_unit="sqrt(1e14*Msun/Mpc**3)*Mpc/(10*Gyr)")

# We now set up the velocities of the two clusters. Assume 1500 km/s
# relative velocity, and then use the M200 of the two clusters to
# set velocity vectors in roughly the CM frame. The velocity is in
# the x-direction only
velocity = (1500.0*u.km/u.s).to_value("kpc/Myr")
velocity1 = np.array([velocity*M200_2/(M200_1+M200_2), 0.0, 0.0])
velocity2 = np.array([-velocity*M200_1/(M200_1+M200_2), 0.0, 0.0])

# Now we set up the cluster initial conditions. use 2e6 DM particles,
# 4e4 star particles. At r_max = 5000.0 kpc, the profiles of each cluster
# are constant
num_particles = {"dm": 2_000_000, "star": 40_000}

ics = cg.ClusterICs("1to3_b0.0", 2, ["profile1.h5", "profile2.h5"],
                    [center1, center2], [velocity1, velocity2],
                    num_particles=num_particles, mag_file="B_IC", r_max=5000.0)

# This writes the GAMER-specific IC files that are needed, generates
# the particles, and prints out the contents of the Input__TestProb
# file which should be used
cg.setup_gamer_ics(ics)
