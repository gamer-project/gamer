#!/usr/bin/env python3.9
################################################################################

# Run Compute_profile.py to call these functions.
# Contains functions called by Compute_profiles.py

################################################################################

import argparse
import csv
import os
import sys

import numpy as np
import yt
from scipy.optimize import curve_fit, ridder
from VelocityDisp import *
from PowerSpec    import *

def soliton( x, xc ,time_a, particle_mass ):
    return 1.9/time_a*((particle_mass/1e-23)**-2)*((xc)**-4)/(1+9.1*1e-2*(x/xc)**2)**8*1e9

def find_virial_mass( r, mass_para, zeta, background_density ):
    mass = 10**np.interp(np.log10(r), np.log10( mass_para[0]), np.log10( mass_para[1] ) )
    return mass-zeta*background_density*(4/3*np.pi*r**3)

def compute_profile( ds, center, _halo_radius, halo_id, path_data ):

    if not os.path.exists(path_data):
        print( "ERROR: %s does not exist !!"%path_data )
        exit(0)

    print( "start computing profile...\n" )
    omega_M0             = ds.omega_matter
    newton_G             = ds.units.newtons_constant.to('(kpc*km**2)/(s**2*Msun)').d     #(kpc*km^2)/(s^2*Msun)
    background_density_0 = (1*ds.units.code_density ).to("Msun/kpc**3").d
    particle_mass        = (ds.parameters[ 'ELBDM_Mass' ]*ds.units.code_mass).to('eV/c**2').d
    zeta_0               = (18*np.pi**2 + 82*( omega_M0 - 1 ) - 39*( omega_M0 - 1 )**2 )/omega_M0

    # plot variable
    nbin             = 256
    max_radius       = 4e2
    coordinates = ds.arr( [center[0], center[1], center[2]], 'code_length' )

    # save file parameter
    halo_parameter_filename = "Halo_Parameter"

    # periodicity
    ds.force_periodicity()

    # find center
    # find the maximum value in a sphere extending halo_radius from center_guess
    halo_radius  = _halo_radius     # halo radius in cMpc/h --> change this value properly
    if ds.cosmological_simulation:
        center_guess       = coordinates
        sphere_guess       = ds.sphere( center_guess, (halo_radius, 'code_length') )
        center_find        = sphere_guess.quantities.max_location( 'density' )
        center_coordinate  = ds.arr( [center_find[1].d, center_find[2].d, center_find[3].d], 'code_length' )
    else:
        center_coordinate = coordinates
    print( "Center is ", center_coordinate.in_units('code_length') )

    # extract halo
    sp         = ds.sphere( center_coordinate, (0.30, 'code_length') ) # extract sphere with max 0.25 Mpc
    max_level  = ds.max_level
    min_radius = ds.domain_width.in_units("kpc")[0].d/2**max_level/ds.domain_dimensions[0]

    prof_mass_accumulate = yt.create_profile( sp, 'radius', fields = 'cell_mass',
                                              weight_field = None, n_bins = nbin ,
                                              units = {'radius': 'kpc', 'cell_mass': 'Msun'},
                                              extrema = {'radius': (min_radius,max_radius)}, accumulation = True )

    prof_dens            = yt.create_profile( sp, 'radius', fields = 'density',
                                              weight_field = 'cell_volume', n_bins = nbin ,
                                              units = {'radius': 'kpc','density': 'Msun/kpc**3'},
                                              extrema = {'radius': (min_radius,max_radius)} )

    prof_volume          = yt.create_profile( sp, 'radius', fields = 'cell_volume',
                                              weight_field = None, n_bins = nbin ,
                                              units = {'radius': 'kpc', 'cell_volume': 'kpc**3'},
                                              extrema = {'radius': (min_radius,max_radius)} )

    radius_o          = prof_dens.x.value
    density_o         = prof_dens['density'].value                       # density at radius
    mass_accumulate_o = prof_mass_accumulate['cell_mass'].value  # all mass within radius
    volume_o          = prof_volume['cell_volume'].value       # volume within radius

    # remove zero
    dens_idx = density_o != 0

    radius            = radius_o[dens_idx]
    density           = density_o[dens_idx]
    mass_accumulate   = mass_accumulate_o[dens_idx]
    circular_velocity = np.sqrt( newton_G * mass_accumulate/radius )

    # output profiles
    if not os.path.exists( path_data + '/prof_dens' ):
        os.makedirs( path_data + '/prof_dens' )

    if not os.path.exists( path_data + '/prof_mass' ):
        os.makedirs( path_data + '/prof_mass' )

    if not os.path.exists( path_data + '/prof_circular_vel' ):
        os.makedirs( path_data + '/prof_circular_vel' )

    with open( '%s/prof_dens/%s_%d_profile_data'%(path_data,ds,halo_id) , 'w' ) as file:
        writer = csv.writer( file, delimiter='\t' )
        writer.writerow( [f"{'#radius(ckpc)':<15}", f"{'density(Msun/ckpc**3)':<15}"] )
        for i in range( len(radius) ):
            writer.writerow( [f"{radius[i]:<15.8f}", f"{density[i]:<15.8f}"] )
#            writer.writerow( [radius[i], density[i]] )

    with open( '%s/prof_mass/%s_%d_mass_accumulate'%(path_data,ds,halo_id) , 'w' ) as file:
        writer = csv.writer( file, delimiter='\t' )
        writer.writerow( [f"{'#radius(ckpc)':<15}", f"{'mass(Msun)':<15}"] )
        for i in range( len(radius) ):
            writer.writerow( [f"{radius[i]:<15.8f}", f"{mass_accumulate[i]:<15.8f}"] )

    with open( '%s/prof_circular_vel/%s_%d_circular_velocity'%(path_data,ds,halo_id) , 'w' ) as file:
        writer = csv.writer( file, delimiter='\t' )
        writer.writerow( [f"{'#radius(ckpc)':<15}", f"{'Vcir(km/s)':<15}"] )
        for i in range( len(radius) ):
            writer.writerow( [f"{radius[i]:<15.8f}", f"{circular_velocity[i]:<15.8f}"] )

    ############# Output Halo's properties (core mass, halo mass, raidus) ####################

    # calculate virial mass to get halo radius
    if ds.cosmological_simulation:
        current_time_z = ds.current_redshift
    else:
        current_time_z = 0
    current_time_a = 1.0/(1+current_time_z)

    # defintion of zeta (halo radius)
    omega_M = (omega_M0*(1 + current_time_z)**3)/(omega_M0*(1 + current_time_z)**3 + (1 - omega_M0))
    zeta    = (18*np.pi**2 + 82*(omega_M - 1) - 39*(omega_M - 1)**2)/omega_M

    # use mass_accumulation directly
    halo_radius = ridder( lambda x:find_virial_mass(x,(radius, mass_accumulate),zeta, background_density_0),min_radius, max_radius )
    halo_mass   = 10**np.interp( np.log10( halo_radius ), np.log10( radius ), np.log10( mass_accumulate ) )

    #core radius 2 : xc = max/2
    try:
        core_radius_2 = ridder( lambda x: 10**np.interp( np.log10( x ), np.log10( radius ), np.log10( density )) - max(density)/2, radius[0], max( radius ) )
    except ValueError as e:
        print("Error while getting core_radius_2: likely unable to locate the center of halo. Please adjust center_first_guess or vicinity in Compute_profiles.py")
        exit(0)
        return
    core_mass_2   = 10**np.interp( np.log10( core_radius_2 ), np.log10( radius ), np.log10( mass_accumulate ) )

    #core radius 1 : curve fit
    avg           = (density > 0.1*max(density))
    popt, pcov    = curve_fit( lambda x, r_c:soliton(x, r_c, current_time_a, particle_mass), radius[avg], density[avg], bounds=(5e-1, 50) )
    core_radius_1 = popt[0]
    core_mass_1   = 10**np.interp( np.log10( core_radius_1 ),np.log10( radius ), np.log10( mass_accumulate ) )

    #core radius 3 : x = 0 solve equation
    a             = (2**(1.0/8) - 1)**(1.0/2)
    core_radius_3 = (max(density)/10**9/1.9*float(current_time_a)*(particle_mass/10**-23)**2)**-0.25
    core_mass_3   = ((4.2*10**9/((particle_mass/10**-23)**2*(float(core_radius_3*current_time_a)*10**3)))*(1/(a**2 + 1)**7)*(3465*a**13 + 23100*a**11 + 65373*a**9 + 101376*a**7 + 92323*a**5 + 48580*a**3 - 3465*a + 3465*(a**2 + 1)**7*np.arctan(a)))

    sto_list = []
    sto_list.append( int(str(ds).split("_")[-1]) )
    sto_list.append( halo_id )
    sto_list.append( particle_mass )
    sto_list.append( center_coordinate.in_units('code_length')[0].d )
    sto_list.append( center_coordinate.in_units('code_length')[1].d )
    sto_list.append( center_coordinate.in_units('code_length')[2].d )
    sto_list.append( current_time_z )
    sto_list.append( halo_radius*current_time_a )
    sto_list.append( halo_mass )
    sto_list.append( max(density)*current_time_a**-3 )
    sto_list.append( core_radius_1*current_time_a )
    sto_list.append( core_radius_2*current_time_a )
    sto_list.append( core_radius_3*current_time_a )
    sto_list.append( core_mass_1 )
    sto_list.append( core_mass_2 )
    sto_list.append( core_mass_3 )

    file_exists = os.path.exists( "%s/%s"%(path_data, halo_parameter_filename) )
    with open( "%s/%s"%(path_data,halo_parameter_filename), 'a+' ) as file:
        writer = csv.writer( file, delimiter='\t' )

        # write header
        if not file_exists:
            writer.writerow( [f"{'#snap':<6}",
                              f"{'halo_id':<8}",
                              f"{'m[eV]':<8}",
                              f"{'x[cMpc/h]':<11}",
                              f"{'y[cMpc/h]':<11}",
                              f"{'z[cMpc/h]':<11}",
                              f"{'redshift':<11}",
                              f"{'r_vir[kpc]':<13}",
                              f"{'Mvir[Msub]':<13}",
                              f"{'rho_max[Msun/kpc**3]':<20}",
                              f"{'r_c1[kpc]':<13}",
                              f"{'r_c2[kpc]':<13}",
                              f"{'r_c3[kpc]':<13}",
                              f"{'M_c1[Msun]':<13}",
                              f"{'M_c2[Msun]':<13}",
                              f"{'M_c3[Msun]':<13}"
                             ] )
        # write data
        writer.writerow( [
        f"{sto_list[0]:<6}",
        f"{sto_list[1]:<8}",
        f"{sto_list[2]:<8.1e}",
        f"{sto_list[3]:<11.8f}",
        f"{sto_list[4]:<11.8f}",
        f"{sto_list[5]:<11.8f}",
        f"{sto_list[6]:<11.8f}",
        f"{sto_list[7]:<13.8f}",
        f"{sto_list[8]:<13.8e}",
        f"{sto_list[9]:<20.8e}",
        f"{sto_list[10]:<13.8f}",
        f"{sto_list[11]:<13.8f}",
        f"{sto_list[12]:<13.8f}",
        f"{sto_list[13]:<13.8e}",
        f"{sto_list[14]:<13.8e}",
        f"{sto_list[15]:<13.8e}"
        ] )

    # Compute Velocity Dispersion & Power spectrum seperately
    Compute_VelocityDispersion( ds, center_coordinate.in_units('code_length').d, halo_id, path_data )
    Compute_PowerSpectrum( ds, center_coordinate.in_units('code_length').d, core_radius_1, halo_id, path_data )
