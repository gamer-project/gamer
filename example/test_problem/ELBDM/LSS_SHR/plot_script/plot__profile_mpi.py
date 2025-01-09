#!/usr/bin/env python3

import argparse
import csv
import os
import sys

import fit_NFW
import numpy as np
import pandas as pd
import yt
from scipy.optimize import curve_fit, ridder

# load the command-line parameters
parser = argparse.ArgumentParser( description='Plot profile and out put Halo_parameter' )

parser.add_argument( '-p', action='store', required=False, type=str, dest='prefix',
                     help='prefix [%(default)s]', default='../' )
parser.add_argument( '-s', action='store', required=True,  type=int, dest='idx_start',
                     help='first data index' )
parser.add_argument( '-e', action='store', required=True,  type=int, dest='idx_end',
                     help='last data index' )
parser.add_argument( '-d', action='store', required=False, type=int, dest='didx',
                     help='delta data index [%(default)d]', default=1 )
parser.add_argument( '-halo', action='store', required=False, type=int, dest='halo',
                     help='which halo [%(default)d]', default=1 )

args=parser.parse_args()

# take note
print( '\nCommand-line arguments:' )
print( '-------------------------------------------------------------------' )
for t in range( len(sys.argv) ):
   print(sys.argv[t], end = ' '),
print( '' )
print( '-------------------------------------------------------------------' )


idx_start   = args.idx_start
idx_end     = args.idx_end
didx        = args.didx
prefix      = args.prefix
halo        = args.halo

yt.enable_parallelism()
ts = yt.DatasetSeries( [ prefix+'/Data_%06d'%idx for idx in range(idx_start, idx_end+1, didx) ] )


# constants
omega_M0             = ts[0].omega_matter
background_density_0 = (1*ts[0].units.code_density).to("Msun/kpc**3").d
particle_mass        = (ts[0].parameters['ELBDM_Mass']*ts[0].units.code_mass).to('eV/c**2').d
zeta_0               = (18*np.pi**2 + 82*(omega_M0 - 1) - 39*(omega_M0 - 1)**2)/omega_M0 


# plot variable
nbin             = 256
max_radius       = 4e2

# save file parameter
halo_parameter_filename = "Halo_Parameter_%d"%halo
writing_mode = 'append' if os.path.exists(halo_parameter_filename) else 'new'
# storage dictionary
storage = {}

# new : change x, y, z to the guessed halo center
# append : according to the last line of halo_parameter_filename
if writing_mode == 'new':
    coordinates = ts[0].arr( [ 0, 2.2, 0.4 ] , 'code_length' ) 
elif writing_mode == 'append':
    df = pd.read_csv( halo_parameter_filename, sep = r'\s+' , header = 0 )
    coordinates = ts[0].arr( [ df['x'].iloc[-1], df['y'].iloc[-1], df['z'].iloc[-1] ] , 'code_length' ) 
    print('%.21f'%coordinates.d[0], '%.21f'%coordinates.d[1], '%.21f'%coordinates.d[2])
else:
    print( 'writing_mode error' )
    sys.exit()

def soliton(x, xc ,time_a):
    return 1.9/time_a*((particle_mass/1e-23)**-2)*((xc)**-4)/(1+9.1*1e-2*(x/xc)**2)**8*1e9

def find_virial(para,zeta,current_time_a, NFW = False):
    # in comoving coordinate
    def enclosed_mass(x):
        if NFW:
            x = x/current_time_a
            rho0 = para[0]
            Rs = para[1]
            mass = 4*np.pi*rho0*Rs**3*(np.log((Rs+x)/Rs)-(x/(Rs+x)))
        else:
            r = para[0]
            mass = para[1]
            mass = 10**np.interp(np.log10(x), np.log10(r), np.log10(mass))
        return mass
    
    def zeta_mass(x):
        return zeta*background_density_0*(4/3*np.pi*x**3)
    
    halo_radius = ridder(lambda x:enclosed_mass(x) - zeta_mass(x), min_radius, max_radius)
    
    return halo_radius

def get_zeta(current_time_z):
    # defintion of zeta (halo radius)
    omega_M = (omega_M0*(1 + current_time_z)**3)/(omega_M0*(1 + current_time_z)**3 + (1 - omega_M0))
    zeta = (18*np.pi**2 + 82*(omega_M - 1) - 39*(omega_M - 1)**2)/omega_M
    return zeta

for sto, ds in ts.piter(storage=storage):
    # periodicity
    ds.force_periodicity()

    # find center 
    # find the maximum value in a sphere extending halo_radius from center_guess
    halo_radius  = 0.1            # halo radius in Mpc/h --> change this value properly
    center_guess = coordinates
    sphere_guess = ds.sphere( center_guess, (halo_radius, 'code_length' ) )
    center_find  = sphere_guess.quantities.max_location( 'density' )
    center_coordinate  = ds.arr( [center_find[1].d, center_find[2].d, center_find[3].d], 'code_length' )

    print(center_coordinate.in_units('code_length')) 

    # extract halo
    sp  = ds.sphere( center_coordinate, (max_radius, 'kpc') )
    max_level = ds.max_level
    min_radius = ds.domain_width.in_units("kpc")[0].d/2**max_level/ds.domain_dimensions[0]
    print('resolution',min_radius)

    prof_mass_accumulate = yt.create_profile( sp, 'radius', fields = 'cell_mass',
                                              weight_field = None, n_bins = nbin ,
                                              units = {'radius': 'kpc', 'cell_mass': 'Msun'}, 
                                              extrema = {'radius': (min_radius,max_radius)}, accumulation = True)

    prof_dens            = yt.create_profile( sp, 'radius', fields = 'density',
                                              weight_field = 'cell_volume', n_bins = nbin ,
                                              units = {'radius': 'kpc','density': 'Msun/kpc**3'}, 
                                              extrema = {'radius': (min_radius,max_radius)})

    prof_volume          = yt.create_profile( sp, 'radius', fields = 'cell_volume',
                                              weight_field = None, n_bins = nbin ,
                                              units = {'radius': 'kpc', 'cell_volume': 'kpc**3'}, 
                                              extrema = {'radius': (min_radius,max_radius)})
    
    radius_o = prof_dens.x.value
    density_o = prof_dens['density'].value                       # density at radius 
    mass_accumulate_o = prof_mass_accumulate['cell_mass'].value  # all mass within radius
    volume_o = prof_volume['cell_volume'].value       # volume within radius

    # remove zero
    radius = []
    density = []
    mass_accumulate = []
    volume = []

    for i in range(len(radius_o)):
        if(density_o[i]!=0):
            radius.append(radius_o[i])
            density.append(density_o[i])
            mass_accumulate.append(mass_accumulate_o[i])
            volume.append(volume_o[i])
            
    radius = np.array(radius)
    density = np.array(density)
    mass_accumulate = np.array(mass_accumulate)
    
    # calculate virial mass to know halo radius
    if ds.cosmological_simulation:
        current_time_z = ds.current_redshift
    else:
        current_time_z = 0
    current_time_a = 1.0/(1+current_time_z)

    zeta = get_zeta(current_time_z)

    # use mass_accumulation directly
    # halo_radius = ridder(lambda x:find_virial_mass(x,(radius, mass_accumulate),zeta),min_radius, max_radius)
    halo_radius = find_virial((radius, mass_accumulate), zeta, current_time_a)
    halo_mass = 10**np.interp(np.log10(halo_radius), np.log10(radius), np.log10(mass_accumulate))
        
    #core radius 1 : curve fit
    #core radius 2 : xc = max/2
    #core radius 3 : x = 0 solve equation 

    #core radius 2
    core_radius_2 = ridder(lambda x: 10**np.interp(np.log10(x),np.log10(radius), np.log10(density)) - max(density)/2, radius[0], max(radius))
    core_mass_2 = 10**np.interp(np.log10(core_radius_2),np.log10(radius), np.log10(mass_accumulate))   
    
    #core radius 1
    avg = (density > 0.1*max(density))
    popt, pcov = curve_fit(lambda x, r_c:soliton(x, r_c, current_time_a), radius[avg], density[avg],bounds=(0, np.inf))
    core_radius_1 = popt[0]
    # a = (2**(1.0/8) - 1)**(1.0/2)
    # core_mass_1 = ((4.2*10**9/((particle_mass/10**-23)**2*((core_radius_1*current_time_a)*10**3)))*(1/(a**2 + 1)**7)*(3465*a**13 + 23100*a**11 + 65373*a**9 + 101376*a**7 + 92323*a**5 + 48580*a**3 - 3465*a + 3465*(a**2 + 1)**7*np.arctan(a)))
    core_mass_1 = 10**np.interp(np.log10(core_radius_1),np.log10(radius), np.log10(mass_accumulate))   

    #core radius 3
    core_radius_3 = (max(density)/10**9/1.9*current_time_a*(particle_mass/10**-23)**2)**-0.25
    a = (2**(1.0/8) - 1)**(1.0/2)	
    core_mass_3 = ((4.2*10**9/((particle_mass/10**-23)**2*((core_radius_3*current_time_a)*10**3)))*(1/(a**2 + 1)**7)*(3465*a**13 + 23100*a**11 + 65373*a**9 + 101376*a**7 + 92323*a**5 + 48580*a**3 - 3465*a + 3465*(a**2 + 1)**7*np.arctan(a)))


    ### NFW
    range_FDM = (radius>core_radius_1*6) & (radius<halo_radius)
    rho0, Rs = fit_NFW.fit_NFW_dens(radius*current_time_a, density/current_time_a**3, range_FDM)
    halo_radius_NFW = find_virial((rho0,Rs), zeta, current_time_a, True)
    halo_mass_NFW = fit_NFW.NFW_mass((rho0, Rs), halo_radius_NFW*current_time_a)
    Ep_fit = fit_NFW.get_Ep_ideal(halo_mass_NFW, halo_radius_NFW, halo_radius_NFW/Rs, 'NFW')

    Ep_sim = fit_NFW.get_Ep_sph_sym((radius,density), halo_radius, halo_mass)

    sto_list = []
    sto_list.append(int(str(ds).split("_")[-1]))
    sto_list.append(particle_mass)
    sto_list.append(center_coordinate.in_units('code_length')[0].d)
    sto_list.append(center_coordinate.in_units('code_length')[1].d)
    sto_list.append(center_coordinate.in_units('code_length')[2].d)
    sto_list.append(current_time_z)
    sto_list.append(halo_radius*current_time_a)
    sto_list.append(halo_mass)
    sto_list.append(max(density)*current_time_a**-3)
    sto_list.append(core_radius_1*current_time_a)
    sto_list.append(core_radius_2*current_time_a)
    sto_list.append(core_radius_3*current_time_a)
    sto_list.append(core_mass_1)
    sto_list.append(core_mass_2)
    sto_list.append(core_mass_3)
    sto_list.append(rho0)
    sto_list.append(Rs)
    sto_list.append(halo_radius_NFW*current_time_a)
    sto_list.append(halo_mass_NFW)
    sto_list.append(Ep_fit)
    sto_list.append(Ep_sim)


    sto.result_id = str(ds)
    sto.result    = sto_list
    
    with open('prof_dens/%s_%d_profile_data'%(ds,halo) , 'w') as file:
        writer = csv.writer(file, delimiter='\t')
        writer.writerow(['radius(kpccm)', 'density(Msun/kpccm**3)' ])
        for i in range(len(radius)):
            writer.writerow([radius[i], density[i]])

    with open('prof_mass/%s_%d_mass_accumulate'%(ds,halo) , 'w') as file:
        writer = csv.writer(file, delimiter='\t')
        writer.writerow(['radius(kpccm)', 'mass(Msun)' ])
        for i in range(len(radius)):
            writer.writerow([radius[i], mass_accumulate[i]])

    with open('volume/%s_%d_volume'%(ds,halo) , 'w') as file:
        writer = csv.writer(file, delimiter='\t')
        writer.writerow(['radius(kpccm)', 'volume(kpccm**3)' ])
        for i in range(len(radius)):
            writer.writerow([radius[i], volume[i]])

if yt.is_root():
   if writing_mode != 'append':
      with open("%s"%halo_parameter_filename , 'w') as f:
         f.write(" #    mass                        x                        y                        z          time   halo_radius     halo_mass  peak_density core_radius_1 core_radius_2 core_radius_3   core_mass_1   core_mass_2   core_mass_3          rho0            Rs    NFW_radius      NFW_mass        fit_Ep        sim_Ep\n")
   with open("%s"%halo_parameter_filename , 'a') as f:
      writer = csv.writer(f, lineterminator='\n')
      for L in sorted(storage.items()):
        sto_list = []
        for i in range(len(L[1])):
            if i == 0:
                sto_list.append(str(L[1][i]))
            elif i == 1:
                sto_list.append(str("{:.1e}".format(L[1][i])))
            elif i <= 4:
                sto_list.append(str("{:.22f}".format(L[1][i])))
            else:
                sto_list.append(str("{:.7e}".format(L[1][i])))

        writer.writerow([' '.join(sto_list)])

