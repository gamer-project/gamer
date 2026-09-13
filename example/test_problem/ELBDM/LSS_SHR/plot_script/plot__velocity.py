#!/usr/bin/env python3

import argparse
import csv
import os
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import yt
from scipy import constants
from scipy.integrate import quad

# load the command-line parameters
parser = argparse.ArgumentParser( description='Profile' )

parser.add_argument( '-p', action='store', required=False, type=str, dest='prefix',
                     help='path prefix [%(default)s]', default='' )
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
   print (str(sys.argv[t]),end=' ')
#   print str(sys.argv[t]), 
print( '' )
print( '-------------------------------------------------------------------\n' )

idx_start   = args.idx_start
idx_end     = args.idx_end
didx        = args.didx
prefix      = args.prefix
halo        = args.halo

ds = yt.load('../Data_%06d'%idx_start)
### constant
kpc2km       = (1*yt.units.kpc).to('km').d
h            = ds.hubble_constant #dimensionless Hubble parameter 0.1km/(s*kpc)
hubble0      = h*0.1 # km/(s*kpc)
omega_M0     = ds.omega_matter
omega_lambda = ds.omega_lambda
newton_G     = ds.units.newtons_constant.to('(kpc**3)/(s**2*Msun)').d  #(kpc^3)/(s^2*Msun)


### read data
df_FDM = pd.read_csv( 'Halo_Parameter_%d'%halo , sep = r'\s+' , header = 0 , index_col='#')

FDM_v_path = ''

core_Both    = np.zeros(idx_end+1-idx_start)
core_QP      = np.zeros(idx_end+1-idx_start)
core_Bulk    = np.zeros(idx_end+1-idx_start)
inner_Both   = np.zeros(idx_end+1-idx_start)
inner_QP     = np.zeros(idx_end+1-idx_start)
inner_Bulk   = np.zeros(idx_end+1-idx_start)
halo_Both    = np.zeros(idx_end+1-idx_start)
halo_QP      = np.zeros(idx_end+1-idx_start)
halo_Bulk    = np.zeros(idx_end+1-idx_start)
average_Both = np.zeros(idx_end+1-idx_start)
average_QP   = np.zeros(idx_end+1-idx_start)
average_Bulk = np.zeros(idx_end+1-idx_start)

def shell_mass(r,dens):
    return 4*np.pi*r**2*dens

def CDM_dens(x, dens_parameter):
    r_CDM = dens_parameter[0]
    dens_CDM = dens_parameter[1]
    return 10**np.interp(np.log10(x), np.log10(r_CDM), np.log10(dens_CDM))

def potential_r(r,shell_mass,dens,dens_parameter,halo_radius,typ):
    def p(s):
        if s<r:
            return shell_mass(s,dens(s,dens_parameter))/r
        else:
            return shell_mass(s,dens(s,dens_parameter))/s
    if typ == 'small':
        potential_r, error = quad(p, 0, halo_radius)
    else:
        potential_r, error = quad(p, 0, halo_radius, epsrel = 0.01)
    return potential_r

def Jeans(r,dens,dens_parameter):
    enclose_mass,error = quad(lambda x:shell_mass(x,dens(r,dens_parameter)), 0,r)
    return -1*newton_G*enclose_mass/r**2*dens(r,dens_parameter)

def NFW_dens(r, dens_parameter):
    rho0 = dens_parameter[0]
    Rs = dens_parameter[1]
    return rho0/(r/Rs*(1+r/Rs)**2)


def density(x, r, d):
    return 10**np.interp(np.log10(x), np.log10(r), np.log10(d))

def enclosed_mass(x, r, d):
    mass, error = quad( lambda a: shell_mass(a, density(a, r, d)), 0, x, epsrel=1e-2)
    return mass

def gforce_mult_dens(x, r, d):
    return enclosed_mass(x, r, d)/x**2*density(x, r, d)*newton_G

def jeans_v(x, r, d):
    inte, error = quad(lambda a:gforce_mult_dens(a, r, d), x, r[-1], epsrel=1e-2)
    return (inte/density(x, r, d))**0.5*kpc2km


count = 0
max_v = 0

def plot_v(path, suffix, name, idx, jeans =False, soliton_v =False):

    df_halo_parameter = pd.read_csv( path+'../../Halo_Parameter_%d'%halo , sep = r'\s+' , header = 0 , index_col='#')
    current_time_z = df_halo_parameter['time'][idx]
    current_time_a = 1/(1+current_time_z)
    halo_radius = df_halo_parameter['halo_radius'][idx]/current_time_a
    core_radius = df_halo_parameter['core_radius_1'][idx]/current_time_a

    ### load data
    df_gamer_dens         = pd.read_csv( path+'AveDens_%06d%s'%(idx,suffix) , sep = r'\s+' , header = None ,skiprows =[0])
    df_gamer_v_bulk_r     = pd.read_csv( path+'AveVr_%06d%s'%(idx,suffix) , sep = r'\s+' , header = None ,skiprows =[0])
    df_gamer_v_bulk_theta = pd.read_csv( path+'AveVtheta_%06d%s'%(idx,suffix) , sep = r'\s+' , header = None ,skiprows =[0])
    df_gamer_v_bulk_phi   = pd.read_csv( path+'AveVphi_%06d%s'%(idx,suffix) , sep = r'\s+' , header = None ,skiprows =[0])
    df_gamer_v_qp_r       = pd.read_csv( path+'AveWr_%06d%s'%(idx,suffix) , sep = r'\s+' , header = None ,skiprows =[0])
    df_gamer_v_qp_theta   = pd.read_csv( path+'AveWtheta_%06d%s'%(idx,suffix) , sep = r'\s+' , header = None ,skiprows =[0])
    df_gamer_v_qp_phi     = pd.read_csv( path+'AveWphi_%06d%s'%(idx,suffix) , sep = r'\s+' , header = None ,skiprows =[0])
    # df_gamer_virial       = pd.read_csv( path+'VirSurf_%06d%s'%(idx,suffix) , sep = r'\s+' , header = None ,skiprows =[0])
    
    ### read data
    gamer_dens            = df_gamer_dens[2] # rho_bg
    gamer_r               = df_gamer_dens[0]*1000/h # kpccm
    gamer_shell_mass      = np.array(df_gamer_dens[6]) # code mass
 
    ### gamer v unit = 100 km/s v:bulk w:qp
    v_ave_r               = df_gamer_v_bulk_r[2]*100/current_time_a
    v_ave_theta           = df_gamer_v_bulk_theta[2]*100/current_time_a
    v_ave_phi             = df_gamer_v_bulk_phi[2]*100/current_time_a
    w_ave_r               = df_gamer_v_qp_r[2]*100/current_time_a
    w_ave_theta           = df_gamer_v_qp_theta[2]*100/current_time_a
    w_ave_phi             = df_gamer_v_qp_phi[2]*100/current_time_a
    v_sigma_r             = df_gamer_v_bulk_r[3]*100/current_time_a
    v_sigma_theta         = df_gamer_v_bulk_theta[3]*100/current_time_a
    v_sigma_phi           = df_gamer_v_bulk_phi[3]*100/current_time_a
    w_sigma_r             = df_gamer_v_qp_r[3]*100/current_time_a
    w_sigma_theta         = df_gamer_v_qp_theta[3]*100/current_time_a
    w_sigma_phi           = df_gamer_v_qp_phi[3]*100/current_time_a

    ### add hubble flow
    ### Recession Velocity = H(z)* distance, H^2(z) = H0^2(omega_M0(1+z)^3+omega_lambda)
    hubble_flow           = (((hubble0)**2*(omega_M0*(1+current_time_z)**3+omega_lambda))**0.5)*(gamer_r*current_time_a)
    v_ave_r               = v_ave_r + hubble_flow

    v_ave_2   = v_ave_r**2 + v_ave_theta**2 + v_ave_phi**2 + w_ave_r**2 + w_ave_theta**2 + w_ave_phi**2             # ave
    v_sigma_2 = v_sigma_r**2 + v_sigma_theta**2 + v_sigma_phi**2 + w_sigma_r**2 + w_sigma_theta**2 + w_sigma_phi**2 # sigma
    v2        = v_ave_2 + v_sigma_2                                                                                 # all
    v2bulk    = v_ave_r**2 + v_ave_theta**2 + v_ave_phi**2 + v_sigma_r**2 + v_sigma_theta**2 + v_sigma_phi**2       # bulk
    v2qp      = w_ave_r**2 + w_ave_theta**2 + w_ave_phi**2 + w_sigma_r**2 + w_sigma_theta**2 + w_sigma_phi**2       # qp

    ### soliton = 3.3 core_radius_1 (95% energy)
    soliton_range = (gamer_r <= core_radius*3.3)
    v2_core_both  = np.sum(v2[soliton_range]*gamer_shell_mass[soliton_range])/np.sum(gamer_shell_mass[soliton_range])
    v2_core_qp    = np.sum(v2qp[soliton_range]*gamer_shell_mass[soliton_range])/np.sum(gamer_shell_mass[soliton_range])
    v2_core_bulk  = np.sum(v2bulk[soliton_range]*gamer_shell_mass[soliton_range])/np.sum(gamer_shell_mass[soliton_range])

    ### inner = core_radius*6 ~ core_radius*8
    inner_range   = [core_radius*6, core_radius*8]
    filter_1      = (gamer_r >= inner_range[0])
    filter_2      = (gamer_r <= inner_range[1])
    v2_inner_both = np.sum(v2[filter_1 & filter_2]*gamer_shell_mass[filter_1 & filter_2])/np.sum(gamer_shell_mass[filter_1 & filter_2])
    v2_inner_qp   = np.sum(v2qp[filter_1 & filter_2]*gamer_shell_mass[filter_1 & filter_2])/np.sum(gamer_shell_mass[filter_1 & filter_2])
    v2_inner_bulk = np.sum(v2bulk[filter_1 & filter_2]*gamer_shell_mass[filter_1 & filter_2])/np.sum(gamer_shell_mass[filter_1 & filter_2])

    ### halo = core_radius*6 ~ halo_radius
    halo_range    = [core_radius*6, halo_radius]
    filter_1      = (gamer_r >= halo_range[0])
    filter_2      = (gamer_r <= halo_range[1])
    v2_halo_both  = np.sum(v2[filter_1 & filter_2]*gamer_shell_mass[filter_1 & filter_2])/np.sum(gamer_shell_mass[filter_1 & filter_2])
    v2_halo_qp    = np.sum(v2qp[filter_1 & filter_2]*gamer_shell_mass[filter_1 & filter_2])/np.sum(gamer_shell_mass[filter_1 & filter_2])
    v2_halo_bulk  = np.sum(v2bulk[filter_1 & filter_2]*gamer_shell_mass[filter_1 & filter_2])/np.sum(gamer_shell_mass[filter_1 & filter_2])

    ### ave = 0 ~ halo_radius
    all_range = (gamer_r<halo_radius)
    v2_ave_both = np.sum(v2[all_range]*gamer_shell_mass[all_range])/np.sum(gamer_shell_mass[all_range])
    v2_ave_QP   = np.sum(v2qp[all_range]*gamer_shell_mass[all_range])/np.sum(gamer_shell_mass[all_range])
    v2_ave_bulk = np.sum(v2bulk[all_range]*gamer_shell_mass[all_range])/np.sum(gamer_shell_mass[all_range])

    ### plot
    global count
    if count <3:
        plt.plot(gamer_r, v2**0.5, color = colors_tab10[count], lw = 2,label = 'v '+name)
        plt.plot(gamer_r, v2bulk**0.5*2**0.5, ':', color = colors_tab10[count*2+1], lw = 1,label = '$\sqrt{2}$ v bulk '+name)
        plt.plot(gamer_r, v2qp**0.5*2**0.5, '--', color = colors_tab10[count*2+2], lw = 1,label = '$\sqrt{2}$ v qp '+name)

    if jeans:

        file_path = "./velocity/jeans/jeans_v_%d_%d"%(halo,idx)

        if os.path.isfile(file_path):
            # File exists, proceed with reading
            df = pd.read_csv(file_path, sep = r'\s+' , header = 0)
            radius = df['radius(kpc)']
            j_v = df['jeans_v(km/s)']            

        else:

            df_dens = pd.read_csv(path+'../../prof_dens/Data_%06d_%d_profile_data'%(idx,halo), sep = r'\s+' , header = 0 )
            dens= np.array(df_dens['density(Msun/kpccm**3)'])/current_time_a**3
            radius = np.array(df_dens['radius(kpccm)'])*current_time_a

            j_v = np.zeros(len(radius))
            for i in range(len(j_v)):
                j_v[i] = jeans_v(radius[i], radius, dens)
            
            with open(file_path, 'w') as file:
                writer = csv.writer(file, delimiter = '\t')
                writer.writerow(['radius(kpc)', 'jeans_v(km/s)'])
                for i in range(len(radius)):
                    writer.writerow([radius[i], j_v[i]])
            

        # plt.plot(radius/current_time_a, sigma_jeans2**0.5*3**0.5, '-.',color = colors_tab10[count+8], lw = 1, label = 'v jeans')
        plt.plot(radius/current_time_a, j_v*3**0.5, '-.',color = colors_tab10[count+8], lw = 1, label = '$\sqrt{3}$ v jeans')

    # if count ==0:
        plt.plot([halo_radius, halo_radius],[0,310], '--', lw = 0.5, label = 'halo radius FDM')
        plt.plot([core_radius, core_radius],[0,310], '--', lw = 0.5, label = 'core radius FDM')
        plt.plot([core_radius*3.3, core_radius*3.3],[0,310],'--',lw = 0.3,color = 'gold', label = '3.3rc')
    count+=1

    core_Both[idx-idx_start]    = v2_core_both**0.5
    core_QP[idx-idx_start]      = v2_core_qp**0.5
    core_Bulk[idx-idx_start]    = v2_core_bulk**0.5
    inner_Both[idx-idx_start]   = v2_inner_both**0.5
    inner_QP[idx-idx_start]     = v2_inner_qp**0.5
    inner_Bulk[idx-idx_start]   = v2_inner_bulk**0.5
    halo_Both[idx-idx_start]    = v2_halo_both**0.5
    halo_QP[idx-idx_start]      = v2_halo_qp**0.5
    halo_Bulk[idx-idx_start]    = v2_halo_bulk**0.5
    average_Both[idx-idx_start] = v2_ave_both**0.5
    average_QP[idx-idx_start]   = v2_ave_QP**0.5
    average_Bulk[idx-idx_start] = v2_ave_bulk**0.5

    global max_v
    max_v = np.max(v2**0.5)


halo_vel_filename = 'halo_velocity_%d'%halo
writing_mode = 'append' if os.path.exists(halo_vel_filename) else 'new'

if writing_mode == 'new':
    with open( halo_vel_filename , 'w') as file:
        # writer = csv.writer(file, delimiter=' ')
        # writer.writerow(['#','a','rho0_FDM','Rs_FDM','Mass_NFW_fit_FDM','Radius_NFW_fit_FDM','Ep_NFW_fit_FDM',
        #                 'rho0_CDM','Rs_CDM','Mass_NFW_fit_CDM','Radius_NFW_fit_CDM','Ep_NFW_fit_CDM',
        #                 'Mass_FDM','Radius_FDM','Ep_FDM','Mass_CDM','Radius_CDM','Ep_CDM'])
        file.write(" #     time_a   core_Both     core_QP   core_Bulk  inner_Both    inner_QP  inner_Bulk   halo_Both    halo_QP    halo_Bulk    ave_Both      ave_QP    ave_Bulk\n")


writing_mode = 'append'

for idx in range(idx_start,idx_end+1,didx):

    group_n = 3
    count = 0

    colors_spring = plt.cm.spring(np.linspace(0, 1, group_n))
    colors_summer = plt.cm.summer(np.linspace(0, 1, group_n))
    colors_autumn = plt.cm.autumn(np.linspace(0, 1, group_n))
    colors_winter = plt.cm.winter(np.linspace(0, 1, group_n))
    colors_tab10 = plt.cm.tab10(np.linspace(0, 1, 10))

    plot_v('velocity/output_%d/'%halo, '', '', idx, jeans=True)
   

    plt.legend(loc = 'upper left', fontsize = 8)
    plt.xlim(2e-1,4e2)
    plt.ylim(0,np.ceil(max_v/50)*50)
    plt.xlabel('radius (kpccm)')
    plt.ylabel('$\sigma$ (km/s)')
    plt.xscale('log')
    plt.savefig('velocity/v2_of_radius/v2_of_r_%2d_%d.png'%(idx,halo), dpi = 150)
    plt.close()
    

time_z = df_FDM['time'][:]
time_a = 1.0/(1.0+time_z)
time_a = time_a[list(range(idx_start,idx_end+1,didx))]

    
with open(halo_vel_filename , 'a') as file:
    writer = csv.writer(file, delimiter='\t')
    for i in range(idx_start,idx_end+1,didx):
        writer.writerow([i,'%.3e'%time_a[i], '%.5f'%core_Both[i-idx_start],'%.5f'%core_QP[i-idx_start],'%.5f'%core_Bulk[i-idx_start],\
            '%.5f'%inner_Both[i-idx_start], '%.5f'%inner_QP[i-idx_start], '%.5f'%inner_Bulk[i-idx_start],\
            '%.5f'%halo_Both[i-idx_start], '%.5f'%halo_QP[i-idx_start], '%.5f'%halo_Bulk[i-idx_start],\
            '%.5f'%average_Both[i-idx_start], '%.5f'%average_QP[i-idx_start], '%.5f'%average_Bulk[i-idx_start]])
