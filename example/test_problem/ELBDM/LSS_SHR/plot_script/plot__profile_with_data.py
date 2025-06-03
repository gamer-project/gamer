#!/usr/bin/env python3

import argparse
import sys

import fit_NFW
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import yt

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
   print(sys.argv[t], end = ' ')
print( '' )
print( '-------------------------------------------------------------------\n' )


idx_start   = args.idx_start
idx_end     = args.idx_end
didx        = args.didx
prefix      = args.prefix
halo        = args.halo

# get background_density_0
ds            = yt.load('../'+prefix+'Data_%06d'%idx_start)
background_density_0 = (1*ds.units.code_density).to("Msun/kpc**3").d

df_Halo_Parameter = pd.read_csv( 'Halo_Parameter_%d'%halo , sep = r'\s+' , header = 0 , index_col='#')


def plot_profile(path,name, core_is_true = True, NFW = True):

    def soliton(x):   
        return ((1.9*(current_time_a**-1)*(particle_mass/10**-23)**-2*(core_radius_1**-4))/((1 + 9.1*10**-2*(x/core_radius_1)**2)**8))*10**9

    # read data
    df = pd.read_csv( path+'/prof_dens/Data_%06d_%d_profile_data'%(idx,halo) , sep = '\t' , header = 0 )
    df_halo_parameter = pd.read_csv( path+'/Halo_Parameter_%d'%halo , sep = r'\s+' , header = 0 , index_col='#')

    current_time_z = df_halo_parameter['time'][idx]
    current_time_a = 1/(current_time_z+1)
    radius = df['radius(kpccm)'][:]
    density = df['density(Msun/kpccm**3)'][:]
    dens = density
    # dens = density/background_density_0
    halo_radius = df_halo_parameter['halo_radius'][idx]/current_time_a
    # plot
    plt.plot( radius, dens, '.', label=name)
    plt.plot( [halo_radius, halo_radius], [1e-1, 1e9], '--', label= name+' halo radius')

    # plot core
    if (core_is_true):
        particle_mass = df_halo_parameter['mass'][idx]
        core_radius_1 = df_halo_parameter['core_radius_1'][idx]/current_time_a
        x = np.logspace(-1, 3, num=50)
        plt.plot( x, soliton(x) )

    # plot NFW
    if (NFW):
        rho0 = df_halo_parameter['rho0'][idx]
        Rs = df_halo_parameter['Rs'][idx]
        x = np.logspace(-2, 3, num = 50)
        plt.plot(x/current_time_a, fit_NFW.NFW_dens(x, (rho0, Rs))*current_time_a**3)

for idx in range(idx_start, idx_end+1, didx):

    current_time_z = df_Halo_Parameter['time'][idx]
    current_time_a = 1/(current_time_z+1)

    plot_profile('./', 'FDM')
    # plot_profile(compare_path, 'N-body', False, False)

    plt.xscale("log")
    plt.yscale("log")
    plt.xlim(1e-1,1e3)
    plt.ylim(1e0,1e10)
    plt.ylabel('$\\rho(r)$ (Msun/kpccm)')
    # plt.ylabel('$\\rho(r)/\\rho_{m0}$')
    plt.xlabel('radius (kpccm)')
    plt.legend(loc = 'upper right')
    plt.title('z = %.2e'%(current_time_z), fontsize=12) 


    FileOut = 'prof_dens/fig_profile_density_%d_%02d'%(halo, idx)+'.png'
    plt.savefig( FileOut, dpi = 200)
    plt.close()
    
# FileOut = 'fig_profile_density'+'.png'
# plt.savefig( FileOut)
