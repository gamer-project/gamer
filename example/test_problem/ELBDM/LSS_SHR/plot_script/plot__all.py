#!/usr/bin/env python3

import argparse
import sys

import matplotlib.pyplot as plt
import pandas as pd
import shr.SHR as SHR
import yt

### load the command-line parameters
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
   print(sys.argv[t], end = ' '),
print( '' )
print( '-------------------------------------------------------------------' )

idx_start   = args.idx_start
idx_end     = args.idx_end
didx        = args.didx
prefix      = args.prefix
halo        = args.halo


### constant
kpc2km                = (1*yt.units.kpc).to('km').d
newton_G              = (1*yt.units.newtons_constant).to('(kpc**3)/(s**2*Msun)').d


shr_calculator = SHR.SHR_calculator('planck18')

df_vel = pd.read_csv('halo_velocity_%d'%halo, sep = r'\s+' , header = 0, index_col = '#').loc[idx_start:idx_end:didx]
df_halo = pd.read_csv('Halo_Parameter_%d'%halo, sep = r'\s+' , header = 0, index_col = '#').loc[idx_start:idx_end:didx]

m22 = df_halo['mass'].iloc[0]/1e-22
halo_mass = df_halo['halo_mass']
halo_radius = df_halo['halo_radius']
core_mass = df_halo['core_mass_1']
time_a = df_vel['time_a']
soliton_m_div_v = SHR.soliton_m_div_v(m22)

### plot
Ep_theo_FDM = -SHR.get_Ep(halo_mass, halo_radius, shr_calculator.concentration_para_FDM(df_halo['halo_mass'], 0, m22), 'NFW')
Ep_fit_FDM = df_halo['fit_Ep']*df_halo['NFW_mass']**2/df_halo['NFW_radius']*newton_G
Ep_sim_FDM = df_halo['sim_Ep']*halo_mass**2/halo_radius*newton_G

ms_Ep_theo = (Ep_theo_FDM/halo_mass)**0.5*soliton_m_div_v
ms_Ep_fit = (Ep_fit_FDM/df_halo['NFW_mass'])**0.5*soliton_m_div_v
ms_Ep_sim = (Ep_sim_FDM/halo_mass)**0.5*soliton_m_div_v

plt.plot(time_a, ms_Ep_fit/ms_Ep_theo,':', color = 'lime', label='Ep fit/theo')
plt.plot(time_a, ms_Ep_sim/ms_Ep_fit,'--', color = 'mediumseagreen', label='Ep sim/fit')
plt.plot(time_a, ms_Ep_sim/ms_Ep_theo, color = 'green', label='Ep sim/theo')


ms_Ek = df_vel['ave_Both']/kpc2km*soliton_m_div_v
plt.plot(time_a, ms_Ek/ms_Ep_sim,'--', color = 'gold', label='Ek/Ep')

w_ave = df_vel['ave_QP']
B_ave = df_vel['ave_Both']
plt.plot(time_a, w_ave/B_ave/0.5**0.5, ':',color = 'pink', label='QP ratio/√2')

w_inner = df_vel['inner_QP']
plt.plot(time_a, (w_inner/w_ave/shr_calculator.nonisothermality(df_halo['halo_radius']/df_halo['Rs'])), color = 'red', label='nonisothermality')

w_s = df_vel['core_QP']
plt.plot(time_a,(w_s/w_inner), color = 'blue', label = 'w_s/w_inner')
B_s = df_vel['core_Both']
B_inner = df_vel['inner_Both']
# plt.plot(time_a,(B_s/B_inner), color = 'skyblue', label = 'Both_s/Both_inner')

ms_ws = df_vel['core_QP']/kpc2km*soliton_m_div_v
plt.plot(time_a, (core_mass/ms_ws),':', color = 'gray', label = 'fidelity')

plt.plot(time_a,(core_mass/shr_calculator.revised_theo_c_FDM_Ms(0,halo_mass,m22)),'-.', color = 'black', label = 'ms sim/theo')

plt.legend(loc = 'lower left', fontsize = 7)
plt.yscale('log')
plt.ylabel('ratio')
plt.xlabel('scale factor')
plt.ylim(0.25, 4)
plt.tight_layout()
plt.savefig('all_factor_%d.png'%halo, bbox_inches='tight', dpi = 150)
plt.close()


