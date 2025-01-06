#!/usr/bin/env python3

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import yt

figure_width = 3.375
font_size_regular = 8
font_size_small = 6
font_size_large = 8
line_width = 0.6

path = '../'
halo = 1
idx = 68

ds = yt.load(path+'/../Data_0000%d'%idx)
### constant
kpc2km                = (1*yt.units.kpc).to('km').d
m_sun2_kg             = (1*yt.units.Msun).to('kg').d
omega_M0              = ds.omega_matter
omega_lambda          = ds.omega_lambda
h                     = ds.hubble_constant #dimensionless Hubble parameter 0.1km/(s*kpc)
hubble0               = h*0.1 # km/(s*kpc)
newton_G              = ds.units.newtons_constant.to('(kpc**3)/(s**2*Msun)').d  #(kpc^3)/(s^2*Msun)
background_density_0  = (1*ds.units.code_density).to("Msun/kpc**3").d
particle_mass         = (ds.parameters['ELBDM_Mass']*ds.units.code_mass).to('eV/c**2').d
code_mass             = (1*ds.units.code_mass).to("kg").d
CDM_particle_mass     = 170399.32174374 # Msun


def plot_v(path, suffix, name):

    df_halo_parameter = pd.read_csv( path+'/../../Halo_Parameter_%d'%halo , sep = r'\s+' , header = 0 , index_col='#')
    current_time_z = df_halo_parameter['time'][idx]
    current_time_a = 1/(1+current_time_z)
    halo_radius = df_halo_parameter['halo_radius'][idx]
    core_radius = df_halo_parameter['core_radius_1'][idx]
    ### load data
    df_gamer_dens = pd.read_csv( path+'AveDens_%06d%s'%(idx,suffix) , sep = r'\s+' , header = None ,skiprows =[0])
    df_gamer_v_bulk_r = pd.read_csv( path+'AveVr_%06d%s'%(idx,suffix) , sep = r'\s+' , header = None ,skiprows =[0])
    df_gamer_v_bulk_theta = pd.read_csv( path+'AveVtheta_%06d%s'%(idx,suffix) , sep = r'\s+' , header = None ,skiprows =[0])
    df_gamer_v_bulk_phi = pd.read_csv( path+'AveVphi_%06d%s'%(idx,suffix) , sep = r'\s+' , header = None ,skiprows =[0])
    df_gamer_v_qp_r = pd.read_csv( path+'AveWr_%06d%s'%(idx,suffix) , sep = r'\s+' , header = None ,skiprows =[0])
    df_gamer_v_qp_theta = pd.read_csv( path+'AveWtheta_%06d%s'%(idx,suffix) , sep = r'\s+' , header = None ,skiprows =[0])
    df_gamer_v_qp_phi = pd.read_csv( path+'AveWphi_%06d%s'%(idx,suffix) , sep = r'\s+' , header = None ,skiprows =[0])
    
    ### read data
    gamer_dens = df_gamer_dens[2] # rho_bg
    gamer_r = df_gamer_dens[0]*1000/h*current_time_a # kpc
    gamer_shell_mass = np.array(df_gamer_dens[6]) # code mass
    
    ### gamer v unit = 100 km/s v:bulk w:qp
    v_ave_r     = df_gamer_v_bulk_r[2]*100/current_time_a
    v_ave_theta = df_gamer_v_bulk_theta[2]*100/current_time_a
    v_ave_phi   = df_gamer_v_bulk_phi[2]*100/current_time_a
    w_ave_r     = df_gamer_v_qp_r[2]*100/current_time_a
    w_ave_theta = df_gamer_v_qp_theta[2]*100/current_time_a
    w_ave_phi   = df_gamer_v_qp_phi[2]*100/current_time_a
    v_sigma_r   = df_gamer_v_bulk_r[3]*100/current_time_a
    v_sigma_theta = df_gamer_v_bulk_theta[3]*100/current_time_a
    v_sigma_phi = df_gamer_v_bulk_phi[3]*100/current_time_a
    w_sigma_r   = df_gamer_v_qp_r[3]*100/current_time_a
    w_sigma_theta = df_gamer_v_qp_theta[3]*100/current_time_a
    w_sigma_phi = df_gamer_v_qp_phi[3]*100/current_time_a

    ### add hubble flow
    ### Recession Velocity = H(z)* distance, H^2(z) = H0^2(omega_M0(1+z)^3+omega_lambda)
    hubble_flow = (((hubble0)**2*(omega_M0*(1+current_time_z)**3+omega_lambda))**0.5)*(gamer_r*current_time_a)
    v_ave_r = v_ave_r + hubble_flow

    v_ave_2 = v_ave_r**2 + v_ave_theta**2 + v_ave_phi**2 + w_ave_r**2 + w_ave_theta**2 + w_ave_phi**2
    v_sigma_2 = v_sigma_r**2 + v_sigma_theta**2 + v_sigma_phi**2 + w_sigma_r**2 + w_sigma_theta**2 + w_sigma_phi**2
    v2 = v_ave_2 + v_sigma_2
    v2bulk = v_ave_r**2 + v_ave_theta**2 + v_ave_phi**2 + v_sigma_r**2 + v_sigma_theta**2 + v_sigma_phi**2
    v2qp = w_ave_r**2 + w_ave_theta**2 + w_ave_phi**2 + w_sigma_r**2 + w_sigma_theta**2 + w_sigma_phi**2

    ### soliton = 3.3 core_radius_1 (95% energy)
    soliton_range = (gamer_r <= core_radius*3.3)
    v2_core_qp = np.sum(v2qp[soliton_range]*gamer_shell_mass[soliton_range])/np.sum(gamer_shell_mass[soliton_range])

    ### plot

    ax.plot(gamer_r, v2qp**0.5, '-', color = 'mediumslateblue', lw = 1,label = 'v qp '+name)
    ax.plot(gamer_r, v2bulk**0.5, '--', color = 'darkorange', lw = 1,label = 'v bulk '+name)
 
    # w_s
    ax.plot([2e-1,core_radius*3.3], [v2_core_qp**0.5,v2_core_qp**0.5], ':', color = 'purple', lw = 0.6,label = 'v core '+name)
    ax.annotate(r'${\langle w \rangle}_{\mathrm{s}}$', xy=(1.5, v2_core_qp**0.5), xycoords='data', xytext=(1.3, v2_core_qp**0.5*1.02), textcoords='data', color = 'purple',fontsize = font_size_regular)

    plt.plot([halo_radius, halo_radius],[0,310], 'seagreen', ls = '--', lw = 0.6, label = 'halo radius FDM')
    ax.annotate('halo radius', xy=(0.85, v2_core_qp**0.5), xycoords='data', xytext=(78, 115), textcoords='data', color = 'seagreen',fontsize = font_size_small)

    ### region
    ax.fill_between([core_radius*0, core_radius*3.3],[19,19],[140,140], color = 'gold', alpha = 0.2)
    ax.annotate(r'soliton: $w \gg v, {\langle w \rangle}_\mathrm{s} \sim w_\mathrm{h, in}$'+'\n(thermal equilibrium)', xy=(1, 120), xycoords='data', xytext=(1.8, 3), ha="center", va="bottom", color = 'darkorange',fontsize = font_size_small )

    ax.fill_between([core_radius*6, 30],[105,105],[135,135], color = 'skyblue', alpha = 0.3)
    ax.annotate('inner halo:\n'+ r'$w_\mathrm{h, in} \sim v_\mathrm{h, in} \sim const$'+'\n(equipartition and isothermal)', xy=(3, 113), xycoords='data', xytext=(13, 83), ha="center", va="bottom", color = 'royalblue',fontsize = font_size_small )


    ### plot Jeans velocity
    df_jeans = pd.read_csv( path+"../jeans/jeans_v_%d_%d"%(halo,idx), sep = r'\s+' , header = 0)
    radius = df_jeans['radius(kpc)']
    j_v = df_jeans['jeans_v(km/s)']
    plt.plot(radius, np.array(j_v)*3**0.5/2**0.5, '-.',color = 'lightskyblue', lw = 1, label = 'v jeans')


fig, ax = plt.subplots(figsize=(figure_width, 2.5))

plot_v(path+'/velocity/output_%d/'%halo, '', 'FDM')

# plot legend 
from matplotlib.legend_handler import HandlerTuple
from matplotlib.lines import Line2D

custom_legend = [
    Line2D([0], [0], ls = '-', color='mediumslateblue',lw = 1),
    Line2D([0], [0], ls='--', color='darkorange',lw = 1),
    # Line2D([0], [0], ls=':', color='violet',lw = 1),
    Line2D([0], [0], ls='-.',color = 'lightskyblue',lw = 1)
]

legend_labels = [ '$w$', '$v$', '$\sqrt{\dfrac{3}{2}}$ Jeans']

leg = ax.legend(custom_legend, legend_labels, loc=(0.65,0.03), fontsize = font_size_small+1,handler_map={tuple: HandlerTuple(ndivide=None)}, handlelength=1.5)
leg.get_frame().set_linewidth(0.3)

### set plot
plt.xscale('log')
plt.xlim(4e-1,4e2)
plt.ylim(0,150)

plt.xlabel('radius $r$ (kpc)', fontsize = font_size_regular, labelpad= 1)
plt.ylabel('velocity (km $\mathrm{s}^{-1}$)', fontsize = font_size_regular, labelpad= 1)
plt.setp(ax.get_xticklabels(), fontsize = font_size_regular)
plt.setp(ax.get_yticklabels(), fontsize = font_size_regular)

ax.tick_params(bottom=True, top=True, left=True, right=True,  which='both',direction ='in')
ax.tick_params('both', length=2, width=0.5, which='major')
ax.tick_params('both', length=1, width=0.2, which='minor')

# spines
for axis in ['top','bottom','left','right']:
    ax.spines[axis].set_linewidth(0.5)

plt.tight_layout(pad=0.1, w_pad=0, h_pad=0)
plt.savefig('fig_velocity.pdf', dpi = 150)
plt.close()

