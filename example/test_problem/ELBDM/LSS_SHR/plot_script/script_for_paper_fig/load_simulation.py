#!/usr/bin/env python3

import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import yt
from numpy import linalg as LA
from scipy.stats import linregress

sys.path.append('../shr/')
import SHR

ds = yt.load('../../Data_000000')
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
zeta_0 = (18*np.pi**2 + 82*(omega_M0 - 1) - 39*(omega_M0 - 1)**2)/omega_M0

shr_calculator = SHR.SHR_calculator('planck18')


### halo            0              1            2              3         4              5                   6               7                  8                 9       10             11              12              13                      14            15          16                 17             18              19      20
# factor_name = ['time_a', 'halo_mass', 'Ep_fit_theo', 'Ep_sim_fit', 'Ek2_Ep', 'v_inner_ave', 'v_core_Both_inner', 'v_core_QP_Both', 'Ep_sim_theo', 'v_core_QP_inner_QP', 'v_ave', 'v_soliton_QP', 'core_mass', 'v_soliton_Both', 'v_inner_QP_avg_QP', 'v_inner_QP_Both', 'v_avg_QP_Both', 'v_ave_fit', 'v_ave_QP', 'v_inner_QP', 'v_inner']
factor_name = {'time_a': 'time a', 'time_z' :'time z', 'halo_mass': 'halo mass ($M_\odot$)', 'halo_mass_CDM': 'halo mass CDM ($M_\odot$)', 'particle_mass':'$m_\psi (eV/c^2)$',\
               'halo_radius':'halo radius (kpc)', 'halo_radius_CDM':'halo radius CDM (kpc)', 'core_mass': 'core mass ($M_\odot$)', 'core_radius': 'core radius (kpc)',\
               'c_theo':'theoretical concentration parmeter', 'c_fit':'fit concentration parmeter','c_fit_c0':'$\dfrac{c_{fit}}{c_{fit,\ z=0}}$', 'c_CDM':'gadget concentration parmeter','c_theo_CDM':'theoretical concentration parmeter CDM',\
               'rho0':'rho0', 'Rs':'Rs',\
               'Ep_theo':'$E_{p,\ theo}$', 'Ep_fit':'$E_{p,\ fit}$', 'Ep_sim':'$E_{p,\ sim}$',\
               'Ep_fit_theo':'$\dfrac{E_{p,\ fit}}{E_{p,\ theo}}$', 'Ep_sim_fit':'$\dfrac{E_{p,\ sim}}{E_{p,\ fit}}$', 'Ep_sim_theo':'$\dfrac{E_{p,\ sim}}{E_{p,\ theo}}$',\
               'Ek2_Ep':'$\dfrac{2E_{k,\ sim}}{E_{p,\ sim}}$',\
               'v_soliton_QP': '$v_{soliton,\ QP}$', 'v_soliton_Both': '$v_{soliton,\ Both}$', 'soliton_QP_ratio': '$\dfrac{v_{soliton,\ QP}}{v_{soliton,\ Both}}$',\
               'v_inner_QP': '$v_{inner,\ QP}$', 'v_inner_Both': '$v_{inner,\ Both}$', 'inner_QP_ratio': '$\dfrac{v_{inner,\ QP}}{v_{inner,\ Both}}$',\
               'v_ave_QP': '$v_{ave,\ QP}$', 'v_ave_Both': '$v_{ave,\ Both}$', 'ave_QP_ratio': '$\dfrac{v_{ave,\ QP}}{v_{ave,\ Both}}$',\
               'v_soliton_QP_inner_QP': '$\dfrac{v_{soliton,\ QP}}{v_{inner,\ QP}}$', 'v_soliton_Both_inner_Both': '$\dfrac{v_{soliton,\ Both}}{v_{inner,\ Both}}$',\
               'v_inner_QP_ave_QP': '$\dfrac{v_{inner,\ QP}}{v_{ave,\ QP}}$', 'v_inner_Both_ave_Both': '$\dfrac{v_{inner,\ Both}}{v_{ave,\ Both}}$',\
                'count': 'count', 'halo_age_Gyr': 'halo age (Gyr)'
               }


def load( path, halo, start, end, name, tag='N' ):


    df_halo = pd.read_csv(path+'/Halo_Parameter_%d'%halo, sep = r'\s+' , header = 0, index_col = '#').loc[start:end]
    df_vel = pd.read_csv(path+'/halo_velocity_%d'%halo, sep = r'\s+' , header = 0, index_col = '#').loc[start:end]
    
    time_a        = df_vel['time_a']
    time_z        = df_halo['time']
    particle_mass = df_halo['mass']
    m22           = particle_mass/1e-22
    halo_mass     = df_halo['halo_mass']
    halo_radius   = df_halo['halo_radius']
    core_mass     = df_halo['core_mass_1']
    core_radius   = df_halo['core_radius_1']
    mcrc          = df_halo['core_radius_2']*df_halo['core_mass_2']

    rho0        = df_halo['rho0']
    Rs          = df_halo['Rs']
    c_theo      = shr_calculator.concentration_para_FDM(halo_mass, time_z, m22)
    c_fit       = halo_radius/Rs

    Ep_theo_FDM = newton_G*SHR.get_Ep(halo_mass, halo_radius, c_theo, 'NFW')
    Ep_fit_FDM  = newton_G*SHR.get_Ep(halo_mass, halo_radius, c_fit, 'NFW')
    Ep_sim_FDM  = newton_G*(df_halo['sim_Ep']*halo_mass**2/halo_radius)


    v_ave_Both   = df_vel['ave_Both']
    v_halo_Both  = df_vel['halo_Both']
    v_inner_Both = df_vel['inner_Both']
    v_core_Both  = df_vel['core_Both']
    v_ave_QP     = df_vel['ave_QP']
    v_halo_QP    = df_vel['halo_QP']
    v_inner_QP   = df_vel['inner_QP']
    v_core_QP    = df_vel['core_QP']

    Ek           = 0.5*v_ave_Both**2*halo_mass/kpc2km**2

    Ep_fit_theo = (Ep_fit_FDM/Ep_theo_FDM)
    Ep_sim_fit  = (Ep_sim_FDM/Ep_fit_FDM)
    Ep_sim_theo = (Ep_sim_FDM/Ep_theo_FDM)
    Ek2_Ep      = (2*Ek/Ep_sim_FDM)

    soliton_QP_ratio = v_core_QP/v_core_Both
    inner_QP_ratio   = v_inner_QP/v_inner_Both
    halo_QP_ratio    = v_halo_QP/v_halo_Both
    ave_QP_ratio     = v_ave_QP/v_ave_Both

    v_soliton_Both_inner_Both = v_core_Both/v_inner_Both
    v_soliton_QP_inner_QP     = v_core_QP/v_inner_QP
    v_inner_QP_ave_QP         = v_inner_QP/v_ave_QP
    v_inner_Both_ave_Both     = v_inner_Both/v_ave_Both

    # make a dataframe
    df = pd.DataFrame()
    df['time_a'] = time_a
    df['time_z'] = time_z
    df['particle_mass'] = particle_mass
    df['m22'] = m22
    df['halo_mass'] = halo_mass
    df['halo_radius'] = halo_radius
    df['core_mass'] = core_mass
    df['core_radius'] = core_radius
    df['mcrc'] = mcrc
    df['c_theo'] = c_theo
    df['c_fit'] = c_fit
    # if len(df[df['c_fit']<1])>0:
    #     print(name)
    #     print(df[df['c_fit']<1])
    df['c_fit_c0'] = c_fit/c_fit[end]
    df['rho0'] = rho0
    df['Rs'] = Rs
    df['Ep_theo'] = Ep_theo_FDM
    df['Ep_fit'] = Ep_fit_FDM
    df['Ep_sim'] = Ep_sim_FDM
    df['Ep_fit_theo'] = Ep_fit_theo
    df['Ep_sim_fit'] = Ep_sim_fit
    df['Ep_sim_theo'] = Ep_sim_theo
    df['Ek2_Ep'] = Ek2_Ep
    df['v_soliton_QP'] = v_core_QP
    df['v_soliton_Both'] = v_core_Both
    df['soliton_QP_ratio'] = soliton_QP_ratio
    df['v_inner_QP'] = v_inner_QP
    df['v_inner_Both'] = v_inner_Both
    df['inner_QP_ratio'] = inner_QP_ratio
    df['v_halo_Both'] = v_halo_Both
    df['v_halo_QP'] = v_halo_QP
    df['halo_QP_ratio'] = halo_QP_ratio
    df['v_ave_QP'] = v_ave_QP
    df['v_ave_Both'] = v_ave_Both
    df['ave_QP_ratio'] = ave_QP_ratio
    df['v_soliton_QP_inner_QP'] = v_soliton_QP_inner_QP
    df['v_soliton_Both_inner_Both'] = v_soliton_Both_inner_Both
    df['v_inner_QP_ave_QP'] = v_inner_QP_ave_QP
    df['v_inner_Both_ave_Both'] = v_inner_Both_ave_Both

    df.name = name
    df = check_virial(df)
    df['count']        = df.index-df.index[0]

    # tag (for fake halo, can exclude if needed)
    # fake halo
    if tag[0] == 'Y':
        df.fake = 1
    else:
        df.fake = 0

    
    df.name = name

    frames.append(df)

    return df


def load_sim():

    global frames
    frames = []

    df = load('../', 1, 59, 68, 'light', 'N')


    return frames, factor_name

def check_virial(df):
    filtered_df = df[df['Ek2_Ep'] > 1.35]
    if not filtered_df.empty:
        start = filtered_df.index[-1]+1
        # print(df.name, start, df.index[0])
        end = df.index[-1]
        df = df.loc[start:end]
    return df

def bin(num_bins, x_axis, scale, df, min=None, max=None):
    # Define the bins for mass

    min = min or df[x_axis].min()
    max = max or df[x_axis].max()

    if scale == 'log':
        interval = np.logspace(np.log10(min), np.log10(max), num_bins+1)
        x_ticks = (interval[:-1]*interval[1:])**0.5
    elif scale == 'lin':
        interval = np.linspace(min, max, num_bins+1)
        x_ticks = (interval[:-1]+interval[1:])*0.5
    else:
        raise ValueError('scale must be "log" or "lin"')
    return interval, x_ticks


colors_tab20 = plt.cm.tab20(np.linspace(0, 1, 20))

global i

def fit_regres(x, y, fit_type, i=0):
    x_log, y_log = fit_type
    plot_x = np.linspace(x.min(), x.max(), 100)
    if x_log and y_log:
        res = linregress(np.log10(x), np.log10(y))
        # plt.plot(10**np.log10(plot_x), 10**(res.intercept + res.slope*np.log10(plot_x)), ':', color = colors_tab20[i])
        # plt.annotate('linregress: ratio = 10**(%.2e*log(x_axis)+%.2e)'%(res.slope, res.intercept), xy = (0.30, 0.05+i*0.04), xycoords='axes fraction', color = colors_tab20[i],size = 8 )
    elif x_log:
        res = linregress(np.log10(x), y)
        # plt.plot(10**np.log10(plot_x), res.intercept + res.slope*np.log10(plot_x), ':', color = colors_tab20[i])
        # plt.annotate('linregress: ratio = %.2e*log(x_axis)+%.2e'%(res.slope, res.intercept), xy = (0.30, 0.05+i*0.04), xycoords='axes fraction',color = colors_tab20[i], size = 8 )
    elif y_log:
        res = linregress(x, np.log10(y))
        # plt.plot(plot_x, 10**(res.intercept + res.slope*plot_x), ':', color = colors_tab20[i])
        # plt.annotate('linregress: log(ratio) = %.2e*x_axis+%.2e'%(res.slope, res.intercept), xy = (0.30, 0.05+i*0.04), xycoords='axes fraction',color = colors_tab20[i], size = 8 )
    else:
        res = linregress(x, y)
        # plt.plot(plot_x, res.intercept + res.slope*plot_x, ':', color = colors_tab20[i])
        # plt.annotate('linregress: ratio = %.2e*x_axis+%.2e'%(res.slope, res.intercept), xy = (0.30, 0.05+i*0.04), xycoords='axes fraction',color = colors_tab20[i], size = 8 )
    return res
    
markers      = ['o', 's', 'v', 'X', '^', 'P', 'd', '*', 'p']
markers_size = [ 4,   4,   4,   5,   5,   5,   4,   6,   5]

colors_jet = plt.cm.jet(np.linspace(0, 1, 5))
red = plt.cm.Reds(np.linspace(0, 1, 40))
blue = plt.cm.Blues(np.linspace(0, 1, 40))
green = plt.cm.Greens(np.linspace(0, 1, 45))

def plot_single(frames, x_axis, fac, ax, group_m22 = True, time_range=None):
    i=0
    for frame in frames:
        name = frame.name
        frame = frame[frame['time_a']>time_range] if time_range else frame
        if frame.empty:
            print(name+' is empty')
            i+=1
            continue
        frame.name = name
        if frame['particle_mass'].iloc[-1] == 8.0e-23:
            if group_m22 == False:
                ax.plot(frame[x_axis], frame[fac], marker=markers[i % len(markers)], ls = 'None', mfc = 'white', mec = colors_jet[(i)%len(colors_jet)],\
                        mew = 1, markersize = 6, label = frame.name)
            else: 
                # ax.plot(frame[x_axis], frame[fac], marker='^', ls = 'None', mfc = green[20], mec = green[30],\
                #         mew = 0.2, markersize = 8, label = frame.name)
                ax.plot(frame[x_axis], frame[fac], marker='^', ls = 'None', mfc = green[20], mec = 'None',\
                        mew = 0.2, markersize = 6, label = frame.name)
            
        elif frame['particle_mass'].iloc[-1] == 1.0e-23:
            if group_m22 == False:
                ax.plot(frame[x_axis], frame[fac], marker=markers[i % len(markers)], ls = 'None', mfc = colors_jet[(i)%len(colors_jet)], mec = 'white',\
                        mew = 0.5, markersize = 6, label = frame.name)
            else:
                # ax.plot(frame[x_axis], frame[fac], marker='o', ls = 'None', mfc = red[20], mec = red[30],\
                #         mew = 0.2, markersize = 6, label = frame.name)
                ax.plot(frame[x_axis], frame[fac], marker='o', ls = 'None', mfc = red[20], mec = 'None',\
                        mew = 0.2, markersize = 4, label = frame.name)
        
        else:
            if group_m22 == False:
                ax.plot(frame[x_axis], frame[fac], marker=markers[i % len(markers)],ls = 'None', mfc = colors_jet[(i+2)%len(colors_jet)], mec = colors_jet[(i)%len(colors_jet)],\
                        mew = 0.4, markersize = markers_size[i%len(markers_size)], label = frame.name)
            else:
                # ax.plot(frame[x_axis], frame[fac], marker='s', ls = 'None', mfc = blue[20], mec = blue[30],\
                #         mew = 0.2, markersize = 6, label = frame.name)
                ax.plot(frame[x_axis], frame[fac], marker='s', ls = 'None', mfc = blue[20], mec = 'None',\
                        mew = 0.2, markersize = 4, label = frame.name)
        i+=1


def plot_m22_redshift(combined_df, x_axis, fac, ax):

    m22_8 = combined_df[combined_df['particle_mass'] == 8.0e-23]
    m22_2 = combined_df[combined_df['particle_mass'] == 2.0e-23]
    m22_1 = combined_df[combined_df['particle_mass'] == 1.0e-23]


    for j in range(len(m22_8)):
        ax.plot(m22_8[x_axis].iloc[j], m22_8[fac].iloc[j], marker='^', ls = 'None', mfc = green[m22_8.index[j]-25], mec = 'None',mew = 0, markersize = 4)
    for j in range(len(m22_2)):
        ax.plot(m22_2[x_axis].iloc[j], m22_2[fac].iloc[j], marker='s', ls = 'None', mfc = blue[m22_2.index[j]-30], mec = 'None',mew = 0, markersize = 4)
    for j in range(len(m22_1)):
        ax.plot(m22_1[x_axis].iloc[j], m22_1[fac].iloc[j], marker='o', ls = 'None', mfc = red[m22_1.index[j]-30], mec = 'None',mew = 0, markersize = 4)

