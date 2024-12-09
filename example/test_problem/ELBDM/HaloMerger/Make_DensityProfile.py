import numpy as np
import matplotlib.pyplot as plt


###################################################################################################################
# code units (in cgs)
UNIT_L     = 4.4366320345075490e+24
UNIT_D     = 2.5758579724476994e-30

# length units
Const_cm   = 1.0;
Const_pc   = 3.08567758149e18;
Const_kpc  = 1.0e3*Const_pc;

# mass units
Const_g    = 1.0;
Const_kg   = 1.0e3*Const_g;
Const_Msun = 1.9885e33;
###################################################################################################################


###################################################################################################################
# density profile
def Halo_fitting_analytical_dens(r, rho_0, r_0, alpha, beta, gamma):
    x = r / r_0
    return rho_0 * ( x**(-gamma) ) * ( ( 1.0+(x**alpha) )**( (gamma-beta)/alpha ) )

def Soliton_rho_c_from_r_c(m22, r_c):
    return 1.945e7/( m22**2 * (r_c*UNIT_L/Const_kpc)**4 )*(Const_Msun/Const_kpc**3)/UNIT_D

def Soliton_fitting_analytical_dens(r, m22, r_c):
    rho_c = Soliton_rho_c_from_r_c(m22, r_c)
    x = r / r_c
    return rho_c*( 1.0+9.06e-2*(x**2) )**(-8)
###################################################################################################################


###################################################################################################################
# parameters for the halo
halo_fitting_rho_0    = 1633769.320 *(Const_Msun/Const_kpc**3)/UNIT_D   # in code_density
halo_fitting_r_0      =       4.461 *(Const_kpc)/UNIT_L                 # in code_length
halo_fitting_alpha    =       4.800
halo_fitting_beta     =       3.098
halo_fitting_gamma    =       0.000

# parameters for the soliton
soliton_fitting_m22   =       1.0
soliton_fitting_r_c   =       0.001                                     # in code_length
soliton_fitting_rho_c = Soliton_rho_c_from_r_c(soliton_fitting_m22, soliton_fitting_r_c)


# parameters for the density profile sampling
halo_r_min            =       0.01  *(Const_kpc)/UNIT_L                 # in code_length
halo_r_max            =     200.0   *(Const_kpc)/UNIT_L                 # in code_length
halo_nbins            =    4096

soliton_r_min         =       0.01  *(Const_kpc)/UNIT_L                 # in code_length
soliton_r_max         =      30.0   *(Const_kpc)/UNIT_L                 # in code_length
soliton_nbins         =    4096
###################################################################################################################


###################################################################################################################
# output the information
print( 'Information'                                                                 )
print( 'UNIT_L                 = {: >16.8e} cm'.format(      UNIT_L                ) )
print( 'UNIT_D                 = {: >16.8e} g/cm**3'.format( UNIT_D                ) )
print( 'halo_fitting_rho_0     = {: >16.8e} UNIT_D'.format(  halo_fitting_rho_0    ) )
print( 'halo_fitting_r_0       = {: >16.8e} UNIT_L'.format(  halo_fitting_r_0      ) )
print( 'halo_fitting_alpha     = {: >16.8e}'.format(         halo_fitting_alpha    ) )
print( 'halo_fitting_beta      = {: >16.8e}'.format(         halo_fitting_beta     ) )
print( 'halo_fitting_gamma     = {: >16.8e}'.format(         halo_fitting_gamma    ) )
print( ''                                                                            )
print( 'soliton_fitting_rho_c  = {: >16.8e} UNIT_D'.format(  soliton_fitting_rho_c ) )
print( 'soliton_fitting_r_c    = {: >16.8e} UNIT_L'.format(  soliton_fitting_r_c   ) )
print( ''                                                                            )
print( 'halo_r_min             = {: >16.8e} UNIT_L'.format(  halo_r_min            ) )
print( 'halo_r_max             = {: >16.8e} UNIT_L'.format(  halo_r_max            ) )
print( 'halo_nbins             = {: >16d}'.format(           halo_nbins            ) )
print( ''                                                                            )
print( 'soliton_r_min          = {: >16.8e} UNIT_L'.format(  soliton_r_min         ) )
print( 'soliton_r_max          = {: >16.8e} UNIT_L'.format(  soliton_r_max         ) )
print( 'soliton_nbins          = {: >16d}'.format(           soliton_nbins         ) )
###################################################################################################################


###################################################################################################################
# create the density profile
halo_densprof_radius     = np.logspace( np.log10(halo_r_min), np.log10(halo_r_max), num=halo_nbins )
halo_densprof_density    = Halo_fitting_analytical_dens( halo_densprof_radius,
                                                         halo_fitting_rho_0, halo_fitting_r_0,
                                                         halo_fitting_alpha, halo_fitting_beta, halo_fitting_gamma )

soliton_densprof_radius  = np.logspace( np.log10(soliton_r_min), np.log10(soliton_r_max), num=soliton_nbins )
soliton_densprof_density = Soliton_fitting_analytical_dens( soliton_densprof_radius,
                                                            soliton_fitting_m22, soliton_fitting_r_c )
###################################################################################################################


###################################################################################################################
# save to file
np.savetxt( 'HaloDensityProfile',
            np.column_stack( (halo_densprof_radius, halo_densprof_density) ),
            fmt='          %9.8e',
            header='                     r                  density' )

np.savetxt( 'SolitonDensityProfile',
            np.column_stack( (soliton_densprof_radius, soliton_densprof_density) ),
            fmt='          %9.8e',
            header='                     r                  density' )
###################################################################################################################


###################################################################################################################
# plot to images
fig = plt.figure()
ax  = fig.add_subplot(111)

# plot some important values for reference
ax.plot( [halo_fitting_r_0,                 halo_fitting_r_0                ], [0.3*np.min(halo_densprof_density), 3.0*np.max(halo_densprof_density)],    '--',  color='grey',  label=r'$r_0$'     )
ax.plot( [0.3*np.min(halo_densprof_radius), 3.0*np.max(halo_densprof_radius)], [halo_fitting_rho_0,                halo_fitting_rho_0               ],    '--',  color='grey',  label=r'$\rho_0$'  )

# plot the density profile
ax.plot( halo_densprof_radius,                                                 halo_densprof_density,                                                     '-',   color='r',     label=r'$\rho(r)$' )

# annotate the information
ax.annotate( r'$\rho_{0}$ = %.8e'%(halo_fitting_rho_0)+'\n'+
             r'$r_{0}$ = %.8e'%(halo_fitting_r_0)+'\n'+
             r'$\alpha$ = %.8e'%(halo_fitting_alpha)+'\n'+
             r'$\beta$ = %.8e'%(halo_fitting_beta)+'\n'+
             r'$\gamma$ = %.8e'%(halo_fitting_gamma)+'\n'+'\n'+
             r'$r_{\rm min}$ = %.8e'%(halo_r_min)+'\n'+
             r'$r_{\rm max}$ = %.8e'%(halo_r_max)+'\n'+
             r'nbins = %d'%(halo_nbins),
             xy=(0.02,0.02), xycoords='axes fraction' )

# setting for the figure
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlim( 0.5*np.min(halo_densprof_radius),  2.0*np.max(halo_densprof_radius)  )
ax.set_ylim( 0.5*np.min(halo_densprof_density), 2.0*np.max(halo_densprof_density) )

# set the labels
ax.set_xlabel( r'$r$'+' (code_length)'     )
ax.set_ylabel( r'$\rho$'+' (code_density)' )
fig.suptitle(   'Density Profile of Halo'  )
ax.legend( loc='upper right' )

# save the figure
fig.subplots_adjust( top=0.93, bottom=0.1, left=0.1, right=0.97 )
fig.savefig( 'fig_HaloDensityProfile.png' )
plt.close()

fig = plt.figure()
ax  = fig.add_subplot(111)

# plot some important values for reference
ax.plot( [soliton_fitting_r_c,              soliton_fitting_r_c                   ], [0.3*np.min(soliton_densprof_density), 3.0*np.max(soliton_densprof_density)],    '--',  color='grey',  label=r'$r_{\rm c}$'     )
ax.plot( [0.3*np.min(soliton_densprof_radius), 3.0*np.max(soliton_densprof_radius)], [soliton_fitting_rho_c,                soliton_fitting_rho_c               ],    '--',  color='grey',  label=r'$\rho_{\rm c}$'  )

# plot the density profile
ax.plot( soliton_densprof_radius,                                                    soliton_densprof_density,                                                        '-',   color='r',     label=r'$\rho(r)$' )

# annotate the information
ax.annotate( r'$\rho_{\rm c}$ = %.8e'%(soliton_fitting_rho_c)+'\n'+
             r'$r_{\rm c}$ = %.8e'%(soliton_fitting_r_c)+'\n'+
             r'$r_{\rm min}$ = %.8e'%(soliton_r_min)+'\n'+
             r'$r_{\rm max}$ = %.8e'%(soliton_r_max)+'\n'+
             r'nbins = %d'%(soliton_nbins),
             xy=(0.02,0.02), xycoords='axes fraction' )

# setting for the figure
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlim(    0.5*np.min(soliton_densprof_radius),  2.0*np.max(soliton_densprof_radius)  )
ax.set_ylim( 0.0001*np.max(soliton_densprof_density), 3.0*np.max(soliton_densprof_density) )

# set the labels
ax.set_xlabel( r'$r$'+' (code_length)'     )
ax.set_ylabel( r'$\rho$'+' (code_density)' )
fig.suptitle(   'Density Profile of Soliton'  )
ax.legend( loc='upper right' )

# save the figure
fig.subplots_adjust( top=0.93, bottom=0.1, left=0.1, right=0.97 )
fig.savefig( 'fig_SolitonDensityProfile.png' )
plt.close()
###################################################################################################################
