import numpy as np
import matplotlib.pyplot as plt


###################################################################################################################
# code units (in cgs)
UNIT_L     = 4.436632034507549e+24
UNIT_D     = 2.5758579724476994e-30

# length units
Const_cm   = 1.0;
Const_m    = 1.0e2*Const_cm;
Const_km   = 1.0e5*Const_cm;
Const_pc   = 3.08567758149e18;
Const_kpc  = 1.0e3*Const_pc;
Const_Mpc  = 1.0e6*Const_pc;
Const_Gpc  = 1.0e9*Const_pc;
Const_au   = 1.495978707e13;
Const_ly   = 9.46053e17;

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
###################################################################################################################


###################################################################################################################
# parameters for the halo
fitting_rho_0 = 1976031.787 *(Const_Msun/Const_kpc**3)/UNIT_D   # in code_density
fitting_r_0   =       6.081 *(Const_kpc)/UNIT_L                 # in code_length
fitting_alpha =       1.815
fitting_beta  =       4.002
fitting_gamma =       0.000

# parameters for the density profile sampling
r_min         =       0.1   *(Const_kpc)/UNIT_L                 # in code_length
r_max         =     200.0   *(Const_kpc)/UNIT_L                 # in code_length
nbins         =    4096
###################################################################################################################


###################################################################################################################
# output the information
print( 'Information'                                                   )
print( 'UNIT_L            = {: >16.8e} cm'.format(      UNIT_L       ) )
print( 'UNIT_D            = {: >16.8e} g/cm**3'.format( UNIT_D       ) )
print( 'fitting_rho_0     = {: >16.8e} UNIT_D'.format( fitting_rho_0 ) )
print( 'fitting_r_0       = {: >16.8e} UNIT_L'.format( fitting_r_0   ) )
print( 'fitting_alpha     = {: >16.8e}'.format(        fitting_alpha ) )
print( 'fitting_beta      = {: >16.8e}'.format(        fitting_beta  ) )
print( 'fitting_gamma     = {: >16.8e}'.format(        fitting_gamma ) )
print( ''                                                              )
print( 'r_min             = {: >16.8e} UNIT_L'.format( r_min         ) )
print( 'r_max             = {: >16.8e} UNIT_L'.format( r_max         ) )
print( 'nbins             = {: >16d}'.format(          nbins         ) )
###################################################################################################################


###################################################################################################################
# create the density profile
halo_densprof_radius  = np.logspace( np.log10(r_min), np.log10(r_max), num=nbins )
halo_densprof_density = Halo_fitting_analytical_dens( halo_densprof_radius,
                                                      fitting_rho_0, fitting_r_0,
                                                      fitting_alpha, fitting_beta, fitting_gamma )
###################################################################################################################


###################################################################################################################
# save to file
np.savetxt( 'HaloDensityProfile',
            np.column_stack( (halo_densprof_radius, halo_densprof_density) ),
            fmt='          %9.8e',
            header='                     r                  density' )
###################################################################################################################


###################################################################################################################
# plot to images
fig = plt.figure()
ax  = fig.add_subplot(111)

# plot some important values for reference
ax.plot( [fitting_r_0,                      fitting_r_0                     ], [0.3*np.min(halo_densprof_density), 3.0*np.max(halo_densprof_density)],    '--',  color='grey',  label=r'$r_0$'     )
ax.plot( [0.3*np.min(halo_densprof_radius), 3.0*np.max(halo_densprof_radius)], [fitting_rho_0,                     fitting_rho_0                    ],    '--',  color='grey',  label=r'$\rho_0$'  )

# plot the density profile
ax.plot( halo_densprof_radius,                                                 halo_densprof_density,                                                     '-',   color='r',     label=r'$\rho(r)$' )

# annotate the information
ax.annotate( r'$\rho_{0}$ = %.8e'%(fitting_rho_0)+'\n'+
             r'$r_{0}$ = %.8e'%(fitting_r_0)+'\n'+
             r'$\alpha$ = %.8e'%(fitting_alpha)+'\n'+
             r'$\beta$ = %.8e'%(fitting_beta)+'\n'+
             r'$\gamma$ = %.8e'%(fitting_gamma)+'\n'+'\n'+
             r'$r_{\rm min}$ = %.8e'%(r_min)+'\n'+
             r'$r_{\rm max}$ = %.8e'%(r_max)+'\n'+
             r'nbins = %d'%(nbins),
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
###################################################################################################################
