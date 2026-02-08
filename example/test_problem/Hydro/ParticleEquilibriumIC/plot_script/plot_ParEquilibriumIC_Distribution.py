import matplotlib.pyplot as plt
import numpy as np

# Model
MODEL = 1 # 0=Table, 1=Plummer, 2=NFW, 3=Burkert, 4=Jaffe, 5=Hernquist, 6=Einasto
ID    = 1


# Parameters for analytical model
G    = 1.0
R0   = 1.0
Rho0 = 1.0/(2*np.pi)


# Analytical solutoins for Hernquist model
def Density_Hernquist( r ):
    return Rho0/( ( r/R0 )*( 1.0+(r/R0) )**3 )

def EnclosedMass_Hernquist( r ):
    M0 = 2*np.pi*Rho0*R0**3
    return M0*( r**2 )/( r+R0 )**2

def Potential_Hernquist( r ):
    M0 = 2*np.pi*Rho0*R0**3
    return ( -G*M0/R0 )/( 1.0+(r/R0) )

def DistributionFunction_Hernquist( E ):
    # The analytical solution is given in Eq. (4.51) of Galactic Dynamics: Second Edition, Binney & Tremaine 2008.
    # Note that we multiply it by M0 here because the distribution function in the textbook is normalized to the total probability: $\int \int f(E) d^{3}v d^{3}x = 1$, while in our calculation, it is normalized to the total mass.
    M0   = 2*np.pi*Rho0*R0**3
    Etil = E*R0/(G*M0)

    return M0 * (2**(-0.5))*(2*np.pi*(G*M0*R0)**0.5)**(-3)*((Etil**0.5)/(1-Etil)**2)*((1-2*Etil)*(8*Etil**2-8*Etil-3)+(3*np.arcsin(Etil**0.5))/(Etil*(1-Etil))**0.5)

# Analytical solutoins for Plummer model
def Density_Plummer( r ):
    return Rho0/( 1.0+(r/R0)**2 )**(2.5)

def EnclosedMass_Plummer( r ):
    M0 = 4.0/3.0*np.pi*Rho0*R0**3
    return M0*( r**3 )/( r**2 + R0**2 )**(1.5)

def Potential_Plummer( r ):
    M0 = 4.0/3.0*np.pi*Rho0*R0**3
    return -G*M0/( r**2 + R0**2 )**0.5

def DistributionFunction_Plummer( E ):
    M0 = 4.0/3.0*np.pi*Rho0*R0**3
    return (24.0*np.sqrt(2.0))/(7.0*np.pi**3)*(R0**2/((G**5)*(M0**4)))*(E**3.5)


# Analytical solutoins
def Density_Analytical( r ):
    if MODEL == 5:
        return Density_Hernquist( r )
    elif MODEL == 1:
        return Density_Plummer( r )
    else:
        return np.full( np.shape(r), np.nan )

def EnclosedMass_Analytical( r ):
    if MODEL == 5:
        return EnclosedMass_Hernquist( r )
    elif MODEL == 1:
        return EnclosedMass_Plummer( r )
    else:
        return np.full( np.shape(r), np.nan )

def Potential_Analytical( r ):
    if MODEL == 5:
        return Potential_Hernquist( r )
    elif MODEL == 1:
        return Potential_Plummer( r )
    else:
        return np.full( np.shape(r), np.nan )

def DistributionFunction_Analytical( E ):
    if MODEL == 5:
        return DistributionFunction_Hernquist( E )
    elif MODEL == 1:
        return DistributionFunction_Plummer( E )
    else:
        return np.full( np.shape(E), np.nan )


# Load the data
R, Rho, M_Enc, Phi, dRho_dPsi = np.loadtxt( '../Record__ParEquilibriumIC_Model_%d_RArray_%d'%(MODEL,ID), unpack=True )
E, IntDFunc, DFunc            = np.loadtxt( '../Record__ParEquilibriumIC_Model_%d_EArray_%d'%(MODEL,ID), unpack=True )


# Sample points of analytical solutions
R_anal = np.logspace( np.log10(R[1])-2, np.log10(R[-1])+2, num=1000 )
E_anal = np.logspace( np.log10(E[1])-2, np.log10(E[-1])+2, num=1000 )


# Plot the data
fig = plt.figure( figsize=(24,12) )

ax1 = fig.add_subplot( 5,2,1      )
ax2 = fig.add_subplot( 5,2,2      )
ax3 = fig.add_subplot( 5,2,3      )
ax4 = fig.add_subplot( 5,2,4      )
ax5 = fig.add_subplot( 5,2,(5,6)  )
ax6 = fig.add_subplot( 5,2,(7,8)  )
ax7 = fig.add_subplot( 5,2,(9,10) )

ax1.plot( R,                               Rho,                                                                         '.-', color='b', label='Numerical'  )
ax1.plot( R_anal,                          Density_Analytical( R_anal ),                                                '--', color='r', label='Analytical' )

ax2.plot( R,                               M_Enc,                                                                       '.-', color='b', label='Numerical'  )
ax2.plot( R_anal,                          EnclosedMass_Analytical( R_anal ),                                           '--', color='r', label='Analytical' )

ax3.plot( R,                               Phi,                                                                         '.-', color='b', label='Numerical'  )
ax3.plot( R_anal,                          Potential_Analytical( R_anal ),                                              '--', color='r', label='Analytical' )

ax4.plot( R,                               dRho_dPsi,                                                                   '.-', color='b', label='Numerical'  )
ax4.plot( R_anal,                          np.gradient(Density_Analytical( R_anal ), -Potential_Analytical( R_anal )),  '--', color='r', label='Analytical' )

ax5.plot( -Phi,                            dRho_dPsi,                                                                   '.-', color='g', label='Numerical'  )
ax5.plot( -Potential_Analytical( R_anal ), np.gradient(Density_Analytical( R_anal ), -Potential_Analytical( R_anal )),  '--', color='r', label='Analytical' )

ax6.plot( E,                               IntDFunc,                                                                    '.-', color='g', label='Numerical'  )

ax7.plot( E,                               DFunc,                                                                       '.-', color='g', label='Numerical'  )
ax7.plot( E_anal,                          DistributionFunction_Analytical( E_anal ),                                   '--', color='r', label='Analytical' )

ax1.set_title( r'$\rho(r)$'             )
ax2.set_title( r'$M_{\rm enc}(r)$'      )
ax3.set_title( r'$\Phi(r)$'             )
ax4.set_title( r'$\frac{d\rho}{d\Psi}$' )
ax5.set_title( r'$\frac{d\rho}{d\Psi}$' )
ax6.set_title( r'Integrated Distribution Function, $\int_0^{\mathcal{E}}\frac{1}{\sqrt{\mathcal{E}-\Psi}}\frac{d\rho}{d\Psi}d\Psi$' )
ax7.set_title( r'Distribution Function, $f(\mathcal{E})$' )

ax1.set_xlabel( r'$r$'           )
ax2.set_xlabel( r'$r$'           )
ax3.set_xlabel( r'$r$'           )
ax4.set_xlabel( r'$r$'           )
ax5.set_xlabel( r'$\Psi$'        )
ax6.set_xlabel( r'$\mathcal{E}$' )
ax7.set_xlabel( r'$\mathcal{E}$' )

ax1.set_xscale( 'symlog', linthresh=1e+0*R[1]                       )
ax2.set_xscale( 'symlog', linthresh=1e+0*R[1]                       )
ax3.set_xscale( 'symlog', linthresh=1e+0*R[1]                       )
ax4.set_xscale( 'symlog', linthresh=1e+0*R[1]                       )
ax5.set_xscale( 'symlog', linthresh=1e+0*E[1]                       )
ax6.set_xscale( 'symlog', linthresh=1e+0*E[1]                       )
ax7.set_xscale( 'symlog', linthresh=1e+0*E[1]                       )
ax1.set_yscale( 'symlog', linthresh=1e-3*Rho[Rho>0][-1]             )
ax2.set_yscale( 'symlog', linthresh=1e-3*M_Enc[1]                   )
ax3.set_yscale( 'symlog', linthresh=1e-3*-Phi[Phi<0][-1]            )
ax4.set_yscale( 'symlog', linthresh=1e-3*dRho_dPsi[dRho_dPsi>0][-1] )
ax5.set_yscale( 'symlog', linthresh=1e-3*dRho_dPsi[dRho_dPsi>0][-1] )
ax6.set_yscale( 'symlog', linthresh=1e-3*IntDFunc[IntDFunc>0][0]    )
ax7.set_yscale( 'symlog', linthresh=1e-3*DFunc[DFunc>0][0]          )

ax1.set_xlim( -2e-2*R[1], 2e+1*R[-1] )
ax2.set_xlim( -2e-2*R[1], 2e+1*R[-1] )
ax3.set_xlim( -2e-2*R[1], 2e+1*R[-1] )
ax4.set_xlim( -2e-2*R[1], 2e+1*R[-1] )
ax5.set_xlim( -2e-2*E[1], 2e+0*E[-1] )
ax6.set_xlim( -2e-2*E[1], 2e+0*E[-1] )
ax7.set_xlim( -2e-2*E[1], 2e+0*E[-1] )

ax1.legend()
ax2.legend()
ax3.legend()
ax4.legend()
ax5.legend()
ax6.legend()
ax7.legend()


# Save to file
plt.tight_layout()
fig.savefig( 'fig_ParEquilibriumIC_Distribution_Model_%d_ID_%d.png'%(MODEL,ID) )
