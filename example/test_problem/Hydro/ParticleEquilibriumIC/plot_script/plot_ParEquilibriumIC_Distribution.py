import matplotlib.pyplot as plt
import numpy as np

# Parameters for Hernquist model
G    = 1.0
R0   = 1.0
Rho0 = 1.0/(2*np.pi)
M0   = 2*np.pi*Rho0*R0**3


# Analytical solutoins for Hernquist model
def Density_Hernquist( r ):
    return Rho0/( ( r/R0 )*( 1.0+(r/R0) )**3 )

def EnclosedMass_Hernquist( r ):
    return M0*( r**2 )/( r+R0 )**2

def Potential_Hernquist( r ):
    return ( -G*M0/R0 )/( 1.0+(r/R0) )

def DistributionFunction_Hernquist( E ):

    Etil = E*R0/(G*M0)

    return (2**(-0.5))*(2*np.pi*(G*M0*R0)**0.5)**(-3)*((Etil**0.5)/(1-Etil)**2)*((1-2*Etil)*(8*Etil**2-8*Etil-3)+(3*np.arcsin(Etil**0.5))/(Etil*(1-Etil))**0.5)

# Load the data
R, Rho, M_Enc, Phi, dRho_dPsi = np.loadtxt( '../Record__ParEquilibriumIC_Model_5_RArray', unpack=True )
E, IntDFunc, DFunc            = np.loadtxt( '../Record__ParEquilibriumIC_Model_5_EArray', unpack=True )

# Plot the data
fig = plt.figure( figsize=(24,12) )

ax1 = fig.add_subplot( 321 )
ax1.plot( R, Rho,                                 '.-', color='b', label='Numerical'  )
ax1.plot( R, Density_Hernquist( R ),              '--', color='r', label='Analytical' )
ax1.set_title( r'$\rho(r)$' )
ax1.set_xlabel( 'r' )
ax1.set_xscale( 'symlog', linthresh=1e-2 )
ax1.set_yscale( 'symlog', linthresh=1e-9 )
ax1.set_xlim( left=-2e-4 )
ax1.legend()

ax2 = fig.add_subplot( 322 )
ax2.plot( R, M_Enc,                               '.-', color='b', label='Numerical'  )
ax2.plot( R, EnclosedMass_Hernquist( R ),         '--', color='r', label='Analytical' )
ax2.set_title( r'$M_{\rm enc}(r)$' )
ax2.set_xlabel( 'r' )
ax2.set_xscale( 'symlog', linthresh=1e-2 )
ax2.set_yscale( 'symlog', linthresh=1e-9 )
ax2.set_xlim( left=-2e-4 )
ax2.legend()

ax3 = fig.add_subplot( 323 )
ax3.plot( R, Phi,                                 '.-', color='b', label='Numerical'  )
ax3.plot( R, Potential_Hernquist( R ),            '--', color='r', label='Analytical' )
ax3.set_title( r'$\Phi(r)$' )
ax3.set_xlabel( 'r' )
ax3.set_xscale( 'symlog', linthresh=1e-2 )
ax3.set_yscale( 'symlog', linthresh=1e-9 )
ax3.set_xlim( left=-2e-4 )
ax3.legend()

ax4 = fig.add_subplot( 324 )
ax4.plot( R, dRho_dPsi,                           '.-', color='b', label='Numerical'  )
ax4.set_title( r'$\frac{d\rho}{d\Psi}$' )
ax4.set_xlabel( 'r' )
ax4.set_xscale( 'symlog', linthresh=1e-2 )
ax4.set_yscale( 'symlog', linthresh=1e-9 )
ax4.set_xlim( left=-2e-4 )
ax4.set_ylim( bottom=1e-5 )
ax4.legend()

ax5 = fig.add_subplot( 325 )
ax5.plot( E, IntDFunc,                            '.-', color='g', label='Numerical'  )
ax5.set_title( 'Integrated Distribution Function' )
ax5.set_xlabel( r'$\mathcal{E}$' )
ax5.set_xscale( 'symlog', linthresh=1e-3 )
ax5.set_yscale( 'symlog', linthresh=1e-9 )
ax5.set_xlim( left=-2e-5 )
ax5.legend()

ax6 = fig.add_subplot( 326 )
ax6.plot( E, DFunc,                               '.-', color='g', label='Numerical'  )
ax6.plot( E, DistributionFunction_Hernquist( E ), '--', color='r', label='Analytical' )
ax6.set_title( 'Distribution Function, $f(\mathcal{E})$' )
ax6.set_xlabel( r'$\mathcal{E}$' )
ax6.set_xscale( 'symlog', linthresh=1e-3 )
ax6.set_yscale( 'symlog', linthresh=1e-9 )
ax6.set_xlim( left=-2e-5 )
ax6.legend()

# Save to file
plt.tight_layout()
fig.savefig( 'fig_ParEquilibriumIC_Distribution.png' )
