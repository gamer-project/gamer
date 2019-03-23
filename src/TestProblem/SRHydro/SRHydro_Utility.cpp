#include <math.h>
#include "../../../include/Macro.h"

#if ( MODEL == SR_HYDRO )

void 
SRHydro_4Velto3Vel_Double (const double In[], double Out[])
{
  double Factor = 1 / SQRT (1 + SQR (In[1]) + SQR (In[2]) + SQR (In[3]));

  Out[0] = In[0];
  Out[1] = In[1] * Factor;
  Out[2] = In[2] * Factor;
  Out[3] = In[3] * Factor;
  Out[4] = In[4];
}				// FUNCTION : SRHydro_4Velto3Vel

void 
SRHydro_3Velto4Vel_Double (const double In[], double Out[])
{
  double Factor = 1 / SQRT (1 - SQR (In[1]) - SQR (In[2]) - SQR (In[3]));

  Out[0] = In[0];
  Out[1] = In[1] * Factor;
  Out[2] = In[2] * Factor;
  Out[3] = In[3] * Factor;
  Out[4] = In[4];
}
				// FUNCTION : SRHydro_4Velto3Vel
void
SRHydro_Pri2Con_Double (const double In[], double Out[], const double Gamma)
{
# if ( EOS == RELATIVISTIC_IDEAL_GAS )
  double nh = 2.5*In[4] + SQRT(2.25*SQR(In[4]) + SQR(In[0])); // approximate enthalpy * proper number density
# elif ( EOS == IDEAL_GAS )
  double Gamma_m1 = Gamma - 1.0;
  double nh = In[0] + ( Gamma / Gamma_m1) * In[4]; // enthalpy * proper number density
# else
# error: unsupported EoS!
# endif

  double Factor0 = 1.0 + SQR (In[1]) + SQR (In[2]) + SQR (In[3]);
  double Factor1 = SQRT(Factor0); // Lorentz factor
  double Factor2 = nh * Factor1;
  
  Out[0] = In[0] * Factor1; // number density in inertial frame
  Out[1] = Factor2 * In[1]; // MomX
  Out[2] = Factor2 * In[2]; // MomX
  Out[3] = Factor2 * In[3]; // MomX
# if   ( CONSERVED_ENERGY == 1 )
  Out[4] = nh * Factor0 - In[4]; // total_energy
# elif ( CONSERVED_ENERGY == 2 )
  Out[4] = nh * Factor0 - In[4] - Out[0]; // ( total_energy ) - ( rest_mass_energy )
# else
# error: CONSERVED_ENERGY must be 1 or 2!
# endif
}				// FUNCTION : SRHydro_Pri2Con_Double

#endif
