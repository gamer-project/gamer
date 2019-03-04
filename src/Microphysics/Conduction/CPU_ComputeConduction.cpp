#ifndef __COMPUTE_CONDUCTION__
#define __COMPUTE_CONDUCTION__

#include "Microphysics.h"

#if ( ( MODEL == HYDRO ) && defined CONDUCTION )

GPU_DEVICE
real Hydro_ComputeConduction( const real fluid[NCOMP_FLUID],
                              const real Gamma_m1, const real MinPres )
{
    const bool CheckMinPres_Yes = true;

    real chi, _Rho;

    _Rho  = (real)1.0 / fluid[DENS];

    if ( VISCOSITY_TYPE == CONSTANT_CONDUCTION ) {

        // Constant conduction
        if ( CONDUCTION_COEFF_TYPE == CONDUCTION_KINETIC_COEFF ) {

            nu = (real)CONDUCTION_COEFF;

        } else if ( VISCOSITY_COEFF_TYPE == CONDUCTION_DYNAMIC_COEFF ) {

            nu = (real)CONDUCTION_COEFF*_Rho;

        }

    } else if ( CONDUCTION_TYPE == SPITZER_CONDUCTION ) { 

        // Spitzer conduction
        real Pres, Temp, Freq_ii;

        Pres = Hydro_GetPressure( fluid[DENS], fluid[MOMX], fluid[MOMY], fluid[MOMZ], 
                                  fluid[ENGY], Gamma_m1, CheckMinPres_Yes, MinPres );
        
        Temp = Pres*_Rho;

        Freq_ii = FreqPrefactor*fluid[DENS]*POW( Temp, (real)-1.5 );

        chi = CONDUCTION_SPITZER_FRACTION*0.96*Pres/Freq_ii;

    }

    chi = FMIN( FMAX( chi, CONDUCTION_COEFF_MIN ), CONDUCTION_COEFF_MAX );

    return nu;

} // FUNCTION : Hydro_ComputeConduction

#endif // #if ( ( MODEL == HYDRO ) && defined CONDUCTION )

#endif // #ifndef __COMPUTE_CONDUCTION__

