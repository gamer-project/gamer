#ifndef __CUFLU_COMPUTECONDUCTION__
#define __CUFLU_COMPUTECONDUCTION__

#include "CUFLU.h"

#if ( ( MODEL == HYDRO ) && defined CONDUCTION )

GPU_DEVICE
void Hydro_ComputeConduction( real &cond_kappa, real &cond_chi, const MicroPhy_t *MicroPhy,
                              const real Dens, const real Temp )
{

    real _Rho  = (real)1.0 / Dens;

    if ( MicroPhy->CondType == CONSTANT_CONDUCTION )
       // Constant conductivity
       cond_kappa = MicroPhy->CondConstCoeff;
    else if ( MicroPhy->CondType == SPITZER_CONDUCTION )
       // Spitzer conductivity, dependent on T
       cond_kappa = MicroPhy->CondPrefactor*POW( (real)1.0e-7*Temp, (real)2.5 );
    cond_kappa = FMIN( cond_kappa, MicroPhy->CondMaxDiffusivity * Dens * MicroPhy->CondSpecificHeat );
    cond_chi = cond_kappa*_Rho / MicroPhy->CondSpecificHeat;

    return;

} // FUNCTION : Hydro_ComputeConduction

#endif // #if ( ( MODEL == HYDRO ) && defined CONDUCTION )

#endif // #ifndef __CUFLU_COMPUTECONDUCTION__

