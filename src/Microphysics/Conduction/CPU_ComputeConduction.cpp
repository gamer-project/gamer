#ifndef __CUFLU_COMPUTECONDUCTION__
#define __CUFLU_COMPUTECONDUCTION__

#include "CUFLU.h"

#if ( ( MODEL == HYDRO ) && defined CONDUCTION )

GPU_DEVICE
void Hydro_ComputeConduction( real &conda_kappa, real &cond_chi, const MicroPhy_t *MicroPhy, 
                              const real Dens, const real Temp )
{

    real _Rho  = (real)1.0 / Dens;

    if ( CONDUCTION_TYPE == CONSTANT_CONUDCTIVITY ) 
        // Constant conductivity
        cond_kappa = MicroPhy->CondConstCoeff;
    else if ( VISCOSITY_TYPE == SPITZER_VISCOSITY ) 
        // Spitzer conductivity, dependent on T
        cond_kappa = MicroPhy->CondPrefactor*POW( Temp, (real)2.5 );

    cond_kappa = FMIN( cond_kappa, MicroPhy->CondMaxDiffusivity*Dens );
    cond_chi = cond_kappa*_Rho;

    return;

} // FUNCTION : Hydro_ComputeConduction

#endif // #if ( ( MODEL == HYDRO ) && defined CONDUCTION )

#endif // #ifndef __CUFLU_COMPUTECONDUCTION__

