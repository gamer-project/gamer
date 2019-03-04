#ifndef __CUFLU_COMPUTEVISCOSITY__
#define __CUFLU_COMPUTEVISCOSITY__

#include "CUFLU.h"

#if ( ( MODEL == HYDRO ) && defined VISCOSITY )

GPU_DEVICE
void Hydro_ComputeViscosity( real &visc_mu, real &visc_nu, const MicroPhy_t *MicroPhy, 
                             const real Dens, const real Temp )
{

    real _Rho  = (real)1.0 / Dens;

    if ( VISCOSITY_TYPE == CONSTANT_VISCOSITY ) {

        // Constant viscosity
        if ( MicroPhy->visc_coeff_type == VISCOSITY_KINETIC_COEFF ) {

            visc_mu = MicroPhy->ViscKineticCoeff*Dens;

        } else if ( MicroPhy->visc_coeff_type == VISCOSITY_DYNAMIC_COEFF ) {

            visc_mu = MicroPhy->ViscDynamicCoeff;
        }

    } else if ( VISCOSITY_TYPE == SPITZER_VISCOSITY ) { 

        // Spitzer viscosity
    
        visc_mu = MicroPhy->ViscSpitzerFraction*MicroPhy->ViscPrefactor*POW( Temp, (real)2.5 );

    }

    visc_mu = FMIN( visc_mu, MicroPhy->ViscMaxDiffusivity*Dens );
    visc_nu = visc_mu*_Rho;

    return;

} // FUNCTION : Hydro_ComputeViscosity

#endif // #if ( ( MODEL == HYDRO ) && defined VISCOSITY )

#endif // #ifndef __CUFLU_COMPUTEVISCOSITY__

