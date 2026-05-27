#ifndef __CUFLU_COMPUTEVISCOSITY__
#define __CUFLU_COMPUTEVISCOSITY__

#include "CUFLU.h"

#if ( ( MODEL == HYDRO ) && defined VISCOSITY )

GPU_DEVICE
void Hydro_ComputeViscosity( real &visc_mu, real &visc_nu, const MicroPhy_t *MicroPhy,
                             const real Dens, const real Temp )
{

    real _Rho  = (real)1.0 / Dens;

    if ( MicroPhy->ViscType == CONSTANT_VISCOSITY ) {

        // Constant viscosity

        if ( MicroPhy->ViscCoeffType == VISCOSITY_KINETIC_COEFF ) {

            // nu is constant

            visc_mu = MicroPhy->ViscConstCoeff*Dens;

        } else if ( MicroPhy->ViscCoeffType == VISCOSITY_DYNAMIC_COEFF ) {

            // mu is constant

            visc_mu = MicroPhy->ViscConstCoeff;
        }

    } else if ( MicroPhy->ViscType == SPITZER_VISCOSITY ) {

        // Spitzer viscosity, dependent on T

        visc_mu = MicroPhy->ViscPrefactor*POW( (real)1.0e-7*Temp, (real)2.5 );

    }

    visc_mu = FMIN( visc_mu, MicroPhy->ViscMaxDiffusivity*Dens );
    visc_nu = visc_mu*_Rho;

    return;

} // FUNCTION : Hydro_ComputeViscosity

#endif // #if ( ( MODEL == HYDRO ) && defined VISCOSITY )

#endif // #ifndef __CUFLU_COMPUTEVISCOSITY__

