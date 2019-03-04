#ifndef __CPU_CONDUCTIVE_FLUXES__
#define __CPU_CONDUCTIVE_FLUXES__

#include "Microphysics.h"

#if ( ( MODEL == HYDRO ) && defined CONDUCTION )

// external functions
#ifdef __CUDACC__
# include "CUFLU_Conduction.cu"
#else // #ifdef __CUDACC__

#endif // #ifdef __CUDACC__ ... else ...

GPU_DEVICE
void Hydro_ComputeConductiveFluxes( const real g_FC_Var [][NCOMP_TOTAL_PLUS_MAG][ CUBE(N_FC_VAR) ],
                                    const real Flux_1Face[NCOMP_TOTAL_PLUS_MAG],
                                    const int i_flux, const int j_flux, const int k_flux,
                                    const int i_fc, const int j_fc, const int k_fc, 
                                    const int d, const real dt, const real dh, const double Time )
{

    real vx, vy, vz, StressMomX, StressMomY, StressMomZ;
    const real one_third = 1./3.;

    chi = Hydro_ComputeConduction( fluid, Gamma_m1, MinPres );

    if ( CONDUCTIVE_FLUX_TYPE == ANISOTROPIC ) {

#ifdef MHD

//      Anisotropic conduction
        switch ( d ) {
            case 0:
                break;
            case 1:
                break;
            case 2:
                break; 
        }

#endif // #ifdef MHD

    } else if ( CONDUCTIVE_FLUX_TYPE == ISOTROPIC ) {

//      Isotropic conduction

        switch ( d ) {
            case 0:
                break;
            case 1:
                break;
            case 2:
                break; 
        }

    } else {

        Aux_Error();

    }
   
    Flux_1Face[4] += vx*StressMomX + vy*StressMomY + vz*StressMomZ;

#  ifdef __CUDACC__
    __syncthreads();
#  endif


} // FUNCTION : Hydro_ComputeConductiveFluxes

#endif // #if ( ( MODEL == HYDRO ) && defined VISCOSITY )

#endif // #ifndef __CPU_CONDUCTIVE_FLUXES__
