#ifndef __SRC_DELEPTONIZATION__
#define __SRC_DELEPTONIZATION__



#include "CUFLU.h"



// external functions and GPU-related set-up
#ifdef __CUDACC__

#include "CUFLU_Shared_FluUtility.cu"
#include "CUDA_ConstMemory.h"

#endif // #ifdef __CUDACC__




//-----------------------------------------------------------------------------------------
// Function    :  Src_Deleptonization
// Description :  Deleptonization for neutron star simulations
//
// Note        :  1. Invoked by CPU/GPU_SrcSolver_IterateAllCells()
//                2. Enabled by the runtime option "SRC_DELEPTONIZATION"
//
// Parameter   :  fluid             : Fluid array storing both the input and updated values
//                                    --> Including both active and passive variables
//                B                 : Cell-centered magnetic field
//                SrcTerms          : Structure storing all source-term variables
//                dt                : Time interval to advance solution
//                dh                : Grid size
//                x/y/z             : Target physical coordinates
//                TimeNew           : Target physical time to reach
//                TimeOld           : Physical time before update
//                                    --> This function updates physical time from TimeOld to TimeNew
//                MinDens/Pres/Eint : Density, pressure, and internal energy floors
//                EoS_DensEint2Pres : EoS routine to compute the gas pressure
//                EoS_DensPres2Eint : EoS routine to compute the gas internal energy
//                EoS_DensPres2CSqr : EoS routine to compute the sound speed square
//                EoS_AuxArray_*    : Auxiliary arrays for the EoS routines
//                EoS_Table         : EoS tables
//
// Return      :  fluid[]
//-----------------------------------------------------------------------------------------
GPU_DEVICE
void Src_Deleptonization( real fluid[], const real B[],
                          const SrcTerms_t SrcTerms, const real dt, const real dh,
                          const double x, const double y, const double z,
                          const double TimeNew, const double TimeOld,
                          const real MinDens, const real MinPres, const real MinEint,
                          const EoS_DE2P_t EoS_DensEint2Pres,
                          const EoS_DP2E_t EoS_DensPres2Eint,
                          const EoS_DP2C_t EoS_DensPres2CSqr,
                          const double EoS_AuxArray_Flt[],
                          const int    EoS_AuxArray_Int[],
                          const real *const EoS_Table[EOS_NTABLE_MAX] )
{

// TBF

} // FUNCTION : Src_Deleptonization



#endif // #ifndef __SRC_DELEPTONIZATION__
