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
//
// Return      :  fluid
//-----------------------------------------------------------------------------------------
GPU_DEVICE
void Src_Deleptonization( real fluid[], const real B[],
                          const SrcTerms_t SrcTerms, const real dt, const real dh,
                          const double x, const double y, const double z,
                          const double TimeNew, const double TimeOld,
                          const real MinDens, const real MinPres, const real MinEint )
{

// TBF

} // FUNCTION : Src_Deleptonization



#endif // #ifndef __SRC_DELEPTONIZATION__
