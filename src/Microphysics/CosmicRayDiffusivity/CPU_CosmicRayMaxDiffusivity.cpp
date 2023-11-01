#ifndef __COMPUTE_COSMICRAYMAXDIFFUSIVITY__
#define __COMPUTE_COSMICRAYMAXDIFFUSIVITY__

#include "CUFLU.h"

#ifdef CR_DIFFUSION

//external functions
#ifdef __CUDACC__

# include "CUFLU_ComputeCosmicRayDiffusivity.cu"

#else // #ifdef __CUDACC__

void CR_ComputeDiffusivity( real &diff_cr_para, real &diff_cr_perp, const MicroPhy_t *Mic );

#endif // #ifdef __CUDACC__ ... else ...


//-----------------------------------------------------------------------------------------
// Function    : Hydro_CosmicRayMaxDiffusivity
//
// Description : Compute the maximum diffusion coefficient of cosmic ray
//
// Note        :
//
// Parameter   : max_diff : Variable to store the maximum diffusion coefficient.
//               Mic      : Microphysics object
//-----------------------------------------------------------------------------------------
GPU_DEVICE
void CR_MaxDiffusivity( real &max_diff, const MicroPhy_t *Mic )
{
   real diff_cr_para, diff_cr_perp;
   CR_ComputeDiffusivity( diff_cr_para, diff_cr_perp, Mic );
   max_diff = MAX(diff_cr_para, diff_cr_perp);

} // FUNCTION : CR_MaxDiffusivity

#endif // #ifdef CR_DIFFUSION

#endif // #ifndef __COMPUTE_COSMICRAYMAXDIFFUSIVITY__
