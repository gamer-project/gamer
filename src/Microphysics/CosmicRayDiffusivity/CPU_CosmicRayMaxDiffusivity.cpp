#ifndef __COMPUTE_COSMICRAYMAXDIFFUSIVITY__
#define __COMPUTE_COSMICRAYMAXDIFFUSIVITY__

#include "CUFLU.h"

#if ( ( MODEL == HYDRO ) && defined CR_DIFFUSION )

//external functions
#ifdef __CUDACC__

# include "CUFLU_ComputeCosmicRayDiffusivity.cu"

#else // #ifdef __CUDACC__

void CR_ComputeDiffusivity( /*dens*/ real &diff_cr_para, real &diff_cr_perp, const MicroPhy_t *Mic );

#endif // #ifdef __CUDACC__ ... else ...


//-----------------------------------------------------------------------------------------
// Function    : Hydro_CosmicRayMaxDiffusivity
// Description : Compute the max diffusivity of cosmic ray
// Note        : This is a constant for now. The coefficient should depend on cell fluid.
// Parameter   :  
// Return      : The max diffusion coefficient.
//-----------------------------------------------------------------------------------------
GPU_DEVICE
void CR_MaxDiffusivity( real &max_diff, const MicroPhy_t *Mic )
{
   real diff_cr_para, diff_cr_perp;
   CR_ComputeDiffusivity( diff_cr_para, diff_cr_perp, Mic );
   max_diff = MAX(diff_cr_para, diff_cr_perp);

} // FUNCTION : CR_MaxDiffusivity

#endif // #if ( ( MODEL == HYDRO ) && defined CR_DIFFUSION )

#endif // #ifndef __COMPUTE_COSMICRAYMAXDIFFUSIVITY__
