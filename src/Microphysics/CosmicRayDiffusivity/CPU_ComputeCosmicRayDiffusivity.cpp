#ifndef __COMPUTE_COSMICRAYDIFFUSIVITY__
#define __COMPUTE_COSMICRAYDIFFUSIVITY__

#include "CUFLU.h"

#ifdef CR_DIFFUSION

//-----------------------------------------------------------------------------------------
// Function    : CR_ComputeDiffusivity
//
// Description : Compute the diffusion coefficients of cosmic ray.
//
// Note        : This is a constant for now. The coefficient should depend on the cell condition.
//
// Parameter   : diff_cr_para : Variable to store the parallel      diffusion coefficient
//               diff_cr_perp : Variable to store the perpendicular diffusion coefficient
//               Mic          : Microphysics object
//-----------------------------------------------------------------------------------------
GPU_DEVICE
void CR_ComputeDiffusivity( real &diff_cr_para, real &diff_cr_perp, const MicroPhy_t *Mic )
{
   diff_cr_para = Mic->CR_diff_coeff_para;
   diff_cr_perp = Mic->CR_diff_coeff_perp;

} // FUNCTION : CR_ComputeDiffusivity

#endif // #ifdef CR_DIFFUSION

#endif // #ifndef __COMPUTE_COSMICRAYDIFFUSIVITY__
