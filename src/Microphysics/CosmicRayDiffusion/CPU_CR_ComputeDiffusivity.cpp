#ifndef __CUFLU_CR_COMPUTEDIFFUSIVITY__
#define __CUFLU_CR_COMPUTEDIFFUSIVITY__



#include "CUFLU.h"

#ifdef CR_DIFFUSION




//-----------------------------------------------------------------------------------------
// Function    : CR_ComputeDiffusivity
//
// Description : Compute the diffusion coefficients of cosmic rays
//
// Note        : This is a constant for now. In general, the coefficients should depend on the
//               gas and CR properties.
//
// Parameter   : diff_cr_para : Variables to store the parallel/perpendicular diffusion coefficients
//               diff_cr_perp   --> Call-by-reference
//               MicroPhy     : Microphysics object
//
// Return      : diff_cr_para, diff_cr_perp
//-----------------------------------------------------------------------------------------
GPU_DEVICE
void CR_ComputeDiffusivity( real &diff_cr_para, real &diff_cr_perp, const MicroPhy_t *MicroPhy )
{

   diff_cr_para = MicroPhy->CR_diff_coeff_para;
   diff_cr_perp = MicroPhy->CR_diff_coeff_perp;

} // FUNCTION : CR_ComputeDiffusivity



#endif // #ifdef CR_DIFFUSION



#endif // #ifndef __CUFLU_CR_COMPUTEDIFFUSIVITY__
