#ifndef __COMPUTE_COSMICRAYDIFFUSIVITY__
#define __COMPUTE_COSMICRAYDIFFUSIVITY__

#include "CUFLU.h"

#if ( ( MODEL == HYDRO ) && defined CR_DIFFUSION )

//-----------------------------------------------------------------------------------------
// Function    : CR_ComputeDiffusivity
// Description : Compute the diffusivity of cosmic ray.
// Note        : This is a constant for now. The coefficient should depend on the cell condition.
// Parameter   :  
//-----------------------------------------------------------------------------------------
GPU_DEVICE
void CR_ComputeDiffusivity( /*dens*/ real &diff_cr_para, real &diff_cr_perp, const MicroPhy_t *Mic )
{
   // ne = dens/sim_mue/sim_mp;
   // diff_cr_para = minval( input_diff_para * (diff_cr_ne0 / ne),
   //                        input_diff_para                      );
   // diff_cr_perp = minval( input_diff_perp * (diff_cr_ne0 / ne),
   //                        input_diff_perp                      );
   
   diff_cr_para = Mic->CR_diff_coeff_para;
   diff_cr_perp = Mic->CR_diff_coeff_perp;
   // diff_cr_para = CR_DIFF_PARA;
   // diff_cr_perp = CR_DIFF_PERP;

} // FUNCTION : CR_ComputeDiffusivity

#endif // #if ( ( MODEL == HYDRO ) && defined CR_DIFFUION )

#endif // #ifndef __COMPUTE_COSMICRAYDIFFUSIVITY__
