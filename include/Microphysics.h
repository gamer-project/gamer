#ifndef __MICROPHYSICS__
#define __MICROPHYSICS__



#include "Macro.h"
#include "Typedef.h"



//-------------------------------------------------------------------------------------------------------
// Structure   :  MicroPhy_t
// Description :  Data structure storing the microphysics variables (e.g. diffusion, heat condution) to be passed to the CPU/GPU solvers
//
// Data Member :  safety          :
//                diff_coeff_para : Diffusion coefficient of cosmic ray along magnetic field.
//                diff_coeff_perp : Diffusion coefficient of cosmic ray perpendicular to magnetic field.
//
// Method      :  None --> It seems that CUDA does not support functions in a struct
//-------------------------------------------------------------------------------------------------------
struct MicroPhy_t
{

#  ifdef CR_DIFFUSION
   double CR_safety;
   double CR_diff_coeff_para;
   double CR_diff_coeff_perp;
#  endif

}; // struct MicroPhy_t

#endif // #ifndef __MICROPHYSICS__
