#ifndef __MICROPHYSICS__
#define __MICROPHYSICS__



#include "Macro.h"
#include "Typedef.h"




//-------------------------------------------------------------------------------------------------------
// Structure   :  MicroPhy_t
// Description :  Data structure storing the microphysics variables (e.g. cosmic-ray diffusion coefficients)
//                to be passed to the CPU/GPU solvers
//
// Data Member :  CR_safety          : CFL safety factor of cosmic-ray diffusion (runtime parameter: DT__CR_DIFFUSION)
//                CR_diff_coeff_para : Cosmic-ray diffusion coefficients parallel/perpendicular to the
//                CR_diff_coeff_perp   magnetic field (runtime parameters: CR_DIFF_PARA, CR_DIFF_PERP)
//                CR_diff_min_b      : minimum magnetic field for enabling diffusion (runtime parameter: CR_DIFF_MIN_B)
//
//                CR_source          : Flag for enabling cosmic-ray source terms (runtime parameter: CR_SOURCE)
//                CR_stream          : Flag for enabling cosmic-ray streaming (runtime parameter: CR_STREAM)
//                CR_Ec_source       : Flag for including the CR energy source term (runtime parameter: CR_EC_SOURCE)
//                CR_vmax            : Maximum velocity (Vm, effective speed of light) (runtime parameter: CR_VMAX)
//                CR_sigma           : Parallel opacity in code-unit = Vm*CR_SIGMA (runtime parameter: CR_SIGMA)
//                CR_sigma_perp      : Perpendicular opacity in code-unit = Vm*CR_SIGMA_PERP (runtime parameter: CR_SIGMA_PERP)
//                CR_max_opacity     : Maximum opacity (runtime parameter: CR_MAX_OPACITY)
//                CR_tau_asym_lim    : Optical depth limit for asymptotic expansion (runtime parameter: CR_TAU_ASYM_LIM)
//                CR_taufact         : Tau factor for optical depth calculation (runtime parameter: CR_TAUFACT)
//                CR_vel_flx_flag    : Flag for adding CR sound speed to v_diff (runtime parameter: CR_VEL_FLX_FLAG)
//                CR_cfl             : CFL safety factor for two-moment cosmic-ray (runtime parameter: CR_CFL)
//
// Method      :  None --> It seems that CUDA does not support functions in a struct
//-------------------------------------------------------------------------------------------------------
struct MicroPhy_t
{

#  ifdef CR_DIFFUSION
   real CR_safety;
   real CR_diff_coeff_para;
   real CR_diff_coeff_perp;
   real CR_diff_min_b;
#  endif

#  ifdef CR_TWOMOMENT
   bool CR_source;
   bool CR_stream;
   bool CR_Ec_source;
   real CR_vmax;
   real CR_sigma;
   real CR_sigma_perp;
   real CR_max_opacity;
   real CR_tau_asym_lim;
   real CR_taufact;
   bool CR_vel_flx_flag;
   real CR_cfl;
#  endif

// somehow the structure itself cannot be empty, so we declare a useless bool variable to avoid the issue
// --> error message from valgrind: Address 0x176aa160 is 0 bytes after a block of size 0 alloc'd
   bool useless;

}; // struct MicroPhy_t



#endif // #ifndef __MICROPHYSICS__
