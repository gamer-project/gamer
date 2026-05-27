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
// Method      :  None --> It seems that CUDA does not support functions in a struct
//-------------------------------------------------------------------------------------------------------
struct MicroPhy_t
{

// somehow the structure itself cannot be empty, so we declare a useless bool variable to avoid the issue
// --> error message from valgrind: Address 0x176aa160 is 0 bytes after a block of size 0 alloc'd
   bool useless;

#  ifdef CR_DIFFUSION
   real CR_safety;
   real CR_diff_coeff_para;
   real CR_diff_coeff_perp;
   real CR_diff_min_b;
#  endif

#  ifdef CONDUCTION
   bool CondSaturation;
   bool CondSatWhistler;
   int  CondType;
   int  CondFluxType;
   real Cond_safety;
   real CondConstCoeff;
   real CondMaxDiffusivity;
   real CondSpitzerFraction;
   real CondPrefactor;
   real CondCoulombLog;
   real CondMFPConst;
   real CondMue;
   real CondSpecificHeat;
   real CondPresConv;
#  endif // #ifdef CONDUCTION

#  ifdef VISCOSITY
   bool ViscSaturation;
   bool ViscBounds;
   int  ViscType;
   int  ViscFluxType;
   int  ViscCoeffType;
   real Visc_safety;
   real ViscConstCoeff;
   real ViscMaxDiffusivity;
   real ViscSpitzerFraction;
   real ViscPrefactor;
   real ViscCoulombLog;
   real ViscMFPConst;
   real ViscMui;
   real ViscThermalSpeedConv;
#  endif // #ifdef VISCOSITY

}; // struct MicroPhy_t



#endif // #ifndef __MICROPHYSICS__
