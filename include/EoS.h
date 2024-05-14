#ifndef __EOS__
#define __EOS__



#include "Macro.h"
#include "Typedef.h"




//-------------------------------------------------------------------------------------------------------
// Structure   :  EoS_t
// Description :  Data structure storing the EoS variables to be passed to the CPU/GPU solvers
//
// Note        :  1. This object stores the GPU addresses and thus should never be used in CPU codes
//                   (except for codes shared by both CPU and GPU with a function parameter EoS_t EoS)
//                   --> For CPU-only codes, use the global variables EoS_AuxArray_Flt/Int, EoS_*_CPUPtr,
//                       and h_EoS_Table[] instead
//                2. When MODEL != HYDRO, we still define EoS_t for compilation but remove all data members
//
// Data Member :  *_AuxArrayDevPtr_* : Auxiliary array pointers
//                *_FuncPtr          : Function pointers to the major EoS functions
//                Table              : Pointer arrays to the EoS tables
//
// Method      :  None --> It seems that CUDA does not support functions in a struct
//-------------------------------------------------------------------------------------------------------
struct EoS_t
{

#  if ( MODEL == HYDRO )
// auxiliary array pointers
   double *AuxArrayDevPtr_Flt;
   int    *AuxArrayDevPtr_Int;

// function pointers
   EoS_GUESS_t   GuessHTilde_FuncPtr;
   EoS_H2TEM_t   HTilde2Temp_FuncPtr;
   EoS_TEM2H_t   Temp2HTilde_FuncPtr;
   EoS_DE2P_t    DensEint2Pres_FuncPtr;
   EoS_DP2E_t    DensPres2Eint_FuncPtr;
   EoS_DP2C_t    DensPres2CSqr_FuncPtr;
   EoS_DE2T_t    DensEint2Temp_FuncPtr;
   EoS_DT2P_t    DensTemp2Pres_FuncPtr;
   EoS_DE2S_t    DensEint2Entr_FuncPtr;
   EoS_GENE_t    General_FuncPtr;
#  ifdef COSMIC_RAY
   EoS_CRE2CRP_t CREint2CRPres_FuncPtr;
#  endif

// table pointers
   real **Table;
#  endif

}; // struct EoS_t



#endif // #ifndef __EOS__
