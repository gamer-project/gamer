#ifndef __EOS__
#define __EOS__



#include "Macro.h"
#include "Typedef.h"



//-------------------------------------------------------------------------------------------------------
// Structure   :  EoS_t
// Description :  Data structure storing the EoS variables to be passed to the CPU/GPU solvers
//
// Data Member :  *_AuxArrayDevPtr_* : Auxiliary array pointers
//                                     --> For GPU, we store the addresses of constant memory arrays, which
//                                         should NOT be used by host
//                *_FuncPtr          : Function pointers to the major EoS functions
//                Table              : Pointer arrays to the EoS tables
//                                     --> For GPU, we store the addresses of global memory arrays, which
//                                         should NOT be used by host
//
// Method      :  None --> It seems that CUDA does not support functions in a struct
//-------------------------------------------------------------------------------------------------------
struct EoS_t
{

// auxiliary array pointers
   double *AuxArrayDevPtr_Flt;
   int    *AuxArrayDevPtr_Int;

// function pointers
   EoS_DE2P_t DensEint2Pres_FuncPtr;
   EoS_DP2E_t DensPres2Eint_FuncPtr;
   EoS_DP2C_t DensPres2CSqr_FuncPtr;
   EoS_DE2T_t DensEint2Temp_FuncPtr;
   EoS_DT2P_t DensTemp2Pres_FuncPtr;
   EoS_DE2S_t DensEint2Entr_FuncPtr;
   EoS_GENE_t General_FuncPtr;

// table pointers
   real **Table;

}; // struct EoS_t



#endif // #ifndef __EOS__
