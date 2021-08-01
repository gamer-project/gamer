#include "CUAPI.h"
#include "CUDA_ConstMemory.h"

#if ( defined GPU  &&  defined GRAVITY )




//-------------------------------------------------------------------------------------------------------
// Function    :  CUAPI_SetConstMemory_ExtAccPot
// Description :  Set the constant memory variables on GPU used by the external acceleration and
//                potential routines
//
// Note        :  1. Adopt the suggested approach for CUDA version >= 5.0
//                2. Invoked by CUAPI_SetConstMemory()
//                3. EXT_ACC_NAUX_MAX and EXT_POT_NAUX_MAX are defined in Macro.h
//
// Parameter   :  None
//
// Return      :  c_ExtAcc_AuxArray[], c_ExtPot_AuxArray_Flt[], c_ExtPot_AuxArray_Int[]
//---------------------------------------------------------------------------------------------------
void CUAPI_SetConstMemory_ExtAccPot()
{

   if ( OPT__EXT_ACC )
      CUDA_CHECK_ERROR(  cudaMemcpyToSymbol( c_ExtAcc_AuxArray,     ExtAcc_AuxArray,     EXT_ACC_NAUX_MAX*sizeof(double) )  );

   if ( OPT__EXT_POT ) {
      CUDA_CHECK_ERROR(  cudaMemcpyToSymbol( c_ExtPot_AuxArray_Flt, ExtPot_AuxArray_Flt, EXT_POT_NAUX_MAX*sizeof(double) )  );
      CUDA_CHECK_ERROR(  cudaMemcpyToSymbol( c_ExtPot_AuxArray_Int, ExtPot_AuxArray_Int, EXT_POT_NAUX_MAX*sizeof(int)    )  );
   }

} // FUNCTION : CUAPI_SetConstMemory_ExtAccPot



#endif // #if ( defined GPU  &&  defined GRAVITY )
