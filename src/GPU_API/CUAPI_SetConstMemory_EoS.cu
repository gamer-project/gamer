#include "CUAPI.h"
#include "CUDA_ConstMemory.h"

#if ( defined GPU  &&  MODEL == HYDRO )

extern real *d_EoS_Table[EOS_NTABLE_MAX];




//-------------------------------------------------------------------------------------------------------
// Function    :  CUAPI_SetConstMemory_EoS
// Description :  Set the EoS constant memory variables on GPU
//
// Note        :  1. Adopt the suggested approach for CUDA version >= 5.0
//                2. Invoked by EoS_Init()
//
// Parameter   :  None
//
// Return      :  c_EoS_AuxArray_Flt[], c_EoS_AuxArray_Int[], c_EoS_Table[]
//                EoS.AuxArrayDevPtr_Flt, EoS.AuxArrayDevPtr_Int, EoS.Table
//---------------------------------------------------------------------------------------------------
void CUAPI_SetConstMemory_EoS()
{

// copy data to constant memory
   CUDA_CHECK_ERROR(  cudaMemcpyToSymbol( c_EoS_AuxArray_Flt, EoS_AuxArray_Flt, EOS_NAUX_MAX  *sizeof(double) )  );
   CUDA_CHECK_ERROR(  cudaMemcpyToSymbol( c_EoS_AuxArray_Int, EoS_AuxArray_Int, EOS_NAUX_MAX  *sizeof(int   ) )  );
   CUDA_CHECK_ERROR(  cudaMemcpyToSymbol( c_EoS_Table,        d_EoS_Table,      EOS_NTABLE_MAX*sizeof(real* ) )  );

// obtain the constant-memory pointers
   CUDA_CHECK_ERROR(  cudaGetSymbolAddress( (void **)&EoS.AuxArrayDevPtr_Flt, c_EoS_AuxArray_Flt )   );
   CUDA_CHECK_ERROR(  cudaGetSymbolAddress( (void **)&EoS.AuxArrayDevPtr_Int, c_EoS_AuxArray_Int )   );
   CUDA_CHECK_ERROR(  cudaGetSymbolAddress( (void **)&EoS.Table,              c_EoS_Table        )   );

} // FUNCTION : CUAPI_SetConstMemory_EoS



#endif // #if ( defined GPU  &&  MODEL == HYDRO )
