#include "Macro.h"
#include "CUFLU.h"

#if ( MODEL == HYDRO  &&  defined GPU )



#ifdef UNSPLIT_GRAVITY
#include "CUPOT.h"
__constant__ double ExtAcc_AuxArray_d_Flu[EXT_ACC_NAUX_MAX];

//-------------------------------------------------------------------------------------------------------
// Function    :  CUFLU_SetConstMem_FluidSolver_ExtAcc
// Description :  Set the constant memory of ExtAcc_AuxArray_d_Flu used by CUFLU_FluidSolver_CTU/MHM
//
// Note        :  1. Adopt the suggested approach for CUDA version >= 5.0
//                2. Invoked by CUAPI_Init_ExternalAccPot()
//
// Parameter   :  None
//
// Return      :  0/-1 : successful/failed
//---------------------------------------------------------------------------------------------------
__host__
int CUFLU_SetConstMem_FluidSolver_ExtAcc( double ExtAcc_AuxArray_h[] )
{

   if (  cudaSuccess != cudaMemcpyToSymbol( ExtAcc_AuxArray_d_Flu, ExtAcc_AuxArray_h, EXT_ACC_NAUX_MAX*sizeof(double),
                                            0, cudaMemcpyHostToDevice)  )
      return -1;

   else
      return 0;

} // FUNCTION : CUFLU_SetConstMem_FluidSolver_ExtAcc
#endif // #ifdef UNSPLIT_GRAVITY


#if ( NCOMP_PASSIVE > 0 )
__constant__ int NormIdx[NCOMP_PASSIVE];

//-------------------------------------------------------------------------------------------------------
// Function    :  CUFLU_SetConstMem_FluidSolver_NormIdx
// Description :  Set the constant memory of NormIdx[] used by CUFLU_FluidSolver_CTU/MHM
//
// Note        :  1. Adopt the suggested approach for CUDA version >= 5.0
//                2. Invoked by CUAPI_Set_Default_GPU_Parameter()
//
// Parameter   :  None
//
// Return      :  0/-1 : successful/failed
//---------------------------------------------------------------------------------------------------
__host__
int CUFLU_SetConstMem_FluidSolver_NormIdx( int NormIdx_h[] )
{

   if (  cudaSuccess != cudaMemcpyToSymbol( NormIdx, NormIdx_h, NCOMP_PASSIVE*sizeof(int),
                                            0, cudaMemcpyHostToDevice)  )
      return -1;

   else
      return 0;

} // FUNCTION : CUFLU_SetConstMem_FluidSolver_NormIdx

#else // #if ( NCOMP_PASSIVE > 0 )
__constant__ int *NormIdx = NULL;

#endif // #if ( NCOMP_PASSIVE > 0 ) ... else ...



#endif // #if ( MODEL == HYDRO  &&  defined GPU )
