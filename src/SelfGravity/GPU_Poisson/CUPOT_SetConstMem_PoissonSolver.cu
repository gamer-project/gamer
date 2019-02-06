#include "Macro.h"
#include "CUPOT.h"

#if ( defined GRAVITY  &&  defined GPU )


__constant__ real c_Mp[3];
__constant__ real c_Mm[3];




//-------------------------------------------------------------------------------------------------------
// Function    :  CUPOT_SetConstMem_PoissonSolver
// Description :  Set the constant memory used by CUPOT_PoissonSolver_SOR_10to14cube()
//
// Note        :  1. Adopt the suggested approach for CUDA version >= 5.0
//                2. Invoked by CUAPI_Set_Default_GPU_Parameter()
//
// Parameter   :  None
//
// Return      :  0/-1 : successful/failed
//---------------------------------------------------------------------------------------------------
__host__
int CUPOT_SetConstMem_PoissonSolver()
{

   const real h_Mp[3] = { -3.0/32.0, +30.0/32.0, +5.0/32.0 };
   const real h_Mm[3] = { +5.0/32.0, +30.0/32.0, -3.0/32.0 };

   if (  cudaSuccess != cudaMemcpyToSymbol( c_Mp, h_Mp, 3*sizeof(real), 0, cudaMemcpyHostToDevice)  )
   return -1;

   if (  cudaSuccess != cudaMemcpyToSymbol( c_Mm, h_Mm, 3*sizeof(real), 0, cudaMemcpyHostToDevice)  )
   return -1;

   return 0;

} // FUNCTION : CUPOT_SetConstMem_PoissonSolver



#endif // #if ( defined GRAVITY  &&  defined GPU )
