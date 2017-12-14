#include "Macro.h"
#include "CUPOT.h"

#if ( defined GRAVITY  &&  defined GPU )


__constant__ real Mp[3];
__constant__ real Mm[3];




//-------------------------------------------------------------------------------------------------------
// Function    :  CUPOT_PoissonSolver_SetConstMem
// Description :  Set the constant memory used by CUPOT_PoissonSolver_SOR_10to14cube
//
// Note        :  Adopt the suggested approach for CUDA version >= 5.0
//
// Parameter   :  None
//
// Return      :  0/-1 : successful/failed
//---------------------------------------------------------------------------------------------------
int CUPOT_PoissonSolver_SetConstMem()
{

   const real Mp_h[3] = { -3.0/32.0, +30.0/32.0, +5.0/32.0 };
   const real Mm_h[3] = { +5.0/32.0, +30.0/32.0, -3.0/32.0 };

   if (  cudaSuccess != cudaMemcpyToSymbol( Mp, Mp_h, 3*sizeof(real), 0, cudaMemcpyHostToDevice)  )
   return -1;

   if (  cudaSuccess != cudaMemcpyToSymbol( Mm, Mm_h, 3*sizeof(real), 0, cudaMemcpyHostToDevice)  )
   return -1;

   return 0;

} // FUNCTION : CUPOT_PoissonSolver_SetConstMem



#endif // #if ( defined GRAVITY  &&  defined GPU )
