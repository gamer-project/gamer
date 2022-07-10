#include "CUAPI.h"
#ifdef GPU

//-------------------------------------------------------------------------------------------------------
// Function    :  CUAPI_Synchronize
// Description :  Block until the device has completed all preceding requested tasks
//
// Note        :  1. Replace the deprecated cudaThreadSynchronize() with cudaDeviceSynchronize()
//-------------------------------------------------------------------------------------------------------
void CUAPI_Synchronize()
{
   CUDA_CHECK_ERROR(  cudaDeviceSynchronize()  );
}

#endif // #ifdef GPU
