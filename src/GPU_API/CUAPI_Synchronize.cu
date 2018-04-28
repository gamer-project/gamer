#include "CUAPI.h"
#ifdef GPU

//-------------------------------------------------------------------------------------------------------
// Function    :  CUAPI_Synchronize
// Description :  Block until the device has completed all preceding requested tasks
//-------------------------------------------------------------------------------------------------------
void CUAPI_Synchronize()
{
   CUDA_CHECK_ERROR(  cudaThreadSynchronize()  );
}

#endif // #ifdef GPU
