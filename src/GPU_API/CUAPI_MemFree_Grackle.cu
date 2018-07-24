#include "CUAPI.h"

#if ( defined GPU  &&  defined SUPPORT_GRACKLE )




//-------------------------------------------------------------------------------------------------------
// Function    :  CUAPI_MemFree_Grackle
// Description :  Free the device and host memory previously allocated by CUAPI_MemAllocate_Grackle()
//
// Parameter   :  None
//-------------------------------------------------------------------------------------------------------
void CUAPI_MemFree_Grackle()
{

// free the device memory
// --> not necessary for now


// free the host memory allocated by CUDA
   for (int t=0; t<2; t++)
   {
      if ( h_Che_Array[t] != NULL )    CUDA_CHECK_ERROR(  cudaFreeHost( h_Che_Array[t] )  );

      h_Che_Array[t] = NULL;
   }

} // FUNCTION : CUAPI_MemFree_Grackle



#endif // #if ( defined GPU  &&  defined SUPPORT_GRACKLE )
