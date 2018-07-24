#include "CUAPI.h"

#if ( defined GPU  &&  defined SUPPORT_GRACKLE )




//-------------------------------------------------------------------------------------------------------
// Function    :  CUAPI_MemAllocate_Grackle
// Description :  Allocate GPU and CPU memory for the Grackle solver
//
// Parameter   :  Che_NPG : Number of patch groups evaluated simultaneously by GPU
//-------------------------------------------------------------------------------------------------------
void CUAPI_MemAllocate_Grackle( const int Che_NPG )
{

// nothing to do if Grackle is disabled
   if ( !GRACKLE_ACTIVATE )   return;


// size of the global memory array(s)
   const long Che_MemSize_In   = sizeof(real)*Che_NPG*CHE_NIN  *CUBE(PS2);
   const long Che_MemSize_Prep = sizeof(real)*Che_NPG*CHE_NPREP*CUBE(PS2);

// output the total memory requirement
   long TotalSize = Che_MemSize_In;

   if ( MPI_Rank == 0 )
      Aux_Message( stdout, "NOTE : total memory requirement in GPU Grackle solver = %ld MB\n", TotalSize/(1<<20) );


// allocate the device memory
// --> not necessary for now


// allocate the host memory by CUDA
   for (int t=0; t<2; t++)
   {
      CUDA_CHECK_ERROR(  cudaMallocHost( (void**) &h_Che_Array[t], Che_MemSize_Prep )  );
   }

} // FUNCTION : CUAPI_MemAllocate_Grackle



#endif // #if ( defined GPU  &&  defined SUPPORT_GRACKLE )
