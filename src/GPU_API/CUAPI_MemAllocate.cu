#include "CUAPI.h"

#ifdef GPU



void CUAPI_SetMemSize( int &GPU_NStream, int &Flu_GPU_NPGroup, int &Pot_GPU_NPGroup, int &Che_GPU_NPGroup,
                       int &Src_GPU_NPGroup );
int CUAPI_MemAllocate_Fluid( const int Flu_NPG, const int Pot_NPG, const int Src_NPG, const int GPU_NStream );
#ifdef GRAVITY
int CUAPI_MemAllocate_PoissonGravity( const int Pot_NPatchGroup );
#endif




//-------------------------------------------------------------------------------------------------------
// Function    :  CUAPI_MemAllocate
// Description :  Allocate GPU memory for several global arrays
//-------------------------------------------------------------------------------------------------------
void CUAPI_MemAllocate()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );


#  ifndef GRAVITY
   int POT_GPU_NPGROUP = NULL_INT;
#  endif
#  ifndef SUPPORT_GRACKLE
   int CHE_GPU_NPGROUP = NULL_INT;
#  endif

// set parameters related to GPU memory allocation
   CUAPI_SetMemSize( GPU_NSTREAM, FLU_GPU_NPGROUP, POT_GPU_NPGROUP, CHE_GPU_NPGROUP, SRC_GPU_NPGROUP );


// allocate memory for all GPU solvers
   CUAPI_MemAllocate_Fluid( FLU_GPU_NPGROUP, POT_GPU_NPGROUP, SRC_GPU_NPGROUP, GPU_NSTREAM );

#  ifdef GRAVITY
   CUAPI_MemAllocate_PoissonGravity( POT_GPU_NPGROUP );
#  endif


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : CUAPI_MemAllocate



#endif // #ifdef GPU
