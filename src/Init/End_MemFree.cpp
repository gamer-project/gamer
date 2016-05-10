#include "Copyright.h"
#include "GAMER.h"

#ifdef LOAD_BALANCE
void LB_GetBufferData_MemFree();
#endif




//-------------------------------------------------------------------------------------------------------
// Function    :  End_MemFree
// Description :  Free memory allocated by the function "Init_MemoryAllocate"
//-------------------------------------------------------------------------------------------------------
void End_MemFree()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... ", __FUNCTION__ );


// a. deallocate the AMR structure
   if ( amr != NULL )   
   {
      delete amr;
      amr = NULL;
   }


// b. deallocate the BaseP
   if ( BaseP != NULL )   
   {
      delete [] BaseP;
      BaseP = NULL;
   }


// c. deallocate arrays for GPU (or CPU) solvers
#  ifdef GPU
      CUAPI_MemFree_Fluid( GPU_NSTREAM );
#  else
      End_MemFree_Fluid();
#  endif

#  ifdef GRAVITY
#     ifdef GPU
         CUAPI_MemFree_PoissonGravity();
#     else
         End_MemFree_PoissonGravity();
#     endif
#  endif


// d. deallocate the dump table
   if ( DumpTable != NULL )   
   {
      delete [] DumpTable;
      DumpTable = NULL;
   }


// e. free the MPI buffers used by LOAD_BALANCE
   #ifdef LOAD_BALANCE
   LB_GetBufferData_MemFree();
   #endif


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "done\n" );
   
} // FUNCTION : End_MemFree
