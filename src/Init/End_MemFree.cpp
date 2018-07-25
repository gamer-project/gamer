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


// 1. AMR structure
   if ( amr != NULL )
   {
//    free particle variables first to avoid warning messages when deleting patches with particles
#     ifdef PARTICLE
      for (int lv=0; lv<NLEVEL; lv++)
      for (int PID=0; PID<amr->num[lv]; PID++)
      {
         amr->patch[0][lv][PID]->NPar = 0;
         free( amr->patch[0][lv][PID]->ParList );
         amr->patch[0][lv][PID]->ParList = NULL;
      }
#     endif

      delete amr;
      amr = NULL;
   }


// 2. BaseP
   if ( BaseP != NULL )
   {
      delete [] BaseP;
      BaseP = NULL;
   }


// 3. arrays for GPU (or CPU) solvers
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

#  ifdef SUPPORT_GRACKLE
   if ( GRACKLE_ACTIVATE )
   {
      End_MemFree_Grackle();

      if ( Che_FieldData != NULL )
      {
         delete Che_FieldData;
         Che_FieldData = NULL;
      }
   }
#  endif


// 4. dump table
   if ( DumpTable != NULL )
   {
      delete [] DumpTable;
      DumpTable = NULL;
   }


// 5. MPI buffers used by LOAD_BALANCE
#  ifdef LOAD_BALANCE
   LB_GetBufferData_MemFree();
#  endif


// 6. star formation random number generator
#  ifdef STAR_FORMATION
   SF_FreeRNG();
#  endif


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "done\n" );

} // FUNCTION : End_MemFree
