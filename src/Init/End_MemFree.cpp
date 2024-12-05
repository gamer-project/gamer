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
         for (int i=0; i<PAR_NTYPE; i++) amr->patch[0][lv][PID]->NParType[i] = 0;
         free( amr->patch[0][lv][PID]->ParList );
         amr->patch[0][lv][PID]->ParList = NULL;
      }
#     endif

      delete amr;    amr = NULL;
   }


// 2. BaseP
   delete [] BaseP;  BaseP = NULL;


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

      delete Che_FieldData;   Che_FieldData = NULL;
   }
#  endif


// 4. dump table
   delete [] DumpTable;    DumpTable = NULL;


// 5. MPI buffers used by LOAD_BALANCE
#  ifdef LOAD_BALANCE
   LB_GetBufferData_MemFree();
#  endif


// 6. star formation random number generator
#  ifdef STAR_FORMATION
   SF_FreeRNG();
#  endif


// 7. user-defined table for grid refinement
   for (int lv=0; lv<NLEVEL-1; lv++)   free( FlagTable_User[lv] );


// 8. user-defined derived fields
   delete [] UserDerField_Label;    UserDerField_Label = NULL;
   delete [] UserDerField_Unit;     UserDerField_Unit  = NULL;


// 9. table of refinement region for OPT__UM_IC_NLEVEL>1
   delete [] UM_IC_RefineRegion;    UM_IC_RefineRegion = NULL;


// 10. global AMR structure
   delete GlobalTree;   GlobalTree = NULL;


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "done\n" );

} // FUNCTION : End_MemFree
