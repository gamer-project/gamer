#include "GAMER.h"





//-------------------------------------------------------------------------------------------------------
// Function    :  Aux_Check_Refinement
// Description :  Verify the refinement result at the input level
//                --> Output grids at level lv that should be refined into level lv+1 but were not refined
//
// Note        :  1. This test will not pass if some patches are not refined due to the proper-nestiing constraint
//                2. This test should always pass for the initial condition constructed by the UM_Init function
//                3. This test only checks density and the option "OPT__FLAG_RHO" should be turned on
//
// Parameter   :  lv       : Target refinement level
//                comment  : You can put the location where this function is invoked in this string
//-------------------------------------------------------------------------------------------------------
void Aux_Check_Refinement( const int lv, const char *comment )
{

// check
#  ifndef DENS
   Aux_Message( stderr, "WARNING : function \"%s\" is supported only if the variable \"DENS\" is defined !!\n",
                __FUNCTION__ );
   OPT__CK_REFINE = false;
   return;
#  else

// nothing to check at the maximum level
   if ( lv == NLEVEL-1 )
   {
      if ( MPI_Rank == 0 )
         Aux_Message( stderr, "WARNING : function \"%s\" should NOT be applied to the finest level !!\n",
                      __FUNCTION__ );
      return;
   }

// this check function must work with the table "FlagTable_Rho"
   if ( !OPT__FLAG_RHO )
   {
      if ( MPI_Rank == 0 )
         Aux_Message( stderr, "WARNING : function \"%s\" must work with the option \"%s\"  !!\n",
                      __FUNCTION__, "OPT__FLAG_RHO == 1" );
      return;
   }


   int   Pass        = true;
   bool  PatchTest   = true;
   int   CountPatch  = 0;
   int   CountGrid   = 0;
   real  Rho;

   for (int TargetRank=0; TargetRank<MPI_NRank; TargetRank++)
   {
      if ( MPI_Rank == TargetRank )
      {
         for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
         {
            if ( amr->patch[0][lv][PID]->son == -1 )
            {
               PatchTest = true;

               for (int k=0; k<PATCH_SIZE; k++)
               for (int j=0; j<PATCH_SIZE; j++)
               for (int i=0; i<PATCH_SIZE; i++)
               {
                  Rho = amr->patch[ amr->FluSg[lv] ][lv][PID]->fluid[DENS][k][j][i];

                  if ( Rho > FlagTable_Rho[lv] )
                  {
                     if ( Pass )
                     {
                        Aux_Message( stderr, "\"%s\" : <%s> FAILED at level %2d, Time = %13.7e, Step = %ld !!\n",
                                     comment, __FUNCTION__, lv, Time[lv], Step );
                        Aux_Message( stderr, "%4s\t%7s\t\t%19s\t\t%10s\t%14s\n",
                                     "Rank", "PatchID", "Patch Corner", "Grid ID", "Density" );

                        Pass = false;
                     }

                     Aux_Message( stderr, "%4d\t%7d\t\t(%10d,%10d,%10d)\t\t(%2d,%2d,%2d)\t%14.7e\n",
                                  MPI_Rank, PID, amr->patch[0][lv][PID]->corner[0],
                                                 amr->patch[0][lv][PID]->corner[1],
                                                 amr->patch[0][lv][PID]->corner[2], i, j, k, Rho );

                     PatchTest = false;
                     CountGrid ++;

                  } // if ( Rho > FlagTable_Rho[lv] )
               } // i, j, k

               if ( !PatchTest )    CountPatch ++;

            } // if ( amr->patch[0][lv][PID]->son == -1 )
         } // for (int PID=0; PID<NPachComma[lv][1]; PID++)
      } // if ( MPI_Rank == TargetRank )

      MPI_Bcast( &Pass, 1, MPI_INT, TargetRank, MPI_COMM_WORLD );

      MPI_Barrier( MPI_COMM_WORLD );

   } // for (int TargetRank=0; TargetRank<MPI_NRank; TargetRank++)

   if ( Pass )
   {
      if ( MPI_Rank == 0 )
         Aux_Message( stdout, "\"%s\" : <%s> PASSED at level %2d, Time = %13.7e, Step = %ld\n",
                      comment, __FUNCTION__, lv, Time[lv], Step );
   }

   else
   {
      int CountPatch_Sum, CountGrid_Sum, NPatch_Sum;

      MPI_Reduce( &CountPatch, &CountPatch_Sum, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD );
      MPI_Reduce( &CountGrid , &CountGrid_Sum , 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD );
      MPI_Reduce( &amr->NPatchComma[lv][1], &NPatch_Sum , 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD );

      if ( MPI_Rank == 0 )
         Aux_Message( stderr, "# of patches = %d, # of failed patches = %d, # of failed grids = %d\n",
                      NPatch_Sum, CountPatch_Sum, CountGrid_Sum );
   }
#  endif // #ifndef DENS ... else ...

} // FUNCTION : Aux_Check_Refinement
