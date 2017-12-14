#include "GAMER.h"




//-------------------------------------------------------------------------------------------------------
// Function    :  Aux_Check_Finite
// Description :  Verify that all variables at level "lv" are neither "NaN" nor "Infinite"
//
// Note        :  It only checks data stored in the current Sg
//
// Parameter   :  lv       : Target refinement level
//                comment  : You can put the location where this function is invoked in this string
//-------------------------------------------------------------------------------------------------------
void Aux_Check_Finite( const int lv, const char *comment )
{

#  ifdef GRAVITY
   const int NVar = NCOMP_TOTAL+1;
#  else
   const int NVar = NCOMP_TOTAL;
#  endif

   int Pass = true;
   real Data[NVar];


   for (int TargetRank=0; TargetRank<MPI_NRank; TargetRank++)
   {
      if ( MPI_Rank == TargetRank )
      {
         for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
         {
            for (int k=0; k<PATCH_SIZE; k++)
            for (int j=0; j<PATCH_SIZE; j++)
            for (int i=0; i<PATCH_SIZE; i++)
            {
               for (int v=0; v<NCOMP_TOTAL; v++)
               Data[          v] = amr->patch[ amr->FluSg[lv] ][lv][PID]->fluid[v][k][j][i];

#              ifdef GRAVITY
               Data[NCOMP_TOTAL] = amr->patch[ amr->PotSg[lv] ][lv][PID]->pot[k][j][i];
#              endif

               for (int v=0; v<NVar; v++)
               {
                  if ( ! isfinite(Data[v]) )
                  {
                     if ( Pass )
                     {
                        Aux_Message( stderr, "\"%s\" : <%s> FAILED at level %2d, Time = %13.7e, Step = %ld !!\n",
                                     comment, __FUNCTION__, lv, Time[lv], Step );
                        Aux_Message( stderr, "%4s\t%7s\t\t%19s\t%8s\n", "Rank", "PatchID", "Coordinate",
                                                                        "Variable" );

                        Pass = false;
                     }

                     Aux_Message( stderr, "%4d\t%7d\t\t(%5d,%5d,%5d)\t%8d\n",
                                  MPI_Rank, PID,
                                  i*amr->scale[lv] + amr->patch[0][lv][PID]->corner[0],
                                  j*amr->scale[lv] + amr->patch[0][lv][PID]->corner[1],
                                  k*amr->scale[lv] + amr->patch[0][lv][PID]->corner[2],
                                  v );
                  } // if ( ! isfinite(Data[v]) )
               } // for (int v=0; v<NVar; v++)
            } // i,j,k
         } // for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
      } // if ( MPI_Rank == TargetRank )

      MPI_Bcast( &Pass, 1, MPI_INT, TargetRank, MPI_COMM_WORLD );

      MPI_Barrier( MPI_COMM_WORLD );

   } // for (int TargetRank=0; TargetRank<MPI_NRank; TargetRank++)


   if ( Pass )
   {
      if ( MPI_Rank == 0 )
         Aux_Message( stdout, "\"%s\" : <%s> PASSED at level %2d, Time = %13.7e, Step = %ld \n",
                      comment, __FUNCTION__, lv, Time[lv], Step );
   }

   else
      MPI_Exit();

} // FUNCTION : Aux_Check_Finite
