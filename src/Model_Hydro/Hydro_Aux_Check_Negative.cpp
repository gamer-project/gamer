#include "Copyright.h"
#include "GAMER.h"

#if ( MODEL == HYDRO )




//-------------------------------------------------------------------------------------------------------
// Function    :  Hydro_Aux_Check_Negative
// Description :  Check if there is any cell with negative density or pressure
//
// Parameter   :  lv       : Targeted refinement level
//                Mode     : 1 : Check negative density
//                           2 : Check negative pressure
//                           3 : Both
//                comment  : You can put the location where this function is invoked in this string
//-------------------------------------------------------------------------------------------------------
void Hydro_Aux_Check_Negative( const int lv, const int Mode, const char *comment )
{

// check
   if ( lv < 0  ||  lv >= NLEVEL )  Aux_Error( ERROR_INFO, "incorrect parameter %s = %d !!\n", "lv", lv );
   if ( Mode < 1  ||  Mode > 3 )    Aux_Error( ERROR_INFO, "incorrect parameter %s = %d !!\n", "Mode", Mode );


   const real Gamma_m1        = GAMMA - (real)1.0;
   const bool CheckMinPres_No = false;

   int  Pass = true;
   real Rho, Pres, Fluid[NCOMP];

   for (int TargetRank=0; TargetRank<MPI_NRank; TargetRank++)
   {
      if ( MPI_Rank == TargetRank )
      {
         for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
         for (int k=0; k<PATCH_SIZE; k++)
         for (int j=0; j<PATCH_SIZE; j++)
         for (int i=0; i<PATCH_SIZE; i++)
         {
            for (int v=0; v<NCOMP; v++)   Fluid[v] = amr->patch[ amr->FluSg[lv] ][lv][PID]->fluid[v][k][j][i];

            Rho  = Fluid[DENS];
            Pres = CPU_GetPressure( Fluid[DENS], Fluid[MOMX], Fluid[MOMY], Fluid[MOMZ], Fluid[ENGY],
                                    Gamma_m1, CheckMinPres_No, NULL_REAL );

            if ( Mode == 1  ||  Mode == 3 )
            {
               if ( Rho <= 0.0 )
               {
                  if ( Pass )
                  {
                     Aux_Message( stderr, "\"%s\" : <%s> FAILED at level %2d, Time = %13.7e, Step = %ld !!\n",
                                  comment, __FUNCTION__, lv, Time[lv], Step );
                     Aux_Message( stderr, "%4s\t%7s\t\t%19s\t%10s\t%14s\t%14s\n",
                                  "Rank", "PID", "Patch Corner", "Grid ID", "Density", "Pressure" );

                     Pass = false;
                  }

                  Aux_Message( stderr, "%4d\t%7d\t\t(%10d,%10d,%10d)\t(%2d,%2d,%2d)\t%14.7e\t%14.7e\n",
                               MPI_Rank, PID, amr->patch[0][lv][PID]->corner[0],
                                              amr->patch[0][lv][PID]->corner[1],
                                              amr->patch[0][lv][PID]->corner[2], i, j, k, Rho, Pres );
               }
            } // if ( Mode == 1  ||  Mode == 3 )

            if ( Mode == 2  ||  Mode == 3 )
            {
               if ( Pres <= 0.0 )
               {
                  if ( Pass )
                  {
                     Aux_Message( stderr, "\"%s\" : <%s> FAILED at level %2d, Time = %13.7e, Step = %ld !!\n",
                                  comment, __FUNCTION__, lv, Time[lv], Step );
                     Aux_Message( stderr, "%4s\t%7s\t\t%19s\t%10s\t%14s\t%14s\n",
                                  "Rank", "PID", "Patch Corner", "Grid ID", "Density", "Pressure" );

                     Pass = false;
                  }

                  Aux_Message( stderr, "%4d\t%7d\t\t(%10d,%10d,%10d)\t(%2d,%2d,%2d)\t%14.7e\t%14.7e\n",
                               MPI_Rank, PID, amr->patch[0][lv][PID]->corner[0],
                                              amr->patch[0][lv][PID]->corner[1],
                                              amr->patch[0][lv][PID]->corner[2], i, j, k, Rho, Pres );
               }
            } // if ( Mode == 2  ||  Mode == 3 )

         } // i,j,k,PID
      } // if ( MPI_Rank == TargetRank )

      MPI_Bcast( &Pass, 1, MPI_INT, TargetRank, MPI_COMM_WORLD );

      MPI_Barrier( MPI_COMM_WORLD );

   } // for (int TargetRank = 0; TargetRank<MPI_NRank; TargetRank++)


   if ( Pass )
   {
      if ( MPI_Rank == 0 )
         Aux_Message( stdout, "\"%s\" : <%s> PASSED at level %2d, Time = %13.7e, Step = %ld\n",
                      comment, __FUNCTION__, lv, Time[lv], Step );
   }

} // FUNCTION : Hydro_Aux_Check_Negative



#endif // #if ( MODEL == HYDRO )
