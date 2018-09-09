#include "GAMER.h"

#if ( MODEL == HYDRO || MODEL == SR_HYDRO  )


// 1: check negative values
// 2: check values close to the floor
// 3: check E is smaller than or equal to |M|
#define CHECK_MODE         1

// check if density/pressure is smaller than CLOSE_FACTOR*floor
#if   ( CHECK_MODE == 1 )
#  define CLOSE_FACTOR     NULL_REAL
#elif ( CHECK_MODE == 2 )
#  define CLOSE_FACTOR     2.0
#else
#  error : ERROR : only support CHECK_MODE == 1/2 !!
#endif



//-------------------------------------------------------------------------------------------------------
// Function    :  SRHydro_Aux_Check_Negative
// Description :  Check if there is any cell with negative density or pressure
//
// Parameter   :  lv       : Target refinement level
//                Mode     : 1 : Check negative density
//                           2 : Check negative pressure (and entropy when DUAL_ENERGY == DE_ENPY)
//                           3 : Check ratio of magnitute of momemtum to energy
//                           4 : Check above all criteria
//                comment  : You can put the location where this function is invoked in this string
//-------------------------------------------------------------------------------------------------------
void SRHydro_Aux_Check_Negative( const int lv, const int Mode, const char *comment )
{

// check
   if ( lv < 0  ||  lv >= NLEVEL )  Aux_Error( ERROR_INFO, "incorrect parameter %s = %d !!\n", "lv", lv );
   if ( Mode < 1  ||  Mode > 4 )    Aux_Error( ERROR_INFO, "incorrect parameter %s = %d !!\n", "Mode", Mode );

   const real Gamma_m1        = GAMMA - (real)1.0;
   const bool CheckMinPres_No = false;

   int  Pass = true;
   real Pres, Fluid[NCOMP_TOTAL];

// set the minimum thresholds for this check
   const real DensCheck  = ( CHECK_MODE == 1 ) ? 0.0 : CLOSE_FACTOR*MIN_DENS;
   const real PresCheck  = ( CHECK_MODE == 1 ) ? 0.0 : CLOSE_FACTOR*MIN_PRES;
   const real RatioCheck = 1.0;

   for (int TargetRank=0; TargetRank<MPI_NRank; TargetRank++)
   {
      if ( MPI_Rank == TargetRank )
      {
         for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
         for (int k=0; k<PATCH_SIZE; k++)
         for (int j=0; j<PATCH_SIZE; j++)
         for (int i=0; i<PATCH_SIZE; i++)
         {
            for (int v=0; v<NCOMP_TOTAL; v++)   Fluid[v] = amr->patch[ amr->FluSg[lv] ][lv][PID]->fluid[v][k][j][i];
/*
            Pres = CPU_GetPressure( Fluid[DENS], Fluid[MOMX], Fluid[MOMY], Fluid[MOMZ], Fluid[ENGY],
                                    Gamma_m1, CheckMinPres_No, NULL_REAL );
*/
            if ( Mode == 1  ||  Mode == 4 ) // check negative density
            {
               if ( Fluid[DENS] <= DensCheck )
               {
                  if ( Pass )
                  {
                     Aux_Message( stderr, "\"%s\" : <%s> FAILED at level %2d, Time = %13.7e, Step = %ld !!\n",
                                  comment, __FUNCTION__, lv, Time[lv], Step );
                     Aux_Message( stderr, "%4s  %7s  %19s  %10s  %21s  %21s  %21s  %21s  %21s\n",
                                  "Rank", "PID", "Patch Corner", "Grid ID", "Dens", "MomX", "MomY", "MomZ", "Engy" );
//                     Aux_Message( stderr, "  %21s\n", "Pres" );

                     Pass = false;
                  }

                  Aux_Message( stderr, "%4d  %7d  (%5d,%5d,%5d)  (%2d,%2d,%2d)  %21.14e  %21.14e  %21.14e  %21.14e  %21.14e",
                               MPI_Rank, PID, amr->patch[0][lv][PID]->corner[0],
                                              amr->patch[0][lv][PID]->corner[1],
                                              amr->patch[0][lv][PID]->corner[2], i, j, k,
                               Fluid[DENS], Fluid[MOMX], Fluid[MOMY], Fluid[MOMZ], Fluid[ENGY] );
//                  Aux_Message( stderr, "  %21.14e\n", Pres );
               }
            } // if ( Mode == 1  ||  Mode == 4 )


            if ( Mode == 2  ||  Mode == 4 ) // check negative pressure
            {
               if ( Pres <= PresCheck )
               {
                  if ( Pass )
                  {
                     Aux_Message( stderr, "\"%s\" : <%s> FAILED at level %2d, Time = %13.7e, Step = %ld !!\n",
                                  comment, __FUNCTION__, lv, Time[lv], Step );
                     Aux_Message( stderr, "%4s  %7s  %19s  %10s  %21s  %21s  %21s  %21s  %21s\n",
                                  "Rank", "PID", "Patch Corner", "Grid ID", "Dens", "MomX", "MomY", "MomZ", "Engy" );
//                     Aux_Message( stderr, "  %21s\n", "Pres" );

                     Pass = false;
                  }

                  Aux_Message( stderr, "%4d  %7d  (%5d,%5d,%5d)  (%2d,%2d,%2d)  %21.14e  %21.14e  %21.14e  %21.14e  %21.14e",
                               MPI_Rank, PID, amr->patch[0][lv][PID]->corner[0],
                                              amr->patch[0][lv][PID]->corner[1],
                                              amr->patch[0][lv][PID]->corner[2], i, j, k,
                               Fluid[DENS], Fluid[MOMX], Fluid[MOMY], Fluid[MOMZ], Fluid[ENGY] );
//                  Aux_Message( stderr, "  %21.14e\n", Pres );
               }
            } // if ( Mode == 2  ||  Mode == 4 )

            if ( Mode == 3  ||  Mode == 4 ) // check ratio of magnitute of momentum to energy
            {
               real M = SQRT(SQR(Fluid[MOMX]) + SQR(MOMY) + SQR(MOMZ));

               if ( M / Fluid[ENGY] >= RatioCheck )
               {
                  if ( Pass )
                  {
                     Aux_Message( stderr, "\"%s\" : <%s> FAILED at level %2d, Time = %13.7e, Step = %ld !!\n",
                                  comment, __FUNCTION__, lv, Time[lv], Step );
                     Aux_Message( stderr, "%4s  %7s  %19s  %10s  %21s  %21s  %21s  %21s  %21s\n",
                                  "Rank", "PID", "Patch Corner", "Grid ID", "Dens", "MomX", "MomY", "MomZ", "Engy" );
//                     Aux_Message( stderr, "  %21s\n", "Pres" );

                     Pass = false;
                  }

                  Aux_Message( stderr, "%4d  %7d  (%5d,%5d,%5d)  (%2d,%2d,%2d)  %21.14e  %21.14e  %21.14e  %21.14e  %21.14e",
                               MPI_Rank, PID, amr->patch[0][lv][PID]->corner[0],
                                              amr->patch[0][lv][PID]->corner[1],
                                              amr->patch[0][lv][PID]->corner[2], i, j, k,
                               Fluid[DENS], Fluid[MOMX], Fluid[MOMY], Fluid[MOMZ], Fluid[ENGY] );
//                  Aux_Message( stderr, "  %21.14e\n", Pres );
               }
            } // if ( Mode == 3  ||  Mode == 4 )

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

} // FUNCTION : SRHydro_Aux_Check_Negative



#endif // #if ( MODEL == SR_HYDRO )
