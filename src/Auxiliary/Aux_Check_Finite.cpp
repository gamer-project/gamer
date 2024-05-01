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

   int NVar = NCOMP_TOTAL;    // fluid
#  ifdef GRAVITY
   NVar ++;                   // potential
#  endif
#  if ( MODEL == HYDRO )
   NVar ++;                   // pressure
#  ifdef MHD
   NVar ++;                   // magnetic energy
#  endif
#  endif // HYDRO

   int Pass = true, NextIdx;
   real *Data = new real [NVar];


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
               NextIdx = 0;

               for (int v=0; v<NCOMP_TOTAL; v++) {
#                 if ( ELBDM_SCHEME == ELBDM_HYBRID )
//                ignore stub field on fluid levels in hybrid scheme
                  if ( !amr->use_wave_flag[lv]  &&  v == STUB )
                     Data[ NextIdx ++ ] = 0.0;
                  else
#                 endif
                     Data[ NextIdx ++ ] = amr->patch[ amr->FluSg[lv] ][lv][PID]->fluid[v][k][j][i];
               }

#              ifdef GRAVITY
               Data[ NextIdx ++ ] = amr->patch[ amr->PotSg[lv] ][lv][PID]->pot[k][j][i];
#              endif

#              if ( MODEL == HYDRO )
#              ifdef MHD
               const real Emag = MHD_GetCellCenteredBEnergyInPatch( lv, PID, i, j, k, amr->MagSg[lv] );
#              else
               const real Emag = NULL_REAL;
#              endif
               const real Pres = Hydro_Con2Pres( Data[DENS], Data[MOMX], Data[MOMY], Data[MOMZ], Data[ENGY], Data+NCOMP_FLUID,
                                                 false, NULL_REAL, Emag,
                                                 EoS_DensEint2Pres_CPUPtr, EoS_GuessHTilde_CPUPtr, EoS_HTilde2Temp_CPUPtr,
                                                 EoS_AuxArray_Flt,
                                                 EoS_AuxArray_Int, h_EoS_Table, NULL );
               Data[ NextIdx ++ ] = Pres;
#              ifdef MHD
               Data[ NextIdx ++ ] = Emag;
#              endif
#              endif // HYDRO

               for (int v=0; v<NVar; v++)
               {
                  if ( ! Aux_IsFinite(Data[v]) )
                  {
                     if ( Pass )
                     {
                        Aux_Message( stderr, "\"%s\" : <%s> FAILED at level %2d, Time = %13.7e, Step = %ld !!\n",
                                     comment, __FUNCTION__, lv, Time[lv], Step );
                        Aux_Message( stderr, "%4s  %7s  %45s  %8s  %14s\n",
                                     "Rank", "PatchID", "Coordinates", "Field", "Value" );

                        Pass = false;
                     }

                     Aux_Message( stderr, "%4d  %7d  (%13.7e, %13.7e, %13.7e)  %8d  %14.7e\n",
                                  MPI_Rank, PID,
                                  amr->patch[0][lv][PID]->EdgeL[0] + (i+0.5)*amr->dh[lv],
                                  amr->patch[0][lv][PID]->EdgeL[1] + (j+0.5)*amr->dh[lv],
                                  amr->patch[0][lv][PID]->EdgeL[2] + (k+0.5)*amr->dh[lv],
                                  v, Data[v] );
                  } // if ( ! Aux_IsFinite(Data[v]) )
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
         Aux_Message( stdout, "\"%s\" : <%s> PASSED at level %2d, Time = %13.7e, Step = %ld\n",
                      comment, __FUNCTION__, lv, Time[lv], Step );
   }

   else
      MPI_Exit();


   delete [] Data;

} // FUNCTION : Aux_Check_Finite
