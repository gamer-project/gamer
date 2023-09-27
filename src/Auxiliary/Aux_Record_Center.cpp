#include "GAMER.h"
#include "TestProb.h"




//-------------------------------------------------------------------------------------------------------
// Function    :  Aux_Record_Center
// Description :  Record various center coordinates
//
// Note        :
//                2. For the center coordinates, it will record the position of maximum density, minimum potential,
//                   and center-of-mass
//                3. Output filename is fixed to "Record__Center"
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Aux_Record_Center()
{

   const char filename_center  [] = "Record__Center";
   const int  CountMPI            = 8;
   const int  DensRecMode         = 0;                               // different modes for recording the peak density (1/2/3: ELBDM/particle/total density)

   double dens, max_dens_loc=-__DBL_MAX__, max_dens_pos_loc[3];
   double pote, min_pote_loc=+__DBL_MAX__, min_pote_pos_loc[3];
   double send[CountMPI], (*recv)[CountMPI]=new double [MPI_NRank][CountMPI];

   const long   DensMode          = ( DensRecMode == 1 ) ? _DENS :
#                                   ifdef PARTICLE
                                    ( DensRecMode == 2 ) ? _PAR_DENS :
#                                   endif
                                                           _TOTAL_DENS;
   const bool   IntPhase_No       = false;
   const real   MinDens_No        = -1.0;
   const real   MinPres_No        = -1.0;
   const real   MinTemp_No        = -1.0;
   const real   MinEntr_No        = -1.0;
   const bool   DE_Consistency_No = false;
#  ifdef PARTICLE
   const bool   TimingSendPar_No  = false;
   const bool   PredictParPos_No  = false;
   const bool   JustCountNPar_No  = false;
#  ifdef LOAD_BALANCE
   const bool   SibBufPatch       = true;
   const bool   FaSibBufPatch     = true;
#  else
   const bool   SibBufPatch       = NULL_BOOL;
   const bool   FaSibBufPatch     = NULL_BOOL;
#  endif
#  endif // #ifdef PARTICLE


// collect local data
   for (int lv=0; lv<NLEVEL; lv++)
   {
//    initialize the particle density array (rho_ext) and collect particles to the target level
#     ifdef PARTICLE
      Prepare_PatchData_InitParticleDensityArray( lv );

      Par_CollectParticle2OneLevel( lv, _PAR_MASS|_PAR_POSX|_PAR_POSY|_PAR_POSZ|_PAR_TYPE, PredictParPos_No, NULL_REAL,
                                    SibBufPatch, FaSibBufPatch, JustCountNPar_No, TimingSendPar_No );
#     endif

//    get the total density on grids
      real (*TotalDens)[PS1][PS1][PS1] = new real [ amr->NPatchComma[lv][1] ][PS1][PS1][PS1];
      int   *PID0List                  = new int  [ amr->NPatchComma[lv][1]/8 ];

      for (int PID0=0, t=0; PID0<amr->NPatchComma[lv][1]; PID0+=8, t++)    PID0List[t] = PID0;

      Prepare_PatchData( lv, Time[lv], TotalDens[0][0][0], NULL, 0, amr->NPatchComma[lv][1]/8, PID0List, DensMode, _NONE,
                         OPT__RHO_INT_SCHEME, INT_NONE, UNIT_PATCH, NSIDE_00, IntPhase_No, OPT__BC_FLU, BC_POT_NONE,
                         MinDens_No, MinPres_No, MinTemp_No, MinEntr_No, DE_Consistency_No );

      delete [] PID0List;


//    free memory for collecting particles from other ranks and levels, and free density arrays with ghost zones (rho_ext)
#     ifdef PARTICLE
      Par_CollectParticle2OneLevel_FreeMemory( lv, SibBufPatch, FaSibBufPatch );

      Prepare_PatchData_FreeParticleDensityArray( lv );
#     endif

      for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
      {
//       skip non-leaf patches
         if ( amr->patch[0][lv][PID]->son != -1 )  continue;

         for (int k=0; k<PS1; k++)  {  const double z = amr->patch[0][lv][PID]->EdgeL[2] + (k+0.5)*amr->dh[lv];
         for (int j=0; j<PS1; j++)  {  const double y = amr->patch[0][lv][PID]->EdgeL[1] + (j+0.5)*amr->dh[lv];
         for (int i=0; i<PS1; i++)  {  const double x = amr->patch[0][lv][PID]->EdgeL[0] + (i+0.5)*amr->dh[lv];

//          dens = amr->patch[ amr->FluSg[lv] ][lv][PID]->fluid[DENS][k][j][i];
            dens = TotalDens[PID][k][j][i];
            pote = amr->patch[ amr->PotSg[lv] ][lv][PID]->pot[k][j][i];

            if ( dens > max_dens_loc )
            {
               max_dens_loc        = dens;
               max_dens_pos_loc[0] = x;
               max_dens_pos_loc[1] = y;
               max_dens_pos_loc[2] = z;
            }

            if ( pote < min_pote_loc )
            {
               min_pote_loc        = pote;
               min_pote_pos_loc[0] = x;
               min_pote_pos_loc[1] = y;
               min_pote_pos_loc[2] = z;
            }
         }}}
      } // for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)

      delete [] TotalDens;
   } // for (int lv=0; lv<NLEVEL; lv++)


// gather data to the root rank
   send[0] = max_dens_loc;
   send[1] = max_dens_pos_loc[0];
   send[2] = max_dens_pos_loc[1];
   send[3] = max_dens_pos_loc[2];
   send[4] = min_pote_loc;
   send[5] = min_pote_pos_loc[0];
   send[6] = min_pote_pos_loc[1];
   send[7] = min_pote_pos_loc[2];

   MPI_Gather( send, CountMPI, MPI_DOUBLE, recv[0], CountMPI, MPI_DOUBLE, 0, MPI_COMM_WORLD );


// record the maximum density and center coordinates
   double max_dens      = -__DBL_MAX__;
   double min_pote      = +__DBL_MAX__;
   int    max_dens_rank = -1;
   int    min_pote_rank = -1;

   if ( MPI_Rank == 0 )
   {
      for (int r=0; r<MPI_NRank; r++)
      {
         if ( recv[r][0] > max_dens )
         {
            max_dens      = recv[r][0];
            max_dens_rank = r;
         }

         if ( recv[r][4] < min_pote )
         {
            min_pote      = recv[r][4];
            min_pote_rank = r;
         }
      }

      if ( max_dens_rank < 0  ||  max_dens_rank >= MPI_NRank )
         Aux_Error( ERROR_INFO, "incorrect max_dens_rank (%d) !!\n", max_dens_rank );

      if ( min_pote_rank < 0  ||  min_pote_rank >= MPI_NRank )
         Aux_Error( ERROR_INFO, "incorrect min_pote_rank (%d) !!\n", min_pote_rank );

      static bool FirstTime = true;

      if ( FirstTime )
      {
         if ( Aux_CheckFileExist(filename_center) )
            Aux_Message( stderr, "WARNING : file \"%s\" already exists !!\n", filename_center );
         else
         {
            FILE *file_center = fopen( filename_center, "w" );
            fprintf( file_center, "#%19s  %10s  %14s  %14s  %14s  %14s  %14s  %14s  %14s  %14s  %10s  %14s  %14s  %14s\n",
                     "Time", "Step", "Dens", "Dens_x", "Dens_y", "Dens_z", "Pote", "Pote_x", "Pote_y", "Pote_z",
                     "NIter", "CM_x", "CM_y", "CM_z" );
            fclose( file_center );
         }

         FirstTime = false;
      }

      FILE *file_center = fopen( filename_center, "a" );
      fprintf( file_center, "%20.14e  %10ld  %14.7e  %14.7e  %14.7e  %14.7e  %14.7e  %14.7e  %14.7e  %14.7e",
               Time[0], Step, recv[max_dens_rank][0], recv[max_dens_rank][1], recv[max_dens_rank][2], recv[max_dens_rank][3],
                              recv[min_pote_rank][4], recv[min_pote_rank][5], recv[min_pote_rank][6], recv[min_pote_rank][7] );
      fclose( file_center );
   } // if ( MPI_Rank == 0 )


// compute the center of mass until convergence
   const double CM_MaxR    = SQRT( SQR(amr->BoxSize[0]) + SQR(amr->BoxSize[1]) + SQR(amr->BoxSize[2]) ); // maximum radius for determining CM
   const double CM_TolErrR = amr->dh[0];                      // maximum allowed errors for determining CM
   const double TolErrR2 = SQR( CM_TolErrR );
   const int    NIterMax = 10;

   double dR2, CM_Old[3], CM_New[3];
   int NIter = 0;

// set an initial guess by the peak density position
   if ( MPI_Rank == 0 )
      for (int d=0; d<3; d++)    CM_Old[d] = recv[max_dens_rank][1+d];

   MPI_Bcast( CM_Old, 3, MPI_DOUBLE, 0, MPI_COMM_WORLD );

   while ( true )
   {
      Aux_FindCenterOfMass( CM_Old, CM_New, CM_MaxR );

      dR2 = SQR( CM_Old[0] - CM_New[0] )
          + SQR( CM_Old[1] - CM_New[1] )
          + SQR( CM_Old[2] - CM_New[2] );
      NIter ++;

      if ( dR2 <= TolErrR2  ||  NIter >= NIterMax )
         break;
      else
         memcpy( CM_Old, CM_New, sizeof(double)*3 );
   }

   if ( MPI_Rank == 0 )
   {
      if ( dR2 > TolErrR2 )
         Aux_Message( stderr, "WARNING : dR (%13.7e) > CM_TolErrR (%13.7e) !!\n", sqrt(dR2), CM_TolErrR );

      FILE *file_center = fopen( filename_center, "a" );
      fprintf( file_center, "  %10d  %14.7e  %14.7e  %14.7e\n", NIter, CM_New[0], CM_New[1], CM_New[2] );
      fclose( file_center );
   }


   delete [] recv;

} // FUNCTION : Aux_Record_Center
