#include "GAMER.h"
#include "TestProb.h"


#  ifdef PARTICLE
static void Aux_FindMaxDensWithPar( double *Value, double Coord[], const long DensField );
#  endif


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

#  ifdef PARTICLE
   double Max_ParDens_Value, Max_ParDens_Coord[3];
   double Max_TotDens_Value, Max_TotDens_Coord[3];

   Aux_FindMaxDensWithPar( &Max_ParDens_Value, Max_ParDens_Coord, _PAR_DENS   );
   Aux_FindMaxDensWithPar( &Max_TotDens_Value, Max_TotDens_Coord, _TOTAL_DENS );
#  endif

   Extrema_t Max_Dens;
   Max_Dens.Field     = _DENS;
   Max_Dens.Radius    = HUGE_NUMBER; // entire domain
   Max_Dens.Center[0] = amr->BoxCenter[0];
   Max_Dens.Center[1] = amr->BoxCenter[1];
   Max_Dens.Center[2] = amr->BoxCenter[2];

   Aux_FindExtrema( &Max_Dens, EXTREMA_MAX, 0, TOP_LEVEL, PATCH_LEAF );

   Extrema_t Min_Pote;
   Min_Pote.Field     = _POTE;
   Min_Pote.Radius    = HUGE_NUMBER; // entire domain
   Min_Pote.Center[0] = amr->BoxCenter[0];
   Min_Pote.Center[1] = amr->BoxCenter[1];
   Min_Pote.Center[2] = amr->BoxCenter[2];

   Aux_FindExtrema( &Min_Pote, EXTREMA_MIN, 0, TOP_LEVEL, PATCH_LEAF );


// compute the center of mass until convergence
   const double CM_MaxR    = SQRT( SQR(amr->BoxSize[0]) + SQR(amr->BoxSize[1]) + SQR(amr->BoxSize[2]) ); // maximum radius for determining CM
   const double CM_MinRho  = 0.0;
   const double CM_TolErrR = amr->dh[0];                      // maximum allowed errors for determining CM
   const double TolErrR2 = SQR( CM_TolErrR );
   const int    NIterMax = 10;

   double dR2, CM_Old[3], CM_New[3];
   int NIter = 0;

// set an initial guess by the peak density position
   if ( MPI_Rank == 0 )
      for (int d=0; d<3; d++)    CM_Old[d] = Max_Dens.Coord[d];

   MPI_Bcast( CM_Old, 3, MPI_DOUBLE, 0, MPI_COMM_WORLD );

   while ( true )
   {
      Aux_FindWeightedCenter( CM_New, CM_Old, CM_MaxR, CM_MinRho, 1, _TOTAL_DENS );

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
   }

   // Output
   if ( MPI_Rank == 0 )
   {

      static bool FirstTime = true;

      if ( FirstTime )
      {
         if ( Aux_CheckFileExist(filename_center) )
            Aux_Message( stderr, "WARNING : file \"%s\" already exists !!\n", filename_center );
         else
         {
            FILE *file_center = fopen( filename_center, "w" );
            fprintf( file_center, "#%19s  %10s  %14s  %14s  %14s  %14s",
                     "Time", "Step", "Dens", "Dens_x", "Dens_y", "Dens_z");
            fprintf( file_center, "  %14s  %14s  %14s  %14s",
                     "Pote", "Pote_x", "Pote_y", "Pote_z");
#           ifdef PARTICLE
            fprintf( file_center, "  %14s  %14s  %14s  %14s",
                     "ParDens", "ParDens_x", "ParDens_y", "ParDens_z");
            fprintf( file_center, "  %14s  %14s  %14s  %14s",
                     "TotalDens", "TotalDens_x", "TotalDens_y", "TotalDens_z");
#           endif
            fprintf( file_center, "  %10s  %14s  %14s  %14s\n",
                     "NIter", "CM_x", "CM_y", "CM_z" );
            fclose( file_center );
         }

         FirstTime = false;
      }

      FILE *file_center = fopen( filename_center, "a" );
      fprintf( file_center, "%20.14e  %10ld  %14.7e  %14.7e  %14.7e  %14.7e",
               Time[0], Step, Max_Dens.Value, Max_Dens.Coord[0], Max_Dens.Coord[1], Max_Dens.Coord[2]);
      fprintf( file_center, "  %14.7e  %14.7e  %14.7e  %14.7e",
                              Min_Pote.Value, Min_Pote.Coord[0], Min_Pote.Coord[1], Min_Pote.Coord[2]);
#     ifdef PARTICLE
      fprintf( file_center, "  %14.7e  %14.7e  %14.7e  %14.7e",
                              Max_ParDens_Value, Max_ParDens_Coord[0], Max_ParDens_Coord[1], Max_ParDens_Coord[2] );
      fprintf( file_center, "  %14.7e  %14.7e  %14.7e  %14.7e",
                              Max_TotDens_Value, Max_TotDens_Coord[0], Max_TotDens_Coord[1], Max_TotDens_Coord[2] );
#     endif
      fprintf( file_center, "  %10d  %14.7e  %14.7e  %14.7e\n",
                              NIter, CM_New[0], CM_New[1], CM_New[2] );
      fclose( file_center );
   }

} // FUNCTION : Aux_Record_Center



#     ifdef PARTICLE
//-------------------------------------------------------------------------------------------------------
// Function    :  Aux_FindMaxDensWithPar
// Description :  Find the value and location of the maximum value of a total density or particle density
//                within the whole box
//
// Note        :
//
//
// Parameter   : Value     : Maximum value of density
//               Coord     : Coordinates of the maximum density
//               DensField : Target density to be calculated (_PAR_DENS or _TOTAL_DENS)
//
// Return      :  Value, Coord[]
//-------------------------------------------------------------------------------------------------------
void Aux_FindMaxDensWithPar( double *Value, double Coord[], const long DensField )
{

   const int  CountMPI            = 4; // value + 3 coord

   double target_dens, max_target_dens_loc=-__DBL_MAX__, max_target_dens_pos_loc[3];
   double send[CountMPI], (*recv)[CountMPI]=new double [MPI_NRank][CountMPI];

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
      real (*TargetDens)[PS1][PS1][PS1] = new real [ amr->NPatchComma[lv][1] ][PS1][PS1][PS1];
      int   *PID0List                  = new int  [ amr->NPatchComma[lv][1]/8 ];

      for (int PID0=0, t=0; PID0<amr->NPatchComma[lv][1]; PID0+=8, t++)    PID0List[t] = PID0;

      Prepare_PatchData( lv, Time[lv], TargetDens[0][0][0], NULL, 0, amr->NPatchComma[lv][1]/8, PID0List, DensField, _NONE,
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

            target_dens = TargetDens[PID][k][j][i];

            if ( target_dens > max_target_dens_loc )
            {
               max_target_dens_loc        = target_dens;
               max_target_dens_pos_loc[0] = x;
               max_target_dens_pos_loc[1] = y;
               max_target_dens_pos_loc[2] = z;
            }

         }}}
      } // for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)

      delete [] TargetDens;
   } // for (int lv=0; lv<NLEVEL; lv++)


// gather data to the root rank
   send[0] = max_target_dens_loc;
   send[1] = max_target_dens_pos_loc[0];
   send[2] = max_target_dens_pos_loc[1];
   send[3] = max_target_dens_pos_loc[2];

   MPI_Gather( send, CountMPI, MPI_DOUBLE, recv[0], CountMPI, MPI_DOUBLE, 0, MPI_COMM_WORLD );


// record the maximum density and center coordinates
   double max_target_dens      = -__DBL_MAX__;
   int    max_target_dens_rank = -1;

   if ( MPI_Rank == 0 )
   {
      for (int r=0; r<MPI_NRank; r++)
      {
         if ( recv[r][0] > max_target_dens )
         {
            max_target_dens      = recv[r][0];
            max_target_dens_rank = r;
         }
      }

      if ( max_target_dens_rank < 0  ||  max_target_dens_rank >= MPI_NRank )
         Aux_Error( ERROR_INFO, "incorrect max_target_dens_rank (%d) !!\n", max_target_dens_rank );
   } // if ( MPI_Rank == 0 )

   *Value = recv[max_target_dens_rank][0];
   for(int d=0; d<3; d++) Coord[d] = recv[max_target_dens_rank][d+1];

   delete [] recv;

}
#     endif
