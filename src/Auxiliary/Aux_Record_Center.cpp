#include "GAMER.h"
#include "TestProb.h"




//-------------------------------------------------------------------------------------------------------
// Function    :  Aux_Record_Center
// Description :  Record various center coordinates
//
// Note        :  1. Invoked by main()
//                2. Enabled by the runtime option "OPT__RECORD_CENTER"
//                3. This function will be called both during the program initialization and after each full update
//                4. It will record the position of maximum density, minimum potential, and center of mass
//                5. Output filename is fixed to "Record__Center"
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Aux_Record_Center()
{

   const char filename_center  [] = "Record__Center";

   // Maximum fluid density in HYDRO/ELBDM
   Extrema_t Max_Dens;
   Max_Dens.Field     = _DENS;
   Max_Dens.Radius    = HUGE_NUMBER; // entire domain
   Max_Dens.Center[0] = amr->BoxCenter[0];
   Max_Dens.Center[1] = amr->BoxCenter[1];
   Max_Dens.Center[2] = amr->BoxCenter[2];

   Aux_FindExtrema( &Max_Dens, EXTREMA_MAX, 0, TOP_LEVEL, PATCH_LEAF );

#  ifdef PARTICLE
   // Maximum particle density
   Extrema_t Max_ParDens;
   Max_ParDens.Field     = _PAR_DENS;
   Max_ParDens.Radius    = HUGE_NUMBER; // entire domain
   Max_ParDens.Center[0] = amr->BoxCenter[0];
   Max_ParDens.Center[1] = amr->BoxCenter[1];
   Max_ParDens.Center[2] = amr->BoxCenter[2];

   Aux_FindExtrema_ParDens( &Max_ParDens, EXTREMA_MAX, 0, TOP_LEVEL, PATCH_LEAF );

   // Maximun total density including fluid density and particle density
   Extrema_t Max_TotDens;
   Max_TotDens.Field     = _TOTAL_DENS;
   Max_TotDens.Radius    = HUGE_NUMBER; // entire domain
   Max_TotDens.Center[0] = amr->BoxCenter[0];
   Max_TotDens.Center[1] = amr->BoxCenter[1];
   Max_TotDens.Center[2] = amr->BoxCenter[2];

   Aux_FindExtrema_ParDens( &Max_TotDens, EXTREMA_MAX, 0, TOP_LEVEL, PATCH_LEAF );
#  endif

   // Minimun gravitational potential
   Extrema_t Min_Pote;
   Min_Pote.Field     = _POTE;
   Min_Pote.Radius    = HUGE_NUMBER; // entire domain
   Min_Pote.Center[0] = amr->BoxCenter[0];
   Min_Pote.Center[1] = amr->BoxCenter[1];
   Min_Pote.Center[2] = amr->BoxCenter[2];

   Aux_FindExtrema_ParDens( &Min_Pote, EXTREMA_MIN, 0, TOP_LEVEL, PATCH_LEAF );


// compute the center of mass until convergence
   const double TolErrR2 = SQR( COM_TOL_ERR_R );

   double dR2, CoM_Old[3], CoM_New[3];
   int NIter = 0;

// set an initial guess by the peak density position or the user-specified center
   if ( MPI_Rank == 0 )
   {
      if ( COM_CEN_X < 0  ||  COM_CEN_Y < 0  ||  COM_CEN_Z < 0 )
      {
#        ifdef PARTICLE
         for (int d=0; d<3; d++) CoM_Old[d] = Max_TotDens.Coord[d];
#        else
         for (int d=0; d<3; d++) CoM_Old[d] = Max_Dens.Coord[d];
#        endif
      }
      else
      {
         CoM_Old[0] = COM_CEN_X;
         CoM_Old[1] = COM_CEN_Y;
         CoM_Old[2] = COM_CEN_Z;
      }
   }

   MPI_Bcast( CoM_Old, 3, MPI_DOUBLE, 0, MPI_COMM_WORLD );

   while ( true )
   {
      Aux_FindWeightedCenter( CoM_New, CoM_Old, COM_MAX_R, COM_MIN_RHO, _TOTAL_DENS );

      dR2 = SQR( CoM_Old[0] - CoM_New[0] )
          + SQR( CoM_Old[1] - CoM_New[1] )
          + SQR( CoM_Old[2] - CoM_New[2] );
      NIter++;

      if ( dR2 <= TolErrR2  ||  NIter >= COM_N_ITER_MAX )
         break;
      else
         memcpy( CoM_Old, CoM_New, sizeof(double)*3 );
   }

   if ( MPI_Rank == 0 )
   {
      if ( dR2 > TolErrR2 )
         Aux_Message( stderr, "WARNING : dR (%13.7e) > COM_TOL_ERR_R (%13.7e) !!\n", sqrt(dR2), COM_TOL_ERR_R );
   }


   // Output the center to file
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
            fprintf( file_center, "#%19s  %10s", "Time", "Step" );
            fprintf( file_center, "  %14s  %14s  %14s  %14s",
                     "MaxDens", "MaxDens_x", "MaxDens_y", "MaxDens_z" );
#           ifdef PARTICLE
            fprintf( file_center, "  %14s  %14s  %14s  %14s",
                     "MaxParDens", "MaxParDens_x", "MaxParDens_y", "MaxParDens_z" );
            fprintf( file_center, "  %14s  %14s  %14s  %14s",
                     "MaxTotalDens", "MaxTotalDens_x", "MaxTotalDens_y", "MaxTotalDens_z" );
#           endif
            fprintf( file_center, "  %14s  %14s  %14s  %14s",
                     "MinPote", "MinPote_x", "MinPote_y", "MinPote_z" );
            fprintf( file_center, "  %10s  %14s  %14s  %14s",
                     "NIter", "CoM_x", "CoM_y", "CoM_z" );
            fprintf( file_center, "\n" );
            fclose( file_center );
         }

         FirstTime = false;
      }

      FILE *file_center = fopen( filename_center, "a" );
      fprintf( file_center, "%20.14e  %10ld", Time[0], Step );
      fprintf( file_center, "  %14.7e  %14.7e  %14.7e  %14.7e",
               Max_Dens.Value, Max_Dens.Coord[0], Max_Dens.Coord[1], Max_Dens.Coord[2] );
#     ifdef PARTICLE
      fprintf( file_center, "  %14.7e  %14.7e  %14.7e  %14.7e",
               Max_ParDens.Value, Max_ParDens.Coord[0], Max_ParDens.Coord[1], Max_ParDens.Coord[2] );
      fprintf( file_center, "  %14.7e  %14.7e  %14.7e  %14.7e",
               Max_TotDens.Value, Max_TotDens.Coord[0], Max_TotDens.Coord[1], Max_TotDens.Coord[2] );
#     endif
      fprintf( file_center, "  %14.7e  %14.7e  %14.7e  %14.7e",
               Min_Pote.Value, Min_Pote.Coord[0], Min_Pote.Coord[1], Min_Pote.Coord[2] );
      fprintf( file_center, "  %10d  %14.7e  %14.7e  %14.7e",
               NIter, CoM_New[0], CoM_New[1], CoM_New[2] );
      fprintf( file_center, "\n" );
      fclose( file_center );
   }

} // FUNCTION : Aux_Record_Center
