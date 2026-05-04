#include "GAMER.h"




//-------------------------------------------------------------------------------------------------------
// Function    :  Aux_Record_Center
// Description :  Record various center coordinates
//
// Note        :  1. Invoked by main()
//                2. Enabled by the runtime option "OPT__RECORD_CENTER"
//                3. This function will be called both during the program initialization and after each global step
//                4. It will record the position of maximum density, minimum potential, and center of mass
//                5. Output filename is fixed to "Record__Center"
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Aux_Record_Center()
{

   static bool FirstTime = true;
   char FileName[2*MAX_STRING];
   sprintf( FileName, "%s/Record__Center", OUTPUT_DIR );


// 1. Maximum fluid density in HYDRO/ELBDM
   Extrema_t Max_Dens;
   Max_Dens.Field     = _DENS;
   Max_Dens.Radius    = __FLT_MAX__; // entire domain
   Max_Dens.Center[0] = amr->BoxCenter[0];
   Max_Dens.Center[1] = amr->BoxCenter[1];
   Max_Dens.Center[2] = amr->BoxCenter[2];

   Aux_FindExtrema( &Max_Dens, EXTREMA_MAX, 0, TOP_LEVEL, PATCH_LEAF );


#  ifdef PARTICLE
// 2. Maximum particle density
   Extrema_t Max_ParDens;
   Max_ParDens.Field     = _PAR_DENS;
   Max_ParDens.Radius    = __FLT_MAX__; // entire domain
   Max_ParDens.Center[0] = amr->BoxCenter[0];
   Max_ParDens.Center[1] = amr->BoxCenter[1];
   Max_ParDens.Center[2] = amr->BoxCenter[2];

   Aux_FindExtrema( &Max_ParDens, EXTREMA_MAX, 0, TOP_LEVEL, PATCH_LEAF );


// 3. Maximum total density including fluid density and particle density
   Extrema_t Max_TotDens;
   Max_TotDens.Field     = _TOTAL_DENS;
   Max_TotDens.Radius    = __FLT_MAX__; // entire domain
   Max_TotDens.Center[0] = amr->BoxCenter[0];
   Max_TotDens.Center[1] = amr->BoxCenter[1];
   Max_TotDens.Center[2] = amr->BoxCenter[2];

   Aux_FindExtrema( &Max_TotDens, EXTREMA_MAX, 0, TOP_LEVEL, PATCH_LEAF );
#  endif


#  ifdef GRAVITY
// 4. Minimum gravitational potential
   Extrema_t Min_Pote;
   Min_Pote.Field     = _POTE;
   Min_Pote.Radius    = __FLT_MAX__; // entire domain
   Min_Pote.Center[0] = amr->BoxCenter[0];
   Min_Pote.Center[1] = amr->BoxCenter[1];
   Min_Pote.Center[2] = amr->BoxCenter[2];

   Aux_FindExtrema( &Min_Pote, EXTREMA_MIN, 0, TOP_LEVEL, PATCH_LEAF );
#  endif


// 5. Center of mass for the total density field
// set an initial guess by the peak density position or the user-specified center
   double CoM_ref[3];
   if ( COM_CEN_X < 0.0  ||  COM_CEN_Y < 0.0  ||  COM_CEN_Z < 0.0 )
   {
#     ifdef PARTICLE
      for (int d=0; d<3; d++) CoM_ref[d] = Max_TotDens.Coord[d];
#     else
      for (int d=0; d<3; d++) CoM_ref[d] = Max_Dens.Coord[d];
#     endif
   }
   else
   {
      CoM_ref[0] = COM_CEN_X;
      CoM_ref[1] = COM_CEN_Y;
      CoM_ref[2] = COM_CEN_Z;
   }

// find the center of mass
   double CoM_Coord[3];
   double FinaldR;
   int    FinalNIter;
   Aux_FindWeightedAverageCenter( CoM_Coord, CoM_ref, COM_MAX_R, COM_MIN_RHO, _TOTAL_DENS, COM_TOLERR_R, COM_MAX_ITER, &FinaldR, &FinalNIter );


// Output the center to file
   if ( MPI_Rank == 0 )
   {
//    Output the header
      if ( FirstTime )
      {
         if ( Aux_CheckFileExist(FileName) )
            Aux_Message( stderr, "WARNING : file \"%s\" already exists !!\n", FileName );
         else
         {
            FILE *File = fopen( FileName, "w" );
            fprintf( File, "#%19s  %10s  %14s  %14s  %14s  %14s",
                           "Time", "Step", "MaxDens", "MaxDens_x", "MaxDens_y", "MaxDens_z" );

#           ifdef PARTICLE
            fprintf( File, "  %14s  %14s  %14s  %14s  %14s  %14s  %14s  %14s",
                           "MaxParDens", "MaxParDens_x", "MaxParDens_y", "MaxParDens_z",
                           "MaxTotalDens", "MaxTotalDens_x", "MaxTotalDens_y", "MaxTotalDens_z" );
#           endif

#           ifdef GRAVITY
            fprintf( File, "  %14s  %14s  %14s  %14s",
                           "MinPote", "MinPote_x", "MinPote_y", "MinPote_z" );
#           endif

            fprintf( File, "  %14s  %14s  %14s  %14s  %14s",
                           "Final_NIter", "Final_dR", "CoM_x", "CoM_y", "CoM_z" );

            fprintf( File, "\n" );
            fclose( File );
         }

         FirstTime = false;
      } // if ( FirstTime )

      FILE *File = fopen( FileName, "a" );
      fprintf( File, "%20.14e  %10ld  %14.7e  %14.7e  %14.7e  %14.7e",
                     Time[0], Step, Max_Dens.Value, Max_Dens.Coord[0], Max_Dens.Coord[1], Max_Dens.Coord[2] );

#     ifdef PARTICLE
      fprintf( File, "  %14.7e  %14.7e  %14.7e  %14.7e  %14.7e  %14.7e  %14.7e  %14.7e",
                     Max_ParDens.Value, Max_ParDens.Coord[0], Max_ParDens.Coord[1], Max_ParDens.Coord[2],
                     Max_TotDens.Value, Max_TotDens.Coord[0], Max_TotDens.Coord[1], Max_TotDens.Coord[2] );
#     endif

#     ifdef GRAVITY
      fprintf( File, "  %14.7e  %14.7e  %14.7e  %14.7e",
                     Min_Pote.Value, Min_Pote.Coord[0], Min_Pote.Coord[1], Min_Pote.Coord[2] );
#     endif

      fprintf( File, "  %14d  %14.7e  %14.7e  %14.7e  %14.7e",
                     FinalNIter, FinaldR, CoM_Coord[0], CoM_Coord[1], CoM_Coord[2] );

      fprintf( File, "\n" );
      fclose( File );
   } // if ( MPI_Rank == 0 )

} // FUNCTION : Aux_Record_Center
