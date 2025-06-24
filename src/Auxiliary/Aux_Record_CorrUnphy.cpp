#include "GAMER.h"




//-------------------------------------------------------------------------------------------------------
// Function    :  Aux_Record_CorrUnphy
// Description :  Record the number of cells with unphysical results and are corrected by Flu_Close()->CorrectUnphysical()
//
// Note        :  1. These cells are corrected by either "OPT__1ST_FLUX_CORR" or "MIN_DENS/PRES"
//                2. The number of corrected cells is recorded in NCorrUnphy in the file "Flu_Close.cpp"
//                3. The total number of cell updates recorded here for the individual time-step integration is
//                   only approximate since the number of patches at each level may change during one global time-step
//-------------------------------------------------------------------------------------------------------
void Aux_Record_CorrUnphy()
{

   static bool FirstTime = true;
   char FileName[2*MAX_STRING];
   sprintf( FileName, "%s/Record__NCorrUnphy", OUTPUT_DIR );

   long NCorrAllRank[NLEVEL];
   FILE *File = NULL;


// collect data from all ranks
   MPI_Reduce( NCorrUnphy, NCorrAllRank, NLEVEL, MPI_LONG, MPI_SUM, 0, MPI_COMM_WORLD );


// only rank 0 needs to take a note
   if ( MPI_Rank == 0 )
   {
//    header
      if ( FirstTime )
      {
         if ( Aux_CheckFileExist(FileName) )
            Aux_Message( stderr, "WARNING : file \"%s\" already exists !!\n", FileName );

         FirstTime = false;

         File = fopen( FileName, "a" );

         fprintf( File, "#%13s %9s %10s", "Time", "Step", "NCorrAllLv" );
         for (int lv=0; lv<NLEVEL; lv++)  fprintf( File, "%16s %2d ", "Level", lv );

         fprintf( File, "\n" );

         fclose( File );
      }


//    count the total number of cells and cell updates
      long   NCorrAllLv=0, NUpdate;
      double Frac;

      for (int lv=0; lv<NLEVEL; lv++)  NCorrAllLv += NCorrAllRank[lv];

      File = fopen( FileName, "a" );

      fprintf( File, "%14.7e %9ld %10ld", Time[0], Step, NCorrAllLv );

      for (int lv=0; lv<NLEVEL; lv++)
      {
         NUpdate = amr->NUpdateLv[lv]*NPatchTotal[lv]*CUBE( PATCH_SIZE );

         if ( NUpdate == 0 )     Frac = 0.0;
         else                    Frac = 100.0 * NCorrAllRank[lv] / NUpdate;

         fprintf( File, " %8ld(%8.2e%%)", NCorrAllRank[lv], Frac );
      }

      fprintf( File, "\n" );

      fclose( File );

   } // if ( MPI_Rank == 0 )


// reset the counter
   for (int lv=0; lv<NLEVEL; lv++)  NCorrUnphy[lv] = 0;

} // FUNCTION : Aux_Record_CorrUnphy


