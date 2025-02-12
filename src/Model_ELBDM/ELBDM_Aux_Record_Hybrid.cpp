#include "GAMER.h"

#if ( ELBDM_SCHEME == ELBDM_HYBRID )




//-------------------------------------------------------------------------------------------------------
// Function    :  ELBDM_Aux_Record_Hybrid
// Description :  Record the ratio of number of wave patches to total patches and the fraction of simulation volume
//                using wave solver
//
// Return      :  Log file "Record__Hybrid"
//-------------------------------------------------------------------------------------------------------
void ELBDM_Aux_Record_Hybrid()
{

   if ( MPI_Rank != 0 )    return;


   static bool FirstTime = true;
   char FileName[2*MAX_STRING];
   sprintf( FileName, "%s/Record__Hybrid", OUTPUT_DIR );
   FILE *File = NULL;

   if ( FirstTime )
   {
      if ( Aux_CheckFileExist(FileName) )
         Aux_Message( stderr, "WARNING : file \"%s\" already exists !!\n", FileName );

      File = fopen( FileName, "a" );
      fprintf( File, "#%13s  %8s  %10s  %9s  %17s  %17s\n", "Time", "Step", "NPatch", "WaveLevel", "NWavePatch/NPatch", "WaveVolume/Volume" );
      fclose( File );
   }


// compute the number of wave patches on all levels
   long WavePatchCount  = 0;
   long TotalPatchCount = 0;
   int  WaveLevel       = -1;

   for (int lv=0; lv<NLEVEL; lv++)
   {
      const long NPatch = NPatchTotal[lv];

      if ( amr->use_wave_flag[lv] ) {
         if ( WaveLevel == -1 )  WaveLevel = lv;   // store the first level using wave scheme
         WavePatchCount += NPatch;
      }

      TotalPatchCount += NPatch;
   }


// output to file
   const double TotalVolume = NPatchTotal[0] * CUBE( amr->dh[0] );
   const double WaveVolume  = ( WaveLevel == -1 ) ? 0.0 : NPatchTotal[WaveLevel] * CUBE( amr->dh[WaveLevel] );

   File = fopen( FileName, "a" );
   fprintf( File, "%14.7e  %8ld  %10ld  %9d  %17.4f  %17.4f\n",
            Time[0], Step, TotalPatchCount, WaveLevel, (double)WavePatchCount/(double)TotalPatchCount, WaveVolume/TotalVolume );
   fclose( File );


   if ( FirstTime )  FirstTime = false;

} // FUNCTION : ELBDM_Aux_Record_Hybrid



#endif // #if ( ELBDM_SCHEME == ELBDM_HYBRID )
