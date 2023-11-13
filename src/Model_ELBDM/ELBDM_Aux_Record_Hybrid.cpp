#include "GAMER.h"


# if ( ELBDM_SCHEME == ELBDM_HYBRID )

//-------------------------------------------------------------------------------------------------------
// Function    :  ELBDM_Aux_Record_Hybrid
// Description :  Record ratio of number of wave patches to total patches and fraction of simulation volume using wave solver
//
//
// Return      :  Log file "Record__Hybrid"
//-------------------------------------------------------------------------------------------------------
void ELBDM_Aux_Record_Hybrid()
{

   if ( MPI_Rank != 0 ) return;

   static bool FirstTime = true;
   const char *FileName  = "Record__Hybrid";


   if ( FirstTime )
   {
      if ( Aux_CheckFileExist(FileName) )
         Aux_Message( stderr, "WARNING : file \"%s\" already exists !!\n", FileName );
   }


// compute number of wave patches on all levels
   long WavePatchCount   = 0;
   long TotalPatchCount  = 0;
   int  WaveLevel        = -1;


   for (int lv=0; lv<NLEVEL; lv++)
   {
      const real dv      = CUBE( amr->dh[lv] );
      const long NPatch  = NPatchTotal[lv];

      if ( amr->use_wave_flag[lv] ) {
//       store first level using wave scheme in WaveLevel
         if ( WaveLevel == -1 ) WaveLevel = lv;
         WavePatchCount   += NPatch;
      }

      TotalPatchCount += NPatch;
   }


   const real TotalVolume = NPatchTotal[0] * CUBE( amr->dh[0] );
   const real WaveVolume  = (WaveLevel == -1) ? 0.0 : NPatchTotal[WaveLevel] * CUBE( amr->dh[WaveLevel] );

   // output to file
   FILE *File = fopen( FileName, "a" );
   if ( WaveLevel == -1 )
      fprintf( File, "Time = %13.7e,  Step = %7ld,  NPatch = %10ld,  No wave scheme used\n", Time[0], Step, TotalPatchCount );
   else
      fprintf( File, "Time = %13.7e,  Step = %7ld,  NPatch = %10ld,  Wave scheme on lv %2d,  NWavePatch/NPatch = %6.2f,  WaveVolume/Volume = %6.2f\n", Time[0], Step, TotalPatchCount, WaveLevel, (real) WavePatchCount / (real) TotalPatchCount, WaveVolume / TotalVolume );
   fclose( File );

   if ( FirstTime )  FirstTime = false;

} // FUNCTION : ELBDM_Aux_Record_Hybrid

#endif // # if ( ELBDM_SCHEME == ELBDM_HYBRID )