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
   real WavePatchVolume  = 0;
   real TotalVolume      = 0;

   for (int lv=0; lv<NLEVEL; lv++)
   {
      const real dv = CUBE( amr->dh[lv] );
      const NPatch  = NPatchTotal[lv];

      if ( amr->use_wave_flag[lv] ) {
         WavePatchCount   += NPatch;
         WavePatchVolume  += NPatch * dv;
      }

      TotalVolume += NPatch * dv;
   }

   // output to file
   FILE *File = fopen( FileName, "a" );
   fprintf( File, "Time = %13.7e,  Step = %7ld,  NPatch = %10ld,  WavePatchFrac = %6.2f,  WaveVolFrac = %6.2f\n\n", Time[0], Step, NPatchAllLevel, WavePatchCount / NPatchAllLevel, WaveVolume / TotalVolume );
   fclose( File );

} // FUNCTION : ELBDM_Aux_Record_Hybrid

#endif // # if ( ELBDM_SCHEME == ELBDM_HYBRID )