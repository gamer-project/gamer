#include "GAMER.h"




//-------------------------------------------------------------------------------------------------------
// Function    :  Aux_Record_Performance
// Description :  Record the code performance
//
// Note        :  1. The code performance is defined as "total number of cell updates per second"
//                   --> Note that for the individual time-step integration cells at higher levels will be updated
//                       more frequently
//                   --> The total number of cell and particle updates recorded here for the individual time-step
//                       integration is only approximate since the number of patches at each level may change
//                       during one global time-step
//                2. When PARTICLE is on, this routine also records the "total number of particle updates per second"
//
// Parameter   :  ElapsedTime : Elapsed time of the current global step
//-------------------------------------------------------------------------------------------------------
void Aux_Record_Performance( const double ElapsedTime )
{

   const char FileName[] = "Record__Performance";
   static bool FirstTime = true;


// get the total number of active particles in each rank
#  ifdef PARTICLE
   long NPar_Lv_AllRank[NLEVEL];
   MPI_Reduce( amr->Par->NPar_Lv, NPar_Lv_AllRank, NLEVEL, MPI_LONG, MPI_SUM, 0, MPI_COMM_WORLD );
#  endif


// only rank 0 needs to take a note
   if ( MPI_Rank == 0 )
   {
//    header
      if ( FirstTime )
      {
         if ( Aux_CheckFileExist(FileName) )
            Aux_Message( stderr, "WARNING : file \"%s\" already exists !!\n", FileName );

         FirstTime = false;

         FILE *File_Record = fopen( FileName, "a" );

         fprintf( File_Record, "#%13s%14s%3s%14s%14s%14s%14s%14s%14s",
                  "Time", "Step", "", "dt", "NCell", "NUpdate_Cell", "ElapsedTime", "Perf_Overall", "Perf_PerRank" );
#        ifdef PARTICLE
         fprintf( File_Record, "%14s%14s%17s%17s",
                  "NParticle", "NUpdate_Par", "ParPerf_Overall", "ParPerf_PerRank" );
#        endif

         for (int lv=0; lv<NLEVEL; lv++)
         {
            char tmp[MAX_STRING];
            sprintf( tmp, "NUpdate_Lv%d", lv );
            fprintf( File_Record, "%14s", tmp );
         }

         fprintf( File_Record, "\n" );
         fclose( File_Record );
      } // if ( FirstTime )


//    count the total number of cells, cell updates, and particle updates
      long NCell=0, NUpdateCell=0, NCellThisLevel;
#     ifdef PARTICLE
      long NUpdatePar=0;
#     endif

      for (int lv=0; lv<NLEVEL; lv++)
      {
         NCellThisLevel = (long)NPatchTotal[lv]*CUBE( PATCH_SIZE );
         NCell         += NCellThisLevel;
         NUpdateCell   += NCellThisLevel     *amr->NUpdateLv[lv];
#        ifdef PARTICLE
         NUpdatePar    += NPar_Lv_AllRank[lv]*amr->NUpdateLv[lv];
#        endif
      }


//    record performance
      const double NUpdateCell_PerSec         = NUpdateCell/ElapsedTime;
#     ifdef GPU
      const double NUpdateCell_PerSec_PerRank = NUpdateCell_PerSec/MPI_NRank;
#     else
      const double NUpdateCell_PerSec_PerRank = NUpdateCell_PerSec/MPI_NRank/OMP_NTHREAD;
#     endif

#     ifdef PARTICLE
      const double NUpdatePar_PerSec          = NUpdatePar/ElapsedTime;
#     ifdef GPU
      const double NUpdatePar_PerSec_PerRank  = NUpdatePar_PerSec/MPI_NRank;
#     else
      const double NUpdatePar_PerSec_PerRank  = NUpdatePar_PerSec/MPI_NRank/OMP_NTHREAD;
#     endif
#     endif

      FILE *File_Record = fopen( FileName, "a" );

      fprintf( File_Record, "%14.7e%14ld%3s%14.7e%14.2e%14.2e%14.2e%14.2e%14.2e",
               Time[0], Step, "", dTime_Base, (double)NCell, (double)NUpdateCell, ElapsedTime, NUpdateCell_PerSec,
               NUpdateCell_PerSec_PerRank );

#     ifdef PARTICLE
      fprintf( File_Record, "%14.2e%14.2e%17.2e%17.2e",
               (double)amr->Par->NPar_Active_AllRank, (double)NUpdatePar, NUpdatePar_PerSec, NUpdatePar_PerSec_PerRank );
#     endif

      for (int lv=0; lv<NLEVEL; lv++)
      fprintf( File_Record, "%14ld", amr->NUpdateLv[lv] );

      fprintf( File_Record, "\n" );

      fclose( File_Record );

   } // if ( MPI_Rank == 0 )

} // FUNCTION : Aux_Record_Performance


