#include "Copyright.h"
#include "GAMER.h"




//-------------------------------------------------------------------------------------------------------
// Function    :  Aux_Record_Performance
// Description :  Record the code performance
//
// Note        :  1. The code performance is defined as "total number of cell updates per second"
//                   --> Note that for the individual time-step integration cells at higher levels will be updated
//                       more frequently
//                   --> The total number of cell updates recorded here for individual time-step integration is
//                       only an approximate number since the grid structure can change during one global time-step
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

#        ifdef GPU
         fprintf( File_Record, "#%13s%14s%3s%14s%14s%14s%14s%14s%14s",
                  "Time", "Step", "", "dt", "NCell", "NUpdate", "ElapsedTime", "Perf_Overall", "Perf_per_GPU" );
#        ifdef PARTICLE
         fprintf( File_Record, "%14s%14s%17s%17s",
                  "NParticle", "ParNUpdate", "ParPerf_Overall", "ParPerf_per_GPU" );
#        endif

#        else // #ifdef GPU

         fprintf( File_Record, "#%13s%14s%3s%14s%14s%14s%14s%14s%14s",
                  "Time", "Step", "", "dt", "NCell", "NUpdate", "ElapsedTime", "Perf_Overall", "Perf_per_CPU" );
#        ifdef PARTICLE
         fprintf( File_Record, "%14s%14s%17s%17s",
                  "NParticle", "ParNUpdate", "ParPerf_Overall", "ParPerf_per_CPU" );
#        endif
#        endif // #ifdef GPU ... else ...

         fprintf( File_Record, "\n" );
         fclose( File_Record );
      } // if ( FirstTime )


//    count the total number of cells, cell updates, and particle updates
      long NCell=0, NUpdate=0, NCellThisLevel, LoadThisLv;
#     ifdef PARTICLE
      long ParNUpdate=0;
#     endif

      for (int lv=0; lv<NLEVEL; lv++)
      {
         NCellThisLevel = (long)NPatchTotal[lv]*CUBE( PATCH_SIZE );
         NCell         += NCellThisLevel;
         LoadThisLv     = ( OPT__DT_LEVEL == DT_LEVEL_DIFF_BY_2 ) ? ( 1<<(lv+1) ) : 1;
         NUpdate       += NCellThisLevel     *LoadThisLv;
#        ifdef PARTICLE
         ParNUpdate    += NPar_Lv_AllRank[lv]*LoadThisLv;
#        endif
      }


//    record performance
      const double NUpdate_PerSec            = NUpdate/ElapsedTime;
#     ifdef GPU
      const double NUpdate_PerSec_PerRank    = NUpdate_PerSec/MPI_NRank;
#     else
      const double NUpdate_PerSec_PerRank    = NUpdate_PerSec/MPI_NRank/OMP_NTHREAD;
#     endif

#     ifdef PARTICLE
      const double ParNUpdate_PerSec         = ParNUpdate/ElapsedTime;
#     ifdef GPU
      const double ParNUpdate_PerSec_PerRank = ParNUpdate_PerSec/MPI_NRank;
#     else
      const double ParNUpdate_PerSec_PerRank = ParNUpdate_PerSec/MPI_NRank/OMP_NTHREAD;
#     endif
#     endif

      FILE *File_Record = fopen( FileName, "a" );
      fprintf( File_Record, "%14.7e%14ld%3s%14.7e%14.2e%14.2e%14.2e%14.2e%14.2e",
               Time[0], Step, "", dTime_Base, (double)NCell, (double)NUpdate, ElapsedTime, NUpdate_PerSec,
               NUpdate_PerSec_PerRank );
#     ifdef PARTICLE
      fprintf( File_Record, "%14.2e%14.2e%17.2e%17.2e",
               (double)amr->Par->NPar_Active_AllRank, (double)ParNUpdate, ParNUpdate_PerSec, ParNUpdate_PerSec_PerRank );
#     endif
      fprintf( File_Record, "\n" );
      fclose( File_Record );

   } // if ( MPI_Rank == 0 )

} // FUNCTION : Aux_Record_Performance


