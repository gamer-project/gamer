#include "Copyright.h"
#include "GAMER.h"




//-------------------------------------------------------------------------------------------------------
// Function    :  Aux_RecordPerformance
// Description :  Record the code performance
//
// Note        :  The code performance is defined as "total number of cell updates per second" 
//                --> Note that for the individual time-step integration cells at higher levels will be updated
//                    more frequently
//                --> The total number of cell updates recorded here for individual time-step integration is
//                    only approximate since it can change during one global time-step
// 
// Parameter   :  ElapsedTime : Elapsed time of the current global step
//-------------------------------------------------------------------------------------------------------
void Aux_RecordPerformance( const double ElapsedTime )
{

   const char FileName[] = "Record__Performance";
   static bool FirstTime = true;


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
         fprintf( File_Record, "#%13s%14s%3s%14s%14s%14s%14s%14s%14s\n", 
                  "Time", "Step", "", "dt", "NCell", "NUpdate", "ElapsedTime", "Perf_Overall", "Perf_per_GPU" );
#        else
         fprintf( File_Record, "#%13s%14s%3s%14s%14s%14s%14s%14s%14s\n", 
                  "Time", "Step", "", "dt", "NCell", "NUpdate", "ElapsedTime", "Perf_Overall", "Perf_per_CPU" );
#        endif
         fclose( File_Record );
      }


//    count the total number of cells and cell updates
      long NCell=0, NUpdate=0, NCellThisLevel;

      for (int lv=0; lv<NLEVEL; lv++)
      {
         NCellThisLevel = NPatchTotal[lv]*CUBE( PATCH_SIZE );

         NCell   += NCellThisLevel;
#        ifdef INDIVIDUAL_TIMESTEP
         NUpdate += NCellThisLevel*(1<<(lv+1));
#        else
         NUpdate += NCellThisLevel;
#        endif
      }


//    record performance
      const double NUpdate_PerSec         = NUpdate/ElapsedTime;
#     ifdef GPU
      const double NUpdate_PerSec_PerRank = NUpdate_PerSec/MPI_NRank;
#     else
      const double NUpdate_PerSec_PerRank = NUpdate_PerSec/MPI_NRank/OMP_NTHREAD;
#     endif

      FILE *File_Record = fopen( FileName, "a" );
      fprintf( File_Record, "%14.7e%14ld%3s%14.7e%14.2e%14.2e%14.2e%14.2e%14.2e\n", 
               Time[0], Step, "", dTime_Base, (double)NCell, (double)NUpdate, ElapsedTime, NUpdate_PerSec, 
               NUpdate_PerSec_PerRank );
      fclose( File_Record );

   } // if ( MPI_Rank == 0 )

} // FUNCTION : Aux_RecordPerformance


