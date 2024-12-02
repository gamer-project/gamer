#include "GAMER.h"




//-------------------------------------------------------------------------------------------------------
// Function    :  Aux_PauseManually
// Description :  Pause a simulation when detecting a file named "PAUSE_GAMER_PAUSE"
//
// Note        :  To resume a simulation, simply delete the above file
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Aux_PauseManually()
{

   const char FileName[] = "PAUSE_GAMER_PAUSE";
   const int  PauseSec   = 10;   // sleep for PauseSec seconds before checking the target file again
                                 // --> set rather arbitrarily here


// enforce NFS to flush the file information
   system( "ls > /dev/null" );


   bool FirstTime   = true;
   int  Pause_local = true;

   while ( Pause_local )
   {
//    check the target file
      Pause_local = Aux_CheckFileExist( FileName );

//    pause a simulation as long as ONE process has detected the target file
      int Pause_global;
      MPI_Allreduce( &Pause_local, &Pause_global, 1, MPI_INT, MPI_BOR, MPI_COMM_WORLD );

//    pause
      if ( Pause_global )
      {
         if ( MPI_Rank == 0  &&  FirstTime )
            Aux_Message( stdout, "\nDetecting file \"PAUSE_GAMER_PAUSE\" --> simulation is paused until this file is deleted ...\n\n" );

         sleep( PauseSec );

         FirstTime = false;
      }
   } // while ( Pause_local )

} // FUNCTION : Aux_PauseManually


