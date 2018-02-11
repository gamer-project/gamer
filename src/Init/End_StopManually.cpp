#include "GAMER.h"




//-------------------------------------------------------------------------------------------------------
// Function    :  End_StopManually
// Description :  Terminate the program if the stop file named "STOP_GAMER_STOP" is found
//
// Parameter   :  Terminate_global : Boolean variable determining whether or not to terminate the run
//-------------------------------------------------------------------------------------------------------
void End_StopManually( int &Terminate_global )
{

   const char FileName[] = "STOP_GAMER_STOP";

   int Terminate_local = false;


// enforce NFS to flush the file information
   system( "ls > /dev/null" );


// check the stop file
   FILE *File = fopen( FileName, "r" );

   if ( File != NULL )
   {
      Terminate_local = true;
      fclose( File );
   }


// the program will be terminated as long as ONE process has detected the stop file
   MPI_Allreduce( &Terminate_local, &Terminate_global, 1, MPI_INT, MPI_BOR, MPI_COMM_WORLD );


   if ( MPI_Rank == 0  &&  Terminate_global )
   {
      Aux_Message( stdout, "\nDetecting file \"STOP_GAMER_STOP\" --> terminating the program ...\n\n" );

//    remove the stop file
      system( "rm -f STOP_GAMER_STOP" );
   }

} // FUNCTION : End_StopManually
