#include "GAMER.h"




//-------------------------------------------------------------------------------------------------------
// Function    :  Output_DumpManually
// Description :  Dump data if the file named "DUMP_GAMER_DUMP" is found 
//
// Parameter   :  Dump_global : Boolean variable determining whether or not to dump data
//-------------------------------------------------------------------------------------------------------
void Output_DumpManually( int &Dump_global )
{

   const char FileName[] = "DUMP_GAMER_DUMP";

   int Dump_local = false;


// enforce NFS to flush the file information
   system( "ls > /dev/null" );


// check the target file
   FILE *File = fopen( FileName, "r" );

   if ( File != NULL )  
   {
      Dump_local = true;
      fclose( File );
   }

// the program will dump data as long as ONE process has detected the target file
   MPI_Allreduce( &Dump_local, &Dump_global, 1, MPI_INT, MPI_BOR, MPI_COMM_WORLD );

   if ( MPI_Rank == 0  &&  Dump_global )  
   {
      Aux_Message( stdout, "\nThe file \"DUMP_GAMER_DUMP\" has been detected --> dump data ...\n\n" );

//    remove the target file
      system( "rm -f DUMP_GAMER_DUMP" );
   }

} // FUNCTION : Output_DumpManually


