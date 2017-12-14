#include "GAMER.h"

#ifndef SERIAL




//-------------------------------------------------------------------------------------------------------
// Function    :  MPI_Exit
// Description :  Force the program to be terminated when any error occurs
//-------------------------------------------------------------------------------------------------------
void MPI_Exit()
{

// flush all previous messages
   fflush( stdout ); fflush( stdout ); fflush( stdout );
   fflush( stderr ); fflush( stderr ); fflush( stderr );

   Aux_Message( stderr, "\nProgram terminated with error ...... rank %d\n\n", MPI_Rank );

   MPI_Abort( MPI_COMM_WORLD, 1 );

   exit( 1 );

} // FUNCTION : MPI_Exit



#endif // #ifndef SERIAL
