#include "ExtractUniform.h"




//-------------------------------------------------------------------------------------------------------
// Function    :  MPI_Exit
// Description :  Force the program to be terminated when any error occurs
//-------------------------------------------------------------------------------------------------------
void MPI_Exit()
{

   cout << flush << flush << flush;

   fprintf( stderr, "\nProgram termination ...... rank %d\n\n", MyRank );

   MPI_Abort( MPI_COMM_WORLD, 0 );

}


