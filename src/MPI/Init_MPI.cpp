#include "GAMER.h"

#ifndef SERIAL




//-------------------------------------------------------------------------------------------------------
// Function    :  Init_MPI
// Description :  Initialize MPI
//
// Parameter   :  argc, argv --> please inherit from the main function
//-------------------------------------------------------------------------------------------------------
void Init_MPI( int *argc, char ***argv )
{

#  ifdef OVERLAP_MPI
   int MPI_Thread_Status = -1;

// MPI_Init_thread( argc, argv, MPI_THREAD_SINGLE,     &MPI_Thread_Status );
// MPI_Init_thread( argc, argv, MPI_THREAD_FUNNELED,   &MPI_Thread_Status );
   MPI_Init_thread( argc, argv, MPI_THREAD_SERIALIZED, &MPI_Thread_Status );
// MPI_Init_thread( argc, argv, MPI_THREAD_MULTIPLE,   &MPI_Thread_Status );
   MPI_Comm_rank( MPI_COMM_WORLD, &MPI_Rank );

   if ( MPI_Rank == 0 )
   {
      Aux_Message( stdout, "Init_MPI ... done\n" );
      Aux_Message( stdout, "   Current level of MPI thread support : " );

      switch ( MPI_Thread_Status )
      {
         case MPI_THREAD_SINGLE:       Aux_Message( stdout, "MPI_THREAD_SINGLE\n" );       break;
         case MPI_THREAD_FUNNELED:     Aux_Message( stdout, "MPI_THREAD_FUNNELED\n" );     break;
         case MPI_THREAD_SERIALIZED:   Aux_Message( stdout, "MPI_THREAD_SERIALIZED\n" );   break;
         case MPI_THREAD_MULTIPLE:     Aux_Message( stdout, "MPI_THREAD_MULTIPLE\n" );     break;

         default:                      Aux_Message( stdout, "UNKNOWN\n" );
      }
   }

#  else

   MPI_Init( argc, argv );
   MPI_Comm_rank( MPI_COMM_WORLD, &MPI_Rank );

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "Init_MPI ... done\n" );
#  endif

} // FUNCTION : Init_MPI



#endif // #ifndef SERIAL
