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

// initialise MPI with OpenMP support for FFTW3 according to (https://www.fftw.org/fftw3_doc/Combining-MPI-and-Threads.html)
#  if (  defined(OVERLAP_MPI) || ( SUPPORT_FFTW == FFTW3 && defined(OPENMP) )  )
   int MPI_Thread_Status = -1;

// MPI_Init_thread( argc, argv, MPI_THREAD_SINGLE,     &MPI_Thread_Status );
// MPI_Init_thread( argc, argv, MPI_THREAD_FUNNELED,   &MPI_Thread_Status );
   MPI_Init_thread( argc, argv, MPI_THREAD_SERIALIZED, &MPI_Thread_Status );
// MPI_Init_thread( argc, argv, MPI_THREAD_MULTIPLE,   &MPI_Thread_Status );
   MPI_Comm_rank( MPI_COMM_WORLD, &MPI_Rank );
   MPI_Comm_size( MPI_COMM_WORLD, &MPI_NRank );

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
   MPI_Comm_size( MPI_COMM_WORLD, &MPI_NRank );

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "Init_MPI ... done\n" );
#  endif

// detect the number of MPI ranks per node and issue a warning if each node has only one rank
// (assume all nodes have the same number of MPI ranks per node)
// reference: https://stackoverflow.com/questions/9022496/how-to-determine-mpi-rank-process-number-local-to-a-socket-node
   int MPI_SizePerNode;
   MPI_Comm shmcomm;
   MPI_Comm_split_type( MPI_COMM_WORLD, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, &shmcomm );
   MPI_Comm_size( shmcomm, &MPI_SizePerNode );

   if ( MPI_SizePerNode == 1  &&  MPI_Rank == 0 )
      Aux_Message( stderr, "WARNING : Each node has only one MPI rank. Using more ranks per node may improve performance !!\n" );

   MPI_Comm_free( &shmcomm );

} // FUNCTION : Init_MPI



#endif // #ifndef SERIAL
