#include "ExtractUniform.h"




//-------------------------------------------------------------------------------------------------------
// Function    :  MPI_ExchangeInfo
// Description :  Exchange the real-patch data stored in the SendBuffer and RecvBuffer between neighbor ranks
//
// Parameter   :  TargetRank : The MPI rank to send and receive data
//                SendSize   : The number of data to be sent
//                RecvSize   : The number of data to be receive
//                SendBuffer : Send buffer
//                RecvBuffer : Receive buffer
//-------------------------------------------------------------------------------------------------------
void MPI_ExchangeInfo( const int TargetRank[2], const int SendSize[2], const int RecvSize[2],
                       real *SendBuffer[2], real *RecvBuffer[2] )
{

// set the target rank == MPI_PROC_NULL if there is nothing to be sent or received
   int SendTarget[2], RecvTarget[2];

   for (int t=0; t<2; t++)
   {
      SendTarget[t] = ( SendSize[t] == 0 ) ? MPI_PROC_NULL : TargetRank[t];
      RecvTarget[t] = ( RecvSize[t] == 0 ) ? MPI_PROC_NULL : TargetRank[t];
   }


// exchange data
   MPI_Request Req[4];

#  ifdef FLOAT8
   MPI_Isend( SendBuffer[0], SendSize[0], MPI_DOUBLE, SendTarget[0], 0, MPI_COMM_WORLD, &Req[0] );
   MPI_Isend( SendBuffer[1], SendSize[1], MPI_DOUBLE, SendTarget[1], 1, MPI_COMM_WORLD, &Req[1] );

   MPI_Irecv( RecvBuffer[0], RecvSize[0], MPI_DOUBLE, RecvTarget[0], 1, MPI_COMM_WORLD, &Req[2] );
   MPI_Irecv( RecvBuffer[1], RecvSize[1], MPI_DOUBLE, RecvTarget[1], 0, MPI_COMM_WORLD, &Req[3] );

#  else

   MPI_Isend( SendBuffer[0], SendSize[0], MPI_FLOAT,  SendTarget[0], 0, MPI_COMM_WORLD, &Req[0] );
   MPI_Isend( SendBuffer[1], SendSize[1], MPI_FLOAT,  SendTarget[1], 1, MPI_COMM_WORLD, &Req[1] );

   MPI_Irecv( RecvBuffer[0], RecvSize[0], MPI_FLOAT,  RecvTarget[0], 1, MPI_COMM_WORLD, &Req[2] );
   MPI_Irecv( RecvBuffer[1], RecvSize[1], MPI_FLOAT,  RecvTarget[1], 0, MPI_COMM_WORLD, &Req[3] );
#  endif

   MPI_Waitall( 4, Req, MPI_STATUSES_IGNORE );

}
