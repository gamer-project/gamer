#include "GAMER.h"

#ifdef GRAVITY




//-------------------------------------------------------------------------------------------------------
// Function    :  MPI_CubeSlice
// Description :  Transfer data for the functions "Cube_to_Slice" and "Slice_to_Cube" 
//
// Parameter   :  Dir     : true  --> Cube_to_Slice
//                          false --> Slice_to_Cube
//                SendBuf : Send buffer
//                RecvBuf : Receive buffer
//-------------------------------------------------------------------------------------------------------
void MPI_CubeSlice( const bool Dir, real *SendBuf, real *RecvBuf )
{

   const int MPI_NRank_XY = MPI_NRank_X[0] * MPI_NRank_X[1];
   const int Size         = NX0[0] * NX0[1] * NX0_TOT[2] / MPI_NRank;
   const int Tag          = MPI_Rank % MPI_NRank_XY;
   const int Target1      = MPI_Rank_X[2] * MPI_NRank_XY;
   const int Target2      = ( MPI_Rank / MPI_NRank_XY ) * MPI_NRank_XY;

   int SendTag, RecvTag, SendTarget0, RecvTarget0;
   MPI_Request Req[MPI_NRank_XY][2];

   SendTag = RecvTag = NULL_INT;

   if ( Dir )
   {
      RecvTag     = Tag;
      SendTarget0 = Target1;
      RecvTarget0 = Target2;
   }

   else
   {
      SendTag     = Tag;
      SendTarget0 = Target2;
      RecvTarget0 = Target1;
   }


   for (int lay=0; lay<MPI_NRank_XY; lay++)
   {
      if ( Dir )  SendTag = lay;
      else        RecvTag = lay;

#     ifdef FLOAT8
      MPI_Isend( SendBuf+lay*Size, Size, MPI_DOUBLE, SendTarget0+lay, SendTag, MPI_COMM_WORLD, &Req[lay][0] );
      MPI_Irecv( RecvBuf+lay*Size, Size, MPI_DOUBLE, RecvTarget0+lay, RecvTag, MPI_COMM_WORLD, &Req[lay][1] );
#     else
      MPI_Isend( SendBuf+lay*Size, Size, MPI_FLOAT,  SendTarget0+lay, SendTag, MPI_COMM_WORLD, &Req[lay][0] );
      MPI_Irecv( RecvBuf+lay*Size, Size, MPI_FLOAT,  RecvTarget0+lay, RecvTag, MPI_COMM_WORLD, &Req[lay][1] );
#     endif
   }

   MPI_Waitall( 2*MPI_NRank_XY, &Req[0][0], MPI_STATUSES_IGNORE );

} // FUNCTION : MPI_CubeSlice



#endif // #ifdef GRAVITY
