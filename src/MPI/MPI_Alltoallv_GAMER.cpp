#include "GAMER.h"

#ifndef SERIAL



//-------------------------------------------------------------------------------------------------------
// Function    :  MPI_Alltoallv_GAMER
// Description :  wrapper for replacing official MPI_Alltoallv() when the number of elements in Send_NDisp/Recv_NDisp exceed __INT_MAX__
//
// Parameter   :  SendBuf:      data to be send by this rank to other ranks via MPI_Alltoallv
//                Send_NCount:  number of elements to be sent by each rank to other ranks in SendBuff; length euqals MPI_NRank
//                Send_NDisp:   displacement indicating the stride where the sent data (to other ranks) starts in SendBuf for each rank;
//                              length equals to MPI_NRank
//                Send_Datatype: sent data type for MPI
//                RecvBuf:      data to be received by this rank from other ranks via MPI_Alltoallv
//                Recv_NCount:  number of elements to be received by each rank from other ranks in RecvBuff; length euqals MPI_NRank
//                Recv_NDisp:   displacement indicating the stride where the received data (from other ranks) starts in RecvdBuf for each rank;
//                              length equals to MPI_NRank
//                Recv_Datatype: received data type for MPI (MPI_GAMER_REAL/MPI_GAMER_REAL_PAR)
//                comm:          MPI communicator
//-------------------------------------------------------------------------------------------------------
template<typename T> 
void MPI_Alltoallv_GAMER(T *SendBuf, int *Send_NCount, long *Send_NDisp, MPI_Datatype Send_Datatype, T *RecvBuf, int *Recv_NCount, long *Recv_NDisp, MPI_Datatype Recv_Datatype, MPI_Comm comm)
{
   bool use_mpi_gamer_flag = false;
   if ( (Send_NDisp[MPI_NRank-1] > __INT_MAX__) || ( Recv_NDisp[MPI_NRank-1] > __INT_MAX__ ) )  use_mpi_gamer_flag = true;
   MPI_Allreduce(MPI_IN_PLACE, &use_mpi_gamer_flag , 1, MPI_C_BOOL, MPI_LOR, MPI_COMM_WORLD);

   // ********************************purely for testing purpose!! Remember to comment this out for distribution!! ********************************
   Aux_Message( stdout, "  Use MPI_Alltoallv_GAMER wrapper..." );
   use_mpi_gamer_flag = true;
   //

   if ( use_mpi_gamer_flag )
   {
      MPI_Request *req_send_and_recv = new MPI_Request[2*MPI_NRank];

      for(int r=0; r<MPI_NRank; r++)
      {
          MPI_Isend(SendBuf+Send_NDisp[r], (int)Send_NCount[r], Send_Datatype, r, MPI_Rank*MPI_NRank + r       , comm, &req_send_and_recv[2*r  ]);
          MPI_Irecv(RecvBuf+Recv_NDisp[r], (int)Recv_NCount[r], Recv_Datatype, r,        r*MPI_NRank + MPI_Rank, comm, &req_send_and_recv[2*r+1]);
      }
      MPI_Waitall(2*MPI_NRank, req_send_and_recv, MPI_STATUSES_IGNORE);

      delete [] req_send_and_recv;
   }
   else
   {
      int *Send_NDisp_int  = new int [MPI_NRank];
      int *Recv_NDisp_int  = new int [MPI_NRank];

      for (int r=0; r<MPI_NRank; r++)
      {
         Send_NDisp_int[r]  = (int)Send_NDisp[r];
         Recv_NDisp_int[r]  = (int)Recv_NDisp[r];
      }

      MPI_Alltoallv( SendBuf, Send_NCount, Send_NDisp_int, Send_Datatype,
                     RecvBuf, Recv_NCount, Recv_NDisp_int, Recv_Datatype, comm );   

      delete [] Send_NDisp_int;
      delete [] Recv_NDisp_int;
   }
       
} // FUNCTION :MPI_Alltoallv_GAMER()


// explicit template instantiation
template void MPI_Alltoallv_GAMER <float>    ( float  *SendBuf, int *Send_NCount, long *Send_NDisp, MPI_Datatype Send_Datatype, float  *RecvBuf, int *Recv_NCount, long *Recv_NDisp, MPI_Datatype Recv_Datatype, MPI_Comm comm );
template void MPI_Alltoallv_GAMER <double>   ( double *SendBuf, int *Send_NCount, long *Send_NDisp, MPI_Datatype Send_Datatype, double *RecvBuf, int *Recv_NCount, long *Recv_NDisp, MPI_Datatype Recv_Datatype, MPI_Comm comm );

#endif // #ifndef SERIAL


