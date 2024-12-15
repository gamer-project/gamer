#include "GAMER.h"

#ifndef SERIAL




//-------------------------------------------------------------------------------------------------------
// Function    :  MPI_Alltoallv_GAMER
// Description :  Wrapper for replacing official MPI_Alltoallv() when the numbers of elements in Send_NDisp/Recv_NDisp exceed __INT_MAX__
//
// Parameter   :  SendBuf:       Data to be sent by this rank to other ranks via MPI_Alltoallv
//                Send_NCount:   Number of elements to be sent by each rank to other ranks in SendBuf; length equals MPI_NRank
//                Send_NDisp:    Displacement indicating the stride where the sent data (to other ranks) starts in SendBuf for each rank;
//                               length equals MPI_NRank
//                Send_Datatype: Sent data type for MPI
//                RecvBuf:       Data to be received by this rank from other ranks via MPI_Alltoallv
//                Recv_NCount:   Number of elements to be received by each rank from other ranks in RecvBuf; length equals MPI_NRank
//                Recv_NDisp:    Displacement indicating the stride where the received data (from other ranks) starts in RecvdBuf for each rank;
//                               length equals MPI_NRank
//                Recv_Datatype: Received data type for MPI (MPI_GAMER_REAL/MPI_GAMER_REAL_PAR/MPI_GAMER_LONG_PAR)
//                comm:          MPI communicator
//
// Return      :  RecvBuf
//-------------------------------------------------------------------------------------------------------
template<typename T>
void MPI_Alltoallv_GAMER( T *SendBuf, long *Send_NCount, long *Send_NDisp, MPI_Datatype Send_Datatype,
                          T *RecvBuf, long *Recv_NCount, long *Recv_NDisp, MPI_Datatype Recv_Datatype, MPI_Comm comm )
{

// sanity check for Send_NCount and Recv_NCount
   for (int r=0; r<MPI_NRank; r++)
   {
      if ( Send_NCount[r] > __INT_MAX__ ) Aux_Error( ERROR_INFO, "Send_NCount[%d] (%ld) > __INT_MAX__ (%ld)!!\n", r, Send_NCount[r], (long)__INT_MAX__ );
      if ( Recv_NCount[r] > __INT_MAX__ ) Aux_Error( ERROR_INFO, "Recv_NCount[%d] (%ld) > __INT_MAX__ (%ld)!!\n", r, Recv_NCount[r], (long)__INT_MAX__ );
   }

   bool use_mpi_gamer_flag = false;
   if (  ( Send_NDisp[MPI_NRank-1] > __INT_MAX__ ) || ( Recv_NDisp[MPI_NRank-1] > __INT_MAX__ )  )    use_mpi_gamer_flag = true;
   MPI_Allreduce( MPI_IN_PLACE, &use_mpi_gamer_flag , 1, MPI_C_BOOL, MPI_LOR, MPI_COMM_WORLD );

   if ( use_mpi_gamer_flag )
   {
      MPI_Request *req_send_and_recv = new MPI_Request[2*MPI_NRank];

// for numbering the MPI_Isend/Irecv tags: since there are total MPI_NRank*MPI_NRank data-transferring tasks, we tag the data-transferring task between rank r_1(send) <--> rank r_2(recv) by number: r_1*MPI_NRank+r_2
      for (int r=0; r<MPI_NRank; r++)
      {
          MPI_Isend( SendBuf+Send_NDisp[r], (int)Send_NCount[r], Send_Datatype, r, MPI_Rank*MPI_NRank + r       , comm, &req_send_and_recv[2*r  ] );
          MPI_Irecv( RecvBuf+Recv_NDisp[r], (int)Recv_NCount[r], Recv_Datatype, r,        r*MPI_NRank + MPI_Rank, comm, &req_send_and_recv[2*r+1] );
      }
      MPI_Waitall( 2*MPI_NRank, req_send_and_recv, MPI_STATUSES_IGNORE );

      delete [] req_send_and_recv;
   }
   else
   {
      int *Send_NCount_int = new int [MPI_NRank];
      int *Recv_NCount_int = new int [MPI_NRank];
      int *Send_NDisp_int  = new int [MPI_NRank];
      int *Recv_NDisp_int  = new int [MPI_NRank];

      for (int r=0; r<MPI_NRank; r++)
      {
         Send_NCount_int[r] = (int)Send_NCount[r];
         Recv_NCount_int[r] = (int)Recv_NCount[r];
         Send_NDisp_int [r] = (int)Send_NDisp [r];
         Recv_NDisp_int [r] = (int)Recv_NDisp [r];
      }

      MPI_Alltoallv( SendBuf, Send_NCount_int, Send_NDisp_int, Send_Datatype,
                     RecvBuf, Recv_NCount_int, Recv_NDisp_int, Recv_Datatype, comm );

      delete [] Send_NCount_int;
      delete [] Recv_NCount_int;
      delete [] Send_NDisp_int;
      delete [] Recv_NDisp_int;
   }

} // FUNCTION : MPI_Alltoallv_GAMER



// explicit template instantiation
template void MPI_Alltoallv_GAMER <float>  ( float  *SendBuf, long *Send_NCount, long *Send_NDisp, MPI_Datatype Send_Datatype, float  *RecvBuf, long *Recv_NCount, long *Recv_NDisp, MPI_Datatype Recv_Datatype, MPI_Comm comm );
template void MPI_Alltoallv_GAMER <double> ( double *SendBuf, long *Send_NCount, long *Send_NDisp, MPI_Datatype Send_Datatype, double *RecvBuf, long *Recv_NCount, long *Recv_NDisp, MPI_Datatype Recv_Datatype, MPI_Comm comm );
template void MPI_Alltoallv_GAMER <int>    ( int    *SendBuf, long *Send_NCount, long *Send_NDisp, MPI_Datatype Send_Datatype, int    *RecvBuf, long *Recv_NCount, long *Recv_NDisp, MPI_Datatype Recv_Datatype, MPI_Comm comm );
template void MPI_Alltoallv_GAMER <long>   ( long   *SendBuf, long *Send_NCount, long *Send_NDisp, MPI_Datatype Send_Datatype, long   *RecvBuf, long *Recv_NCount, long *Recv_NDisp, MPI_Datatype Recv_Datatype, MPI_Comm comm );



#endif // #ifndef SERIAL
