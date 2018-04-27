#include "GAMER.h"

#ifndef SERIAL




//-------------------------------------------------------------------------------------------------------
// Function    :  MPI_ExchangeBufferPosition
// Description :  Exchange the positions of buffer patches
//
// Note        :  Currently it only works for Buf_AllocateBufferPatch()
//
// Parameter   :  NSend        : Number of data to be sent in 26 sibling directions
//                NRecv        : Number of data to be received in 26 sibling directions
//                Send_PosList : Send buffer
//                Recv_PosList : Receive buffer
//-------------------------------------------------------------------------------------------------------
void MPI_ExchangeBufferPosition( int NSend[26], int NRecv[26], int *Send_PosList[26], int *Recv_PosList[26] )
{

   const int v[26] = { 0,1,2,3,4,5,6,9,7,8,10,13,11,12,14,17,16,15,18,25,19,24,20,23,21,22 };
   int SendTarget[2], RecvTarget[2];
   MPI_Request Req[4];


// a. get the NRecv
   for (int s=0; s<26; s+=2)
   {
//    properly deal with the non-periodic B.C.
      for (int t=0; t<2; t++)
      {
         SendTarget[t] = ( MPI_SibRank[ v[s+t] ] < 0 ) ? MPI_PROC_NULL : MPI_SibRank[ v[s+t] ];
         RecvTarget[t] = ( MPI_SibRank[ v[s+t] ] < 0 ) ? MPI_PROC_NULL : MPI_SibRank[ v[s+t] ];
      }

      MPI_Isend( &NSend[ v[s  ] ], 1, MPI_INT, SendTarget[0], 0, MPI_COMM_WORLD, &Req[0] );
      MPI_Isend( &NSend[ v[s+1] ], 1, MPI_INT, SendTarget[1], 1, MPI_COMM_WORLD, &Req[1] );

      MPI_Irecv( &NRecv[ v[s  ] ], 1, MPI_INT, RecvTarget[0], 1, MPI_COMM_WORLD, &Req[2] );
      MPI_Irecv( &NRecv[ v[s+1] ], 1, MPI_INT, RecvTarget[1], 0, MPI_COMM_WORLD, &Req[3] );

      MPI_Waitall( 4, Req, MPI_STATUSES_IGNORE );

//    properly deal with the non-periodic B.C.
      for (int t=0; t<2; t++)
         if ( MPI_SibRank[ v[s+t] ] < 0 )    NRecv[ v[s+t] ] = 0;

   } // for (int s=0; s<26; s+=2)


// b. allocate memory for Recv_PosList (which will be deallocated in the function "Buf_AllocateBufferPatch")
   for (int s=0; s<26; s++)   Recv_PosList[s] = new int [ NRecv[s] ];


// c. get the Recv_PosList
   for (int s=0; s<26; s+=2)
   {
      for (int t=0; t<2; t++)
      {
         SendTarget[t] = ( NSend[ v[s+t] ] == 0 ) ? MPI_PROC_NULL : MPI_SibRank[ v[s+t] ];
         RecvTarget[t] = ( NRecv[ v[s+t] ] == 0 ) ? MPI_PROC_NULL : MPI_SibRank[ v[s+t] ];

#        ifdef GAMER_DEBUG
         if (  NSend[ v[s+t] ] != 0  &&  ( SendTarget[t] < 0 || SendTarget[t] >= MPI_NRank )  )
            Aux_Error( ERROR_INFO, "incorrect SendTarget[%d] = %d !!\n", t, SendTarget[t] );
         if (  NRecv[ v[s+t] ] != 0  &&  ( RecvTarget[t] < 0 || RecvTarget[t] >= MPI_NRank )  )
            Aux_Error( ERROR_INFO, "incorrect RecvTarget[%d] = %d !!\n", t, RecvTarget[t] );
#        endif
      }

      MPI_Isend( Send_PosList[ v[s  ] ], NSend[ v[s  ] ], MPI_INT, SendTarget[0], 2, MPI_COMM_WORLD, &Req[0] );
      MPI_Isend( Send_PosList[ v[s+1] ], NSend[ v[s+1] ], MPI_INT, SendTarget[1], 3, MPI_COMM_WORLD, &Req[1] );

      MPI_Irecv( Recv_PosList[ v[s  ] ], NRecv[ v[s  ] ], MPI_INT, RecvTarget[0], 3, MPI_COMM_WORLD, &Req[2] );
      MPI_Irecv( Recv_PosList[ v[s+1] ], NRecv[ v[s+1] ], MPI_INT, RecvTarget[1], 2, MPI_COMM_WORLD, &Req[3] );

      MPI_Waitall( 4, Req, MPI_STATUSES_IGNORE );

   } // for (int s=0; s<26; s+=2)

} // FUNCTION : MPI_ExchangeBufferPosition



#endif // #ifndef SERIAL
