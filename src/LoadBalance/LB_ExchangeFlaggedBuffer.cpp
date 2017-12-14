#include "GAMER.h"

#ifdef LOAD_BALANCE




//-------------------------------------------------------------------------------------------------------
// Function    :  LB_ExchangeFlaggedBuffer 
// Description :  Flag real patches flagged by either flag-buffer check or grandson check in other 
//                MPI ranks
//
// Note        :  1. After calling this function, all ranks should have correct flag results for 
//                   all real patches
//                2. Flag list is sorted by the recv rank instead of the send rank
//                   --> Sorting and matching only need to be done once
//
// Parameter   :  lv : Target refinement level
//-------------------------------------------------------------------------------------------------------
void LB_ExchangeFlaggedBuffer( const int lv )
{

   const int NSibBuf = amr->NPatchComma[lv][2] - amr->NPatchComma[lv][1];
   const int MemUnit = 1 + NSibBuf/MPI_NRank;      // set arbitrarily

   int   TRank, MemSize[MPI_NRank], NSend[MPI_NRank];
   long  LBIdx;
   long *Send_Temp[MPI_NRank];

// initialize arrays
   for (int r=0; r<MPI_NRank; r++)
   {
      MemSize  [r] = MemUnit;
      Send_Temp[r] = (long*)malloc( MemSize[r]*sizeof(long) );
      NSend    [r] = 0;
   }


// 1. collect the LB_Idx of all flagged sibling-buffer patches
// ==========================================================================================
   for (int PID=amr->NPatchComma[lv][1]; PID<amr->NPatchComma[lv][2]; PID++)
   {
      if ( amr->patch[0][lv][PID]->flag )
      {
         LBIdx = amr->patch[0][lv][PID]->LB_Idx;
         TRank = LB_Index2Rank( lv, LBIdx, CHECK_ON );

//       allocate enough memory
         if ( NSend[TRank] >= MemSize[TRank] )  
         {
            MemSize  [TRank] += MemUnit;
            Send_Temp[TRank]  = (long*)realloc( Send_Temp[TRank], MemSize[TRank]*sizeof(long) );
         }

//       record list
         Send_Temp[TRank][ NSend[TRank] ++ ] = LBIdx;

      } // if ( amr->patch[0][lv][PID]->flag )
   } // for (int PID=amr->NPatchComma[lv][1]; PID<amr->NPatchComma[lv][2]; PID++)


// 2. broadcast the unsorted send list to all other ranks
// ==========================================================================================
   int   NRecv[MPI_NRank], Send_Disp[MPI_NRank], Recv_Disp[MPI_NRank], NSend_Total, NRecv_Total, Counter;
   long *SendBuf=NULL, *RecvBuf=NULL; 

// 2.1 broadcast the number of elements sent to different ranks
   MPI_Alltoall( NSend, 1, MPI_INT, NRecv, 1, MPI_INT, MPI_COMM_WORLD );

// 2.2 prepare the MPI buffers
   Send_Disp[0] = 0;
   Recv_Disp[0] = 0;
   for (int r=1; r<MPI_NRank; r++)
   {
      Send_Disp[r] = Send_Disp[r-1] + NSend[r-1];
      Recv_Disp[r] = Recv_Disp[r-1] + NRecv[r-1];
   }
   NSend_Total = Send_Disp[MPI_NRank-1] + NSend[MPI_NRank-1];
   NRecv_Total = Recv_Disp[MPI_NRank-1] + NRecv[MPI_NRank-1];

   SendBuf = new long [NSend_Total];
   RecvBuf = new long [NRecv_Total];

   Counter = 0;
   for (int r=0; r<MPI_NRank; r++)
   for (int t=0; t<NSend[r]; t++)
      SendBuf[ Counter ++ ] = Send_Temp[r][t];

// 2.3 broadcast the send list
   MPI_Alltoallv( SendBuf, NSend, Send_Disp, MPI_LONG, 
                  RecvBuf, NRecv, Recv_Disp, MPI_LONG, MPI_COMM_WORLD );


// 3. flag real patches according to the received flag results
// ============================================================================================================
   int  TPID; 
   int *Match = new int [NRecv_Total];

// 3.1 sort the received list
   Mis_Heapsort( NRecv_Total, RecvBuf, NULL );

// 3.2 match the corresponding PID
   Mis_Matching_int( amr->NPatchComma[lv][1], amr->LB->IdxList_Real[lv], NRecv_Total, RecvBuf, Match );

// 3.3 flag
   for (int t=0; t<NRecv_Total; t++)
   {
//    all target real patches must be found
#     ifdef GAMER_DEBUG
      if ( Match[t] == -1 )
         Aux_Error( ERROR_INFO, "lv %d, LB_Idx %ld found no matching patches !!\n", lv, RecvBuf[t] );
#     endif

      TPID = amr->LB->IdxList_Real_IdxTable[lv][ Match[t] ];

      amr->patch[0][lv][TPID]->flag = true;
   }


// free memory
   for (int r=0; r<MPI_NRank; r++)  free( Send_Temp[r] );

   delete [] SendBuf;
   delete [] RecvBuf;
   delete [] Match;

} // FUNCTION : LB_ExchangeFlaggedBuffer



#endif // #ifdef LOAD_BALANCE
