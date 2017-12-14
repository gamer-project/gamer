#include "GAMER.h"

#ifdef LOAD_BALANCE



void Flag_Grandson( const int lv, const int PID, const int LocalID );




//-------------------------------------------------------------------------------------------------------
// Function    :  LB_GrandsonCheck 
// Description :  Grandson check for the flag operation 
//
// Note        :  Invoked by the function "Flag_Real" 
//
// Parameter   :  FaLv : Target refinement level to be flagged
//-------------------------------------------------------------------------------------------------------
void LB_GrandsonCheck( const int FaLv )
{

// check
   if ( FaLv > NLEVEL-3 )
   {
      Aux_Message( stderr, "WARNING : function <%s> should not be applied to level %d (not > %d) !!\n",
                 __FUNCTION__, FaLv, NLEVEL-3 );
      return;
   }


   const int SonLv   = FaLv + 1;
   const int FaNReal = amr->NPatchComma[ FaLv][1];
   const int MemUnit = 1 + FaNReal/MPI_NRank;         // set arbitrarily

   int   SonPID, TRank, MemSize[MPI_NRank], NQuery[MPI_NRank];
   long  SonLBIdx;
   int  *FaPID_List[MPI_NRank], *FaPID_IdxTable[MPI_NRank];
   long *Query_Temp[MPI_NRank];


// initialize arrays
   for (int r=0; r<MPI_NRank; r++)     
   {
      MemSize   [r] = MemUnit;
      Query_Temp[r] = (long*)malloc( MemSize[r]*sizeof(long) );
      FaPID_List[r] = (int* )malloc( MemSize[r]*sizeof(int ) );
      NQuery    [r] = 0;
   }


// loop over all real patches
   for (int FaPID=0; FaPID<FaNReal; FaPID++)
   {
      SonPID = amr->patch[0][FaLv][FaPID]->son;

      if ( SonPID != -1 )
      {
//       1. son at home --> flag directly
// ==========================================================================================
         if ( SonPID >= 0 )
         {
            for (int LocalID=0; LocalID<8; LocalID++)
            {
               if ( amr->patch[0][SonLv][SonPID+LocalID]->son != -1 )   // if grandson exists
               {
//                flag the corresponding siblings of patch PID
                  Flag_Grandson( FaLv, FaPID, LocalID );
   
//                flag the patch PID
                  amr->patch[0][FaLv][FaPID]->flag = true;
               }
            }
         }


//       2. son not home --> record the query list
// ==========================================================================================
         else
         {
//          determine the target rank and LB_Idx
#           if ( LOAD_BALANCE == HILBERT )
            SonLBIdx = 8*amr->patch[0][FaLv][FaPID]->LB_Idx;   // faster, but may not correspond to LocalID==0
#           else
            SonLBIdx = LB_Corner2Index( SonLv, amr->patch[0][FaLv][FaPID]->corner, CHECK_OFF );
#           endif
            TRank    = SON_OFFSET_LB - SonPID;

//          allocate enough memory
            if ( NQuery[TRank] >= MemSize[TRank] )  
            {
               MemSize   [TRank] += MemUnit;
               Query_Temp[TRank]  = (long*)realloc( Query_Temp[TRank], MemSize[TRank]*sizeof(long) );
               FaPID_List[TRank]  = (int* )realloc( FaPID_List[TRank], MemSize[TRank]*sizeof(int ) );
            }

//          record list
            Query_Temp[TRank][ NQuery[TRank] ] = SonLBIdx;
            FaPID_List[TRank][ NQuery[TRank] ] = FaPID;
            NQuery    [TRank] ++;

         } // if ( SonPID >= 0 ) ... else ...
      } // if ( SonPID != -1 )
   } // for (int FaPID=0; FaPID<FaNReal; FaPID++)


// 3. sort the query list for different ranks
// ==========================================================================================
   for (int r=0; r<MPI_NRank; r++)
   {
      FaPID_IdxTable[r] = new int [ NQuery[r] ];

      Mis_Heapsort( NQuery[r], Query_Temp[r], FaPID_IdxTable[r] );
   }


// 4. transfer data to get reply : (SendBuf_Query --> RecvBuf_Query --> SendBuf_Reply --> RecvBuf_Reply)
// ==========================================================================================   
   int   Query_Disp[MPI_NRank], Reply_Disp[MPI_NRank], NReply[MPI_NRank], NQuery_Total, NReply_Total; 
   int   Counter, SonPID0;
   long *SendBuf_Query=NULL, *RecvBuf_Query=NULL, *QueryPtr=NULL;
   int  *SendBuf_Reply=NULL, *RecvBuf_Reply=NULL, *ReplyPtr=NULL;
   int  *Match=NULL;

// 4.1 send the number of queries
   MPI_Alltoall( NQuery, 1, MPI_INT, NReply, 1, MPI_INT, MPI_COMM_WORLD );

// 4.2 prepare the query and reply arrays
   Query_Disp[0] = 0;
   Reply_Disp[0] = 0;
   for (int r=1; r<MPI_NRank; r++)
   {
      Query_Disp[r] = Query_Disp[r-1] + NQuery[r-1]; 
      Reply_Disp[r] = Reply_Disp[r-1] + NReply[r-1]; 
   }
   NQuery_Total = Query_Disp[MPI_NRank-1] + NQuery[MPI_NRank-1];
   NReply_Total = Reply_Disp[MPI_NRank-1] + NReply[MPI_NRank-1];

   SendBuf_Query = new long [NQuery_Total];
   RecvBuf_Query = new long [NReply_Total];
   SendBuf_Reply = new int  [NReply_Total];
   RecvBuf_Reply = new int  [NQuery_Total];

   Counter = 0;
   for (int r=0; r<MPI_NRank; r++)
   for (int t=0; t<NQuery[r]; t++)
      SendBuf_Query[ Counter ++ ] = Query_Temp[r][t];

   for (int t=0; t<NReply_Total; t++)   SendBuf_Reply[t] = 0;

// 4.3 send queries
   MPI_Alltoallv( SendBuf_Query, NQuery, Query_Disp, MPI_LONG, 
                  RecvBuf_Query, NReply, Reply_Disp, MPI_LONG, MPI_COMM_WORLD );

// 4.4 prepare replies
   for (int r=0; r<MPI_NRank; r++)
   {
      QueryPtr = RecvBuf_Query + Reply_Disp[r];
      ReplyPtr = SendBuf_Reply + Reply_Disp[r];
      Match    = new int [ NReply[r] ];

      Mis_Matching_int( amr->NPatchComma[SonLv][1], amr->LB->IdxList_Real[SonLv], NReply[r], QueryPtr, Match );

      for (int t=0; t<NReply[r]; t++)
      {
#        ifdef GAMER_DEBUG
         if ( Match[t] == -1 )
            Aux_Error( ERROR_INFO, "SonLBIdx = %ld found no matching in Rank %d !!\n", QueryPtr[t], MPI_Rank );
#        endif

         SonPID0 = amr->LB->IdxList_Real_IdxTable[SonLv][ Match[t] ];
         SonPID0 = SonPID0 - SonPID0%8;

         for (int LocalID=0; LocalID<8; LocalID++)
         {
            SonPID = SonPID0 + LocalID;

//          check if grandson exists
            if ( amr->patch[0][SonLv][SonPID]->son != -1 )  ReplyPtr[t] |= ( 1 << LocalID );
         }
      } // for (int t=0; t<NReply[r]; t++)

      delete [] Match;

   } // for (int r=0; r<MPI_NRank; r++)

// 4.5 send replies
   MPI_Alltoallv( SendBuf_Reply, NReply, Reply_Disp, MPI_INT, 
                  RecvBuf_Reply, NQuery, Query_Disp, MPI_INT, MPI_COMM_WORLD );


// 5. flag father patches according to the reply list
// ==========================================================================================
   int FaPID;
   Counter=0;

   for (int r=0; r<MPI_NRank; r++)
   for (int t=0; t<NQuery[r]; t++)
   {
      if ( RecvBuf_Reply[Counter] != 0 )
      {
         FaPID = FaPID_List[r][ FaPID_IdxTable[r][t] ];

//       flag the patch FaPID
         amr->patch[0][FaLv][FaPID]->flag = true;

//       flag the corresponding siblings of patch FaPID
         for (int LocalID=0; LocalID<8; LocalID++)
            if (  RecvBuf_Reply[Counter]  &  ( 1 << LocalID )  )     Flag_Grandson( FaLv, FaPID, LocalID );
      }

      Counter ++;
   }


// free memory
   for (int r=0; r<MPI_NRank; r++)
   {
      free( Query_Temp[r] );
      free( FaPID_List[r] );
      delete [] FaPID_IdxTable[r];
   }

   delete [] SendBuf_Query;
   delete [] RecvBuf_Query;
   delete [] SendBuf_Reply;
   delete [] RecvBuf_Reply;

} // FUNCTION : LB_GrandsonCheck



#endif // #ifdef LOAD_BALANCE
