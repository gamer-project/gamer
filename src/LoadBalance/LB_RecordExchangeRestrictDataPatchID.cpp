#include "GAMER.h"

#ifdef LOAD_BALANCE




//-------------------------------------------------------------------------------------------------------
// Function    :  LB_RecordExchangeRestrictDataPatchID_AllLv
// Description :  Construct the MPI sending and receiving data lists for exchanging the restricted data
//
// Note        :  1. The restricted data will be stored in father patches. If any of these father patches
//                   are NOT real patches, we need to send these restricted data back to their homes
//                2. This function will NOT deallocate any fluid array
//                   --> We have assumed that all fluid arrays have been deallocated by
//                       LB_RecordExchangeDataPatchID()
//                   --> Otherwise some fluid arrays may be useless if they are allocated by
//                       LB_RecordExchangeRestrictDataPatchID() previously but later on the corresponding
//                       real son patches are deallocated
//                       --> Actually it is the case in the current version
//
// Parameter   :  FaLv : Target refinement level for constructing MPI lists
//-------------------------------------------------------------------------------------------------------
void LB_RecordExchangeRestrictDataPatchID( const int FaLv )
{

// nothing to do for the maximum level
   if ( FaLv == NLEVEL-1 )    return;


   const int SonLv    = FaLv + 1;
   const int FaNReal  = amr->NPatchComma[FaLv ][1];
   const int SonNReal = amr->NPatchComma[SonLv][1];
   const int MemUnit  = 1 + ( amr->NPatchComma[FaLv][3] - FaNReal ) / MPI_NRank;    // set arbitrarily

   int   TRank, FaPID, MemSize[MPI_NRank];
   long  FaLBIdx;
   long *LB_SendR_LBIdx[MPI_NRank];

   int  *LB_SendR_NList           = amr->LB->SendR_NList          [FaLv];
   int **LB_SendR_IDList          = amr->LB->SendR_IDList         [FaLv];
   int **LB_SendR_IDList_IdxTable = amr->LB->SendR_IDList_IdxTable[FaLv];
   int  *LB_RecvR_NList           = amr->LB->RecvR_NList          [FaLv];
   int **LB_RecvR_IDList          = amr->LB->RecvR_IDList         [FaLv];


// 1. get the unsorted send list
// ============================================================================================================
   for (int r=0; r<MPI_NRank; r++)
   {
      if ( LB_SendR_IDList[r] != NULL )   free( LB_SendR_IDList[r] );

      MemSize        [r] = MemUnit;
      LB_SendR_IDList[r] = (int* )malloc( MemSize[r]*sizeof(int ) );
      LB_SendR_LBIdx [r] = (long*)malloc( MemSize[r]*sizeof(long) );
      LB_SendR_NList [r] = 0;
   }

// loop over all "real" patches at SonLv with LocalID == 0
   for (int SonPID0=0; SonPID0<SonNReal; SonPID0+=8)
   {
      FaPID = amr->patch[0][SonLv][SonPID0]->father;

#     ifdef GAMER_DEBUG
      if ( FaPID < 0 )
         Aux_Error( ERROR_INFO, "SonLv %d, SonPID0 %d has no father patch (FaPID = %d) !!\n",
                    SonLv, SonPID0, FaPID );
#     endif

//    if father patch is not a real patch --> the restricted data need to be sent back home
      if ( FaPID >= FaNReal )
      {
         FaLBIdx = amr->patch[0][FaLv][FaPID]->LB_Idx;
         TRank   = LB_Index2Rank( FaLv, FaLBIdx, CHECK_ON );

#        ifdef GAMER_DEBUG
         if ( TRank == MPI_Rank )
            Aux_Error( ERROR_INFO, "FaLv %d, FaPID %d, FaLBIdx %ld, TRank == MPI_Rank !!\n",
                       FaLv, FaPID, FaLBIdx );
#        endif

//       allocate enough memory
         if ( LB_SendR_NList[TRank] >= MemSize[TRank] )
         {
            MemSize        [TRank] += MemUnit;
            LB_SendR_IDList[TRank] = (int* )realloc( LB_SendR_IDList[TRank], MemSize[TRank]*sizeof(int ) );
            LB_SendR_LBIdx [TRank] = (long*)realloc( LB_SendR_LBIdx [TRank], MemSize[TRank]*sizeof(long) );
         }

//       record the send list
         LB_SendR_IDList[TRank][ LB_SendR_NList[TRank] ] = FaPID;
         LB_SendR_LBIdx [TRank][ LB_SendR_NList[TRank] ] = FaLBIdx;

//       allocate fluid and pot arrays
         for (int Sg=0; Sg<2; Sg++)    amr->patch[Sg][FaLv][FaPID]->hnew();

#        ifdef GRAVITY // so that the XXX_R lists can also be applied to the potential data
         for (int Sg=0; Sg<2; Sg++)    amr->patch[Sg][FaLv][FaPID]->gnew();
#        endif

         LB_SendR_NList[TRank] ++;

      } // if ( FaPID > FaNReal )
   } // for (int SonPID0=0; SonPID0<SonNReal; SonPID0+=8)


// 2. sort the send list
// ============================================================================================================
   for (int r=0; r<MPI_NRank; r++)
   {
      if ( LB_SendR_IDList_IdxTable[r] != NULL )   delete [] LB_SendR_IDList_IdxTable[r];

      LB_SendR_IDList_IdxTable[r] = new int [ LB_SendR_NList[r] ];

      Mis_Heapsort( LB_SendR_NList[r], LB_SendR_LBIdx[r], LB_SendR_IDList_IdxTable[r] );

//    send PID should not repeat
#     ifdef GAMER_DEBUG
      int *TempIDList = new int [ LB_SendR_NList[r] ];
      memcpy( TempIDList, LB_SendR_IDList[r], LB_SendR_NList[r]*sizeof(int) );

      Mis_Heapsort( LB_SendR_NList[r], TempIDList, NULL );

      for (int t=1; t<LB_SendR_NList[r]; t++)
         if ( TempIDList[t] == TempIDList[t-1] )
            Aux_Error( ERROR_INFO, "FaLv %d, Rank %d, repeated send PID %d !!\n",
                       FaLv, r, TempIDList[t] );

      delete [] TempIDList;
#     endif
   } // for (int r=0; r<MPI_NRank; r++)


// 3. broadcast the send list to all other ranks
// ============================================================================================================
   int   Send_Disp_R[MPI_NRank], Recv_Disp_R[MPI_NRank], NSend_Total_R, NRecv_Total_R, Counter;
   long *SendBuf_R=NULL, *RecvBuf_R=NULL;

// 3.1 broadcast the number of elements sent to different ranks
   MPI_Alltoall( LB_SendR_NList, 1, MPI_INT, LB_RecvR_NList, 1, MPI_INT, MPI_COMM_WORLD );

// 3.2 prepare the MPI buffers
   Send_Disp_R[0] = 0;
   Recv_Disp_R[0] = 0;
   for (int r=1; r<MPI_NRank; r++)
   {
      Send_Disp_R[r] = Send_Disp_R[r-1] + LB_SendR_NList[r-1];
      Recv_Disp_R[r] = Recv_Disp_R[r-1] + LB_RecvR_NList[r-1];
   }
   NSend_Total_R = Send_Disp_R[MPI_NRank-1] + LB_SendR_NList[MPI_NRank-1];
   NRecv_Total_R = Recv_Disp_R[MPI_NRank-1] + LB_RecvR_NList[MPI_NRank-1];

   SendBuf_R = new long [NSend_Total_R];
   RecvBuf_R = new long [NRecv_Total_R];

   Counter = 0;
   for (int r=0; r<MPI_NRank; r++)
   for (int t=0; t<LB_SendR_NList[r]; t++)
      SendBuf_R[ Counter ++ ] = LB_SendR_LBIdx[r][t];

// 3.3 broadcast the send list
   MPI_Alltoallv( SendBuf_R, LB_SendR_NList, Send_Disp_R, MPI_LONG,
                  RecvBuf_R, LB_RecvR_NList, Recv_Disp_R, MPI_LONG, MPI_COMM_WORLD );


// 4. construct the recv list
// ============================================================================================================
   int *Match_R=NULL;

   for (int r=0; r<MPI_NRank; r++)
   {
//    4.1 allocate the recv and matching array
      if ( LB_RecvR_IDList[r] != NULL )   delete [] LB_RecvR_IDList[r];

      LB_RecvR_IDList[r] = new int [ LB_RecvR_NList[r] ];
      Match_R            = new int [ LB_RecvR_NList[r] ];

//    4.2 match the corresponding PID
      Mis_Matching_int( amr->NPatchComma[FaLv][1], amr->LB->IdxList_Real[FaLv], LB_RecvR_NList[r],
                        RecvBuf_R+Recv_Disp_R[r], Match_R );

//    4.3 check: all target real patches must be found
#     ifdef GAMER_DEBUG
      for (int t=0; t<LB_RecvR_NList[r]; t++)
         if ( Match_R[t] == -1 )
            Aux_Error( ERROR_INFO, "FaLv %d, TRank %d, LB_Idx %ld found no matching patches !!\n",
                       FaLv, r, *(RecvBuf_R+Recv_Disp_R[r]+t) );
#     endif

//    4.4 store the recv patch indices
      for (int t=0; t<LB_RecvR_NList[r]; t++)
         LB_RecvR_IDList[r][t] = amr->LB->IdxList_Real_IdxTable[FaLv][ Match_R[t] ];

#     ifdef GAMER_DEBUG
//    4.5 check: every patch in the recv list should have son not home
      int SonPID;

      for (int t=0; t<LB_RecvR_NList[r]; t++)
      {
         FaPID  = LB_RecvR_IDList[r][t];
         SonPID = amr->patch[0][FaLv][FaPID]->son;

         if ( SonPID >= -1 )
            Aux_Error( ERROR_INFO, "TRank %d, FaLv %d, FaPID %d, SonNReal %d, incorrect SonPID = %d !!\n",
                       r, FaLv, FaPID, SonNReal, SonPID );
      }

//    4.6 check: verify the recv list
      int  *TempList_PID            = new int  [FaNReal];
      long *TempList_LBIdx          = new long [FaNReal];
      int  *TempList_LBIdx_IdxTable = new int  [FaNReal];

      Counter = 0;

      for (FaPID=0; FaPID<FaNReal; FaPID++)
      {
         SonPID = amr->patch[0][FaLv][FaPID]->son;

         if ( SonPID < -1 )
         {
            TRank = SON_OFFSET_LB - SonPID;

            if ( TRank == r )
            {
               TempList_PID  [Counter] = FaPID;
               TempList_LBIdx[Counter] = amr->patch[0][FaLv][FaPID]->LB_Idx;
               Counter ++;
            }
         }
      }

      Mis_Heapsort( Counter, TempList_LBIdx, TempList_LBIdx_IdxTable );

      if ( Counter != LB_RecvR_NList[r] )
         Aux_Error( ERROR_INFO, "TRank %d, FaLv %d, Counter (%d) != NList (%d) !!\n",
                    r, FaLv, Counter, LB_RecvR_NList[r] );

      for (int t=0; t<Counter; t++)
         if ( TempList_PID[ TempList_LBIdx_IdxTable[t] ] != LB_RecvR_IDList[r][t] )
            Aux_Error( ERROR_INFO, "TRank %d, FaLv %d, t %d, CheckPID (%d) != RecvPID (%d) !!\n",
                       r, FaLv, t, TempList_PID[ TempList_LBIdx_IdxTable[t] ], LB_RecvR_IDList[r][t] );

      delete [] TempList_PID;
      delete [] TempList_LBIdx;
      delete [] TempList_LBIdx_IdxTable;
#     endif


      delete [] Match_R;

   } // for (int r=0; r<MPI_NRank; r++)


// free memory
   for (int r=0; r<MPI_NRank; r++)  free( LB_SendR_LBIdx[r] );
   delete [] SendBuf_R;
   delete [] RecvBuf_R;

} // FUNCTION : LB_RecordExchangeRestrictDataPatchID



#endif // #ifdef LOAD_BALANCE
