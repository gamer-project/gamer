#include "GAMER.h"

#ifdef LOAD_BALANCE




//-------------------------------------------------------------------------------------------------------
// Function    :  LB_AllocateFluxArray
// Description :  1. Allocate flux arrays for the real patches at FaLv adjacent to the coarse-fine boundaries
//                2. Construct the MPI send and recv lists for exchanging fluxes at FaLv
//                3. According to the send list, allocate flux arrays for the buffer patches at FaLv
//
// Parameter   :  FaLv : Coarse-grid level
//-------------------------------------------------------------------------------------------------------
void LB_AllocateFluxArray( const int FaLv )
{

// check
   if ( !amr->WithFlux )
      Aux_Message( stderr, "WARNING : why invoking %s when amr->WithFlux is off ??\n", __FUNCTION__ );


   const int SonLv   = FaLv + 1;
   const int FaNReal = amr->NPatchComma[FaLv ][1];
   const int MemUnit = 1 + FaNReal/MPI_NRank;         // set arbitrarily

   int   MemSize[MPI_NRank], TRank, SonPID, SibPID, SibSonPID;
   long  SibSonLBIdx;
   int  *Temp_SibList[MPI_NRank];
   long *LB_RecvF_SibSonLBIdx[MPI_NRank];

   int  *LB_SendF_NList           = amr->LB->SendF_NList          [FaLv];
   int **LB_SendF_IDList          = amr->LB->SendF_IDList         [FaLv];
   int **LB_SendF_SibList         = amr->LB->SendF_SibList        [FaLv];
   int  *LB_RecvF_NList           = amr->LB->RecvF_NList          [FaLv];
   int **LB_RecvF_IDList          = amr->LB->RecvF_IDList         [FaLv];
   int **LB_RecvF_IDList_IdxTable = amr->LB->RecvF_IDList_IdxTable[FaLv];
   int **LB_RecvF_SibList         = amr->LB->RecvF_SibList        [FaLv];


// initialize arrays
   for (int r=0; r<MPI_NRank; r++)
   {
      if ( LB_RecvF_IDList[r] != NULL )   free( LB_RecvF_IDList[r] );

      MemSize             [r] = MemUnit;
      LB_RecvF_SibSonLBIdx[r] = (long*)malloc( MemSize[r]*sizeof(long) );
      Temp_SibList        [r] = (int* )malloc( MemSize[r]*sizeof(int ) );
      LB_RecvF_IDList     [r] = (int* )malloc( MemSize[r]*sizeof(int ) );
      LB_RecvF_NList      [r] = 0;
   }


// 1. deallocate the flux arrays allocated previously
// ============================================================================================================
#  pragma omp parallel for schedule( runtime )
   for (int FaPID=0; FaPID<amr->NPatchComma[FaLv][3]; FaPID++)  amr->patch[0][FaLv][FaPID]->fdelete();


// 2. allocate flux arrays for the real patches and record the unsorted recv list
// ============================================================================================================
   if ( NPatchTotal[SonLv] != 0 )
   {
      for (int FaPID=0; FaPID<FaNReal; FaPID++)
      {
         SonPID = amr->patch[0][FaLv][FaPID]->son;

         if ( SonPID == -1 )
         {
            for (int Sib=0; Sib<6; Sib++)
            {
               SibPID = amr->patch[0][FaLv][FaPID]->sibling[Sib];

               if ( SibPID >= 0 )   // work for the non-periodic B.C. as well
               {
                  SibSonPID = amr->patch[0][FaLv][SibPID]->son;

                  if ( SibSonPID != -1 )
                  {
//                   allocate flux array
                     amr->patch[0][FaLv][FaPID]->fnew( Sib, AUTO_REDUCE_DT );

//                   record the MPI recv list
                     if ( SibSonPID < -1 )   // son is not home
                     {
//                      determine the target rank and LB_Idx
#                       if ( LOAD_BALANCE == HILBERT )
                        SibSonLBIdx = 8*amr->patch[0][FaLv][SibPID]->LB_Idx;  // faster
#                       else
                        SibSonLBIdx = LB_Corner2Index( SonLv, amr->patch[0][FaLv][SibPID]->corner, CHECK_OFF );
#                       endif
                        TRank       = SON_OFFSET_LB - SibSonPID;

//                      allocate enough memory
                        if ( LB_RecvF_NList[TRank] >= MemSize[TRank] )
                        {
                           MemSize             [TRank] += MemUnit;
                           LB_RecvF_IDList     [TRank] = (int* )realloc( LB_RecvF_IDList     [TRank],
                                                                         MemSize[TRank]*sizeof(int ) );
                           Temp_SibList        [TRank] = (int* )realloc( Temp_SibList        [TRank],
                                                                         MemSize[TRank]*sizeof(int ) );
                           LB_RecvF_SibSonLBIdx[TRank] = (long*)realloc( LB_RecvF_SibSonLBIdx[TRank],
                                                                         MemSize[TRank]*sizeof(long) );
                        }

//                      record the recv list
                        LB_RecvF_IDList     [TRank][ LB_RecvF_NList[TRank] ] = FaPID;
                        Temp_SibList        [TRank][ LB_RecvF_NList[TRank] ] = Sib;
                        LB_RecvF_SibSonLBIdx[TRank][ LB_RecvF_NList[TRank] ] = SibSonLBIdx;

                        LB_RecvF_NList      [TRank] ++;

                     } // if ( SibSonPID < -1 )
                  } // if ( SibSonPID != -1 )
               } // if ( SibPID >= 0 )
            } // for (int Sib=0; Sib<6; Sib++)
         } // if ( SonPID == -1 )
      } // for (int FaPID=0; FaPID<FaNReal; FaPID++)
   } // if ( NPatchTotal[SonLv] != 0 )


// 3. sort the recv list
// ============================================================================================================
   for (int r=0; r<MPI_NRank; r++)
   {
      if ( LB_RecvF_IDList_IdxTable[r] != NULL )   delete [] LB_RecvF_IDList_IdxTable[r];
      if ( LB_RecvF_SibList        [r] != NULL )   delete [] LB_RecvF_SibList        [r];

      LB_RecvF_IDList_IdxTable[r] = new int [ LB_RecvF_NList[r] ];
      LB_RecvF_SibList        [r] = new int [ LB_RecvF_NList[r] ];

      Mis_Heapsort( LB_RecvF_NList[r], LB_RecvF_SibSonLBIdx[r], LB_RecvF_IDList_IdxTable[r] );

      for (int t=0; t<LB_RecvF_NList[r]; t++)
         LB_RecvF_SibList[r][t] = Temp_SibList[r][ LB_RecvF_IDList_IdxTable[r][t] ];

//    recv PID and Sib should not repeat
#     ifdef GAMER_DEBUG
      int *TempIDList   = new int [ LB_RecvF_NList[r] ];
      int *TempIdxTable = new int [ LB_RecvF_NList[r] ];
      memcpy( TempIDList, LB_RecvF_IDList[r], LB_RecvF_NList[r]*sizeof(int) );

      Mis_Heapsort( LB_RecvF_NList[r], TempIDList, TempIdxTable );

      for (int t=1; t<LB_RecvF_NList[r]; t++)
         if ( TempIDList[t] == TempIDList[t-1]  &&
              Temp_SibList[r][ TempIdxTable[t] ] == Temp_SibList[r][ TempIdxTable[t-1] ])
            Aux_Error( ERROR_INFO, "FaLv %d, Rank %d, repeated recv PID %d and Sib %d !!\n",
                       FaLv, r, TempIDList[t], Temp_SibList[r][ TempIdxTable[t] ] );

      delete [] TempIDList;
      delete [] TempIdxTable;
#     endif
   } // for (int r=0; r<MPI_NRank; r++)


// 4. broadcast the recv list to all other ranks
// ============================================================================================================
   int   Send_Disp_F[MPI_NRank], Recv_Disp_F[MPI_NRank], NSend_Total_F, NRecv_Total_F, Counter;
   int  *SendBuf_SibID=NULL, *RecvBuf_SibID=NULL;
   long *SendBuf_LBIdx=NULL, *RecvBuf_LBIdx=NULL;

// 4.1 broadcast the number of elements received from different ranks
   MPI_Alltoall( LB_RecvF_NList, 1, MPI_INT, LB_SendF_NList, 1, MPI_INT, MPI_COMM_WORLD );

// 4.2 prepare the MPI buffers
   Send_Disp_F[0] = 0;
   Recv_Disp_F[0] = 0;
   for (int r=1; r<MPI_NRank; r++)
   {
      Send_Disp_F[r] = Send_Disp_F[r-1] + LB_RecvF_NList[r-1];
      Recv_Disp_F[r] = Recv_Disp_F[r-1] + LB_SendF_NList[r-1];
   }
   NSend_Total_F = Send_Disp_F[MPI_NRank-1] + LB_RecvF_NList[MPI_NRank-1];
   NRecv_Total_F = Recv_Disp_F[MPI_NRank-1] + LB_SendF_NList[MPI_NRank-1];

   SendBuf_SibID = new int  [NSend_Total_F];
   RecvBuf_SibID = new int  [NRecv_Total_F];
   SendBuf_LBIdx = new long [NSend_Total_F];
   RecvBuf_LBIdx = new long [NRecv_Total_F];

   Counter = 0;
   for (int r=0; r<MPI_NRank; r++)
   for (int t=0; t<LB_RecvF_NList[r]; t++)
   {
      SendBuf_SibID[ Counter ] = LB_RecvF_SibList    [r][t];
      SendBuf_LBIdx[ Counter ] = LB_RecvF_SibSonLBIdx[r][t];
      Counter ++;
   }

// 4.3 broadcast the recv list
   MPI_Alltoallv( SendBuf_SibID, LB_RecvF_NList, Send_Disp_F, MPI_INT,
                  RecvBuf_SibID, LB_SendF_NList, Recv_Disp_F, MPI_INT,  MPI_COMM_WORLD );

   MPI_Alltoallv( SendBuf_LBIdx, LB_RecvF_NList, Send_Disp_F, MPI_LONG,
                  RecvBuf_LBIdx, LB_SendF_NList, Recv_Disp_F, MPI_LONG, MPI_COMM_WORLD );


// 5. construct the send list and allocate flux arrays for the buffer patches
// ============================================================================================================
   const int MirrorSib[6] = { 1,0,3,2,5,4 };
   int TPID, TSib, MSib, *Match_F=NULL, *RecvPtr=NULL;

   for (int r=0; r<MPI_NRank; r++)
   {
//    5.1 allocate the send and matching array
      if ( LB_SendF_IDList [r] != NULL )  delete [] LB_SendF_IDList [r];
      if ( LB_SendF_SibList[r] != NULL )  delete [] LB_SendF_SibList[r];

      LB_SendF_IDList [r] = new int [ LB_SendF_NList[r] ];
      LB_SendF_SibList[r] = new int [ LB_SendF_NList[r] ];
      Match_F             = new int [ LB_SendF_NList[r] ];

//    5.2 store the send sibling directions
      RecvPtr = RecvBuf_SibID + Recv_Disp_F[r];
      for (int t=0; t<LB_SendF_NList[r]; t++)   LB_SendF_SibList[r][t] = RecvPtr[t];

//    5.3 match the corresponding PID
      Mis_Matching_int( amr->NPatchComma[SonLv][1], amr->LB->IdxList_Real[SonLv], LB_SendF_NList[r],
                        RecvBuf_LBIdx+Recv_Disp_F[r], Match_F );

//    check : all target patches must be found
      #ifdef GAMER_DEBUG
      for (int t=0; t<LB_SendF_NList[r]; t++)
         if ( Match_F[t] == -1 )
            Aux_Error( ERROR_INFO, "FaLv %d, TRank %d, LB_Idx %ld found no matching patches !!\n",
                       FaLv, r, *(RecvBuf_LBIdx+Recv_Disp_F[r]+t) );
#     endif

//    5.4 store the send patch indices and allocate flux arrays
      for (int t=0; t<LB_SendF_NList[r]; t++)
      {
         TSib      = LB_SendF_SibList[r][t];
         MSib      = MirrorSib[ TSib ];
         SibSonPID = amr->LB->IdxList_Real_IdxTable[SonLv][ Match_F[t] ];
         SibPID    = amr->patch[0][SonLv][SibSonPID]->father;

#        ifdef GAMER_DEBUG
         if ( SibPID < 0 )
            Aux_Error( ERROR_INFO, "TRank %d, SonLv %d, SibSonPID %d, SibPID (%d) < 0 !!\n", r, SonLv, SibSonPID, SibPID );
#        endif

         TPID = amr->patch[0][FaLv][SibPID]->sibling[MSib];

#        ifdef GAMER_DEBUG
         if ( TPID < 0 )
            Aux_Error( ERROR_INFO, "TRank %d, FaLv %d, SibPID %d, TPID (%d) < 0 !!\n", r, FaLv, SibPID, TPID );

         if ( TPID < FaNReal )
            Aux_Error( ERROR_INFO, "TRank %d, FaLv %d, SibPID %d, TPID %d is not a buffer patch !!\n",
                       r, FaLv, SibPID, TPID );

         if ( amr->patch[0][FaLv][TPID]->son != -1 )
            Aux_Error( ERROR_INFO, "TRank %d, FaLv %d, TPID %d has son (SonPID = %d) !!\n",
                       r, FaLv, TPID, amr->patch[0][FaLv][TPID]->son );
#        endif

//       record send PID
         LB_SendF_IDList[r][t] = TPID;

//       allocate flux array
         amr->patch[0][FaLv][TPID]->fnew( TSib, AUTO_REDUCE_DT );

      } // for (int t=0; t<LB_SendF_NList[FaLv][r]; t++)


      delete [] Match_F;

   } // for (int r=0; r<MPI_NRank; r++)


// free memory
   for (int r=0; r<MPI_NRank; r++)
   {
      free( Temp_SibList        [r] );
      free( LB_RecvF_SibSonLBIdx[r] );
   }
   delete [] SendBuf_SibID;
   delete [] RecvBuf_SibID;
   delete [] SendBuf_LBIdx;
   delete [] RecvBuf_LBIdx;

} // LB_AllocateFluxArray



#endif // #ifdef LOAD_BALANCE
