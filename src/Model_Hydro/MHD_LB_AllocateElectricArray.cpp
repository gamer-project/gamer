#include "GAMER.h"

#if ( MODEL == HYDRO  &&  defined MHD  &&  defined LOAD_BALANCE )




//-------------------------------------------------------------------------------------------------------
// Function    :  MHD_LB_AllocateElectricArray
// Description :  Prepare for transferring electric field on the coarse-fine boundaries
//
// Note        :  1. Procedure:
//                   (1) Allocate electric field arrays for the real patches at FaLv adjacent to the coarse-fine boundaries
//                   (2) Construct the MPI send and recv lists for exchanging electric field at FaLv
//                   (3) According to the send list, allocate electric field arrays for the buffer patches at FaLv
//                2. Invoked by LB_Init_LoadBalance() and LB_Refine()
//
// Parameter   :  FaLv : Coarse-grid level
//-------------------------------------------------------------------------------------------------------
void MHD_LB_AllocateElectricArray( const int FaLv )
{

// check
   if ( !amr->WithElectric )
      Aux_Message( stderr, "WARNING : why invoking %s when amr->WithElectric is off ??\n", __FUNCTION__ );


   const int SonLv      = FaLv + 1;
   const int FaNReal    = amr->NPatchComma[FaLv][1];
   const int MemUnit    = 1 + FaNReal/MPI_NRank;      // set arbitrarily
   const int MirSib[18] = { 1,0,3,2,5,4,9,8,7,6,13,12,11,10,17,16,15,14 };

   int   MemSize         [MPI_NRank];
   int  *Temp_SibEList   [MPI_NRank];  // sibling direction of the target electric field
   int  *Temp_SibF2CList [MPI_NRank];  // sibling direction from fine to coarse patches
   int  *Recv_SibF2CList [MPI_NRank];
   int  *Send_SibF2CList [MPI_NRank];
   long *Recv_SibSonLBIdx[MPI_NRank];

   int  *LB_SendE_NList           = amr->LB->SendE_NList          [FaLv];
   int **LB_SendE_IDList          = amr->LB->SendE_IDList         [FaLv];
   int **LB_SendE_SibList         = amr->LB->SendE_SibList        [FaLv];
   int  *LB_RecvE_NList           = amr->LB->RecvE_NList          [FaLv];
   int **LB_RecvE_IDList          = amr->LB->RecvE_IDList         [FaLv];
   int **LB_RecvE_IDList_IdxTable = amr->LB->RecvE_IDList_IdxTable[FaLv];
   int **LB_RecvE_SibList         = amr->LB->RecvE_SibList        [FaLv];


// initialize arrays
   for (int r=0; r<MPI_NRank; r++)
   {
      if ( LB_RecvE_IDList[r] != NULL )   free( LB_RecvE_IDList[r] );

      MemSize         [r] = MemUnit;
      Temp_SibEList   [r] = (int* )malloc( MemSize[r]*sizeof(int ) );
      Temp_SibF2CList [r] = (int* )malloc( MemSize[r]*sizeof(int ) );
      Recv_SibSonLBIdx[r] = (long*)malloc( MemSize[r]*sizeof(long) );
      LB_RecvE_IDList [r] = (int* )malloc( MemSize[r]*sizeof(int ) );
      LB_RecvE_NList  [r] = 0;
   }



// 1. deallocate the electric field arrays allocated previously
// ============================================================================================================
#  pragma omp parallel for schedule( runtime )
   for (int FaPID=0; FaPID<amr->NPatchComma[FaLv][3]; FaPID++)  amr->patch[0][FaLv][FaPID]->edelete();



// 2. allocate electric field arrays for the real patches and record the unsorted recv list
// ============================================================================================================
   if ( NPatchTotal[SonLv] != 0 )
   {
      for (int FaPID=0; FaPID<FaNReal; FaPID++)
      {
         const int SonPID = amr->patch[0][FaLv][FaPID]->son;

         if ( SonPID == -1 )
         {
//          2-1 work on faces 0 ~ 5
            for (int FaceID=0; FaceID<6; FaceID++)
            {
               const int SibPID = amr->patch[0][FaLv][FaPID]->sibling[FaceID];

               if ( SibPID >= 0 )   // work for the non-periodic B.C. as well
               {
                  const int SibSonPID = amr->patch[0][FaLv][SibPID]->son;

                  if ( SibSonPID != -1 )
                  {
//                   allocate electrid field array
                     amr->patch[0][FaLv][FaPID]->enew( FaceID, AUTO_REDUCE_DT );

//                   record the MPI recv list
                     if ( SibSonPID < -1 )   // son is not home
                     {
//                      determine the target rank and LB_Idx
#                       if ( LOAD_BALANCE == HILBERT )
                        const int SibSonLBIdx = 8*amr->patch[0][FaLv][SibPID]->LB_Idx;  // faster
#                       else
                        const int SibSonLBIdx = LB_Corner2Index( SonLv, amr->patch[0][FaLv][SibPID]->corner, CHECK_OFF );
#                       endif
                        const int TRank       = SON_OFFSET_LB - SibSonPID;

//                      allocate enough memory
                        if ( LB_RecvE_NList[TRank] >= MemSize[TRank] )
                        {
                           MemSize         [TRank] += MemUnit;
                           Temp_SibEList   [TRank] = (int* )realloc( Temp_SibEList   [TRank],
                                                                     MemSize[TRank]*sizeof(int ) );
                           Temp_SibF2CList [TRank] = (int* )realloc( Temp_SibF2CList [TRank],
                                                                     MemSize[TRank]*sizeof(int ) );
                           Recv_SibSonLBIdx[TRank] = (long*)realloc( Recv_SibSonLBIdx[TRank],
                                                                     MemSize[TRank]*sizeof(long) );
                           LB_RecvE_IDList [TRank] = (int* )realloc( LB_RecvE_IDList [TRank],
                                                                     MemSize[TRank]*sizeof(int ) );
                        }

//                      record the recv list
                        Temp_SibEList   [TRank][ LB_RecvE_NList[TRank] ] = FaceID;
                        Temp_SibF2CList [TRank][ LB_RecvE_NList[TRank] ] = MirSib[FaceID];
                        Recv_SibSonLBIdx[TRank][ LB_RecvE_NList[TRank] ] = SibSonLBIdx;
                        LB_RecvE_IDList [TRank][ LB_RecvE_NList[TRank] ] = FaPID;

                        LB_RecvE_NList  [TRank] ++;

                     } // if ( SibSonPID < -1 )
                  } // if ( SibSonPID != -1 )
               } // if ( SibPID >= 0 )
            } // for (int FaceID=0; FaceID<6; FaceID++)


//          2-2 work on edges 6 ~ 17
            for (int EdgeID=6; EdgeID<18; EdgeID++)
            {
               int  SibC2F[3], MPISibC2F=-1, MPISibPID=-1, MPISibSonPID=-1;
               bool NeedMPI=false, Allocated=false;

               TABLE_SiblingSharingSameEdge( EdgeID, SibC2F, NULL );

               for (int s=0; s<3; s++)
               {
                  const int SibPID    = amr->patch[0][FaLv][FaPID]->sibling[ SibC2F[s] ];
                  const int SibSonPID = ( SibPID >= 0 ) ? amr->patch[0][FaLv][SibPID]->son : -1;

                  if ( SibSonPID != -1 )
                  {
//                   allocate electric field array
                     if ( Allocated == false )
                     {
                        amr->patch[0][FaLv][FaPID]->enew( EdgeID, AUTO_REDUCE_DT );
                        Allocated = true; // allocate it just once
                     }

//                   check if MPI communication is required
                     if ( SibSonPID >= 0 )
                     {
                        NeedMPI = false;
                        break;   // no MPI if any of the sibling-son patches is real
                     }

                     else
                     {
                        NeedMPI      = true;
                        MPISibC2F    = SibC2F[s];
                        MPISibPID    = SibPID;
                        MPISibSonPID = SibSonPID;
                     }
                  } // if ( SibSonPID != -1 )
               } // for (int s=0; s<3; s++)

//             record the MPI recv list
               if ( NeedMPI )
               {
//                determine the target rank and LB_Idx
#                 if ( LOAD_BALANCE == HILBERT )
                  const int SibSonLBIdx = 8*amr->patch[0][FaLv][ MPISibPID ]->LB_Idx;
#                 else
                  const int SibSonLBIdx = LB_Corner2Index( SonLv, amr->patch[0][FaLv][ MPISibPID ]->corner, CHECK_OFF );
#                 endif
                  const int TRank       = SON_OFFSET_LB - MPISibSonPID;

//                allocate enough memory
                  if ( LB_RecvE_NList[TRank] >= MemSize[TRank] )
                  {
                       MemSize         [TRank] += MemUnit;
                       Temp_SibEList   [TRank] = (int* )realloc( Temp_SibEList   [TRank],
                                                                 MemSize[TRank]*sizeof(int ) );
                       Temp_SibF2CList [TRank] = (int* )realloc( Temp_SibF2CList [TRank],
                                                                 MemSize[TRank]*sizeof(int ) );
                       Recv_SibSonLBIdx[TRank] = (long*)realloc( Recv_SibSonLBIdx[TRank],
                                                                 MemSize[TRank]*sizeof(long) );
                       LB_RecvE_IDList [TRank] = (int* )realloc( LB_RecvE_IDList [TRank],
                                                                 MemSize[TRank]*sizeof(int ) );
                  }

//                record the recv list
                  Temp_SibEList   [TRank][ LB_RecvE_NList[TRank] ] = EdgeID;
                  Temp_SibF2CList [TRank][ LB_RecvE_NList[TRank] ] = MirSib[ MPISibC2F ];
                  Recv_SibSonLBIdx[TRank][ LB_RecvE_NList[TRank] ] = SibSonLBIdx;
                  LB_RecvE_IDList [TRank][ LB_RecvE_NList[TRank] ] = FaPID;

                  LB_RecvE_NList  [TRank] ++;

               } // if ( NeedMPI )
            } // for (int EdgeID=6; EdgeID<18; EdgeID++)
         } // if ( SonPID == -1 )
      } // for (int FaPID=0; FaPID<FaNReal; FaPID++)
   } // if ( NPatchTotal[SonLv] != 0 )



// 3. sort the recv list
// ============================================================================================================
   for (int r=0; r<MPI_NRank; r++)
   {
      if ( LB_RecvE_IDList_IdxTable[r] != NULL )   delete [] LB_RecvE_IDList_IdxTable[r];
      if ( LB_RecvE_SibList        [r] != NULL )   delete [] LB_RecvE_SibList        [r];

      Recv_SibF2CList         [r] = new int [ LB_RecvE_NList[r] ];
      LB_RecvE_IDList_IdxTable[r] = new int [ LB_RecvE_NList[r] ];
      LB_RecvE_SibList        [r] = new int [ LB_RecvE_NList[r] ];

      Mis_Heapsort( LB_RecvE_NList[r], Recv_SibSonLBIdx[r], LB_RecvE_IDList_IdxTable[r] );

      for (int t=0; t<LB_RecvE_NList[r]; t++)
      {
         Recv_SibF2CList [r][t] = Temp_SibF2CList[r][ LB_RecvE_IDList_IdxTable[r][t] ];
         LB_RecvE_SibList[r][t] = Temp_SibEList  [r][ LB_RecvE_IDList_IdxTable[r][t] ];
      }

//    recv PID and Sib should not repeat
#     ifdef GAMER_DEBUG
      int *TempIDList   = new int [ LB_RecvE_NList[r] ];
      int *TempIdxTable = new int [ LB_RecvE_NList[r] ];
      memcpy( TempIDList, LB_RecvE_IDList[r], LB_RecvE_NList[r]*sizeof(int) );

      Mis_Heapsort( LB_RecvE_NList[r], TempIDList, TempIdxTable );

      for (int t=1; t<LB_RecvE_NList[r]; t++)
         if ( TempIDList[t] == TempIDList[t-1]  &&
              Temp_SibEList[r][ TempIdxTable[t] ] == Temp_SibEList[r][ TempIdxTable[t-1] ])
            Aux_Error( ERROR_INFO, "FaLv %d, Rank %d, repeated recv PID %d and Sib %d !!\n",
                       FaLv, r, TempIDList[t], Temp_SibEList[r][ TempIdxTable[t] ] );

      delete [] TempIDList;
      delete [] TempIdxTable;
#     endif
   } // for (int r=0; r<MPI_NRank; r++)



// 4. broadcast the recv list to all other ranks
// ============================================================================================================
   int   Send_Disp_E[MPI_NRank], Recv_Disp_E[MPI_NRank], NSend_Total_E, NRecv_Total_E, Counter;
   int  *SendBuf_SibE=NULL, *RecvBuf_SibE=NULL;
   int  *SendBuf_SibF2C=NULL, *RecvBuf_SibF2C=NULL;
   long *SendBuf_LBIdx=NULL, *RecvBuf_LBIdx=NULL;

// 4.1 broadcast the number of elements received from different ranks
   MPI_Alltoall( LB_RecvE_NList, 1, MPI_INT, LB_SendE_NList, 1, MPI_INT, MPI_COMM_WORLD );


// 4.2 prepare the MPI buffers
   Send_Disp_E[0] = 0;
   Recv_Disp_E[0] = 0;
   for (int r=1; r<MPI_NRank; r++)
   {
      Send_Disp_E[r] = Send_Disp_E[r-1] + LB_RecvE_NList[r-1];
      Recv_Disp_E[r] = Recv_Disp_E[r-1] + LB_SendE_NList[r-1];
   }
   NSend_Total_E = Send_Disp_E[MPI_NRank-1] + LB_RecvE_NList[MPI_NRank-1];
   NRecv_Total_E = Recv_Disp_E[MPI_NRank-1] + LB_SendE_NList[MPI_NRank-1];

   SendBuf_SibE   = new int  [NSend_Total_E];
   RecvBuf_SibE   = new int  [NRecv_Total_E];
   SendBuf_SibF2C = new int  [NSend_Total_E];
   RecvBuf_SibF2C = new int  [NRecv_Total_E];
   SendBuf_LBIdx  = new long [NSend_Total_E];
   RecvBuf_LBIdx  = new long [NRecv_Total_E];

   Counter = 0;
   for (int r=0; r<MPI_NRank; r++)
   for (int t=0; t<LB_RecvE_NList[r]; t++)
   {
      SendBuf_SibE  [Counter] = LB_RecvE_SibList[r][t];
      SendBuf_SibF2C[Counter] = Recv_SibF2CList [r][t];
      SendBuf_LBIdx [Counter] = Recv_SibSonLBIdx[r][t];
      Counter ++;
   }


// 4.3 broadcast the recv list
   MPI_Alltoallv( SendBuf_SibE,   LB_RecvE_NList, Send_Disp_E, MPI_INT,
                  RecvBuf_SibE,   LB_SendE_NList, Recv_Disp_E, MPI_INT,  MPI_COMM_WORLD );

   MPI_Alltoallv( SendBuf_SibF2C, LB_RecvE_NList, Send_Disp_E, MPI_INT,
                  RecvBuf_SibF2C, LB_SendE_NList, Recv_Disp_E, MPI_INT,  MPI_COMM_WORLD );

   MPI_Alltoallv( SendBuf_LBIdx,  LB_RecvE_NList, Send_Disp_E, MPI_LONG,
                  RecvBuf_LBIdx,  LB_SendE_NList, Recv_Disp_E, MPI_LONG, MPI_COMM_WORLD );



// 5. construct the send list and allocate electric arrays for the buffer patches
// ============================================================================================================
   int *Match_E=NULL, *RecvPtr=NULL;

   for (int r=0; r<MPI_NRank; r++)
   {
//    5.1 allocate the send and matching arrays
      if ( LB_SendE_IDList [r] != NULL )  delete [] LB_SendE_IDList [r];
      if ( LB_SendE_SibList[r] != NULL )  delete [] LB_SendE_SibList[r];

      Send_SibF2CList [r] = new int [ LB_SendE_NList[r] ];
      LB_SendE_IDList [r] = new int [ LB_SendE_NList[r] ];
      LB_SendE_SibList[r] = new int [ LB_SendE_NList[r] ];
      Match_E             = new int [ LB_SendE_NList[r] ];


//    5.2 store the send sibling directions
      RecvPtr = RecvBuf_SibF2C + Recv_Disp_E[r];
      for (int t=0; t<LB_SendE_NList[r]; t++)   Send_SibF2CList [r][t] = RecvPtr[t];

      RecvPtr = RecvBuf_SibE   + Recv_Disp_E[r];
      for (int t=0; t<LB_SendE_NList[r]; t++)   LB_SendE_SibList[r][t] = RecvPtr[t];


//    5.3 match the corresponding PID
      Mis_Matching_int( amr->NPatchComma[SonLv][1], amr->LB->IdxList_Real[SonLv], LB_SendE_NList[r],
                        RecvBuf_LBIdx+Recv_Disp_E[r], Match_E );

//    check : all target patches must be found
#     ifdef GAMER_DEBUG
      for (int t=0; t<LB_SendE_NList[r]; t++)
         if ( Match_E[t] == -1 )
            Aux_Error( ERROR_INFO, "FaLv %d, TRank %d, LB_Idx %ld found no matching patches !!\n",
                       FaLv, r, *(RecvBuf_LBIdx+Recv_Disp_E[r]+t) );
#     endif


//    5.4 store the send patch indices and allocate the electric field arrays
      for (int t=0; t<LB_SendE_NList[r]; t++)
      {
         const int SibE      = LB_SendE_SibList[r][t];
         const int SibF2C    = Send_SibF2CList [r][t];
         const int SibSonPID = amr->LB->IdxList_Real_IdxTable[SonLv][ Match_E[t] ];
         const int SibPID    = amr->patch[0][SonLv][SibSonPID]->father;

#        ifdef GAMER_DEBUG
         if ( SibPID < 0 )
            Aux_Error( ERROR_INFO, "TRank %d, SonLv %d, SibSonPID %d, SibPID (%d) < 0 !!\n", r, SonLv, SibSonPID, SibPID );
#        endif

         const int TPID = amr->patch[0][FaLv][SibPID]->sibling[SibF2C];

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
         LB_SendE_IDList[r][t] = TPID;

//       allocate electric field
         amr->patch[0][FaLv][TPID]->enew( SibE, AUTO_REDUCE_DT );

      } // for (int t=0; t<LB_SendF_NList[FaLv][r]; t++)

      delete [] Match_E;

   } // for (int r=0; r<MPI_NRank; r++)


// free memory
   for (int r=0; r<MPI_NRank; r++)
   {
      free( Temp_SibEList   [r] );
      free( Temp_SibF2CList [r] );
      free( Recv_SibSonLBIdx[r] );

      delete [] Recv_SibF2CList[r];
      delete [] Send_SibF2CList[r];
   }

   delete [] SendBuf_SibE;
   delete [] RecvBuf_SibE;
   delete [] SendBuf_SibF2C;
   delete [] RecvBuf_SibF2C;
   delete [] SendBuf_LBIdx;
   delete [] RecvBuf_LBIdx;

} // FUNCTION : MHD_LB_AllocateElectricArray



#endif // #if ( MODEL == HYDRO  &&  defined MHD  &&  defined LOAD_BALANCE )
