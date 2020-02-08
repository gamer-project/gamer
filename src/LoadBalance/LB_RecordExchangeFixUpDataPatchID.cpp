#include "GAMER.h"

#ifdef LOAD_BALANCE




//-------------------------------------------------------------------------------------------------------
// Function    :  LB_RecordExchangeFixUpDataPatchID
// Description :  Construct the MPI sending and receiving data lists for exchanging hydro/MHD data after
//                the fix-up operations
//
// Note        :  1. Must be invoked AFTER LB_RecordExchangeDataPatchID()
//                2. All real and buffer patches must know whether or not they have sons
//                3. Lists constructed here (SendX/RecvX/SendY/RecvY) are subsets of the lists (SendH/RecvH)
//                   --> To reduce the amount of data to be transferred after the fix-up operations
//                4. Aassuming Flu_ParaBuf < PATCH_SIZE
//                   --> Because we only check if the interfaces between real and sibling-buffer patches
//                       coincide with coarse-fine boundaaries
//                       --> But for Flu_ParaBuf == PATCH_SIZE, we may need to exchange data even if
//                           the interfaces between real and sibling-buffer patches do NOT coincide with
//                           coarse-fine boundaaries
//                       --> For example: real_patch(lv 1, GID 0) | buffer_patch(lv 1, GID 1) | buffer_patch(lv 2, GID 2)
//                           In this case, even though GID 0 and GID 1 are at the same level, we still need to exchange
//                           data for GID 1 because it's right interface is a coarse-fine boundary. But this corner case
//                           is currently not taken care of.
//                       --> One quick solution is to replace DATA_AFTER_FIXUP by DATA_GENERAL when calling Buf_GetBufferData()
//                           after Flu_FixUp() in EvolveLevel()
//
// Parameter   :  Lv : Target refinement level for recording MPI lists
//-------------------------------------------------------------------------------------------------------
void LB_RecordExchangeFixUpDataPatchID( const int Lv )
{

   int  *LB_SendH_NList           = amr->LB->SendH_NList          [Lv];
   int **LB_SendH_IDList          = amr->LB->SendH_IDList         [Lv];
   int **LB_SendH_SibList         = amr->LB->SendH_SibList        [Lv];
   int  *LB_RecvH_NList           = amr->LB->RecvH_NList          [Lv];
   int **LB_RecvH_IDList          = amr->LB->RecvH_IDList         [Lv];
   int **LB_RecvH_IDList_IdxTable = amr->LB->RecvH_IDList_IdxTable[Lv];
   int **LB_RecvH_SibList         = amr->LB->RecvH_SibList        [Lv];

   int  *LB_SendX_NList           = amr->LB->SendX_NList          [Lv];
   int  *LB_SendX_NResList        = amr->LB->SendX_NResList       [Lv];
   int **LB_SendX_IDList          = amr->LB->SendX_IDList         [Lv];
   int **LB_SendX_SibList         = amr->LB->SendX_SibList        [Lv];
   int  *LB_RecvX_NList           = amr->LB->RecvX_NList          [Lv];
   int  *LB_RecvX_NResList        = amr->LB->RecvX_NResList       [Lv];
   int **LB_RecvX_IDList          = amr->LB->RecvX_IDList         [Lv];
   int **LB_RecvX_SibList         = amr->LB->RecvX_SibList        [Lv];

#  ifdef MHD
   int  *LB_SendY_NList           = amr->LB->SendY_NList          [Lv];
   int **LB_SendY_IDList          = amr->LB->SendY_IDList         [Lv];
   int **LB_SendY_SibList         = amr->LB->SendY_SibList        [Lv];
   int  *LB_RecvY_NList           = amr->LB->RecvY_NList          [Lv];
   int **LB_RecvY_IDList          = amr->LB->RecvY_IDList         [Lv];
   int **LB_RecvY_SibList         = amr->LB->RecvY_SibList        [Lv];

   int **LB_SendY_SubID           = new int* [MPI_NRank];
#  endif

   const int SibMask[6] =
   {   ~(  (1<<19) | (1<<16) | (1<<21) | (1<< 7) | (1<< 1) | (1<< 9) | (1<<23) | (1<<17) | (1<<25)  ),
       ~(  (1<<18) | (1<<14) | (1<<20) | (1<< 6) | (1<< 0) | (1<< 8) | (1<<22) | (1<<15) | (1<<24)  ),
       ~(  (1<<20) | (1<<11) | (1<<21) | (1<< 8) | (1<< 3) | (1<< 9) | (1<<24) | (1<<13) | (1<<25)  ),
       ~(  (1<<18) | (1<<10) | (1<<19) | (1<< 6) | (1<< 2) | (1<< 7) | (1<<22) | (1<<12) | (1<<23)  ),
       ~(  (1<<22) | (1<<12) | (1<<23) | (1<<15) | (1<< 5) | (1<<17) | (1<<24) | (1<<13) | (1<<25)  ),
       ~(  (1<<18) | (1<<10) | (1<<19) | (1<<14) | (1<< 4) | (1<<16) | (1<<20) | (1<<11) | (1<<21)  )   };


// 1. initialize arrays
// ============================================================================================================
   for (int r=0; r<MPI_NRank; r++)
   {
      if ( LB_SendX_IDList [r] != NULL )  delete [] LB_SendX_IDList [r];
      if ( LB_SendX_SibList[r] != NULL )  delete [] LB_SendX_SibList[r];
      if ( LB_RecvX_IDList [r] != NULL )  delete [] LB_RecvX_IDList [r];
      if ( LB_RecvX_SibList[r] != NULL )  delete [] LB_RecvX_SibList[r];

#     ifdef MHD
      if ( LB_SendY_IDList [r] != NULL )  delete [] LB_SendY_IDList [r];
      if ( LB_SendY_SibList[r] != NULL )  delete [] LB_SendY_SibList[r];
      if ( LB_RecvY_IDList [r] != NULL )  delete [] LB_RecvY_IDList [r];
      if ( LB_RecvY_SibList[r] != NULL )  delete [] LB_RecvY_SibList[r];
#     endif

      LB_SendX_NList   [r] = 0;
      LB_SendX_NResList[r] = 0;
      LB_SendX_IDList  [r] = new int [ LB_SendH_NList[r] ];
      LB_SendX_SibList [r] = new int [ LB_SendH_NList[r] ];
      LB_RecvX_NList   [r] = 0;
      LB_RecvX_NResList[r] = 0;
      LB_RecvX_IDList  [r] = new int [ LB_RecvH_NList[r] ];
      LB_RecvX_SibList [r] = new int [ LB_RecvH_NList[r] ];

#     ifdef MHD
      LB_SendY_NList   [r] = 0;
      LB_SendY_IDList  [r] = new int [ LB_SendH_NList[r] ];
      LB_SendY_SibList [r] = new int [ LB_SendH_NList[r] ];
      LB_RecvY_NList   [r] = 0;
      LB_RecvY_IDList  [r] = new int [ LB_RecvH_NList[r] ];
      LB_RecvY_SibList [r] = new int [ LB_RecvH_NList[r] ];

      LB_SendY_SubID   [r] = new int [ LB_SendH_NList[r] ];
#     endif
   } // for (int r=0; r<MPI_NRank; r++)



// 2. data to be sent and received after the restriction fix-up operation
// ============================================================================================================
// note that even when OPT__FIXUP_RESTRICT is off we still need to do data restriction in several places
// (e.g., restart and OPT__CORR_AFTER_ALL_SYNC)
// --> for simplicity and sustainability, we always prepare the following data transfer list (even when
//     OPT__FIXUP_RESTRICT is off)
// if ( OPT__FIXUP_RESTRICT )
   if ( true )
   {
//    2.1 send
      for (int r=0; r<MPI_NRank; r++)
      {
         for (int t=0; t<LB_SendH_NList[r]; t++)
         {
            const int SPID = LB_SendH_IDList[r][t];

            if ( amr->patch[0][Lv][SPID]->son != -1 ) // sons may not be home --> need LB_RecordExchangeRestrictDataPatchID()
            {
               LB_SendX_IDList [r][ LB_SendX_NList[r] ] = SPID;
               LB_SendX_SibList[r][ LB_SendX_NList[r] ] = LB_SendH_SibList[r][t];
               LB_SendX_NList  [r] ++;
            }
         }

         LB_SendX_NResList[r] = LB_SendX_NList[r];
      } // for (int r=0; r<MPI_NRank; r++)


//    2.2 recv
      for (int r=0; r<MPI_NRank; r++)
      {
         for (int t=0; t<LB_RecvH_NList[r]; t++)
         {
            const int RPID = LB_RecvH_IDList[r][ LB_RecvH_IDList_IdxTable[r][t] ];

            if ( amr->patch[0][Lv][RPID]->son != -1 ) // assuming that all buffer patches know whether or not they have sons
                                                      // --> LB_FindSonNotHome() applies to both real and buffer patches
            {
               LB_RecvX_IDList [r][ LB_RecvX_NList[r] ] = RPID;
               LB_RecvX_SibList[r][ LB_RecvX_NList[r] ] = LB_RecvH_SibList[r][t];
               LB_RecvX_NList  [r] ++;
            }
         }

         LB_RecvX_NResList[r] = LB_RecvX_NList[r];
      } // for (int r=0; r<MPI_NRank; r++)
   } // if ( true )



// 3. data to be sent and received after the flux fix-up operation
// ============================================================================================================
   if ( OPT__FIXUP_FLUX )
   {
//    check
      if ( !amr->WithFlux )
         Aux_Error( ERROR_INFO, "amr->WithFlux is off --> no flux array is required !!\n" );


//    3.1 send
      for (int r=0; r<MPI_NRank; r++)
      for (int t=0; t<LB_SendH_NList[r]; t++)
      {
         const int SPID = LB_SendH_IDList[r][t];
         bool GotYou = false;

         for (int s=0; s<6; s++)
         {
            if ( amr->patch[0][Lv][SPID]->flux[s] != NULL )
            {
#              ifdef GAMER_DEBUG
               if ( amr->patch[0][Lv][SPID]->son != -1 )
                  Aux_Error( ERROR_INFO, "Lv %d, SPID %d, SonPID = %d != -1 <-> flux[%d] != NULL !!\n",
                             Lv, SPID, amr->patch[0][Lv][SPID]->son, s );
#              endif

               if ( LB_SendH_SibList[r][t] & SibMask[s] )
               {
                  if ( GotYou )
                     LB_SendX_SibList[r][ LB_SendX_NList[r] ] |= 1 << s;

                  else
                  {
                     LB_SendX_SibList[r][ LB_SendX_NList[r] ]  = 1 << s;
                     LB_SendX_IDList [r][ LB_SendX_NList[r] ]  = SPID;
                     GotYou = true;
                  }
               }
            } // if ( amr->patch[0][Lv][SPID]->flux[s] != NULL )
         } // for (int s=0; s<6; s++)

         if ( GotYou )  LB_SendX_NList[r] ++;
      } // for (int r=0; r<MPI_NRank; r++) for (int t=0; t<LB_SendH_NList[r]; t++)


//    3.2 recv
      for (int r=0; r<MPI_NRank; r++)
      for (int t=0; t<LB_RecvH_NList[r]; t++)
      {
         const int RPID = LB_RecvH_IDList[r][ LB_RecvH_IDList_IdxTable[r][t] ];
         bool GotYou = false;

         if ( amr->patch[0][Lv][RPID]->son == -1 )
         for (int s=0; s<6; s++)
         {
            const int SibPID = amr->patch[0][Lv][RPID]->sibling[s];

            if ( SibPID >= 0  &&  amr->patch[0][Lv][SibPID]->son != -1 )   // work for non-periodic B.C. as well
            {
               if ( LB_RecvH_SibList[r][t] & SibMask[s] )
               {
                  if ( GotYou )
                     LB_RecvX_SibList[r][ LB_RecvX_NList[r] ] |= 1 << s;

                  else
                  {
                     LB_RecvX_SibList[r][ LB_RecvX_NList[r] ]  = 1 << s;
                     LB_RecvX_IDList [r][ LB_RecvX_NList[r] ]  = RPID;
                     GotYou = true;
                  }
               }
            } // if ( SibPID >= 0  &&  amr->patch[0][Lv][SibPID]->son != -1 )
         } // for (int s=0; s<6; s++)

         if ( GotYou )  LB_RecvX_NList[r] ++;
      } // for (int r=0; r<MPI_NRank; r++) for (int t=0; t<LB_RecvH_NList[r]; t++)

   } // if ( OPT__FIXUP_FLUX )



// 4. B field data to be sent and received after the electric field fix-up operation
// ============================================================================================================
#  ifdef MHD
   if ( OPT__FIXUP_ELECTRIC )
   {
//    check
      if ( !amr->WithElectric )
         Aux_Error( ERROR_INFO, "amr->WithElectric is off --> no electric field array is required !!\n" );


//    4.1 send
      for (int r=0; r<MPI_NRank; r++)
      for (int t=0; t<LB_SendH_NList[r]; t++)
      {
         const int SPID = LB_SendH_IDList[r][t];

         for (int s=0; s<18; s++)
         {
//          for simplicity, just re-send the data of patches with electric field arrays being allocated
            if ( amr->patch[0][Lv][SPID]->electric[s] != NULL )
            {
               LB_SendY_IDList [r][ LB_SendY_NList[r] ] = SPID;
               LB_SendY_SibList[r][ LB_SendY_NList[r] ] = LB_SendH_SibList[r][t];
               LB_SendY_SubID  [r][ LB_SendY_NList[r] ] = t;
               LB_SendY_NList  [r] ++;
               break;
            }
         }
      }


//    4.2 recv
      int  SendY_Disp[MPI_NRank], RecvY_Disp[MPI_NRank], NSendY_Total, NRecvY_Total, CounterY;
      int *SendY_Buf=NULL, *RecvY_Buf=NULL, *RecvY_Ptr=NULL;

//    4.2.1. set the N list
      MPI_Alltoall( LB_SendY_NList, 1, MPI_INT, LB_RecvY_NList, 1, MPI_INT, MPI_COMM_WORLD );

//    4.2.2. exchange the subset ID list
      SendY_Disp[0] = 0;
      RecvY_Disp[0] = 0;
      for (int r=1; r<MPI_NRank; r++)
      {
         SendY_Disp[r] = SendY_Disp[r-1] + LB_SendY_NList[r-1];
         RecvY_Disp[r] = RecvY_Disp[r-1] + LB_RecvY_NList[r-1];
      }
      NSendY_Total = SendY_Disp[MPI_NRank-1] + LB_SendY_NList[MPI_NRank-1];
      NRecvY_Total = RecvY_Disp[MPI_NRank-1] + LB_RecvY_NList[MPI_NRank-1];

      SendY_Buf = new int [NSendY_Total];
      RecvY_Buf = new int [NRecvY_Total];

      CounterY = 0;
      for (int r=0; r<MPI_NRank; r++)
      for (int t=0; t<LB_SendY_NList[r]; t++)
         SendY_Buf[ CounterY ++ ] = LB_SendY_SubID[r][t];

      MPI_Alltoallv( SendY_Buf, LB_SendY_NList, SendY_Disp, MPI_INT,
                     RecvY_Buf, LB_RecvY_NList, RecvY_Disp, MPI_INT,  MPI_COMM_WORLD );

//    4.2.3. set the recv lists
      for (int r=0; r<MPI_NRank; r++)
      {
         RecvY_Ptr = RecvY_Buf + RecvY_Disp[r];

         for (int t=0; t<LB_RecvY_NList[r]; t++)
         {
            const int SubID = RecvY_Ptr[t];

            LB_RecvY_IDList [r][t] = LB_RecvH_IDList [r][ LB_RecvH_IDList_IdxTable[r][SubID] ];
            LB_RecvY_SibList[r][t] = LB_RecvH_SibList[r][SubID];
         }
      }

      delete [] SendY_Buf;
      delete [] RecvY_Buf;
   } // if ( OPT__FIXUP_ELECTRIC )

   for (int r=0; r<MPI_NRank; r++)  delete [] LB_SendY_SubID[r];
   delete [] LB_SendY_SubID;
#  endif // #ifdef MHD



// 5. check
// ============================================================================================================
#  ifdef GAMER_DEBUG
   int  LB_RecvX_NList_Ck[MPI_NRank], LB_RecvX_NResList_Ck[MPI_NRank], Send_Disp[MPI_NRank], Recv_Disp[MPI_NRank];
   int  NSend_Total, NRecv_Total, Counter;
   int *SendBuf, *RecvBuf, *RecvPtr;

// 5.1. check the N lists
   MPI_Alltoall( LB_SendX_NList,    1, MPI_INT, LB_RecvX_NList_Ck,    1, MPI_INT, MPI_COMM_WORLD );
   MPI_Alltoall( LB_SendX_NResList, 1, MPI_INT, LB_RecvX_NResList_Ck, 1, MPI_INT, MPI_COMM_WORLD );

   for (int r=0; r<MPI_NRank; r++)
   {
      if ( LB_RecvX_NList[r] != LB_RecvX_NList_Ck[r] )
         Aux_Error( ERROR_INFO, "Lv %d, TRank %d, LB_RecvX_NList (%d) != Expect (%d) !!\n",
                    Lv, r, LB_RecvX_NList[r], LB_RecvX_NList_Ck[r] );

      if ( LB_RecvX_NResList[r] != LB_RecvX_NResList_Ck[r] )
         Aux_Error( ERROR_INFO, "Lv %d, TRank %d, LB_RecvX_NResList (%d) != Expect (%d) !!\n",
                    Lv, r, LB_RecvX_NResList[r], LB_RecvX_NResList_Ck[r] );
   }


// 5.2. check the recv sibling lists
   Send_Disp[0] = 0;
   Recv_Disp[0] = 0;
   for (int r=1; r<MPI_NRank; r++)
   {
      Send_Disp[r] = Send_Disp[r-1] + LB_SendX_NList[r-1];
      Recv_Disp[r] = Recv_Disp[r-1] + LB_RecvX_NList[r-1];
   }
   NSend_Total = Send_Disp[MPI_NRank-1] + LB_SendX_NList[MPI_NRank-1];
   NRecv_Total = Recv_Disp[MPI_NRank-1] + LB_RecvX_NList[MPI_NRank-1];

   SendBuf = new int [NSend_Total];
   RecvBuf = new int [NRecv_Total];

   Counter = 0;
   for (int r=0; r<MPI_NRank; r++)
   for (int t=0; t<LB_SendX_NList[r]; t++)
      SendBuf[ Counter ++ ] = LB_SendX_SibList[r][t];

   MPI_Alltoallv( SendBuf, LB_SendX_NList, Send_Disp, MPI_INT,
                  RecvBuf, LB_RecvX_NList, Recv_Disp, MPI_INT,  MPI_COMM_WORLD );

   for (int r=0; r<MPI_NRank; r++)
   {
      RecvPtr = RecvBuf + Recv_Disp[r];

      for (int t=0; t<LB_RecvX_NList[r]; t++)
      {
         if ( RecvPtr[t] != LB_RecvX_SibList[r][t] )
         {
            Aux_Message( stderr, "ERROR : lv %d, MyRank %d: RecvX_SibList[%d][%d] (%d) != Expect (%d) !!\n",
                         Lv, MPI_Rank, r, t, LB_RecvX_SibList[r][t], RecvPtr[t] );

            Aux_Message( stderr, "        LB_RecvX_NList = %d, LB_RecvX_NResList = %d\n",
                         LB_RecvX_NList[r], LB_RecvX_NResList[r] );

            Aux_Message( stderr, "        RecvX_SibList = " );
            for (int b=0; b<32; b++)
               Aux_Message( stderr, "%d", ( LB_RecvX_SibList[r][t] & (1<<b) ) ? 1 : 0 );
            Aux_Message( stderr, "\n" );

            Aux_Message( stderr, "        Expect        = " );
            for (int b=0; b<32; b++)
               Aux_Message( stderr, "%d", ( RecvPtr[t] & (1<<b) ) ? 1 : 0 );
            Aux_Message( stderr, "\n" );

            Aux_Error( ERROR_INFO, "Program is going to be terminated ...\n" );
         }
      } // for (int t=0; t<LB_RecvX_NList[r]; t++)
   } // for (int r=0; r<MPI_NRank; r++)

   delete [] SendBuf;
   delete [] RecvBuf;
#  endif // #ifdef GAMER_DEBUG

} // FUNCTION : LB_RecordExchangeFixUpDataPatchID



#endif // #ifdef LOAD_BALANCE
