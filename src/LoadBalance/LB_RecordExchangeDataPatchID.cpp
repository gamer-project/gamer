#include "GAMER.h"

#ifdef LOAD_BALANCE



static void SetTargetSibling( int NTSib[], int* TSib_List[] );
static void SetSiblingMask( int SibMask_Check[], int SibMask_Clear[], int SibMask_Duplicate[] );
static void SetReceiveSibling( int* RSib_List[] );




//-------------------------------------------------------------------------------------------------------
// Function    :  LB_RecordExchangeDataPatchID
// Description :  Construct the MPI sending and receiving data lists for exchanging fluid, magnetic field,
//                and potential data
//
// Note        :  1. LB_RecvH_IDList[] is unsorted --> use LB_RecvH_IDList_Idxtable[] to obtain the correct order
//                   <--> All other lists are sorted
//                2. This function will NOT deallocate any fluid/magnetic/pot arrays allocated previously
//
// Parameter   :  Lv          : Target refinement level for recording MPI lists
//                AfterRefine : Record the difference between old and new MPI lists after grid refinement
//                              --> Minimizing the MPI time after grid refinement by only exchanging the
//                                  buffer data not existing in the old MPI lists
//-------------------------------------------------------------------------------------------------------
void LB_RecordExchangeDataPatchID( const int Lv, const bool AfterRefine )
{

//###OPTIMIZATION: NSib_C = 6 for some interpolation schemes
   const int MirSib[27] = { 1,0,3,2,5,4,9,8,7,6,13,12,11,10,17,16,15,14,25,24,23,22,21,20,19,18,26 };
   const int SonLv      = Lv + 1;
   const int NSib_F     = 26;
   const int NSib_C     = 26;
   const int NReal      = amr->NPatchComma[Lv][1];
   const int NBuff      = amr->NPatchComma[Lv][3] - amr->NPatchComma[Lv][1];
   const int MemUnit_H  = 1 + NBuff/MPI_NRank;           // set arbitrarily
#  ifdef GRAVITY
   const int MemUnit_G  = 1 + (NReal+NBuff)/MPI_NRank;   // set arbitrarily
#  endif

   int  FaPID, FaSibPID, FaSibPID0, SibPID=-1, SibPID0, TRank=-1, SibPID0_List[26], NSibPID_Delta[26], *SibPID_Delta[26];
   int  NFaBuff, FaBuff[27], TPID, NTSib[26], *TSib_List[26], Side, RSib, *RSib_List[26], SibIdx, NID, OID;
   int  SibMask_Check[27], SibMask_Clear[27], SibMask_Duplicate[27];
   long SibLBIdx;

   int   MemSize_H[MPI_NRank], *LB_RecvH_SibList_Unsorted[MPI_NRank], *SibList_H, *Match_H;
   int   Old_RecvH_NList[MPI_NRank], *Old_RecvH_SibList[MPI_NRank], *Old_RecvH_PCr1D_IdxTable[MPI_NRank];
  ulong *Old_RecvH_PCr1D[MPI_NRank];

   int   *LB_SendH_NList           = amr->LB->SendH_NList          [Lv];
   int  **LB_SendH_IDList          = amr->LB->SendH_IDList         [Lv];
   int  **LB_SendH_SibList         = amr->LB->SendH_SibList        [Lv];
   int  **LB_SendH_SibDiffList     = amr->LB->SendH_SibDiffList    [Lv];
   long **LB_SendH_LBIdxList       = amr->LB->SendH_LBIdxList      [Lv];
   int   *LB_RecvH_NList           = amr->LB->RecvH_NList          [Lv];
   int  **LB_RecvH_IDList          = amr->LB->RecvH_IDList         [Lv];
   int  **LB_RecvH_IDList_IdxTable = amr->LB->RecvH_IDList_IdxTable[Lv];
   int  **LB_RecvH_SibList         = amr->LB->RecvH_SibList        [Lv];
   int  **LB_RecvH_SibDiffList     = amr->LB->RecvH_SibDiffList    [Lv];
   long **LB_RecvH_LBIdxList       = amr->LB->RecvH_LBIdxList      [Lv];
  ulong **LB_RecvH_PCr1D           = amr->LB->RecvH_PCr1D          [Lv];
   int  **LB_RecvH_PCr1D_IdxTable  = amr->LB->RecvH_PCr1D_IdxTable [Lv];

#  ifdef GRAVITY
   int   MemSize_G[MPI_NRank], *LB_RecvG_SibList_Unsorted[MPI_NRank], *SibList_G, *Match_G;
   int   Old_RecvG_NList[MPI_NRank], *Old_RecvG_SibList[MPI_NRank], *Old_RecvG_PCr1D_IdxTable[MPI_NRank];
  ulong *Old_RecvG_PCr1D[MPI_NRank];

   int   *LB_SendG_NList           = amr->LB->SendG_NList          [Lv];
   int  **LB_SendG_IDList          = amr->LB->SendG_IDList         [Lv];
   int  **LB_SendG_SibList         = amr->LB->SendG_SibList        [Lv];
   int  **LB_SendG_SibDiffList     = amr->LB->SendG_SibDiffList    [Lv];
   long **LB_SendG_LBIdxList       = amr->LB->SendG_LBIdxList      [Lv];
   int   *LB_RecvG_NList           = amr->LB->RecvG_NList          [Lv];
   int  **LB_RecvG_IDList          = amr->LB->RecvG_IDList         [Lv];
   int  **LB_RecvG_IDList_IdxTable = amr->LB->RecvG_IDList_IdxTable[Lv];
   int  **LB_RecvG_SibList         = amr->LB->RecvG_SibList        [Lv];
   int  **LB_RecvG_SibDiffList     = amr->LB->RecvG_SibDiffList    [Lv];
   long **LB_RecvG_LBIdxList       = amr->LB->RecvG_LBIdxList      [Lv];
  ulong **LB_RecvG_PCr1D           = amr->LB->RecvG_PCr1D          [Lv];
   int  **LB_RecvG_PCr1D_IdxTable  = amr->LB->RecvG_PCr1D_IdxTable [Lv];
#  endif


// 1. initialize arrays
// ============================================================================================================
   for (int r=0; r<MPI_NRank; r++)
   {
      Old_RecvH_NList          [r] = LB_RecvH_NList         [r];
      Old_RecvH_SibList        [r] = LB_RecvH_SibList       [r];
      Old_RecvH_PCr1D          [r] = LB_RecvH_PCr1D         [r];
      Old_RecvH_PCr1D_IdxTable [r] = LB_RecvH_PCr1D_IdxTable[r];

      MemSize_H                [r] = MemUnit_H;
      LB_RecvH_NList           [r] = 0;
      LB_RecvH_IDList          [r] = (int*)realloc( LB_RecvH_IDList[r], MemSize_H[r]*sizeof(int) );
      LB_RecvH_SibList_Unsorted[r] = (int*)malloc( MemSize_H[r]*sizeof(int) );

#     ifdef GRAVITY
      Old_RecvG_NList          [r] = LB_RecvG_NList         [r];
      Old_RecvG_SibList        [r] = LB_RecvG_SibList       [r];
      Old_RecvG_PCr1D          [r] = LB_RecvG_PCr1D         [r];
      Old_RecvG_PCr1D_IdxTable [r] = LB_RecvG_PCr1D_IdxTable[r];

      MemSize_G                [r] = MemUnit_G;
      LB_RecvG_NList           [r] = 0;
      LB_RecvG_IDList          [r] = (int*)realloc( LB_RecvG_IDList[r], MemSize_G[r]*sizeof(int) );
      LB_RecvG_SibList_Unsorted[r] = (int*)malloc( MemSize_G[r]*sizeof(int) );
#     endif
   } // for (int r=0; r<MPI_NRank; r++)

// get the sibling index differences along different directions
   TABLE_GetSibPID_Delta( NSibPID_Delta, SibPID_Delta );

// set up the target sibling indices at SonLv
   SetTargetSibling( NTSib, TSib_List );

// set up the bit mask arrays
   SetSiblingMask( SibMask_Check, SibMask_Clear, SibMask_Duplicate );

// set up the receiving sibling indices at SonLv
   SetReceiveSibling( RSib_List );

// initialize the sibling-indices list
   SibList_H = (int*)calloc( NBuff, sizeof(int) );
#  ifdef GRAVITY
   SibList_G = (int*)calloc( NBuff, sizeof(int) );
#  endif



// 2. get the sibling indices of all data to be received
// ============================================================================================================
// 2.1 sibling patches at Lv (for fluid solver, gravity solver, and estimating time-step by gravity)
   for (int PID0=0; PID0<NReal; PID0+=8)  // loop over all "real" patches at Lv with LocalID == 0
   {
      TABLE_GetSibPID_Based( Lv, PID0, SibPID0_List );

      for (int s=0; s<NSib_F; s++)
      {
         SibPID0 = SibPID0_List[s];

         if ( SibPID0 >= NReal ) // work for both periodic and non-periodic boundary conditions
         {
            RSib = MirSib[s];    // sibling direction to receive data

            for (int Count=0; Count<NSibPID_Delta[s]; Count++)
            {
               SibPID = SibPID0 + SibPID_Delta[s][Count];
               SibIdx = SibPID - NReal;

//             2.1.1 hydro and magnetic field
//             record sibling indices
               if (  ( SibList_H[SibIdx] & SibMask_Check[RSib] ) == false  )
               {
                  SibList_H[SibIdx] |= ( 1 << RSib );
                  SibList_H[SibIdx] &= SibMask_Clear[RSib];
               }

//             allocate memory for the buffer patches that will receive data
               for (int Sg=0; Sg<2; Sg++)    amr->patch[Sg][Lv][SibPID]->hnew();

#              ifdef MHD
               for (int Sg=0; Sg<2; Sg++)    amr->patch[Sg][Lv][SibPID]->mnew();
#              endif

#              ifdef GRAVITY // so that the XXX_H lists can also be applied to the potential data
               for (int Sg=0; Sg<2; Sg++)    amr->patch[Sg][Lv][SibPID]->gnew();
#              endif


//             2.1.2 potential
//###OPTIMIZATION: these buffer patches are allocated for the gravity solver at the base level and the
//                 time-step estimation --> only need to exchange "GRA_GHOST_SIZE" cells
#              ifdef GRAVITY
//             record sibling indices
               if (  ( SibList_G[SibIdx] & SibMask_Check[RSib] ) == false  )
               {
                  SibList_G[SibIdx] |= ( 1 << RSib );
                  SibList_G[SibIdx] &= SibMask_Clear[RSib];
               }
#              endif

            } // for (int Count=0; Count<NSibPID_Delta[s]; Count++)
         } // if ( SibPID0 >= NReal )
      } // for (int s=0; s<NSib_F; s++)
   } // for (int PID0=0; PID0<NReal; PID0+=8)


// 2.2 father-sibling patches at Lv (for fluid solver only)
   if ( SonLv < NLEVEL )
   for (int SonPID0=0; SonPID0<amr->NPatchComma[SonLv][1]; SonPID0+=8)  // loop over all "real" patches at SonLv
   {
      FaPID = amr->patch[0][SonLv][SonPID0]->father;

#     ifdef GAMER_DEBUG
      if ( FaPID == -1 )   Aux_Error( ERROR_INFO, "SonLv %d, SonPID0 %d has no father patch !!\n", SonLv, SonPID0 );
#     endif

//    exchange complete father patch if option OPT__LB_EXCHANGE_FATHER is set
      if ( OPT__LB_EXCHANGE_FATHER ) {

         TPID = FaPID;

         if ( TPID >= NReal )
         {
            RSib = 26;
            SibIdx = TPID - NReal;

            if (  ( SibList_H[SibIdx] & SibMask_Check[RSib] ) == false  )
            {
               SibList_H[SibIdx] |= ( 1 << RSib );
               SibList_H[SibIdx] &= SibMask_Clear[RSib];
            }

//          allocate memory for the buffer patches that will receive data
            for (int Sg=0; Sg<2; Sg++)    amr->patch[Sg][Lv][TPID]->hnew();

#           ifdef MHD
            for (int Sg=0; Sg<2; Sg++)    amr->patch[Sg][Lv][TPID]->mnew();
#           endif

#           ifdef GRAVITY // so that the XXX_H lists can also be applied to the potential data
            for (int Sg=0; Sg<2; Sg++)    amr->patch[Sg][Lv][TPID]->gnew();
#           endif
         } // if ( TPID >= NReal )
      } // if ( OPT__LB_EXCHANGE_FATHER )

      TABLE_GetSibPID_Based( SonLv, SonPID0, SibPID0_List );

      for (int s=0; s<NSib_F; s++)
      {
         SibPID0 = SibPID0_List[s];

         if ( SibPID0 == -1 ) // work for both periodic and non-periodic boundary conditions
         {
            NFaBuff   = 0;
            FaSibPID0 = amr->patch[0][Lv][FaPID]->sibling[s];
            FaBuff[ NFaBuff ++ ] = FaSibPID0;

#           ifdef GAMER_DEBUG
//          FaSibPID0 < -1 is disallowed even for non-periodic B.C. (since SibPID0 == -1 instead of < -1)
            if ( FaSibPID0 < 0 )
               Aux_Error( ERROR_INFO, "Lv %d, SonPID0 %d, FaPID %d has no sibling [%d] !!\n",
                          Lv, SonPID0, FaPID, s );
#           endif

            for (int CSib=0; CSib<NTSib[s]; CSib++)
            {
               Side     = TSib_List[s][CSib];
               FaSibPID = amr->patch[0][Lv][FaSibPID0]->sibling[Side];
               FaBuff[ NFaBuff ++ ] = FaSibPID;

#              ifdef GAMER_DEBUG
//             FaSibPID < -1 is allowed for non-periodic B.C. even though FaSibPID0 < -1 is disallowed
               if ( FaSibPID == -1 )
                  Aux_Error( ERROR_INFO, "Lv %d, FaSibPID0 %d has no sibling [%d] !!\n",
                             Lv, FaSibPID0, Side );
#              endif
            }

            for (int f=0; f<NFaBuff; f++)
            {
               TPID = FaBuff[f];

               if ( TPID >= NReal ) // work for both periodic and non-periodic boundary conditions
               {
//                record sibling indices
                  RSib   = RSib_List[s][f];  // sibling direction to receive data
                  SibIdx = TPID - NReal;

                  if (  ( SibList_H[SibIdx] & SibMask_Check[RSib] ) == false  )
                  {
                     SibList_H[SibIdx] |= ( 1 << RSib );
                     SibList_H[SibIdx] &= SibMask_Clear[RSib];
                  }

//                allocate memory for the buffer patches that will receive data
                  for (int Sg=0; Sg<2; Sg++)    amr->patch[Sg][Lv][TPID]->hnew();

#                 ifdef MHD
                  for (int Sg=0; Sg<2; Sg++)    amr->patch[Sg][Lv][TPID]->mnew();
#                 endif

#                 ifdef GRAVITY // so that the XXX_H lists can also be applied to the potential data
                  for (int Sg=0; Sg<2; Sg++)    amr->patch[Sg][Lv][TPID]->gnew();
#                 endif
               } // if ( TPID >= NReal )
            } // for (int f=0; f<NFaBuff; f++)
         } // if ( SibPID0 == -1 )
      } // for (int s=0; s<NSib_F; s++)
   } // for (int SonPID0=0; SonPID0<amr->NPatchComma[SonLv][1]; SonPID0+=8)


// 2.3 father-sibling patches at Lv (for Poisson solver only)
#  ifdef GRAVITY
   if ( SonLv < NLEVEL )
   for (int SonPID0=0; SonPID0<amr->NPatchComma[SonLv][1]; SonPID0+=8)  // loop over all "real" patches at SonLv
   {
      FaPID = amr->patch[0][SonLv][SonPID0]->father;

#     ifdef GAMER_DEBUG
      if ( FaPID == -1 )   Aux_Error( ERROR_INFO, "SonLv %d, SonPID0 %d has no father patch !!\n", SonLv, SonPID0 );
#     endif

      for (int s=0; s<NSib_C; s++)  FaBuff[s] = amr->patch[0][Lv][FaPID]->sibling[s];
      FaBuff[NSib_C] = FaPID;    // here we assume NSib_C == 26

#     ifdef GAMER_DEBUG
      for (int s=0; s<NSib_C+1; s++)
      {
         if ( FaBuff[s] == -1 )
            Aux_Error( ERROR_INFO, "Lv %d, SonPID0 %d, FaPID %d, FaBuff[%d] == -1 !!\n", Lv, SonPID0, FaPID, s );
      }
#     endif

      for (int f=0; f<NSib_C+1; f++)
      {
         TPID = FaBuff[f];

         if ( TPID >= NReal ) // work for both periodic and non-periodic boundary conditions
         {
//          record sibling indices
            RSib   = MirSib[f];     // sibling direction to receive data
            SibIdx = TPID - NReal;

            if (  ( SibList_G[SibIdx] & SibMask_Check[RSib] ) == false  )
            {
               SibList_G[SibIdx] |= ( 1 << RSib );
               SibList_G[SibIdx] &= SibMask_Clear[RSib];
            }

//          allocate memory for the buffer patches that will receive data
            for (int Sg=0; Sg<2; Sg++)    amr->patch[Sg][Lv][TPID]->gnew();

         } // if ( TPID >= NReal )
      } // for (int f=0; f<NFaBuff; f++)
   } // for (int SonPID0=0; SonPID0<amr->NPatchComma[SonLv][1]; SonPID0+=8)
#  endif // #ifdef GRAVITY


// check if there are duplicate target data in SibList
#  ifdef GAMER_DEBUG
   for (int t=0; t<NBuff; t++)
   for (int s=0; s<27; s++)
   {
      if (  ( SibList_H[t] & (1<<s) )  &&  ( SibList_H[t] & SibMask_Duplicate[s] )  )
      {
         Aux_Message( stderr, "ERROR : duplicate data in SibList_H[%d], Lv %d, Rank %d, PID %d, s %d !!\n",
                      t, Lv, MPI_Rank, t+NReal, s );

         Aux_Message( stderr, "        SibList_H   = " );
         for (int b=0; b<32; b++)
            Aux_Message( stderr, "%d", ( SibList_H[t] & (1<<b) ) ? 1 : 0 );
         Aux_Message( stderr, "\n" );

         Aux_Message( stderr, "        SibMask_Dup = " );
         for (int b=0; b<32; b++)
            Aux_Message( stderr, "%d", ( SibMask_Duplicate[s] & (1<<b) ) ? 1 : 0 );
         Aux_Message( stderr, "\n" );

         Aux_Error( ERROR_INFO, "Program is going to be terminated ...\n" );
      }

#     ifdef GRAVITY
      if (  ( SibList_G[t] & (1<<s) )  &&  ( SibList_G[t] & SibMask_Duplicate[s] )  )
      {
         Aux_Message( stderr, "ERROR : duplicate data in SibList_G[%d], Lv %d, Rank %d, PID %d, s %d !!\n",
                      t, Lv, MPI_Rank, t+NReal, s );

         Aux_Message( stderr, "        SibList_G   = " );
         for (int b=0; b<32; b++)
            Aux_Message( stderr, "%d", ( SibList_G[t] & (1<<b) ) ? 1 : 0 );
         Aux_Message( stderr, "\n" );

         Aux_Message( stderr, "        SibMask_Dup = " );
         for (int b=0; b<32; b++)
            Aux_Message( stderr, "%d", ( SibMask_Duplicate[s] & (1<<b) ) ? 1 : 0 );
         Aux_Message( stderr, "\n" );

         Aux_Error( ERROR_INFO, "Program is going to be terminated ...\n" );
      }
#     endif
   } // for (int t=0; t<NBuff; t++) for (int s=0; s<27; s++)
#  endif // #ifdef GAMER_DEBUG



// 3. get the unsorted recv list
// ============================================================================================================
   for (int t=0; t<NBuff; t++)
   {
//    3.1 hydro and magnetic field
      if ( SibList_H[t] )
      {
         SibPID   = t + NReal;
         SibLBIdx = amr->patch[0][Lv][SibPID]->LB_Idx;
         TRank    = LB_Index2Rank( Lv, SibLBIdx, CHECK_ON );

//       allocate enough memory
         if ( LB_RecvH_NList[TRank] >= MemSize_H[TRank] )
         {
            MemSize_H                [TRank] += MemUnit_H;
            LB_RecvH_IDList          [TRank]  = (int*)realloc( LB_RecvH_IDList          [TRank],
                                                               MemSize_H[TRank]*sizeof(int) );
            LB_RecvH_SibList_Unsorted[TRank]  = (int*)realloc( LB_RecvH_SibList_Unsorted[TRank],
                                                               MemSize_H[TRank]*sizeof(int) );
         }

         LB_RecvH_IDList          [TRank][ LB_RecvH_NList[TRank] ] = SibPID;
         LB_RecvH_SibList_Unsorted[TRank][ LB_RecvH_NList[TRank] ] = SibList_H[t];
         LB_RecvH_NList[TRank] ++;
      }


//    3.2 potential
#     ifdef GRAVITY
      if ( SibList_G[t] )
      {
//       get TRank only if it is not calculated yet
         if ( SibList_H[t] == 0 )
         {
            SibPID   = t + NReal;
            SibLBIdx = amr->patch[0][Lv][SibPID]->LB_Idx;
            TRank    = LB_Index2Rank( Lv, SibLBIdx, CHECK_ON );
         }

//       allocate enough memory
         if ( LB_RecvG_NList[TRank] >= MemSize_G[TRank] )
         {
            MemSize_G                [TRank] += MemUnit_G;
            LB_RecvG_IDList          [TRank]  = (int*)realloc( LB_RecvG_IDList          [TRank],
                                                               MemSize_G[TRank]*sizeof(int) );
            LB_RecvG_SibList_Unsorted[TRank]  = (int*)realloc( LB_RecvG_SibList_Unsorted[TRank],
                                                               MemSize_G[TRank]*sizeof(int) );
         }

         LB_RecvG_IDList          [TRank][ LB_RecvG_NList[TRank] ] = SibPID;
         LB_RecvG_SibList_Unsorted[TRank][ LB_RecvG_NList[TRank] ] = SibList_G[t];
         LB_RecvG_NList[TRank] ++;
      }
#     endif // #ifdef GRAVITY
   } // for (int t=0; t<NBuff; t++)



// 4. get the MPI displacement arrays
// ============================================================================================================
// 4.1 hydro and magnetic field
   int Send_Disp_H[MPI_NRank], Recv_Disp_H[MPI_NRank], NSend_Total_H, NRecv_Total_H;

// 4.1.1 broadcast the number of elements received from different ranks
   MPI_Alltoall( LB_RecvH_NList, 1, MPI_INT, LB_SendH_NList, 1, MPI_INT, MPI_COMM_WORLD );

// 4.1.2 get the displacement and total count
   Send_Disp_H[0] = 0;
   Recv_Disp_H[0] = 0;
   for (int r=1; r<MPI_NRank; r++)
   {
      Send_Disp_H[r] = Send_Disp_H[r-1] + LB_RecvH_NList[r-1];
      Recv_Disp_H[r] = Recv_Disp_H[r-1] + LB_SendH_NList[r-1];
   }
   NSend_Total_H = Send_Disp_H[MPI_NRank-1] + LB_RecvH_NList[MPI_NRank-1];
   NRecv_Total_H = Recv_Disp_H[MPI_NRank-1] + LB_SendH_NList[MPI_NRank-1];


// 4.2 potential
#  ifdef GRAVITY
   int Send_Disp_G[MPI_NRank], Recv_Disp_G[MPI_NRank], NSend_Total_G, NRecv_Total_G;

// 4.2.1 broadcast the number of elements received from different ranks
   MPI_Alltoall( LB_RecvG_NList, 1, MPI_INT, LB_SendG_NList, 1, MPI_INT, MPI_COMM_WORLD );

// 4.2.2 prepare the MPI buffers
   Send_Disp_G[0] = 0;
   Recv_Disp_G[0] = 0;
   for (int r=1; r<MPI_NRank; r++)
   {
      Send_Disp_G[r] = Send_Disp_G[r-1] + LB_RecvG_NList[r-1];
      Recv_Disp_G[r] = Recv_Disp_G[r-1] + LB_SendG_NList[r-1];
   }
   NSend_Total_G = Send_Disp_G[MPI_NRank-1] + LB_RecvG_NList[MPI_NRank-1];
   NRecv_Total_G = Recv_Disp_G[MPI_NRank-1] + LB_SendG_NList[MPI_NRank-1];
#  endif // #ifdef GRAVITY



// 5. sort the recv list
// ============================================================================================================
// 5.1 hydro and magnetic field
// 5.1.1 allocate memory
   if ( LB_RecvH_IDList_IdxTable[0] != NULL )   delete [] LB_RecvH_IDList_IdxTable[0];
   if ( LB_RecvH_SibDiffList    [0] != NULL )   delete [] LB_RecvH_SibDiffList    [0];
   if ( LB_RecvH_LBIdxList      [0] != NULL )   delete [] LB_RecvH_LBIdxList      [0];
   if ( LB_SendH_IDList         [0] != NULL )   delete [] LB_SendH_IDList         [0];
   if ( LB_SendH_SibList        [0] != NULL )   delete [] LB_SendH_SibList        [0];
   if ( LB_SendH_SibDiffList    [0] != NULL )   delete [] LB_SendH_SibDiffList    [0];
   if ( LB_SendH_LBIdxList      [0] != NULL )   delete [] LB_SendH_LBIdxList      [0];

   LB_RecvH_IDList_IdxTable[0] = new int  [NSend_Total_H];
   LB_RecvH_SibList        [0] = new int  [NSend_Total_H];
   LB_RecvH_SibDiffList    [0] = new int  [NSend_Total_H];
   LB_RecvH_LBIdxList      [0] = new long [NSend_Total_H];
   LB_RecvH_PCr1D          [0] = new ulong[NSend_Total_H];
   LB_RecvH_PCr1D_IdxTable [0] = new int  [NSend_Total_H];
   LB_SendH_IDList         [0] = new int  [NRecv_Total_H];
   LB_SendH_SibList        [0] = new int  [NRecv_Total_H];
   LB_SendH_SibDiffList    [0] = new int  [NRecv_Total_H];
   LB_SendH_LBIdxList      [0] = new long [NRecv_Total_H];

   for (int r=0; r<MPI_NRank; r++)
   {
//    5.1.2 assign pointers
      LB_RecvH_IDList_IdxTable[r] = LB_RecvH_IDList_IdxTable[0] + Send_Disp_H[r];
      LB_RecvH_SibList        [r] = LB_RecvH_SibList        [0] + Send_Disp_H[r];
      LB_RecvH_SibDiffList    [r] = LB_RecvH_SibDiffList    [0] + Send_Disp_H[r];
      LB_RecvH_LBIdxList      [r] = LB_RecvH_LBIdxList      [0] + Send_Disp_H[r];
      LB_RecvH_PCr1D          [r] = LB_RecvH_PCr1D          [0] + Send_Disp_H[r];
      LB_RecvH_PCr1D_IdxTable [r] = LB_RecvH_PCr1D_IdxTable [0] + Send_Disp_H[r];
      LB_SendH_IDList         [r] = LB_SendH_IDList         [0] + Recv_Disp_H[r];
      LB_SendH_SibList        [r] = LB_SendH_SibList        [0] + Recv_Disp_H[r];
      LB_SendH_SibDiffList    [r] = LB_SendH_SibDiffList    [0] + Recv_Disp_H[r];
      LB_SendH_LBIdxList      [r] = LB_SendH_LBIdxList      [0] + Recv_Disp_H[r];

//    5.1.3 get the sorted LBIdx list
      for (int t=0; t<LB_RecvH_NList[r]; t++)
         LB_RecvH_LBIdxList[r][t] = amr->patch[0][Lv][ LB_RecvH_IDList[r][t] ]->LB_Idx;

      Mis_Heapsort( LB_RecvH_NList[r], LB_RecvH_LBIdxList[r], LB_RecvH_IDList_IdxTable[r] );

//    5.1.4 get the sorted SibList and PaddedCr1D
      for (int t=0; t<LB_RecvH_NList[r]; t++)
      {
         TPID = LB_RecvH_IDList[r][ LB_RecvH_IDList_IdxTable[r][t] ];

         LB_RecvH_SibList[r][t] = LB_RecvH_SibList_Unsorted[r][ LB_RecvH_IDList_IdxTable[r][t] ];
         LB_RecvH_PCr1D  [r][t] = amr->patch[0][Lv][TPID]->PaddedCr1D;
      }

//    5.1.5 sort PaddedCr1D again and record the index table (for constructing the SibDiff lists later)
      Mis_Heapsort( LB_RecvH_NList[r], LB_RecvH_PCr1D[r], LB_RecvH_PCr1D_IdxTable[r] );

//    5.1.6 get the recv SibDiff lists
      if ( AfterRefine )
      {
         Match_H = new int [ LB_RecvH_NList[r] ];

//       match patches with the same LB_Idx
         Mis_Matching_int( Old_RecvH_NList[r], Old_RecvH_PCr1D[r], LB_RecvH_NList[r], LB_RecvH_PCr1D[r], Match_H );

//       remove sibling directions already exist in the previous MPI lists
         for (int t=0; t<LB_RecvH_NList[r]; t++)
         {
            NID = LB_RecvH_PCr1D_IdxTable[r][t];

            if ( Match_H[t] == -1 )
               LB_RecvH_SibDiffList[r][NID] = LB_RecvH_SibList[r][NID];

            else
            {
               OID = Old_RecvH_PCr1D_IdxTable[r][ Match_H[t] ];

//             add duplicate sibling directions to the old MPI lists (to minimize the amount of data to be resent)
               if ( Old_RecvH_SibList[r][OID] & (1<<26) )   Old_RecvH_SibList[r][OID] |= SibMask_Duplicate[26];
               else
               for (int s=0; s<18; s++)
               if ( Old_RecvH_SibList[r][OID] & (1<< s) )   Old_RecvH_SibList[r][OID] |= SibMask_Duplicate[ s];

//             record SibDiff List
               LB_RecvH_SibDiffList[r][NID] = LB_RecvH_SibList[r][NID] & ( ~Old_RecvH_SibList[r][OID] );
            }
         }

         delete [] Match_H;
      } // if ( AfterRefine )

//    recv PID should not repeat
#     ifdef GAMER_DEBUG
      int *TempIDList = new int [ LB_RecvH_NList[r] ];
      memcpy( TempIDList, LB_RecvH_IDList[r], LB_RecvH_NList[r]*sizeof(int) );

      Mis_Heapsort<int,int>( LB_RecvH_NList[r], TempIDList, NULL );

      for (int t=1; t<LB_RecvH_NList[r]; t++)
         if ( TempIDList[t] == TempIDList[t-1] )
            Aux_Error( ERROR_INFO, "Lv %d, Rank %d, repeated recv fluid PID %d !!\n",
                       Lv, r, TempIDList[t] );

      delete [] TempIDList;
#     endif
   } // for (int r=0; r<MPI_NRank; r++)


// 5.2 potential
#  ifdef GRAVITY
// 5.2.1 allocate memory
   if ( LB_RecvG_IDList_IdxTable[0] != NULL )   delete [] LB_RecvG_IDList_IdxTable[0];
   if ( LB_RecvG_SibDiffList    [0] != NULL )   delete [] LB_RecvG_SibDiffList    [0];
   if ( LB_RecvG_LBIdxList      [0] != NULL )   delete [] LB_RecvG_LBIdxList      [0];
   if ( LB_SendG_IDList         [0] != NULL )   delete [] LB_SendG_IDList         [0];
   if ( LB_SendG_SibList        [0] != NULL )   delete [] LB_SendG_SibList        [0];
   if ( LB_SendG_SibDiffList    [0] != NULL )   delete [] LB_SendG_SibDiffList    [0];
   if ( LB_SendG_LBIdxList      [0] != NULL )   delete [] LB_SendG_LBIdxList      [0];

   LB_RecvG_IDList_IdxTable[0] = new int  [NSend_Total_G];
   LB_RecvG_SibList        [0] = new int  [NSend_Total_G];
   LB_RecvG_SibDiffList    [0] = new int  [NSend_Total_G];
   LB_RecvG_LBIdxList      [0] = new long [NSend_Total_G];
   LB_RecvG_PCr1D          [0] = new ulong[NSend_Total_G];
   LB_RecvG_PCr1D_IdxTable [0] = new int  [NSend_Total_G];
   LB_SendG_IDList         [0] = new int  [NRecv_Total_G];
   LB_SendG_SibList        [0] = new int  [NRecv_Total_G];
   LB_SendG_SibDiffList    [0] = new int  [NRecv_Total_G];
   LB_SendG_LBIdxList      [0] = new long [NRecv_Total_G];

   for (int r=0; r<MPI_NRank; r++)
   {
//    5.2.2 assign pointers
      LB_RecvG_IDList_IdxTable[r] = LB_RecvG_IDList_IdxTable[0] + Send_Disp_G[r];
      LB_RecvG_SibList        [r] = LB_RecvG_SibList        [0] + Send_Disp_G[r];
      LB_RecvG_SibDiffList    [r] = LB_RecvG_SibDiffList    [0] + Send_Disp_G[r];
      LB_RecvG_LBIdxList      [r] = LB_RecvG_LBIdxList      [0] + Send_Disp_G[r];
      LB_RecvG_PCr1D          [r] = LB_RecvG_PCr1D          [0] + Send_Disp_G[r];
      LB_RecvG_PCr1D_IdxTable [r] = LB_RecvG_PCr1D_IdxTable [0] + Send_Disp_G[r];
      LB_SendG_IDList         [r] = LB_SendG_IDList         [0] + Recv_Disp_G[r];
      LB_SendG_SibList        [r] = LB_SendG_SibList        [0] + Recv_Disp_G[r];
      LB_SendG_SibDiffList    [r] = LB_SendG_SibDiffList    [0] + Recv_Disp_G[r];
      LB_SendG_LBIdxList      [r] = LB_SendG_LBIdxList      [0] + Recv_Disp_G[r];

//    5.2.3 get the sorted LBIdx list
      for (int t=0; t<LB_RecvG_NList[r]; t++)
         LB_RecvG_LBIdxList[r][t] = amr->patch[0][Lv][ LB_RecvG_IDList[r][t] ]->LB_Idx;

      Mis_Heapsort( LB_RecvG_NList[r], LB_RecvG_LBIdxList[r], LB_RecvG_IDList_IdxTable[r] );

//    5.2.4 get the sorted SibList and PaddedCr1D
      for (int t=0; t<LB_RecvG_NList[r]; t++)
      {
         TPID = LB_RecvG_IDList[r][ LB_RecvG_IDList_IdxTable[r][t] ];

         LB_RecvG_SibList[r][t] = LB_RecvG_SibList_Unsorted[r][ LB_RecvG_IDList_IdxTable[r][t] ];
         LB_RecvG_PCr1D  [r][t] = amr->patch[0][Lv][TPID]->PaddedCr1D;
      }

//    5.2.5 sort PaddedCr1D again and record the index table (for constructing the SibDiff lists later)
      Mis_Heapsort( LB_RecvG_NList[r], LB_RecvG_PCr1D[r], LB_RecvG_PCr1D_IdxTable[r] );

//    5.2.6 get the recv SibDiff lists
      if ( AfterRefine )
      {
         Match_G = new int [ LB_RecvG_NList[r] ];

//       match patches with the same LB_Idx
         Mis_Matching_int( Old_RecvG_NList[r], Old_RecvG_PCr1D[r], LB_RecvG_NList[r], LB_RecvG_PCr1D[r], Match_G );

//       remove sibling directions already exist in the previous MPI lists
         for (int t=0; t<LB_RecvG_NList[r]; t++)
         {
            NID = LB_RecvG_PCr1D_IdxTable[r][t];

            if ( Match_G[t] == -1 )
               LB_RecvG_SibDiffList[r][NID] = LB_RecvG_SibList[r][NID];

            else
            {
               OID = Old_RecvG_PCr1D_IdxTable[r][ Match_G[t] ];

//             add duplicate sibling directions to the old MPI lists (to minimize the amount of data to be resent)
               if ( Old_RecvG_SibList[r][OID] & (1<<26) )   Old_RecvG_SibList[r][OID] |= SibMask_Duplicate[26];
               else
               for (int s=0; s<18; s++)
               if ( Old_RecvG_SibList[r][OID] & (1<< s) )   Old_RecvG_SibList[r][OID] |= SibMask_Duplicate[ s];

//             record SibDiff List
               LB_RecvG_SibDiffList[r][NID] = LB_RecvG_SibList[r][NID] & ( ~Old_RecvG_SibList[r][OID] );
            }
         }

         delete [] Match_G;
      } // if ( AfterRefine )

//    recv PID should not repeat
#     ifdef GAMER_DEBUG
      int *TempIDList = new int [ LB_RecvG_NList[r] ];
      memcpy( TempIDList, LB_RecvG_IDList[r], LB_RecvG_NList[r]*sizeof(int) );

      Mis_Heapsort<int,int>( LB_RecvG_NList[r], TempIDList, NULL );

      for (int t=1; t<LB_RecvG_NList[r]; t++)
         if ( TempIDList[t] == TempIDList[t-1] )
            Aux_Error( ERROR_INFO, "Lv %d, Rank %d, repeated recv potential PID %d !!\n",
                       Lv, r, TempIDList[t] );

      delete [] TempIDList;
#     endif
   } // for (int r=0; r<MPI_NRank; r++)
#  endif // #ifdef GRAVITY



// 6. broadcast the recv list to all other ranks
// ============================================================================================================
// 6.1 hydro and magnetic field
   int  *SendBuf_SibList_H     = LB_RecvH_SibList    [0];
   int  *SendBuf_SibDiffList_H = LB_RecvH_SibDiffList[0];
   long *SendBuf_LBIdx_H       = LB_RecvH_LBIdxList  [0];
   int  *RecvBuf_SibList_H     = LB_SendH_SibList    [0];
   int  *RecvBuf_SibDiffList_H = LB_SendH_SibDiffList[0];
   long *RecvBuf_LBIdx_H       = LB_SendH_LBIdxList  [0];

   MPI_Alltoallv( SendBuf_SibList_H,     LB_RecvH_NList, Send_Disp_H, MPI_INT,
                  RecvBuf_SibList_H,     LB_SendH_NList, Recv_Disp_H, MPI_INT,  MPI_COMM_WORLD );

   if ( AfterRefine )
   MPI_Alltoallv( SendBuf_SibDiffList_H, LB_RecvH_NList, Send_Disp_H, MPI_INT,
                  RecvBuf_SibDiffList_H, LB_SendH_NList, Recv_Disp_H, MPI_INT,  MPI_COMM_WORLD );

   MPI_Alltoallv( SendBuf_LBIdx_H,       LB_RecvH_NList, Send_Disp_H, MPI_LONG,
                  RecvBuf_LBIdx_H,       LB_SendH_NList, Recv_Disp_H, MPI_LONG, MPI_COMM_WORLD );


// 6.2 potential
#  ifdef GRAVITY
   int  *SendBuf_SibList_G     = LB_RecvG_SibList    [0];
   int  *SendBuf_SibDiffList_G = LB_RecvG_SibDiffList[0];
   long *SendBuf_LBIdx_G       = LB_RecvG_LBIdxList  [0];
   int  *RecvBuf_SibList_G     = LB_SendG_SibList    [0];
   int  *RecvBuf_SibDiffList_G = LB_SendG_SibDiffList[0];
   long *RecvBuf_LBIdx_G       = LB_SendG_LBIdxList  [0];

   MPI_Alltoallv( SendBuf_SibList_G,     LB_RecvG_NList, Send_Disp_G, MPI_INT,
                  RecvBuf_SibList_G,     LB_SendG_NList, Recv_Disp_G, MPI_INT,  MPI_COMM_WORLD );

   if ( AfterRefine )
   MPI_Alltoallv( SendBuf_SibDiffList_G, LB_RecvG_NList, Send_Disp_G, MPI_INT,
                  RecvBuf_SibDiffList_G, LB_SendG_NList, Recv_Disp_G, MPI_INT,  MPI_COMM_WORLD );

   MPI_Alltoallv( SendBuf_LBIdx_G,       LB_RecvG_NList, Send_Disp_G, MPI_LONG,
                  RecvBuf_LBIdx_G,       LB_SendG_NList, Recv_Disp_G, MPI_LONG, MPI_COMM_WORLD );
#  endif // #ifdef GRAVITY



// 7. construct the send list
// ============================================================================================================
// 7.1 hydro and magnetic field
   for (int r=0; r<MPI_NRank; r++)
   {
//    7.1.1 allocate the matching array
      Match_H = new int [ LB_SendH_NList[r] ];

//    7.1.2 match the corresponding PID
      Mis_Matching_int( amr->NPatchComma[Lv][1], amr->LB->IdxList_Real[Lv], LB_SendH_NList[r],
                        LB_SendH_LBIdxList[r], Match_H );

//    7.1.3 check: all target patches must be found
#     ifdef GAMER_DEBUG
      for (int t=0; t<LB_SendH_NList[r]; t++)
         if ( Match_H[t] == -1 )
            Aux_Error( ERROR_INFO, "Lv %d, TRank %d, LB_Idx %ld found no matching patches (hydro) !!\n",
                       Lv, r, LB_SendH_LBIdxList[r][t] );
#     endif

//    7.1.4 store the patch indices to send data
      for (int t=0; t<LB_SendH_NList[r]; t++)
         LB_SendH_IDList[r][t] = amr->LB->IdxList_Real_IdxTable[Lv][ Match_H[t] ];

      delete [] Match_H;
   } // for (int r=0; r<MPI_NRank; r++)


// 7.2 potential
#  ifdef GRAVITY
   for (int r=0; r<MPI_NRank; r++)
   {
//    7.2.1 allocate the matching array
      Match_G = new int [ LB_SendG_NList[r] ];

//    7.2.2 match the corresponding PID
      Mis_Matching_int( amr->NPatchComma[Lv][1], amr->LB->IdxList_Real[Lv], LB_SendG_NList[r],
                        LB_SendG_LBIdxList[r], Match_G );

//    7.2.3 check: all target patches must be found
#     ifdef GAMER_DEBUG
      for (int t=0; t<LB_SendG_NList[r]; t++)
         if ( Match_G[t] == -1 )
            Aux_Error( ERROR_INFO, "Lv %d, TRank %d, LB_Idx %ld found no matching patches (potential) !!\n",
                       Lv, r, LB_SendG_LBIdxList[r][t] );
#     endif

//    7.2.4 store the patch indices to send data
      for (int t=0; t<LB_SendG_NList[r]; t++)
         LB_SendG_IDList[r][t] = amr->LB->IdxList_Real_IdxTable[Lv][ Match_G[t] ];

      delete [] Match_G;
   }
#  endif // #ifdef GRAVITY



// free memory
   for (int s=0; s<26; s++)
   {
      delete [] TSib_List   [s];
      delete [] RSib_List   [s];
      delete [] SibPID_Delta[s];
   }

   free( SibList_H );
   for (int r=0; r<MPI_NRank; r++)  free( LB_RecvH_SibList_Unsorted[r] );

   if ( Old_RecvH_SibList       [0] != NULL )    delete [] Old_RecvH_SibList       [0];
   if ( Old_RecvH_PCr1D         [0] != NULL )    delete [] Old_RecvH_PCr1D         [0];
   if ( Old_RecvH_PCr1D_IdxTable[0] != NULL )    delete [] Old_RecvH_PCr1D_IdxTable[0];

#  ifdef GRAVITY
   free( SibList_G );
   for (int r=0; r<MPI_NRank; r++)  free( LB_RecvG_SibList_Unsorted[r] );

   if ( Old_RecvG_SibList       [0] != NULL )    delete [] Old_RecvG_SibList       [0];
   if ( Old_RecvG_PCr1D         [0] != NULL )    delete [] Old_RecvG_PCr1D         [0];
   if ( Old_RecvG_PCr1D_IdxTable[0] != NULL )    delete [] Old_RecvG_PCr1D_IdxTable[0];
#  endif

} // FUNCTION : LB_RecordExchangeDataPatchID



//-------------------------------------------------------------------------------------------------------
// Function    :  SetSiblingMask
// Description :  Set the mask for determining which part of patch data to be sent
//
// Note        :  Element 0 ~ 25 : data in 26 sibling directions
//                Element 26     : entire patch data
//
// Parameter   :  SibMask_Check     : Mask for checking whether the target sibling direction has already
//                                    been included
//                SibMask_Clear     : Mask for clearing the duplicated sibling directions
//                SibMask_Duplicate : Mask for checking if there are duplicate target sibling directions
//-------------------------------------------------------------------------------------------------------
void SetSiblingMask( int SibMask_Check[], int SibMask_Clear[], int SibMask_Duplicate[] )
{

// check
   SibMask_Check[ 0] = (1<<26) | (1<< 0);
   SibMask_Check[ 1] = (1<<26) | (1<< 1);
   SibMask_Check[ 2] = (1<<26) | (1<< 2);
   SibMask_Check[ 3] = (1<<26) | (1<< 3);
   SibMask_Check[ 4] = (1<<26) | (1<< 4);
   SibMask_Check[ 5] = (1<<26) | (1<< 5);

   SibMask_Check[ 6] = (1<<26) | (1<< 6) | (1<<0) | (1<<2);
   SibMask_Check[ 7] = (1<<26) | (1<< 7) | (1<<1) | (1<<2);
   SibMask_Check[ 8] = (1<<26) | (1<< 8) | (1<<0) | (1<<3);
   SibMask_Check[ 9] = (1<<26) | (1<< 9) | (1<<1) | (1<<3);
   SibMask_Check[10] = (1<<26) | (1<<10) | (1<<2) | (1<<4);
   SibMask_Check[11] = (1<<26) | (1<<11) | (1<<3) | (1<<4);
   SibMask_Check[12] = (1<<26) | (1<<12) | (1<<2) | (1<<5);
   SibMask_Check[13] = (1<<26) | (1<<13) | (1<<3) | (1<<5);
   SibMask_Check[14] = (1<<26) | (1<<14) | (1<<0) | (1<<4);
   SibMask_Check[15] = (1<<26) | (1<<15) | (1<<0) | (1<<5);
   SibMask_Check[16] = (1<<26) | (1<<16) | (1<<1) | (1<<4);
   SibMask_Check[17] = (1<<26) | (1<<17) | (1<<1) | (1<<5);

   SibMask_Check[18] = (1<<26) | (1<<18) | (1<<0) | (1<<2) | (1<<4) | (1<<6) | (1<<10) | (1<<14);
   SibMask_Check[19] = (1<<26) | (1<<19) | (1<<1) | (1<<2) | (1<<4) | (1<<7) | (1<<10) | (1<<16);
   SibMask_Check[20] = (1<<26) | (1<<20) | (1<<0) | (1<<3) | (1<<4) | (1<<8) | (1<<11) | (1<<14);
   SibMask_Check[21] = (1<<26) | (1<<21) | (1<<1) | (1<<3) | (1<<4) | (1<<9) | (1<<11) | (1<<16);
   SibMask_Check[22] = (1<<26) | (1<<22) | (1<<0) | (1<<2) | (1<<5) | (1<<6) | (1<<12) | (1<<15);
   SibMask_Check[23] = (1<<26) | (1<<23) | (1<<1) | (1<<2) | (1<<5) | (1<<7) | (1<<12) | (1<<17);
   SibMask_Check[24] = (1<<26) | (1<<24) | (1<<0) | (1<<3) | (1<<5) | (1<<8) | (1<<13) | (1<<15);
   SibMask_Check[25] = (1<<26) | (1<<25) | (1<<1) | (1<<3) | (1<<5) | (1<<9) | (1<<13) | (1<<17);

   SibMask_Check[26] = (1<<26);


// duplicate
   SibMask_Duplicate[ 0] = (1<<18) | (1<<14) | (1<<20) | (1<< 6) | (1<< 8) | (1<<22) | (1<<15) | (1<<24);
   SibMask_Duplicate[ 1] = (1<<19) | (1<<16) | (1<<21) | (1<< 7) | (1<< 9) | (1<<23) | (1<<17) | (1<<25);
   SibMask_Duplicate[ 2] = (1<<18) | (1<<10) | (1<<19) | (1<< 6) | (1<< 7) | (1<<22) | (1<<12) | (1<<23);
   SibMask_Duplicate[ 3] = (1<<20) | (1<<11) | (1<<21) | (1<< 8) | (1<< 9) | (1<<24) | (1<<13) | (1<<25);
   SibMask_Duplicate[ 4] = (1<<18) | (1<<10) | (1<<19) | (1<<14) | (1<<16) | (1<<20) | (1<<11) | (1<<21);
   SibMask_Duplicate[ 5] = (1<<22) | (1<<12) | (1<<23) | (1<<15) | (1<<17) | (1<<24) | (1<<13) | (1<<25);

   SibMask_Duplicate[ 6] = (1<<18) | (1<<22);
   SibMask_Duplicate[ 7] = (1<<19) | (1<<23);
   SibMask_Duplicate[ 8] = (1<<20) | (1<<24);
   SibMask_Duplicate[ 9] = (1<<21) | (1<<25);
   SibMask_Duplicate[10] = (1<<18) | (1<<19);
   SibMask_Duplicate[11] = (1<<20) | (1<<21);
   SibMask_Duplicate[12] = (1<<22) | (1<<23);
   SibMask_Duplicate[13] = (1<<24) | (1<<25);
   SibMask_Duplicate[14] = (1<<18) | (1<<20);
   SibMask_Duplicate[15] = (1<<22) | (1<<24);
   SibMask_Duplicate[16] = (1<<19) | (1<<21);
   SibMask_Duplicate[17] = (1<<23) | (1<<25);

   SibMask_Duplicate[18] = 0;
   SibMask_Duplicate[19] = 0;
   SibMask_Duplicate[20] = 0;
   SibMask_Duplicate[21] = 0;
   SibMask_Duplicate[22] = 0;
   SibMask_Duplicate[23] = 0;
   SibMask_Duplicate[24] = 0;
   SibMask_Duplicate[25] = 0;

   SibMask_Duplicate[26] = ~(1<<26);


// clear
   for (int s=0; s<27; s++)   SibMask_Clear[s] = ~SibMask_Duplicate[s];

} // FUNCTION : SetSiblingMask



//-------------------------------------------------------------------------------------------------------
// Function    :  SetReceiveSibling
// Description :  Set the sibling directions for receiving ghost-zone data at the coarse-grid level
//
// Note        :  1. RSib_List[] needs to be deallocated manually
//                2. The order of sibling indices recorded in RSib_List[] must be defined consistently with those
//                   defined in SetTargetSibling()
//
// Parameter   :  RSib_List : Target sibling indices along different sibling directions
//-------------------------------------------------------------------------------------------------------
void SetReceiveSibling( int* RSib_List[] )
{

   int NRSib[26];

   for (int t= 0; t< 6; t++)  NRSib[t] = 18;
   for (int t= 6; t<18; t++)  NRSib[t] = 12;
   for (int t=18; t<26; t++)  NRSib[t] =  8;

   for (int s=0; s<26; s++)   RSib_List[s] = new int [ NRSib[s] ];

   RSib_List[ 0][ 0] =  1;
   RSib_List[ 0][ 1] = 25;
   RSib_List[ 0][ 2] = 24;
   RSib_List[ 0][ 3] = 17;
   RSib_List[ 0][ 4] = 15;
   RSib_List[ 0][ 5] = 23;
   RSib_List[ 0][ 6] = 22;
   RSib_List[ 0][ 7] =  9;
   RSib_List[ 0][ 8] =  8;
   RSib_List[ 0][ 9] =  0;
   RSib_List[ 0][10] =  7;
   RSib_List[ 0][11] =  6;
   RSib_List[ 0][12] = 21;
   RSib_List[ 0][13] = 20;
   RSib_List[ 0][14] = 16;
   RSib_List[ 0][15] = 14;
   RSib_List[ 0][16] = 19;
   RSib_List[ 0][17] = 18;

   RSib_List[ 1][ 0] =  0;
   RSib_List[ 1][ 1] = 25;
   RSib_List[ 1][ 2] = 24;
   RSib_List[ 1][ 3] = 17;
   RSib_List[ 1][ 4] = 15;
   RSib_List[ 1][ 5] = 23;
   RSib_List[ 1][ 6] = 22;
   RSib_List[ 1][ 7] =  9;
   RSib_List[ 1][ 8] =  8;
   RSib_List[ 1][ 9] =  1;
   RSib_List[ 1][10] =  7;
   RSib_List[ 1][11] =  6;
   RSib_List[ 1][12] = 21;
   RSib_List[ 1][13] = 20;
   RSib_List[ 1][14] = 16;
   RSib_List[ 1][15] = 14;
   RSib_List[ 1][16] = 19;
   RSib_List[ 1][17] = 18;

   RSib_List[ 2][ 0] =  3;
   RSib_List[ 2][ 1] = 25;
   RSib_List[ 2][ 2] = 13;
   RSib_List[ 2][ 3] = 24;
   RSib_List[ 2][ 4] = 23;
   RSib_List[ 2][ 5] = 12;
   RSib_List[ 2][ 6] = 22;
   RSib_List[ 2][ 7] =  9;
   RSib_List[ 2][ 8] =  8;
   RSib_List[ 2][ 9] =  7;
   RSib_List[ 2][10] =  2;
   RSib_List[ 2][11] =  6;
   RSib_List[ 2][12] = 21;
   RSib_List[ 2][13] = 11;
   RSib_List[ 2][14] = 20;
   RSib_List[ 2][15] = 19;
   RSib_List[ 2][16] = 10;
   RSib_List[ 2][17] = 18;

   RSib_List[ 3][ 0] =  2;
   RSib_List[ 3][ 1] = 25;
   RSib_List[ 3][ 2] = 13;
   RSib_List[ 3][ 3] = 24;
   RSib_List[ 3][ 4] = 23;
   RSib_List[ 3][ 5] = 12;
   RSib_List[ 3][ 6] = 22;
   RSib_List[ 3][ 7] =  9;
   RSib_List[ 3][ 8] =  3;
   RSib_List[ 3][ 9] =  8;
   RSib_List[ 3][10] =  7;
   RSib_List[ 3][11] =  6;
   RSib_List[ 3][12] = 21;
   RSib_List[ 3][13] = 11;
   RSib_List[ 3][14] = 20;
   RSib_List[ 3][15] = 19;
   RSib_List[ 3][16] = 10;
   RSib_List[ 3][17] = 18;

   RSib_List[ 4][ 0] =  5;
   RSib_List[ 4][ 1] = 25;
   RSib_List[ 4][ 2] = 13;
   RSib_List[ 4][ 3] = 24;
   RSib_List[ 4][ 4] = 17;
   RSib_List[ 4][ 5] = 15;
   RSib_List[ 4][ 6] = 23;
   RSib_List[ 4][ 7] = 12;
   RSib_List[ 4][ 8] = 22;
   RSib_List[ 4][ 9] = 21;
   RSib_List[ 4][10] = 11;
   RSib_List[ 4][11] = 20;
   RSib_List[ 4][12] = 16;
   RSib_List[ 4][13] =  4;
   RSib_List[ 4][14] = 14;
   RSib_List[ 4][15] = 19;
   RSib_List[ 4][16] = 10;
   RSib_List[ 4][17] = 18;

   RSib_List[ 5][ 0] =  4;
   RSib_List[ 5][ 1] = 25;
   RSib_List[ 5][ 2] = 13;
   RSib_List[ 5][ 3] = 24;
   RSib_List[ 5][ 4] = 17;
   RSib_List[ 5][ 5] =  5;
   RSib_List[ 5][ 6] = 15;
   RSib_List[ 5][ 7] = 23;
   RSib_List[ 5][ 8] = 12;
   RSib_List[ 5][ 9] = 22;
   RSib_List[ 5][10] = 21;
   RSib_List[ 5][11] = 11;
   RSib_List[ 5][12] = 20;
   RSib_List[ 5][13] = 16;
   RSib_List[ 5][14] = 14;
   RSib_List[ 5][15] = 19;
   RSib_List[ 5][16] = 10;
   RSib_List[ 5][17] = 18;

   RSib_List[ 6][ 0] =  9;
   RSib_List[ 6][ 1] = 25;
   RSib_List[ 6][ 2] = 24;
   RSib_List[ 6][ 3] = 23;
   RSib_List[ 6][ 4] = 22;
   RSib_List[ 6][ 5] =  8;
   RSib_List[ 6][ 6] =  7;
   RSib_List[ 6][ 7] =  6;
   RSib_List[ 6][ 8] = 21;
   RSib_List[ 6][ 9] = 20;
   RSib_List[ 6][10] = 19;
   RSib_List[ 6][11] = 18;

   RSib_List[ 7][ 0] =  8;
   RSib_List[ 7][ 1] = 25;
   RSib_List[ 7][ 2] = 24;
   RSib_List[ 7][ 3] = 23;
   RSib_List[ 7][ 4] = 22;
   RSib_List[ 7][ 5] =  9;
   RSib_List[ 7][ 6] =  7;
   RSib_List[ 7][ 7] =  6;
   RSib_List[ 7][ 8] = 21;
   RSib_List[ 7][ 9] = 20;
   RSib_List[ 7][10] = 19;
   RSib_List[ 7][11] = 18;

   RSib_List[ 8][ 0] =  7;
   RSib_List[ 8][ 1] = 25;
   RSib_List[ 8][ 2] = 24;
   RSib_List[ 8][ 3] = 23;
   RSib_List[ 8][ 4] = 22;
   RSib_List[ 8][ 5] =  9;
   RSib_List[ 8][ 6] =  8;
   RSib_List[ 8][ 7] =  6;
   RSib_List[ 8][ 8] = 21;
   RSib_List[ 8][ 9] = 20;
   RSib_List[ 8][10] = 19;
   RSib_List[ 8][11] = 18;

   RSib_List[ 9][ 0] =  6;
   RSib_List[ 9][ 1] = 25;
   RSib_List[ 9][ 2] = 24;
   RSib_List[ 9][ 3] = 23;
   RSib_List[ 9][ 4] = 22;
   RSib_List[ 9][ 5] =  9;
   RSib_List[ 9][ 6] =  8;
   RSib_List[ 9][ 7] =  7;
   RSib_List[ 9][ 8] = 21;
   RSib_List[ 9][ 9] = 20;
   RSib_List[ 9][10] = 19;
   RSib_List[ 9][11] = 18;

   RSib_List[10][ 0] = 13;
   RSib_List[10][ 1] = 25;
   RSib_List[10][ 2] = 24;
   RSib_List[10][ 3] = 23;
   RSib_List[10][ 4] = 12;
   RSib_List[10][ 5] = 22;
   RSib_List[10][ 6] = 21;
   RSib_List[10][ 7] = 11;
   RSib_List[10][ 8] = 20;
   RSib_List[10][ 9] = 19;
   RSib_List[10][10] = 10;
   RSib_List[10][11] = 18;

   RSib_List[11][ 0] = 12;
   RSib_List[11][ 1] = 25;
   RSib_List[11][ 2] = 13;
   RSib_List[11][ 3] = 24;
   RSib_List[11][ 4] = 23;
   RSib_List[11][ 5] = 22;
   RSib_List[11][ 6] = 21;
   RSib_List[11][ 7] = 11;
   RSib_List[11][ 8] = 20;
   RSib_List[11][ 9] = 19;
   RSib_List[11][10] = 10;
   RSib_List[11][11] = 18;

   RSib_List[12][ 0] = 11;
   RSib_List[12][ 1] = 25;
   RSib_List[12][ 2] = 13;
   RSib_List[12][ 3] = 24;
   RSib_List[12][ 4] = 23;
   RSib_List[12][ 5] = 12;
   RSib_List[12][ 6] = 22;
   RSib_List[12][ 7] = 21;
   RSib_List[12][ 8] = 20;
   RSib_List[12][ 9] = 19;
   RSib_List[12][10] = 10;
   RSib_List[12][11] = 18;

   RSib_List[13][ 0] = 10;
   RSib_List[13][ 1] = 25;
   RSib_List[13][ 2] = 13;
   RSib_List[13][ 3] = 24;
   RSib_List[13][ 4] = 23;
   RSib_List[13][ 5] = 12;
   RSib_List[13][ 6] = 22;
   RSib_List[13][ 7] = 21;
   RSib_List[13][ 8] = 11;
   RSib_List[13][ 9] = 20;
   RSib_List[13][10] = 19;
   RSib_List[13][11] = 18;

   RSib_List[14][ 0] = 17;
   RSib_List[14][ 1] = 25;
   RSib_List[14][ 2] = 24;
   RSib_List[14][ 3] = 15;
   RSib_List[14][ 4] = 23;
   RSib_List[14][ 5] = 22;
   RSib_List[14][ 6] = 21;
   RSib_List[14][ 7] = 20;
   RSib_List[14][ 8] = 16;
   RSib_List[14][ 9] = 14;
   RSib_List[14][10] = 19;
   RSib_List[14][11] = 18;

   RSib_List[15][ 0] = 16;
   RSib_List[15][ 1] = 25;
   RSib_List[15][ 2] = 24;
   RSib_List[15][ 3] = 17;
   RSib_List[15][ 4] = 15;
   RSib_List[15][ 5] = 23;
   RSib_List[15][ 6] = 22;
   RSib_List[15][ 7] = 21;
   RSib_List[15][ 8] = 20;
   RSib_List[15][ 9] = 14;
   RSib_List[15][10] = 19;
   RSib_List[15][11] = 18;

   RSib_List[16][ 0] = 15;
   RSib_List[16][ 1] = 25;
   RSib_List[16][ 2] = 24;
   RSib_List[16][ 3] = 17;
   RSib_List[16][ 4] = 23;
   RSib_List[16][ 5] = 22;
   RSib_List[16][ 6] = 21;
   RSib_List[16][ 7] = 20;
   RSib_List[16][ 8] = 16;
   RSib_List[16][ 9] = 14;
   RSib_List[16][10] = 19;
   RSib_List[16][11] = 18;

   RSib_List[17][ 0] = 14;
   RSib_List[17][ 1] = 25;
   RSib_List[17][ 2] = 24;
   RSib_List[17][ 3] = 17;
   RSib_List[17][ 4] = 15;
   RSib_List[17][ 5] = 23;
   RSib_List[17][ 6] = 22;
   RSib_List[17][ 7] = 21;
   RSib_List[17][ 8] = 20;
   RSib_List[17][ 9] = 16;
   RSib_List[17][10] = 19;
   RSib_List[17][11] = 18;

   RSib_List[18][ 0] = 25;
   RSib_List[18][ 1] = 24;
   RSib_List[18][ 2] = 23;
   RSib_List[18][ 3] = 22;
   RSib_List[18][ 4] = 21;
   RSib_List[18][ 5] = 20;
   RSib_List[18][ 6] = 19;
   RSib_List[18][ 7] = 18;

   RSib_List[19][ 0] = 24;
   RSib_List[19][ 1] = 25;
   RSib_List[19][ 2] = 23;
   RSib_List[19][ 3] = 22;
   RSib_List[19][ 4] = 21;
   RSib_List[19][ 5] = 20;
   RSib_List[19][ 6] = 19;
   RSib_List[19][ 7] = 18;

   RSib_List[20][ 0] = 23;
   RSib_List[20][ 1] = 25;
   RSib_List[20][ 2] = 24;
   RSib_List[20][ 3] = 22;
   RSib_List[20][ 4] = 21;
   RSib_List[20][ 5] = 20;
   RSib_List[20][ 6] = 19;
   RSib_List[20][ 7] = 18;

   RSib_List[21][ 0] = 22;
   RSib_List[21][ 1] = 25;
   RSib_List[21][ 2] = 24;
   RSib_List[21][ 3] = 23;
   RSib_List[21][ 4] = 21;
   RSib_List[21][ 5] = 20;
   RSib_List[21][ 6] = 19;
   RSib_List[21][ 7] = 18;

   RSib_List[22][ 0] = 21;
   RSib_List[22][ 1] = 25;
   RSib_List[22][ 2] = 24;
   RSib_List[22][ 3] = 23;
   RSib_List[22][ 4] = 22;
   RSib_List[22][ 5] = 20;
   RSib_List[22][ 6] = 19;
   RSib_List[22][ 7] = 18;

   RSib_List[23][ 0] = 20;
   RSib_List[23][ 1] = 25;
   RSib_List[23][ 2] = 24;
   RSib_List[23][ 3] = 23;
   RSib_List[23][ 4] = 22;
   RSib_List[23][ 5] = 21;
   RSib_List[23][ 6] = 19;
   RSib_List[23][ 7] = 18;

   RSib_List[24][ 0] = 19;
   RSib_List[24][ 1] = 25;
   RSib_List[24][ 2] = 24;
   RSib_List[24][ 3] = 23;
   RSib_List[24][ 4] = 22;
   RSib_List[24][ 5] = 21;
   RSib_List[24][ 6] = 20;
   RSib_List[24][ 7] = 18;

   RSib_List[25][ 0] = 18;
   RSib_List[25][ 1] = 25;
   RSib_List[25][ 2] = 24;
   RSib_List[25][ 3] = 23;
   RSib_List[25][ 4] = 22;
   RSib_List[25][ 5] = 21;
   RSib_List[25][ 6] = 20;
   RSib_List[25][ 7] = 19;

} // FUNCTION : SetReceiveSibling



//-------------------------------------------------------------------------------------------------------
// Function    :  SetTargetSibling
// Description :  Set the target sibling directions for looping over all father-sibling patches
//
// Note        :  1. TSib needs to be deallocated manually
//                2. The order of sibling indices recorded in TSib must be defined consistently with those
//                   defined in "SetReceiveSibling"
//
// Parameter   :  NTSib : Number of target sibling patches along different sibling directions
//                TSib  : Target sibling indices along different sibling directions
//-------------------------------------------------------------------------------------------------------
void SetTargetSibling( int NTSib[], int* TSib[] )
{

   for (int t= 0; t< 6; t++)  NTSib[t] = 17;
   for (int t= 6; t<18; t++)  NTSib[t] = 11;
   for (int t=18; t<26; t++)  NTSib[t] =  7;

   for (int s=0; s<26; s++)   TSib[s] = new int [ NTSib[s] ];

   TSib[ 0][ 0] = 10;
   TSib[ 0][ 1] = 19;
   TSib[ 0][ 2] =  4;
   TSib[ 0][ 3] = 16;
   TSib[ 0][ 4] = 11;
   TSib[ 0][ 5] = 21;
   TSib[ 0][ 6] =  2;
   TSib[ 0][ 7] =  7;
   TSib[ 0][ 8] =  1;
   TSib[ 0][ 9] =  3;
   TSib[ 0][10] =  9;
   TSib[ 0][11] = 12;
   TSib[ 0][12] = 23;
   TSib[ 0][13] =  5;
   TSib[ 0][14] = 17;
   TSib[ 0][15] = 13;
   TSib[ 0][16] = 25;

   TSib[ 1][ 0] = 18;
   TSib[ 1][ 1] = 10;
   TSib[ 1][ 2] = 14;
   TSib[ 1][ 3] =  4;
   TSib[ 1][ 4] = 20;
   TSib[ 1][ 5] = 11;
   TSib[ 1][ 6] =  6;
   TSib[ 1][ 7] =  2;
   TSib[ 1][ 8] =  0;
   TSib[ 1][ 9] =  8;
   TSib[ 1][10] =  3;
   TSib[ 1][11] = 22;
   TSib[ 1][12] = 12;
   TSib[ 1][13] = 15;
   TSib[ 1][14] =  5;
   TSib[ 1][15] = 24;
   TSib[ 1][16] = 13;

   TSib[ 2][ 0] = 14;
   TSib[ 2][ 1] =  4;
   TSib[ 2][ 2] = 16;
   TSib[ 2][ 3] = 20;
   TSib[ 2][ 4] = 11;
   TSib[ 2][ 5] = 21;
   TSib[ 2][ 6] =  0;
   TSib[ 2][ 7] =  1;
   TSib[ 2][ 8] =  8;
   TSib[ 2][ 9] =  3;
   TSib[ 2][10] =  9;
   TSib[ 2][11] = 15;
   TSib[ 2][12] =  5;
   TSib[ 2][13] = 17;
   TSib[ 2][14] = 24;
   TSib[ 2][15] = 13;
   TSib[ 2][16] = 25;

   TSib[ 3][ 0] = 18;
   TSib[ 3][ 1] = 10;
   TSib[ 3][ 2] = 19;
   TSib[ 3][ 3] = 14;
   TSib[ 3][ 4] =  4;
   TSib[ 3][ 5] = 16;
   TSib[ 3][ 6] =  6;
   TSib[ 3][ 7] =  2;
   TSib[ 3][ 8] =  7;
   TSib[ 3][ 9] =  0;
   TSib[ 3][10] =  1;
   TSib[ 3][11] = 22;
   TSib[ 3][12] = 12;
   TSib[ 3][13] = 23;
   TSib[ 3][14] = 15;
   TSib[ 3][15] =  5;
   TSib[ 3][16] = 17;

   TSib[ 4][ 0] =  6;
   TSib[ 4][ 1] =  2;
   TSib[ 4][ 2] =  7;
   TSib[ 4][ 3] =  0;
   TSib[ 4][ 4] =  1;
   TSib[ 4][ 5] =  8;
   TSib[ 4][ 6] =  3;
   TSib[ 4][ 7] =  9;
   TSib[ 4][ 8] = 22;
   TSib[ 4][ 9] = 12;
   TSib[ 4][10] = 23;
   TSib[ 4][11] = 15;
   TSib[ 4][12] =  5;
   TSib[ 4][13] = 17;
   TSib[ 4][14] = 24;
   TSib[ 4][15] = 13;
   TSib[ 4][16] = 25;

   TSib[ 5][ 0] = 18;
   TSib[ 5][ 1] = 10;
   TSib[ 5][ 2] = 19;
   TSib[ 5][ 3] = 14;
   TSib[ 5][ 4] =  4;
   TSib[ 5][ 5] = 16;
   TSib[ 5][ 6] = 20;
   TSib[ 5][ 7] = 11;
   TSib[ 5][ 8] = 21;
   TSib[ 5][ 9] =  6;
   TSib[ 5][10] =  2;
   TSib[ 5][11] =  7;
   TSib[ 5][12] =  0;
   TSib[ 5][13] =  1;
   TSib[ 5][14] =  8;
   TSib[ 5][15] =  3;
   TSib[ 5][16] =  9;

   TSib[ 6][ 0] =  4;
   TSib[ 6][ 1] = 16;
   TSib[ 6][ 2] = 11;
   TSib[ 6][ 3] = 21;
   TSib[ 6][ 4] =  1;
   TSib[ 6][ 5] =  3;
   TSib[ 6][ 6] =  9;
   TSib[ 6][ 7] =  5;
   TSib[ 6][ 8] = 17;
   TSib[ 6][ 9] = 13;
   TSib[ 6][10] = 25;

   TSib[ 7][ 0] = 14;
   TSib[ 7][ 1] =  4;
   TSib[ 7][ 2] = 20;
   TSib[ 7][ 3] = 11;
   TSib[ 7][ 4] =  0;
   TSib[ 7][ 5] =  8;
   TSib[ 7][ 6] =  3;
   TSib[ 7][ 7] = 15;
   TSib[ 7][ 8] =  5;
   TSib[ 7][ 9] = 24;
   TSib[ 7][10] = 13;

   TSib[ 8][ 0] = 10;
   TSib[ 8][ 1] = 19;
   TSib[ 8][ 2] =  4;
   TSib[ 8][ 3] = 16;
   TSib[ 8][ 4] =  2;
   TSib[ 8][ 5] =  7;
   TSib[ 8][ 6] =  1;
   TSib[ 8][ 7] = 12;
   TSib[ 8][ 8] = 23;
   TSib[ 8][ 9] =  5;
   TSib[ 8][10] = 17;

   TSib[ 9][ 0] = 18;
   TSib[ 9][ 1] = 10;
   TSib[ 9][ 2] = 14;
   TSib[ 9][ 3] =  4;
   TSib[ 9][ 4] =  6;
   TSib[ 9][ 5] =  2;
   TSib[ 9][ 6] =  0;
   TSib[ 9][ 7] = 22;
   TSib[ 9][ 8] = 12;
   TSib[ 9][ 9] = 15;
   TSib[ 9][10] =  5;

   TSib[10][ 0] =  0;
   TSib[10][ 1] =  1;
   TSib[10][ 2] =  8;
   TSib[10][ 3] =  3;
   TSib[10][ 4] =  9;
   TSib[10][ 5] = 15;
   TSib[10][ 6] =  5;
   TSib[10][ 7] = 17;
   TSib[10][ 8] = 24;
   TSib[10][ 9] = 13;
   TSib[10][10] = 25;

   TSib[11][ 0] =  6;
   TSib[11][ 1] =  2;
   TSib[11][ 2] =  7;
   TSib[11][ 3] =  0;
   TSib[11][ 4] =  1;
   TSib[11][ 5] = 22;
   TSib[11][ 6] = 12;
   TSib[11][ 7] = 23;
   TSib[11][ 8] = 15;
   TSib[11][ 9] =  5;
   TSib[11][10] = 17;

   TSib[12][ 0] = 14;
   TSib[12][ 1] =  4;
   TSib[12][ 2] = 16;
   TSib[12][ 3] = 20;
   TSib[12][ 4] = 11;
   TSib[12][ 5] = 21;
   TSib[12][ 6] =  0;
   TSib[12][ 7] =  1;
   TSib[12][ 8] =  8;
   TSib[12][ 9] =  3;
   TSib[12][10] =  9;

   TSib[13][ 0] = 18;
   TSib[13][ 1] = 10;
   TSib[13][ 2] = 19;
   TSib[13][ 3] = 14;
   TSib[13][ 4] =  4;
   TSib[13][ 5] = 16;
   TSib[13][ 6] =  6;
   TSib[13][ 7] =  2;
   TSib[13][ 8] =  7;
   TSib[13][ 9] =  0;
   TSib[13][10] =  1;

   TSib[14][ 0] =  2;
   TSib[14][ 1] =  7;
   TSib[14][ 2] =  1;
   TSib[14][ 3] =  3;
   TSib[14][ 4] =  9;
   TSib[14][ 5] = 12;
   TSib[14][ 6] = 23;
   TSib[14][ 7] =  5;
   TSib[14][ 8] = 17;
   TSib[14][ 9] = 13;
   TSib[14][10] = 25;

   TSib[15][ 0] = 10;
   TSib[15][ 1] = 19;
   TSib[15][ 2] =  4;
   TSib[15][ 3] = 16;
   TSib[15][ 4] = 11;
   TSib[15][ 5] = 21;
   TSib[15][ 6] =  2;
   TSib[15][ 7] =  7;
   TSib[15][ 8] =  1;
   TSib[15][ 9] =  3;
   TSib[15][10] =  9;

   TSib[16][ 0] =  6;
   TSib[16][ 1] =  2;
   TSib[16][ 2] =  0;
   TSib[16][ 3] =  8;
   TSib[16][ 4] =  3;
   TSib[16][ 5] = 22;
   TSib[16][ 6] = 12;
   TSib[16][ 7] = 15;
   TSib[16][ 8] =  5;
   TSib[16][ 9] = 24;
   TSib[16][10] = 13;

   TSib[17][ 0] = 18;
   TSib[17][ 1] = 10;
   TSib[17][ 2] = 14;
   TSib[17][ 3] =  4;
   TSib[17][ 4] = 20;
   TSib[17][ 5] = 11;
   TSib[17][ 6] =  6;
   TSib[17][ 7] =  2;
   TSib[17][ 8] =  0;
   TSib[17][ 9] =  8;
   TSib[17][10] =  3;

   TSib[18][ 0] =  1;
   TSib[18][ 1] =  3;
   TSib[18][ 2] =  9;
   TSib[18][ 3] =  5;
   TSib[18][ 4] = 17;
   TSib[18][ 5] = 13;
   TSib[18][ 6] = 25;

   TSib[19][ 0] =  0;
   TSib[19][ 1] =  8;
   TSib[19][ 2] =  3;
   TSib[19][ 3] = 15;
   TSib[19][ 4] =  5;
   TSib[19][ 5] = 24;
   TSib[19][ 6] = 13;

   TSib[20][ 0] =  2;
   TSib[20][ 1] =  7;
   TSib[20][ 2] =  1;
   TSib[20][ 3] = 12;
   TSib[20][ 4] = 23;
   TSib[20][ 5] =  5;
   TSib[20][ 6] = 17;

   TSib[21][ 0] =  6;
   TSib[21][ 1] =  2;
   TSib[21][ 2] =  0;
   TSib[21][ 3] = 22;
   TSib[21][ 4] = 12;
   TSib[21][ 5] = 15;
   TSib[21][ 6] =  5;

   TSib[22][ 0] =  4;
   TSib[22][ 1] = 16;
   TSib[22][ 2] = 11;
   TSib[22][ 3] = 21;
   TSib[22][ 4] =  1;
   TSib[22][ 5] =  3;
   TSib[22][ 6] =  9;

   TSib[23][ 0] = 14;
   TSib[23][ 1] =  4;
   TSib[23][ 2] = 20;
   TSib[23][ 3] = 11;
   TSib[23][ 4] =  0;
   TSib[23][ 5] =  8;
   TSib[23][ 6] =  3;

   TSib[24][ 0] = 10;
   TSib[24][ 1] = 19;
   TSib[24][ 2] =  4;
   TSib[24][ 3] = 16;
   TSib[24][ 4] =  2;
   TSib[24][ 5] =  7;
   TSib[24][ 6] =  1;

   TSib[25][ 0] = 18;
   TSib[25][ 1] = 10;
   TSib[25][ 2] = 14;
   TSib[25][ 3] =  4;
   TSib[25][ 4] =  6;
   TSib[25][ 5] =  2;
   TSib[25][ 6] =  0;

} // FUNCTION : SetTargetSibling



#endif // #ifdef LOAD_BALANCE
