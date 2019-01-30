#include "GAMER.h"

#if ( defined PARTICLE  &&  defined LOAD_BALANCE )



// defined in LB_RecordExchangeDataPatchID.cpp
extern void SetTargetLocalID( int NTLocalID[], int *TLocalID[] );
extern void SetTargetSibPID0( const int lv, const int PID0, int SibPID0_List[] );




//-------------------------------------------------------------------------------------------------------
// Function    :  Par_LB_RecordExchangeParticlePatchID
// Description :  Record the patch indices for exchanging particles between different ranks
//
// Note        :  1. R2B : send particles from real to buffer patches
//                   1-1. Include buffer patches at both MainLv and MainLv-1 surrounding real patches at MainLv
//                   1-2. Construct R2B_Real/Buff_NPatchTotal[lv][0/1], R2B_Real/Buff_NPatchEachRank[lv][0/1],
//                        R2B_Real/Buff_PIDList[lv][0/1]
//                        --> They all have the dimension [NLEVEL][2], where [MainLv][0/1] is for receiving particles
//                            at MainLv/MainLv-1 buffer patches adjacent to real patches at MainLv
//                        --> Mainly for the Poisson Solver at MainLv (i.e., for calculating the total density field at MainLv)
//                        --> More specific,
//                            [MainLv][0] is for receiving the particles of sibling-buffer patches at MainLv
//                            adjacent to real patches at MainLv
//                            [MainLv][1] is for receiving the particles of father-sibling-buffer patches at MainLv-1
//                            adjacent to real patches at MainLv
//
//                2. B2R : send particles from buffer to real patches
//                   2-1. Include buffer patches at both MainLv and MainLv-1 surrounding real patches at MainLv
//                   2-2. Construct B2R_Send/Buff_NPatchTotal[lv][0/1], B2R_Send/Buff_NPatchEachRank[lv][0/1],
//                        B2R_Send/Buff_PIDList[lv][0/1]
//                        --> Mainly for sending particles temporarily stored in the buffer patches just after updating
//                            their position to their corresponding real patches
//                        --> Similar to the R2B case. But instead of considering **all** real patches, it only
//                            considers **leaf** real patches (since only these patches can have particles)
//                        --> Therefore, B2R list is a subset of the R2B list
//
//                3. F2S : send particles from fathers (at MainLv-1) to sons (at MainLv)
//                   3-1. Exact procedure is to send particles from real father patches at MainLv-1 to the
//                        corresponding father-buffer patches first, and then call Par_PassParticle2Son_MultiPatch() to transfer
//                        particles from father-buffer patches to their real son patches at MainLv in the same rank
//                        --> This function only records the real father patches at MainLv-1 (those to send particles)
//                            and the corresponding father-buffer patches at MainLv-1 (those to receive particles)
//
// Parameter   :  MainLv : Main target refinement level
//
// Return      :  R2B_Real/Buff_NPatchTotal[MainLv][0/1], R2B_Real/Buff_NPatchEachRank[MainLv][0/1], R2B_Real/Buff_PIDList[MainLv][0/1]
//                B2R_Real/Buff_NPatchTotal[MainLv][0/1], B2R_Real/Buff_NPatchEachRank[MainLv][0/1], B2R_Real/Buff_PIDList[MainLv][0/1]
//                F2S_Send/Recv_NPatchTotal[MainLv-1], F2S_Send/Recv_NPatchEachRank[MainLv-1], F2S_Send/Recv_PIDList[MainLv-1]
//-------------------------------------------------------------------------------------------------------
void Par_LB_RecordExchangeParticlePatchID( const int MainLv )
{

// nothing to do for levels above MAX_LEVEL
   if ( MainLv > MAX_LEVEL )  return;


   const int FaLv         = MainLv - 1;
   const int RelatedLv[2] = { MainLv, FaLv };
   const int NLv          = ( MainLv > 0 ) ? 2 : 1;

   int lv, NReal[2], NBuff[2], MemUnit_R2B[2], MemUnit_B2R[2], MemUnit_F2S, MemSize_R2B[2], MemSize_B2R[2], MemSize_F2S;
   int FaPID, FaSibPID, SibPID, SibPID0, SibPID0_List[26], NTLocalID[26], *TLocalID[26], Buff_NPatchTotal_Dup;


// 1. initialize arrays
   for (int t=0; t<NLv; t++)
   {
      lv             = RelatedLv[t];
      NReal      [t] = amr->NPatchComma[lv][1];
      NBuff      [t] = amr->NPatchComma[lv][3] - amr->NPatchComma[lv][1];
      MemUnit_R2B[t] = MAX( 1, NBuff[t]/4  );   // set arbitrarily (but must > 0)
      MemUnit_B2R[t] = MAX( 1, NBuff[t]/16 );   // set arbitrarily (but must > 0)
      MemSize_R2B[t] = MemUnit_R2B[t];
      MemSize_B2R[t] = MemUnit_B2R[t];

      if ( amr->Par->R2B_Buff_PIDList[MainLv][t] != NULL )  free( amr->Par->R2B_Buff_PIDList[MainLv][t] );

      amr->Par->R2B_Buff_NPatchTotal[MainLv][t] = 0;
      amr->Par->R2B_Buff_PIDList    [MainLv][t] = (int*)malloc( MemSize_R2B[t]*sizeof(int) );

      if ( amr->Par->B2R_Buff_PIDList[MainLv][t] != NULL )  free( amr->Par->B2R_Buff_PIDList[MainLv][t] );

      amr->Par->B2R_Buff_NPatchTotal[MainLv][t] = 0;
      amr->Par->B2R_Buff_PIDList    [MainLv][t] = (int*)malloc( MemSize_B2R[t]*sizeof(int) );

      if ( lv == FaLv ) {
      MemUnit_F2S = MAX( 1, NReal[t]/8 );    // set arbitrarily (but must > 0)
      MemSize_F2S = MemUnit_F2S;

      if ( amr->Par->F2S_Send_PIDList[FaLv] != NULL )    free( amr->Par->F2S_Send_PIDList[FaLv] );

      amr->Par->F2S_Send_NPatchTotal[FaLv] = 0;
      amr->Par->F2S_Send_PIDList    [FaLv] = (int*)malloc( MemSize_F2S*sizeof(int) );
      }
   } // for (int t=0; t<NLv; t++)

// set up the target local indices
   SetTargetLocalID( NTLocalID, TLocalID );


// 2. get the buffer patches at MainLv and MainLv-1 to RECEIVE particles
//    and the buffer patches at MainLv-1 to SEND particles
// loop over all "real (both leaf and non-leaf)" patches at MainLv with LocalID == 0
   for (int PID0=0; PID0<NReal[0]; PID0+=8)
   {
      SetTargetSibPID0( MainLv, PID0, SibPID0_List );

      for (int s=0; s<26; s++)
      {
         SibPID0 = SibPID0_List[s];

//       check if SibPID0 exists and is a buffer patch
         if ( SibPID0 >= NReal[0] )    // work for both periodic and non-periodic boundary conditions
         {
            for (int Count=0; Count<NTLocalID[s]; Count++)
            {
               SibPID = SibPID0 + TLocalID[s][Count];

//             allocate enough memory for the PID array
               if ( amr->Par->R2B_Buff_NPatchTotal[MainLv][0] >= MemSize_R2B[0] )
               {
                  MemSize_R2B[0] += MemUnit_R2B[0];
                  amr->Par->R2B_Buff_PIDList[MainLv][0] = (int*)realloc( amr->Par->R2B_Buff_PIDList[MainLv][0],
                                                                         MemSize_R2B[0]*sizeof(int) );
               }

//             2-1. store the target sibling-buffer patch index into the R2B list (note that there may be duplicate PID)
               amr->Par->R2B_Buff_PIDList[MainLv][0][ amr->Par->R2B_Buff_NPatchTotal[MainLv][0] ++ ] = SibPID;
            } // for (int Count=0; Count<NTLocalID[s]; Count++)
         } // if ( SibPID0 >= NReal[0] )

         else if ( SibPID0 == -1 )  // work for both periodic and non-periodic boundary conditions
         {
#           ifdef DEBUG_PARTICLE
            if ( MainLv == 0 )
               Aux_Error( ERROR_INFO, "Root-level PID0 %d has no sibling patch along sib %d !!\n", PID0, s );
#           endif

            FaPID = amr->patch[0][MainLv][PID0]->father;

#           ifdef DEBUG_PARTICLE
            if ( FaPID < 0 )
               Aux_Error( ERROR_INFO, "Lv %d, PID0 %d has no father patch (FaPID %d) !!\n", MainLv, PID0, FaPID );
#           endif

            FaSibPID = amr->patch[0][FaLv][FaPID]->sibling[s];

#           ifdef DEBUG_PARTICLE
            if ( FaSibPID < 0 )
               Aux_Error( ERROR_INFO, "Lv %d, PID0 %d, FaPID %d has no sibling [%d] (FaSibPID = %d) !!\n",
                          MainLv, PID0, FaPID, s, FaSibPID );
#           endif

//          check if FaSibPID is a buffer patch
            if ( FaSibPID >= NReal[1] )   // work for both periodic and non-periodic boundary conditions
            {
//             2-2. store the target father-sibling-buffer patch index into the R2B list (note that there may be duplicate PID)
//             allocate enough memory for the PID array
               if ( amr->Par->R2B_Buff_NPatchTotal[MainLv][1] >= MemSize_R2B[1] )
               {
                  MemSize_R2B[1] += MemUnit_R2B[1];
                  amr->Par->R2B_Buff_PIDList[MainLv][1] = (int*)realloc( amr->Par->R2B_Buff_PIDList[MainLv][1],
                                                                         MemSize_R2B[1]*sizeof(int) );
               }

               amr->Par->R2B_Buff_PIDList[MainLv][1][ amr->Par->R2B_Buff_NPatchTotal[MainLv][1] ++ ] = FaSibPID;

//             2-3. store the same patch index into the B2R list (note that there may be duplicate PID)
//                  --> because of the proper-nesting condition, PID0 must be a leaf patch
               if ( amr->Par->B2R_Buff_NPatchTotal[MainLv][1] >= MemSize_B2R[1] )
               {
                  MemSize_B2R[1] += MemUnit_B2R[1];
                  amr->Par->B2R_Buff_PIDList[MainLv][1] = (int*)realloc( amr->Par->B2R_Buff_PIDList[MainLv][1],
                                                                         MemSize_B2R[1]*sizeof(int) );
               }

               amr->Par->B2R_Buff_PIDList[MainLv][1][ amr->Par->B2R_Buff_NPatchTotal[MainLv][1] ++ ] = FaSibPID;
            }
         } // if ( SibPID0 >= NReal[0] ) ... else if ( SibPID0 == -1 )
      } // for (int s=0; s<26; s++)
   } // for (int PID0=0; PID0<NReal[0]; PID0+=8)


// 3. get the buffer patches at MainLv to SEND particles
// loop over all "real leaf" patches at MainLv
   for (int PID=0; PID<NReal[0]; PID++)
   {
//    skip patches with son
      if ( amr->patch[0][MainLv][PID]->son != -1 )    continue;

      for (int s=0; s<26; s++)
      {
         SibPID = amr->patch[0][MainLv][PID]->sibling[s];

//       check if SibPID exists and is a buffer patch
         if ( SibPID >= NReal[0] )  // work for both periodic and non-periodic boundary conditions
         {
//          allocate enough memory for the PID array
            if ( amr->Par->B2R_Buff_NPatchTotal[MainLv][0] >= MemSize_B2R[0] )
            {
               MemSize_B2R[0] += MemUnit_B2R[0];
               amr->Par->B2R_Buff_PIDList[MainLv][0] = (int*)realloc( amr->Par->B2R_Buff_PIDList[MainLv][0],
                                                                      MemSize_B2R[0]*sizeof(int) );
            }

//          store the target sibling-buffer patch index (note that there may be duplicate PID)
            amr->Par->B2R_Buff_PIDList[MainLv][0][ amr->Par->B2R_Buff_NPatchTotal[MainLv][0] ++ ] = SibPID;
         }
      } // for (int s=0; s<26; s++)
   } // for (int PID=0; PID<NReal[0]; PID++)


// 4. sort the candidate lists and remove duplicates (defined as the patches with the same PID)
//    --> note that patches with different PID may still have the same physical coordinates (i.e., PaddedCr1D)
   for (int t=0; t<NLv; t++)
   {
//    4-1. R2B list
      Buff_NPatchTotal_Dup = amr->Par->R2B_Buff_NPatchTotal[MainLv][t];

      Mis_Heapsort( Buff_NPatchTotal_Dup, amr->Par->R2B_Buff_PIDList[MainLv][t], NULL );

      amr->Par->R2B_Buff_NPatchTotal[MainLv][t] = ( Buff_NPatchTotal_Dup > 0 ) ? 1 : 0;

      for (int p=1; p<Buff_NPatchTotal_Dup; p++)
         if ( amr->Par->R2B_Buff_PIDList[MainLv][t][p] != amr->Par->R2B_Buff_PIDList[MainLv][t][p-1] )
            amr->Par->R2B_Buff_PIDList[MainLv][t][ amr->Par->R2B_Buff_NPatchTotal[MainLv][t] ++ ]
               = amr->Par->R2B_Buff_PIDList[MainLv][t][p];

#     ifdef DEBUG_PARTICLE
      for (int p=1; p<amr->Par->R2B_Buff_NPatchTotal[MainLv][t]; p++)
         if ( amr->Par->R2B_Buff_PIDList[MainLv][t][p] <= amr->Par->R2B_Buff_PIDList[MainLv][t][p-1] )
            Aux_Error( ERROR_INFO, "Duplicate PID (Lv %d, t %d, p %d, PID %d, next PID %d) !!\n",
                       MainLv, t, p, amr->Par->R2B_Buff_PIDList[MainLv][t][p-1], amr->Par->R2B_Buff_PIDList[MainLv][t][p] );
#     endif


//    4-2. B2R list
      Buff_NPatchTotal_Dup = amr->Par->B2R_Buff_NPatchTotal[MainLv][t];

      Mis_Heapsort( Buff_NPatchTotal_Dup, amr->Par->B2R_Buff_PIDList[MainLv][t], NULL );

      amr->Par->B2R_Buff_NPatchTotal[MainLv][t] = ( Buff_NPatchTotal_Dup > 0 ) ? 1 : 0;

      for (int p=1; p<Buff_NPatchTotal_Dup; p++)
         if ( amr->Par->B2R_Buff_PIDList[MainLv][t][p] != amr->Par->B2R_Buff_PIDList[MainLv][t][p-1] )
            amr->Par->B2R_Buff_PIDList[MainLv][t][ amr->Par->B2R_Buff_NPatchTotal[MainLv][t] ++ ]
               = amr->Par->B2R_Buff_PIDList[MainLv][t][p];

#     ifdef DEBUG_PARTICLE
      for (int p=1; p<amr->Par->B2R_Buff_NPatchTotal[MainLv][t]; p++)
         if ( amr->Par->B2R_Buff_PIDList[MainLv][t][p] <= amr->Par->B2R_Buff_PIDList[MainLv][t][p-1] )
            Aux_Error( ERROR_INFO, "Duplicate PID (Lv %d, t %d, p %d, PID %d, next PID %d) !!\n",
                       MainLv, t, p, amr->Par->B2R_Buff_PIDList[MainLv][t][p-1], amr->Par->B2R_Buff_PIDList[MainLv][t][p] );
#     endif
   } // for (int t=0; t<NLv; t++)


// 5. get the father patches at MainLv-1 to send particles to their sons
// loop over all real patches at MainLv-1
   if ( FaLv >= 0 )
   for (FaPID=0; FaPID<amr->NPatchComma[FaLv][1]; FaPID++)
   {
//    skip patches without son (son==-1) or with sons at home (son>=0)
      if ( amr->patch[0][FaLv][FaPID]->son >= -1 )    continue;

      for (int s=0; s<26; s++)
      {
         FaSibPID = amr->patch[0][FaLv][FaPID]->sibling[s];

#        ifdef DEBUG_PARTICLE
         if ( FaSibPID == -1 )
            Aux_Error( ERROR_INFO, "FaLv %d, FaPID %d, sib %d, FaSibPID == -1 !!\n", FaLv, FaPID, s );
#        endif

//       store FaPID if any of its sibling is a leaf patch (which can be either real or buffer)
//       --> note that FaSibPID can be <-1 for non-periodic B.C.
         if ( FaSibPID >= 0  &&  amr->patch[0][FaLv][FaSibPID]->son == -1 )
         {
//          allocate enough memory for the PID array
            if ( amr->Par->F2S_Send_NPatchTotal[FaLv] >= MemSize_F2S )
            {
               MemSize_F2S += MemUnit_F2S;
               amr->Par->F2S_Send_PIDList[FaLv] = (int*)realloc( amr->Par->F2S_Send_PIDList[FaLv],
                                                                 MemSize_F2S*sizeof(int) );
            }

//          store the target real father patch index
            amr->Par->F2S_Send_PIDList[FaLv][ amr->Par->F2S_Send_NPatchTotal[FaLv] ++ ] = FaPID;
            break;
         }
      } // for (int s=0; s<26; s++)
   } // for (FaPID=0; FaPID<amr->NPatchComma[FaLv][1]; FaPID++)


// 6. map buffer patches to real patches
   const bool UseInputLBIdx_Yes = true;
   const bool UseInputLBIdx_No  = false;

   for (int t=0; t<NLv; t++)
   {
      lv = RelatedLv[t];

//    6-1. R2B list
      Par_LB_MapBuffer2RealPatch( lv, amr->Par->R2B_Buff_NPatchTotal   [MainLv][t],
                                      amr->Par->R2B_Buff_PIDList       [MainLv][t],
                                      amr->Par->R2B_Buff_NPatchEachRank[MainLv][t],
                                      amr->Par->R2B_Real_NPatchTotal   [MainLv][t],
                                      amr->Par->R2B_Real_PIDList       [MainLv][t],
                                      amr->Par->R2B_Real_NPatchEachRank[MainLv][t],
                                      UseInputLBIdx_No, NULL );
//    6-2. B2R list
      Par_LB_MapBuffer2RealPatch( lv, amr->Par->B2R_Buff_NPatchTotal   [MainLv][t],
                                      amr->Par->B2R_Buff_PIDList       [MainLv][t],
                                      amr->Par->B2R_Buff_NPatchEachRank[MainLv][t],
                                      amr->Par->B2R_Real_NPatchTotal   [MainLv][t],
                                      amr->Par->B2R_Real_PIDList       [MainLv][t],
                                      amr->Par->B2R_Real_NPatchEachRank[MainLv][t],
                                      UseInputLBIdx_No, NULL );
   } // for (int t=0; t<NLv; t++)

// 6-3. F2S list
   if ( FaLv >= 0 )
   {
      int   SonPID;
      long  SonLBIdx;
      long *SonLBIdxList = new long [ amr->Par->F2S_Send_NPatchTotal[FaLv] ];

//    6-3-1. get the LBIdx of sons living abroad
//###NOTE: faster version can only be applied to the Hilbert space-filling curve
      for (int p=0; p<amr->Par->F2S_Send_NPatchTotal[FaLv]; p++)
      {
         FaPID    = amr->Par->F2S_Send_PIDList[FaLv][p];
#        if ( LOAD_BALANCE == HILBERT )
         SonLBIdx = 8*amr->patch[0][FaLv][FaPID]->LB_Idx;   // faster, LB_Idx of one of the eight sons
#        else
         SonLBIdx = LB_Corner2Index( MainLv, amr->patch[0][FaLv][FaPID]->corner, CHECK_ON );    // LB_Idx of son 0
#        endif

         SonLBIdxList[p] = SonLBIdx;
      }

//    6-3-2. find the real patches at MainLv corresponding to SonLBIdxList
      Par_LB_MapBuffer2RealPatch( MainLv, amr->Par->F2S_Send_NPatchTotal   [FaLv],
                                          amr->Par->F2S_Send_PIDList       [FaLv],
                                          amr->Par->F2S_Send_NPatchEachRank[FaLv],
                                          amr->Par->F2S_Recv_NPatchTotal   [FaLv],
                                          amr->Par->F2S_Recv_PIDList       [FaLv],
                                          amr->Par->F2S_Recv_NPatchEachRank[FaLv],
                                          UseInputLBIdx_Yes, SonLBIdxList );

//    check: there should be no particles exchanged within the same rank
#     ifdef DEBUG_PARTICLE
      if ( amr->Par->F2S_Send_NPatchEachRank[FaLv][MPI_Rank] != 0 )
         Aux_Error( ERROR_INFO, "Send %d patches to itself (FaLv %d) !!\n",
                    amr->Par->F2S_Send_NPatchEachRank[FaLv][MPI_Rank], FaLv );

      if ( amr->Par->F2S_Recv_NPatchEachRank[FaLv][MPI_Rank] != 0 )
         Aux_Error( ERROR_INFO, "Recv %d patches from itself (FaLv %d) !!\n",
                    amr->Par->F2S_Recv_NPatchEachRank[FaLv][MPI_Rank], FaLv );
#     endif

//    6-3-3. record the father-buffer patches to actually receive particles
      for (int p=0; p<amr->Par->F2S_Recv_NPatchTotal[FaLv]; p++)
      {
         SonPID = amr->Par->F2S_Recv_PIDList[FaLv][p];
         FaPID  = amr->patch[0][MainLv][SonPID]->father;

//       check: father patches must be buffer patches
#        ifdef DEBUG_PARTICLE
         if ( FaPID < amr->NPatchComma[FaLv][1] )
            Aux_Error( ERROR_INFO, "This is NOT a buffer patch (FaLv %d, SonPID %d, FaPID %d, NReal %d) !!\n",
                       FaLv, SonPID, FaPID, amr->NPatchComma[FaLv][1] );
#        endif

         amr->Par->F2S_Recv_PIDList[FaLv][p] = FaPID;
      }

      delete [] SonLBIdxList;
   } // if ( FaLv >= 0 )


// 7. free memory
   for (int s=0; s<26; s++)   delete [] TLocalID[s];

} // FUNCTION : Par_LB_RecordExchangeParticlePatchID



#endif // #if ( defined PARTICLE  &&  defined LOAD_BALANCE )
