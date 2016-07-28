#include "Copyright.h"
#include "GAMER.h"

#if ( defined PARTICLE  &&  defined LOAD_BALANCE )




//-------------------------------------------------------------------------------------------------------
// Function    :  Par_LB_CollectParticleFromRealPatch
// Description :  Collect particles for the specified buffer patches from the corresponding real patches
//                in all MPI ranks
//
// Note        :  1. Target patches (those in PID_List) must be buffer patches in level "lv"
//                2. Currently this function only collects particle mass and position
//                   --> For particle mass assignment only
//                3. This function is called by Par_LB_CollectParticleFromDescendant
//                4. Array ParMassPos_Away will be allocated for all target patches with particles in the
//                   corresponding real patches
//                   --> Must be deallocated afterward by calling Par_LB_CollectParticleFromDescendant_FreeMemory
//
// Parameter   :  lv                : Target refinement level
//                NSendPatchTotal   : Total number of patches in the PID_List
//                PID_List          : Target PID list
//                PredictPos        : Predict particle position, which is useful for particle mass assignement
//                                    --> We send particle position **after** prediction so that we don't have to
//                                        send particle velocity
//                TargetTime        : Target time for predicting the particle position
//
// Return      :  NPar_Away and ParMassPos_Away (if NPar_Away > 0) for all patches specified in PID_List
//-------------------------------------------------------------------------------------------------------
void Par_LB_CollectParticleFromRealPatch( const int lv, const int NSendPatchTotal, const int *PID_List, const bool PredictPos,
                                          const double TargetTime )
{

// check
#  if ( PAR_MASS >= 4  ||  PAR_POSX >= 4  ||  PAR_POSY >= 4  ||  PAR_POSZ >= 4 )
#     error : ERROR : PAR_MASS, PAR_POSX/Y/Z must be < 4 !!
#  endif

#  ifdef DEBUG_PARTICLE
   if ( lv < 0  ||  lv > MAX_LEVEL )   Aux_Error( ERROR_INFO, "incorrect target level (%d) !!\n", lv );
   if ( NSendPatchTotal < 0 )    Aux_Error( ERROR_INFO, "NSendPatchTotal = %d < 0 !!\n", NSendPatchTotal );
   if ( NSendPatchTotal > 0  &&  PID_List == NULL )
      Aux_Error( ERROR_INFO, "PID_LIst == NULL (NPatch = %d) !!\n", NSendPatchTotal );

   for (int t=0; t<NSendPatchTotal; t++)
   {
      if ( PID_List[t] < amr->NPatchComma[lv][1] )
         Aux_Error( ERROR_INFO, "This function should only be applied to buffer patches (t %d, PID %d, NReal %d) !!\n",
                    t, PID_List[t], amr->NPatchComma[lv][1] );

      if ( amr->patch[0][lv][ PID_List[t] ]->NPar > 0 )
         Aux_Error( ERROR_INFO, "lv %d, PID %d, NPar = %d > 0 !!\n",
                    lv, PID_List[t], amr->patch[0][lv][ PID_List[t] ]->NPar );

      if ( amr->patch[0][lv][ PID_List[t] ]->NPar_Away >= 0 )
         Aux_Error( ERROR_INFO, "lv %d, PID %d, NPar_Away = %d >= 0 !!\n",
                    lv, PID_List[t], amr->patch[0][lv][ PID_List[t] ]->NPar_Away );

      for (int v=0; v<4; v++)
      if ( amr->patch[0][lv][ PID_List[t] ]->ParMassPos_Away[v] != NULL )
         Aux_Error( ERROR_INFO, "lv %d, PID %d, NPar_Away = %d, ParMassPos_Away[%d] != NULL !!\n",
                    lv, PID_List[t], amr->patch[0][lv][ PID_List[t] ]->NPar_Away, v );
   } // for (int t=0; t<NSendPatchTotal; t++)
#  endif // #ifdef DEBUG_PARTICLE


// nothing to do if there is no target patch
   if ( NSendPatchTotal == 0 )   return;


// 1. exchange the target load-balance indices between all ranks
// 1-1. get the LB_Idx of all target patches
// --> note that the LB_Idx stored in each patch always assumes periodicity
// --> so external buffer patches can still find the corresponding real patches easily
   long *LBIdx_List = new long [NSendPatchTotal];
   for (int t=0; t<NSendPatchTotal; t++)  LBIdx_List[t] = amr->patch[0][lv][ PID_List[t] ]->LB_Idx;

// 1-2. get the number of patches exchanged between different ranks
   int *SendBuf_NPatchEachRank = new int [MPI_NRank];
   int *RecvBuf_NPatchEachRank = new int [MPI_NRank];
   int  TRank, NRecvPatchTotal;

   for (int r=0; r<MPI_NRank; r++)  SendBuf_NPatchEachRank[r] = 0;

   for (int t=0; t<NSendPatchTotal; t++)
   {
      TRank = LB_Index2Rank( lv, LBIdx_List[t], CHECK_ON );

      SendBuf_NPatchEachRank[TRank] ++;
   }

   MPI_Alltoall( SendBuf_NPatchEachRank, 1, MPI_INT, RecvBuf_NPatchEachRank, 1, MPI_INT, MPI_COMM_WORLD );

   NRecvPatchTotal = 0;
   for (int r=0; r<MPI_NRank; r++)  NRecvPatchTotal += RecvBuf_NPatchEachRank[r];

// 1-3. set MPI send/recv counts and displacements
   int  *SendDisp_LBIdxEachPatch  = new int [MPI_NRank];
   int  *RecvDisp_LBIdxEachPatch  = new int [MPI_NRank];
   int  *SendCount_LBIdxEachPatch = SendBuf_NPatchEachRank;
   int  *RecvCount_LBIdxEachPatch = RecvBuf_NPatchEachRank;

   SendDisp_LBIdxEachPatch[0] = 0;
   RecvDisp_LBIdxEachPatch[0] = 0;
   for (int r=1; r<MPI_NRank; r++)
   {
      SendDisp_LBIdxEachPatch[r] = SendDisp_LBIdxEachPatch[r-1] + SendCount_LBIdxEachPatch[r-1];
      RecvDisp_LBIdxEachPatch[r] = RecvDisp_LBIdxEachPatch[r-1] + RecvCount_LBIdxEachPatch[r-1];
   }

// 1-4. prepare the send buffer of LBIdx
   long *SendBuf_LBIdxEachPatch   = new long [NSendPatchTotal];
   long *RecvBuf_LBIdxEachPatch   = new long [NRecvPatchTotal];
   int  *Offset_EachRank          = new int  [MPI_NRank];
   int  *PID_List_SortBySendOrder = new int  [NSendPatchTotal];

#  ifdef DEBUG_PARTICLE
   for (int t=0; t<NSendPatchTotal; t++)  PID_List_SortBySendOrder[t] = -1;
#  endif

   for (int r=0; r<MPI_NRank; r++)  Offset_EachRank[r] = SendDisp_LBIdxEachPatch[r];

   for (int t=0; t<NSendPatchTotal; t++)
   {
      TRank = LB_Index2Rank( lv, LBIdx_List[t], CHECK_ON );

      SendBuf_LBIdxEachPatch  [ Offset_EachRank[TRank] ] = LBIdx_List[t];
      PID_List_SortBySendOrder[ Offset_EachRank[TRank] ] = PID_List  [t];

      Offset_EachRank[TRank] ++;
   }

#  ifdef DEBUG_PARTICLE
   for (int t=0; t<NSendPatchTotal; t++)
      if ( PID_List_SortBySendOrder[t] == -1 )  Aux_Error( ERROR_INFO, "PID_List_SortBySendOrder[%d] == -1 !!\n", t );
#  endif

// 1-5. collect LBIdx from all ranks
   MPI_Alltoallv( SendBuf_LBIdxEachPatch, SendCount_LBIdxEachPatch, SendDisp_LBIdxEachPatch, MPI_LONG,
                  RecvBuf_LBIdxEachPatch, RecvCount_LBIdxEachPatch, RecvDisp_LBIdxEachPatch, MPI_LONG, MPI_COMM_WORLD );


// 2. get the PID list corresponding to the received LBIdx list and also the total number of particles to be sent back
   int *RecvBuf_LBIdxEachPatch_IdxTable = new int [NRecvPatchTotal];
   int *Match_LBIdxEachPatch            = new int [NRecvPatchTotal];
   int *PID_List_Map2RecvBuf_LBIdx      = new int [NRecvPatchTotal];
   int *SendBackBuf_NParEachPatch       = new int [NRecvPatchTotal];

   int PID_Match, RecvBuf_Idx, NParThisPatch, NSendBackParTotal = 0;

#  ifdef DEBUG_PARTICLE
   for (int t=0; t<NRecvPatchTotal; t++)  PID_List_Map2RecvBuf_LBIdx[t] = -1;
#  endif

// 2-1. sort and match
   Mis_Heapsort( NRecvPatchTotal, RecvBuf_LBIdxEachPatch, RecvBuf_LBIdxEachPatch_IdxTable );

   Mis_Matching_int( amr->NPatchComma[lv][1], amr->LB->IdxList_Real[lv], NRecvPatchTotal, RecvBuf_LBIdxEachPatch,
                     Match_LBIdxEachPatch );

// loop over all target patches
   for (int t=0; t<NRecvPatchTotal; t++)
   {
#     ifdef DEBUG_PARTICLE
      if ( Match_LBIdxEachPatch[t] == -1 )
         Aux_Error( ERROR_INFO, "LBIdx (%ld) found no match (lv %d) !!\n", RecvBuf_LBIdxEachPatch[t], lv );
#     endif

      PID_Match   = amr->LB->IdxList_Real_IdxTable[lv][ Match_LBIdxEachPatch[t] ];
      RecvBuf_Idx = RecvBuf_LBIdxEachPatch_IdxTable[t];

//    2-2. store the PID mapped to the received LBIdx list
      PID_List_Map2RecvBuf_LBIdx[RecvBuf_Idx] = PID_Match;

//    2-3. count the total number of particles to be sent back
      if ( amr->patch[0][lv][PID_Match]->son == -1 )  NParThisPatch = amr->patch[0][lv][PID_Match]->NPar;
      else                                            NParThisPatch = amr->patch[0][lv][PID_Match]->NPar_Away;

#     ifdef DEBUG_PARTICLE
      if ( NParThisPatch < 0 )
         Aux_Error( ERROR_INFO, "NParThisPatch (%d) has not been calculated (lv %d, PID_Match %d) !!\n",
                    NParThisPatch, lv, PID_Match );
#     endif

      NSendBackParTotal                      += NParThisPatch;
      SendBackBuf_NParEachPatch[RecvBuf_Idx]  = NParThisPatch;
   } // for (int t=0; t<NRecvPatchTotal; t++)

#  ifdef DEBUG_PARTICLE
   for (int t=0; t<NRecvPatchTotal; t++)
      if ( PID_List_Map2RecvBuf_LBIdx[t] < 0  ||  PID_List_Map2RecvBuf_LBIdx[t] >= amr->NPatchComma[lv][1] )
         Aux_Error( ERROR_INFO, "incorrect PID (lv %d, t %d, PID %d, NReal %d) !!\n",
                    lv, t, PID_List_Map2RecvBuf_LBIdx[t], amr->NPatchComma[lv][1] );
#  endif


// 3. prepare the particle data to be sent back
   const int NParVar = 4;  // mass*1 + position*3

   real *SendBackBuf_ParDataEachPatch = new real [NSendBackParTotal*NParVar];

   real  *SendPtr         = SendBackBuf_ParDataEachPatch;
   long  *ParList         = NULL;
   real **ParMassPos_Away = NULL;
   long   ParID;

   for (int t=0; t<NRecvPatchTotal; t++)
   {
      PID_Match     = PID_List_Map2RecvBuf_LBIdx[t];
      NParThisPatch = SendBackBuf_NParEachPatch [t];

//    skip patches with no particles
      if ( NParThisPatch == 0 )  continue;

      if ( amr->patch[0][lv][PID_Match]->son == -1 )
      {
         ParList = amr->patch[0][lv][PID_Match]->ParList;

#        ifdef DEBUG_PARTICLE
         if ( ParList == NULL )
            Aux_Error( ERROR_INFO, "ParList == NULL for NParThisPatch (%d) > 0 (lv %d, PID_Match %d) !!\n",
                       NParThisPatch, lv, PID_Match );
#        endif

         for (int p=0; p<NParThisPatch; p++)
         {
            ParID = ParList[p];

//          here we have assumed that both PAR_MASS, PAR_POSX/Y/Z < NParVar
            SendPtr[PAR_MASS] = amr->Par->ParVar[PAR_MASS][ParID];
            SendPtr[PAR_POSX] = amr->Par->ParVar[PAR_POSX][ParID];
            SendPtr[PAR_POSY] = amr->Par->ParVar[PAR_POSY][ParID];
            SendPtr[PAR_POSZ] = amr->Par->ParVar[PAR_POSZ][ParID];

//          predict particle position to TargetTime
            if ( PredictPos )
            {
//             there should be no particles waiting for velocity correction since we are collecting particles from **higher** levels
//             if ( amr->Par->Time[ParID] < (real)0.0 )  continue;
#              ifdef DEBUG_PARTICLE
               if ( amr->Par->Time[ParID] < (real)0.0 )  Aux_Error( ERROR_INFO, "ParTime[%ld] = %21.14e < 0.0 !!\n",
                    ParID, amr->Par->Time[ParID] );
#              endif

//             note that we don't have to worry about the periodic BC here (in other word, Pos can lie outside the box)
               Par_PredictPos( 1, &ParID, SendPtr+PAR_POSX, SendPtr+PAR_POSY, SendPtr+PAR_POSZ, TargetTime );
            }

            SendPtr += NParVar;
         } // for (int p=0; p<NParThisPatch; p++)
      } // if ( amr->patch[0][lv][PID_Match]->son == -1 )

      else
      {
         ParMassPos_Away = amr->patch[0][lv][PID_Match]->ParMassPos_Away;

#        ifdef DEBUG_PARTICLE
         for (int v=0; v<4; v++)
            if ( ParMassPos_Away[v] == NULL )
               Aux_Error( ERROR_INFO, "ParMassPos_Away[%d] == NULL for NParThisPatch (%d) > 0 (lv %d, PID_Match %d) !!\n",
                          v, NParThisPatch, lv, PID_Match );
#        endif

         for (int p=0; p<NParThisPatch; p++)
         {
//          here we have assumed that both PAR_MASS, PAR_POSX/Y/Z < NParVar
//          (also note that these particle position should have already been predicted to TargetTime
//          by Par_LB_CollectParticleFromDescendant)
            SendPtr[PAR_MASS] = ParMassPos_Away[PAR_MASS][p];
            SendPtr[PAR_POSX] = ParMassPos_Away[PAR_POSX][p];
            SendPtr[PAR_POSY] = ParMassPos_Away[PAR_POSY][p];
            SendPtr[PAR_POSZ] = ParMassPos_Away[PAR_POSZ][p];

            SendPtr += NParVar;
         }
      } // if ( amr->patch[0][lv][PID_Match]->son == -1 ) ... else ...
   } // for (int t=0; t<NRecvPatchTotal; t++)


// 4. send the number of particles and their attributes back
   const bool Exchange_NPatchEachRank_No = false;
   const bool Exchange_LBIdxEachRank_No  = false;

   int  *SendBackBuf_NPatchEachRank   = RecvBuf_NPatchEachRank;
   int  *RecvBackBuf_NPatchEachRank   = SendBuf_NPatchEachRank;
   int  *RecvBackBuf_NParEachPatch    = NULL;   // will be allocated by Par_LB_SendParticle and must be free'd later
   real *RecvBackBuf_ParDataEachPatch = NULL;   // will be allocated by Par_LB_SendParticle and must be free'd later

   long *SendBackBuf_LBIdxEachRank    = NULL;   // useless and does not need to be allocated
   long *RecvBackBuf_LBIdxEachRank    = NULL;   // useless and will not be allocated by Par_LB_SendParticle

   int NRecvBackPatchTotal, NRecvBackParTotal;  // returned from Par_LB_SendParticle

// note that we don't exchange NPatchEachRank (which is already known) and LBIdxEachRank (which is useless)
   Par_LB_SendParticle( NParVar, SendBackBuf_NPatchEachRank, SendBackBuf_NParEachPatch, SendBackBuf_LBIdxEachRank,
                        SendBackBuf_ParDataEachPatch, RecvBackBuf_NPatchEachRank, RecvBackBuf_NParEachPatch,
                        RecvBackBuf_LBIdxEachRank, RecvBackBuf_ParDataEachPatch, NRecvBackPatchTotal, NRecvBackParTotal,
                        Exchange_NPatchEachRank_No, Exchange_LBIdxEachRank_No );

#  ifdef DEBUG_PARTICLE
   if ( NRecvBackPatchTotal != NSendPatchTotal )
      Aux_Error( ERROR_INFO, "Total number of received patches (%d) != expected (%d) !!\n",
                 NRecvBackPatchTotal, NSendPatchTotal );
#  endif

// free the send buffer of particle data in advance to save memory
   delete [] SendBackBuf_ParDataEachPatch;


// 5. store the received particle data to each patch
   const real *RecvPtr = RecvBackBuf_ParDataEachPatch;

   for (int t=0; t<NRecvBackPatchTotal; t++)
   {
      PID_Match     = PID_List_SortBySendOrder [t];
      NParThisPatch = RecvBackBuf_NParEachPatch[t];

//    5-1. set the number of particles in this patch
      amr->patch[0][lv][PID_Match]->NPar_Away = NParThisPatch;

      if ( NParThisPatch > 0 )
      {
//       5-2. allocate the ParMassPos_Away array
         for (int v=0; v<NParVar; v++)
            amr->patch[0][lv][PID_Match]->ParMassPos_Away[v] = new real [NParThisPatch];

         for (int p=0; p<NParThisPatch; p++)
         {
//          5-3. store the particle mass and position
            for (int v=0; v<NParVar; v++)
               amr->patch[0][lv][PID_Match]->ParMassPos_Away[v][p] = *RecvPtr++;

//          5-4. check
#           ifdef DEBUG_PARTICLE
//          we do not transfer inactive particles
            if ( amr->patch[0][lv][PID_Match]->ParMassPos_Away[PAR_MASS][p] < (real)0.0 )
               Aux_Error( ERROR_INFO, "found inactive particle (lv %d, PID_Match %d, Mass %14.7e, particle %d) !!\n",
                          lv, PID_Match, amr->patch[0][lv][PID_Match]->ParMassPos_Away[PAR_MASS][p], p );

//          check if the received particle lies within the target patch (may not when PredictPos is on)
            if ( !PredictPos )
            {
//             always assume periodic B.C. in this check since we don't allocate buffer patches lying outside
//             the simulation domain for non-periodic B.C.
//             --> do NOT use the EdgeL/R stored in each patch for this check since they do not consider periodicity
               const double dh_min     = amr->dh[TOP_LEVEL];
               const int    PatchScale = PS1*amr->scale[lv];
               const real   ParPos[3]  = { amr->patch[0][lv][PID_Match]->ParMassPos_Away[PAR_POSX][p],
                                           amr->patch[0][lv][PID_Match]->ParMassPos_Away[PAR_POSY][p],
                                           amr->patch[0][lv][PID_Match]->ParMassPos_Away[PAR_POSZ][p] };

               double EdgeL, EdgeR;
               int    Cr_Periodic;

               for (int d=0; d<3; d++)
               {
                  Cr_Periodic = ( amr->patch[0][lv][PID_Match]->corner[d] + amr->BoxScale[d] ) % amr->BoxScale[d];
                  EdgeL       = (double)( Cr_Periodic              )*dh_min;
                  EdgeR       = (double)( Cr_Periodic + PatchScale )*dh_min;

                  if ( ParPos[d] < EdgeL  ||  ParPos[d] >= EdgeR )
                     Aux_Error( ERROR_INFO, "wrong home patch (L/R edge = %13.6e/%13.6e, pos[%d] = %13.6e, particle %d, lv %d, PID %d) !!\n",
                                EdgeL, EdgeR, d, ParPos[d], p, lv, PID_Match );
               }
            } // if ( !PredictPos )
#           endif // #ifdef DEBUG_PARTICLE

         } // for (int p=0; p<NParThisPatch; p++)
      } // if ( NParThisPatch > 0 )
   } // for (int t=0; t<NRecvPatchTotal; t++)


// 6. free memory
   delete [] LBIdx_List;
   delete [] SendBuf_NPatchEachRank;
   delete [] RecvBuf_NPatchEachRank;
   delete [] SendDisp_LBIdxEachPatch;
   delete [] RecvDisp_LBIdxEachPatch;
   delete [] SendBuf_LBIdxEachPatch;
   delete [] RecvBuf_LBIdxEachPatch;
   delete [] Offset_EachRank;
   delete [] PID_List_SortBySendOrder;
   delete [] RecvBuf_LBIdxEachPatch_IdxTable;
   delete [] Match_LBIdxEachPatch;
   delete [] PID_List_Map2RecvBuf_LBIdx;
   delete [] SendBackBuf_NParEachPatch;
   delete [] RecvBackBuf_NParEachPatch;
   delete [] RecvBackBuf_ParDataEachPatch;

} // FUNCTION : Par_LB_CollectParticleFromRealPatch



#endif // #if ( defined PARTICLE  &&  defined LOAD_BALANCE )
