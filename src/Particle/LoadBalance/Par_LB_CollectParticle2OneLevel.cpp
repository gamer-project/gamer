#include "Copyright.h"
#include "GAMER.h"

#if ( defined PARTICLE  &&  defined LOAD_BALANCE )




//-------------------------------------------------------------------------------------------------------
// Function    :  Par_LB_CollectParticle2OneLevel
// Description :  Parallel version of Par_CollectParticle2OneLevel for collecting particles from all
//                descendants (sons, grandsons, ...) to patches at a target level
//
// Note        :  1. ParList in all descendants will NOT be changeed after calling this function
//                2. Array ParMassPos_Away will be allocated for patches at FaLv with particles in descendant patches
//                   --> Must be deallocated afterward by calling Par_LB_CollectParticle2OneLevel_FreeMemory
//                3. Do not take into account it's own particles (particles at FaLv)
//                   --> Do nothing if this patch is a leaf patch (ParMassPos_Away array will NOT be allocated)
//                4. This function assumes that all particles live in the leaf patch
//                5. Invoked by Par_CollectParticle2OneLevel
//                6. This function only collects particle mass and position
//                   --> Mainly for depositing particle mass onto grids
//                   --> Particle position are predicted before sending to other ranks so that we don't have to
//                       send particle velocity
//                7. When SibBufPatch is on, this function will also collect particles for sibling-buffer patchesat FaLv
//                   --> Moreover, if FaSibBufPatch is also on, it will also collect particles for
//                       father-sibling-buffer patches at FaLv-1 (if FaLv > 0)
//                       --> Useful for constructing the density field at FaLv for the Poisson solver at FaLv
//
// Parameter   :  FaLv           : Father's refinement level
//                PredictPos     : Predict particle position, which is useful for particle mass assignement
//                                 --> We send particle position **after** prediction so that we don't have to
//                                     send particle velocity
//                TargetTime     : Target time for predicting the particle position
//                SibBufPatch    : true --> Collect particles for sibling-buffer patches at FaLv as well
//                FaSibBufPatch  : true --> Collect particles for father-sibling-buffer patches at FaLv-1 as well
//                                          (do nothing if FaLv==0)
//
// Return      :  NPar_Away and ParMassPos_Away array for all non-leaf patches at FaLv
//-------------------------------------------------------------------------------------------------------
void Par_LB_CollectParticle2OneLevel( const int FaLv, const bool PredictPos, const double TargetTime,
                                      const bool SibBufPatch, const bool FaSibBufPatch )
{

// nothing to do for levels above the max level
// (but note that if SibBufPatch or FaSibBufPatch is on, we need to collect particles for buffer patches
// even when FaLv == MAX_LEVEL)
   if ( FaLv > MAX_LEVEL )    return;


   const int NParVar = 4;  // mass*1 + position*3


// check
#  if ( PAR_MASS >= 4  ||  PAR_POSX >= 4  ||  PAR_POSY >= 4  ||  PAR_POSZ >= 4 )
#     error : ERROR : PAR_MASS, PAR_POSX/Y/Z must be < 4 !!
#  endif

#  ifdef DEBUG_PARTICLE
   for (int FaPID=0; FaPID<amr->NPatchComma[FaLv][1]; FaPID++)
   {
      if ( amr->patch[0][FaLv][FaPID]->NPar_Away != -1 )
         Aux_Error( ERROR_INFO, "particle parameters have been initialized already (FaLv %d, FaPID %d, NPar_Away %d) !!\n",
                    FaLv, FaPID, amr->patch[0][FaLv][FaPID]->NPar_Away );

      for (int v=0; v<NParVar; v++)
      {
         if ( amr->patch[0][FaLv][FaPID]->ParMassPos_Away[v] != NULL )
            Aux_Error( ERROR_INFO, "particle parameters have been initialized already (FaLv %d, FaPID %d, NPar_Away %d, v %d) !!\n",
                       FaLv, FaPID, amr->patch[0][FaLv][FaPID]->NPar_Away, v );
      }
   }
#  endif // #ifdef DEBUG_PARTICLE


// 0. jump to step 5 if the target level is the maximum level --> just collect particles for buffer patches
   if ( FaLv == MAX_LEVEL )
   {
//    0-1. sibling-buffer patches at FaLv
      if ( SibBufPatch )
         Par_LB_CollectParticleFromRealPatch(
            FaLv,
            amr->Par->R2B_Buff_NPatchTotal[FaLv][0], amr->Par->R2B_Buff_PIDList[FaLv][0], amr->Par->R2B_Buff_NPatchEachRank[FaLv][0],
            amr->Par->R2B_Real_NPatchTotal[FaLv][0], amr->Par->R2B_Real_PIDList[FaLv][0], amr->Par->R2B_Real_NPatchEachRank[FaLv][0],
            PredictPos, TargetTime );

//    0-2. father-sibling-buffer patches at FaLv-1
      if ( FaSibBufPatch  &&  FaLv > 0 )
         Par_LB_CollectParticleFromRealPatch(
            FaLv-1,
            amr->Par->R2B_Buff_NPatchTotal[FaLv][1], amr->Par->R2B_Buff_PIDList[FaLv][1], amr->Par->R2B_Buff_NPatchEachRank[FaLv][1],
            amr->Par->R2B_Real_NPatchTotal[FaLv][1], amr->Par->R2B_Real_PIDList[FaLv][1], amr->Par->R2B_Real_NPatchEachRank[FaLv][1],
            PredictPos, TargetTime );

      return;
   }


// 1. prepare the send buffers
// these arrays (except NParForEachRank) will be deallocated by Par_LB_SendParticleData
   int  *NParForEachRank          = new int [MPI_NRank];
   int  *NPatchForEachRank        = new int [MPI_NRank];
   int  *SendBuf_NPatchEachRank   = NPatchForEachRank;
   int  *SendBuf_NParEachPatch    = NULL;
   long *SendBuf_LBIdxEachPatch   = NULL;
   real *SendBuf_ParDataEachPatch = NULL;

   long LB_Idx;
   int  TRank;

#  if ( LOAD_BALANCE != HILBERT )
   const int PatchScaleFaLv = PS1 * amr->scale[FaLv];
   int FaCr[3];
#  endif


// 1-1. get the number of patches and particles sent to each rank
   for (int r=0; r<MPI_NRank; r++)
   {
      NParForEachRank  [r] = 0;
      NPatchForEachRank[r] = 0;
   }

// loop over all real patches at higher levels
   for (int lv=FaLv+1; lv<=MAX_LEVEL; lv++)
   for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
   {
//    skip patches with no particles
      if ( amr->patch[0][lv][PID]->NPar <= 0 )  continue;

#     ifdef DEBUG_PARTICLE
      if ( amr->patch[0][lv][PID]->son != -1 )
         Aux_Error( ERROR_INFO, "non-leaf patch has particles (lv %d, PID %d, SonPID %d, NPar %d) !!\n",
                    lv, PID, amr->patch[0][lv][PID]->son, amr->patch[0][lv][PID]->NPar );

      if ( amr->patch[0][lv][PID]->NPar < 0 )
         Aux_Error( ERROR_INFO, "lv %d, PID %d, NPar %d < 0 !!\n", lv, PID, amr->patch[0][lv][PID]->NPar );
#     endif

//###NOTE: faster version can only be applied to the Hilbert space-filling curve
#     if ( LOAD_BALANCE == HILBERT )
      LB_Idx = amr->patch[0][lv][PID]->LB_Idx / ( 1 << (3*(lv-FaLv)) );
#     else
      for (int d=0; d<3; d++)    FaCr[d] = amr->patch[0][lv][PID]->corner[d] - amr->patch[0][lv][PID]->corner[d]%PatchScaleFaLv;
      LB_Idx = LB_Corner2Index( FaLv, FaCr, CHECK_ON );
#     endif
      TRank  = LB_Index2Rank( FaLv, LB_Idx, CHECK_ON );

      NParForEachRank  [TRank] += amr->patch[0][lv][PID]->NPar;
      NPatchForEachRank[TRank] ++;
   } // for lv, PID


// 1-2. allocate the send buffers
   int NSendPatchTotal=0, NSendParTotal=0;

   for (int r=0; r<MPI_NRank; r++)
   {
      NSendPatchTotal += NPatchForEachRank[r];
      NSendParTotal   += NParForEachRank  [r];
   }

   SendBuf_NParEachPatch    = new int  [NSendPatchTotal];
   SendBuf_LBIdxEachPatch   = new long [NSendPatchTotal];
   SendBuf_ParDataEachPatch = new real [ NSendParTotal*NParVar ];


// 1-3. set the array offsets of each send buffer
   int *Offset_NParEachPatch    = new int [MPI_NRank];
   int *Offset_LBIdxEachPatch   = Offset_NParEachPatch;
   int *Offset_ParDataEachPatch = new int [MPI_NRank];

   Offset_NParEachPatch   [0] = 0;
   Offset_ParDataEachPatch[0] = 0;

   for (int r=1; r<MPI_NRank; r++)
   {
      Offset_NParEachPatch   [r] = Offset_NParEachPatch   [r-1] + NPatchForEachRank[r-1];
      Offset_ParDataEachPatch[r] = Offset_ParDataEachPatch[r-1] + NParForEachRank  [r-1]*NParVar;
   }


// 1-4. fill the send buffers
   long ParID;
   int  NParThisPatch;

   for (int lv=FaLv+1; lv<=MAX_LEVEL; lv++)
   for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
   {
      NParThisPatch = amr->patch[0][lv][PID]->NPar;

//    skip patches with no particles
      if ( NParThisPatch <= 0 )  continue;

//###NOTE: faster version can only be applied to the Hilbert space-filling curve
#     if ( LOAD_BALANCE == HILBERT )
      LB_Idx = amr->patch[0][lv][PID]->LB_Idx / ( 1 << (3*(lv-FaLv)) );
#     else
      for (int d=0; d<3; d++)    FaCr[d] = amr->patch[0][lv][PID]->corner[d] - amr->patch[0][lv][PID]->corner[d]%PatchScaleFaLv;
      LB_Idx = LB_Corner2Index( FaLv, FaCr, CHECK_ON );
#     endif
      TRank  = LB_Index2Rank( FaLv, LB_Idx, CHECK_ON );

//    copy data into the send buffer
#     ifdef DEBUG_PARTICLE
      if ( Offset_NParEachPatch[TRank] >= NSendPatchTotal )
         Aux_Error( ERROR_INFO, "Offset_NParEachPatch[%d] (%d) >= max (%d) !!\n",
                    TRank, Offset_NParEachPatch[TRank], NSendPatchTotal );

      if ( Offset_NParEachPatch[TRank] + PAR_POSZ >= NSendParTotal*NParVar )
         Aux_Error( ERROR_INFO, "Offset_NParEachPatch[%d] + PAR_POSZ (%d) >= max (%d) !!\n",
                    TRank, Offset_NParEachPatch[TRank] + PAR_POSZ, NSendParTotal*NParVar );
#     endif

      SendBuf_NParEachPatch [ Offset_NParEachPatch[TRank] ] = NParThisPatch;
      SendBuf_LBIdxEachPatch[ Offset_NParEachPatch[TRank] ] = LB_Idx;

      for (int p=0; p<NParThisPatch; p++)
      {
         ParID = amr->patch[0][lv][PID]->ParList[p];

         SendBuf_ParDataEachPatch[ Offset_ParDataEachPatch[TRank] + PAR_MASS ] = amr->Par->Mass[ParID];
         SendBuf_ParDataEachPatch[ Offset_ParDataEachPatch[TRank] + PAR_POSX ] = amr->Par->PosX[ParID];
         SendBuf_ParDataEachPatch[ Offset_ParDataEachPatch[TRank] + PAR_POSY ] = amr->Par->PosY[ParID];
         SendBuf_ParDataEachPatch[ Offset_ParDataEachPatch[TRank] + PAR_POSZ ] = amr->Par->PosZ[ParID];

//       predict particle position to TargetTime
         if ( PredictPos )
         {
//          there should be no particles waiting for velocity correction since we are collecting particles from **higher** levels
//          if ( amr->Par->Time[ParID] < (real)0.0 )  continue;
#           ifdef DEBUG_PARTICLE
            if ( amr->Par->Time[ParID] < (real)0.0 )  Aux_Error( ERROR_INFO, "ParTime[%ld] = %21.14e < 0.0 !!\n",
                 ParID, amr->Par->Time[ParID] );
#           endif

//          note that we don't have to worry about the periodic BC here (in other word, Pos can lie outside the box)
            Par_PredictPos( 1, &ParID, SendBuf_ParDataEachPatch + Offset_ParDataEachPatch[TRank] + PAR_POSX,
                                       SendBuf_ParDataEachPatch + Offset_ParDataEachPatch[TRank] + PAR_POSY,
                                       SendBuf_ParDataEachPatch + Offset_ParDataEachPatch[TRank] + PAR_POSZ,
                            TargetTime );
         }

//       update array offset
         Offset_ParDataEachPatch[TRank] += NParVar;
      } // for (int p=0; p<NParThisPatch; p++)

//    update array offset
      Offset_NParEachPatch[TRank] ++;
   } // for lv, PID



// 2. send data to all ranks
// these arrays will be allocated by Par_LB_SendParticleData (using call by reference) and must be free'd later
   int  *RecvBuf_NPatchEachRank   = NULL;
   int  *RecvBuf_NParEachPatch    = NULL;
   long *RecvBuf_LBIdxEachPatch   = NULL;
   real *RecvBuf_ParDataEachPatch = NULL;

// 2-1. exchange data
   const bool Exchange_NPatchEachRank_Yes = true;
   const bool Exchange_LBIdxEachRank_Yes  = true;
   int NRecvPatchTotal, NRecvParTotal;

// note that Par_LB_SendParticleData will also return the total number of patches and particles received (using call by reference)
   Par_LB_SendParticleData( NParVar, SendBuf_NPatchEachRank, SendBuf_NParEachPatch, SendBuf_LBIdxEachPatch,
                            SendBuf_ParDataEachPatch, RecvBuf_NPatchEachRank, RecvBuf_NParEachPatch,
                            RecvBuf_LBIdxEachPatch, RecvBuf_ParDataEachPatch, NRecvPatchTotal, NRecvParTotal,
                            Exchange_NPatchEachRank_Yes, Exchange_LBIdxEachRank_Yes );

// 2-2. free memory
   delete [] SendBuf_NPatchEachRank;
   delete [] SendBuf_NParEachPatch;
   delete [] SendBuf_LBIdxEachPatch;
   delete [] SendBuf_ParDataEachPatch;



// 3. store the received particle data to each patch
// 3-1. LBIdx --> PID
   int *RecvBuf_LBIdxEachPatch_IdxTable = new int [NRecvPatchTotal];
   int *Match_LBIdxEachPatch            = new int [NRecvPatchTotal];
   int  FaPID_Match;

   Mis_Heapsort( NRecvPatchTotal, RecvBuf_LBIdxEachPatch, RecvBuf_LBIdxEachPatch_IdxTable );

   Mis_Matching_int( amr->NPatchComma[FaLv][1], amr->LB->IdxList_Real[FaLv], NRecvPatchTotal, RecvBuf_LBIdxEachPatch,
                     Match_LBIdxEachPatch );


// 3-2. get the number of particles to be allocated for each patch (note that there may be duplicate LBIdx)
//      (also construct the list mapping array index to PID)
   int *FaPIDList = new int [NRecvPatchTotal];
   int  RecvBuf_Idx;

// only set NPar_Away for non-leaf real patches --> leaf real patches will always have NPar_Away == -1
   for (int FaPID=0; FaPID<amr->NPatchComma[FaLv][1]; FaPID++)
      if ( amr->patch[0][FaLv][FaPID]->son != -1 )    amr->patch[0][FaLv][FaPID]->NPar_Away = 0;

#  ifdef DEBUG_PARTICLE
   for (int t=0; t<NRecvPatchTotal; t++)  FaPIDList[t] = -1;
#  endif

   for (int t=0; t<NRecvPatchTotal; t++)
   {
#     ifdef DEBUG_PARTICLE
      if ( Match_LBIdxEachPatch[t] == -1 )
         Aux_Error( ERROR_INFO, "LBIdx (%ld) found no match (FaLv %d) !!\n", RecvBuf_LBIdxEachPatch[t], FaLv );
#     endif

      FaPID_Match = amr->LB->IdxList_Real_IdxTable[FaLv][ Match_LBIdxEachPatch[t] ];
      RecvBuf_Idx = RecvBuf_LBIdxEachPatch_IdxTable[t];

#     ifdef DEBUG_PARTICLE
      if ( amr->patch[0][FaLv][FaPID_Match]->NPar != 0 )
         Aux_Error( ERROR_INFO, "FaLv %d, FaPID_Match %d, NPar = %d != 0 !!\n",
                    FaLv, FaPID_Match, amr->patch[0][FaLv][FaPID_Match]->NPar );

      if ( amr->patch[0][FaLv][FaPID_Match]->son == -1 )
         Aux_Error( ERROR_INFO, "FaLv %d, FaPID_Match %d, SonPID == -1 !!\n", FaLv, FaPID_Match );
#     endif

//    accumulate the number of particles in descendant patches
      amr->patch[0][FaLv][FaPID_Match]->NPar_Away += RecvBuf_NParEachPatch[RecvBuf_Idx];

//    construct the list mapping array index to PID
      FaPIDList[RecvBuf_Idx] = FaPID_Match;
   } // for (int t=0; t<NRecvPatchTotal; t++)

#  ifdef DEBUG_PARTICLE
   for (int t=0; t<NRecvPatchTotal; t++)
      if ( FaPIDList[t] < 0  ||  FaPIDList[t] >= amr->NPatchComma[FaLv][1] )
         Aux_Error( ERROR_INFO, "incorrect PID (FaLv %d, t %d, FaPID %d, NReal %d) !!\n",
                    FaLv, t, FaPIDList[t], amr->NPatchComma[FaLv][1] );
#  endif


// 3-3. allocate the ParMassPos_Away array for each patch
   for (int FaPID=0; FaPID<amr->NPatchComma[FaLv][1]; FaPID++)
   {
      if ( amr->patch[0][FaLv][FaPID]->NPar_Away > 0 )
      {
         for (int v=0; v<NParVar; v++)
            amr->patch[0][FaLv][FaPID]->ParMassPos_Away[v] = new real [ amr->patch[0][FaLv][FaPID]->NPar_Away ];

//       reset to zero since it will be updated again when storing particle data
         amr->patch[0][FaLv][FaPID]->NPar_Away = 0;
      }
   }


// 3-4. store the received particle data
   const real *RecvPtr = RecvBuf_ParDataEachPatch;
   int NPar_Away_Old;

   for (int t=0; t<NRecvPatchTotal; t++)
   {
      FaPID_Match   = FaPIDList[t];
      NPar_Away_Old = amr->patch[0][FaLv][FaPID_Match]->NPar_Away;

      amr->patch[0][FaLv][FaPID_Match]->NPar_Away += RecvBuf_NParEachPatch[t];

#     ifdef DEBUG_PARTICLE
      if ( RecvBuf_NParEachPatch[t] <= 0 )
         Aux_Error( ERROR_INFO, "RecvBuf_NParEachPatch[%d] = %d <= 0 !!\n", t, RecvBuf_NParEachPatch[t] );
#     endif

      for (int p=NPar_Away_Old; p<amr->patch[0][FaLv][FaPID_Match]->NPar_Away; p++)
      {
#        ifdef DEBUG_PARTICLE
         for (int v=0; v<NParVar; v++)
         {
            if ( amr->patch[0][FaLv][FaPID_Match]->ParMassPos_Away[v] == NULL )
               Aux_Error( ERROR_INFO, "particle parameters have NOT been initialized (FaLv %d, FaPID %d, NPar_Away %d, v %d) !!\n",
                          FaLv, FaPID_Match, amr->patch[0][FaLv][FaPID_Match]->NPar_Away, v );
         }
#        endif

         for (int v=0; v<NParVar; v++)
            amr->patch[0][FaLv][FaPID_Match]->ParMassPos_Away[v][p] = *RecvPtr++;

#        ifdef DEBUG_PARTICLE
//       we do not transfer inactive particles
         if ( amr->patch[0][FaLv][FaPID_Match]->ParMassPos_Away[PAR_MASS][p] < (real)0.0 )
            Aux_Error( ERROR_INFO, "found inactive particle (FaLv %d, FaPID %d, Mass %14.7e, particle %d) !!\n",
                       FaLv, FaPID_Match, amr->patch[0][FaLv][FaPID_Match]->ParMassPos_Away[PAR_MASS][p], p );

//       check if the received particle lies within the target patch (may not when PredictPos is on)
         if ( !PredictPos )
         {
            const double *EdgeL     = amr->patch[0][FaLv][FaPID_Match]->EdgeL;
            const double *EdgeR     = amr->patch[0][FaLv][FaPID_Match]->EdgeR;
            const real    ParPos[3] = { amr->patch[0][FaLv][FaPID_Match]->ParMassPos_Away[PAR_POSX][p],
                                        amr->patch[0][FaLv][FaPID_Match]->ParMassPos_Away[PAR_POSY][p],
                                        amr->patch[0][FaLv][FaPID_Match]->ParMassPos_Away[PAR_POSZ][p] };

            for (int d=0; d<3; d++)
            {
               if ( ParPos[d] < EdgeL[d]  ||  ParPos[d] >= EdgeR[d] )
                  Aux_Error( ERROR_INFO, "wrong home patch (L/R edge = %13.6e/%13.6e, pos[%d] = %13.6e, particle %d, FaLv %d, FaPID %d) !!\n",
                             EdgeL[d], EdgeR[d], d, ParPos[d], p, FaLv, FaPID_Match );
            }
         }
#        endif // #ifdef DEBUG_PARTICLE
      } // for (int p=NPar_Away_Old; p<amr->patch[0][FaLv][FaPID_Match]->NPar_Away; p++)
   } // for (int t=0; t<NRecvPatchTotal; t++)


// 4. check
#  ifdef DEBUG_PARTICLE
   long NParLocal_Get=0, NParLocal_Check=0, NParAllRank_Get, NParAllRank_Check;

   for (int FaPID=0; FaPID<amr->NPatchComma[FaLv][1]; FaPID++)
      NParLocal_Get += ( amr->patch[0][FaLv][FaPID]->son == -1 ) ? amr->patch[0][FaLv][FaPID]->NPar :
                                                                   amr->patch[0][FaLv][FaPID]->NPar_Away;

   for (int lv=FaLv; lv<=MAX_LEVEL; lv++)    NParLocal_Check += amr->Par->NPar_Lv[lv];

   MPI_Reduce( &NParLocal_Get,   &NParAllRank_Get,   1, MPI_LONG, MPI_SUM, 0, MPI_COMM_WORLD );
   MPI_Reduce( &NParLocal_Check, &NParAllRank_Check, 1, MPI_LONG, MPI_SUM, 0, MPI_COMM_WORLD );

   if ( MPI_Rank == 0 )
   {
      if ( NParAllRank_Get != NParAllRank_Check )
         Aux_Error( ERROR_INFO, "Total number of active particles >= lv %d = %ld != expected = %ld !!\n",
                    FaLv, NParAllRank_Get, NParAllRank_Check );
   }
#  endif // #ifdef DEBUG_PARTICLE


// 5. collect particles for buffer patches
// 5-1. sibling-buffer patches at FaLv
   if ( SibBufPatch )
      Par_LB_CollectParticleFromRealPatch(
         FaLv,
         amr->Par->R2B_Buff_NPatchTotal[FaLv][0], amr->Par->R2B_Buff_PIDList[FaLv][0], amr->Par->R2B_Buff_NPatchEachRank[FaLv][0],
         amr->Par->R2B_Real_NPatchTotal[FaLv][0], amr->Par->R2B_Real_PIDList[FaLv][0], amr->Par->R2B_Real_NPatchEachRank[FaLv][0],
         PredictPos, TargetTime );

// 5-2. father-sibling-buffer patches at FaLv-1
   if ( FaSibBufPatch  &&  FaLv > 0 )
      Par_LB_CollectParticleFromRealPatch(
         FaLv-1,
         amr->Par->R2B_Buff_NPatchTotal[FaLv][1], amr->Par->R2B_Buff_PIDList[FaLv][1], amr->Par->R2B_Buff_NPatchEachRank[FaLv][1],
         amr->Par->R2B_Real_NPatchTotal[FaLv][1], amr->Par->R2B_Real_PIDList[FaLv][1], amr->Par->R2B_Real_NPatchEachRank[FaLv][1],
         PredictPos, TargetTime );


// 6. free memory
   delete [] RecvBuf_NPatchEachRank;
   delete [] RecvBuf_NParEachPatch;
   delete [] RecvBuf_LBIdxEachPatch;
   delete [] RecvBuf_ParDataEachPatch;

   delete [] NParForEachRank;
   delete [] Offset_NParEachPatch;
   delete [] Offset_ParDataEachPatch;
   delete [] RecvBuf_LBIdxEachPatch_IdxTable;
   delete [] Match_LBIdxEachPatch;
   delete [] FaPIDList;

} // FUNCTION : Par_LB_CollectParticle2OneLevel



//-------------------------------------------------------------------------------------------------------
// Function    :  Par_LB_CollectParticle2OneLevel_FreeMemory
// Description :  Release the memory allocated by Par_LB_CollectParticle2OneLevel
//
// Note        :  1. Invoded by Par_CollectParticle2OneLevel_FreeMemory
//
// Parameter   :  lv           : Target refinement level
//                SibBufPatch    : true --> Release memory for sibling-buffer patches at lv as well
//                FaSibBufPatch  : true --> Release memory for father-sibling-buffer patches at lv-1 as well
//                                          (do nothing if lv==0)
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Par_LB_CollectParticle2OneLevel_FreeMemory( const int lv, const bool SibBufPatch, const bool FaSibBufPatch )
{

   const int NParVar = 4;


// 1. real patches at lv
   for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
   {
      for (int v=0; v<NParVar; v++)
      {
         if ( amr->patch[0][lv][PID]->ParMassPos_Away[v] != NULL )
         {
            delete [] amr->patch[0][lv][PID]->ParMassPos_Away[v];
            amr->patch[0][lv][PID]->ParMassPos_Away[v] = NULL;
         }
      }

//    -1 : indicating that NPar_Away is not calculated yet
      amr->patch[0][lv][PID]->NPar_Away = -1;
   }


// 2. sibling-buffer patches at lv
   if ( SibBufPatch )
   for (int p=0; p<amr->Par->R2B_Buff_NPatchTotal[lv][0]; p++)
   {
      const int PID = amr->Par->R2B_Buff_PIDList[lv][0][p];

      for (int v=0; v<NParVar; v++)
      {
         if ( amr->patch[0][lv][PID]->ParMassPos_Away[v] != NULL )
         {
            delete [] amr->patch[0][lv][PID]->ParMassPos_Away[v];
            amr->patch[0][lv][PID]->ParMassPos_Away[v] = NULL;
         }
      }

//    -1 : indicating that NPar_Away is not calculated yet
      amr->patch[0][lv][PID]->NPar_Away = -1;
   }


// 3. father-sibling-buffer patches at lv-1
   const int FaLv = lv - 1;

   if ( FaSibBufPatch  &&  FaLv >= 0 )
   for (int p=0; p<amr->Par->R2B_Buff_NPatchTotal[lv][1]; p++)
   {
      const int FaPID = amr->Par->R2B_Buff_PIDList[lv][1][p];

      for (int v=0; v<NParVar; v++)
      {
         if ( amr->patch[0][FaLv][FaPID]->ParMassPos_Away[v] != NULL )
         {
            delete [] amr->patch[0][FaLv][FaPID]->ParMassPos_Away[v];
            amr->patch[0][FaLv][FaPID]->ParMassPos_Away[v] = NULL;
         }
      }

//    -1 : indicating that NPar_Away is not calculated yet
      amr->patch[0][FaLv][FaPID]->NPar_Away = -1;
   }


// check: if we do everthing correctly, no patches (either real or buffer patches) at lv and lv-1
//        should have ParMassPos_Away allocated
#  ifdef DEBUG_PARTICLE
   for (int TLv=lv; (TLv>=lv-1 && TLv>=0); TLv--)
   {
//    loop over all real and buffer patches
      for (int PID=0; PID<amr->num[TLv]; PID++)
      {
         for (int v=0; v<NParVar; v++)
         if ( amr->patch[0][TLv][PID]->ParMassPos_Away[v] != NULL )
            Aux_Error( ERROR_INFO, "lv %d, PID %d, v %d, ParMassPos_Away != NULL !!\n",
                       TLv, PID, v );

         if ( amr->patch[0][TLv][PID]->NPar_Away != -1 )
            Aux_Error( ERROR_INFO, "lv %d, PID %d, NPar_Away = %d != -1 !!\n",
                       TLv, PID, amr->patch[0][TLv][PID]->NPar_Away );
      }
   }
#  endif // #ifdef DEBUG_PARTICLE

} // FUNCTION : Par_LB_CollectParticle2OneLevel_FreeMemory



#endif // #if ( defined PARTICLE  &&  defined LOAD_BALANCE )
