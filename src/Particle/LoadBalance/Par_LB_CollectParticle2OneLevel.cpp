#include "GAMER.h"

#if ( defined PARTICLE  &&  defined LOAD_BALANCE )

#ifdef TIMING
extern Timer_t *Timer_Par_MPI[NLEVEL][6];
#endif




//-------------------------------------------------------------------------------------------------------
// Function    :  Par_LB_CollectParticle2OneLevel
// Description :  Parallel version of Par_CollectParticle2OneLevel for collecting particles from all
//                descendants (sons, grandsons, ...) to patches at a target level
//
// Note        :  1. ParList[] in all descendants will NOT be changeed after calling this function
//                2. ParMassPos_Copy[] will be allocated for the target patches at FaLv
//                   --> Must be deallocated afterward by calling Par_LB_CollectParticle2OneLevel_FreeMemory()
//                3. Only work on non-leaf real patches and buffer patches (which can be either real or buffer)
//                   --> For leaf real patches the ParList_Copy[] will NOT be allocated
//                   --> However, it does take into account the corner case where non-leaf real patches may
//                       have NPar>0 temporarily after updating particle position.
//                       --> It's because particles travelling from coarse to fine grids will stay in coarse grids
//                           temporarily until the velocity correction is done.
//                       --> For these patches, NPar_Copy will be **the sum of NPar and the number of particles
//                           collected from other patches**, and ParList_Copy[] (or ParMassPos_Copy[]) will contain
//                           information of particles belonging to NPar as well.
//                       --> It makes implementation simplier. **For leaf real patches, one only needs to consider
//                           NPar and ParList[]. While for all other patches, one only needs to consider NPar_Copy,
//                           ParList_Copy[] (or ParMassPos_Copy[]). One never needs to consider both.**
//                4. Invoked by Par_CollectParticle2OneLevel()
//                5. This function only collects particle mass and position
//                   --> Mainly for depositing particle mass onto grids
//                   --> Particle position are predicted before sending to other ranks so that we don't have to
//                       send particle velocity
//                       --> Accordingly, if Par->PredictPos is on, ParMassPos_Copy[] always stores the data
//                           after position prediction
//                6. When SibBufPatch is on, this function will also collect particles for sibling-buffer patches at FaLv
//                   --> Moreover, if FaSibBufPatch is also on, it will also collect particles for
//                       father-sibling-buffer patches at FaLv-1 (if FaLv > 0)
//                       --> Useful for constructing the density field at FaLv for the Poisson solver at FaLv
//                       --> These data will be stored in NPar_Copy and ParMassPos_Copy[] as well
//                7. Option "JustCountNPar" can be used to count the number of particles in each real patch at FaLv
//                   --> Do NOT collect other particle data (e.g., particle mass and position)
//                   --> Particle count is stored in NPar_Copy
//                   --> ParMassPos_Copy[] will NOT be allocated
//                   --> Currently it does NOT work with the options SibBufPatch and FaSibBufPatch
//                       --> It only counts particles for **real** patches at FaLv
//                   --> Does NOT work with "PredictPos"
//
// Parameter   :  FaLv          : Father's refinement level
//                PredictPos    : Predict particle position, which is useful for particle mass assignement
//                                --> We send particle position **after** prediction so that we don't have to
//                                    send particle velocity
//                TargetTime    : Target time for predicting the particle position
//                SibBufPatch   : true --> Collect particles for sibling-buffer patches at FaLv as well
//                FaSibBufPatch : true --> Collect particles for father-sibling-buffer patches at FaLv-1 as well
//                                         (do nothing if FaLv==0)
//                JustCountNPar : Just count the number of particles in each real patch at FaLv. Don't collect
//                                other particle data (e.g., particle mass and position)
//                TimingSendPar : Measure the elapsed time of the routine "Par_LB_SendParticleData"
//
// Return      :  NPar_Copy and ParMassPos_Copy[] (if JustCountNPar == false) for all non-leaf patches at FaLv
//                (and for sibling-buffer patches at FaLv if SibBufPatch is on, and for father-sibling-buffer
//                patches at FaLv-1 if FaSibBufPatch is on and FaLv>0)
//-------------------------------------------------------------------------------------------------------
void Par_LB_CollectParticle2OneLevel( const int FaLv, const bool PredictPos, const double TargetTime,
                                      const bool SibBufPatch, const bool FaSibBufPatch, const bool JustCountNPar,
                                      const bool TimingSendPar )
{

// nothing to do for levels above the max level
// (but note that if SibBufPatch or FaSibBufPatch is on, we need to collect particles for buffer patches
// even when FaLv == MAX_LEVEL)
   if ( FaLv > MAX_LEVEL )    return;


   const int NParAtt = 4;  // mass*1 + position*3

   Timer_t *Timer[3] = { NULL, NULL, NULL };

#  ifdef TIMING
   if ( TimingSendPar )
   {
      Timer[0] = Timer_Par_MPI[FaLv][3];
      Timer[1] = Timer_Par_MPI[FaLv][4];
      Timer[2] = Timer_Par_MPI[FaLv][5];
   }
#  endif


// check
#  if ( PAR_MASS >= 4  ||  PAR_POSX >= 4  ||  PAR_POSY >= 4  ||  PAR_POSZ >= 4 )
#     error : ERROR : PAR_MASS, PAR_POSX/Y/Z must be < 4 !!
#  endif

#  ifdef DEBUG_PARTICLE
   if ( JustCountNPar )
   {
      if ( PredictPos )       Aux_Error( ERROR_INFO, "JustCountNPar does NOT work with PredictPos !!\n" );
      if ( SibBufPatch )      Aux_Error( ERROR_INFO, "JustCountNPar does NOT work with SibBufPatch !!\n" );
      if ( FaSibBufPatch )    Aux_Error( ERROR_INFO, "JustCountNPar does NOT work with FaSibBufPatch !!\n" );
   }

   for (int FaPID=0; FaPID<amr->NPatchComma[FaLv][1]; FaPID++)
   {
      if ( amr->patch[0][FaLv][FaPID]->NPar_Copy != -1 )
         Aux_Error( ERROR_INFO, "particle parameters have been initialized already (FaLv %d, FaPID %d, NPar_Copy %d) !!\n",
                    FaLv, FaPID, amr->patch[0][FaLv][FaPID]->NPar_Copy );

      for (int v=0; v<NParAtt; v++)
      {
         if ( amr->patch[0][FaLv][FaPID]->ParMassPos_Copy[v] != NULL )
            Aux_Error( ERROR_INFO, "particle parameters have been initialized already (FaLv %d, FaPID %d, NPar_Copy %d, v %d) !!\n",
                       FaLv, FaPID, amr->patch[0][FaLv][FaPID]->NPar_Copy, v );
      }
   }
#  endif // #ifdef DEBUG_PARTICLE


// 0. jump to step 5 if FaLv+1 is guaranteed to have no particle
//    --> just collect particles for buffer patches
//    --> note that we don't have to set NPar_Copy here since leaf real patches always have NPar_Copy == -1
   if ( FaLv == TOP_LEVEL  ||  NPatchTotal[FaLv+1] == 0 )
   {
//    0-1. sibling-buffer patches at FaLv
      if ( SibBufPatch  &&  !JustCountNPar )
      {
         char Timer_Comment[50];
         sprintf( Timer_Comment, "%3d %15s", FaLv, "Par_Coll2Sib" );

         Par_LB_CollectParticleFromRealPatch(
            FaLv,
            amr->Par->R2B_Buff_NPatchTotal[FaLv][0], amr->Par->R2B_Buff_PIDList[FaLv][0], amr->Par->R2B_Buff_NPatchEachRank[FaLv][0],
            amr->Par->R2B_Real_NPatchTotal[FaLv][0], amr->Par->R2B_Real_PIDList[FaLv][0], amr->Par->R2B_Real_NPatchEachRank[FaLv][0],
            PredictPos, TargetTime, Timer[1], Timer_Comment );
      }

//    0-2. father-sibling-buffer patches at FaLv-1
      if ( FaSibBufPatch  &&  FaLv > 0  &&  !JustCountNPar )
      {
         char Timer_Comment[50];
         sprintf( Timer_Comment, "%3d %15s", FaLv, "Par_Coll2FaSib" );

         Par_LB_CollectParticleFromRealPatch(
            FaLv-1,
            amr->Par->R2B_Buff_NPatchTotal[FaLv][1], amr->Par->R2B_Buff_PIDList[FaLv][1], amr->Par->R2B_Buff_NPatchEachRank[FaLv][1],
            amr->Par->R2B_Real_NPatchTotal[FaLv][1], amr->Par->R2B_Real_PIDList[FaLv][1], amr->Par->R2B_Real_NPatchEachRank[FaLv][1],
            PredictPos, TargetTime, Timer[2], Timer_Comment );
      }

      return;
   } // if ( FaLv == MAX_LEVEL )


// 1. prepare the send buffers
   int  *NParForEachRank          = new int [MPI_NRank];
   int  *NPatchForEachRank        = new int [MPI_NRank];
   int  *SendBuf_NPatchEachRank   = NPatchForEachRank;
   int  *SendBuf_NParEachPatch    = NULL;
   long *SendBuf_LBIdxEachPatch   = NULL;
   real *SendBuf_ParDataEachPatch = NULL;

#  if ( LOAD_BALANCE != HILBERT )
   const int PatchScaleFaLv = PS1 * amr->scale[FaLv];
#  endif


// 1-1. get the number of patches and particles sent to each rank
   int (*NPatchAcc) = new int [MAX_LEVEL-FaLv];    // accumulative number of patches at each target level
   int NPatchAll;

   NPatchAll    = 0;    // note that it includes patches **without** particles --> different from NSendPatchTotal declared later
   NPatchAcc[0] = 0;    // note that it includes patches **without** particles

   for (int lv=FaLv+1, q=0; lv<=MAX_LEVEL; lv++, q++)
   {
      NPatchAll += amr->NPatchComma[lv][1];

      if ( q > 0 )   NPatchAcc[q] = NPatchAcc[q-1] + amr->NPatchComma[lv-1][1];
   }

   int  *TRank_AllPatch = new int  [NPatchAll];
   long *LBIdx_AllPatch = new long [NPatchAll];

   for (int r=0; r<MPI_NRank; r++)
   {
      NParForEachRank  [r] = 0;
      NPatchForEachRank[r] = 0;
   }

// loop over all real patches at higher levels
   for (int lv=FaLv+1, q=0; lv<=MAX_LEVEL; lv++, q++)
   for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
   {
//    skip patches without particles
      if ( amr->patch[0][lv][PID]->NPar == 0 )  continue;

#     ifdef DEBUG_PARTICLE
      if ( amr->patch[0][lv][PID]->son != -1 )
         Aux_Error( ERROR_INFO, "non-leaf patch has particles (lv %d, PID %d, SonPID %d, NPar %d) !!\n",
                    lv, PID, amr->patch[0][lv][PID]->son, amr->patch[0][lv][PID]->NPar );

      if ( amr->patch[0][lv][PID]->NPar < 0 )
         Aux_Error( ERROR_INFO, "lv %d, PID %d, NPar %d < 0 !!\n", lv, PID, amr->patch[0][lv][PID]->NPar );
#     endif

//###NOTE: faster version can only be applied to the Hilbert space-filling curve
#     if ( LOAD_BALANCE == HILBERT )
      const long LB_Idx = amr->patch[0][lv][PID]->LB_Idx / ( 1 << (3*(lv-FaLv)) );
#     else
      int FaCr[3];
      for (int d=0; d<3; d++)    FaCr[d] = amr->patch[0][lv][PID]->corner[d] - amr->patch[0][lv][PID]->corner[d]%PatchScaleFaLv;
      const long LB_Idx = LB_Corner2Index( FaLv, FaCr, CHECK_ON );
#     endif
      const int TRank   = LB_Index2Rank( FaLv, LB_Idx, CHECK_ON );
      const int AccIdx  = NPatchAcc[q] + PID;

//    both NParForEachRank[] and NPatchForEachRank[] do NOT include patches without particles
      NParForEachRank  [TRank ] += amr->patch[0][lv][PID]->NPar;
      NPatchForEachRank[TRank ] ++;
      TRank_AllPatch   [AccIdx] = TRank;
      LBIdx_AllPatch   [AccIdx] = LB_Idx;
   } // for lv, PID


// 1-2. allocate the send buffers
// both NSendPatchTotal and NSendParTotal do NOT include patches without particles
   int NSendPatchTotal=0, NSendParTotal=0;

   for (int r=0; r<MPI_NRank; r++)
   {
      NSendPatchTotal += NPatchForEachRank[r];
      NSendParTotal   += NParForEachRank  [r];
   }

   SendBuf_NParEachPatch  = new int  [NSendPatchTotal];
   SendBuf_LBIdxEachPatch = new long [NSendPatchTotal];

// reuse the MPI send buffer declared in LB_GetBufferData for better MPI performance
   if ( !JustCountNPar )   SendBuf_ParDataEachPatch = LB_GetBufferData_MemAllocate_Send( NSendParTotal*NParAtt );


// 1-3. set the array offsets of the send buffer of each rank
   int *OffsetEachRank_Patch   = new int [MPI_NRank];
   int *OffsetEachRank_ParData = new int [MPI_NRank];    // actually useless in the JustCountNPar mode

   OffsetEachRank_Patch  [0] = 0;
   OffsetEachRank_ParData[0] = 0;

   for (int r=1; r<MPI_NRank; r++)
   {
      OffsetEachRank_Patch  [r] = OffsetEachRank_Patch  [r-1] + NPatchForEachRank[r-1];
      OffsetEachRank_ParData[r] = OffsetEachRank_ParData[r-1] + NParForEachRank  [r-1]*NParAtt;
   }


// 1-4. set the array offsets of the send buffer of each patch (mainly for the OpenMP parallelization)
   int *OffsetEachPatch_Patch   = new int [NPatchAll];
   int *OffsetEachPatch_ParData = new int [NPatchAll];   // actually useless in the JustCountNPar mode

   for (int lv=FaLv+1, q=0; lv<=MAX_LEVEL; lv++, q++)
   for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
   {
      const int NParThisPatch = amr->patch[0][lv][PID]->NPar;
      const int AccIdx        = NPatchAcc[q] + PID;
      const int TRank         = TRank_AllPatch[AccIdx];

//    skip patches without particles (must skip since OffsetEachRank_* are not set for patches without particles)
      if ( NParThisPatch == 0 )  continue;

      OffsetEachPatch_Patch  [AccIdx] = OffsetEachRank_Patch  [TRank];
      OffsetEachPatch_ParData[AccIdx] = OffsetEachRank_ParData[TRank];

      OffsetEachRank_Patch  [TRank] += 1;
      OffsetEachRank_ParData[TRank] += NParThisPatch*NParAtt;
   }

   delete [] OffsetEachRank_Patch;
   delete [] OffsetEachRank_ParData;


// 1-5. fill the send buffers
   for (int lv=FaLv+1, q=0; lv<=MAX_LEVEL; lv++, q++)
   {
#     pragma omp parallel for schedule( PAR_OMP_SCHED, PAR_OMP_SCHED_CHUNK )
      for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
      {
         const int  NParThisPatch = amr->patch[0][lv][PID]->NPar;
         const int  AccIdx        = NPatchAcc[q] + PID;
         const long LBIdx         = LBIdx_AllPatch[AccIdx];

//       skip patches without particles (must skip since many variables are not set for patches without particles)
         if ( NParThisPatch == 0 )  continue;

//       copy data into the send buffer
#        ifdef DEBUG_PARTICLE
         if ( OffsetEachPatch_Patch[AccIdx] >= NSendPatchTotal )
            Aux_Error( ERROR_INFO, "OffsetEachPatch_Patch[%d] (%d) >= max (%d) !!\n",
                       AccIdx, OffsetEachPatch_Patch[AccIdx], NSendPatchTotal );

         if ( OffsetEachPatch_ParData[AccIdx] + PAR_POSZ >= NSendParTotal*NParAtt )
            Aux_Error( ERROR_INFO, "OffsetEachPatch_ParData[%d] + PAR_POSZ (%d) >= max (%d) !!\n",
                       AccIdx, OffsetEachPatch_ParData[AccIdx] + PAR_POSZ, NSendParTotal*NParAtt );
#        endif

         SendBuf_NParEachPatch [ OffsetEachPatch_Patch[AccIdx] ] = NParThisPatch;
         SendBuf_LBIdxEachPatch[ OffsetEachPatch_Patch[AccIdx] ] = LBIdx;

         if ( !JustCountNPar )
         {
            real *SendPtr = SendBuf_ParDataEachPatch + OffsetEachPatch_ParData[AccIdx];

            for (int p=0; p<NParThisPatch; p++)
            {
               const long ParID = amr->patch[0][lv][PID]->ParList[p];

               *( SendPtr + PAR_MASS ) = amr->Par->Mass[ParID];
               *( SendPtr + PAR_POSX ) = amr->Par->PosX[ParID];
               *( SendPtr + PAR_POSY ) = amr->Par->PosY[ParID];
               *( SendPtr + PAR_POSZ ) = amr->Par->PosZ[ParID];

//             predict particle position to TargetTime
               if ( PredictPos )
               {
//                there should be no particles waiting for velocity correction since these particles are collected from **higher** levels
//                if ( amr->Par->Time[ParID] < (real)0.0 )  continue;
#                 ifdef DEBUG_PARTICLE
                  if ( amr->Par->Time[ParID] < (real)0.0 )  Aux_Error( ERROR_INFO, "ParTime[%ld] = %21.14e < 0.0 !!\n",
                       ParID, amr->Par->Time[ParID] );
#                 endif

//                note that we don't have to worry about the periodic BC here (in other words, Pos can lie outside the box)
                  Par_PredictPos( 1, &ParID, SendPtr+PAR_POSX, SendPtr+PAR_POSY, SendPtr+PAR_POSZ, TargetTime );
               }

//             update array offset
               SendPtr += NParAtt;
            } // for (int p=0; p<NParThisPatch; p++)
         } // if ( !JustCountNPar )
      } // for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
   } // for (int lv=FaLv+1, q=0; lv<=MAX_LEVEL; lv++, q++)


// 1-6. free memory
   delete [] OffsetEachPatch_Patch;
   delete [] OffsetEachPatch_ParData;
   delete [] NParForEachRank;
   delete [] TRank_AllPatch;
   delete [] LBIdx_AllPatch;
   delete [] NPatchAcc;



// 2. send data to all ranks
// these arrays will be allocated by Par_LB_SendParticleData (using call by reference) and must be free'd later
// --> except for RecvBuf_ParDataEachPatch, which is just a pointer to the MPI recv buffer declared in LB_GetBufferData
   int  *RecvBuf_NPatchEachRank   = NULL;
   int  *RecvBuf_NParEachPatch    = NULL;
   long *RecvBuf_LBIdxEachPatch   = NULL;
   real *RecvBuf_ParDataEachPatch = NULL;

// 2-1. exchange data
   const bool Exchange_NPatchEachRank_Yes = true;
   const bool Exchange_LBIdxEachRank_Yes  = true;
   const bool Exchange_ParDataEachRank    = !JustCountNPar;
   int  NRecvPatchTotal, NRecvParTotal;

   char Timer_Comment[50];
   sprintf( Timer_Comment, "%3d %15s", FaLv, "Par_Collect" );

// note that Par_LB_SendParticleData will also return the total number of patches and particles received (using call by reference)
   Par_LB_SendParticleData( NParAtt, SendBuf_NPatchEachRank, SendBuf_NParEachPatch, SendBuf_LBIdxEachPatch,
                            SendBuf_ParDataEachPatch, NSendParTotal, RecvBuf_NPatchEachRank, RecvBuf_NParEachPatch,
                            RecvBuf_LBIdxEachPatch, RecvBuf_ParDataEachPatch, NRecvPatchTotal, NRecvParTotal,
                            Exchange_NPatchEachRank_Yes, Exchange_LBIdxEachRank_Yes, Exchange_ParDataEachRank,
                            Timer[0], Timer_Comment );

// 2-2. free memory
   delete [] SendBuf_NPatchEachRank;
   delete [] SendBuf_NParEachPatch;
   delete [] SendBuf_LBIdxEachPatch;



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
   int *FaPIDList = new int [NRecvPatchTotal];  // actually useless in the JustCountNPar mode
   int  RecvBuf_Idx;

// initialize NPar_Copy as NPar (instead of 0) for non-leaf real patches
// --> add the number of particles (i.e., NPar) temporarily residing in this patch waiting for the velocity correction in KDK
// --> so only leaf real patches will always have NPar_Copy == -1 after calling this function
   for (int FaPID=0; FaPID<amr->NPatchComma[FaLv][1]; FaPID++)
      if ( amr->patch[0][FaLv][FaPID]->son != -1 )    amr->patch[0][FaLv][FaPID]->NPar_Copy = amr->patch[0][FaLv][FaPID]->NPar;

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
      if ( amr->patch[0][FaLv][FaPID_Match]->son == -1 )
         Aux_Error( ERROR_INFO, "FaLv %d, FaPID_Match %d, SonPID == -1 !!\n", FaLv, FaPID_Match );
#     endif

//    accumulate the number of particles in descendant patches
      amr->patch[0][FaLv][FaPID_Match]->NPar_Copy += RecvBuf_NParEachPatch[RecvBuf_Idx];

//    construct the list mapping array index to PID
      FaPIDList[RecvBuf_Idx] = FaPID_Match;
   } // for (int t=0; t<NRecvPatchTotal; t++)

#  ifdef DEBUG_PARTICLE
   for (int t=0; t<NRecvPatchTotal; t++)
      if ( FaPIDList[t] < 0  ||  FaPIDList[t] >= amr->NPatchComma[FaLv][1] )
         Aux_Error( ERROR_INFO, "incorrect PID (FaLv %d, t %d, FaPID %d, NReal %d) !!\n",
                    FaLv, t, FaPIDList[t], amr->NPatchComma[FaLv][1] );
#  endif


// 3-3. allocate the ParMassPos_Copy array for each patch
   if ( !JustCountNPar )
   for (int FaPID=0; FaPID<amr->NPatchComma[FaLv][1]; FaPID++)
   {
      if ( amr->patch[0][FaLv][FaPID]->NPar_Copy > 0 )
      {
         for (int v=0; v<NParAtt; v++)
            amr->patch[0][FaLv][FaPID]->ParMassPos_Copy[v] = new real [ amr->patch[0][FaLv][FaPID]->NPar_Copy ];

//       reset to zero (instead of NPar) since we will use NPar_Copy to record the number of particles that has been
//       added to the ParMassPos_Copy array
         amr->patch[0][FaLv][FaPID]->NPar_Copy = 0;
      }
   }


// 3-4. store the received particle data
// --> currently we do not implement OpenMP here since different received patches at lv>FaLv may map to the same patch at FaLv
// --> one cannot naively parallelize the "for (int t=0; t<NRecvPatchTotal; t++)" loop
   const real *RecvPtr = RecvBuf_ParDataEachPatch;
   int NPar_Copy_Old;

   if ( !JustCountNPar )
   for (int t=0; t<NRecvPatchTotal; t++)
   {
//    different "t" may map to the same "FaPID_Match"
      FaPID_Match   = FaPIDList[t];
      NPar_Copy_Old = amr->patch[0][FaLv][FaPID_Match]->NPar_Copy;

      amr->patch[0][FaLv][FaPID_Match]->NPar_Copy += RecvBuf_NParEachPatch[t];

#     ifdef DEBUG_PARTICLE
      if ( RecvBuf_NParEachPatch[t] <= 0 )
         Aux_Error( ERROR_INFO, "RecvBuf_NParEachPatch[%d] = %d <= 0 !!\n", t, RecvBuf_NParEachPatch[t] );
#     endif

      for (int p=NPar_Copy_Old; p<amr->patch[0][FaLv][FaPID_Match]->NPar_Copy; p++)
      {
#        ifdef DEBUG_PARTICLE
         for (int v=0; v<NParAtt; v++)
         {
            if ( amr->patch[0][FaLv][FaPID_Match]->ParMassPos_Copy[v] == NULL )
               Aux_Error( ERROR_INFO, "particle parameters have NOT been initialized (FaLv %d, FaPID %d, NPar_Copy %d, v %d) !!\n",
                          FaLv, FaPID_Match, amr->patch[0][FaLv][FaPID_Match]->NPar_Copy, v );
         }
#        endif

         for (int v=0; v<NParAtt; v++)
            amr->patch[0][FaLv][FaPID_Match]->ParMassPos_Copy[v][p] = *RecvPtr++;

#        ifdef DEBUG_PARTICLE
//       we do not transfer inactive particles
         if ( amr->patch[0][FaLv][FaPID_Match]->ParMassPos_Copy[PAR_MASS][p] < (real)0.0 )
            Aux_Error( ERROR_INFO, "found inactive particle (FaLv %d, FaPID %d, Mass %14.7e, particle %d) !!\n",
                       FaLv, FaPID_Match, amr->patch[0][FaLv][FaPID_Match]->ParMassPos_Copy[PAR_MASS][p], p );

//       check if the received particle lies within the target patch (may not when PredictPos is on)
         if ( !PredictPos )
         {
            const double *EdgeL     = amr->patch[0][FaLv][FaPID_Match]->EdgeL;
            const double *EdgeR     = amr->patch[0][FaLv][FaPID_Match]->EdgeR;
            const real    ParPos[3] = { amr->patch[0][FaLv][FaPID_Match]->ParMassPos_Copy[PAR_POSX][p],
                                        amr->patch[0][FaLv][FaPID_Match]->ParMassPos_Copy[PAR_POSY][p],
                                        amr->patch[0][FaLv][FaPID_Match]->ParMassPos_Copy[PAR_POSZ][p] };

            for (int d=0; d<3; d++)
            {
               if ( ParPos[d] < EdgeL[d]  ||  ParPos[d] >= EdgeR[d] )
                  Aux_Error( ERROR_INFO, "wrong home patch (L/R edge = %13.6e/%13.6e, pos[%d] = %13.6e, particle %d, FaLv %d, FaPID %d) !!\n",
                             EdgeL[d], EdgeR[d], d, ParPos[d], p, FaLv, FaPID_Match );
            }
         }
#        endif // #ifdef DEBUG_PARTICLE
      } // for (int p=NPar_Copy_Old; p<amr->patch[0][FaLv][FaPID_Match]->NPar_Copy; p++)
   } // for (int t=0; t<NRecvPatchTotal; t++)


// 3-5. free memory
   delete [] RecvBuf_NPatchEachRank;
   delete [] RecvBuf_NParEachPatch;
   delete [] RecvBuf_LBIdxEachPatch;
   delete [] RecvBuf_LBIdxEachPatch_IdxTable;
   delete [] Match_LBIdxEachPatch;
   delete [] FaPIDList;



// 4. add particles temporarily residing in this patch to the ParMassPos_Copy array
   if ( !JustCountNPar )
   for (int FaPID=0; FaPID<amr->NPatchComma[FaLv][1]; FaPID++)
   {
      if ( amr->patch[0][FaLv][FaPID]->son != -1  &&  amr->patch[0][FaLv][FaPID]->NPar > 0 )
      {
         for (int p=0; p<amr->patch[0][FaLv][FaPID]->NPar; p++)
         {
            const long ParID = amr->patch[0][FaLv][FaPID]->ParList[p];
            const int  idx   = amr->patch[0][FaLv][FaPID]->NPar_Copy + p;

//          4-1. check if this particle is indeed waiting for the velocity correction (i.e., ParTime = -dt_half < 0.0 for KDK)
#           ifdef DEBUG_PARTICLE
            if ( amr->Par->Integ == PAR_INTEG_KDK  &&  amr->Par->Time[ParID] >= (real)0.0 )
               Aux_Error( ERROR_INFO, "This particle shouldn't be here (FaLv %d, FaPID %d, ParID %ld, ParTime %21.14e) !!\n",
                          FaLv, FaPID, ParID, amr->Par->Time[ParID] );
#           endif

//          4-2. add particle data
//          --> no need for position prediction here since these particles are all waiting for velocity correction
//          and should already be synchronized with TargetTime
            amr->patch[0][FaLv][FaPID]->ParMassPos_Copy[PAR_MASS][idx] = amr->Par->Attribute[PAR_MASS][ParID];
            amr->patch[0][FaLv][FaPID]->ParMassPos_Copy[PAR_POSX][idx] = amr->Par->Attribute[PAR_POSX][ParID];
            amr->patch[0][FaLv][FaPID]->ParMassPos_Copy[PAR_POSY][idx] = amr->Par->Attribute[PAR_POSY][ParID];
            amr->patch[0][FaLv][FaPID]->ParMassPos_Copy[PAR_POSZ][idx] = amr->Par->Attribute[PAR_POSZ][ParID];
         } // for (int p=0; p<amr->patch[0][FaLv][FaPID]->NPar; p++)

//       4-3. update NPar_Copy
         amr->patch[0][FaLv][FaPID]->NPar_Copy += amr->patch[0][FaLv][FaPID]->NPar;
      } // if ( amr->patch[0][FaLv][FaPID]->son != -1  &&  amr->patch[0][FaLv][FaPID]->NPar > 0 )
   } // for (int FaPID=0; FaPID<amr->NPatchComma[FaLv][1]; FaPID++)



// 5. check if we do collect all particles at levels >= FaLv
#  ifdef DEBUG_PARTICLE
   long NParLocal_Get=0, NParLocal_Check=0, NParAllRank_Get, NParAllRank_Check;

   for (int FaPID=0; FaPID<amr->NPatchComma[FaLv][1]; FaPID++)
      NParLocal_Get += ( amr->patch[0][FaLv][FaPID]->son == -1 ) ? amr->patch[0][FaLv][FaPID]->NPar :
                                                                   amr->patch[0][FaLv][FaPID]->NPar_Copy;

   for (int lv=FaLv; lv<=MAX_LEVEL; lv++)    NParLocal_Check += amr->Par->NPar_Lv[lv];

   MPI_Reduce( &NParLocal_Get,   &NParAllRank_Get,   1, MPI_LONG, MPI_SUM, 0, MPI_COMM_WORLD );
   MPI_Reduce( &NParLocal_Check, &NParAllRank_Check, 1, MPI_LONG, MPI_SUM, 0, MPI_COMM_WORLD );

   if ( MPI_Rank == 0 )
   {
      if ( NParAllRank_Get != NParAllRank_Check )
         Aux_Error( ERROR_INFO, "Total number of active particles >= level %d (%ld) != expected (%ld) !!\n",
                    FaLv, NParAllRank_Get, NParAllRank_Check );
   }
#  endif // #ifdef DEBUG_PARTICLE



// 6. collect particles for buffer patches
// 6-1. sibling-buffer patches at FaLv
   if ( SibBufPatch  &&  !JustCountNPar )
   {
      char Timer_Comment[50];
      sprintf( Timer_Comment, "%3d %15s", FaLv, "Par_Coll2Sib" );

      Par_LB_CollectParticleFromRealPatch(
         FaLv,
         amr->Par->R2B_Buff_NPatchTotal[FaLv][0], amr->Par->R2B_Buff_PIDList[FaLv][0], amr->Par->R2B_Buff_NPatchEachRank[FaLv][0],
         amr->Par->R2B_Real_NPatchTotal[FaLv][0], amr->Par->R2B_Real_PIDList[FaLv][0], amr->Par->R2B_Real_NPatchEachRank[FaLv][0],
         PredictPos, TargetTime, Timer[1], Timer_Comment );
   }

// 6-2. father-sibling-buffer patches at FaLv-1
   if ( FaSibBufPatch  &&  FaLv > 0  &&  !JustCountNPar )
   {
      char Timer_Comment[50];
      sprintf( Timer_Comment, "%3d %15s", FaLv, "Par_Coll2FaSib" );

      Par_LB_CollectParticleFromRealPatch(
         FaLv-1,
         amr->Par->R2B_Buff_NPatchTotal[FaLv][1], amr->Par->R2B_Buff_PIDList[FaLv][1], amr->Par->R2B_Buff_NPatchEachRank[FaLv][1],
         amr->Par->R2B_Real_NPatchTotal[FaLv][1], amr->Par->R2B_Real_PIDList[FaLv][1], amr->Par->R2B_Real_NPatchEachRank[FaLv][1],
         PredictPos, TargetTime, Timer[2], Timer_Comment );
   }

} // FUNCTION : Par_LB_CollectParticle2OneLevel



//-------------------------------------------------------------------------------------------------------
// Function    :  Par_LB_CollectParticle2OneLevel_FreeMemory
// Description :  Release the memory allocated by Par_LB_CollectParticle2OneLevel
//
// Note        :  1. Invoded by Par_CollectParticle2OneLevel_FreeMemory
//
// Parameter   :  lv            : Target refinement level
//                SibBufPatch   : true --> Release memory for sibling-buffer patches at lv as well
//                FaSibBufPatch : true --> Release memory for father-sibling-buffer patches at lv-1 as well
//                                         (do nothing if lv==0)
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Par_LB_CollectParticle2OneLevel_FreeMemory( const int lv, const bool SibBufPatch, const bool FaSibBufPatch )
{

   const int NParAtt = 4;


// 1. real patches at lv
   for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
   {
      for (int v=0; v<NParAtt; v++)
      {
         if ( amr->patch[0][lv][PID]->ParMassPos_Copy[v] != NULL )
         {
            delete [] amr->patch[0][lv][PID]->ParMassPos_Copy[v];
            amr->patch[0][lv][PID]->ParMassPos_Copy[v] = NULL;
         }
      }

//    -1 : indicating that NPar_Copy is not calculated yet
      amr->patch[0][lv][PID]->NPar_Copy = -1;
   }


// 2. sibling-buffer patches at lv
   if ( SibBufPatch )
   for (int p=0; p<amr->Par->R2B_Buff_NPatchTotal[lv][0]; p++)
   {
      const int PID = amr->Par->R2B_Buff_PIDList[lv][0][p];

      for (int v=0; v<NParAtt; v++)
      {
         if ( amr->patch[0][lv][PID]->ParMassPos_Copy[v] != NULL )
         {
            delete [] amr->patch[0][lv][PID]->ParMassPos_Copy[v];
            amr->patch[0][lv][PID]->ParMassPos_Copy[v] = NULL;
         }
      }

//    -1 : indicating that NPar_Copy is not calculated yet
      amr->patch[0][lv][PID]->NPar_Copy = -1;
   }


// 3. father-sibling-buffer patches at lv-1
   const int FaLv = lv - 1;

   if ( FaSibBufPatch  &&  FaLv >= 0 )
   for (int p=0; p<amr->Par->R2B_Buff_NPatchTotal[lv][1]; p++)
   {
      const int FaPID = amr->Par->R2B_Buff_PIDList[lv][1][p];

      for (int v=0; v<NParAtt; v++)
      {
         if ( amr->patch[0][FaLv][FaPID]->ParMassPos_Copy[v] != NULL )
         {
            delete [] amr->patch[0][FaLv][FaPID]->ParMassPos_Copy[v];
            amr->patch[0][FaLv][FaPID]->ParMassPos_Copy[v] = NULL;
         }
      }

//    -1 : indicating that NPar_Copy is not calculated yet
      amr->patch[0][FaLv][FaPID]->NPar_Copy = -1;
   }


// check: if we do everthing correctly, no patches (either real or buffer patches) at lv and lv-1
//        should have ParMassPos_Copy allocated
#  ifdef DEBUG_PARTICLE
   for (int TLv=lv; (TLv>=lv-1 && TLv>=0); TLv--)
   {
//    loop over all real and buffer patches
      for (int PID=0; PID<amr->num[TLv]; PID++)
      {
         for (int v=0; v<NParAtt; v++)
         if ( amr->patch[0][TLv][PID]->ParMassPos_Copy[v] != NULL )
            Aux_Error( ERROR_INFO, "lv %d, PID %d, v %d, ParMassPos_Copy != NULL !!\n",
                       TLv, PID, v );

         if ( amr->patch[0][TLv][PID]->NPar_Copy != -1 )
            Aux_Error( ERROR_INFO, "lv %d, PID %d, NPar_Copy = %d != -1 !!\n",
                       TLv, PID, amr->patch[0][TLv][PID]->NPar_Copy );
      }
   }
#  endif // #ifdef DEBUG_PARTICLE

} // FUNCTION : Par_LB_CollectParticle2OneLevel_FreeMemory



#endif // #if ( defined PARTICLE  &&  defined LOAD_BALANCE )
