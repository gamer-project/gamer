#include "GAMER.h"

#if ( defined PARTICLE  &&  defined LOAD_BALANCE )




//-------------------------------------------------------------------------------------------------------
// Function    :  Par_LB_CollectParticleFromRealPatch
// Description :  Collect particles for the specified buffer patches from the corresponding real patches
//                in all MPI ranks
//
// Note        :  1. Information of both the target buffer patches "Buff_NPatchTotal, Buff_PIDList, Buff_NPatchEachRank"
//                   and the corresponding real patches "Real_NPatchTotal, Real_PIDList, Real_NPatchEachRank" must be
//                   provided. The information of real patches can be calculated in advance by using Par_LB_MapBuffer2RealPatch()
//                2. All Target patches (those in Buff_PIDList[] and Real_PIDList[]) must be patches at the same level "lv"
//                3. This function is called by Par_LB_CollectParticle2OneLevel()
//                4. ParAtt_Copy[] will be allocated for all target buffer patches with particles in the
//                   corresponding real patches
//                   --> Must be deallocated afterward by calling Par_LB_CollectParticle2OneLevel_FreeMemory()
//
// Parameter   :  lv                  : Target refinement level
//                AttBitIdx           : Bitwise indices of the target particle attributes (e.g., _PAR_MASS | _PAR_VELX)
//                                      --> A user-defined attribute with an integer index AttIntIdx returned by
//                                          AddParticleAttribute() can be converted to a bitwise index by BIDX(AttIntIdx)
//                Buff_NPatchTotal    : Total number of buffer patches in Buff_PIDList
//                Buff_PIDList        : Target buffer patch indices
//                Buff_NPatchEachRank : Number of buffer patches to receive particles from each rank
//                Real_NPatchTotal    : Total number of real patches in Real_PIDList
//                Real_PIDList        : Target real patch indices
//                Real_NPatchEachRank : Number of real patches to send particles to each rank
//                PredictPos          : Predict particle position, which is useful for particle mass assignement
//                                      --> We send particle position **after** prediction so that we don't have to
//                                          send particle velocity
//                TargetTime          : Target time for predicting the particle position
//                Timer               : Timer used by Par_LB_SendParticleData()
//                Timer_Comment       : String used by Par_LB_SendParticleData()
//
// Return      :  NPar_Copy and ParAtt_Copy[] (if NPar_Copy > 0) for all buffer patches specified in Buff_PIDList[]
//-------------------------------------------------------------------------------------------------------
void Par_LB_CollectParticleFromRealPatch( const int lv, const long AttBitIdx,
                                          const int Buff_NPatchTotal, const int *Buff_PIDList, int *Buff_NPatchEachRank,
                                          const int Real_NPatchTotal, const int *Real_PIDList, int *Real_NPatchEachRank,
                                          const bool PredictPos, const double TargetTime,
                                          Timer_t *Timer, const char *Timer_Comment )
{

// nothing to do for levels above MAX_LEVEL
   if ( lv > MAX_LEVEL )  return;



// 0. determine the target particle attributes
// --> assuming _VAR_NAME = 1L<<VAR_NAME (e.g., _PAR_MASS == 1L<<PAR_MASS == BIDX(PAR_MASS))
// --> PosSendIdx[] is used by Par_PredictPos()
   int NAtt=0, AttIntIdx[PAR_NATT_TOTAL], PosSendIdx[3]={-1, -1, -1};

   for (int v=0; v<PAR_NATT_TOTAL; v++)
      if ( AttBitIdx & (1L<<v) )    AttIntIdx[ NAtt ++ ] = v;

   if ( PredictPos )
   {
      for (int v=0; v<NAtt; v++)
      {
         if      ( AttIntIdx[v] == PAR_POSX )   PosSendIdx[0] = v;
         else if ( AttIntIdx[v] == PAR_POSY )   PosSendIdx[1] = v;
         else if ( AttIntIdx[v] == PAR_POSZ )   PosSendIdx[2] = v;
      }

#     ifdef DEBUG_PARTICLE
      for (int d=0; d<3; d++)
         if ( PosSendIdx[d] == -1 )    Aux_Error( ERROR_INFO, "PosSendIdx[%d] == -1 for PredictPos !!\n", d );
#     endif
   }


// check
#  ifdef DEBUG_PARTICLE
   if ( lv < 0  ||  lv >= NLEVEL )     Aux_Error( ERROR_INFO, "incorrect target level (%d) !!\n", lv );

   if ( NAtt == 0  &&  MPI_Rank == 0 )    Aux_Message( stderr, "WARNING : NAtt == 0 !!\n" );

   if      ( Buff_NPatchTotal < 0 )    Aux_Error( ERROR_INFO, "Buff_NPatchTotal = %d < 0 !!\n", Buff_NPatchTotal );
   else if ( Buff_NPatchTotal > 0 )
   {
      if ( Buff_PIDList == NULL )
         Aux_Error( ERROR_INFO, "Buff_PIDList == NULL (Buff_NPatchTotal = %d) !!\n", Buff_NPatchTotal );

      if ( Buff_NPatchEachRank == NULL )
         Aux_Error( ERROR_INFO, "Buff_NPatchEachRank == NULL (Buff_NPatchTotal = %d) !!\n", Buff_NPatchTotal );
   }

   if      ( Real_NPatchTotal < 0 )    Aux_Error( ERROR_INFO, "Real_NPatchTotal = %d < 0 !!\n", Real_NPatchTotal );
   else if ( Real_NPatchTotal > 0 )
   {
      if ( Real_PIDList == NULL )
         Aux_Error( ERROR_INFO, "Real_PIDList == NULL (Real_NPatchTotal = %d) !!\n", Real_NPatchTotal );

      if ( Real_NPatchEachRank == NULL )
         Aux_Error( ERROR_INFO, "Real_NPatchEachRank == NULL (Real_NPatchTotal = %d) !!\n", Real_NPatchTotal );
   }

   for (int t=0; t<Buff_NPatchTotal; t++)
   {
      const int PID = Buff_PIDList[t];

      if ( PID < amr->NPatchComma[lv][1]  ||  PID >= amr->num[lv] )
         Aux_Error( ERROR_INFO, "This is NOT a buffer patch (t %d, PID %d, NReal %d, NTotal %d) !!\n",
                    t, PID, amr->NPatchComma[lv][1], amr->num[lv] );

      if ( amr->patch[0][lv][PID]->NPar > 0 )
         Aux_Error( ERROR_INFO, "lv %d, PID %d, NPar = %d > 0 !!\n", lv, PID, amr->patch[0][lv][PID]->NPar );

      if ( amr->patch[0][lv][PID]->NPar_Copy >= 0 )
         Aux_Error( ERROR_INFO, "lv %d, PID %d, NPar_Copy = %d >= 0 !!\n", lv, PID, amr->patch[0][lv][PID]->NPar_Copy );

      for (int v=0; v<PAR_NATT_TOTAL; v++)
      if ( amr->patch[0][lv][PID]->ParAtt_Copy[v] != NULL )
         Aux_Error( ERROR_INFO, "lv %d, PID %d, NPar_Copy = %d, ParAtt_Copy[%d] != NULL !!\n",
                    lv, PID, amr->patch[0][lv][PID]->NPar_Copy, v );
   } // for (int t=0; t<Buff_NPatchTotal; t++)

   for (int t=0; t<Real_NPatchTotal; t++)
   {
      const int PID = Real_PIDList[t];

      if ( PID < 0  ||  PID >= amr->NPatchComma[lv][1] )
         Aux_Error( ERROR_INFO, "This is NOT a real patch (t %d, PID %d, NReal %d, NTotal %d) !!\n",
                    t, PID, amr->NPatchComma[lv][1], amr->num[lv] );

      if ( amr->patch[0][lv][PID]->NPar < 0 )
         Aux_Error( ERROR_INFO, "lv %d, PID %d, NPar = %d < 0 !!\n", lv, PID, amr->patch[0][lv][PID]->NPar );

//    all leaf real patches should have NPar_Copy == -1
      if ( amr->patch[0][lv][PID]->son == -1  &&  amr->patch[0][lv][PID]->NPar_Copy >= 0 )
         Aux_Error( ERROR_INFO, "lv %d, PID %d has no son, NPar_Copy = %d >= 0 !!\n",
                    lv, PID, amr->patch[0][lv][PID]->NPar_Copy );

//    all non-leaf real patches should have NPar_Copy >= 0 set by Par_LB_CollectParticle2OneLevel()
      if ( amr->patch[0][lv][PID]->son != -1  &&  amr->patch[0][lv][PID]->NPar_Copy < 0 )
         Aux_Error( ERROR_INFO, "lv %d, PID %d, SonPID %d, NPar_Copy = %d < 0 !!\n",
                    lv, PID, amr->patch[0][lv][PID]->son, amr->patch[0][lv][PID]->NPar_Copy );
   }

// check Real_NPatchEachRank and Real_NPatchTotal
   int Real_NPatchEachRank_Check[MPI_NRank], Real_NPatchTotal_Check=0;

   MPI_Alltoall( Buff_NPatchEachRank, 1, MPI_INT, Real_NPatchEachRank_Check, 1, MPI_INT, MPI_COMM_WORLD );

   for (int r=0; r<MPI_NRank; r++)
   if ( Real_NPatchEachRank[r] != Real_NPatchEachRank_Check[r] )
      Aux_Error( ERROR_INFO, "Real_NPatchEachRank[%d] (%d) != expected (%d) !!\n",
                 r, Real_NPatchEachRank[r], Real_NPatchEachRank_Check[r] );

   for (int r=0; r<MPI_NRank; r++)  Real_NPatchTotal_Check += Real_NPatchEachRank_Check[r];

   if ( Real_NPatchTotal != Real_NPatchTotal_Check )
      Aux_Error( ERROR_INFO, "Real_NPatchTotal (%d) != expected (%d) !!\n",
                 Real_NPatchTotal, Real_NPatchTotal_Check );
#  endif // #ifdef DEBUG_PARTICLE


// must NOT call "return" here even if Buff/Real_NPatchTotal==0 since this rank still needs to call Par_LB_SendParticleData()
// if ( Real_NPatchTotal == 0 )   return;
// if ( Buff_NPatchTotal == 0 )   return;


// 1. get the number of particles to be sent
   int  *SendBuf_NParEachPatch = new int  [Real_NPatchTotal];
   long *SendBuf_Offset        = new long [Real_NPatchTotal];

   int  PID, NParThisPatch;
   long NSendParTotal = 0L;

// loop over all target real patches
#  pragma omp parallel for private( PID, NParThisPatch ) reduction( +:NSendParTotal ) schedule( PAR_OMP_SCHED, PAR_OMP_SCHED_CHUNK )
   for (int t=0; t<Real_NPatchTotal; t++)
   {
      PID = Real_PIDList[t];

//    we assume that non-leaf real patches have already collected particles from their descendant patches
      if ( amr->patch[0][lv][PID]->son == -1 )  NParThisPatch = amr->patch[0][lv][PID]->NPar;
      else                                      NParThisPatch = amr->patch[0][lv][PID]->NPar_Copy;

#     ifdef DEBUG_PARTICLE
      if ( NParThisPatch < 0 )
         Aux_Error( ERROR_INFO, "NParThisPatch (%d) has not been calculated (lv %d, PID %d) !!\n",
                    NParThisPatch, lv, PID );
#     endif

      NSendParTotal            += (long)NParThisPatch;
      SendBuf_NParEachPatch[t]  =       NParThisPatch;
   } // for (int t=0; t<Real_NPatchTotal; t++)

// get the array offset of each patch (mainly for the OpenMP parallelization)
   if ( Real_NPatchTotal > 0 )   SendBuf_Offset[0] = 0L;
   for (int t=0; t<Real_NPatchTotal-1; t++)  SendBuf_Offset[t+1] = SendBuf_Offset[t] + long(SendBuf_NParEachPatch[t]*NAtt);


// 2. prepare the particle data to be sent
// reuse the MPI send buffer declared in LB_GetBufferData() for better MPI performance
   real_par  *SendBuf_ParDataEachPatch = (real_par *)LB_GetBufferData_MemAllocate_Send( NSendParTotal*(long)NAtt*sizeof(real_par) );

   real_par  *SendPtr     = NULL;
   long      *ParList     = NULL;
   real_par **ParAtt_Copy = NULL;
   long       ParID;

#  pragma omp parallel for private( PID, NParThisPatch, SendPtr, ParList, ParAtt_Copy, ParID ) \
                           schedule( PAR_OMP_SCHED, PAR_OMP_SCHED_CHUNK )
   for (int t=0; t<Real_NPatchTotal; t++)
   {
      PID           = Real_PIDList         [t];
      NParThisPatch = SendBuf_NParEachPatch[t];
      SendPtr       = SendBuf_ParDataEachPatch + SendBuf_Offset[t];

//    skip patches with no particles
      if ( NParThisPatch == 0 )  continue;

      if ( amr->patch[0][lv][PID]->son == -1 )
      {
         ParList = amr->patch[0][lv][PID]->ParList;

#        ifdef DEBUG_PARTICLE
         if ( ParList == NULL )
            Aux_Error( ERROR_INFO, "ParList == NULL for NParThisPatch (%d) > 0 (lv %d, PID %d) !!\n",
                       NParThisPatch, lv, PID );
#        endif

         for (int p=0; p<NParThisPatch; p++)
         {
            ParID = ParList[p];

            for (int v=0; v<NAtt; v++)    SendPtr[v] = amr->Par->Attribute[ AttIntIdx[v] ][ParID];

//          predict particle position to TargetTime
//          --> note that we need to skip particles waiting for velocity correction since these are leaf real patches
//              which may have particles just been updated
//              --> done in Par_PredictPos()
//          --> also note that we don't have to worry about the periodic BC here (in other words, Pos can lie outside the box)
            if ( PredictPos )    Par_PredictPos( 1, &ParID, SendPtr+PosSendIdx[0], SendPtr+PosSendIdx[1], SendPtr+PosSendIdx[2],
                                                 TargetTime );

            SendPtr += NAtt;
         } // for (int p=0; p<NParThisPatch; p++)
      } // if ( amr->patch[0][lv][PID]->son == -1 )

      else
      {
         ParAtt_Copy = amr->patch[0][lv][PID]->ParAtt_Copy;

#        ifdef DEBUG_PARTICLE
         for (int v=0; v<NAtt; v++)
            if ( ParAtt_Copy[ AttIntIdx[v] ] == NULL )
               Aux_Error( ERROR_INFO, "ParAtt_Copy[%d] == NULL for NParThisPatch (%d) > 0 (lv %d, PID %d) !!\n",
                          AttIntIdx[v], NParThisPatch, lv, PID );
#        endif

         for (int p=0; p<NParThisPatch; p++)
         {
//          note that the particle position should have already been predicted to TargetTime
//          by Par_LB_CollectParticle2OneLevel()
            for (int v=0; v<NAtt; v++)    SendPtr[v] = ParAtt_Copy[ AttIntIdx[v] ][p];

            SendPtr += NAtt;
         }
      } // if ( amr->patch[0][lv][PID]->son == -1 ) ... else ...
   } // for (int t=0; t<Real_NPatchTotal; t++)


// 3. send the number of particles and their attributes
   const bool Exchange_NPatchEachRank_No   = false;
   const bool Exchange_LBIdxEachRank_No    = false;
   const bool Exchange_ParDataEachRank_Yes = true;

   int      *SendBuf_NPatchEachRank   = Real_NPatchEachRank;
   int      *RecvBuf_NPatchEachRank   = Buff_NPatchEachRank;
   int      *RecvBuf_NParEachPatch    = NULL;   // will be allocated by Par_LB_SendParticleData and must be free'd later
   real_par *RecvBuf_ParDataEachPatch = NULL;   // a pointer to the MPI recv buffer declared in LB_GetBufferData
                                                // --> don't have to be free'd here

   long     *SendBuf_LBIdxEachRank    = NULL;   // useless and does not need to be allocated
   long     *RecvBuf_LBIdxEachRank    = NULL;   // useless and will not be allocated by Par_LB_SendParticleData

   int       NRecvPatchTotal;                   // returned from Par_LB_SendParticleData
   long      NRecvParTotal;                     // returned from Par_LB_SendParticleData

// note that we don't exchange NPatchEachRank (which is already known) and LBIdxEachRank (which is useless here)
   Par_LB_SendParticleData(
      NAtt,
      SendBuf_NPatchEachRank, SendBuf_NParEachPatch, SendBuf_LBIdxEachRank, SendBuf_ParDataEachPatch, NSendParTotal,
      RecvBuf_NPatchEachRank, RecvBuf_NParEachPatch, RecvBuf_LBIdxEachRank, RecvBuf_ParDataEachPatch,
      NRecvPatchTotal, NRecvParTotal, Exchange_NPatchEachRank_No, Exchange_LBIdxEachRank_No, Exchange_ParDataEachRank_Yes,
      Timer, Timer_Comment );

#  ifdef DEBUG_PARTICLE
   if ( NRecvPatchTotal != Buff_NPatchTotal )
      Aux_Error( ERROR_INFO, "Total number of received patches (%d) != expected (%d) !!\n",
                 NRecvPatchTotal, Buff_NPatchTotal );
#  endif

// free the send buffer in advance to save memory
   delete [] SendBuf_NParEachPatch;
   delete [] SendBuf_Offset;


// 4. store the received particle data to each patch
   const real_par *RecvPtr = NULL;

// 4-0. get the array offset of each patch (mainly for the OpenMP parallelization)
   long *RecvBuf_Offset = new long [Buff_NPatchTotal];
   if ( Buff_NPatchTotal > 0 )   RecvBuf_Offset[0] = 0L;
   for (int t=0; t<Buff_NPatchTotal-1; t++)  RecvBuf_Offset[t+1] = RecvBuf_Offset[t] + long(RecvBuf_NParEachPatch[t]*NAtt);

#  pragma omp parallel for private( PID, NParThisPatch, RecvPtr ) schedule( PAR_OMP_SCHED, PAR_OMP_SCHED_CHUNK )
   for (int t=0; t<Buff_NPatchTotal; t++)
   {
      PID           = Buff_PIDList         [t];
      NParThisPatch = RecvBuf_NParEachPatch[t];
      RecvPtr       = RecvBuf_ParDataEachPatch + RecvBuf_Offset[t];

//    4-1. set the number of particles in this patch
      amr->patch[0][lv][PID]->NPar_Copy = NParThisPatch;

      if ( NParThisPatch > 0 )
      {
//       4-2. allocate ParAtt_Copy[]
         for (int v=0; v<NAtt; v++)
            amr->patch[0][lv][PID]->ParAtt_Copy[ AttIntIdx[v] ] = new real_par [NParThisPatch];

         for (int p=0; p<NParThisPatch; p++)
         {
//          4-3. store the received data
            for (int v=0; v<NAtt; v++)
               amr->patch[0][lv][PID]->ParAtt_Copy[ AttIntIdx[v] ][p] = *RecvPtr++;

//          4-4. check
#           ifdef DEBUG_PARTICLE
//          we do not transfer inactive particles
            if ( AttBitIdx & _PAR_MASS )
            if ( amr->patch[0][lv][PID]->ParAtt_Copy[PAR_MASS][p] < (real_par)0.0 )
               Aux_Error( ERROR_INFO, "found inactive particle (lv %d, PID %d, Mass %14.7e, particle %d) !!\n",
                          lv, PID, amr->patch[0][lv][PID]->ParAtt_Copy[PAR_MASS][p], p );

//          check if the received particle lies within the target patch (may not when PredictPos is on)
            if ( !PredictPos  &&  ( AttBitIdx & _PAR_POSX )  &&  ( AttBitIdx & _PAR_POSY )  &&  ( AttBitIdx & _PAR_POSZ ) )
            {
//             always assume periodic B.C. in this check since we don't allocate buffer patches lying outside
//             the simulation domain for non-periodic B.C.
//             --> we can use EdgeL/R stored in each patch directly since they assume periodicity as well
               const double     *EdgeL     = amr->patch[0][lv][PID]->EdgeL;
               const double     *EdgeR     = amr->patch[0][lv][PID]->EdgeR;
               const real_par    ParPos[3] = { amr->patch[0][lv][PID]->ParAtt_Copy[PAR_POSX][p],
                                               amr->patch[0][lv][PID]->ParAtt_Copy[PAR_POSY][p],
                                               amr->patch[0][lv][PID]->ParAtt_Copy[PAR_POSZ][p] };

               for (int d=0; d<3; d++)
               {
                  if ( ParPos[d] < EdgeL[d]  ||  ParPos[d] >= EdgeR[d] )
                     Aux_Error( ERROR_INFO, "wrong home patch (L/R edge = %13.6e/%13.6e, pos[%d] = %13.6e, particle %d, lv %d, PID %d) !!\n",
                                EdgeL[d], EdgeR[d], d, ParPos[d], p, lv, PID );
               }
            } // if ( !PredictPos )
#           endif // #ifdef DEBUG_PARTICLE

         } // for (int p=0; p<NParThisPatch; p++)
      } // if ( NParThisPatch > 0 )
   } // for (int t=0; t<Real_NPatchTotal; t++)


// 5. free memory
   delete [] RecvBuf_NParEachPatch;
   delete [] RecvBuf_Offset;

} // FUNCTION : Par_LB_CollectParticleFromRealPatch



#endif // #if ( defined PARTICLE  &&  defined LOAD_BALANCE )
