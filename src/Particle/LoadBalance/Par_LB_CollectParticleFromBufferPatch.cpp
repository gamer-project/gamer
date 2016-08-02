#include "Copyright.h"
#include "GAMER.h"

#if ( defined PARTICLE  &&  defined LOAD_BALANCE )




//-------------------------------------------------------------------------------------------------------
// Function    :  Par_LB_CollectParticleFromBufferPatch
// Description :  Collect particles for the specified buffer patches from the corresponding real patches
//                in all MPI ranks
//
// Note        :  1. Information of both the target buffer patches "Buff_NPatchTotal, Buff_PIDList, Buff_NPatchEachRank"
//                   and the corresponding real patches "Real_NPatchTotal, Real_PIDList, Real_NPatchEachRank" must be
//                   provided. The information of real patches can be calculated in advance by using "Par_LB_MapBuffer2RealPatch"
//                2. All Target patches (those in Buff_PIDList and Real_PIDList) must be patches at the same level "lv"
//                3. Currently this function exchange all particles attributes (both ParVar and Passive arrays)
//                   --> But it can be generalized to work with arbitrary particle attributes
//                4. This function is called by Par_PassParticle2Sibling
//                5. Note that after calling this function, some particles may reside in **non-leaf** real patches
//                   --> They will be sent again to leaf real patches after the velocity correction operation
//                       (by the function Par_PassParticle2Son_AllPatch)
//
// Parameter   :  lv                   : Target refinement level
//                Buff_NPatchTotal     : Total number of buffer patches in Buff_PIDList
//                Buff_PIDList         : Target buffer patch indices
//                Buff_NPatchEachRank  : Number of buffer patches to receive particles from each rank
//                Real_NPatchTotal     : Total number of real patches in Real_PIDList
//                Real_PIDList         : Target real patch indices
//                Real_NPatchEachRank  : Number of real patches to send particles to each rank
//
// Return      :  New particles will be added to the particle repository of this rank and linked to the
//                target real patches (which may NOT be leaf patches)
//-------------------------------------------------------------------------------------------------------
void Par_LB_CollectParticleFromBufferPatch( const int lv,
                                            const int Buff_NPatchTotal, const int *Buff_PIDList, int *Buff_NPatchEachRank,
                                            const int Real_NPatchTotal, const int *Real_PIDList, int *Real_NPatchEachRank )
{

// nothing to do for levels above MAX_LEVEL
   if ( lv > MAX_LEVEL )  return;


// check
#  ifdef DEBUG_PARTICLE
   if ( lv < 0  ||  lv >= NLEVEL )     Aux_Error( ERROR_INFO, "incorrect target level (%d) !!\n", lv );

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

      if ( amr->patch[0][lv][PID]->NPar < 0 )
         Aux_Error( ERROR_INFO, "lv %d, PID %d, NPar = %d < 0 !!\n", lv, PID, amr->patch[0][lv][PID]->NPar );

      if ( amr->patch[0][lv][PID]->NPar_Away >= 0 )
         Aux_Error( ERROR_INFO, "lv %d, PID %d, NPar_Away = %d >= 0 !!\n", lv, PID, amr->patch[0][lv][PID]->NPar_Away );

      for (int v=0; v<4; v++)
      if ( amr->patch[0][lv][PID]->ParMassPos_Away[v] != NULL )
         Aux_Error( ERROR_INFO, "lv %d, PID %d, NPar_Away = %d, ParMassPos_Away[%d] != NULL !!\n",
                    lv, PID, amr->patch[0][lv][PID]->NPar_Away, v );
   } // for (int t=0; t<Buff_NPatchTotal; t++)

   for (int t=0; t<Real_NPatchTotal; t++)
   {
      const int PID = Real_PIDList[t];

      if ( PID < 0  ||  PID >= amr->NPatchComma[lv][1] )
         Aux_Error( ERROR_INFO, "This is NOT a real patch (t %d, PID %d, NReal %d, NTotal %d) !!\n",
                    t, PID, amr->NPatchComma[lv][1], amr->num[lv] );

      if ( amr->patch[0][lv][PID]->NPar < 0 )
         Aux_Error( ERROR_INFO, "lv %d, PID %d, NPar = %d < 0 !!\n", lv, PID, amr->patch[0][lv][PID]->NPar );

//    no real patches should have NPar_Away != -1 at this point
      if ( amr->patch[0][lv][PID]->NPar_Away >= 0 )
         Aux_Error( ERROR_INFO, "lv %d, PID %d, NPar_Away = %d >= 0 !!\n",
                    lv, PID, amr->patch[0][lv][PID]->NPar_Away );
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


// must NOT call "return" here even if Buff/Real_NPatchTotal==0 since this rank still needs to call Par_LB_ExchangeParticle
// if ( Real_NPatchTotal == 0 )   return;
// if ( Buff_NPatchTotal == 0 )   return;


// 1. get the number of particles to be sent
   int *SendBuf_NParEachPatch = new int [Buff_NPatchTotal];

   int PID, NParThisPatch, NSendParTotal = 0;

// loop over all target buffer patches
   for (int t=0; t<Buff_NPatchTotal; t++)
   {
      PID = Buff_PIDList[t];

//    we assume that some buffer patches just received particles from nearby real patches in the same rank
//    (in the function Par_PassParticle2Sibling)
      NParThisPatch = amr->patch[0][lv][PID]->NPar;

#     ifdef DEBUG_PARTICLE
      if ( NParThisPatch < 0 )
         Aux_Error( ERROR_INFO, "NParThisPatch (%d) has not been calculated (lv %d, PID %d) !!\n",
                    NParThisPatch, lv, PID );
#     endif

      NSendParTotal            += NParThisPatch;
      SendBuf_NParEachPatch[t]  = NParThisPatch;
   } // for (int t=0; t<Buff_NPatchTotal; t++)


// 2. prepare the particle data to be sent (and then remove these particles from this rank)
   const bool RemoveAllPar_Yes = true;
   const int  NParVar          = NPAR_VAR + NPAR_PASSIVE;

   real *SendBuf_ParDataEachPatch = new real [ NSendParTotal*NParVar ];

   real *SendPtr = SendBuf_ParDataEachPatch;
   long *ParList = NULL;
   long  ParID;

   for (int t=0; t<Buff_NPatchTotal; t++)
   {
      PID           = Buff_PIDList         [t];
      NParThisPatch = SendBuf_NParEachPatch[t];

//    skip patches with no particles
      if ( NParThisPatch == 0 )  continue;

      ParList = amr->patch[0][lv][PID]->ParList;

#     ifdef DEBUG_PARTICLE
      if ( ParList == NULL )
         Aux_Error( ERROR_INFO, "ParList == NULL for NParThisPatch (%d) > 0 (lv %d, PID %d) !!\n",
                    NParThisPatch, lv, PID );
#     endif

      for (int p=0; p<NParThisPatch; p++)
      {
         ParID = ParList[p];

//       2-1. store particle data into the MPI send buffer
         for (int v=0; v<NPAR_VAR;     v++)  *SendPtr++ = amr->Par->ParVar [v][ParID];
         for (int v=0; v<NPAR_PASSIVE; v++)  *SendPtr++ = amr->Par->Passive[v][ParID];

//       2-2. remove this particle from the particle repository of this rank
//       (don't have to update the average density since particles are just sent to other ranks)
         amr->Par->RemoveOneParticle( ParID, PAR_INACTIVE_MPI, NULL, NULL_REAL );
      }

//    2-3. remove all particles in this buffer patch
      amr->patch[0][lv][PID]->RemoveParticle( NULL_INT, NULL, &amr->Par->NPar_Lv[lv], RemoveAllPar_Yes );
   } // for (int t=0; t<Buff_NPatchTotal; t++)


// 3. send the number of particles and their attributes
   const bool Exchange_NPatchEachRank_No = false;
   const bool Exchange_LBIdxEachRank_No  = false;

   int  *SendBuf_NPatchEachRank   = Buff_NPatchEachRank;
   int  *RecvBuf_NPatchEachRank   = Real_NPatchEachRank;
   int  *RecvBuf_NParEachPatch    = NULL;   // will be allocated by Par_LB_ExchangeParticle and must be free'd later
   real *RecvBuf_ParDataEachPatch = NULL;   // will be allocated by Par_LB_ExchangeParticle and must be free'd later

   long *SendBuf_LBIdxEachRank    = NULL;   // useless and does not need to be allocated
   long *RecvBuf_LBIdxEachRank    = NULL;   // useless and will not be allocated by Par_LB_ExchangeParticle

   int NRecvPatchTotal, NRecvParTotal;  // returned from Par_LB_ExchangeParticle

// note that we don't exchange NPatchEachRank (which is already known) and LBIdxEachRank (which is useless here)
   Par_LB_ExchangeParticle(
      NParVar,
      SendBuf_NPatchEachRank, SendBuf_NParEachPatch, SendBuf_LBIdxEachRank, SendBuf_ParDataEachPatch,
      RecvBuf_NPatchEachRank, RecvBuf_NParEachPatch, RecvBuf_LBIdxEachRank, RecvBuf_ParDataEachPatch,
      NRecvPatchTotal, NRecvParTotal, Exchange_NPatchEachRank_No, Exchange_LBIdxEachRank_No );

#  ifdef DEBUG_PARTICLE
   if ( NRecvPatchTotal != Real_NPatchTotal )
      Aux_Error( ERROR_INFO, "Total number of received patches (%d) != expected (%d) !!\n",
                 NRecvPatchTotal, Real_NPatchTotal );
#  endif

// free the send buffer in advance to save memory
   delete [] SendBuf_NParEachPatch;
   delete [] SendBuf_ParDataEachPatch;


// 4. store the received particle data to the particle repository and link to each real patch
   const real *RecvPtr = RecvBuf_ParDataEachPatch;

// 4-1. get the maximum number of particles in one patch
   int   NParThisPatch_Max;
   long *NewParIDList = NULL;

   NParThisPatch_Max = 0;
   for (int t=0; t<Real_NPatchTotal; t++)    NParThisPatch_Max = MAX( NParThisPatch_Max, RecvBuf_NParEachPatch[t] );

   NewParIDList = new long [NParThisPatch_Max];


   for (int t=0; t<Real_NPatchTotal; t++)
   {
      NParThisPatch = RecvBuf_NParEachPatch[t];

//    4-2. add particles to the particle repository
      for (int p=0; p<NParThisPatch; p++)
      {
//       don't have to update the average density since particles are just received from other ranks
         ParID    = amr->Par->AddOneParticle( RecvPtr, RecvPtr+NPAR_VAR, NULL, NULL_REAL );
         RecvPtr += NParVar;

//       store the new particle index
         NewParIDList[p] = ParID;

//       we do not transfer inactive particles
#        ifdef DEBUG_PARTICLE
         if ( amr->Par->ParVar[PAR_MASS][ParID] < (real)0.0 )
            Aux_Error( ERROR_INFO, "found inactive particle (ParID %d, Mass %14.7e) !!\n",
                       ParID, amr->Par->ParVar[PAR_MASS][ParID] );
#        endif
      }

//    4-3. add particles to the real patch
      PID = Real_PIDList[t];

#     ifdef DEBUG_PARTICLE
//    do not set ParPos too early since pointers to the particle repository (e.g., amr->Par->PosX)
//    may change after calling amr->Par->AddOneParticle
      const real *ParPos[3] = { amr->Par->PosX, amr->Par->PosY, amr->Par->PosZ };
      char Comment[100];
      sprintf( Comment, "%s", __FUNCTION__ );

      amr->patch[0][lv][PID]->AddParticle( NParThisPatch, NewParIDList, &amr->Par->NPar_Lv[lv],
                                           ParPos, amr->Par->NPar_AcPlusInac, Comment );
#     else
      amr->patch[0][lv][PID]->AddParticle( NParThisPatch, NewParIDList, &amr->Par->NPar_Lv[lv] );
#     endif
   } // for (int t=0; t<Real_NPatchTotal; t++)


// 5. free memory
   delete [] RecvBuf_NParEachPatch;
   delete [] RecvBuf_ParDataEachPatch;
   delete [] NewParIDList;

} // FUNCTION : Par_LB_CollectParticleFromBufferPatch



#endif // #if ( defined PARTICLE  &&  defined LOAD_BALANCE )
