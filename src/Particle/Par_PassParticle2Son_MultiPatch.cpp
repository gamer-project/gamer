#include "GAMER.h"

#ifdef PARTICLE



#ifdef TIMING
extern Timer_t *Timer_Par_MPI[NLEVEL][6];
#endif




//-------------------------------------------------------------------------------------------------------
// Function    :  Par_PassParticle2Son_MultiPatch
// Description :  Send particles from real father patches at FaLv to their real son patches at FaLv+1
//
// Note        :  1. For father patches with sons at home, it simply invokes Par_PassParticle2Son_SinlgePatch().
//                   For father patches with sons living abroad, it first sends particles to the corresponding
//                   father-buffer patches by calling Par_LB_ExchangeParticleBetweenPatch(), and then pass these
//                   particles to the real son patches by calling Par_PassParticle2Son_SinglePatch().
//                2. Invoked by
//                   (1) EvolveLevel()
//                       a. Use "Mode = PAR_PASS2SON_EVOLVE"
//                       b. Invoked after the velocity correction in KDK
//                   (2) Flu_CorrAfterAllSync()
//                       a. Use "Mode = PAR_PASS2SON_EVOLVE"
//                   (3) LB_Refine()
//                       a. Use "Mode = PAR_PASS2SON_GENERAL"
//                       b. After allocating new real patches at FaLv+1, their particles may still reside in
//                          **real father** patches temporarily if sons and fathers are NOT in the same rank.
//                          This function send these particles back to their real son patches in other ranks.
//                   (4) Par_AddParticleAfterInit()
//                       a. Use "Mode = PAR_PASS2SON_GENERAL"
//                3. This function does NOT check **buffer** patches at FaLv
//
// Parameter   :  FaLv          : Target father level
//                Mode          : PAR_PASS2SON_GENERAL
//                                --> Check patches in PIDList[], whose sons can be either home or living abroad
//                                PAR_PASS2SON_EVOLVE
//                                --> For non-leaf real patches at FaLv **with sons at home**
//                                    --> Check all of them
//                                    For non-leaf real patches at FaLv **with sons living abroad**
//                                    --> Only check patches adjacent to coarse-fine boundaries
//                                    --> Use the list amr->Par->F2S_* for that
//                                --> No need to set NFaPatch and FaPIDList[] since the target patches have been
//                                    pre-calculated and stored in amr->Par->F2S_*
//                TimingSendPar : Measure the elapsed time of Par_LB_SendParticleData(), which is
//                                called by Par_LB_ExchangeParticleBetweenPatch()
//                                --> LOAD_BALANCE only
//                NFaPatch      : Total number of real patches at FaLv for sending particles
//                                --> Useless for PAR_PASS2SON_EVOLVE
//                FaPIDList     : Patch indices at FaLv for sending particles
//                                --> They must be real instead of buffer patches
//                                --> Useless for PAR_PASS2SON_EVOLVE
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Par_PassParticle2Son_MultiPatch( const int FaLv, const ParPass2Son_t Mode, const bool TimingSendPar,
                                      const int NFaPatch, const int *FaPIDList )
{

   const int SonLv = FaLv + 1;


// nothing to do if there is no patch at FaLv+1
   if ( FaLv == TOP_LEVEL  ||  NPatchTotal[SonLv] == 0 )    return;


#  ifdef DEBUG_PARTICLE
   if ( Mode != PAR_PASS2SON_EVOLVE  &&  Mode != PAR_PASS2SON_GENERAL )
      Aux_Error( ERROR_INFO, "unsupported mode = %d !!\n", Mode );

   if ( Mode == PAR_PASS2SON_GENERAL )
   {
      if ( NFaPatch > 0  &&  FaPIDList == NULL )
         Aux_Error( ERROR_INFO, "FaPIDList == NULL for NFaPatch = %d > 0 !!\n", NFaPatch );

      for (int t=0; t<NFaPatch; t++)
      {
         const int FaPID = FaPIDList[t];

         if ( FaPID < 0  ||  FaPID >= amr->NPatchComma[FaLv][1] )
            Aux_Error( ERROR_INFO, "invalid FaPID %d (FaLv %d, t %d, NReal %d) !!\n",
                       FaPID, FaLv, t, amr->NPatchComma[FaLv][1] );
      }
   } // if ( Mode == PAR_PASS2SON_GENERAL )
#  endif


// 1. pass particles from fathers to sons if they are in the same rank
   const int loop_size = ( Mode == PAR_PASS2SON_EVOLVE ) ? amr->NPatchComma[FaLv][1] : NFaPatch;

   for (int t=0; t<loop_size; t++)
   {
      const int FaPID = ( Mode == PAR_PASS2SON_EVOLVE ) ? t : FaPIDList[t];

      if ( amr->patch[0][FaLv][FaPID]->son >= 0 )  Par_PassParticle2Son_SinglePatch( FaLv, FaPID );
   }


// 3. pass particles from fathers to sons if they are NOT in the same rank
#  ifdef LOAD_BALANCE

// 3-1. set timer
   Timer_t *Timer = NULL;
   char Timer_Comment[20];

#  ifdef TIMING
   if ( TimingSendPar )
   {
      Timer = Timer_Par_MPI[FaLv][2];
      sprintf( Timer_Comment, "%3d %15s", FaLv, "Par_2Son" );
   }
#  endif


// 3-2. find the corresponding father-buffer patches at FaLv to receive particles
//      --> used by Par_LB_ExchangeParticleBetweenPatch()
   int Send_NPatchTotal, *Send_PIDList=NULL, *Send_NPatchEachRank=NULL;
   int Recv_NPatchTotal, *Recv_PIDList=NULL, *Recv_NPatchEachRank=NULL;

// 3-2-a. for PAR_PASS2SON_EVOLVE, just link to the pre-computed info
   if ( Mode == PAR_PASS2SON_EVOLVE )
   {
      Send_NPatchTotal    = amr->Par->F2S_Send_NPatchTotal   [FaLv];
      Send_PIDList        = amr->Par->F2S_Send_PIDList       [FaLv];
      Send_NPatchEachRank = amr->Par->F2S_Send_NPatchEachRank[FaLv];
      Recv_NPatchTotal    = amr->Par->F2S_Recv_NPatchTotal   [FaLv];
      Recv_PIDList        = amr->Par->F2S_Recv_PIDList       [FaLv];
      Recv_NPatchEachRank = amr->Par->F2S_Recv_NPatchEachRank[FaLv];
   }


// 3-2-b. for PAR_PASS2SON_GENERAL, we need to compute the required info
   else
   {
//    b-1. get the number of father patches with sons not home and with particles
      Send_NPatchTotal = 0;
      Send_PIDList     = new int [NFaPatch];    // allocate with the maximum number possible

      for (int t=0; t<NFaPatch; t++)
      {
         const int FaPID = FaPIDList[t];

         if ( amr->patch[0][FaLv][FaPID]->son < -1  &&  amr->patch[0][FaLv][FaPID]->NPar > 0 )
            Send_PIDList[ Send_NPatchTotal ++ ] = FaPID;
      }


//    b-2. compute the load-balance indices of the son patches corresponding to Send_PIDList[]
//         --> note that patches in Send_PIDList[] are at FaLv while Send_LBIdxList[] should
//             store the indices of son patches at FaLv+1
      long *Send_LBIdxList = new long [Send_NPatchTotal];

//###NOTE: faster version can only be applied to the Hilbert space-filling curve
      for (int t=0; t<Send_NPatchTotal; t++)
      {
         const int FaPID     = Send_PIDList[t];
#        if ( LOAD_BALANCE == HILBERT )
         const long SonLBIdx = 8*amr->patch[0][FaLv][FaPID]->LB_Idx;    // faster, LB_Idx of one of the eight sons
#        else
         const long SonLBIdx = LB_Corner2Index( SonLv, amr->patch[0][FaLv][FaPID]->corner, CHECK_ON );   // LB_Idx of son 0
#        endif

         Send_LBIdxList[t] = SonLBIdx;
      }


//    b-3. find the real son patches corresponding to Send_LBIdxList[]
      const bool UseInputLBIdx_Yes = true;

      Send_NPatchEachRank = new int [MPI_NRank];
      Recv_NPatchEachRank = new int [MPI_NRank];

      Par_LB_MapBuffer2RealPatch(
         SonLv,
         Send_NPatchTotal, Send_PIDList, Send_NPatchEachRank,
         Recv_NPatchTotal, Recv_PIDList, Recv_NPatchEachRank,
         UseInputLBIdx_Yes, Send_LBIdxList );

#     ifdef DEBUG_PARTICLE
//    no particles should be exchanged within the same rank
      if ( Send_NPatchEachRank[MPI_Rank] != 0 )
         Aux_Error( ERROR_INFO, "%d patches are sent to itself !!\n", Send_NPatchEachRank[MPI_Rank] );

      if ( Recv_NPatchEachRank[MPI_Rank] != 0 )
         Aux_Error( ERROR_INFO, "%d patches are received from itself !!\n", Recv_NPatchEachRank[MPI_Rank] );
#     endif


//    b-4. record the father-buffer patches to receive particles
//         --> overwrite the original Recv_PIDList[]
      for (int t=0; t<Recv_NPatchTotal; t++)
      {
         const int SonPID = Recv_PIDList[t];

//       check: son patches must be real patches
#        ifdef DEBUG_PARTICLE
         if ( SonPID < 0  ||  SonPID >= amr->NPatchComma[SonLv][1] )
            Aux_Error( ERROR_INFO, "invalid SonPID = %d (FaLv %d, NReal %d) !!\n",
                       SonPID, FaLv, amr->NPatchComma[SonLv][1] );
#        endif

         const int FaPID = amr->patch[0][SonLv][SonPID]->father;

//       check: father patches must be buffer patches
#        ifdef DEBUG_PARTICLE
         if ( FaPID < amr->NPatchComma[FaLv][1] )
            Aux_Error( ERROR_INFO, "this is NOT a buffer patch (FaLv %d, SonPID %d, FaPID %d, NReal %d) !!\n",
                       FaLv, SonPID, FaPID, amr->NPatchComma[FaLv][1] );
#        endif

         Recv_PIDList[t] = FaPID;
      } // for (int t=0; t<Recv_NPatchTotal; t++)

      delete [] Send_LBIdxList;
   } // if ( Mode == PAR_PASS2SON_EVOLVE ) ... else ...


// 3.3. collect particles from other ranks to the father-buffer patches in this rank
   Par_LB_ExchangeParticleBetweenPatch(
      FaLv,
      Send_NPatchTotal, Send_PIDList, Send_NPatchEachRank,
      Recv_NPatchTotal, Recv_PIDList, Recv_NPatchEachRank,
      Timer, Timer_Comment );


// 3.4. pass particles from father-buffer patches to their real son patches in the same rank
   for (int t=0; t<Recv_NPatchTotal; t++)
      Par_PassParticle2Son_SinglePatch( FaLv, Recv_PIDList[t] );


// 3-5. free memory
   if ( Mode == PAR_PASS2SON_GENERAL )
   {
      delete [] Send_PIDList;
      delete [] Recv_PIDList;
      delete [] Send_NPatchEachRank;
      delete [] Recv_NPatchEachRank;
   }

#  endif // #  ifdef LOAD_BALANCE


// 4. check
// neither non-leaf real patches nor buffer patches (either leaf of non-leaf) at FaLv can have particles at this point
#  ifdef DEBUG_PARTICLE
   for (int FaPID=0; FaPID<amr->NPatchComma[FaLv][1]; FaPID++)
      if ( amr->patch[0][FaLv][FaPID]->son != -1  &&  amr->patch[0][FaLv][FaPID]->NPar != 0 )
         Aux_Error( ERROR_INFO, "non-leaf real patch has particles (FaLv %d, FaPID %d, NPar %d) !!\n",
                    FaLv, FaPID, amr->patch[0][FaLv][FaPID]->NPar );

   for (int FaPID=amr->NPatchComma[FaLv][1]; FaPID<amr->NPatchComma[FaLv][3]; FaPID++)
      if ( amr->patch[0][FaLv][FaPID]->NPar != 0 )
         Aux_Error( ERROR_INFO, "buffer patch has particles (FaLv %d, FaPID %d, NPar %d) !!\n",
                    FaLv, FaPID, amr->patch[0][FaLv][FaPID]->NPar );
#  endif

} // FUNCTION : Par_PassParticle2Son_MultiPatch



#endif // #ifdef PARTICLE
