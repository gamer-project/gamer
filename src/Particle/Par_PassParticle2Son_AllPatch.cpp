#include "GAMER.h"

#ifdef PARTICLE

#ifdef TIMING
extern Timer_t *Timer_Par_MPI[NLEVEL][6];
#endif




//-------------------------------------------------------------------------------------------------------
// Function    :  Par_PassParticle2Son_AllPatch
// Description :  Pass particles from father to sons for all patches at the target level
//
// Note        :  1. It simply invokes Par_PassParticle2Son() for all patches at the target level
//                   (and with sons at home)
//                   --> For patches with sons living abroad, this function will first collect
//                       particles from other ranks to the father-buffer patches in this rank
//                       by calling Par_LB_ExchangeParticleBetweenPatch(), and then pass
//                       these particles to the real son patches in the same rank by calling
//                       Par_PassParticle2Son() again
//                2. It is invoked by EvolveLevel() after the velocity correction in KDK
//
// Parameter   :  FaLv          : Father's refinement level
//                TimingSendPar : Measure the elapsed time of Par_LB_SendParticleData(), which is
//                                called by Par_LB_ExchangeParticleBetweenPatch()
//                                --> LOAD_BALANCE only
//-------------------------------------------------------------------------------------------------------
void Par_PassParticle2Son_AllPatch( const int FaLv, const bool TimingSendPar )
{

// nothing to do if there is no patch at FaLv+1
   if ( FaLv == TOP_LEVEL  ||  NPatchTotal[FaLv+1] == 0 )   return;


// pass particles from fathers to sons if they are in the same rank
   for (int FaPID=0; FaPID<amr->NPatchComma[FaLv][1]; FaPID++)
   {
//    patches with sons in other ranks (those with SonPID<-1) will be dealt with later
      if ( amr->patch[0][FaLv][FaPID]->son >= 0 )  Par_PassParticle2Son( FaLv, FaPID );
   }


// deal with the case when fathers and sons are NOT in the same rank
#  ifdef LOAD_BALANCE
   int FaBufPID;

   Timer_t *Timer = NULL;
   char Timer_Comment[20];
#  ifdef TIMING
   if ( TimingSendPar )
   {
      Timer = Timer_Par_MPI[FaLv][2];
      sprintf( Timer_Comment, "%3d %15s", FaLv, "Par_2Son" );
   }
#  endif

// collect particles from other ranks to the father-buffer patches in this rank
   Par_LB_ExchangeParticleBetweenPatch(
      FaLv,
      amr->Par->F2S_Send_NPatchTotal[FaLv], amr->Par->F2S_Send_PIDList[FaLv], amr->Par->F2S_Send_NPatchEachRank[FaLv],
      amr->Par->F2S_Recv_NPatchTotal[FaLv], amr->Par->F2S_Recv_PIDList[FaLv], amr->Par->F2S_Recv_NPatchEachRank[FaLv],
      Timer, Timer_Comment );

// pass particles from father-buffer patches to their real son patches in the same rank
   for (int t=0; t<amr->Par->F2S_Recv_NPatchTotal[FaLv]; t++)
   {
      FaBufPID = amr->Par->F2S_Recv_PIDList[FaLv][t];

#     ifdef DEBUG_PARTICLE
      if ( FaBufPID < amr->NPatchComma[FaLv][1] )
         Aux_Error( ERROR_INFO, "This is NOT a father-buffer patch (FaLv %d, FaPID %d, NReal %d) !!\n",
                    FaLv, FaBufPID, amr->NPatchComma[FaLv][1] );

      if ( amr->patch[0][FaLv][FaBufPID]->son < 0 )
         Aux_Error( ERROR_INFO, "Father-buffer patch has no son at home (FaLv %d, FaPID %d, SonPID %d) !!\n",
                    FaLv, FaBufPID, amr->patch[0][FaLv][FaBufPID]->son );
#     endif

      Par_PassParticle2Son( FaLv, FaBufPID );
   }

// check: neither non-leaf real patches nor buffer patches (either leaf of non-leaf) at FaLv can have particles at this point
#  ifdef DEBUG_PARTICLE
   for (int FaPID=0; FaPID<amr->NPatchComma[FaLv][1]; FaPID++)
      if ( amr->patch[0][FaLv][FaPID]->son != -1  &&  amr->patch[0][FaLv][FaPID]->NPar != 0 )
         Aux_Error( ERROR_INFO, "Non-leaf real patch has particles (FaLv %d, FaPID %d, NPar %d) !!\n",
                    FaLv, FaPID, amr->patch[0][FaLv][FaPID]->NPar );

   for (int FaPID=amr->NPatchComma[FaLv][1]; FaPID<amr->NPatchComma[FaLv][3]; FaPID++)
      if ( amr->patch[0][FaLv][FaPID]->NPar != 0 )
         Aux_Error( ERROR_INFO, "Buffer patch has particles (FaLv %d, FaPID %d, NPar %d) !!\n",
                    FaLv, FaPID, amr->patch[0][FaLv][FaPID]->NPar );
#  endif
#  endif // #ifdef LOAD_BALANCE

} // FUNCTION : Par_PassParticle2Son_AllPatch



#endif // #ifdef PARTICLE
