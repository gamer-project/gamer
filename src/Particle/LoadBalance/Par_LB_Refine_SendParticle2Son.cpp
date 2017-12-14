#include "GAMER.h"

#if ( defined PARTICLE  &&  defined LOAD_BALANCE )




//-------------------------------------------------------------------------------------------------------
// Function    :  Par_LB_Refine_SendParticle2Son
// Description :  Send particles to leaf real patches at FaLv+1 after grid refinement
//
// Note        :  1. Called by LB_Refine
//                2. After allocating new real patches at FaLv+1, their particles may still reside in **real**
//                   father patches temporarily if sons and fathers are NOT in the same rank. This function
//                   send these particles from real father patches at one rank to their real son patches in
//                   a different rank.
//                   --> Similar to the function "Par_PassParticle2Son_AllPatch"
//                3. The input arrays RefineF2S_Send_PIDList and RefineF2S_Send_LBIdxList must be free'd manually
//                   after calling this function
//
// Parameter   :  FaLv                       : Target father level
//                RefineF2S_Send_NPatchTotal : Total number of real patches at FaLv for sending particles
//                RefineF2S_Send_PIDList     : Patch indices at FaLv for sending particles
//                RefineF2S_Send_LBIdxList   : Load-balance indices of the sons of patches in RefineF2S_Send_PIDList
//                                             --> Used by Par_LB_MapBuffer2RealPatch to find the rank of these sons
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Par_LB_Refine_SendParticle2Son( const int FaLv, const int RefineF2S_Send_NPatchTotal, int *RefineF2S_Send_PIDList,
                                     long *RefineF2S_Send_LBIdxList )
{

// nothing to do if there is no higher level
   if ( FaLv == TOP_LEVEL )   return;


#  ifdef DEBUG_PARTICLE
   if ( RefineF2S_Send_NPatchTotal > 0  &&  RefineF2S_Send_PIDList == NULL )
      Aux_Error( ERROR_INFO, "RefineF2S_Send_PIDList == NULL for RefineF2S_Send_NPatchTotal = %d > 0 !!\n",
                 RefineF2S_Send_NPatchTotal );

   if ( RefineF2S_Send_NPatchTotal > 0  &&  RefineF2S_Send_LBIdxList == NULL )
      Aux_Error( ERROR_INFO, "RefineF2S_Send_LBIdxList == NULL for RefineF2S_Send_NPatchTotal = %d > 0 !!\n",
                 RefineF2S_Send_NPatchTotal );
#  endif


// 1. find the corresponding buffer patches at FaLv to receive particles
   const bool UseInputLBIdx_Yes = true;
   const int  SonLv             = FaLv + 1;

   int  FaPID, SonPID;
   int  RefineF2S_Send_NPatchEachRank[MPI_NRank];
   int  RefineF2S_Recv_NPatchEachRank[MPI_NRank];
   int  RefineF2S_Recv_NPatchTotal;
   int *RefineF2S_Recv_PIDList = NULL;

// 1-1. find the real patches at SonLv corresponding to RefineF2S_Send_LBIdxList
   Par_LB_MapBuffer2RealPatch( SonLv, RefineF2S_Send_NPatchTotal,
                                      RefineF2S_Send_PIDList,
                                      RefineF2S_Send_NPatchEachRank,
                                      RefineF2S_Recv_NPatchTotal,
                                      RefineF2S_Recv_PIDList,
                                      RefineF2S_Recv_NPatchEachRank,
                                      UseInputLBIdx_Yes, RefineF2S_Send_LBIdxList );

#  ifdef DEBUG_PARTICLE
// no particles should be exchanged within the same rank
   if ( RefineF2S_Send_NPatchEachRank[MPI_Rank] != 0 )
      Aux_Error( ERROR_INFO, "%d patches are sent to itself !!\n", RefineF2S_Send_NPatchEachRank[MPI_Rank] );

   if ( RefineF2S_Recv_NPatchEachRank[MPI_Rank] != 0 )
      Aux_Error( ERROR_INFO, "%d patches are received from itself !!\n", RefineF2S_Recv_NPatchEachRank[MPI_Rank] );

// all recv patches at SonLv should have no particles
   for (int t=0; t<RefineF2S_Recv_NPatchTotal; t++)
   {
      SonPID = RefineF2S_Recv_PIDList[t];

      if ( amr->patch[0][SonLv][SonPID]->NPar != 0 )
         Aux_Error( ERROR_INFO, "SonLv %d, SonPID %d, NPar = %d != 0 !!\n",
                    SonLv, SonPID, amr->patch[0][SonLv][SonPID]->NPar );
   }
#  endif


// 1-2. record the father-buffer patches to actually receive particles
   for (int t=0; t<RefineF2S_Recv_NPatchTotal; t++)
   {
      SonPID = RefineF2S_Recv_PIDList[t];
      FaPID  = amr->patch[0][SonLv][SonPID]->father;

//    check: father patches must be buffer patches
#     ifdef DEBUG_PARTICLE
      if ( FaPID < amr->NPatchComma[FaLv][1] )
         Aux_Error( ERROR_INFO, "This is NOT a buffer patch (FaLv %d, SonPID %d, FaPID %d, NReal %d) !!\n",
                    FaLv, SonPID, FaPID, amr->NPatchComma[FaLv][1] );
#     endif

      RefineF2S_Recv_PIDList[t] = FaPID;
   }


// 2. exchange particles between different ranks
   Par_LB_ExchangeParticleBetweenPatch(
      FaLv,
      RefineF2S_Send_NPatchTotal, RefineF2S_Send_PIDList, RefineF2S_Send_NPatchEachRank,
      RefineF2S_Recv_NPatchTotal, RefineF2S_Recv_PIDList, RefineF2S_Recv_NPatchEachRank,
      NULL, NULL );


// 3. pass particles from father-buffer patches to their real son patches in the same rank
   for (int t=0; t<RefineF2S_Recv_NPatchTotal; t++)
      Par_PassParticle2Son( FaLv, RefineF2S_Recv_PIDList[t] );


// 4. free memory (RefineF2S_Send_PIDList and RefineF2S_Send_LBIdxList will be free'd by LB_Refine)
   delete [] RefineF2S_Recv_PIDList;

} // FUNCTION : Par_LB_Refine_SendParticle2Son



#endif // #if ( defined PARTICLE  &&  defined LOAD_BALANCE )
