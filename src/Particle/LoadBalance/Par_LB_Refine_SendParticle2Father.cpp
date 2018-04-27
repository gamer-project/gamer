#include "GAMER.h"

#if ( defined PARTICLE  &&  defined LOAD_BALANCE )




//-------------------------------------------------------------------------------------------------------
// Function    :  Par_LB_Refine_SendParticle2Father
// Description :  Send particles to leaf real patches at FaLv after grid refinement
//
// Note        :  1. Called by LB_Refine()
//                2. After deleting patches at FaLv+1, their particles may temporarily reside in the
//                   buffer-father patches at FaLv. This function send these particles to their corresponding
//                   real patches at Falv.
//                3. Input array RefineS2F_Send_PIDList[] must be free'd manually after calling this function
//
// Parameter   :  FaLv                       : Target father level
//                RefineS2F_Send_NPatchTotal : Total number of buffer patches at FaLv for sending particles
//                RefineS2F_Send_PIDList     : Patch indices at FaLv for sending particles
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Par_LB_Refine_SendParticle2Father( const int FaLv, const int RefineS2F_Send_NPatchTotal, int *RefineS2F_Send_PIDList )
{

#  ifdef DEBUG_PARTICLE
   if ( RefineS2F_Send_NPatchTotal > 0  &&  RefineS2F_Send_PIDList == NULL )
      Aux_Error( ERROR_INFO, "RefineS2F_Send_PIDList == NULL for RefineS2F_Send_NPatchTotal = %d > 0 !!\n",
                 RefineS2F_Send_NPatchTotal );
#  endif


// 1. find the corresponding real patches at FaLv
   const bool UseInputLBIdx_No = false;

   int  RefineS2F_Send_NPatchEachRank[MPI_NRank];
   int  RefineS2F_Recv_NPatchEachRank[MPI_NRank];
   int  RefineS2F_Recv_NPatchTotal;
   int *RefineS2F_Recv_PIDList = NULL;

   Par_LB_MapBuffer2RealPatch( FaLv, RefineS2F_Send_NPatchTotal,
                                     RefineS2F_Send_PIDList,
                                     RefineS2F_Send_NPatchEachRank,
                                     RefineS2F_Recv_NPatchTotal,
                                     RefineS2F_Recv_PIDList,
                                     RefineS2F_Recv_NPatchEachRank,
                                     UseInputLBIdx_No, NULL );
// check
#  ifdef DEBUG_PARTICLE
// no particles should be exchanged within the same rank
   if ( RefineS2F_Send_NPatchEachRank[MPI_Rank] != 0 )
      Aux_Error( ERROR_INFO, "%d patches are sent to itself !!\n", RefineS2F_Send_NPatchEachRank[MPI_Rank] );

   if ( RefineS2F_Recv_NPatchEachRank[MPI_Rank] != 0 )
      Aux_Error( ERROR_INFO, "%d patches are received from itself !!\n", RefineS2F_Recv_NPatchEachRank[MPI_Rank] );

// all recv patches should have no particles and sons
   for (int t=0; t<RefineS2F_Recv_NPatchTotal; t++)
   {
      const int FaPID = RefineS2F_Recv_PIDList[t];

      if ( amr->patch[0][FaLv][FaPID]->son != -1 )
         Aux_Error( ERROR_INFO, "FaLv %d, FaPID %d, SonPID = %d != -1 !!\n", FaLv, FaPID, amr->patch[0][FaLv][FaPID]->son );

      if ( amr->patch[0][FaLv][FaPID]->NPar != 0 )
         Aux_Error( ERROR_INFO, "FaLv %d, FaPID %d, NPar = %d != 0 !!\n", FaLv, FaPID, amr->patch[0][FaLv][FaPID]->NPar );
   }
#  endif // #ifdef DEBUG_PARTICLE


// 2. exchange particles
   Par_LB_ExchangeParticleBetweenPatch(
      FaLv,
      RefineS2F_Send_NPatchTotal, RefineS2F_Send_PIDList, RefineS2F_Send_NPatchEachRank,
      RefineS2F_Recv_NPatchTotal, RefineS2F_Recv_PIDList, RefineS2F_Recv_NPatchEachRank,
      NULL, NULL );


// 3. free memory (RefineS2F_Send_PIDList[] will be free'd by LB_Refine())
   delete [] RefineS2F_Recv_PIDList;

} // FUNCTION : Par_LB_Refine_SendParticle2Father



#endif // #if ( defined PARTICLE  &&  defined LOAD_BALANCE )
