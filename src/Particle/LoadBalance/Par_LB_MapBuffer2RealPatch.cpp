#include "GAMER.h"

#if ( defined PARTICLE  &&  defined LOAD_BALANCE )




//-------------------------------------------------------------------------------------------------------
// Function    :  Par_LB_MapBuffer2RealPatch
// Description :  Obtain the indices of real patches corresponding to the input buffer patches
//
// Note        :  1. Target patches (those in Buff_PIDList) must be buffer patches at level "lv"
//                   --> Unless UseInputLBIdx==true, in which case we don't care whether patches stored in
//                       Buff_PIDList are buffer or real patches
//                       --> Useful for constructing the PID lists for exchanging particles between fathers
//                           and son (the F2S lists constructed in Par_LB_RecordExchangeParticlePatchID)
//                           --> Since in this special case Buff_PIDList actually stores the PID of real father
//                               patches at lv-1
//                2. Arrays "Buff_PIDList, Buff_NPatchEachRank, and Real_NPatchEachRank" must be preallocated
//                   Array "Real_PIDList" will be allocated (using call-by-reference) here and must be free'd
//                   manually by calling "delete []"
//                3. Input array "Buff_PIDList" will be sorted according to their target MPI rank
//                4. Both "Buff_PIDList, Real_NPatchTotal, and Real_PIDList" are inputted as "call-by-reference"
//
// Parameter   :  lv                   : Target refinement level
//                Buff_NPatchTotal     : Total number of buffer patches in the Buff_PIDList
//                Buff_PIDList         : Input array storing the target buffer patch indices
//                                       --> It will be sorted after calling this function
//                Buff_NPatchEachRank  : Output array to store the number of buffer patches belonging to each rank
//                                       (after they are mapped to real patches)
//                Real_NPatchTotal     : Total number of real patches in this rank mapping to the input buffer patches
//                                       (call-by-reference)
//                Real_PIDList         : Output array to store the real patch indices mapping to the input buffer patches
//                                       --> It will be allocated (using call-by-reference) by this function and must be
//                                           free'd manulally later
//                Real_NPatchEachRank  : Output array to store the number of real patches belonging to each rank
//                UseInputLBIdx        : true --> Use the input load-balance indices instead of recalculating it
//                Buff_LBIdxList_Input : Load-balance indices used by UseInputLBIdx
//
// Return      :  Buff_PIDList, Buff_NPatchEachRank, Real_NPatchTotal, Real_PIDList, Real_NPatchEachRank
//-------------------------------------------------------------------------------------------------------
void Par_LB_MapBuffer2RealPatch( const int lv, const int  Buff_NPatchTotal, int *&Buff_PIDList, int *Buff_NPatchEachRank,
                                                     int &Real_NPatchTotal, int *&Real_PIDList, int *Real_NPatchEachRank,
                                 const bool UseInputLBIdx, long *Buff_LBIdxList_Input )
{

#  ifdef DEBUG_PARTICLE
   if ( lv < 0  ||  lv >= NLEVEL )     Aux_Error( ERROR_INFO, "incorrect target level (%d) !!\n", lv );
   if ( Buff_NPatchTotal < 0 )         Aux_Error( ERROR_INFO, "Buff_NPatchTotal = %d < 0 !!\n", Buff_NPatchTotal );
   else if ( Buff_NPatchTotal > 0 )
   {
      if ( Buff_PIDList == NULL )
         Aux_Error( ERROR_INFO, "Buff_PIDList == NULL (Buff_NPatchTotal = %d) !!\n", Buff_NPatchTotal );

      if ( Buff_NPatchEachRank == NULL )
         Aux_Error( ERROR_INFO, "Buff_NPatchEachRank == NULL (Buff_NPatchTotal = %d) !!\n", Buff_NPatchTotal );

      if ( Real_NPatchEachRank == NULL )
         Aux_Error( ERROR_INFO, "Real_NPatchEachRank == NULL (Buff_NPatchTotal = %d) !!\n", Buff_NPatchTotal );

      if ( UseInputLBIdx  &&  Buff_LBIdxList_Input == NULL )
         Aux_Error( ERROR_INFO, "Buff_LBIdxList_Input == NULL (Buff_NPatchTotal = %d) !!\n", Buff_NPatchTotal );
   }

// PID stored in Buff_PIDList needs to be buffer patches only if UseInputLBIdx == false
   if ( !UseInputLBIdx )
   for (int t=0; t<Buff_NPatchTotal; t++)
   {
      if ( Buff_PIDList[t] < amr->NPatchComma[lv][1]  ||  Buff_PIDList[t] >= amr->num[lv] )
         Aux_Error( ERROR_INFO, "This function should only be applied to buffer patches (t %d, PID %d, NReal %d, NTotal %d) !!\n",
                    t, Buff_PIDList[t], amr->NPatchComma[lv][1], amr->num[lv] );
   } // for (int t=0; t<Buff_NPatchTotal; t++)
#  endif // #ifdef DEBUG_PARTICLE


// must NOT call "return" here even if Buff_NPatchTotal==0 since this rank still needs to SEND data to other ranks
// if ( Buff_NPatchTotal == 0 )  return;


// 1. exchange the target load-balance indices between all ranks
// 1-1. get the LB_Idx of all target buffer patches
// --> note that the LB_Idx stored in each patch always assumes periodicity
// --> so external buffer patches can still find the corresponding real patches correctly
   long *Buff_LBIdxList = NULL;

   if ( UseInputLBIdx )
      Buff_LBIdxList = Buff_LBIdxList_Input;

   else
   {
      Buff_LBIdxList = new long [Buff_NPatchTotal];
      for (int t=0; t<Buff_NPatchTotal; t++)    Buff_LBIdxList[t] = amr->patch[0][lv][ Buff_PIDList[t] ]->LB_Idx;
   }

// 1-2. get the number of patches to be exchanged between different ranks
   int TRank;

   for (int r=0; r<MPI_NRank; r++)  Buff_NPatchEachRank[r] = 0;

   for (int t=0; t<Buff_NPatchTotal; t++)
   {
      TRank = LB_Index2Rank( lv, Buff_LBIdxList[t], CHECK_ON );

      Buff_NPatchEachRank[TRank] ++;
   }

   MPI_Alltoall( Buff_NPatchEachRank, 1, MPI_INT, Real_NPatchEachRank, 1, MPI_INT, MPI_COMM_WORLD );

   Real_NPatchTotal = 0;
   for (int r=0; r<MPI_NRank; r++)  Real_NPatchTotal += Real_NPatchEachRank[r];

// 1-3. set MPI send/recv counts and displacements for exchanging LBIdx
   int *SendCount_LBIdxList = Buff_NPatchEachRank;
   int *RecvCount_LBIdxList = Real_NPatchEachRank;
   int *SendDisp_LBIdxList  = new int [MPI_NRank];
   int *RecvDisp_LBIdxList  = new int [MPI_NRank];

   SendDisp_LBIdxList[0] = 0;
   RecvDisp_LBIdxList[0] = 0;
   for (int r=1; r<MPI_NRank; r++)
   {
      SendDisp_LBIdxList[r] = SendDisp_LBIdxList[r-1] + SendCount_LBIdxList[r-1];
      RecvDisp_LBIdxList[r] = RecvDisp_LBIdxList[r-1] + RecvCount_LBIdxList[r-1];
   }

// 1-4. prepare the send buffer of LBIdx
   long *Buff_LBIdxList_Sort = new long [Buff_NPatchTotal];
   long *Real_LBIdxList_Sort = new long [Real_NPatchTotal];
   int  *Offset_EachRank     = new int  [MPI_NRank];
   int  *Buff_PIDList_Sort   = new int  [Buff_NPatchTotal];

#  ifdef DEBUG_PARTICLE
   for (int t=0; t<Buff_NPatchTotal; t++)    Buff_PIDList_Sort[t] = -1;
#  endif

   for (int r=0; r<MPI_NRank; r++)  Offset_EachRank[r] = SendDisp_LBIdxList[r];

   for (int t=0; t<Buff_NPatchTotal; t++)
   {
      TRank = LB_Index2Rank( lv, Buff_LBIdxList[t], CHECK_ON );

      Buff_LBIdxList_Sort[ Offset_EachRank[TRank] ] = Buff_LBIdxList[t];
      Buff_PIDList_Sort  [ Offset_EachRank[TRank] ] = Buff_PIDList  [t];

      Offset_EachRank[TRank] ++;
   }

#  ifdef DEBUG_PARTICLE
   for (int t=0; t<Buff_NPatchTotal; t++)
      if ( Buff_PIDList_Sort[t] == -1 )   Aux_Error( ERROR_INFO, "Buff_PIDList_Sort[%d] == -1 !!\n", t );
#  endif

// 1-5. collect LBIdx from all ranks
   MPI_Alltoallv( Buff_LBIdxList_Sort, SendCount_LBIdxList, SendDisp_LBIdxList, MPI_LONG,
                  Real_LBIdxList_Sort, RecvCount_LBIdxList, RecvDisp_LBIdxList, MPI_LONG, MPI_COMM_WORLD );

// 1-6. store the sorted PID list
   memcpy( Buff_PIDList, Buff_PIDList_Sort, Buff_NPatchTotal*sizeof(int) );
   delete [] Buff_PIDList_Sort;
// do NOT use the following codes since we don't know whether Buff_PIDList is allocated by malloc or new
// delete [] Buff_PIDList;
// Buff_PIDList = Buff_PIDList_Sort


// 2. get the PID list corresponding to the received LBIdx list
   if ( Real_PIDList != NULL )   delete [] Real_PIDList;

   int *Real_LBIdxList_Sort_IdxTable = new int [Real_NPatchTotal];
   int *Match_LBIdxList              = new int [Real_NPatchTotal];
        Real_PIDList                 = new int [Real_NPatchTotal];

#  ifdef DEBUG_PARTICLE
   for (int t=0; t<Real_NPatchTotal; t++)    Real_PIDList[t] = -1;
#  endif

// 2-1. sort the received LBIdxlist again and find the matching real patch indices
   Mis_Heapsort( Real_NPatchTotal, Real_LBIdxList_Sort, Real_LBIdxList_Sort_IdxTable );

   Mis_Matching_int( amr->NPatchComma[lv][1], amr->LB->IdxList_Real[lv], Real_NPatchTotal, Real_LBIdxList_Sort,
                     Match_LBIdxList );

// 2-2. store the PID mapped to the received LBIdx list
   for (int t=0; t<Real_NPatchTotal; t++)
   {
#     ifdef DEBUG_PARTICLE
      if ( Match_LBIdxList[t] == -1 )
         Aux_Error( ERROR_INFO, "LBIdx (%ld) found no match (lv %d) !!\n", Real_LBIdxList_Sort[t], lv );
#     endif

      Real_PIDList[ Real_LBIdxList_Sort_IdxTable[t] ] = amr->LB->IdxList_Real_IdxTable[lv][ Match_LBIdxList[t] ];
   }

#  ifdef DEBUG_PARTICLE
   for (int t=0; t<Real_NPatchTotal; t++)
      if ( Real_PIDList[t] < 0  ||  Real_PIDList[t] >= amr->NPatchComma[lv][1] )
         Aux_Error( ERROR_INFO, "incorrect PID (lv %d, t %d, PID %d, NReal %d) !!\n",
                    lv, t, Real_PIDList[t], amr->NPatchComma[lv][1] );
#  endif


// 3. free memory
   if ( !UseInputLBIdx )
   delete [] Buff_LBIdxList;
   delete [] SendDisp_LBIdxList;
   delete [] RecvDisp_LBIdxList;
   delete [] Buff_LBIdxList_Sort;
   delete [] Real_LBIdxList_Sort;
   delete [] Offset_EachRank;
   delete [] Real_LBIdxList_Sort_IdxTable;
   delete [] Match_LBIdxList;

} // FUNCTION : Par_LB_MapBuffer2RealPatch



#endif // #if ( defined PARTICLE  &&  defined LOAD_BALANCE )
