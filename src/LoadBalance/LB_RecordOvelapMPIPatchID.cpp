#include "GAMER.h"

#ifdef LOAD_BALANCE




//-------------------------------------------------------------------------------------------------------
// Function    :  LB_RecordOverlapMPIPatchID
// Description :  Construct the patch indices for overlapping MPI time with CPU/GPU computation
//
// Note        :  This function must be invoked AFTER LB_RecordExchangeDataPatchID()
//
// Parameter   :  Lv : Target refinement level for recording lists
//-------------------------------------------------------------------------------------------------------
void LB_RecordOverlapMPIPatchID( const int Lv )
{

// check
#  ifdef MHD
   Aux_Error( ERROR_INFO, "%s does NOT support MHD !!\n", __FUNCTION__ );
#  endif


   const int NReal0 = amr->NPatchComma[Lv][1]/8;

   int  *LB_SendH_NList  = amr->LB->SendH_NList [Lv];
   int **LB_SendH_IDList = amr->LB->SendH_IDList[Lv];
   bool *FluSyncList     = new bool [NReal0];
#  ifdef GRAVITY
   int  *LB_SendG_NList  = amr->LB->SendG_NList [Lv];
   int **LB_SendG_IDList = amr->LB->SendG_IDList[Lv];
   bool *PotSyncList     = new bool [NReal0];
#  endif


// 1. initialize arrays
// ============================================================================================================
   amr->LB->OverlapMPI_FluSyncN [Lv] = 0;
   amr->LB->OverlapMPI_FluAsyncN[Lv] = 0;

   if ( amr->LB->OverlapMPI_FluSyncPID0 [Lv] != NULL )   delete [] amr->LB->OverlapMPI_FluSyncPID0 [Lv];
   if ( amr->LB->OverlapMPI_FluAsyncPID0[Lv] != NULL )   delete [] amr->LB->OverlapMPI_FluAsyncPID0[Lv];

   amr->LB->OverlapMPI_FluSyncPID0 [Lv] = new int [NReal0];
   amr->LB->OverlapMPI_FluAsyncPID0[Lv] = new int [NReal0];

   for (int t=0; t<NReal0; t++)  FluSyncList[t] = false;

#  ifdef GRAVITY
   amr->LB->OverlapMPI_PotSyncN [Lv] = 0;
   amr->LB->OverlapMPI_PotAsyncN[Lv] = 0;

   if ( amr->LB->OverlapMPI_PotSyncPID0 [Lv] != NULL )   delete [] amr->LB->OverlapMPI_PotSyncPID0 [Lv];
   if ( amr->LB->OverlapMPI_PotAsyncPID0[Lv] != NULL )   delete [] amr->LB->OverlapMPI_PotAsyncPID0[Lv];

   amr->LB->OverlapMPI_PotSyncPID0 [Lv] = new int [NReal0];
   amr->LB->OverlapMPI_PotAsyncPID0[Lv] = new int [NReal0];

   for (int t=0; t<NReal0; t++)  PotSyncList[t] = false;
#  endif


// 2. construct the SyncList (true/false --> sync/async)
// ============================================================================================================
   int ID0;

   for (int r=0; r<MPI_NRank; r++)
   for (int t=0; t<LB_SendH_NList[r]; t++)
   {
      ID0 = LB_SendH_IDList[r][t] / 8;

#     ifdef GAMER_DEBUG
      if ( ID0 >= NReal0 )    Aux_Error( ERROR_INFO, "ID0 (%d) >= NReal0 (%d) (fluid) !!\n", ID0, NReal0 );
#     endif

      FluSyncList[ID0] = true;
#     ifdef GRAVITY
      PotSyncList[ID0] = true;
#     endif
   }

#  ifdef GRAVITY
   for (int r=0; r<MPI_NRank; r++)
   for (int t=0; t<LB_SendG_NList[r]; t++)
   {
      ID0 = LB_SendG_IDList[r][t] / 8;

#     ifdef GAMER_DEBUG
      if ( ID0 >= NReal0 )    Aux_Error( ERROR_INFO, "ID0 (%d) >= NReal0 (%d) (gravity) !!\n", ID0, NReal0 );
#     endif

      PotSyncList[ID0] = true;
   }
#  endif


// 3. record the Sync and Async lists
// ============================================================================================================
   int PID0;

   for (int t=0; t<NReal0; t++)
   {
      PID0 = t*8;

      if ( FluSyncList[t] )   amr->LB->OverlapMPI_FluSyncPID0 [Lv][ amr->LB->OverlapMPI_FluSyncN [Lv] ++ ] = PID0;
      else                    amr->LB->OverlapMPI_FluAsyncPID0[Lv][ amr->LB->OverlapMPI_FluAsyncN[Lv] ++ ] = PID0;

#     ifdef GRAVITY
      if ( PotSyncList[t] )   amr->LB->OverlapMPI_PotSyncPID0 [Lv][ amr->LB->OverlapMPI_PotSyncN [Lv] ++ ] = PID0;
      else                    amr->LB->OverlapMPI_PotAsyncPID0[Lv][ amr->LB->OverlapMPI_PotAsyncN[Lv] ++ ] = PID0;
#     endif
   }


// free memory
   delete [] FluSyncList;
#  ifdef GRAVITY
   delete [] PotSyncList;
#  endif

} // FUNCTION : LB_RecordOverlapMPIPatchID



#endif // #ifdef LOAD_BALANCE
