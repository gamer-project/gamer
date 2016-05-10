#include "Copyright.h"
#include "GAMER.h"

#ifdef LOAD_BALANCE



static void LB_RedistributeRealPatch( const int lv );




//-------------------------------------------------------------------------------------------------------
// Function    :  LB_Init_LoadBalance
// Description :  Initialize the load-balance process 
//
// Note        :  1. Patches at all levels will be redistributed, and all patch relations will be reconstructed
//                2. All parameters in the data structure "LB_t LB" will be reconstructed 
//                3. Data structure "ParaVar_t ParaVar" will be removed
//                4. This function is used in both initialization phase and run-time data redistribution
//                5. Data in the buffer patches will be filled up
//                6. Arrays "NPatchTotal" and "NPatchComma" must be prepared in advance
//
// Parameter   :  DuringRestart  : true --> This function is invoked during the RESTART process
//                                      --> In this case, no "LB_SetCutPoint" and "LB_RedistributeRealPatch" are required 
//                                          since these tasks are already done in the function "Init_Reload"
//-------------------------------------------------------------------------------------------------------
void LB_Init_LoadBalance( const bool DuringRestart )
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );


// check   
   if ( amr->LB == NULL )  Aux_Error( ERROR_INFO, "amr->LB has not been allocated !!\n" );


// check the synchronization
   for (int lv=1; lv<NLEVEL; lv++)
      if ( NPatchTotal[lv] != 0 )   Mis_CompareRealValue( Time[0], Time[lv], __FUNCTION__, true );


// 0. delete ParaVar which is no longer useful
   if ( amr->ParaVar != NULL )
   {
      delete amr->ParaVar;
      amr->ParaVar = NULL;
   }


// 1. re-distribute and allocate all patches   
   const bool InputLBIdxList_No = false;

   for (int lv=0; lv<NLEVEL; lv++)
   {
      if ( OPT__VERBOSE  &&  MPI_Rank == 0 )    Aux_Message( stdout, "   Re-distributing patches at Lv %2d ... \n", lv );

      if ( !DuringRestart )
      {
//       1.1 set the load-balance cut points
         LB_SetCutPoint( lv, amr->LB->CutPoint[lv], InputLBIdxList_No, NULL );

//       1.2 re-distribute real patches
         LB_RedistributeRealPatch( lv );
      }

//    1.3 allocate sibling-buffer patches at lv
      LB_AllocateBufferPatch_Sibling( lv );

//    1.4 allocate father-buffer patches at lv-1
      if ( lv > 0 )
      LB_AllocateBufferPatch_Father( lv, true, NULL_INT, NULL, false, NULL, NULL );
   }


// 2. contruct the patch relation
   for (int lv=0; lv<NLEVEL; lv++)
   {
      if ( OPT__VERBOSE  &&  MPI_Rank == 0 )    Aux_Message( stdout, "   Constructing patch relation at Lv %2d ... \n", lv );

      LB_SiblingSearch( lv, true, NULL_INT, NULL );

      if ( lv > 0 )
      {
         LB_FindFather    ( lv,   true, NULL_INT, NULL );
         LB_FindSonNotHome( lv-1, true, NULL_INT, NULL );
      }
   }


// 3. construct the MPI send and recv data list 
   for (int lv=0; lv<NLEVEL; lv++)
   {
//    3.1 list for exchanging hydro and potential data
      LB_RecordExchangeDataPatchID( lv, false );

//    3.2 list for exchanging restricted hydro data
#     ifndef GAMER_DEBUG
      if ( OPT__FIXUP_RESTRICT )    
#     endif
      LB_RecordExchangeRestrictDataPatchID( lv );

//    3.3 list for exchanging hydro fluxes (also allocate flux arrays)
      if ( amr->WithFlux )        
      LB_AllocateFluxArray( lv ); 

//    3.4 list for exchanging hydro data after the fix-up operation
      if ( OPT__FIXUP_RESTRICT  ||  OPT__FIXUP_FLUX )
      LB_RecordExchangeFixUpDataPatchID( lv );

//    3.5 list for overlapping MPI time with CPU/GPU computation
      if ( OPT__OVERLAP_MPI )
      LB_RecordOverlapMPIPatchID( lv );
   }


// 4. get the buffer data   
   for (int lv=0; lv<NLEVEL; lv++)  
   {
      Buf_GetBufferData( lv, amr->FluSg[lv], NULL_INT, DATA_GENERAL,    _FLU,  Flu_ParaBuf, USELB_YES );

#     ifdef GRAVITY
      Buf_GetBufferData( lv, NULL_INT, amr->PotSg[lv], POT_FOR_POISSON, _POTE, Pot_ParaBuf, USELB_YES );
#     endif
   }


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : LB_Init_LoadBalance



//-------------------------------------------------------------------------------------------------------
// Function    :  LB_SetCutPoint
// Description :  Set the range of LB_Idx for distributing patches to different ranks 
//
// Note        :  1. This function assumes that "NPatchTotal[lv]" has already been set by 
//                   invoking the function "Mis_GetTotalPatchNumber( lv )"
//                2. Input array "CutPoint" will be set in this function
//                3. Real patches with LB_Idx in the range "CutPoint[r] <= LB_Idx < CutPoint[r+1]"
//                   will be calculated by rank "r"
//
// Parameter   :  lv             : Targeted refinement level
//                CutPoint       : Cut point array to be set
//                InputLBIdxList : Input "LBIdx_AllRank"   --> useful during RESTART
//                LBIdx_AllRank  : LBIdx list of all ranks --> useful during RESTART
//                                 (only rank 0 needs to provide this list, which can be unsorted)
//
// Return      :  CutPoint
//-------------------------------------------------------------------------------------------------------
void LB_SetCutPoint( const int lv, long *CutPoint, const bool InputLBIdxList, long *LBIdx_AllRank )
{

// check
   if ( NPatchTotal[lv]%8 != 0 )    
      Aux_Error( ERROR_INFO, "NPatchTotal[%d] = %d is NOT a multiple of 8 !!\n", lv, NPatchTotal[lv] );

   if ( MPI_Rank == 0  &&  InputLBIdxList  &&  LBIdx_AllRank == NULL )
      Aux_Error( ERROR_INFO, "LBIdx_AllRank == NULL when the option InputLBIdxList is on !!\n" );


// 1. determine the load of each MPI rank   
   const int NPG            = NPatchTotal[lv]/8;
   const int NPG_per_Rank   = NPG/MPI_NRank; 
   const int Rank_with_more = NPG%MPI_NRank;
   int *LoadList = NULL;
   
   if ( MPI_Rank == 0 )
   {
      LoadList= new int [MPI_NRank];

      for (int r=0; r<MPI_NRank; r++)
      {
         LoadList[r] = NPG_per_Rank*8;

         if ( r < Rank_with_more )  LoadList[r] += 8;
      }

#     ifdef GAMER_DEBUG
      int LoadSum = 0;
      for (int r=0; r<MPI_NRank; r++)  LoadSum += LoadList[r];
      if ( LoadSum != NPatchTotal[lv] )
         Aux_Error( ERROR_INFO, "LoadSum (%d) != NPatchTotal (%d) at level %d !!\n", 
                    LoadSum, NPatchTotal[lv], lv );
#     endif
   }


// 2. collect LB_Idx from all ranks 
   long *SendBuf_LBIdx = NULL;
   long *RecvBuf_LBIdx = NULL; 
   int  *NPatch_each_Rank = NULL, *Recv_Disp = NULL;

   if ( InputLBIdxList )
   {
      if ( MPI_Rank == 0 )    RecvBuf_LBIdx = LBIdx_AllRank; 
   }

   else
   {
      SendBuf_LBIdx = new long [ amr->NPatchComma[lv][1] ];

      if ( MPI_Rank == 0 )    
      {
         RecvBuf_LBIdx    = new long [ NPatchTotal[lv] ];
         NPatch_each_Rank = new int  [ MPI_NRank ];
         Recv_Disp        = new int  [ MPI_NRank ];
      }

      for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)   SendBuf_LBIdx[PID] = amr->patch[0][lv][PID]->LB_Idx;

      MPI_Gather( &amr->NPatchComma[lv][1], 1, MPI_INT, NPatch_each_Rank, 1, MPI_INT, 0, MPI_COMM_WORLD );

      if ( MPI_Rank == 0 )
      {
         Recv_Disp[0] = 0;
         for (int r=0; r<MPI_NRank-1; r++)   Recv_Disp[r+1] = Recv_Disp[r] + NPatch_each_Rank[r]; 
      }

      MPI_Gatherv( SendBuf_LBIdx, amr->NPatchComma[lv][1], MPI_LONG, RecvBuf_LBIdx, NPatch_each_Rank, Recv_Disp,
                   MPI_LONG, 0, MPI_COMM_WORLD );
   } // if ( InputLBIdxList ) ... else ...


   if ( MPI_Rank == 0 )
   {
//    3. sort LB_Idx
      Mis_Heapsort( NPatchTotal[lv], RecvBuf_LBIdx, NULL );


//    4. set the cut points   
      if ( NPatchTotal[lv] == 0 )
         for (int t=0; t<MPI_NRank+1; t++)   CutPoint[t] = -1;

      else
      {
         const long LBIdx_Max = RecvBuf_LBIdx[ NPatchTotal[lv] - 1 ];
         int Counter = 0;    

         for (int r=0; r<MPI_NRank; r++)   
         {
            CutPoint[r] = ( Counter < NPatchTotal[lv] ) ? RecvBuf_LBIdx[Counter] : LBIdx_Max+1;
            Counter += LoadList[r];
         }

         CutPoint[MPI_NRank] = LBIdx_Max + 1;
      }
   }
    

// 5. broadcast the cut points
   MPI_Bcast( CutPoint, MPI_NRank+1, MPI_LONG, 0, MPI_COMM_WORLD );


// 6. output the cut points and load in each MPI rank
   if ( OPT__VERBOSE  &&  MPI_Rank == 0 )
   {
      const double Load_Ave = (double)NPatchTotal[lv] / MPI_NRank;
      int Load_Max = -1;

      for (int r=0; r<MPI_NRank; r++)
      {
         Aux_Message( stdout, "   Lv %2d: Rank %4d, Cut %10ld -> %10ld, NPatch %10d\n", 
                      lv, r, CutPoint[r], CutPoint[r+1], LoadList[r] );

         if ( LoadList[r] > Load_Max )   Load_Max = LoadList[r];
      }

      Aux_Message( stdout, "   Load_Ave %9.3e, Load_Max %8d --> Load Imbalance = %6.2f%%\n", 
                   Load_Ave, Load_Max, (Load_Max==0) ? 0.0 : 100.0*(Load_Max-Load_Ave)/Load_Ave );
      Aux_Message( stdout, "   =============================================================================\n" );
   }


   if ( !InputLBIdxList )  delete [] SendBuf_LBIdx;
   if ( MPI_Rank == 0 )
   {
      delete [] LoadList;

      if ( !InputLBIdxList )
      {
         delete [] RecvBuf_LBIdx;
         delete [] NPatch_each_Rank;
         delete [] Recv_Disp;
      }
   }

} // FUNCTION : LB_SetCutPoint



//-------------------------------------------------------------------------------------------------------
// Function    :  LB_RedistrubteRealPatch
// Description :  Redistribute real patches to different ranks according to the cut point array
//                "amr->LB->CutPoint[lv]"
//
// Note        :  1. All ranks must have the array "LB_CutPoint" prepared
//                2. This function assumes that the "patch group" is adopted as the basic unit for data 
//                   redistribution
//                3. Real patches with LB_Idx in the range "CutPoint[r] <= LB_Idx < CutPoint[r+1]"
//                   will be sent to rank "r"
//
// Parameter   :  lv : Targeted refinement level
//-------------------------------------------------------------------------------------------------------
void LB_RedistributeRealPatch( const int lv )
{

// 1. count the number of real patches to be sent and received   
// ==========================================================================================
   const int PatchSize1v = CUBE( PATCH_SIZE );

   int  NSend_Total_Patch, NRecv_Total_Patch, TRank;
   long LB_Idx;

   int *Send_NCount_Patch  = new int [MPI_NRank];
   int *Recv_NCount_Patch  = new int [MPI_NRank];
   int *Send_NDisp_Patch   = new int [MPI_NRank];
   int *Recv_NDisp_Patch   = new int [MPI_NRank];
   int *Send_NCount_Data1v = new int [MPI_NRank];
   int *Recv_NCount_Data1v = new int [MPI_NRank];
   int *Send_NDisp_Data1v  = new int [MPI_NRank];
   int *Recv_NDisp_Data1v  = new int [MPI_NRank];
   int *Counter            = new int [MPI_NRank];

   for (int r=0; r<MPI_NRank; r++)  
   {
      Send_NCount_Patch[r] = 0;
      Recv_NCount_Patch[r] = 0;
   }
   Send_NDisp_Patch[0] = 0;
   Recv_NDisp_Patch[0] = 0;

// 1.1 send count 
   for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
   {
      LB_Idx = amr->patch[0][lv][PID]->LB_Idx;
      TRank  = LB_Index2Rank( lv, LB_Idx, CHECK_ON );

      Send_NCount_Patch[TRank] ++;
   }

// 1.2 receive count  
   MPI_Alltoall( Send_NCount_Patch, 1, MPI_INT, Recv_NCount_Patch, 1, MPI_INT, MPI_COMM_WORLD );

// 1.3 send/recv displacement
   for (int r=1; r<MPI_NRank; r++)     
   {
      Send_NDisp_Patch[r] = Send_NDisp_Patch[r-1] + Send_NCount_Patch[r-1];
      Recv_NDisp_Patch[r] = Recv_NDisp_Patch[r-1] + Recv_NCount_Patch[r-1];
   }

// 1.4 send/recv data displacement
   for (int r=0; r<MPI_NRank; r++)     
   {
      Send_NCount_Data1v[r] = PatchSize1v*Send_NCount_Patch[r];
      Recv_NCount_Data1v[r] = PatchSize1v*Recv_NCount_Patch[r];
      Send_NDisp_Data1v [r] = PatchSize1v*Send_NDisp_Patch [r];
      Recv_NDisp_Data1v [r] = PatchSize1v*Recv_NDisp_Patch [r];
   }

// 1.5 total number of patches to be sent and received
   NSend_Total_Patch = Send_NDisp_Patch[ MPI_NRank-1 ] + Send_NCount_Patch[ MPI_NRank-1 ];
   NRecv_Total_Patch = Recv_NDisp_Patch[ MPI_NRank-1 ] + Recv_NCount_Patch[ MPI_NRank-1 ];


// 2. prepare the MPI send buffers
// ==========================================================================================
   const int SendDataSize1v = NSend_Total_Patch*PatchSize1v;
   const int RecvDataSize1v = NRecv_Total_Patch*PatchSize1v;
   const int FluSg          = amr->FluSg[lv];
#  ifdef GRAVITY
   const int PotSg          = amr->PotSg[lv];
#  endif

   real *SendPtr       = NULL;
   long *SendBuf_LBIdx = new long [ NSend_Total_Patch ]; 
   real *SendBuf_Flu   = new real [ SendDataSize1v*NCOMP ]; 
#  ifdef GRAVITY
   real *SendBuf_Pot   = new real [ SendDataSize1v ]; 
#  endif

   for (int r=0; r<MPI_NRank; r++)  Counter[r] = 0;

   for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
   {
      LB_Idx = amr->patch[0][lv][PID]->LB_Idx;
      TRank  = LB_Index2Rank( lv, LB_Idx, CHECK_ON );

//    2.1 LB_Idx
      SendBuf_LBIdx[ Send_NDisp_Patch[TRank] + Counter[TRank] ] = LB_Idx;

//    2.2 fluid
      for (int v=0; v<NCOMP; v++)
      {
         SendPtr = SendBuf_Flu + v*SendDataSize1v + Send_NDisp_Data1v[TRank] + Counter[TRank]*PatchSize1v;
         memcpy( SendPtr, &amr->patch[FluSg][lv][PID]->fluid[v][0][0][0], PatchSize1v*sizeof(real) );
      }

//    2.3 potential
#     ifdef GRAVITY
      SendPtr = SendBuf_Pot + Send_NDisp_Data1v[TRank] + Counter[TRank]*PatchSize1v;
      memcpy( SendPtr, &amr->patch[PotSg][lv][PID]->pot[0][0][0], PatchSize1v*sizeof(real) );
#     endif

      Counter[TRank] ++;
   } // for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)


// 4. delete old patches and allocate the MPI recv buffers
// ==========================================================================================
// free memory first to reduce the memory consumption
   amr->Lvdelete( lv );  // NPatchComma is also reset to 0 here

// allocate recv buffers AFTER deleting old patches
   long *RecvBuf_LBIdx = new long [ NRecv_Total_Patch ];
   real *RecvBuf_Flu   = new real [ RecvDataSize1v*NCOMP ];
#  ifdef GRAVITY
   real *RecvBuf_Pot   = new real [ RecvDataSize1v ];
#  endif


// 5. transfer data by MPI_Alltoallv
// ==========================================================================================
// 5.1 LB_Idx
   MPI_Alltoallv( SendBuf_LBIdx, Send_NCount_Patch, Send_NDisp_Patch, MPI_LONG, 
                  RecvBuf_LBIdx, Recv_NCount_Patch, Recv_NDisp_Patch, MPI_LONG, MPI_COMM_WORLD );

// 5.2 fluid (transfer one component at a time to avoid exceeding the maximum allowed transferred size in MPI)
   for (int v=0; v<NCOMP; v++)
   {
#     ifdef FLOAT8   
      MPI_Alltoallv( SendBuf_Flu + v*SendDataSize1v, Send_NCount_Data1v, Send_NDisp_Data1v, MPI_DOUBLE, 
                     RecvBuf_Flu + v*RecvDataSize1v, Recv_NCount_Data1v, Recv_NDisp_Data1v, MPI_DOUBLE, MPI_COMM_WORLD );
#     else
      MPI_Alltoallv( SendBuf_Flu + v*SendDataSize1v, Send_NCount_Data1v, Send_NDisp_Data1v, MPI_FLOAT, 
                     RecvBuf_Flu + v*RecvDataSize1v, Recv_NCount_Data1v, Recv_NDisp_Data1v, MPI_FLOAT,  MPI_COMM_WORLD );
#     endif
   }

// 5.3 potential
#  ifdef GRAVITY
#  ifdef FLOAT8   
   MPI_Alltoallv( SendBuf_Pot, Send_NCount_Data1v, Send_NDisp_Data1v, MPI_DOUBLE, 
                  RecvBuf_Pot, Recv_NCount_Data1v, Recv_NDisp_Data1v, MPI_DOUBLE, MPI_COMM_WORLD );
#  else
   MPI_Alltoallv( SendBuf_Pot, Send_NCount_Data1v, Send_NDisp_Data1v, MPI_FLOAT, 
                  RecvBuf_Pot, Recv_NCount_Data1v, Recv_NDisp_Data1v, MPI_FLOAT,  MPI_COMM_WORLD );
#  endif
#  endif


// 6. deallocate the MPI send buffers (BEFORE creating new patches to reduce the memory consumption)
// ==========================================================================================
   delete [] Send_NCount_Patch;
   delete [] Send_NDisp_Patch;
   delete [] Send_NCount_Data1v;
   delete [] Send_NDisp_Data1v;
   delete [] Counter;
   delete [] SendBuf_LBIdx;
   delete [] SendBuf_Flu;
#  ifdef GRAVITY
   delete [] SendBuf_Pot;
#  endif


// 7. allocate new patches with the data just received (use "patch group" as the basic unit)
// ==========================================================================================
   const int PScale  = PATCH_SIZE*amr->scale[lv];
   const int PGScale = 2*PScale;
   int TPID, TPID0, Cr0[3];
   real *RecvPtr = NULL;

   for (int t=0; t<NRecv_Total_Patch; t+=8)
   {
      LB_Idx = RecvBuf_LBIdx[t]; 

      LB_Index2Corner( lv, LB_Idx, Cr0, CHECK_ON );
      
      for (int d=0; d<3; d++)    Cr0[d] -= Cr0[d]%PGScale; // currently this line has no effect 

//    father patch is still unkown ...      
      amr->pnew( lv, Cr0[0],        Cr0[1],        Cr0[2],        -1, true, true );
      amr->pnew( lv, Cr0[0]+PScale, Cr0[1],        Cr0[2],        -1, true, true );
      amr->pnew( lv, Cr0[0],        Cr0[1]+PScale, Cr0[2],        -1, true, true );
      amr->pnew( lv, Cr0[0],        Cr0[1],        Cr0[2]+PScale, -1, true, true );
      amr->pnew( lv, Cr0[0]+PScale, Cr0[1]+PScale, Cr0[2],        -1, true, true );
      amr->pnew( lv, Cr0[0],        Cr0[1]+PScale, Cr0[2]+PScale, -1, true, true );
      amr->pnew( lv, Cr0[0]+PScale, Cr0[1],        Cr0[2]+PScale, -1, true, true );
      amr->pnew( lv, Cr0[0]+PScale, Cr0[1]+PScale, Cr0[2]+PScale, -1, true, true );

//    assign data      
      TPID0 = amr->num[lv] - 8;

      for (int m=0; m<8; m++)
      {
         TPID = TPID0 + m;

//       fluid
         for (int v=0; v<NCOMP; v++)
         {
            RecvPtr = RecvBuf_Flu + v*RecvDataSize1v + (t+m)*PatchSize1v;
            memcpy( &amr->patch[FluSg][lv][TPID]->fluid[v][0][0][0], RecvPtr, PatchSize1v*sizeof(real) );
         }

//       potential
#        ifdef GRAVITY
         RecvPtr = RecvBuf_Pot + (t+m)*PatchSize1v;
         memcpy( &amr->patch[PotSg][lv][TPID]->pot[0][0][0], RecvPtr, PatchSize1v*sizeof(real) );
#        endif
      }
   } // for (int t=0; t<NRecv_Total_Patch; t+=8)

// reset NPatchComma
   for (int m=1; m<28; m++)   amr->NPatchComma[lv][m] = NRecv_Total_Patch;

// check the amr->NPatchComma recording
   if ( amr->NPatchComma[lv][1] != amr->num[lv] )
      Aux_Error( ERROR_INFO, "amr->NPatchComma[%d][1] (%d) != amr->num[%d] (%d) !!\n",
                 lv, amr->NPatchComma[lv][1], lv, amr->num[lv] );


// 8. record LB_IdxList_Real
// ==========================================================================================
   if ( amr->LB->IdxList_Real         [lv] != NULL )  delete [] amr->LB->IdxList_Real         [lv];
   if ( amr->LB->IdxList_Real_IdxTable[lv] != NULL )  delete [] amr->LB->IdxList_Real_IdxTable[lv];

   amr->LB->IdxList_Real         [lv] = new long [NRecv_Total_Patch];
   amr->LB->IdxList_Real_IdxTable[lv] = new int  [NRecv_Total_Patch];

   for (int PID=0; PID<NRecv_Total_Patch; PID++)   amr->LB->IdxList_Real[lv][PID] = amr->patch[0][lv][PID]->LB_Idx;

   Mis_Heapsort( NRecv_Total_Patch, amr->LB->IdxList_Real[lv], amr->LB->IdxList_Real_IdxTable[lv] );


// 9. deallocate the MPI recv buffers
// ==========================================================================================
   delete [] Recv_NCount_Patch;
   delete [] Recv_NDisp_Patch;
   delete [] Recv_NCount_Data1v;
   delete [] Recv_NDisp_Data1v;
   delete [] RecvBuf_LBIdx;
   delete [] RecvBuf_Flu;
#  ifdef GRAVITY
   delete [] RecvBuf_Pot;
#  endif

} // FUNCTION : LB_RedistributePatch



#endif // #ifdef LOAD_BALANCE
