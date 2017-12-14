#include "GAMER.h"

#ifdef LOAD_BALANCE



static void LB_RedistributeRealPatch( const int lv, real **ParVar_Old, real **Passive_Old );
#ifdef PARTICLE
static void LB_RedistributeParticle_Init( real **ParVar_Old, real **Passive_Old );
static void LB_RedistributeParticle_End( real **ParVar_Old, real **Passive_Old );
#endif




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
//                7. Particles will also be redistributed
//
// Parameter   :  Redistribute : true  --> Redistribute all real patches according to the load-balance weighting of each patch
//                                         and initialize all load-balance related set-up
//                               false --> Initialize all load-balance related set-up, but do NOT invoke "LB_SetCutPoint"
//                                         and "LB_RedistributeRealPatch" to redistribute all real patches
//                                     --> Currently it is used only during the RESTART process since we already call
//                                         LB_SetCutPoint() and load real patches accordingly when calling Init_Reload()
//                ParWeight    : Relative load-balance weighting of particles
//                               --> Weighting of each patch is estimated as "PATCH_SIZE^3 + NParThisPatch*ParWeight"
//                               --> <= 0.0 : do not consider particle weighting
//                                            --> Currently we force ParWeight==0.0 when calling LB_Init_LoadBalance()
//                                                for the first time in Init_GAMER() and main() since we don't have enough
//                                                information for calculating particle weighting at that time
//                Reset        : Call LB->reset() to reset the load-balance variables
//                               --> Note that CutPoint[] will NOT be reset even when "Reset == true"
//-------------------------------------------------------------------------------------------------------
void LB_Init_LoadBalance( const bool Redistribute, const double ParWeight, const bool Reset )
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );


// check
   if ( amr->LB == NULL )  Aux_Error( ERROR_INFO, "amr->LB has not been allocated !!\n" );


// check the synchronization
   for (int lv=1; lv<NLEVEL; lv++)
      if ( NPatchTotal[lv] != 0 )   Mis_CompareRealValue( Time[0], Time[lv], __FUNCTION__, true );


// delete ParaVar which is no longer useful
   if ( amr->ParaVar != NULL )
   {
      delete amr->ParaVar;
      amr->ParaVar = NULL;
   }


// 1. set up the load-balance cut points (must do this before calling LB_RedistributeParticle_Init)
   const bool InputLBIdxAndLoad_No = false;

   if ( Redistribute )
   for (int lv=0; lv<NLEVEL; lv++)
   {
      if ( OPT__VERBOSE  &&  MPI_Rank == 0 )    Aux_Message( stdout, "   Calculating load-balance indices at Lv %2d ... \n", lv );

      LB_SetCutPoint( lv, amr->LB->CutPoint[lv], InputLBIdxAndLoad_No, NULL, NULL, ParWeight );
   }


// 2. reinitialize arrays used by the load-balance routines
//    --> must do this AFTER calling LB_SetCutPoint() since it still needs to access load-balance information when PARTICLE is on
//        --> for example, LB_EstimateWorkload_AllPatchGroup()->Par_CollectParticle2OneLevel()->Par_LB_CollectParticle2OneLevel(),
//            which will access amr->LB->IdxList_Real. But amr->LB->IdxList_Real will be reset when calling amr->LB->reset()
//    --> amr->LB->reset() will NOT reset CutPoint[] (otherwise it will just overwrite the cut points set by LB_SetCutPoint())
   if ( Reset )   amr->LB->reset();


// 3. re-distribute and allocate all patches (and particles)
#  ifdef PARTICLE
   real  *ParVar_Old [PAR_NVAR    ];
   real  *Passive_Old[PAR_NPASSIVE];
#  else
   real **ParVar_Old  = NULL;
   real **Passive_Old = NULL;
#  endif

#  ifdef PARTICLE
   if ( Redistribute )  LB_RedistributeParticle_Init( ParVar_Old, Passive_Old );
#  endif

   for (int lv=0; lv<NLEVEL; lv++)
   {
      if ( OPT__VERBOSE  &&  MPI_Rank == 0 )    Aux_Message( stdout, "   Re-distributing patches at Lv %2d ... \n", lv );

//    3.1 re-distribute real patches
      if ( Redistribute )
      LB_RedistributeRealPatch( lv, ParVar_Old, Passive_Old );

//    3.2 allocate sibling-buffer patches at lv
      LB_AllocateBufferPatch_Sibling( lv );

//    3.3 allocate father-buffer patches at lv-1
      if ( lv > 0 )
      LB_AllocateBufferPatch_Father( lv, true, NULL_INT, NULL, false, NULL, NULL );
   }

#  ifdef PARTICLE
   if ( Redistribute )  LB_RedistributeParticle_End( ParVar_Old, Passive_Old );
#  endif


// 4. contruct the patch relation
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


// 5. construct the MPI send and recv data list
   for (int lv=0; lv<NLEVEL; lv++)
   {
//    5.1 list for exchanging hydro and potential data
      LB_RecordExchangeDataPatchID( lv, false );

//    5.2 list for exchanging restricted hydro data
//        --> note that even when OPT__FIXUP_RESTRICT is off we still need to do data restriction in several places
//            (e.g., restart, and OPT__CORR_AFTER_ALL_SYNC)
//        --> for simplicity and sustainability, we always invoke LB_RecordExchangeRestrictDataPatchID()
      LB_RecordExchangeRestrictDataPatchID( lv );

//    5.3 list for exchanging hydro fluxes (also allocate flux arrays)
      if ( amr->WithFlux )
      LB_AllocateFluxArray( lv );

//    5.4 list for exchanging hydro data after the fix-up operation
//        --> for simplicity and sustainability, we always invoke LB_RecordExchangeFixUpDataPatchID()
//        --> see the comments 3.2 above
      LB_RecordExchangeFixUpDataPatchID( lv );

//    5.5 list for overlapping MPI time with CPU/GPU computation
      if ( OPT__OVERLAP_MPI )
      LB_RecordOverlapMPIPatchID( lv );

//    5.6 list for exchanging particles
#     ifdef PARTICLE
      Par_LB_RecordExchangeParticlePatchID( lv );
#     endif
   } // for (int lv=0; lv<NLEVEL; lv++)


// 6. get the buffer data
   for (int lv=0; lv<NLEVEL; lv++)
   {
      Buf_GetBufferData( lv, amr->FluSg[lv], NULL_INT, DATA_GENERAL,    _TOTAL, Flu_ParaBuf, USELB_YES );

#     ifdef GRAVITY
      Buf_GetBufferData( lv, NULL_INT, amr->PotSg[lv], POT_FOR_POISSON, _POTE,  Pot_ParaBuf, USELB_YES );
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
//                   will be set to rank "r"
//                4. The input option "InputLBIdx0AndLoad" is useful during RESTART where we have very
//                   limited information (e.g., we don't know the number of patches in each rank,
//                   amr->NPatchComma, and any particle information yet ...)
//                   --> See the description of "InputLBIdx0AndLoad, LBIdx0_AllRank_Input, and
//                       Load_AllRank_Input" below
//
// Parameter   :  lv                   : Target refinement level
//                CutPoint             : Cut point array to be set
//                InputLBIdx0AndLoad   : Provide both "LBIdx0_AllRank_Input" and "Load_AllRank_Input" directly
//                                       so that they don't have to be collected from all ranks again
//                                       --> Useful during RESTART
//                LBIdx0_AllRank_Input : LBIdx of all patch groups in all ranks
//                                       --> Useful only when InputLBIdx0AndLoad == true
//                                       --> Only need the **minimum** LBIdx in each patch group
//                                       --> Only rank 0 needs to provide this list
//                                       --> Can be unsorted
//                Load_AllRank_Input   : Load-balance weighting of all patch groups in all ranks
//                                       --> Useful only when InputLBIdx0AndLoad == true
//                                       --> Please provide the **sum** of all patches within each patch group
//                                       --> Only rank 0 needs to provide this list
//                                       --> Must be in the same order as LBIdx0_AllRank_Input
//                ParWeight            : Relative load-balance weighting of particles
//                                       --> Weighting of each patch is estimated as "PATCH_SIZE^3 + NParThisPatch*ParWeight"
//                                       --> <= 0.0 : do not consider particle weighting
//
// Return      :  CutPoint
//-------------------------------------------------------------------------------------------------------
void LB_SetCutPoint( const int lv, long *CutPoint, const bool InputLBIdx0AndLoad, long *LBIdx0_AllRank_Input,
                     double *Load_AllRank_Input, const double ParWeight )
{

// check
   if ( NPatchTotal[lv]%8 != 0 )
      Aux_Error( ERROR_INFO, "NPatchTotal[%d] = %d is NOT a multiple of 8 !!\n", lv, NPatchTotal[lv] );

   if ( MPI_Rank == 0  &&  InputLBIdx0AndLoad  &&  ( LBIdx0_AllRank_Input == NULL || Load_AllRank_Input == NULL )  )
      Aux_Error( ERROR_INFO, "LBIdx0_AllRank_Input/Load_AllRank_Input == NULL when InputLBIdx0AndLoad is on !!\n" );


// 1. collect the load-balance weighting and LB_Idx of all patch groups from all ranks
   const int NPG_Total = NPatchTotal[lv] / 8;

   long   *LBIdx0_ThisRank = NULL;
   long   *LBIdx0_AllRank  = NULL;
   double *Load_ThisRank   = NULL;
   double *Load_AllRank    = NULL;
   int    *NPG_AllRank     = NULL;
   int    *Recv_Disp       = NULL;
   int    *IdxTable        = ( MPI_Rank == 0 ) ? new int [NPG_Total] : NULL;

// use the input tables directly
// --> useful during RESTART, where we have very limited information
//     (e.g., we don't know the number of patches in each rank, amr->NPatchComma, and any particle information yet ...)
   if ( InputLBIdx0AndLoad )
   {
      if ( MPI_Rank == 0 )
      {
         LBIdx0_AllRank = LBIdx0_AllRank_Input;
         Load_AllRank   = Load_AllRank_Input;
      }
   }

   else
   {
//    allocate memory
      int NPG_ThisRank = amr->NPatchComma[lv][1] / 8;

      LBIdx0_ThisRank = new long   [ NPG_ThisRank ];
      Load_ThisRank   = new double [ NPG_ThisRank ];

      if ( MPI_Rank == 0 )
      {
         LBIdx0_AllRank = new long   [ NPG_Total ];
         Load_AllRank   = new double [ NPG_Total ];
         NPG_AllRank    = new int    [ MPI_NRank ];
         Recv_Disp      = new int    [ MPI_NRank ];
      }


//    collect the number of patch groups in each rank
      MPI_Gather( &NPG_ThisRank, 1, MPI_INT, NPG_AllRank, 1, MPI_INT, 0, MPI_COMM_WORLD );

      if ( MPI_Rank == 0 )
      {
         Recv_Disp[0] = 0;
         for (int r=0; r<MPI_NRank-1; r++)   Recv_Disp[r+1] = Recv_Disp[r] + NPG_AllRank[r];
      }


//    collect the minimum LBIdx in each patch group
//    --> assuming patches within the same patch group have consecutive LBIdx
      for (int t=0; t<NPG_ThisRank; t++)
      {
         const int PID0 = t*8;

         LBIdx0_ThisRank[t]  = amr->patch[0][lv][PID0]->LB_Idx;
         LBIdx0_ThisRank[t] -= LBIdx0_ThisRank[t] % 8;         // get the **minimum** LBIdx in this patch group
      }

      MPI_Gatherv( LBIdx0_ThisRank, NPG_ThisRank, MPI_LONG, LBIdx0_AllRank, NPG_AllRank, Recv_Disp,
                   MPI_LONG, 0, MPI_COMM_WORLD );


//    collect the load-balance weighting in each patch group
      LB_EstimateWorkload_AllPatchGroup( lv, ParWeight, Load_ThisRank );

      MPI_Gatherv( Load_ThisRank, NPG_ThisRank, MPI_DOUBLE, Load_AllRank, NPG_AllRank, Recv_Disp,
                   MPI_DOUBLE, 0, MPI_COMM_WORLD );

   } // if ( InputLBIdx0AndLoad ) ... else ...


   if ( MPI_Rank == 0 )
   {
      double *Load_Record = ( OPT__VERBOSE ) ? new double [MPI_NRank] : NULL;
      double  Load_Ave;

//    3. sort LB_Idx
//    --> after sorting, we must use IdxTable to access the Load_AllRank[] array
      Mis_Heapsort( NPG_Total, LBIdx0_AllRank, IdxTable );


//    4. set the cut points
//    4-1. take care of the case with no patches at all
      if ( NPG_Total == 0 )
      {
         Load_Ave = 0.0;

         for (int t=0; t<MPI_NRank+1; t++)   CutPoint[t] = -1;

         if ( OPT__VERBOSE )
            for (int r=0; r<MPI_NRank; r++)  Load_Record[r] = 0.0;
      }

      else
      {
//       4-2. get the average workload for each rank
         Load_Ave = 0.0;
         for (int t=0; t<NPG_Total; t++)  Load_Ave += Load_AllRank[t];
         Load_Ave /= (double)MPI_NRank;

//       4-3. set the min and max cut points
         const long LBIdx0_Min = LBIdx0_AllRank[             0 ];
         const long LBIdx0_Max = LBIdx0_AllRank[ NPG_Total - 1 ];
         CutPoint[        0] = LBIdx0_Min;
         CutPoint[MPI_NRank] = LBIdx0_Max + 8;  // +8 since the maximum LBIdx in all patches is LBIdx0_Max + 7

//       4-4. find the LBIdx with an accumulated workload (LoadAcc) closest to the average workload of each rank (LoadTarget)
         int    CutIdx     = 1;                 // target array index for CutPoint[]
                                                // --> note that CutPoint[CutIdx] is the **exclusive** upper bound of rank "CutIdx-1"
         double LoadAcc    = 0.0;               // accumulated workload
         double LoadTarget = CutIdx*Load_Ave;   // target accumulated workload for the rank "CutIdx-1"
         double LoadThisPG;                     // workload of the target patch group

         for (int PG=0; PG<NPG_Total; PG++)
         {
//          nothing to do if all cut points have been set already
            if ( CutIdx == MPI_NRank )    break;

//          remember to use IdxTable to access Load_AllRank
            LoadThisPG = Load_AllRank[ IdxTable[PG] ];

//          check if adding a new patch group will exceed the target accumulated workload
            if ( LoadAcc+LoadThisPG >= LoadTarget )
            {
//             determine the cut point with an accumulated workload **closest** to the target accumulated workload
//             (a) if adding a new patch group will exceed the target accumulated workload too much
//                 --> exclude this patch group from the rank "CutIdx-1"
//             note that both "LoadAcc > LoadTarget" and "LoadAcc <= LoadTaget" can happen
               if ( fabs(LoadAcc-LoadTarget) < LoadAcc+LoadThisPG-LoadTarget )
               {
                  CutPoint[CutIdx] = LBIdx0_AllRank[PG];

                  PG --;   // because this patch group has been **excluded** from this cut point
               }

//             (b) if adding a new patch group will NOT exceed the target accumulated workload too much
//                 --> include this patch group in the rank "CutIdx-1"
               else
               {
//                be careful about the special case "PG == NPG_Total-1"
                  CutPoint[CutIdx] = ( PG == NPG_Total-1 ) ? CutPoint[MPI_NRank] : LBIdx0_AllRank[PG+1];

                  LoadAcc += LoadThisPG;
               }

//             record the **accumulated** wordload of each rank
               if ( OPT__VERBOSE )  Load_Record[ CutIdx - 1 ] = LoadAcc;

               CutIdx ++;
               LoadTarget = CutIdx*Load_Ave;
            } // if ( LoadAcc+LoadThisPG >= LoadTarget )

            else
            {
               LoadAcc += LoadThisPG;
            } // if ( LoadAcc+LoadThisPG >= LoadTarget ) ... else ...

         } // for (int PG=0; PG<NPG_Total; PG++)
      } // if ( NPG_Total == 0 ) ... else ...


//    5. output the cut points and workload of each MPI rank
      if ( OPT__VERBOSE )
      {
//       convert the accumulated workload to the actual workload of each rank
         Load_Record[ MPI_NRank - 1 ] = Load_Ave*MPI_NRank;
         for (int r=MPI_NRank-1; r>=1; r--)  Load_Record[r] -= Load_Record[r-1];

         double Load_Max = -1.0;

         for (int r=0; r<MPI_NRank; r++)
         {
            Aux_Message( stdout, "   Lv %2d: Rank %4d, Cut %15ld -> %15ld, Load_Weighted %9.3e\n",
                         lv, r, CutPoint[r], CutPoint[r+1], Load_Record[r] );

            if ( Load_Record[r] > Load_Max )    Load_Max = Load_Record[r];
         }

         Aux_Message( stdout, "   Load_Ave %9.3e, Load_Max %9.3e --> Load_Imbalance = %6.2f%%\n",
                      Load_Ave, Load_Max, (NPG_Total == 0) ? 0.0 : 100.0*(Load_Max-Load_Ave)/Load_Ave );
         Aux_Message( stdout, "   =============================================================================\n" );

         delete [] Load_Record;
      }
   } // if ( MPI_Rank == 0 )


// 6. broadcast the cut points
   MPI_Bcast( CutPoint, MPI_NRank+1, MPI_LONG, 0, MPI_COMM_WORLD );


// free memory
   if ( MPI_Rank == 0 )    delete [] IdxTable;

   if ( !InputLBIdx0AndLoad )
   {
      delete [] LBIdx0_ThisRank;
      delete [] Load_ThisRank;

      if ( MPI_Rank == 0 )
      {
         delete [] LBIdx0_AllRank;
         delete [] Load_AllRank;
         delete [] NPG_AllRank;
         delete [] Recv_Disp;
      }
   }

} // FUNCTION : LB_SetCutPoint



//-------------------------------------------------------------------------------------------------------
// Function    :  LB_RedistrubteRealPatch
// Description :  Redistribute real patches (and particles) to different ranks according to the cut point
//                array "amr->LB->CutPoint[lv]"
//
// Note        :  1. All ranks must have the array "LB_CutPoint" prepared
//                2. This function assumes that the "patch group" is adopted as the basic unit for data
//                   redistribution
//                3. Real patches with LB_Idx in the range "CutPoint[r] <= LB_Idx < CutPoint[r+1]"
//                   will be sent to rank "r"
//                4. Particles will be redistributed along with the leaf patches as well
//
// Parameter   :  lv : Target refinement level
//-------------------------------------------------------------------------------------------------------
void LB_RedistributeRealPatch( const int lv, real **ParVar_Old, real **Passive_Old )
{

// 1. count the number of real patches (and particles) to be sent and received
// ==========================================================================================
   const int PatchSize1v = CUBE( PATCH_SIZE );
#  ifdef STORE_POT_GHOST
   const int GraNxtSize  = CUBE( GRA_NXT );
#  endif

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

#  ifdef STORE_POT_GHOST
   int *Send_NCount_PotExt = new int [MPI_NRank];
   int *Recv_NCount_PotExt = new int [MPI_NRank];
   int *Send_NDisp_PotExt  = new int [MPI_NRank];
   int *Recv_NDisp_PotExt  = new int [MPI_NRank];
#  endif

#  ifdef PARTICLE
   const int  NParVar           = PAR_NVAR + PAR_NPASSIVE;
   const bool RemoveAllParticle = true;

   int  NSend_Total_ParData, NRecv_Total_ParData;
   long ParID;

   int *Send_NCount_ParData = new int [MPI_NRank];
   int *Recv_NCount_ParData = new int [MPI_NRank];
   int *Send_NDisp_ParData  = new int [MPI_NRank];
   int *Recv_NDisp_ParData  = new int [MPI_NRank];
   int *Counter_ParData     = new int [MPI_NRank];

#  ifdef DEBUG_PARTICLE
   if ( ParVar_Old  == NULL )    Aux_Error( ERROR_INFO, "ParVar_Old == NULL !!\n" );
   if ( Passive_Old == NULL )    Aux_Error( ERROR_INFO, "Passive_Old == NULL !!\n" );
#  endif
#  endif // #ifdef PARTICLE

   for (int r=0; r<MPI_NRank; r++)
   {
      Send_NCount_Patch  [r] = 0;
#     ifdef PARTICLE
      Send_NCount_ParData[r] = 0;
#     endif
   }
   Send_NDisp_Patch  [0] = 0;
   Recv_NDisp_Patch  [0] = 0;
#  ifdef PARTICLE
   Send_NDisp_ParData[0] = 0;
   Recv_NDisp_ParData[0] = 0;
#  endif

// 1.1 send count
   for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
   {
      LB_Idx = amr->patch[0][lv][PID]->LB_Idx;
      TRank  = LB_Index2Rank( lv, LB_Idx, CHECK_ON );

      Send_NCount_Patch  [TRank] ++;
#     ifdef PARTICLE
      Send_NCount_ParData[TRank] += amr->patch[0][lv][PID]->NPar;
#     endif
   }
#  ifdef PARTICLE
   for (int r=0; r<MPI_NRank; r++)  Send_NCount_ParData[r] *= NParVar;
#  endif

// 1.2 receive count
   MPI_Alltoall( Send_NCount_Patch,   1, MPI_INT, Recv_NCount_Patch,   1, MPI_INT, MPI_COMM_WORLD );
#  ifdef PARTICLE
   MPI_Alltoall( Send_NCount_ParData, 1, MPI_INT, Recv_NCount_ParData, 1, MPI_INT, MPI_COMM_WORLD );
#  endif

// 1.3 send/recv displacement
   for (int r=1; r<MPI_NRank; r++)
   {
      Send_NDisp_Patch  [r] = Send_NDisp_Patch  [r-1] + Send_NCount_Patch  [r-1];
      Recv_NDisp_Patch  [r] = Recv_NDisp_Patch  [r-1] + Recv_NCount_Patch  [r-1];
#     ifdef PARTICLE
      Send_NDisp_ParData[r] = Send_NDisp_ParData[r-1] + Send_NCount_ParData[r-1];
      Recv_NDisp_ParData[r] = Recv_NDisp_ParData[r-1] + Recv_NCount_ParData[r-1];
#     endif
   }

// 1.4 send/recv data displacement
   for (int r=0; r<MPI_NRank; r++)
   {
      Send_NCount_Data1v[r] = PatchSize1v*Send_NCount_Patch[r];
      Recv_NCount_Data1v[r] = PatchSize1v*Recv_NCount_Patch[r];
      Send_NDisp_Data1v [r] = PatchSize1v*Send_NDisp_Patch [r];
      Recv_NDisp_Data1v [r] = PatchSize1v*Recv_NDisp_Patch [r];
#     ifdef STORE_POT_GHOST
      Send_NCount_PotExt[r] = GraNxtSize*Send_NCount_Patch[r];
      Recv_NCount_PotExt[r] = GraNxtSize*Recv_NCount_Patch[r];
      Send_NDisp_PotExt [r] = GraNxtSize*Send_NDisp_Patch [r];
      Recv_NDisp_PotExt [r] = GraNxtSize*Recv_NDisp_Patch [r];
#     endif
   }

// 1.5 total number of patches (and particle data) to be sent and received
   NSend_Total_Patch   = Send_NDisp_Patch  [ MPI_NRank-1 ] + Send_NCount_Patch  [ MPI_NRank-1 ];
   NRecv_Total_Patch   = Recv_NDisp_Patch  [ MPI_NRank-1 ] + Recv_NCount_Patch  [ MPI_NRank-1 ];
#  ifdef PARTICLE
   NSend_Total_ParData = Send_NDisp_ParData[ MPI_NRank-1 ] + Send_NCount_ParData[ MPI_NRank-1 ];
   NRecv_Total_ParData = Recv_NDisp_ParData[ MPI_NRank-1 ] + Recv_NCount_ParData[ MPI_NRank-1 ];
#  endif

// 1.6 check
#  ifdef GAMER_DEBUG
   if ( NSend_Total_Patch != amr->NPatchComma[lv][1] )
      Aux_Error( ERROR_INFO, "NSend_Total_Patch (%d) != expected (%d) !!\n", NSend_Total_Patch, amr->NPatchComma[lv][1] );
#  endif
#  ifdef DEBUG_PARTICLE
   if ( NSend_Total_ParData != amr->Par->NPar_Lv[lv]*NParVar )
      Aux_Error( ERROR_INFO, "NSend_Total_ParData (%d) != expected (%ld) !!\n", NSend_Total_ParData, amr->Par->NPar_Lv[lv]*NParVar );
#  endif


// 2. prepare the MPI send buffers
// ==========================================================================================
   const int SendDataSize1v     = NSend_Total_Patch*PatchSize1v;
   const int RecvDataSize1v     = NRecv_Total_Patch*PatchSize1v;
   const int FluSg              = amr->FluSg[lv];
#  ifdef GRAVITY
   const int PotSg              = amr->PotSg[lv];
#  ifdef STORE_POT_GHOST
   const int SendDataSizePotExt = NSend_Total_Patch*GraNxtSize;
   const int RecvDataSizePotExt = NRecv_Total_Patch*GraNxtSize;
#  endif
#  endif

   real *SendPtr         = NULL;
   long *SendBuf_LBIdx   = new long [ NSend_Total_Patch ];
   real *SendBuf_Flu     = new real [ SendDataSize1v*NCOMP_TOTAL ];
#  ifdef GRAVITY
   real *SendBuf_Pot     = new real [ SendDataSize1v ];
#  ifdef STORE_POT_GHOST
   real *SendBuf_PotExt  = new real [ SendDataSizePotExt ];
#  endif
#  endif
#  ifdef PARTICLE
   real *SendBuf_ParData = new real [ NSend_Total_ParData ];
   int  *SendBuf_NPar    = new int  [ NSend_Total_Patch ];
#  endif

   for (int r=0; r<MPI_NRank; r++)
   {
      Counter        [r] = 0;
#     ifdef PARTICLE
      Counter_ParData[r] = 0;
#     endif
   }

   for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
   {
      LB_Idx = amr->patch[0][lv][PID]->LB_Idx;
      TRank  = LB_Index2Rank( lv, LB_Idx, CHECK_ON );

//    2.1 LB_Idx
      SendBuf_LBIdx[ Send_NDisp_Patch[TRank] + Counter[TRank] ] = LB_Idx;

//    2.2 fluid
      for (int v=0; v<NCOMP_TOTAL; v++)
      {
         SendPtr = SendBuf_Flu + v*SendDataSize1v + Send_NDisp_Data1v[TRank] + Counter[TRank]*PatchSize1v;
         memcpy( SendPtr, &amr->patch[FluSg][lv][PID]->fluid[v][0][0][0], PatchSize1v*sizeof(real) );
      }

#     ifdef GRAVITY
//    2.3 potential
      SendPtr = SendBuf_Pot + Send_NDisp_Data1v[TRank] + Counter[TRank]*PatchSize1v;
      memcpy( SendPtr, &amr->patch[PotSg][lv][PID]->pot[0][0][0], PatchSize1v*sizeof(real) );

#     ifdef STORE_POT_GHOST
//    2.4 potential with ghost zones
      SendPtr = SendBuf_PotExt + Send_NDisp_PotExt[TRank] + Counter[TRank]*GraNxtSize;
      memcpy( SendPtr, &amr->patch[PotSg][lv][PID]->pot_ext[0][0][0], GraNxtSize*sizeof(real) );
#     endif
#     endif

//    2.5 particle
#     ifdef PARTICLE
      SendBuf_NPar[ Send_NDisp_Patch[TRank] + Counter[TRank] ] = amr->patch[0][lv][PID]->NPar;

      SendPtr = SendBuf_ParData + Send_NDisp_ParData[TRank] + Counter_ParData[TRank];

      for (int p=0; p<amr->patch[0][lv][PID]->NPar; p++)
      {
         ParID = amr->patch[0][lv][PID]->ParList[p];

//       there should be no inactive particles
#        ifdef DEBUG_PARTICLE
         if ( ParVar_Old[PAR_MASS][ParID] < (real)0.0 )
            Aux_Error( ERROR_INFO, "Mass[%ld] = %14.7e < 0.0 !!\n", ParID, ParVar_Old[PAR_MASS][ParID] );
#        endif

         for (int v=0; v<PAR_NVAR; v++)      *SendPtr++ = ParVar_Old [v][ParID];
         for (int v=0; v<PAR_NPASSIVE; v++)  *SendPtr++ = Passive_Old[v][ParID];
      }
#     endif // #ifdef PARTICLE

      Counter        [TRank] ++;
#     ifdef PARTICLE
      Counter_ParData[TRank] += amr->patch[0][lv][PID]->NPar*NParVar;

//    detach particles from patches to avoid warning messages when deleting patches with particles
      amr->patch[0][lv][PID]->RemoveParticle( NULL_INT, NULL, &amr->Par->NPar_Lv[lv], RemoveAllParticle );
#     endif
   } // for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)

// check if all particles are detached from patches at lv
#  ifdef DEBUG_PARTICLE
   if ( amr->Par->NPar_Lv[lv] != 0 )
      Aux_Error( ERROR_INFO, "NPar_Lv[%d] = %ld != 0 !!\n", lv, amr->Par->NPar_Lv[lv] );
#  endif


// 4. delete old patches and allocate the MPI recv buffers
// ==========================================================================================
// free memory first to reduce the memory consumption
// --> for OPT__REUSE_MEMORY == 2 (aggressive mode), we only mark patches as inactive but do not deallocate memory
// --> NPatchComma is also reset to 0 here
   amr->Lvdelete( lv, OPT__REUSE_MEMORY==2 );

// allocate recv buffers AFTER deleting old patches
   long *RecvBuf_LBIdx   = new long [ NRecv_Total_Patch ];
   real *RecvBuf_Flu     = new real [ RecvDataSize1v*NCOMP_TOTAL ];
#  ifdef GRAVITY
   real *RecvBuf_Pot     = new real [ RecvDataSize1v ];
#  ifdef STORE_POT_GHOST
   real *RecvBuf_PotExt  = new real [ RecvDataSizePotExt ];
#  endif
#  endif
#  ifdef PARTICLE
   real *RecvBuf_ParData = new real [ NRecv_Total_ParData ];
   int  *RecvBuf_NPar    = new int  [ NRecv_Total_Patch ];
#  endif


// 5. transfer data by MPI_Alltoallv
// ==========================================================================================
// 5.1 LB_Idx
   MPI_Alltoallv( SendBuf_LBIdx, Send_NCount_Patch, Send_NDisp_Patch, MPI_LONG,
                  RecvBuf_LBIdx, Recv_NCount_Patch, Recv_NDisp_Patch, MPI_LONG, MPI_COMM_WORLD );

// 5.2 fluid (transfer one component at a time to avoid exceeding the maximum allowed transferred size in MPI)
   for (int v=0; v<NCOMP_TOTAL; v++)
   {
#     ifdef FLOAT8
      MPI_Alltoallv( SendBuf_Flu + v*SendDataSize1v, Send_NCount_Data1v, Send_NDisp_Data1v, MPI_DOUBLE,
                     RecvBuf_Flu + v*RecvDataSize1v, Recv_NCount_Data1v, Recv_NDisp_Data1v, MPI_DOUBLE, MPI_COMM_WORLD );
#     else
      MPI_Alltoallv( SendBuf_Flu + v*SendDataSize1v, Send_NCount_Data1v, Send_NDisp_Data1v, MPI_FLOAT,
                     RecvBuf_Flu + v*RecvDataSize1v, Recv_NCount_Data1v, Recv_NDisp_Data1v, MPI_FLOAT,  MPI_COMM_WORLD );
#     endif
   }

#  ifdef GRAVITY
// 5.3 potential
// --> debugger may report that the potential data are NOT initialized when calling LB_Init_LoadBalance() during initialization
// --> it's fine since we calculate potential AFTER invoking LB_Init_LoadBalance() in Init_GAMER()
#  ifdef FLOAT8
   MPI_Alltoallv( SendBuf_Pot, Send_NCount_Data1v, Send_NDisp_Data1v, MPI_DOUBLE,
                  RecvBuf_Pot, Recv_NCount_Data1v, Recv_NDisp_Data1v, MPI_DOUBLE, MPI_COMM_WORLD );
#  else
   MPI_Alltoallv( SendBuf_Pot, Send_NCount_Data1v, Send_NDisp_Data1v, MPI_FLOAT,
                  RecvBuf_Pot, Recv_NCount_Data1v, Recv_NDisp_Data1v, MPI_FLOAT,  MPI_COMM_WORLD );
#  endif

// 5.4 potential with ghost zones
#  ifdef STORE_POT_GHOST
#  ifdef FLOAT8
   MPI_Alltoallv( SendBuf_PotExt, Send_NCount_PotExt, Send_NDisp_PotExt, MPI_DOUBLE,
                  RecvBuf_PotExt, Recv_NCount_PotExt, Recv_NDisp_PotExt, MPI_DOUBLE, MPI_COMM_WORLD );
#  else
   MPI_Alltoallv( SendBuf_PotExt, Send_NCount_PotExt, Send_NDisp_PotExt, MPI_FLOAT,
                  RecvBuf_PotExt, Recv_NCount_PotExt, Recv_NDisp_PotExt, MPI_FLOAT,  MPI_COMM_WORLD );
#  endif
#  endif // STORE_POT_GHOST
#  endif // GRAVITY

#  ifdef PARTICLE
// 5.5 particle count
   MPI_Alltoallv( SendBuf_NPar, Send_NCount_Patch, Send_NDisp_Patch, MPI_INT,
                  RecvBuf_NPar, Recv_NCount_Patch, Recv_NDisp_Patch, MPI_INT, MPI_COMM_WORLD );

// 5.6 particle data
#  ifdef FLOAT8
   MPI_Alltoallv( SendBuf_ParData, Send_NCount_ParData, Send_NDisp_ParData, MPI_DOUBLE,
                  RecvBuf_ParData, Recv_NCount_ParData, Recv_NDisp_ParData, MPI_DOUBLE, MPI_COMM_WORLD );
#  else
   MPI_Alltoallv( SendBuf_ParData, Send_NCount_ParData, Send_NDisp_ParData, MPI_FLOAT,
                  RecvBuf_ParData, Recv_NCount_ParData, Recv_NDisp_ParData, MPI_FLOAT,  MPI_COMM_WORLD );
#  endif
#  endif // #ifdef PARTICLE


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
#  ifdef STORE_POT_GHOST
   delete [] Send_NCount_PotExt;
   delete [] Recv_NCount_PotExt;
   delete [] Send_NDisp_PotExt;
   delete [] Recv_NDisp_PotExt;
   delete [] SendBuf_PotExt;
#  endif
#  endif // GRAVITY
#  ifdef PARTICLE
   delete [] Send_NCount_ParData;
   delete [] Send_NDisp_ParData;
   delete [] Counter_ParData;
   delete [] SendBuf_ParData;
   delete [] SendBuf_NPar;
#  endif


// 7. reset particle parameters and add particles to the particle repository
// ==========================================================================================
   const real *RecvPtr = NULL;

#  ifdef PARTICLE
   const long NParToBeAdded = NRecv_Total_ParData / NParVar;
   const long NParPrevious  = amr->Par->NPar_AcPlusInac;

#  ifdef DEBUG_PARTICLE
   if ( NParPrevious + NParToBeAdded > amr->Par->ParListSize )
      Aux_Error( ERROR_INFO, "NParExpect (%ld) >= ParListSize (%ld) !!\n",
                 NParPrevious + NParToBeAdded, amr->Par->ParListSize );
#  endif

// add particles to the repository
   RecvPtr = RecvBuf_ParData;
   for (long p=0; p<NParToBeAdded; p++)
   {
      amr->Par->AddOneParticle( RecvPtr, RecvPtr+PAR_NVAR );
      RecvPtr += NParVar;
   }

// free memory
   delete [] Recv_NCount_ParData;
   delete [] Recv_NDisp_ParData;
   delete [] RecvBuf_ParData;
#  endif // #ifdef PARTICLE


// 8. allocate new patches with the data just received (use "patch group" as the basic unit)
// ==========================================================================================
   const int PScale  = PATCH_SIZE*amr->scale[lv];
   const int PGScale = 2*PScale;
   int PID, Cr0[3];

#  ifdef PARTICLE
#  ifdef DEBUG_PARTICLE
   const real *ParPos[3] = { amr->Par->PosX, amr->Par->PosY, amr->Par->PosZ };
#  endif
   long *ParList         = NULL;
   int   ParListSizeMax  = 0;    // must NOT be negative to deal with the case NRecv_Total_Patch == 0

   for (int t=0; t<NRecv_Total_Patch; t++)   ParListSizeMax = MAX( ParListSizeMax, RecvBuf_NPar[t] );

   ParList = new long [ParListSizeMax];
   ParID   = NParPrevious;
#  endif

   for (int PID0=0; PID0<NRecv_Total_Patch; PID0+=8)
   {
      LB_Idx = RecvBuf_LBIdx[PID0];

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
      for (int LocalID=0; LocalID<8; LocalID++)
      {
         PID = PID0 + LocalID;

//       fluid
         for (int v=0; v<NCOMP_TOTAL; v++)
         {
            RecvPtr = RecvBuf_Flu + v*RecvDataSize1v + PID*PatchSize1v;
            memcpy( &amr->patch[FluSg][lv][PID]->fluid[v][0][0][0], RecvPtr, PatchSize1v*sizeof(real) );
         }

#        ifdef GRAVITY
//       potential
         RecvPtr = RecvBuf_Pot + PID*PatchSize1v;
         memcpy( &amr->patch[PotSg][lv][PID]->pot[0][0][0], RecvPtr, PatchSize1v*sizeof(real) );

#        ifdef STORE_POT_GHOST
//       potential with ghost zones
         RecvPtr = RecvBuf_PotExt + PID*GraNxtSize;
         memcpy( &amr->patch[PotSg][lv][PID]->pot_ext[0][0][0], RecvPtr, GraNxtSize*sizeof(real) );
#        endif
#        endif // GRAVITY

//       particle
#        ifdef PARTICLE
         for (int p=0; p<RecvBuf_NPar[PID]; p++)   ParList[p] = ParID++;

#        ifdef DEBUG_PARTICLE
         char Comment[100];
         sprintf( Comment, "%s, PID %d, NPar %d", __FUNCTION__, PID, RecvBuf_NPar[PID] );
         amr->patch[0][lv][PID]->AddParticle( RecvBuf_NPar[PID], ParList, &amr->Par->NPar_Lv[lv],
                                              ParPos, amr->Par->NPar_AcPlusInac, Comment );
#        else
         amr->patch[0][lv][PID]->AddParticle( RecvBuf_NPar[PID], ParList, &amr->Par->NPar_Lv[lv] );
#        endif
#        endif // #ifdef PARTICLE
      } // for (int LocalID=0; LocalID<8; LocalID++)
   } // for (int PID0=0; PID0<NRecv_Total_Patch; PID0+=8)

// reset NPatchComma
   for (int m=1; m<28; m++)   amr->NPatchComma[lv][m] = NRecv_Total_Patch;

// check the amr->NPatchComma recording
   if ( amr->NPatchComma[lv][1] != amr->num[lv] )
      Aux_Error( ERROR_INFO, "amr->NPatchComma[%d][1] (%d) != amr->num[%d] (%d) !!\n",
                 lv, amr->NPatchComma[lv][1], lv, amr->num[lv] );


// 9. record LB_IdxList_Real
// ==========================================================================================
   if ( amr->LB->IdxList_Real         [lv] != NULL )  delete [] amr->LB->IdxList_Real         [lv];
   if ( amr->LB->IdxList_Real_IdxTable[lv] != NULL )  delete [] amr->LB->IdxList_Real_IdxTable[lv];

   amr->LB->IdxList_Real         [lv] = new long [NRecv_Total_Patch];
   amr->LB->IdxList_Real_IdxTable[lv] = new int  [NRecv_Total_Patch];

   for (int PID=0; PID<NRecv_Total_Patch; PID++)   amr->LB->IdxList_Real[lv][PID] = amr->patch[0][lv][PID]->LB_Idx;

   Mis_Heapsort( NRecv_Total_Patch, amr->LB->IdxList_Real[lv], amr->LB->IdxList_Real_IdxTable[lv] );


// 10. deallocate the MPI recv buffers
// ==========================================================================================
   delete [] Recv_NCount_Patch;
   delete [] Recv_NDisp_Patch;
   delete [] Recv_NCount_Data1v;
   delete [] Recv_NDisp_Data1v;
   delete [] RecvBuf_LBIdx;
   delete [] RecvBuf_Flu;
#  ifdef GRAVITY
   delete [] RecvBuf_Pot;
#  ifdef STORE_POT_GHOST
   delete [] RecvBuf_PotExt;
#  endif
#  endif
#  ifdef PARTICLE
   delete [] RecvBuf_NPar;
   delete [] ParList;
#  endif

} // FUNCTION : LB_RedistributePatch



#ifdef PARTICLE
//-------------------------------------------------------------------------------------------------------
// Function    :  LB_RedistributeParticle_Init
// Description :  Initialize the procedure for redistributing particles
//
// Note        :  1. This function will get the total number of particles AFTER data redistribution and
//                   then allocate the new particle attribute arrays
//                2. This function will also reallocate particle repository by calling amr->Par->InitRepo.
//                   However, we reset NPar_AcPlusInac to zero since we will update it level by level
//                   when calling LB_RedistributeRealPatch later
//                3. One must call LB_SetCutPoint for all levels in advance
//
// Parameter   :  ParVar_Old  : Pointers for backing up the old particle attribute arrays (amr->Par->ParVar)
//                PassiveOld  : Pointers for backing up the old particle attribute arrays (amr->Par->Passive)
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void LB_RedistributeParticle_Init( real **ParVar_Old, real **Passive_Old )
{

// backup the old particle attribute arrays
// remember to reset ParVar and Passive to NULL so that amr->Par->InitRepo will NOT delete these arrays
   for (int v=0; v<PAR_NVAR; v++)
   {
      ParVar_Old      [v] = amr->Par->ParVar[v];
      amr->Par->ParVar[v] = NULL;
   }

   for (int v=0; v<PAR_NPASSIVE; v++)
   {
      Passive_Old      [v] = amr->Par->Passive[v];
      amr->Par->Passive[v] = NULL;
   }


// get the total number of particles at each rank after data redistribution
   int  TRank, Send_NPar[MPI_NRank], Recv_NPar[MPI_NRank];
   long LB_Idx, Recv_NPar_Sum;

   for (int r=0; r<MPI_NRank; r++)  Send_NPar[r] = 0;
   Recv_NPar_Sum = 0;

   for (int lv=0; lv<NLEVEL; lv++)
   for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
   {
      LB_Idx            = amr->patch[0][lv][PID]->LB_Idx;
      TRank             = LB_Index2Rank( lv, LB_Idx, CHECK_ON );
      Send_NPar[TRank] += amr->patch[0][lv][PID]->NPar;
   }

   MPI_Alltoall( Send_NPar, 1, MPI_INT, Recv_NPar, 1, MPI_INT, MPI_COMM_WORLD );

   for (int r=0; r<MPI_NRank; r++)  Recv_NPar_Sum += Recv_NPar[r];


// reset particle variables (do not reset NPar_Lv since we will need it for debug in LB_RedistributeRealPatch)
   amr->Par->InitRepo( Recv_NPar_Sum, MPI_NRank );

// reset the total number of particles to be zero
// --> so particle repository is pre-allocated, but it contains no active particle yet
// --> we will add active particles in LB_RedistributeRealPatch
   amr->Par->NPar_AcPlusInac = 0;
   amr->Par->NPar_Active     = 0;

} // FUNCTION : LB_RedistributeParticle_Init



//-------------------------------------------------------------------------------------------------------
// Function    :  LB_RedistributeParticle_End
// Description :  End the procedure for redistributing particles
//
// Note        :  1. Free old particle attribute arrays
//
// Parameter   :  ParVar_Old  : Pointers for backing up the old particle attribute arrays (amr->Par->ParVar)
//                PassiveOld  : Pointers for backing up the old particle attribute arrays (amr->Par->Passive)
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void LB_RedistributeParticle_End( real **ParVar_Old, real **Passive_Old )
{

// remove old particle attribute arrays
   for (int v=0; v<PAR_NVAR;     v++)  free( ParVar_Old [v] );
   for (int v=0; v<PAR_NPASSIVE; v++)  free( Passive_Old[v] );


// check the total number of particles
   if ( amr->Par->NPar_AcPlusInac != amr->Par->NPar_Active )
      Aux_Error( ERROR_INFO, "NPar_AcPlusInac (%ld) != NPar_Active (%ld) !!\n", amr->Par->NPar_AcPlusInac, amr->Par->NPar_Active );

   long NPar_Lv_Sum=0;
   for (int lv=0; lv<NLEVEL; lv++)  NPar_Lv_Sum += amr->Par->NPar_Lv[lv];

   if ( NPar_Lv_Sum != amr->Par->NPar_Active )
      Aux_Error( ERROR_INFO, "NPar_Lv_Sum (%ld) != expect (%ld) !!\n", NPar_Lv_Sum, amr->Par->NPar_Active );

#  ifdef DEBUG_PARTICLE
   long NPar_Active_AllRank_Check;

   MPI_Reduce( &amr->Par->NPar_Active, &NPar_Active_AllRank_Check, 1, MPI_LONG, MPI_SUM, 0, MPI_COMM_WORLD );

   if ( MPI_Rank == 0  &&  NPar_Active_AllRank_Check != amr->Par->NPar_Active_AllRank )
      Aux_Error( ERROR_INFO, "NPar_Active_AllRank (%ld) != expected (%ld) !!\n",
                 NPar_Active_AllRank_Check, amr->Par->NPar_Active_AllRank );
#  endif

} // FUNCTION : LB_RedistributeParticle_End
#endif // #ifdef PARTICLE



#endif // #ifdef LOAD_BALANCE
