#include "GAMER.h"

#ifdef LOAD_BALANCE



static void LB_RedistributeRealPatch( const int lv, real_par **ParAtt_Old, const bool RemoveParFromRepo, const bool SendGridData );
static void LB_SortRealPatch( const int lv );
#ifdef PARTICLE
static void LB_RedistributeParticle_Init( real_par **ParAtt_Old );
static void LB_RedistributeParticle_End ( real_par **ParAtt_Old );
#endif




//-------------------------------------------------------------------------------------------------------
// Function    :  LB_Init_LoadBalance
// Description :  Initialize the load-balance process
//
// Note        :  1. Redistribute patches and reconstruct patch relation on the target level(s)
//                2. All parameters in the data structure "LB_t LB" will be reconstructed on the target level(s)
//                3. Data structure "ParaVar_t ParaVar" will be removed
//                4. This function is used in both initialization phase and run-time data redistribution
//                5. Data in the buffer patches will be filled up
//                6. NPatchTotal[] and NPatchComma[] must be prepared in advance
//                7. Particles will also be redistributed
//
// Parameter   :  Redistribute  : true  --> Redistribute all real patches according to the load-balance weighting of
//                                          each patch and initialize all load-balance related set-up
//                                false --> Initialize all load-balance related set-up, but do NOT invoke LB_SetCutPoint()
//                                          and LB_RedistributeRealPatch() to redistribute all real patches
//                                      --> Currently it is used only during the RESTART process since we already call
//                                          LB_SetCutPoint() and load real patches accordingly when calling Init_ByRestart_*()
//                SendGridData  : Transfer grid data
//                                --> Set to false by Init_ByFile() to reduce memory consumption
//                                --> Useless when Redistribute==false
//                ParWeight     : Relative load-balance weighting of particles
//                                --> Weighting of each patch is estimated as "PATCH_SIZE^3 + NParThisPatch*ParWeight"
//                                --> <= 0.0 : do not consider particle weighting
//                                             --> Currently we force ParWeight==0.0 when calling LB_Init_LoadBalance()
//                                                 for the first time during the restart process since we don't have enough
//                                                 information for calculating particle weighting at that time
//                                             --> For example, Par_LB_CollectParticle2OneLevel() invoked by
//                                                 LB_EstimateWorkload_AllPatchGroup() needs amr->LB->IdxList_Real[], which
//                                                 will be constructed only AFTER calling LB_Init_LoadBalance()
//                Reset         : Call LB->reset() to reset the load-balance variables on the target level(s)
//                                --> Note that CutPoint[] will NOT be reset even when "Reset == true"
//                SortRealPatch : Sort real patches by load-balance indices to improve bitwise reproducibility during restart
//                                --> Ensure the order of real patches in the snapshot and after restart is the same when using
//                                    the same number of MPI processes
//                                --> Controlled by the runtime parameter OPT__SORT_PATCH_BY_LBIDX
//                                --> Work even for Redistribute==false
//                TLv           : Target refinement level(s)
//                                --> 0~TOP_LEVEL : only apply to a specific level
//                                    <0          : apply to all levels
//-------------------------------------------------------------------------------------------------------
void LB_Init_LoadBalance( const bool Redistribute, const bool SendGridData, const double ParWeight, const bool Reset,
                          const bool SortRealPatch, const int TLv )
{

   if ( MPI_Rank == 0 )
   {
      char lv_str[MAX_STRING];
      if ( TLv < 0 )    sprintf( lv_str, "%s", "all levels" );
      else              sprintf( lv_str, "Lv %d", TLv );

      Aux_Message( stdout, "   %s at %s ...\n", __FUNCTION__, lv_str );
   }


// check
   if ( amr->LB == NULL )  Aux_Error( ERROR_INFO, "amr->LB has not been allocated !!\n" );

   if ( TLv > TOP_LEVEL )  Aux_Error( ERROR_INFO, "TLv (%d) > TOP_LEVEL (%d) !!\n", TLv, TOP_LEVEL );

// check the synchronization
   if ( TLv < 0 )
   for (int lv=1; lv<NLEVEL; lv++)
      if ( NPatchTotal[lv] != 0 )   Mis_CompareRealValue( Time[0], Time[lv], __FUNCTION__, true );


// delete ParaVar which is no longer useful
   if ( amr->ParaVar != NULL )
   {
      delete amr->ParaVar;
      amr->ParaVar = NULL;
   }


// 0. set the target level(s)
   const int lv_min = ( TLv < 0 ) ?         0 : TLv;
   const int lv_max = ( TLv < 0 ) ? TOP_LEVEL : TLv;


// 1. set up the load-balance cut points (must do this before calling LB_RedistributeParticle_Init())
   const bool InputLBIdxAndLoad_No = false;

   if ( Redistribute )
   for (int lv=lv_min; lv<=lv_max; lv++)
      LB_SetCutPoint( lv, NPatchTotal[lv]/8, amr->LB->CutPoint[lv], InputLBIdxAndLoad_No, NULL, NULL, ParWeight );


// 2. reinitialize arrays used by the load-balance routines
//    --> must do this AFTER calling LB_SetCutPoint() since it still needs to access load-balance information when
//        PARTICLE is on
//        --> for example, Par_LB_CollectParticle2OneLevel() invoked by LB_EstimateWorkload_AllPatchGroup()
//            will access amr->LB->IdxList_Real[], which will be reset when calling amr->LB->reset()
//    --> amr->LB->reset() will NOT reset CutPoint[]
//        --> otherwise it will just overwrite the cut points set by LB_SetCutPoint()
   if ( Reset )
   for (int lv=lv_min; lv<=lv_max; lv++)  amr->LB->reset( lv );


// 3. re-distribute and allocate all patches (and their associated particles)
   const bool RemoveParFromRepo_Yes = true;
   const bool RemoveParFromRepo_No  = false;

#  ifdef PARTICLE
   real_par  *ParAtt_Old[PAR_NATT_TOTAL];
#  else
   real_par **ParAtt_Old = NULL;
#  endif

#  ifdef PARTICLE
   if ( Redistribute )
   {
      if ( TLv < 0 )
         LB_RedistributeParticle_Init( ParAtt_Old );

      else
      {
         for (int v=0; v<PAR_NATT_TOTAL; v++)   ParAtt_Old[v] = amr->Par->Attribute[v];
      }
   }

   else
   {
      for (int v=0; v<PAR_NATT_TOTAL; v++)   ParAtt_Old[v] = NULL;
   }
#  endif

   for (int lv=lv_min; lv<=lv_max; lv++)
   {
      if ( OPT__VERBOSE  &&  MPI_Rank == 0 )    Aux_Message( stdout, "      Re-distributing patches at Lv %2d ... ", lv );

//    3.1 re-distribute real patches (and particles)
      if ( Redistribute )
      LB_RedistributeRealPatch( lv, ParAtt_Old, (TLv<0)?RemoveParFromRepo_No:RemoveParFromRepo_Yes, SendGridData );

//    3.2 sort real patches
      if ( SortRealPatch )
      LB_SortRealPatch( lv );

//    3.3 allocate sibling-buffer patches at lv
      LB_AllocateBufferPatch_Sibling( lv );

//    3.4 allocate father-buffer patches at lv
//        --> only necessary when applying LB_Init_LoadBalance() to a single level
      if ( TLv >= 0  &&  lv < TOP_LEVEL )
      LB_AllocateBufferPatch_Father( lv+1, true, NULL_INT, NULL, false, NULL, NULL );

//    3.5 allocate father-buffer patches at lv-1
      if ( lv > 0 )
      LB_AllocateBufferPatch_Father( lv,   true, NULL_INT, NULL, false, NULL, NULL );

      if ( OPT__VERBOSE  &&  MPI_Rank == 0 )    Aux_Message( stdout, "done\n" );
   } // for (int lv=lv_min; lv<=lv_max; lv++)

#  ifdef PARTICLE
   if ( Redistribute  &&  TLv < 0 )    LB_RedistributeParticle_End( ParAtt_Old );
#  endif


// 4. contruct the patch relation
   const bool ResetSonID_Yes = true;

   for (int lv=lv_min; lv<=lv_max; lv++)
   {
      if ( OPT__VERBOSE  &&  MPI_Rank == 0 )    Aux_Message( stdout, "      Constructing patch relation at Lv %2d ... ", lv );

//    4.1 sibling relation at lv
      LB_SiblingSearch( lv, true, NULL_INT, NULL );

//    4.2 father-son relation between lv and lv-1
//        --> must reset the son indices on lv-1 as -1 (i.e., using ResetSonID_Yes) when applying
//            LB_Init_LoadBalance() to a single level since, after redistributing patches on lv,
//            the sons of some patches on lv-1 may move abroad (and thus their son indices should
//            be changed from >=0 to <=SON_OFFSET_LB)
//        --> but LB_FindSonNotHome() only checks patches with son indices <= -1
      if ( lv > 0 )
      {
         LB_FindFather    ( lv,   true, NULL_INT, NULL, ResetSonID_Yes );
         LB_FindSonNotHome( lv-1, true, NULL_INT, NULL );
      }

//    reconstruct the following relation only when applying LB_Init_LoadBalance() to a single level
      if ( TLv >= 0 )
      {
//       4.3 sibling relation at lv-1
         if ( lv > 0 )  LB_SiblingSearch( lv-1, true, NULL_INT, NULL );

//       4.4 father-son relation between lv+1 and lv
         if ( lv < TOP_LEVEL )
         {
            LB_FindFather    ( lv+1, true, NULL_INT, NULL, ResetSonID_Yes );
            LB_FindSonNotHome( lv,   true, NULL_INT, NULL );
         }
      }

      if ( OPT__VERBOSE  &&  MPI_Rank == 0 )    Aux_Message( stdout, "done\n" );
   } // for (int lv=lv_min; lv<=lv_max; lv++)


// 5. construct the MPI send and recv data list
   const int lv_min_mpi = ( TLv < 0 ) ?         0 : MAX(0,TLv-1);
   const int lv_max_mpi = ( TLv < 0 ) ? TOP_LEVEL : TLv;

   for (int lv=lv_min_mpi; lv<=lv_max_mpi; lv++)
   {
      if ( OPT__VERBOSE  &&  MPI_Rank == 0 )    Aux_Message( stdout, "      Constructing MPI lists at Lv %2d ... ", lv );

//    5.1 list for exchanging hydro and potential data
      LB_RecordExchangeDataPatchID( lv, false );

//    5.2 list for exchanging restricted hydro data
//        --> note that even when OPT__FIXUP_RESTRICT is off we still need to do data restriction in several places
//            (e.g., restart and OPT__CORR_AFTER_ALL_SYNC)
//        --> for simplicity and sustainability, we always invoke LB_RecordExchangeRestrictDataPatchID()
      LB_RecordExchangeRestrictDataPatchID( lv );

//    5.3 list for exchanging hydro fluxes (and also allocate flux arrays)
      if ( amr->WithFlux )
      LB_AllocateFluxArray( lv );

//    5.4 list for exchanging MHD electric field (and also allocate electric field arrays)
#     ifdef MHD
      if ( amr->WithElectric )
      MHD_LB_AllocateElectricArray( lv );
#     endif

//    5.5 list for exchanging hydro data after the fix-up operation
//        --> for simplicity and sustainability, we always invoke LB_RecordExchangeFixUpDataPatchID()
//        --> see the comments 5.2 above
      LB_RecordExchangeFixUpDataPatchID( lv );

//    5.6 list for overlapping MPI time with CPU/GPU computation
      if ( OPT__OVERLAP_MPI )
      LB_RecordOverlapMPIPatchID( lv );

//    5.7 list for exchanging particles
#     ifdef PARTICLE
      Par_LB_RecordExchangeParticlePatchID( lv );
#     endif

      if ( OPT__VERBOSE  &&  MPI_Rank == 0 )    Aux_Message( stdout, "done\n" );
   } // for (int lv=lv_min_mpi; lv<=lv_max_mpi; lv++)

// 5.8 list for exchanging particles on TLv+1
#  ifdef PARTICLE
   if ( TLv >= 0  &&  TLv < TOP_LEVEL )
   Par_LB_RecordExchangeParticlePatchID( TLv+1 );
#  endif


// 6. get the buffer data
   for (int lv=lv_min_mpi; lv<=lv_max_mpi; lv++)
   {
      if ( OPT__VERBOSE  &&  MPI_Rank == 0 )    Aux_Message( stdout, "      Transferring buffer data at Lv %2d ... ", lv );

      Buf_GetBufferData( lv, amr->FluSg[lv], amr->MagSg[lv], NULL_INT, DATA_GENERAL, _TOTAL, _MAG, Flu_ParaBuf, USELB_YES );

#     ifdef GRAVITY
      Buf_GetBufferData( lv, NULL_INT, NULL_INT, amr->PotSg[lv], POT_FOR_POISSON, _POTE, _NONE, Pot_ParaBuf, USELB_YES );
#     endif

      if ( OPT__VERBOSE  &&  MPI_Rank == 0 )    Aux_Message( stdout, "done\n" );
   }


// 7. reset *SgTime[lv][ 1-*Sg[lv] ] to an arbitrary "negative" number to indicate that the
//    data at 1-Sg are no longer available
   if ( Redistribute )
   for (int lv=lv_min; lv<=lv_max; lv++)
   {
      amr->FluSgTime[lv][ 1-amr->FluSg[lv] ] = -__FLT_MAX__;
#     ifdef MHD
      amr->MagSgTime[lv][ 1-amr->MagSg[lv] ] = -__FLT_MAX__;
#     endif
#     ifdef GRAVITY
      amr->PotSgTime[lv][ 1-amr->PotSg[lv] ] = -__FLT_MAX__;
#     endif
   }


   if ( MPI_Rank == 0 )
   {
      char lv_str[MAX_STRING];
      if ( TLv < 0 )    sprintf( lv_str, "%s", "all levels" );
      else              sprintf( lv_str, "Lv %d", TLv );

      Aux_Message( stdout, "   %s at %s ... done\n", __FUNCTION__, lv_str );
   }

} // FUNCTION : LB_Init_LoadBalance



//-------------------------------------------------------------------------------------------------------
// Function    :  LB_RedistributeRealPatch
// Description :  Redistribute real patches (and particles) to different ranks according to the cut point
//                array amr->LB->CutPoint[]
//
// Note        :  1. All ranks must have LB_CutPoint[] prepared in advance
//                2. This function adopts the "patch group" as the basic unit for data redistribution
//                3. Real patches with LB_Idx in the range "CutPoint[lv][r] <= LB_Idx < CutPoint[lv][r+1]"
//                   will be sent to rank "r"
//                4. Particles will be redistributed along with the leaf patches as well
//
// Parameter   :  lv                : Target refinement level
//                ParAtt_Old        : Pointers pointing to the particle attribute arrays (amr->Par->Attribute[])
//                RemoveParFromRepo : Remove particles on lv from the particle repository (amr->Par)
//                                    --> Useful when applying LB_Init_LoadBalance() to a single level (i.e., TLv>=0)
//                SendGridData      : Transfer grid data
//                                    --> Particle data will always be transferred
//-------------------------------------------------------------------------------------------------------
void LB_RedistributeRealPatch( const int lv, real_par **ParAtt_Old, const bool RemoveParFromRepo, const bool SendGridData )
{

// 1. count the number of real patches (and particles) to be sent and received
// ==========================================================================================
   const int FluSize1v  = CUBE( PS1 );
#  ifdef STORE_POT_GHOST
   const int GraNxtSize = CUBE( GRA_NXT );
#  endif
#  ifdef MHD
   const int MagSize1v  = PS1P1*SQR( PS1 );
#  endif

   int  NSend_Total_Patch, NRecv_Total_Patch, TRank;
   long LB_Idx;

   int *Send_NCount_Patch    = new int  [MPI_NRank];
   int *Recv_NCount_Patch    = new int  [MPI_NRank];
   int *Send_NDisp_Patch     = new int  [MPI_NRank];
   int *Recv_NDisp_Patch     = new int  [MPI_NRank];
   int *NDone_Patch          = new int  [MPI_NRank];
   long *Send_NCount_Flu1v   = new long [MPI_NRank];
   long *Recv_NCount_Flu1v   = new long [MPI_NRank];
   long *Send_NDisp_Flu1v    = new long [MPI_NRank];
   long *Recv_NDisp_Flu1v    = new long [MPI_NRank];

#  ifdef STORE_POT_GHOST
   long *Send_NCount_PotExt  = new long [MPI_NRank];
   long *Recv_NCount_PotExt  = new long [MPI_NRank];
   long *Send_NDisp_PotExt   = new long [MPI_NRank];
   long *Recv_NDisp_PotExt   = new long [MPI_NRank];
#  endif

#  ifdef MHD
   long *Send_NCount_Mag1v   = new long [MPI_NRank];
   long *Recv_NCount_Mag1v   = new long [MPI_NRank];
   long *Send_NDisp_Mag1v    = new long [MPI_NRank];
   long *Recv_NDisp_Mag1v    = new long [MPI_NRank];
#  endif

#  ifdef PARTICLE
   const bool RemoveAllParticle = true;

   long NSend_Total_ParData, NRecv_Total_ParData;
   long ParID;

   long *NDone_ParData       = new long [MPI_NRank];
   long *Send_NCount_ParData = new long [MPI_NRank];
   long *Recv_NCount_ParData = new long [MPI_NRank];
   long *Send_NDisp_ParData  = new long [MPI_NRank];
   long *Recv_NDisp_ParData  = new long [MPI_NRank];

#  ifdef DEBUG_PARTICLE
   if ( ParAtt_Old  == NULL )    Aux_Error( ERROR_INFO, "ParAtt_Old == NULL !!\n" );
#  endif
#  endif // #ifdef PARTICLE

   for (int r=0; r<MPI_NRank; r++)
   {
      Send_NCount_Patch  [r] = 0;
#     ifdef PARTICLE
      Send_NCount_ParData[r] = 0L;
#     endif
   }
   Send_NDisp_Patch  [0] = 0;
   Recv_NDisp_Patch  [0] = 0;
#  ifdef PARTICLE
   Send_NDisp_ParData[0] = 0L;
   Recv_NDisp_ParData[0] = 0L;
#  endif

// 1.1 send count
   for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
   {
      LB_Idx = amr->patch[0][lv][PID]->LB_Idx;
      TRank  = LB_Index2Rank( lv, LB_Idx, CHECK_ON );

      Send_NCount_Patch  [TRank] ++;
#     ifdef PARTICLE
      Send_NCount_ParData[TRank] += (long)amr->patch[0][lv][PID]->NPar;
#     endif
   }
#  ifdef PARTICLE
   for (int r=0; r<MPI_NRank; r++)  Send_NCount_ParData[r] *= (long)PAR_NATT_TOTAL;
#  endif

// 1.2 receive count
   MPI_Alltoall( Send_NCount_Patch,   1, MPI_INT,  Recv_NCount_Patch,   1, MPI_INT,  MPI_COMM_WORLD );
#  ifdef PARTICLE
   MPI_Alltoall( Send_NCount_ParData, 1, MPI_LONG, Recv_NCount_ParData, 1, MPI_LONG, MPI_COMM_WORLD );
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
      Send_NCount_Flu1v [r] = (long)FluSize1v  * (long)Send_NCount_Patch[r];
      Recv_NCount_Flu1v [r] = (long)FluSize1v  * (long)Recv_NCount_Patch[r];
      Send_NDisp_Flu1v  [r] = (long)FluSize1v  * (long)Send_NDisp_Patch [r];
      Recv_NDisp_Flu1v  [r] = (long)FluSize1v  * (long)Recv_NDisp_Patch [r];
#     ifdef STORE_POT_GHOST
      Send_NCount_PotExt[r] = (long)GraNxtSize * (long)Send_NCount_Patch[r];
      Recv_NCount_PotExt[r] = (long)GraNxtSize * (long)Recv_NCount_Patch[r];
      Send_NDisp_PotExt [r] = (long)GraNxtSize * (long)Send_NDisp_Patch [r];
      Recv_NDisp_PotExt [r] = (long)GraNxtSize * (long)Recv_NDisp_Patch [r];
#     endif
#     ifdef MHD
      Send_NCount_Mag1v [r] = (long)MagSize1v  * (long)Send_NCount_Patch[r];
      Recv_NCount_Mag1v [r] = (long)MagSize1v  * (long)Recv_NCount_Patch[r];
      Send_NDisp_Mag1v  [r] = (long)MagSize1v  * (long)Send_NDisp_Patch [r];
      Recv_NDisp_Mag1v  [r] = (long)MagSize1v  * (long)Recv_NDisp_Patch [r];
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
      Aux_Error( ERROR_INFO, "NSend_Total_Patch (%d) != expected (%d) !!\n",
                 NSend_Total_Patch, amr->NPatchComma[lv][1] );
#  endif
#  ifdef DEBUG_PARTICLE
   if ( NSend_Total_ParData != (long)amr->Par->NPar_Lv[lv]*(long)PAR_NATT_TOTAL )
      Aux_Error( ERROR_INFO, "NSend_Total_ParData (%ld) != expected (%ld) !!\n",
                 NSend_Total_ParData, (long)amr->Par->NPar_Lv[lv]*(long)PAR_NATT_TOTAL );
#  endif


// 2. prepare the MPI send buffers
// ==========================================================================================
   const long SendDataSizeFlu1v  = NSend_Total_Patch*FluSize1v;
   const long RecvDataSizeFlu1v  = NRecv_Total_Patch*FluSize1v;
   const int  FluSg              = amr->FluSg[lv];
#  ifdef GRAVITY
   const int  PotSg              = amr->PotSg[lv];
#  ifdef STORE_POT_GHOST
   const long SendDataSizePotExt = NSend_Total_Patch*GraNxtSize;
   const long RecvDataSizePotExt = NRecv_Total_Patch*GraNxtSize;
#  endif
#  endif // GRAVITY
#  ifdef MHD
   const int  MagSg              = amr->MagSg[lv];
   const long SendDataSizeMag1v  = NSend_Total_Patch*MagSize1v;
   const long RecvDataSizeMag1v  = NRecv_Total_Patch*MagSize1v;
#  endif

   real     *SendPtr         = NULL;
   real_par *SendPtr_Par     = NULL;
   long     *SendBuf_LBIdx   = new long [ NSend_Total_Patch ];
   real     *SendBuf_Flu     = ( SendGridData ) ? new real [ SendDataSizeFlu1v*NCOMP_TOTAL ] : NULL;
#  ifdef GRAVITY
   real     *SendBuf_Pot     = ( SendGridData ) ? new real [ SendDataSizeFlu1v ]             : NULL;
#  ifdef STORE_POT_GHOST
   real     *SendBuf_PotExt  = ( SendGridData ) ? new real [ SendDataSizePotExt ]            : NULL;
#  endif
#  endif // GRAVITY
#  ifdef MHD
   real     *SendBuf_Mag     = ( SendGridData ) ? new real [ SendDataSizeMag1v*NCOMP_MAG ]   : NULL;
#  endif
#  ifdef PARTICLE
   real_par *SendBuf_ParData = new real_par [ NSend_Total_ParData ];
   int      *SendBuf_NPar    = new int      [ NSend_Total_Patch ];
#  endif

   for (int r=0; r<MPI_NRank; r++)
   {
      NDone_Patch  [r] = 0;
#     ifdef PARTICLE
      NDone_ParData[r] = 0L;
#     endif
   }

   for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
   {
      LB_Idx = amr->patch[0][lv][PID]->LB_Idx;
      TRank  = LB_Index2Rank( lv, LB_Idx, CHECK_ON );

//    2.1 LB_Idx
      SendBuf_LBIdx[ Send_NDisp_Patch[TRank] + NDone_Patch[TRank] ] = LB_Idx;

      if ( SendGridData )
      {
//       2.2 fluid
         for (int v=0; v<NCOMP_TOTAL; v++)
         {
            SendPtr = SendBuf_Flu + v*SendDataSizeFlu1v + Send_NDisp_Flu1v[TRank] + (long)NDone_Patch[TRank]*FluSize1v;
            memcpy( SendPtr, &amr->patch[FluSg][lv][PID]->fluid[v][0][0][0], FluSize1v*sizeof(real) );
         }

#        ifdef GRAVITY
//       2.3 potential
         SendPtr = SendBuf_Pot + Send_NDisp_Flu1v[TRank] + (long)NDone_Patch[TRank]*FluSize1v;
         memcpy( SendPtr, &amr->patch[PotSg][lv][PID]->pot[0][0][0], FluSize1v*sizeof(real) );

//       2.4 potential with ghost zones
#        ifdef STORE_POT_GHOST
         SendPtr = SendBuf_PotExt + Send_NDisp_PotExt[TRank] + (long)NDone_Patch[TRank]*GraNxtSize;
         memcpy( SendPtr, &amr->patch[PotSg][lv][PID]->pot_ext[0][0][0], GraNxtSize*sizeof(real) );
#        endif
#        endif

//       2.5 magnetic field
#        ifdef MHD
         for (int v=0; v<NCOMP_MAG; v++)
         {
            SendPtr = SendBuf_Mag + v*SendDataSizeMag1v + Send_NDisp_Mag1v[TRank] + (long)NDone_Patch[TRank]*MagSize1v;
            memcpy( SendPtr, &amr->patch[MagSg][lv][PID]->magnetic[v][0], MagSize1v*sizeof(real) );
         }
#        endif
      } // if ( SendGridData )

//    2.6 particle
#     ifdef PARTICLE
      SendBuf_NPar[ Send_NDisp_Patch[TRank] + NDone_Patch[TRank] ] = amr->patch[0][lv][PID]->NPar;

      SendPtr_Par = SendBuf_ParData + Send_NDisp_ParData[TRank] + NDone_ParData[TRank];

      for (int p=0; p<amr->patch[0][lv][PID]->NPar; p++)
      {
         ParID = amr->patch[0][lv][PID]->ParList[p];

//       there should be no inactive particles associated with patches
#        ifdef DEBUG_PARTICLE
         if ( ParAtt_Old[PAR_MASS][ParID] < (real_par)0.0 )
            Aux_Error( ERROR_INFO, "Mass[%ld] = %14.7e < 0.0 !!\n", ParID, ParAtt_Old[PAR_MASS][ParID] );
#        endif

         for (int v=0; v<PAR_NATT_TOTAL; v++)   *SendPtr_Par++ = ParAtt_Old[v][ParID];

//       remove this particle from the particle repository
         if ( RemoveParFromRepo )   amr->Par->RemoveOneParticle( ParID, PAR_INACTIVE_MPI );
      }
#     endif // #ifdef PARTICLE

      NDone_Patch  [TRank] ++;
#     ifdef PARTICLE
      NDone_ParData[TRank] += (long)amr->patch[0][lv][PID]->NPar*(long)PAR_NATT_TOTAL;

//    detach particles from patches to avoid warning messages when deleting
//    patches with particles
      const real_par *PType = amr->Par->Type;
      amr->patch[0][lv][PID]->RemoveParticle( NULL_INT, NULL, &amr->Par->NPar_Lv[lv], RemoveAllParticle, PType );
#     endif
   } // for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)

// check if all particles are detached from patches at lv
#  ifdef DEBUG_PARTICLE
   if ( amr->Par->NPar_Lv[lv] != 0 )
      Aux_Error( ERROR_INFO, "NPar_Lv[%d] = %ld != 0 !!\n", lv, amr->Par->NPar_Lv[lv] );
#  endif


// 3. delete old patches and allocate the MPI recv buffers
// ==========================================================================================
// free memory first to reduce the memory consumption
// --> for OPT__REUSE_MEMORY == 2 (aggressive mode), we only mark patches as inactive but do not deallocate memory
// --> NPatchComma is also reset to 0 here
   amr->Lvdelete( lv, OPT__REUSE_MEMORY==2 );

// allocate recv buffers AFTER deleting old patches
   long *RecvBuf_LBIdx   = new long [ NRecv_Total_Patch ];
   real *RecvBuf_Flu     = ( SendGridData ) ? new real [ RecvDataSizeFlu1v*NCOMP_TOTAL ] : NULL;
#  ifdef GRAVITY
   real *RecvBuf_Pot     = ( SendGridData ) ? new real [ RecvDataSizeFlu1v ]             : NULL;
#  ifdef STORE_POT_GHOST
   real *RecvBuf_PotExt  = ( SendGridData ) ? new real [ RecvDataSizePotExt ]            : NULL;
#  endif
#  endif // GRAVITY
#  ifdef MHD
   real *RecvBuf_Mag     = ( SendGridData ) ? new real [ RecvDataSizeMag1v*NCOMP_MAG ]   : NULL;
#  endif
#  ifdef PARTICLE
   real_par *RecvBuf_ParData = new real_par [ NRecv_Total_ParData ];
   int      *RecvBuf_NPar    = new int      [ NRecv_Total_Patch ];
#  endif


// 4. transfer data by MPI_Alltoallv
// ==========================================================================================
// 4.1 LB_Idx
   MPI_Alltoallv( SendBuf_LBIdx, Send_NCount_Patch, Send_NDisp_Patch, MPI_LONG,
                  RecvBuf_LBIdx, Recv_NCount_Patch, Recv_NDisp_Patch, MPI_LONG, MPI_COMM_WORLD );

   if ( SendGridData )
   {
//    4.2 fluid (transfer one component at a time to avoid exceeding the maximum allowed transfer size in MPI)
      for (int v=0; v<NCOMP_TOTAL; v++)
      {
         MPI_Alltoallv_GAMER( SendBuf_Flu + v*SendDataSizeFlu1v, Send_NCount_Flu1v, Send_NDisp_Flu1v, MPI_GAMER_REAL,
                              RecvBuf_Flu + v*RecvDataSizeFlu1v, Recv_NCount_Flu1v, Recv_NDisp_Flu1v, MPI_GAMER_REAL, MPI_COMM_WORLD );
      }

#     ifdef GRAVITY
//    4.3 potential
//    --> debugger may report that the potential data are NOT initialized when calling LB_Init_LoadBalance()
//        during initialization
//    --> it's fine since we will calculate potential AFTER invoking LB_Init_LoadBalance() in Init_GAMER()
      MPI_Alltoallv_GAMER( SendBuf_Pot, Send_NCount_Flu1v, Send_NDisp_Flu1v, MPI_GAMER_REAL,
                           RecvBuf_Pot, Recv_NCount_Flu1v, Recv_NDisp_Flu1v, MPI_GAMER_REAL, MPI_COMM_WORLD );

//    4.4 potential with ghost zones
#     ifdef STORE_POT_GHOST
      MPI_Alltoallv_GAMER( SendBuf_PotExt, Send_NCount_PotExt, Send_NDisp_PotExt, MPI_GAMER_REAL,
                           RecvBuf_PotExt, Recv_NCount_PotExt, Recv_NDisp_PotExt, MPI_GAMER_REAL, MPI_COMM_WORLD );
#     endif // STORE_POT_GHOST
#     endif // GRAVITY

//    4.5 magnetic field (transfer one component at a time to avoid exceeding the maximum allowed transfer size in MPI)
#     ifdef MHD
      for (int v=0; v<NCOMP_MAG; v++)
      {
         MPI_Alltoallv_GAMER( SendBuf_Mag + v*SendDataSizeMag1v, Send_NCount_Mag1v, Send_NDisp_Mag1v, MPI_GAMER_REAL,
                              RecvBuf_Mag + v*RecvDataSizeMag1v, Recv_NCount_Mag1v, Recv_NDisp_Mag1v, MPI_GAMER_REAL, MPI_COMM_WORLD );
      }
#     endif
   } // if ( SendGridData )

#  ifdef PARTICLE
// 4.6 particle count
   MPI_Alltoallv( SendBuf_NPar, Send_NCount_Patch, Send_NDisp_Patch, MPI_INT,
                  RecvBuf_NPar, Recv_NCount_Patch, Recv_NDisp_Patch, MPI_INT, MPI_COMM_WORLD );

// 4.7 particle data
   MPI_Alltoallv_GAMER( SendBuf_ParData, Send_NCount_ParData, Send_NDisp_ParData, MPI_GAMER_REAL_PAR,
                        RecvBuf_ParData, Recv_NCount_ParData, Recv_NDisp_ParData, MPI_GAMER_REAL_PAR, MPI_COMM_WORLD );
#  endif // #ifdef PARTICLE


// 5. deallocate the MPI send buffers (BEFORE creating new patches to reduce the memory consumption)
// ==========================================================================================
   delete [] Send_NCount_Patch;
   delete [] Send_NDisp_Patch;
   delete [] Send_NCount_Flu1v;
   delete [] Send_NDisp_Flu1v;
   delete [] NDone_Patch;
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
#  ifdef MHD
   delete [] Send_NCount_Mag1v;
   delete [] Recv_NCount_Mag1v;
   delete [] Send_NDisp_Mag1v;
   delete [] Recv_NDisp_Mag1v;
   delete [] SendBuf_Mag;
#  endif
#  ifdef PARTICLE
   delete [] Send_NCount_ParData;
   delete [] Send_NDisp_ParData;
   delete [] NDone_ParData;
   delete [] SendBuf_ParData;
   delete [] SendBuf_NPar;
#  endif


// 6. allocate new patches with the data just received (use "patch group" as the basic unit)
//    --> also add particles to the particle repository and associate them with home patches
// ==========================================================================================
   const real *RecvPtr_Grid = NULL;
   const int   PScale       = PATCH_SIZE*amr->scale[lv];
   const int   PGScale      = 2*PScale;
   int PID, Cr0[3];

#  ifdef PARTICLE
// check: for RemoveParFromRepo == false, the size of particle repository should be exactly equal to the received particles
// --> see LB_RedistributeParticle_Init()
#  ifdef DEBUG_PARTICLE
   const long NParExpect = amr->Par->NPar_AcPlusInac + NRecv_Total_ParData/(long)PAR_NATT_TOTAL;
   if ( !RemoveParFromRepo  &&  NParExpect > amr->Par->ParListSize )
      Aux_Error( ERROR_INFO, "NParExpect (%ld) > ParListSize (%ld) !!\n", NParExpect, amr->Par->ParListSize );
#  endif

   const real_par *RecvPtr_Par = RecvBuf_ParData;
   long *ParList               = NULL;
   int   ParListSizeMax        = 0;    // must NOT be negative to deal with the case NRecv_Total_Patch == 0

   for (int t=0; t<NRecv_Total_Patch; t++)   ParListSizeMax = MAX( ParListSizeMax, RecvBuf_NPar[t] );

   ParList = new long [ParListSizeMax];
#  endif // #ifdef PARTICLE

   for (int PID0=0; PID0<NRecv_Total_Patch; PID0+=8)
   {
      LB_Idx = RecvBuf_LBIdx[PID0];

      LB_Index2Corner( lv, LB_Idx, Cr0, CHECK_ON );

      for (int d=0; d<3; d++)    Cr0[d] -= Cr0[d]%PGScale; // currently this line has no effect

//    6.1 allocate patches
//    father patch is still unkown ...
      amr->pnew( lv, Cr0[0],        Cr0[1],        Cr0[2],        -1, true, true, true );
      amr->pnew( lv, Cr0[0]+PScale, Cr0[1],        Cr0[2],        -1, true, true, true );
      amr->pnew( lv, Cr0[0],        Cr0[1]+PScale, Cr0[2],        -1, true, true, true );
      amr->pnew( lv, Cr0[0],        Cr0[1],        Cr0[2]+PScale, -1, true, true, true );
      amr->pnew( lv, Cr0[0]+PScale, Cr0[1]+PScale, Cr0[2],        -1, true, true, true );
      amr->pnew( lv, Cr0[0],        Cr0[1]+PScale, Cr0[2]+PScale, -1, true, true, true );
      amr->pnew( lv, Cr0[0]+PScale, Cr0[1],        Cr0[2]+PScale, -1, true, true, true );
      amr->pnew( lv, Cr0[0]+PScale, Cr0[1]+PScale, Cr0[2]+PScale, -1, true, true, true );

//    6.2 assign data
      for (int LocalID=0; LocalID<8; LocalID++)
      {
         PID = PID0 + LocalID;

         if ( SendGridData )
         {
//          fluid
            for (int v=0; v<NCOMP_TOTAL; v++)
            {
               RecvPtr_Grid = RecvBuf_Flu + v*RecvDataSizeFlu1v + PID*FluSize1v;
               memcpy( &amr->patch[FluSg][lv][PID]->fluid[v][0][0][0], RecvPtr_Grid, FluSize1v*sizeof(real) );
            }

#           ifdef GRAVITY
//          potential
            RecvPtr_Grid = RecvBuf_Pot + PID*FluSize1v;
            memcpy( &amr->patch[PotSg][lv][PID]->pot[0][0][0], RecvPtr_Grid, FluSize1v*sizeof(real) );

//          potential with ghost zones
#           ifdef STORE_POT_GHOST
            RecvPtr_Grid = RecvBuf_PotExt + PID*GraNxtSize;
            memcpy( &amr->patch[PotSg][lv][PID]->pot_ext[0][0][0], RecvPtr_Grid, GraNxtSize*sizeof(real) );
#           endif
#           endif // GRAVITY

//          magnetic field
#           ifdef MHD
            for (int v=0; v<NCOMP_MAG; v++)
            {
               RecvPtr_Grid = RecvBuf_Mag + v*RecvDataSizeMag1v + PID*MagSize1v;
               memcpy( &amr->patch[MagSg][lv][PID]->magnetic[v][0], RecvPtr_Grid, MagSize1v*sizeof(real) );
            }
#           endif
         } // if ( SendGridData )

//       particle
#        ifdef PARTICLE
         for (int p=0; p<RecvBuf_NPar[PID]; p++)
         {
//          add a single particle to the particle repository
            ParID        = amr->Par->AddOneParticle( RecvPtr_Par );
            RecvPtr_Par += PAR_NATT_TOTAL;

//          store the new particle index
            ParList[p] = ParID;

//          we do not transfer inactive particles
#           ifdef DEBUG_PARTICLE
            if ( amr->Par->Attribute[PAR_MASS][ParID] < (real)0.0 )
               Aux_Error( ERROR_INFO, "Transferring inactive particle (ParID %d, Mass %14.7e) !!\n",
                          ParID, amr->Par->Attribute[PAR_MASS][ParID] );
#           endif
         }

//       6.3 associate particles with their home patches
         const real_par *PType = amr->Par->Type;
#        ifdef DEBUG_PARTICLE
//       do not set ParPos too early since pointers to the particle repository (e.g., amr->Par->PosX)
//       may change after calling amr->Par->AddOneParticle()
         const real_par *ParPos[3] = { amr->Par->PosX, amr->Par->PosY, amr->Par->PosZ };
         char Comment[100];
         sprintf( Comment, "%s, PID %d, NPar %d", __FUNCTION__, PID, RecvBuf_NPar[PID] );
         amr->patch[0][lv][PID]->AddParticle( RecvBuf_NPar[PID], ParList, &amr->Par->NPar_Lv[lv],
                                              PType, ParPos, amr->Par->NPar_AcPlusInac, Comment );
#        else
         amr->patch[0][lv][PID]->AddParticle( RecvBuf_NPar[PID], ParList, &amr->Par->NPar_Lv[lv],
                                              PType );
#        endif
#        endif // #ifdef PARTICLE
      } // for (int LocalID=0; LocalID<8; LocalID++)
   } // for (int PID0=0; PID0<NRecv_Total_Patch; PID0+=8)

// 6.4 reset NPatchComma
   for (int m=1; m<28; m++)   amr->NPatchComma[lv][m] = NRecv_Total_Patch;

// check the amr->NPatchComma recording
   if ( amr->NPatchComma[lv][1] != amr->num[lv] )
      Aux_Error( ERROR_INFO, "amr->NPatchComma[%d][1] (%d) != amr->num[%d] (%d) !!\n",
                 lv, amr->NPatchComma[lv][1], lv, amr->num[lv] );


// 7. record LB_IdxList_Real
// ==========================================================================================
   delete [] amr->LB->IdxList_Real         [lv];
   delete [] amr->LB->IdxList_Real_IdxTable[lv];

   amr->LB->IdxList_Real         [lv] = new long [NRecv_Total_Patch];
   amr->LB->IdxList_Real_IdxTable[lv] = new int  [NRecv_Total_Patch];

   for (int PID=0; PID<NRecv_Total_Patch; PID++)   amr->LB->IdxList_Real[lv][PID] = amr->patch[0][lv][PID]->LB_Idx;

   Mis_Heapsort( NRecv_Total_Patch, amr->LB->IdxList_Real[lv], amr->LB->IdxList_Real_IdxTable[lv] );


// 8. deallocate the MPI recv buffers
// ==========================================================================================
   delete [] Recv_NCount_Patch;
   delete [] Recv_NDisp_Patch;
   delete [] Recv_NCount_Flu1v;
   delete [] Recv_NDisp_Flu1v;
   delete [] RecvBuf_LBIdx;
   delete [] RecvBuf_Flu;
#  ifdef GRAVITY
   delete [] RecvBuf_Pot;
#  ifdef STORE_POT_GHOST
   delete [] RecvBuf_PotExt;
#  endif
#  endif // GRAVITY
#  ifdef MHD
   delete [] RecvBuf_Mag;
#  endif
#  ifdef PARTICLE
   delete [] Recv_NCount_ParData;
   delete [] Recv_NDisp_ParData;
   delete [] RecvBuf_ParData;
   delete [] RecvBuf_NPar;
   delete [] ParList;
#  endif

} // FUNCTION : LB_RedistributeRealPatch



//-------------------------------------------------------------------------------------------------------
// Function    :  LB_SortRealPatch
// Description :  Sort real patches on the target level according to their load-balance indices
//
// Note        :  1. Must be called before allocating any buffer patches
//                2. It will reset amr->LB->IdxList_Real[lv] and amr->LB->IdxList_Real_IdxTable[lv]
//
// Parameter   :  lv : Target refinement level
//-------------------------------------------------------------------------------------------------------
void LB_SortRealPatch( const int lv )
{

   const int NReal = amr->NPatchComma[lv][1];

// check: there must be no buffer patches
   if ( NReal != amr->num[lv] )
      Aux_Error( ERROR_INFO, "NReal (%d) != amr->num[%d] (%d) !!\n", NReal, lv, amr->num[lv] );


// 1. back up patch pointers
   patch_t *(*patch_ptr)[2] = new patch_t* [NReal][2];

   for (int PID=0; PID<NReal; PID++)
   for (int Sg=0; Sg<2; Sg++)
      patch_ptr[PID][Sg] = amr->patch[Sg][lv][PID];


// 2. sort real patches by LB_Idx
//    --> use patch groups instead of patches as the sorting unit since we must fix the order of patches
//        within the same patch group
//    --> HILBERT guarantees patches within the same patch group have consecutive LB_Idx
#  if ( defined LOAD_BALANCE  &&  LOAD_BALANCE != HILBERT )
#  warning : validate the adopted LOAD_BALANCE scheme supports the following code !!
#  endif
   for (int NewPID0=0; NewPID0<NReal; NewPID0+=8)
   {
      int OldPID0 = amr->LB->IdxList_Real_IdxTable[lv][NewPID0];
      OldPID0 = OldPID0 - OldPID0%8;

      for (int LocalID=0; LocalID<8; LocalID++)
      for (int Sg=0; Sg<2; Sg++)
         amr->patch[Sg][lv][NewPID0+LocalID] = patch_ptr[OldPID0+LocalID][Sg];
   }

   delete [] patch_ptr;


// 3. reset amr->LB->IdxList_Real[lv] and amr->LB->IdxList_Real_IdxTable[lv]
// --> better still reallocate memory since LB_RedistributeRealPatch() is not called when Redistribute==false
   delete [] amr->LB->IdxList_Real         [lv];
   delete [] amr->LB->IdxList_Real_IdxTable[lv];

   amr->LB->IdxList_Real         [lv] = new long [NReal];
   amr->LB->IdxList_Real_IdxTable[lv] = new int  [NReal];

   for (int PID=0; PID<NReal; PID++)   amr->LB->IdxList_Real[lv][PID] = amr->patch[0][lv][PID]->LB_Idx;

   Mis_Heapsort( NReal, amr->LB->IdxList_Real[lv], amr->LB->IdxList_Real_IdxTable[lv] );


// 4. check
#  ifdef GAMER_DEBUG
// 4-1. order of patches with LocalID==0
   for (int PID0=8; PID0<NReal; PID0+=8)
   {
      const long LBIdx1 = amr->patch[0][lv][PID0  ]->LB_Idx;
      const long LBIdx2 = amr->patch[0][lv][PID0-8]->LB_Idx;

      if ( LBIdx1 <= LBIdx2 )    Aux_Error( ERROR_INFO, "LB_Idx %ld <= %ld (lv %d, PID1 %d, PID2 %d) !!\n",
                                            LBIdx1, LBIdx2, lv, PID0, PID0-8 );
   }

// 4-2. order of patches within each patch group
   const int PScale = PS1*amr->scale[lv];

   for (int PID0=0; PID0<NReal; PID0+=8)
   for (int LocalID=1; LocalID<8; LocalID++)
   {
      const int PID1 = PID0 + LocalID;

      for (int d=0; d<3; d++)
      {
         const int corner0 = amr->patch[0][lv][PID0]->corner[d];
         const int corner1 = amr->patch[0][lv][PID1]->corner[d];
         const int offset = TABLE_02( LocalID, 'x'+d, 0, PScale );

         if ( corner1 != corner0 + offset )
         {
            Output_Patch( lv, PID0, 0, 0, 0, "WrongOrder" );
            Output_Patch( lv, PID1, 0, 0, 0, "WrongOrder" );
            Aux_Error( ERROR_INFO, "incorrect order of local patches (lv %d, PID0 %d, PID1 %d, d %d, corner0 %d, corner1 %d, offset %d)\n"
                       "        --> Check the files Patch_*_WrongOrder !!\n", lv, PID0, PID1, d, corner0, corner1, offset );
         }
      }
   } // for (int LocalID=1; LocalID<8; LocalID++)
#  endif // #ifdef GAMER_DEBUG

} // FUNCTION : LB_SortRealPatch



#ifdef PARTICLE
//-------------------------------------------------------------------------------------------------------
// Function    :  LB_RedistributeParticle_Init
// Description :  Initialize the procedure for redistributing particles
//
// Note        :  1. This function will get the total number of particles AFTER data redistribution and
//                   then allocate the new particle attribute arrays
//                2. This function will also reallocate particle repository by calling amr->Par->InitRepo().
//                   However, we reset NPar_AcPlusInac to zero since we will update it level by level
//                   when calling LB_RedistributeRealPatch() later.
//                3. One must call LB_SetCutPoint() for all levels in advance
//
// Parameter   :  ParAtt_Old : Pointers for backing up the old particle attribute arrays (amr->Par->Attribute[])
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void LB_RedistributeParticle_Init( real_par **ParAtt_Old )
{

// backup the old particle attribute arrays
// remember to reset Attribute[] to NULL so that amr->Par->InitRepo will NOT delete these arrays
   for (int v=0; v<PAR_NATT_TOTAL; v++)
   {
      ParAtt_Old         [v] = amr->Par->Attribute[v];
      amr->Par->Attribute[v] = NULL;
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


// reset particle variables (do not reset NPar_Lv since we will need it for debug in LB_RedistributeRealPatch())
   amr->Par->InitRepo( Recv_NPar_Sum, MPI_NRank );

// reset the total number of particles to be zero
// --> so particle repository is pre-allocated, but it contains no active particle yet
// --> we will add active particles in LB_RedistributeRealPatch()
   amr->Par->NPar_AcPlusInac = 0;
   amr->Par->NPar_Active     = 0;

} // FUNCTION : LB_RedistributeParticle_Init



//-------------------------------------------------------------------------------------------------------
// Function    :  LB_RedistributeParticle_End
// Description :  End the procedure for redistributing particles
//
// Note        :  1. Free old particle attribute arrays
//
// Parameter   :  ParAtt_Old : Pointers for backing up the old particle attribute arrays (amr->Par->Attribute[])
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void LB_RedistributeParticle_End( real_par **ParAtt_Old )
{

// remove old particle attribute arrays
   for (int v=0; v<PAR_NATT_TOTAL; v++)   free( ParAtt_Old [v] );


// check the total number of particles
   if ( amr->Par->NPar_AcPlusInac != amr->Par->NPar_Active )
      Aux_Error( ERROR_INFO, "NPar_AcPlusInac (%ld) != NPar_Active (%ld) !!\n",
                 amr->Par->NPar_AcPlusInac, amr->Par->NPar_Active );

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
