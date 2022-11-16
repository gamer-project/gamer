#include "GAMER.h"



/*
void GID2PID(const long gid, int &level, int &PID) {

#ifdef GAMER_DEBUG
    long NPatchAllLv = 0;
    for (int lv=0; lv<NLEVEL; lv++)  NPatchAllLv += NPatchTotal[lv];
    if ( gid < 0  ||  gid >= NPatchAllLv )  Aux_Error( ERROR_INFO, "incorrect gid %ld (max = %ld) !!\n", gid, NPatchAllLv-1 );
#endif

    level = 0;
    for(int lv=1; lv<NLEVEL; lv++){
        if( gid <  GID_Offset[lv] )      break;
        level = lv;
    }
    PID = gid - GID_Offset[level];
}*/

void GID2PID(const long gid, int &level, int &PID){

#ifdef GAMER_DEBUG
    long NPatchAllLv = 0;
    for (int lv=0; lv<NLEVEL; lv++)  NPatchAllLv += NPatchTotal[lv];
    if ( gid < 0  ||  gid >= NPatchAllLv )  Aux_Error( ERROR_INFO, "incorrect gid %ld (max = %ld) !!\n", gid, NPatchAllLv-1 );
#endif

    level = 0;
    for(int lv=1; lv<NLEVEL; lv++){
        if( gid <  YT_GID_Offset[lv] )      break;
        level = lv;
    }
    PID = gid - YT_GID_Offset[*level];
}

// should also do automatic conversion for negative PIDs with the help of the information gathered from all ranks 
void PID2GID(const int PID, const int level, long &gid) {

#ifdef GAMER_DEBUG
    if ( gid < 0  ||  PID >= NPatchAll[level])  Aux_Error( ERROR_INFO, "incorrect gid %ld (max = %ld) !!\n", gid, NPatchAllLv-1 );
#endif

   if ( PID < 0 ) {
   }

   // patch present on AMR rank
   else if ( PID < (amr->NPatchComma[level][1]) ) {
      gid = (long) (PID + YT_GID_Offset[level]);

   // father patch that is absent
   } else {
      long LBIdx = amr->patch[0][level][PID]->LB_Idx;
      int MatchIdx;
      Mis_Matching_int( NPatchTotal[level], LBIdxList_Sort[level], 1, &FaLBIdx, &MatchIdx );

#     ifdef DEBUG_HDF5
      if ( MatchIdx < 0 ){
          Aux_Error( ERROR_INFO, "Lv %d, PID %d, LBIdx %ld, couldn't find a matching patch !!\n",
                     level, PID, FaLBIdx );
      }
#     endif

      gid = LBIdxList_Sort_IdxTable[level][MatchIdx] + GID_LvStart[level];
   }
}


//-------------------------------------------------------------------------------------------------------
// Function    :  LB_SyncTree
// Description :  Sync AMR structure excluding the field data to all MPI ranks
//
// Note        :  
//
// Parameter   :  
//-------------------------------------------------------------------------------------------------------
void LB_SyncTree( )
{

//***********************************************************************************
// THE AMR structure 
//***********************************************************************************

/*
// data members
// ===================================================================================
   patch_t    *patch[2][NLEVEL][MAX_PATCH];

#  ifdef PARTICLE
   Particle_t *Par;
#  endif

#  ifndef SERIAL
   ParaVar_t  *ParaVar;
#  ifdef LOAD_BALANCE
   LB_t       *LB;
#  endif
#  endif

   int    num         [NLEVEL];
   int    scale       [NLEVEL];
   int    FluSg       [NLEVEL];
   double FluSgTime   [NLEVEL][2];
   int    MagSg       [NLEVEL];     // for convenience, it is defined even when MHD is disabled
#  ifdef MHD
   double MagSgTime   [NLEVEL][2];
#  endif
#  ifdef GRAVITY
   int    PotSg       [NLEVEL];
   double PotSgTime   [NLEVEL][2];
#  endif
   int    NPatchComma [NLEVEL][28];
   double dh          [NLEVEL];
   int    ResPower2   [NLEVEL];
   double BoxEdgeL    [3];
   double BoxEdgeR    [3];
   double BoxCenter   [3];
   double BoxSize     [3];
   int    BoxScale    [3];
   bool   WithFlux;
#  ifdef MHD
   bool   WithElectric;
#  endif
   long   NUpdateLv   [NLEVEL];

#if ( MODEL == ELBDM  && ELBDM_SCHEME == HYBRID )
   bool   use_wave_flag[NLEVEL];
#endif // #if ( MODEL == ELBDM  && ELBDM_SCHEME == HYBRID )
*/

//***********************************************************************************
// The patch structure
//***********************************************************************************

/*
// data members
// ===================================================================================
   real (*fluid)[PS1][PS1][PS1];


#  ifdef MHD
   real (*magnetic)[ PS1P1*SQR(PS1) ];
#  endif

#  ifdef GRAVITY
   real (*pot)[PS1][PS1];
#  ifdef STORE_POT_GHOST
   real (*pot_ext)[GRA_NXT][GRA_NXT];
#  endif
#  endif // GRAVITY

#  ifdef DUAL_ENERGY
   char (*de_status)[PS1][PS1];
#  endif

#  ifdef PARTICLE
   real (*rho_ext)[RHOEXT_NXT][RHOEXT_NXT];
#  endif

   real (*flux       [6])[PS1][PS1];
   real (*flux_tmp   [6])[PS1][PS1];
#  ifdef BIT_REP_FLUX
   real (*flux_bitrep[6])[PS1][PS1];
#  endif

#  ifdef MHD
   real (*electric       [18]);
   real (*electric_tmp   [18]);
#  ifdef BIT_REP_ELECTRIC
   real (*electric_bitrep[18]);
#  endif
   bool ele_corrected[12];
#  endif

   int    corner[3];
   int    sibling[26];
   int    father;
   int    son;
   bool   flag;
   
   bool   Active;
   double EdgeL[3];
   double EdgeR[3];

   ulong  PaddedCr1D;
   long   LB_Idx;


#  ifdef PARTICLE
   int    NPar;
   int    ParListSize;
   long  *ParList;

   int    NPar_Copy;
#  ifdef LOAD_BALANCE
   real  *ParAtt_Copy[PAR_NATT_TOTAL];
#  else
   long  *ParList_Copy;
#  endif

   int    NPar_Escp[26];
   long  *ParList_Escp[26];

#  endif

#  if ( MODEL == ELBDM  && ELBDM_SCHEME == HYBRID )
   bool use_wave_flag;
#  endif // #  if ( MODEL == ELBDM  && ELBDM_SCHEME == HYBRID )
*/

// Setting up GIDs


// 1. gather the number of patches at different MPI ranks, calculate number of local patches
//    and set the corresponding GID offset
   int (*NPatchAllRank)[NLEVEL] = new int [MPI_NRank][NLEVEL];
   int NPatchLocal[NLEVEL], NPatchAllLv=0, NPatchLocalLv=0, GID_LvStart[NLEVEL];

   for (int lv=0; lv<NLEVEL; lv++)
   {
      NPatchLocal[lv] = amr->NPatchComma[lv][1];
      NPatchLocalLv = NPatchLocalLv + NPatchLocal[lv];
   }

   MPI_Allgather( NPatchLocal, NLEVEL, MPI_INT, NPatchAllRank[0], NLEVEL, MPI_INT, MPI_COMM_WORLD );

   for (int lv=0; lv<NLEVEL; lv++)
   {
      // set YT_GID_Offset for searching GID in derived function and particle get attribute function.
      YT_GID_Offset[lv] = 0;

      for (int r=0; r<MPI_Rank; r++)      YT_GID_Offset[lv] += NPatchAllRank[r][lv];

      for (int FaLv=0; FaLv<lv; FaLv++)   YT_GID_Offset[lv] += NPatchTotal[FaLv];

      NPatchAllLv += NPatchTotal[lv];

      GID_LvStart[lv] = ( lv == 0 ) ? 0 : GID_LvStart[lv-1] + NPatchTotal[lv-1];
   }



// get the search table for converting from negative PID to GID + knowledge which rank, needed by parent_id
   int   RecvCount_LBIdx[MPI_NRank], RecvDisp_LBIdx[MPI_NRank];
   long *LBIdxList_Sort[NLEVEL], *LBIdxList_Local[NLEVEL];
   int  *LBIdxList_Sort_IdxTable[NLEVEL];

   for (int lv=0; lv<NLEVEL; lv++)
   {
      LBIdxList_Local        [lv] = new long [ amr->NPatchComma[lv][1] ];
      LBIdxList_Sort         [lv] = new long [ NPatchTotal[lv] ];
      LBIdxList_Sort_IdxTable[lv] = new int  [ NPatchTotal[lv] ];
   }

   for (int lv=0; lv<NLEVEL; lv++)
   {
      for (int r=0; r<MPI_NRank; r++)
      {
         RecvCount_LBIdx[r] = NPatchAllRank[r][lv];
         RecvDisp_LBIdx [r] = ( r == 0 ) ? 0 : RecvDisp_LBIdx[r-1] + RecvCount_LBIdx[r-1];
      }

      for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
         LBIdxList_Local[lv][PID] = amr->patch[0][lv][PID]->LB_Idx;

      MPI_Allgatherv( LBIdxList_Local[lv], amr->NPatchComma[lv][1], MPI_LONG,
                      LBIdxList_Sort[lv], RecvCount_LBIdx, RecvDisp_LBIdx, MPI_LONG,
                      MPI_COMM_WORLD );
   }

   for (int lv=0; lv<NLEVEL; lv++)
      Mis_Heapsort( NPatchTotal[lv], LBIdxList_Sort[lv], LBIdxList_Sort_IdxTable[lv] );



 // Set up patches for all levels
   struct GlobalPatch {
      long GID; 
      long son;
      long father;
      int  lv;
      int  rank; 
      int  son_rank; 
      int  father_rank; 

      GlobalPatch() : GID(-1), son(-1), father(-1), lv(-1), rank(-1), son_rank(-1), father_rank(-1) {
      }
   };

   GlobalPatch  *global_patch = new MyPatch[NPatchAllLv];

// Iterate over real patches
   for (int lv=0; lv<NLEVEL; lv++)
   {
      for (int PID = 0; PID < amr->NPatchComma[lv][1] ; PID++)
      {
         const long GID         = PID + YT_GID_Offset[lv];
         global_patch[GID].gid  = gid;
         global_patch[GID].lv   = lv;
         global_patch[GID].rank = MPI_Rank; 
         

//       Get son's GID and rank
         {
            const int SonPID = amr->patch[0][lv][PID]->son; 

            const int FaLv  = lv; 
            const int SonLv = lv + 1; 

//          No son
            if ( SonPID == -1 ) {
               global_patch[GID].son = -1;
            } 
//          Son at home
            else if ( SonPID > -1 ) {
               const long SonGID           = SonPID + YT_GID_Offset[SonLv];
               global_patch[GID].son       = SonGID;
               global_patch[GID].son_rank  = MPI_Rank;
            } 
//          Son abroad (only possible in LOAD_BALANCE)
            else {
//             Determine the target rank and LB_Idx
#              if ( LOAD_BALANCE == HILBERT )
               const long SonLBIdx = 8*amr->patch[0][FaLv][PID]->LB_Idx;  // faster
#              else
               const long SonLBIdx = LB_Corner2Index( SonLv, amr->patch[0][FaLv][PID]->corner, CHECK_OFF );
#              endif

               int MatchIdx;
               Mis_Matching_int( NPatchTotal[SonLv], LBIdxList_Sort[SonLv], 1, &SonLBIdx, &MatchIdx );

#              ifdef DEBUG_HDF5
               if ( MatchIdx < 0 ){
                   Aux_Error( ERROR_INFO, "Lv %d, PID %d, SonLBIdx %ld, couldn't find a matching patch !!\n",
                              SonLv, SonPID, SonLBIdx );
               }
#              endif

               const long SonGID = LBIdxList_Sort_IdxTable[SonLv][MatchIdx] + GID_LvStart[SonLv];

               global_patch[GID].son       = SonGID;
               global_patch[GID].son_rank  = SON_OFFSET_LB - SonPID;

            }
         }

//       Get father's GID and rank
         {
            const int FaPID = amr->patch[0][lv][PID]->father;
            const int FaLv  = lv - 1;
            const int SonLv = lv; 

//          No father (only possible for the root patches)
            if ( FaPID < 0 ) {
#              ifdef DEBUG_HDF5
               if ( lv != 0 )       Aux_Error( ERROR_INFO, "Lv %d, PID %d, FaPID %d < 0 !!\n", lv, PID, FaPID );
               if ( FaPID != -1 )   Aux_Error( ERROR_INFO, "Lv %d, PID %d, FaPID %d < 0 but != -1 !!\n", lv, PID, FaPID );
#              endif
               global_patch[GID].father = -1; 
            }
//          Father at home 
            else if ( FaPID < amr->NPatchComma[FaLv][1] ){
               const long FaGID = (long) (FaPID + YT_GID_Offset[FaLv]);
               global_patch[GID].father        = FaGID;
               global_patch[GID].father_rank   = MPI_Rank;
            }
//          Father abroad (only possible in LOAD_BALANCE)
            else {
#              ifdef DEBUG_HDF5
#              ifndef LOAD_BALANCE
               Aux_Error( ERROR_INFO, "Lv %d, PID %d, FaPID %d >= NRealFaPatch %d (only possible in LOAD_BALANCE) !!\n",
                          lv, PID, FaPID, amr->NPatchComma[FaLv][1] );
#              endif
               if ( FaPID >= (amr->num[FaLv]) ){
                   Aux_Error( ERROR_INFO, "Lv %d, PID %d, FaPID %d >= total number of patches %d !!\n",
                              lv, PID, FaPID, amr->num[FaLv] );
               }
#              endif // DEBUG_HDF5

               long FaLBIdx = amr->patch[0][FaLv][FaPID]->LB_Idx;
               int MatchIdx;
               Mis_Matching_int( NPatchTotal[FaLv], LBIdxList_Sort[FaLv], 1, &FaLBIdx, &MatchIdx );

#              ifdef DEBUG_HDF5
               if ( MatchIdx < 0 ){
                   Aux_Error( ERROR_INFO, "Lv %d, PID %d, FaPID %d, FaLBIdx %ld, couldn't find a matching patch !!\n",
                              lv, PID, FaPID, FaLBIdx );
               }
#              endif
               const long FaGID = (long) (LBIdxList_Sort_IdxTable[FaLv][MatchIdx] + GID_LvStart[FaLv]);
               global_patch[GID].father        = FaGID;
               global_patch[GID].father_rank   = ;

            }
         }

      }
   }



// Synchronise global patches between different MPI ranks

// -> allgather and then appropriate sorting 


   if (MPI_NRank == 0) {
   rand_nums = create_rand_nums(elements_per_proc * world_size);
   }

   // Create a buffer that will hold a subset of the random numbers
   float *sub_rand_nums = malloc(sizeof(float) * elements_per_proc);

   // Scatter the random numbers to all processes
   MPI_Scatter(rand_nums, elements_per_proc, MPI_FLOAT, sub_rand_nums,
               elements_per_proc, MPI_FLOAT, 0, MPI_COMM_WORLD);

   // Compute the average of your subset
   float sub_avg = compute_avg(sub_rand_nums, elements_per_proc);
   // Gather all partial averages down to the root process
   AMR_t *trees = NULL;
   if (world_rank == 0) {
      trees = malloc(sizeof(AMR_t) * world_size);
   }
   MPI_Gather(&sub_avg, 1, MPI_FLOAT, trees, 1, MPI_FLOAT, 0,
            MPI_COMM_WORLD);

   // Compute the total average of all numbers.
   if (MPI_Rank == 0) {
   float avg = compute_avg(sub_avgs, world_size);
   }

} // FUNCTION : Aux_SwapPointer


// get the search table, needed by parent_id
   int   RecvCount_LBIdx[MPI_NRank], RecvDisp_LBIdx[MPI_NRank];
   long *LBIdxList_Sort[NLEVEL], *LBIdxList_Local[NLEVEL];
   int  *LBIdxList_Sort_IdxTable[NLEVEL];

   for (int lv=0; lv<NLEVEL; lv++)
   {
      LBIdxList_Local        [lv] = new long [ amr->NPatchComma[lv][1] ];
      LBIdxList_Sort         [lv] = new long [ NPatchTotal[lv] ];
      LBIdxList_Sort_IdxTable[lv] = new int  [ NPatchTotal[lv] ];
   }

   for (int lv=0; lv<NLEVEL; lv++)
   {
      for (int r=0; r<MPI_NRank; r++)
      {
         RecvCount_LBIdx[r] = NPatchAllRank[r][lv];
         RecvDisp_LBIdx [r] = ( r == 0 ) ? 0 : RecvDisp_LBIdx[r-1] + RecvCount_LBIdx[r-1];
      }

      for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
         LBIdxList_Local[lv][PID] = amr->patch[0][lv][PID]->LB_Idx;

      MPI_Allgatherv( LBIdxList_Local[lv], amr->NPatchComma[lv][1], MPI_LONG,
                      LBIdxList_Sort[lv], RecvCount_LBIdx, RecvDisp_LBIdx, MPI_LONG,
                      MPI_COMM_WORLD );
   }

   for (int lv=0; lv<NLEVEL; lv++)
      Mis_Heapsort( NPatchTotal[lv], LBIdxList_Sort[lv], LBIdxList_Sort_IdxTable[lv] );

// loop over local patches at all levels
   for (int lv=0; lv<NLEVEL; lv++)
   {
      const int FaLv  = lv - 1;
      const int FluSg = amr->FluSg[lv];
#     ifdef GRAVITY
      const int PotSg = amr->PotSg[lv];
#     endif

#     ifdef MHD
      const int MagSg = amr->MagSg[lv];
#     endif

#     ifdef LIBYT_USE_PATCH_GROUP
      for (int PID=0; PID<(amr->NPatchComma[lv][1]); PID+=8)
#     else
      for (int PID=0; PID<(amr->NPatchComma[lv][1]); PID++)
#     endif // #ifdef LIBYT_USE_PATCH_GROUP
      {
         const int GID = PID + YT_GID_Offset[lv];

         for (int d=0; d<3; d++)
         {
#           ifdef LIBYT_USE_PATCH_GROUP
            YT_Grids[LID].left_edge [d] = amr->patch[0][lv][PID    ]->EdgeL[d];
            YT_Grids[LID].right_edge[d] = amr->patch[0][lv][PID + 7]->EdgeR[d];
            YT_Grids[LID].grid_dimensions[d] = PATCH_SIZE * 2;
#           else
            YT_Grids[LID].left_edge [d] = amr->patch[0][lv][PID]->EdgeL[d];
            YT_Grids[LID].right_edge[d] = amr->patch[0][lv][PID]->EdgeR[d];
            YT_Grids[LID].grid_dimensions[d] = PATCH_SIZE;
#           endif // #ifdef LIBYT_USE_PATCH_GROUP
         }

#        ifdef PARTICLE
#        ifdef LIBYT_USE_PATCH_GROUP
         // input particle num in this patch group.
         long particle_count = 0;
         for(int i=PID; i<PID+8; i++){
             particle_count += (long) amr->patch[0][lv][i]->NPar;
         }
         YT_Grids[LID].particle_count_list[0] = particle_count;
#        else
         // input particle num in this grid
         YT_Grids[LID].particle_count_list[0] = (long) amr->patch[0][lv][PID]->NPar;
#        endif // #ifdef LIBYT_USE_PATCH_GROUP
#        endif // #ifdef PARTICLE

#        ifdef LIBYT_USE_PATCH_GROUP
         YT_Grids[LID].id     = (long) GID / 8;
#        else
         YT_Grids[LID].id     = (long) GID;
#        endif // #ifdef LIBYT_USE_PATCH_GROUP

         YT_Grids[LID].level  = lv;

         // getting parent's id
         int FaPID = amr->patch[0][lv][PID]->father;
         int FaLv = lv - 1;

         if ( FaPID < 0 ){
            // no father (only possible for the root patches)
#           ifdef DEBUG_HDF5
            if ( lv != 0 )       Aux_Error( ERROR_INFO, "Lv %d, PID %d, FaPID %d < 0 !!\n", lv, PID, FaPID );
            if ( FaPID != -1 )   Aux_Error( ERROR_INFO, "Lv %d, PID %d, FaPID %d < 0 but != -1 !!\n", lv, PID, FaPID );
#           endif
            YT_Grids[LID].parent_id = -1;
         }
         else if ( FaPID < (amr->NPatchComma[FaLv][1]) ){
            // has father patch
#           ifdef LIBYT_USE_PATCH_GROUP
            YT_Grids[LID].parent_id = (long) (FaPID + YT_GID_Offset[FaLv]) / 8;
#           else
            YT_Grids[LID].parent_id = (long) (FaPID + YT_GID_Offset[FaLv]);
#           endif // #ifdef LIBYT_USE_PATCH_GROUP
         }
         else{
            // father patch is a buffer patch (only possible in LOAD_BALANCE)
#           ifdef DEBUG_HDF5
#           ifndef LOAD_BALANCE
            Aux_Error( ERROR_INFO, "Lv %d, PID %d, FaPID %d >= NRealFaPatch %d (only possible in LOAD_BALANCE) !!\n",
                       lv, PID, FaPID, amr->NPatchComma[FaLv][1] );
#           endif
            if ( FaPID >= (amr->num[FaLv]) ){
                Aux_Error( ERROR_INFO, "Lv %d, PID %d, FaPID %d >= total number of patches %d !!\n",
                           lv, PID, FaPID, amr->num[FaLv] );
            }
#           endif // DEBUG_HDF5

            long FaLBIdx = amr->patch[0][FaLv][FaPID]->LB_Idx;
            int MatchIdx;
            Mis_Matching_int( NPatchTotal[FaLv], LBIdxList_Sort[FaLv], 1, &FaLBIdx, &MatchIdx );

#           ifdef DEBUG_HDF5
            if ( MatchIdx < 0 ){
                Aux_Error( ERROR_INFO, "Lv %d, PID %d, FaPID %d, FaLBIdx %ld, couldn't find a matching patch !!\n",
                           lv, PID, FaPID, FaLBIdx );
            }
#           endif

#           ifdef LIBYT_USE_PATCH_GROUP
            YT_Grids[LID].parent_id = (long) (LBIdxList_Sort_IdxTable[FaLv][MatchIdx] + GID_LvStart[FaLv]) / 8;
#           else
            YT_Grids[LID].parent_id = (long) (LBIdxList_Sort_IdxTable[FaLv][MatchIdx] + GID_LvStart[FaLv]);
#           endif // #ifdef LIBYT_USE_PATCH_GROUP
         }

#        ifndef LIBYT_USE_PATCH_GROUP
         // load patch data to libyt if not use LIBYT_USE_PATCH_GROUP
         for (int v = 0; v < NCOMP_TOTAL; v++){
             YT_Grids[LID].field_data[v].data_ptr = amr->patch[FluSg][lv][PID]->fluid[v];
         }

#        ifdef GRAVITY
         // find field index of GRAVITY
         int PotIdx = 0;
         for ( int v = 0; v < NField; v++ ){
             if ( strcmp(FieldList[v].field_name, PotLabel) == 0 ){
                 PotIdx = v;
                 break;
             }
         }
         // load Pote patch data to libyt
         YT_Grids[LID].field_data[PotIdx].data_ptr = amr->patch[PotSg][lv][PID]->pot;
#        endif // #ifdef GRAVITY

#        ifdef MHD
         // find field index of CCMagX
         int MHDIdx = 0;
         for ( int v = 0; v < NField; v++ ){
             if ( strcmp(FieldList[v].field_name, "CCMagX") == 0 ){
                 MHDIdx = v;
                 break;
             }
         }

         for (int v = 0; v < NCOMP_MAG; v++){
             // input the data pointer
             YT_Grids[LID].field_data[ MHDIdx + v ].data_ptr = amr->patch[MagSg][lv][PID]->magnetic[v];

             // input the field dimension, since MHD has different dimension.
             for (int d = 0; d < 3; d++){
                 YT_Grids[LID].field_data[ MHDIdx + v ].data_dimensions[d] = PATCH_SIZE;

                 if ( strcmp(FieldList[ MHDIdx + v ].field_name, "CCMagX") == 0 && d == 2) {
                     YT_Grids[LID].field_data[ MHDIdx + v ].data_dimensions[d] = PATCH_SIZE + 1;
                 }
                 else if ( strcmp(FieldList[ MHDIdx + v ].field_name, "CCMagY") == 0 && d == 1) {
                     YT_Grids[LID].field_data[ MHDIdx + v ].data_dimensions[d] = PATCH_SIZE + 1;
                 }
                 else if ( strcmp(FieldList[ MHDIdx + v ].field_name, "CCMagZ") == 0 && d == 0) {
                     YT_Grids[LID].field_data[ MHDIdx + v ].data_dimensions[d] = PATCH_SIZE + 1;
                 }
             }
         }
#        endif // #ifdef MHD

#        endif // #ifndef LIBYT_USE_PATCH_GROUP
         
         LID = LID + 1;
      }
   }

   // free resources
   for (int lv=0; lv<NLEVEL; lv++)
   {
      delete [] LBIdxList_Local        [lv];
      delete [] LBIdxList_Sort         [lv];
      delete [] LBIdxList_Sort_IdxTable[lv];
   }

   if ( yt_commit_grids( ) != YT_SUCCESS )  Aux_Error( ERROR_INFO, "yt_commit_grids() failed !!\n" );

   if ( OPT__VERBOSE  &&  MPI_Rank == 0 )   Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : YT_AddLocalGrid

#endif // #ifdef SUPPORT_LIBYT


#include "GAMER.h"

#ifdef SUPPORT_LIBYT

//-------------------------------------------------------------------------------------------------------
// Function    :  YT_GetPID
// Description :  Get PID from YT passed in gid.
//
// Note        :  1. Get YT passed in gid, and then write its corresponding level and PID at *level and *PID.
//                2. Whether or not using LIBYT_USE_PATCH_GROUP would work.
//
// Parameter   :  gid      : gid YT passed in.
//                level    : fill in level of this gid.
//                PID      : fill in PID of this gid.
//
// Return      :  *level, *PID
//-------------------------------------------------------------------------------------------------------
void YT_GetPID(const long gid, int *level, int *PID){
#ifdef GAMER_DEBUG
    long NPatchAllLv = 0;
    for (int lv=0; lv<NLEVEL; lv++)  NPatchAllLv += NPatchTotal[lv];
    if ( gid < 0  ||  gid >= NPatchAllLv )  Aux_Error( ERROR_INFO, "incorrect gid %ld (max = %ld) !!\n", gid, NPatchAllLv-1 );
#endif
    *level = 0;
    for(int lv=1; lv<NLEVEL; lv++){
#ifdef  LIBYT_USE_PATCH_GROUP
        if( gid < (YT_GID_Offset[lv] / 8) ) break;
#else
        if( gid <  YT_GID_Offset[lv] )      break;
#endif // #ifdef  LIBYT_USE_PATCH_GROUP
        *level = lv;
    }
#ifdef LIBYT_USE_PATCH_GROUP
    *PID = gid * 8 - YT_GID_Offset[*level];
#else
    *PID = gid - YT_GID_Offset[*level];
#endif // #ifdef  LIBYT_USE_PATCH_GROUP
}

#endif // #ifdef SUPPORT_LIBYT


#include "GAMER.h"

#ifdef SUPPORT_LIBYT




//-------------------------------------------------------------------------------------------------------
// Function    :  YT_SetParameter
// Description :  Set YT-specific parameters for the inline analysis
//
// Note        :  1. This function must be called in advance **every time** we invoke the inline analysis
//                2. Invoked by YT_Inline().
//                3. Set up num_species, species_list for supporting PARTICLE.
//
// Parameter   :  NPatchAllLv : Total number of patches at all levels
//                NField      : Total number of fields
//                NPatchLocal : Number of local patches at all levels
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void YT_SetParameter( const int NPatchAllLv, const int NField, const int NPatchLocalLv )
{

   if ( OPT__VERBOSE  &&  MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );


// 1. prepare the simulation information for libyt
   yt_param_yt param_yt;

   param_yt.frontend                = "gamer";           // simulation frontend
   if ( strcmp(YT_FIG_BASENAME, "") != 0 )
       param_yt.fig_basename = YT_FIG_BASENAME;          // figure base name, use default if not set (default=Fig%09d)

   param_yt.length_unit             = UNIT_L;            // units are in cgs
   param_yt.mass_unit               = UNIT_M;
   param_yt.time_unit               = UNIT_T;

#  ifdef MHD
   param_yt.magnetic_unit           = UNIT_B;
#  endif

   param_yt.current_time            = Time[0];
   param_yt.dimensionality          = 3;
   param_yt.refine_by               = 2;
   param_yt.num_fields              = NField;

#  ifdef LIBYT_USE_PATCH_GROUP
   if ( NPatchAllLv % 8 != 0 || NPatchLocalLv % 8 != 0 ) Aux_Error( ERROR_INFO, "Using patch group in libyt failed !!\n" );
   param_yt.num_grids               = NPatchAllLv / 8;
   param_yt.num_grids_local         = NPatchLocalLv / 8;
#  else
   param_yt.num_grids               = NPatchAllLv;
   param_yt.num_grids_local         = NPatchLocalLv;
#  endif

#  ifdef PARTICLE
   yt_species *species_list         = new yt_species [1];
   species_list[0].species_name     = "io";
   species_list[0].num_attr         = PAR_NATT_TOTAL;

   param_yt.num_species             = 1;
   param_yt.species_list            = species_list;
#  endif

   for (int d=0; d<3; d++)
   {
      param_yt.domain_dimensions[d] = NX0_TOT[d];
      param_yt.domain_left_edge [d] = 0.0;
      param_yt.domain_right_edge[d] = amr->BoxSize[d];
      param_yt.periodicity      [d] = ( OPT__BC_FLU[0] == BC_FLU_PERIODIC ) ? 1 : 0;
   }

#  ifdef COMOVING
   param_yt.cosmological_simulation = 1;
   param_yt.current_redshift        = 1.0/Time[0] - 1.0;
   param_yt.omega_lambda            = 1.0 - OMEGA_M0;
   param_yt.omega_matter            = OMEGA_M0;
   param_yt.hubble_constant         = HUBBLE0;
#  else
   param_yt.cosmological_simulation = 0;
   param_yt.current_redshift        = 0.0;
   param_yt.omega_lambda            = 0.0;
   param_yt.omega_matter            = 0.0;
   param_yt.hubble_constant         = 0.0;
#  endif


// 2. transfer simulation information to libyt
   if ( yt_set_parameter( &param_yt ) != YT_SUCCESS )    Aux_Error( ERROR_INFO, "yt_set_parameter() failed !!\n" );

// 2-1. free no longer used resource
#  ifdef PARTICLE
   delete [] species_list;
#  endif

// 3. set code specific parameter
#  ifdef MHD
   const int mhd = 1;
#  else
   const int mhd = 0;
#  endif
   if (yt_add_user_parameter_int("mhd", 1, &mhd) != YT_SUCCESS)  Aux_Error( ERROR_INFO, "yt_add_user_parameter() add mhd failed !!\n" );

#  if ( MODEL == HYDRO )
   const double gamma = (double) GAMMA;
   const double mu = (double) MOLECULAR_WEIGHT;
#  ifdef SRHD
   const int srhd = 1;
#  else
   const int srhd = 0;
#  endif
   if (yt_add_user_parameter_double("gamma", 1, &gamma) != YT_SUCCESS )  Aux_Error( ERROR_INFO, "yt_add_user_parameter() add GAMMA failed !!\n" );
   if (yt_add_user_parameter_double("mu", 1, &mu) != YT_SUCCESS )  Aux_Error( ERROR_INFO, "yt_add_user_parameter() add MOLECULAR_WEIGHT failed !!\n" );
   if (yt_add_user_parameter_int("srhd", 1, &srhd) != YT_SUCCESS ) Aux_Error( ERROR_INFO, "yt_add_user_parameter() add srhd failed !!\n" );

#  elif ( MODEL == ELBDM )
   const int srhd = 0;
   if (yt_add_user_parameter_int("srhd", 1, &srhd) != YT_SUCCESS ) Aux_Error( ERROR_INFO, "yt_add_user_parameter() add srhd failed !!\n" );
#  endif

   if ( OPT__VERBOSE  &&  MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );


} // FUNCTION : YT_SetParameter



#endif // #ifdef SUPPORT_LIBYT
