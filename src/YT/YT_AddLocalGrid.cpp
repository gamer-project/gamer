#include "GAMER.h"

#ifdef SUPPORT_LIBYT

//-------------------------------------------------------------------------------------------------------
// Function    :  YT_AddLocalGrid
// Description :  Send the hierarchy information and data of local patches to libyt.
//
// Note        :  1. One must call YT_SetParameter() before invoking this function.
//                2. Invoked by YT_Inline().
//                3. FieldList is used by MHD field, since it needs to load dimensions to yt_field.
//
// Parameter   :  GID_LvStart   : Glocal patch index that this level starts at
//                NPatchAllRank : Number of patches in [MPI rank][level]
//                NField        : Number of fields loaded to YT.
//                FieldList     : List of field_name, field_define_type.
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void YT_AddLocalGrid( const int *GID_LvStart, const int (*NPatchAllRank)[NLEVEL], int NField, yt_field *FieldList)
{

   if ( OPT__VERBOSE  &&  MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );

// get the libyt local grids array pointer
   yt_grid *YT_Grids;
   yt_get_gridsPtr( &YT_Grids );

// record local grids index and patched grids index if LIBYT_USE_PATCH_GROUP
   int LID = 0;

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
