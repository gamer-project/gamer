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
// Parameter   :  NField        : Number of fields loaded to YT.
//                FieldList     : List of field_name, field_define_type.
//                pc            : Patch count object with information about the number of patches on all ranks
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void YT_AddLocalGrid( int NField, yt_field *FieldList, LB_PatchCount& pc)
{

   if ( OPT__VERBOSE  &&  MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );

// get the libyt local grids array pointer
   yt_grid *YT_Grids;
   yt_get_GridsPtr( &YT_Grids );

// record local grids index and patched grids index if LIBYT_USE_PATCH_GROUP
   int LID = 0;

   LB_LocalPatchExchangeList lel;

// sync load balance ids
   LB_AllgatherLBIdx( pc, lel );
   LB_FillLocalPatchExchangeList( pc, lel );

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
         const long GID = PID + pc.GID_Offset[lv];

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
//       input particle num in this patch group.
         long particle_count = 0;
         for(int i=PID; i<PID+8; i++){
             particle_count += (long) amr->patch[0][lv][i]->NPar;
         }
         YT_Grids[LID].par_count_list[0] = particle_count;
#        else
//       input particle num in this grid
         YT_Grids[LID].par_count_list[0] = (long) amr->patch[0][lv][PID]->NPar;
#        endif // #ifdef LIBYT_USE_PATCH_GROUP
#        endif // #ifdef PARTICLE

#        ifdef LIBYT_USE_PATCH_GROUP
         YT_Grids[LID].id     = (long) GID / 8;
#        else
         YT_Grids[LID].id     = (long) GID;
#        endif // #ifdef LIBYT_USE_PATCH_GROUP

         YT_Grids[LID].level  = lv;

//       getting parent's id
         long FaGID = lel.FaList_Local[lv][PID];

#        ifdef LIBYT_USE_PATCH_GROUP
         if ( FaGID != -1 ) FaGID = (long) (FaGID / 8);
#        endif // #ifdef LIBYT_USE_PATCH_GROUP

         YT_Grids[LID].parent_id = FaGID;

#        ifndef LIBYT_USE_PATCH_GROUP
//       load patch data to libyt if not use LIBYT_USE_PATCH_GROUP
         for (int v=0; v<NCOMP_TOTAL; v++)
         {
             YT_Grids[LID].field_data[v].data_ptr = amr->patch[FluSg][lv][PID]->fluid[v];
         }

#        ifdef GRAVITY
//       find field index of GRAVITY
         int PotIdx = 0;
         for (int v=0; v<NField; v++)
         {
             if ( strcmp( FieldList[v].field_name, PotLabel ) == 0 )
             {
                 PotIdx = v;
                 break;
             }
         }
//       load Pote patch data to libyt
         YT_Grids[LID].field_data[PotIdx].data_ptr = amr->patch[PotSg][lv][PID]->pot;
#        endif // #ifdef GRAVITY

#        ifdef MHD
//       find field index of CCMagX
         int MHDIdx = 0;
         for (int v=0; v<NField; v++)
         {
             if ( strcmp( FieldList[v].field_name, "CCMagX" ) == 0 )
             {
                 MHDIdx = v;
                 break;
             }
         }

         for (int v=0; v<NCOMP_MAG; v++)
         {
//           input the data pointer
             YT_Grids[LID].field_data[ MHDIdx + v ].data_ptr = amr->patch[MagSg][lv][PID]->magnetic[v];

//           input the field dimension, since MHD has different dimension.
             for (int d=0; d<3; d++)
             {
                 YT_Grids[LID].field_data[ MHDIdx + v ].data_dimensions[d] = PATCH_SIZE;

                 if      ( strcmp( FieldList[ MHDIdx + v ].field_name, "CCMagX" ) == 0  &&  d == 2 )
                 {
                     YT_Grids[LID].field_data[ MHDIdx + v ].data_dimensions[d] = PATCH_SIZE + 1;
                 }
                 else if ( strcmp( FieldList[ MHDIdx + v ].field_name, "CCMagY" ) == 0  &&  d == 1 )
                 {
                     YT_Grids[LID].field_data[ MHDIdx + v ].data_dimensions[d] = PATCH_SIZE + 1;
                 }
                 else if ( strcmp( FieldList[ MHDIdx + v ].field_name, "CCMagZ" ) == 0  &&  d == 0 )
                 {
                     YT_Grids[LID].field_data[ MHDIdx + v ].data_dimensions[d] = PATCH_SIZE + 1;
                 }
             }
         }
#        endif // #ifdef MHD

#        endif // #ifndef LIBYT_USE_PATCH_GROUP

         LID = LID + 1;
      }
   } // for (int lv=0; lv<NLEVEL; lv++)

   if ( yt_commit( ) != YT_SUCCESS )  Aux_Error( ERROR_INFO, "yt_commit() failed !!\n" );

   if ( OPT__VERBOSE  &&  MPI_Rank == 0 )   Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : YT_AddLocalGrid

#endif // #ifdef SUPPORT_LIBYT
