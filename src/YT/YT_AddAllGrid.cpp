#include "GAMER.h"

#ifdef SUPPORT_LIBYT




//-------------------------------------------------------------------------------------------------------
// Function    :  YT_AddAllGrid
// Description :  Send the hierarchy information and data of all patches to libyt
//
// Note        :  1. One must call YT_SetParameter() before invoking this function
//                2. Invoked by YT_Inline()
//
// Parameter   :  GID_Offset : Global patch index offset at each refinement level for this rank
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void YT_AddAllGrid( const int *GID_Offset )
{

   if ( OPT__VERBOSE  &&  MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );

// get the libyt local grids array pointer
   yt_grid *YT_Grids;
   yt_get_gridsPtr( &YT_Grids );

// loop over local patches at all levels
   for (int lv=0; lv<NLEVEL; lv++)
   {
      const int FaLv  = lv - 1;
      const int FluSg = amr->FluSg[lv];
#     ifdef GRAVITY
      const int PotSg = amr->PotSg[lv];
#     endif

      for (int PID=0; PID<(amr->NPatchComma[lv][1]); PID++)
      {
         const int GID = PID + GID_Offset[lv];

         for (int d=0; d<3; d++)
         {
            YT_Grids[PID].left_edge [d] = amr->patch[0][lv][PID]->EdgeL[d];
            YT_Grids[PID].right_edge[d] = amr->patch[0][lv][PID]->EdgeR[d];
            YT_Grids[PID].dimensions[d] = PATCH_SIZE;
         }

#        ifdef PARTICLE
         YT_Grids[PID].particle_count = amr->patch[0][lv][PID]->NPar;
#        else
         YT_Grids[PID].particle_count = 0;
#        endif

         YT_Grids[PID].id             = GID;
         YT_Grids[PID].parent_id      = ( amr->patch[0][lv][PID]->father < 0 ) ? -1 : amr->patch[0][lv][PID]->father + GID_Offset[FaLv];
         YT_Grids[PID].level          = lv;

         for (int v = 0; v < NCOMP_TOTAL; v++){
            YT_Grids[PID].field_data[v]           = amr->patch[FluSg][lv][PID]->fluid[v];
         }

#        ifdef GRAVITY
         YT_Grids[PID].field_data[NCOMP_TOTAL] = amr->patch[PotSg][lv][PID]->pot;
#        endif
      }
   }

   if ( yt_add_grids( ) != YT_SUCCESS )  Aux_Error( ERROR_INFO, "yt_add_grid() failed !!\n" );

   if ( OPT__VERBOSE  &&  MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : YT_AddAllGrid



#endif // #ifdef SUPPORT_LIBYT
