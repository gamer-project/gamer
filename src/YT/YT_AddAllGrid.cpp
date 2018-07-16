#include "GAMER.h"

#ifdef SUPPORT_LIBYT




//-------------------------------------------------------------------------------------------------------
// Function    :  YT_AddAllGrid
// Description :  Send the hierarchy information and data of all patches to libyt
//
// Note        :  1. One must call YT_SetParameter() before invoking this function
//                2. Invoked by YT_Inline()
//
// Parameter   :  Grid       : libyt grid object
//                GID_Offset : Global patch index offset at each refinement level for this rank
//                NField     : Number of fields (e.g., density, momentum-x, ...)
//                FieldLabel : Labels of all fields
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void YT_AddAllGrid( yt_grid *Grid, const int *GID_Offset, const int NField, char **FieldLabel )
{

   if ( OPT__VERBOSE  &&  MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );


// loop over all patches at all levels
   for (int lv=0; lv<NLEVEL; lv++)
   {
      const int FaLv  = lv - 1;
      const int FluSg = amr->FluSg[lv];
#     ifdef GRAVITY
      const int PotSg = amr->PotSg[lv];
#     endif

      for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
      {
//       1. collect the hierarchy information of this patch
         const int GID = PID + GID_Offset[lv];

         for (int d=0; d<3; d++)
         {
            Grid[GID].left_edge [d] = amr->patch[0][lv][PID]->EdgeL[d];
            Grid[GID].right_edge[d] = amr->patch[0][lv][PID]->EdgeR[d];
            Grid[GID].dimensions[d] = PATCH_SIZE;
         }

#        ifdef PARTICLE
         Grid[GID].particle_count = amr->patch[0][lv][PID]->NPar;
#        else
         Grid[GID].particle_count = 0;
#        endif
         Grid[GID].id             = GID;
//###ISSUE: do not support parallelism
//          --> see Output/Output_DumpData_Total_HDF5.cpp for the parallel implementation
         Grid[GID].parent_id      = ( amr->patch[0][lv][PID]->father < 0 ) ? -1 : amr->patch[0][lv][PID]->father + GID_Offset[FaLv];
         Grid[GID].level          = lv;


//       2. set pointers pointing to different field data
//          --> "field_data" pointer array must be pre-allocated
         for (int v=0; v<NCOMP_TOTAL; v++)
         Grid[GID].field_data[v]           = amr->patch[FluSg][lv][PID]->fluid[v];
#        ifdef GRAVITY
         Grid[GID].field_data[NCOMP_TOTAL] = amr->patch[PotSg][lv][PID]->pot;
#        endif


//       3. set other field parameters
         Grid[GID].num_fields   = NField;
         Grid[GID].field_labels = (const char **)FieldLabel;
#        ifdef FLOAT8
         Grid[GID].field_ftype  = YT_DOUBLE;
#        else
         Grid[GID].field_ftype  = YT_FLOAT;
#        endif


//       4. send this patch to libyt
         if ( yt_add_grid( &Grid[GID] ) != YT_SUCCESS )  Aux_Error( ERROR_INFO, "yt_add_grid() failed !!\n" );

      } // for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
   } // for (int lv=0; lv<NLEVEL; lv++)


   if ( OPT__VERBOSE  &&  MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : YT_AddAllGrid



#endif // #ifdef SUPPORT_LIBYT
