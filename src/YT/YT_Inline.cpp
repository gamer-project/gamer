#include "GAMER.h"

#ifdef SUPPORT_LIBYT

void YT_SetParameter( const int NPatchAllLv, const int NField, const int NPatchLocalLv, yt_field *FieldList );
void YT_AddLocalGrid( const int *GID_Offset, const int *GID_LvStart, const int (*NPatchAllRank)[NLEVEL], int NField, yt_field *FieldList);



//-------------------------------------------------------------------------------------------------------
// Function    :  YT_Inline
// Description :  Invoke the yt inline analysis
//
// Note        :  1. This function conducts the following three basic steps for performing the yt inline analysis
//                   1-1. YT_SetParameter   --> invoke yt_set_parameter()
//                   1-2. YT_AddLocalGrid   --> invoke yt_get_gridsPtr(), yt_add_grids() for local patches
//                   1-3. yt_inline(), yt_inline_argument
//                   1-4. yt_free_gridsPtr()
//                2. This function is invoked by main() directly
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void YT_Inline()
{

   if ( OPT__VERBOSE  &&  MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );


// 1. gather the number of patches at different MPI ranks, calculate number of local patches
//    and set the corresponding GID offset
   int (*NPatchAllRank)[NLEVEL] = new int [MPI_NRank][NLEVEL];
   int NPatchLocal[NLEVEL], NPatchAllLv=0, NPatchLocalLv=0, GID_Offset[NLEVEL], GID_LvStart[NLEVEL];

   for (int lv=0; lv<NLEVEL; lv++)  
   {
      NPatchLocal[lv] = amr->NPatchComma[lv][1];
      NPatchLocalLv = NPatchLocalLv + NPatchLocal[lv];
   }

   MPI_Allgather( NPatchLocal, NLEVEL, MPI_INT, NPatchAllRank[0], NLEVEL, MPI_INT, MPI_COMM_WORLD );

   for (int lv=0; lv<NLEVEL; lv++)
   {
      GID_Offset[lv] = 0;

      for (int r=0; r<MPI_Rank; r++)      GID_Offset[lv] += NPatchAllRank[r][lv];

      for (int FaLv=0; FaLv<lv; FaLv++)   GID_Offset[lv] += NPatchTotal[FaLv];

      NPatchAllLv += NPatchTotal[lv];

      GID_LvStart[lv] = ( lv == 0 ) ? 0 : GID_LvStart[lv-1] + NPatchTotal[lv-1];
   }


// 2. prepare YT-specific parameters
// 2-1. determine the number of fields
   int NField = NCOMP_TOTAL;
#  ifdef GRAVITY
   int PotIdx = NField;
   NField = NField + 1;
#  endif

#ifdef MHD
   int MHDIdx = NField;
   NField = NField + NCOMP_MAG;
#endif

// 2-2. determine the field labels, and declare a yt_field type array
//      which stores the field labels and field define type (ex: cell-centered, face-centered)
   yt_field *FieldList = new yt_field [NField];
   for (int v=0; v<NCOMP_TOTAL; v++){
       FieldList[v].field_name = FieldLabel[v];
   }

#ifdef GRAVITY
   FieldList[PotIdx].field_name = PotLabel;
#endif

#ifdef MHD
   for (int v=0; v<NCOMP_MAG; v++){
       FieldList[v + MHDIdx].field_name        = MagLabel[v];
       FieldList[v + MHDIdx].field_define_type = "face-centered";
       FieldList[v + MHDIdx].field_unit        = "code_magnetic";
   }
   FieldList[ MHDIdx ].field_display_name = "B_x";
   FieldList[ MHDIdx + 1 ].field_display_name = "B_y";
   FieldList[ MHDIdx + 2 ].field_display_name = "B_z";

   // Add alias names to "MAGX"
   char *MHDX_name_alias[2] = {"test_alias_name", "magnetic_x"};
   FieldList[ MHDIdx ].num_field_name_alias = 2;
   FieldList[ MHDIdx ].field_name_alias = MHDX_name_alias;

#endif

// 2-3. Call YT_SetParameter
   YT_SetParameter( NPatchAllLv, NField, NPatchLocalLv, FieldList );


// 3. prepare local patches for libyt
   YT_AddLocalGrid( GID_Offset, GID_LvStart, NPatchAllRank, NField, FieldList);


// 4. perform yt inline analysis
   if ( yt_inline_argument( "yt_inline_inputArg", 1, "\'Dens\'" ) != YT_SUCCESS )    Aux_Error( ERROR_INFO, "yt_inline_inputArg() failed !!\n" );
   if ( yt_inline( "yt_inline" ) != YT_SUCCESS )     Aux_Error( ERROR_INFO, "yt_inline() failed !!\n" );

// 5. free resource
   if ( yt_free_gridsPtr() != YT_SUCCESS )    Aux_Error( ERROR_INFO, "yt_free_gridsPtr() failed !!\n" );
   delete [] NPatchAllRank;
   delete [] FieldList;

   if ( OPT__VERBOSE  &&  MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : YT_Inline



#endif // #ifdef SUPPORT_LIBYT
