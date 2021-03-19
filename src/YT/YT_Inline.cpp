#include "GAMER.h"

#ifdef SUPPORT_LIBYT

void YT_SetParameter( const int NPatchAllLv, const int NField, const int NPatchLocalLv, char **FieldLabel );
void YT_AddLocalGrid( const int *GID_Offset, const int *GID_LvStart, const int (*NPatchAllRank)[NLEVEL]);




//-------------------------------------------------------------------------------------------------------
// Function    :  YT_Inline
// Description :  Invoke the yt inline analysis
//
// Note        :  1. This function conducts the following three basic steps for performing the yt inline analysis
//                   1-1. YT_SetParameter --> invoke yt_set_parameter()
//                   1-2. YT_AddLocalGrid   --> invoke yt_get_gridsPtr(), yt_add_grids() for local patches
//                   1-3. yt_inline()
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
   NField = NField + 1;
#  endif

// 2-2. determine the field labels
   char **FieldLabelForYT = new char* [NField];
   for (int v=0; v<NField; v++)  { FieldLabelForYT[v] = new char [MAX_STRING]; }
   for (int v=0; v<NCOMP_TOTAL; v++)   { sprintf( FieldLabelForYT[v], FieldLabel[v] ); }

#  ifdef GRAVITY
   sprintf( FieldLabelForYT[NCOMP_TOTAL], PotLabel );
#  endif

// 2-3. Call YT_SetParameter
   YT_SetParameter( NPatchAllLv, NField, NPatchLocalLv, FieldLabelForYT );


// 3. prepare local patches for libyt
   YT_AddLocalGrid( GID_Offset, GID_LvStart, NPatchAllRank );


// 4. perform yt inline analysis
   if ( yt_inline( "yt_inline_inputArg", 1, "\'Dens\'" ) != YT_SUCCESS )    Aux_Error( ERROR_INFO, "yt_inline() failed !!\n" );


// 5. free resource
   if ( yt_free_gridsPtr() != YT_SUCCESS )    Aux_Error( ERROR_INFO, "yt_free_gridsPtr() failed !!\n" );
   delete [] NPatchAllRank;
   for (int v=0; v<NField; v++)  delete [] FieldLabelForYT[v];
   delete [] FieldLabelForYT;

   if ( OPT__VERBOSE  &&  MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : YT_Inline



#endif // #ifdef SUPPORT_LIBYT
