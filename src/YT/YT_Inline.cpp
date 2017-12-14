#include "GAMER.h"

#ifdef SUPPORT_LIBYT

void YT_SetParameter( const int NPatchAllLv );
void YT_AddAllGrid( yt_grid *Grid, const int *GID_Offset, const int NField, char **FieldLabel );




//-------------------------------------------------------------------------------------------------------
// Function    :  YT_Inline
// Description :  Invoke the yt inline analysis
//
// Note        :  1. This function conducts the following three basic steps for performing the yt inline analysis
//                   1-1. YT_SetParameter --> invoke yt_set_parameter()
//                   1-2. YT_AddAllGrid   --> invoke yt_add_grid() for all patches
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


// 1. gather the number of patches at different MPI ranks and set the corresponding GID offset
   int (*NPatchAllRank)[NLEVEL] = new int [MPI_NRank][NLEVEL];
   int NPatchLocal[NLEVEL], NPatchAllLv=0, GID_Offset[NLEVEL];

   for (int lv=0; lv<NLEVEL; lv++)  NPatchLocal[lv] = amr->NPatchComma[lv][1];

   MPI_Allgather( NPatchLocal, NLEVEL, MPI_INT, NPatchAllRank[0], NLEVEL, MPI_INT, MPI_COMM_WORLD );

   for (int lv=0; lv<NLEVEL; lv++)
   {
      GID_Offset[lv] = 0;

      for (int r=0; r<MPI_Rank; r++)      GID_Offset[lv] += NPatchAllRank[r][lv];

      for (int FaLv=0; FaLv<lv; FaLv++)   GID_Offset[lv] += NPatchTotal[FaLv];

      NPatchAllLv += NPatchTotal[lv];
   }


// 2. prepare YT-specific parameters
   YT_SetParameter( NPatchAllLv );


// 3. prepare the hierarchy information and data of all patches
// 3-1. determine the number of fields
   int NField = NCOMP_TOTAL;
#  ifdef GRAVITY
   NField ++;
#  endif

// 3-2. set the field labels
   char **FieldLabel = new char* [NField];
   for (int v=0; v<NField; v++)  FieldLabel[v] = new char [MAX_STRING];

#  if ( MODEL == HYDRO )
   sprintf( FieldLabel[DENS], "Dens" );
   sprintf( FieldLabel[MOMX], "MomX" );
   sprintf( FieldLabel[MOMY], "MomY" );
   sprintf( FieldLabel[MOMZ], "MomZ" );
   sprintf( FieldLabel[ENGY], "Engy" );

#  elif ( MODEL == MHD )
#  warning : WAIT MHD !!!

#  elif ( MODEL == ELBDM )
   sprintf( FieldLabel[DENS], "Dens" );
   sprintf( FieldLabel[REAL], "Real" );
   sprintf( FieldLabel[IMAG], "Imag" );

#  else
#  error : unsupported MODEL !!
#  endif // MODEL

   for (int v=0; v<NCOMP_PASSIVE; v++)    sprintf( FieldLabel[NCOMP_FLUID+v], "%s", PassiveFieldName_Grid[v] );

#  ifdef GRAVITY
   sprintf( FieldLabel[NCOMP_TOTAL], "Pote" );
#  endif

// 3-3. prepare all patches for libyt
   yt_grid *Grid = new yt_grid [NPatchAllLv];

   for (int GID=0; GID<NPatchAllLv; GID++)   Grid[GID].field_data = new void* [NField];

   YT_AddAllGrid( Grid, GID_Offset, NField, FieldLabel );


// 4. perform yt inline analysis
   if ( yt_inline() != YT_SUCCESS )    Aux_Error( ERROR_INFO, "yt_inline() failed !!\n" );


// 5. free resource
   delete [] NPatchAllRank;
   for (int v=0; v<NField; v++)  delete [] FieldLabel[v];
   delete [] FieldLabel;
   for (int GID=0; GID<NPatchAllLv; GID++)   delete [] Grid[GID].field_data;
   delete [] Grid;


   if ( OPT__VERBOSE  &&  MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : YT_Inline



#endif // #ifdef SUPPORT_LIBYT
