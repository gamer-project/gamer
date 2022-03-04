#include "GAMER.h"

#ifdef SUPPORT_LIBYT



void YT_SetParameter( const int NPatchAllLv, const int NField, const int NPatchLocalLv);
void YT_AddLocalGrid( const int *GID_Offset, const int *GID_LvStart, const int (*NPatchAllRank)[NLEVEL], int NField, yt_field *FieldList);

#ifdef MHD
// derived function for Mag to CCMag
void MagX_DerivedFunc(long gid, double *Converted_MagX);
void MagY_DerivedFunc(long gid, double *Converted_MagY);
void MagZ_DerivedFunc(long gid, double *Converted_MagZ);
#endif

#if ( MODEL == HYDRO )
void Temperature_DerivedFunc(long gid, double *TempData);
#endif

#ifdef PARTICLE
// get the particle attribute, since we only have one type of particle "io"
// we only need one function.
void Get_ParticleAttribute(long gid, char *attr, void *raw_data);
#endif




//-------------------------------------------------------------------------------------------------------
// Function    :  YT_Inline
// Description :  Invoke the yt inline analysis
//
// Note        :  1. This function conducts the following three basic steps for performing the yt inline analysis
//                   1-1. YT_SetParameter   --> invoke yt_set_parameter()
//                   1-2. yt_get_fieldsPtr, yt_get_particlesPtr  --> get the yt_field array pointer and
//                        yt_particle array pointer, then fill in the info.
//                   1-3. YT_AddLocalGrid   --> invoke yt_get_gridsPtr(), yt_add_grids() for local patches
//                   1-4. yt_inline(), yt_inline_argument()
//                   1-5. yt_free_gridsPtr()
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

#  ifdef MHD
   int MHDIdx = NField;
   NField = NField + NCOMP_MAG;
#  endif

#  if ( MODEL == HYDRO )
   int EoSTempIdx = NField;
   NField = NField + 1;
#  endif

// 2-2. Call YT_SetParameter and set particle info if need.
   YT_SetParameter( NPatchAllLv, NField, NPatchLocalLv);

// 3.   Get FieldList and ParticleList, fill the info if needed
// 3-1. get yt_field array FieldList, and filled in field info
//      FieldList :
//      +----------------------------------------
//      +       0            |                  +
//      +       :            |   cell-centered  +
//      + (NCOMP_TOTAL - 1)  |                  +
//      +  GRAVITY (PotIdx)  |   cell-centered  +
//      +  MHD     (MHDIdx)  |   face-centered  +
//      +.......................................+
//      +         Other Derived Fields          +
//      +---------------------------------------+
   yt_field *FieldList;
   yt_get_fieldsPtr( &FieldList );

   for (int v=0; v<NCOMP_TOTAL; v++){
       FieldList[v].field_name = FieldLabel[v];
   }

#  ifdef GRAVITY
   FieldList[PotIdx].field_name = const_cast<char*> (PotLabel);
#  endif

#  ifdef MHD
   char *CCMagLabel[] = {"CCMagX", "CCMagY", "CCMagZ"};
   for (int v=0; v<NCOMP_MAG; v++){
       FieldList[v + MHDIdx].field_name        = CCMagLabel[v];
       FieldList[v + MHDIdx].field_define_type = "face-centered";
       FieldList[v + MHDIdx].field_unit        = "code_magnetic";
   }

   // Add field display name
   FieldList[ MHDIdx     ].field_display_name = "B_x";
   FieldList[ MHDIdx + 1 ].field_display_name = "B_y";
   FieldList[ MHDIdx + 2 ].field_display_name = "B_z";

   // Add field derived function pointer
   FieldList[ MHDIdx     ].derived_func = MagX_DerivedFunc;
   FieldList[ MHDIdx + 1 ].derived_func = MagY_DerivedFunc;
   FieldList[ MHDIdx + 2 ].derived_func = MagZ_DerivedFunc;
#  endif

#  if ( MODEL == HYDRO )
   FieldList[EoSTempIdx].field_name = "Temp";
   FieldList[EoSTempIdx].field_define_type = "derived_func";
   FieldList[EoSTempIdx].field_unit = "code_temperature";
   FieldList[EoSTempIdx].field_display_name = "Temperature";
   FieldList[EoSTempIdx].derived_func = Temperature_DerivedFunc;
#  endif

   // Set field's data type
   for (int v=0; v<NField; v++){
#  ifdef FLOAT8
       FieldList[v].field_dtype = YT_DOUBLE;
#  else
       FieldList[v].field_dtype = YT_FLOAT;
#  endif
   }

// 3-2 Get the ParticleList
#  ifdef PARTICLE
   yt_particle *ParticleList;
   yt_get_particlesPtr( &ParticleList );

   // Set attributes
   for (int v=0; v<ParticleList[0].num_attr; v++){
       // set attribute name
       ParticleList[0].attr_list[v].attr_name  = ParAttLabel[v];
       // set attribute data type
#      ifdef FLOAT8
       ParticleList[0].attr_list[v].attr_dtype = YT_DOUBLE;
#      else
       ParticleList[0].attr_list[v].attr_dtype = YT_FLOAT;
#      endif
       // set attribute unit
   }

   // Set label (attribute name) of coordinate x/y/z
   ParticleList[0].coor_x   = "ParPosX";
   ParticleList[0].coor_y   = "ParPosY";
   ParticleList[0].coor_z   = "ParPosZ";

   // Set get attribute function
   ParticleList[0].get_attr = Get_ParticleAttribute;
#  endif // #ifdef PARTICLE

// 4. prepare local patches for libyt
   YT_AddLocalGrid( GID_Offset, GID_LvStart, NPatchAllRank, NField, FieldList);

// 5. perform yt inline analysis
   if ( yt_inline_argument( "yt_inline_inputArg", 1, "\'Dens\'" ) != YT_SUCCESS )    Aux_Error( ERROR_INFO, "yt_inline_inputArg() failed !!\n" );
   if ( yt_inline( "yt_inline" ) != YT_SUCCESS )     Aux_Error( ERROR_INFO, "yt_inline() failed !!\n" );

// 6. free resource
   if ( yt_free_gridsPtr() != YT_SUCCESS )    Aux_Error( ERROR_INFO, "yt_free_gridsPtr() failed !!\n" );
   delete [] NPatchAllRank;

   if ( OPT__VERBOSE  &&  MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : YT_Inline



#endif // #ifdef SUPPORT_LIBYT
