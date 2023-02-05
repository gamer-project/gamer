#include "GAMER.h"

#ifdef SUPPORT_LIBYT

// call libyt API
void YT_SetParameter( const int NPatchAllLv, const int NField, const int NPatchLocalLv);
void YT_AddLocalGrid( const int *GID_LvStart, const int (*NPatchAllRank)[NLEVEL], int NField, yt_field *FieldList);

#ifdef LIBYT_USE_PATCH_GROUP

void DerivedFuncWithName_PatchGroup(const int list_len, const long *list_gid, const char *field, yt_array *data_array);

#else  // #ifdef LIBYT_USE_PATCH_GROUP

#ifdef MHD
// derived function for Mag to CCMag
void MagX_DerivedFunc(const int list_len, const long *list_gid, const char *field_name, yt_array *data_array);
void MagY_DerivedFunc(const int list_len, const long *list_gid, const char *field_name, yt_array *data_array);
void MagZ_DerivedFunc(const int list_len, const long *list_gid, const char *field_name, yt_array *data_array);
#endif

#if ( MODEL == HYDRO )
void Temperature_DerivedFunc(const int list_len, const long *list_gid, const char *field_name, yt_array *data_array);
#endif

#endif  // #ifdef LIBYT_USE_PATCH_GROUP

#ifdef PARTICLE
// get the particle attribute, since we only have one type of particle "io"
// we only need one function.
void Get_ParticleAttribute(const int list_len, const long *list_gid, const char *ptype, const char *attr, yt_array *data_array);
#endif



//-------------------------------------------------------------------------------------------------------
// Function    :  YT_Inline
// Description :  Invoke the yt inline analysis
//
// Note        :  1. This function conducts the following to perform the yt inline analysis
//                   1-1. set yt parameters
//                   1-2. set fields and particles info
//                   1-3. set local grids info
//                   1-4. run in situ analysis
//                   1-5. free resources allocated for in situ analysis
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
   int EoSIdx = NField;
#  ifdef LIBYT_USE_PATCH_GROUP
   // Add field : _TEMP, _PRES, _ENTR
   NField = NField + 3;
#  else
   // Add field : _TEMP
   NField = NField + 1;
#  endif // #ifdef LIBYT_USE_PATCH_GROUP
#  endif

// 2-2. Call YT_SetParameter and set particle info if need.
   YT_SetParameter( NPatchAllLv, NField, NPatchLocalLv);

// 3.   Get FieldList and ParticleList, fill the info if needed
// 3-1. get yt_field array FieldList, and filled in field info
//      FieldList :
//      |                    |                  |   LIBYT_USE_PATCH_GROUP  |
//      +-------------------------------------------------------------------
//      +       0            |                  +                          +
//      +       :            |   cell-centered  +                          +
//      + (NCOMP_TOTAL - 1)  |                  +       derived_func       +
//      +  GRAVITY (PotIdx)  |   cell-centered  +                          +
//      +  MHD     (MHDIdx)  |   face-centered  +                          +
//      +.......................................+..........................+
//      +                      Other Derived Fields                        +
//      +---------------------------------------+--------------------------+
   yt_field *FieldList;
   yt_get_FieldsPtr( &FieldList );

#  ifdef LIBYT_USE_PATCH_GROUP
   for (int v=0; v<NCOMP_TOTAL; v++){
       FieldList[v].field_name   = FieldLabel[v];
       FieldList[v].field_type   = "derived_func";
       FieldList[v].derived_func = DerivedFuncWithName_PatchGroup;
   }

#  ifdef GRAVITY
   FieldList[PotIdx].field_name   = PotLabel;
   FieldList[PotIdx].field_type   = "derived_func";
   FieldList[PotIdx].derived_func = DerivedFuncWithName_PatchGroup;
#  endif

#  ifdef MHD
   const char *CCMagLabel[] = {"CCMagX", "CCMagY", "CCMagZ"};
   for (int v=0; v<NCOMP_MAG; v++){
       FieldList[v + MHDIdx].field_name   = CCMagLabel[v];
       FieldList[v + MHDIdx].field_type   = "derived_func";
       FieldList[v + MHDIdx].field_unit   = "code_magnetic";
       FieldList[v + MHDIdx].derived_func = DerivedFuncWithName_PatchGroup;
   }

   // Add field display name
   FieldList[ MHDIdx     ].field_display_name = "B_x";
   FieldList[ MHDIdx + 1 ].field_display_name = "B_y";
   FieldList[ MHDIdx + 2 ].field_display_name = "B_z";
#  endif

#  if ( MODEL == HYDRO )
   const char *AddFieldLabel[] = {"Temp", "Pres", "Entr"};
   for (int v=0; v<3; v++){
       FieldList[v + EoSIdx].field_name   = AddFieldLabel[v];
       FieldList[v + EoSIdx].field_type   = "derived_func";
       FieldList[v + EoSIdx].derived_func = DerivedFuncWithName_PatchGroup;
   }

   FieldList[EoSIdx].field_unit             = "code_temperature";
   FieldList[EoSIdx].field_display_name     = "Temperature";

   FieldList[EoSIdx + 1].field_unit         = "code_mass / (code_length*code_time**2)";
   FieldList[EoSIdx + 1].field_display_name = "Pressure";

#  if ( EOS == EOS_NUCLEAR )
   FieldList[EoSIdx + 2].field_unit         = "code_mass*code_length**(2) / (code_temperature*code_time**(2))";
#  endif // #if ( EOS == EOS_NUCLEAR )
#  if ( EOS == EOS_GAMMA )
   char EntropyUnit[100];
   real gamma_m1 = (real) GAMMA - 1.0;
   sprintf(EntropyUnit, "code_mass**(1-%.2f) / (code_length**(1-3*%.2f)*code_time**2)", gamma_m1, gamma_m1);
   FieldList[EoSIdx + 2].field_unit         = EntropyUnit;
#  endif // #if ( EOS == EOS_GAMMA )
   FieldList[EoSIdx + 2].field_display_name = "Entropy";

#  endif // #if ( MODEL == HYDRO )

#  else  // #ifdef LIBYT_USE_PATCH_GROUP
   for (int v=0; v<NCOMP_TOTAL; v++){
       FieldList[v].field_name = FieldLabel[v];
   }

#  ifdef GRAVITY
   FieldList[PotIdx].field_name = const_cast<char*> (PotLabel);
#  endif

#  ifdef MHD
   char *CCMagLabel[] = {"CCMagX", "CCMagY", "CCMagZ"};
   for (int v=0; v<NCOMP_MAG; v++){
       FieldList[v + MHDIdx].field_name  = CCMagLabel[v];
       FieldList[v + MHDIdx].field_type  = "face-centered";
       FieldList[v + MHDIdx].field_unit  = "code_magnetic";
   }

   // Add field display name
   FieldList[ MHDIdx     ].field_display_name = "B_x";
   FieldList[ MHDIdx + 1 ].field_display_name = "B_y";
   FieldList[ MHDIdx + 2 ].field_display_name = "B_z";

   // Add field derived function pointer
   // if you wish to use, set field_define_type = "derived_field"
   FieldList[ MHDIdx     ].derived_func = MagX_DerivedFunc;
   FieldList[ MHDIdx + 1 ].derived_func = MagY_DerivedFunc;
   FieldList[ MHDIdx + 2 ].derived_func = MagZ_DerivedFunc;
#  endif

#  if ( MODEL == HYDRO )
   FieldList[EoSIdx].field_name         = "Temp";
   FieldList[EoSIdx].field_type         = "derived_func";
   FieldList[EoSIdx].field_unit         = "code_temperature";
   FieldList[EoSIdx].field_display_name = "Temperature";
   FieldList[EoSIdx].derived_func       = Temperature_DerivedFunc;
#  endif

#  endif // #ifdef LIBYT_USE_PATCH_GROUP

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
   yt_get_ParticlesPtr( &ParticleList );

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
   ParticleList[0].get_par_attr = Get_ParticleAttribute;

#  endif // #ifdef PARTICLE

// 4. prepare local patches for libyt
   YT_AddLocalGrid( GID_LvStart, NPatchAllRank, NField, FieldList);

// 5. perform yt inline analysis
   if ( yt_run_FunctionArguments( "yt_inline_inputArg", 1, "\'Dens\'" ) != YT_SUCCESS )
       Aux_Error( ERROR_INFO, "yt_run_FunctionArguments() failed !!\n" );
   if ( yt_run_Function( "yt_inline" ) != YT_SUCCESS )
       Aux_Error( ERROR_INFO, "yt_run_Function() failed !!\n" );

// 6. free resource
   if ( yt_free() != YT_SUCCESS )    Aux_Error( ERROR_INFO, "yt_free() failed !!\n" );
   delete [] NPatchAllRank;

   if ( OPT__VERBOSE  &&  MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : YT_Inline



#endif // #ifdef SUPPORT_LIBYT
