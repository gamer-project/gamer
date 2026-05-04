#include "GAMER.h"

#ifdef SUPPORT_LIBYT

// call libyt API
void YT_SetParameter( const int NPatchAllLv, const int NField, const int NPatchLocalAllLv );
void YT_AddLocalGrid( int NField, yt_field *FieldList, LB_PatchCount& pc );

#ifdef LIBYT_USE_PATCH_GROUP

void DerivedFuncWithName_PatchGroup( const int list_len, const long *list_gid, const char *field, yt_array *data_array );

#else  // #ifdef LIBYT_USE_PATCH_GROUP

#ifdef MHD
// derived function for Mag to CCMag
void MagX_DerivedFunc( const int list_len, const long *list_gid, const char *field_name, yt_array *data_array );
void MagY_DerivedFunc( const int list_len, const long *list_gid, const char *field_name, yt_array *data_array );
void MagZ_DerivedFunc( const int list_len, const long *list_gid, const char *field_name, yt_array *data_array );
#endif

#if ( MODEL == HYDRO )
void Temperature_DerivedFunc( const int list_len, const long *list_gid, const char *field_name, yt_array *data_array );
#endif

#endif  // #ifdef LIBYT_USE_PATCH_GROUP

#ifdef PARTICLE
// get the particle attribute, since we only have one type of particle "io"
// we only need one function.
void Get_ParticleAttribute( const int list_len, const long *list_gid, const char *ptype, const char *attr, yt_array *data_array );
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

   LB_PatchCount pc;

// 1. get patch counts per level and per rank from all ranks
   LB_AllgatherPatchCount( pc );


// set YT_GID_Offset for searching GID in derived function and particle get attribute function.
   for (int lv=0; lv<NLEVEL; lv++)
   {
      YT_GID_Offset[lv] = pc.GID_Offset[lv];
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
// add field : _TEMP, _PRES, _ENTR
   NField = NField + 3;
#  else
// add field : _TEMP
   NField = NField + 1;
#  endif // #ifdef LIBYT_USE_PATCH_GROUP
#  endif

// 2-2. call YT_SetParameter and set particle info if need.
   YT_SetParameter( pc.NPatchAllLv, NField, pc.NPatchLocalAllLv );

// 3.   get FieldList and ParticleList, fill the info if needed
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
   for (int v=0; v<NCOMP_TOTAL; v++)
   {
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
   for (int v=0; v<NCOMP_MAG; v++)
   {
      FieldList[v + MHDIdx].field_name   = CCMagLabel[v];
      FieldList[v + MHDIdx].field_type   = "derived_func";
      FieldList[v + MHDIdx].field_unit   = "code_magnetic";
      FieldList[v + MHDIdx].derived_func = DerivedFuncWithName_PatchGroup;
   }

// add field display name
   FieldList[ MHDIdx     ].field_display_name = "B_x";
   FieldList[ MHDIdx + 1 ].field_display_name = "B_y";
   FieldList[ MHDIdx + 2 ].field_display_name = "B_z";
#  endif

#  if ( MODEL == HYDRO )
   const char *AddFieldLabel[] = {"Temp", "Pres", "Entr"};
   for (int v=0; v<3; v++)
   {
      FieldList[v + EoSIdx].field_name   = AddFieldLabel[v];
      FieldList[v + EoSIdx].field_type   = "derived_func";
      FieldList[v + EoSIdx].derived_func = DerivedFuncWithName_PatchGroup;
   }

   FieldList[EoSIdx].field_unit             = "code_temperature";
   FieldList[EoSIdx].field_display_name     = "Temperature";

   FieldList[EoSIdx + 1].field_unit         = "code_mass / (code_length*code_time**2)";
   FieldList[EoSIdx + 1].field_display_name = "Pressure";

#  if   ( EOS == EOS_NUCLEAR )
   FieldList[EoSIdx + 2].field_unit         = "code_mass*code_length**(2) / (code_temperature*code_time**(2))";
#  elif ( EOS == EOS_GAMMA  ||  EOS == EOS_COSMIC_RAY )
   char EntropyUnit[100];
   real gamma_m1 = (real) GAMMA - 1.0;
   sprintf( EntropyUnit, "code_mass**(1-%.2f) / (code_length**(1-3*%.2f)*code_time**2)", gamma_m1, gamma_m1 );
   FieldList[EoSIdx + 2].field_unit         = EntropyUnit;
#  else
   Aux_Message( stderr, "WARNING : unknown entropy unit in %s() !!\n", __FUNCTION__ );
#  endif // #if ( EOS == EOS_GAMMA )
   FieldList[EoSIdx + 2].field_display_name = "Entropy";

#  endif // #if ( MODEL == HYDRO )

#  else  // #ifdef LIBYT_USE_PATCH_GROUP
   for (int v=0; v<NCOMP_TOTAL; v++)   FieldList[v].field_name = FieldLabel[v];

#  ifdef GRAVITY
   FieldList[PotIdx].field_name = const_cast<char*> (PotLabel);
#  endif

#  ifdef MHD
   char *CCMagLabel[] = {"CCMagX", "CCMagY", "CCMagZ"};
   for (int v=0; v<NCOMP_MAG; v++)
   {
       FieldList[v + MHDIdx].field_name  = CCMagLabel[v];
       FieldList[v + MHDIdx].field_type  = "face-centered";
       FieldList[v + MHDIdx].field_unit  = "code_magnetic";
   }

// add field display name
   FieldList[ MHDIdx     ].field_display_name = "B_x";
   FieldList[ MHDIdx + 1 ].field_display_name = "B_y";
   FieldList[ MHDIdx + 2 ].field_display_name = "B_z";

// add field derived function pointer
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

// set field's data type
   for (int v=0; v<NField; v++)
   {
#  ifdef FLOAT8
       FieldList[v].field_dtype = YT_DOUBLE;
#  else
       FieldList[v].field_dtype = YT_FLOAT;
#  endif
   }

// 3-2 get the ParticleList
#  ifdef PARTICLE
   yt_particle *ParticleList;
   yt_get_ParticlesPtr( &ParticleList );

// set floating-point attributes
   for (int v=0; v<PAR_NATT_FLT_TOTAL; v++)
   {
//     set attribute name
       ParticleList[0].attr_list[v].attr_name  = ParAttFltLabel[v];
//     set attribute data type
#      ifdef FLOAT8_PAR
       ParticleList[0].attr_list[v].attr_dtype = YT_DOUBLE;
#      else
       ParticleList[0].attr_list[v].attr_dtype = YT_FLOAT;
#      endif
//     set attribute unit
   }

// set integer attributes
   for (int v=PAR_NATT_FLT_TOTAL; v<ParticleList[0].num_attr; v++)
   {
//     set attribute name
       ParticleList[0].attr_list[v].attr_name  = ParAttIntLabel[v-PAR_NATT_FLT_TOTAL];
//     set attribute data type
#      ifdef INT8_PAR
       ParticleList[0].attr_list[v].attr_dtype = YT_LONG;
#      else
       ParticleList[0].attr_list[v].attr_dtype = YT_INT;
#      endif
//     set attribute unit
   }

// set label (attribute name) of coordinate x/y/z
   ParticleList[0].coor_x   = "ParPosX";
   ParticleList[0].coor_y   = "ParPosY";
   ParticleList[0].coor_z   = "ParPosZ";

// set get attribute function
   ParticleList[0].get_par_attr = Get_ParticleAttribute;

#  endif // #ifdef PARTICLE

// 4. prepare local patches for libyt
   YT_AddLocalGrid( NField, FieldList, pc);

// 5-1. perform yt inline analysis
   if ( yt_run_FunctionArguments( "yt_inline_inputArg", 1, "\'Dens\'" ) != YT_SUCCESS )
       Aux_Error( ERROR_INFO, "yt_run_FunctionArguments(\"yt_inline_inputArg\", 1, \"'Dens'\") failed !!\n" );
   if ( yt_run_Function( "yt_inline" ) != YT_SUCCESS )
       Aux_Error( ERROR_INFO, "yt_run_Function(\"yt_inline\") failed !!\n" );

#  ifdef LIBYT_INTERACTIVE
// 5-2. activate libyt interactive python prompt
   if ( yt_run_InteractiveMode( "LIBYT_STOP" ) != YT_SUCCESS )
       Aux_Error( ERROR_INFO, "yt_run_InteractiveMode(\"LIBYT_STOP\") failed !!\n" );
#  endif

#  ifdef LIBYT_RELOAD
// 5-3. activate libyt reloading script feature
   if ( yt_run_ReloadScript( "LIBYT_STOP", "RELOAD", "reload.py" ) != YT_SUCCESS )
       Aux_Error( ERROR_INFO, "yt_run_ReloadScript(\"LIBYT_STOP\", \"RELOAD\", \"reload.py\") failed !!\n" );
#  endif

#  ifdef LIBYT_JUPYTER
// 5-4. activate jupyter kernel
   if ( yt_run_JupyterKernel( "LIBYT_STOP", YT_JUPYTER_USE_CONNECTION_FILE ) != YT_SUCCESS )
       Aux_Error( ERROR_INFO, "yt_run_JupyterKernel(\"LIBYT_STOP\", %s) failed !!\n", YT_JUPYTER_USE_CONNECTION_FILE ? "true" : "false" );
#  endif

// 6. free resource
   if ( yt_free() != YT_SUCCESS )    Aux_Error( ERROR_INFO, "yt_free() failed !!\n" );

   if ( OPT__VERBOSE  &&  MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : YT_Inline

#endif // #ifdef SUPPORT_LIBYT
