#include "GAMER.h"

#ifdef SUPPORT_LIBYT

void YT_GetPID( const long gid, int *level, int *PID );



#ifdef LIBYT_USE_PATCH_GROUP
//-------------------------------------------------------------------------------------------------------
// Function    :  Fields_DerivedFuncWithName_PatchGroup
// Description :  Derived function for fields known in NCOMP_TOTAL.
//
// Note        :  1. This function's pointer will be passed into libyt.
//                2. The argument should be declared like this, in order to match the libyt API.
//                3. yt_getGridInfo_Dimensions() gets the grid_dimensions[0][1][2] in [x][y][z] coordinate.
//                4. libyt asks field through yt_field *FieldList you passed in, so we should convert
//                   to the one gamer understands. Try in this order :
//                   (1) GetFieldIndex: Will get fields defined in FieldLabel.
//                   (2) PotLabel     : Try PotLabel if supports GRAVITY.
//                   (3) CCMagLabel   : Try CCMagLabel if supports MHD.
//                   (3) Hydro        : Try find these fields if supports MODEL == HYDRO.
//
// Parameter   :  list_len    : length of list_gid
//                list_gid    : a list of grid id to prepare.
//                field       : target field.
//                data_array  : store data here, will be returned and wrapped by libyt.
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void DerivedFuncWithName_PatchGroup( const int list_len, const long *list_gid, const char *field, yt_array *data_array )
{

// get gamer field index through field name
   int gamer_fieldIdx = GetFieldIndex( field, CHECK_OFF );
   long gamer_fieldBIdx = -100;

// look through other options if cannot get FieldLabel defined, and convert to bitwise index
   if ( gamer_fieldIdx == Idx_Undefined )
   {
#     ifdef GRAVITY
      if ( strcmp( PotLabel, field ) == 0 )   gamer_fieldBIdx = _POTE;
#     endif

#     ifdef MHD
      const char *CCMagLabel[] = {"CCMagX", "CCMagY", "CCMagZ"};
      for (int v=0; v<NCOMP_MAG; v++)
      {
         if ( strcmp( CCMagLabel[v], field ) == 0 )
         {
            if      ( v == 0 )   gamer_fieldBIdx = _MAGX_CC;
            else if ( v == 1 )   gamer_fieldBIdx = _MAGY_CC;
            else if ( v == 2 )   gamer_fieldBIdx = _MAGZ_CC;
            break;
         }
      }
#     endif

#     if ( MODEL == HYDRO )
      if ( strcmp( "Temp", field ) == 0 )   gamer_fieldBIdx = _TEMP;
      if ( strcmp( "Pres", field ) == 0 )   gamer_fieldBIdx = _PRES;
      if ( strcmp( "Entr", field ) == 0 )   gamer_fieldBIdx = _ENTR;
#     endif
   }
   else
   {
      gamer_fieldBIdx = BIDX( gamer_fieldIdx );
   } // if ( gamer_fieldIdx == Idx_Undefined ) ... else ...

// if cannot find proper matching gamer bitwise index, raise an error
   if ( gamer_fieldBIdx == -100 )
      Aux_Error( ERROR_INFO, "cannot find the matching gamer field bitwise index for libyt field \"%s\" !!\n", field );

// loop through list_gid and fill in data
   for(int lid=0; lid<list_len; lid++)
   {
//    parse level and PID0
      int level, PID0;
      YT_GetPID( list_gid[lid], &level, &PID0 );

//    generate data in patch
      Prepare_PatchData( level, Time[0], (real*) data_array[lid].data_ptr, NULL, 0, 1, &PID0,
                         gamer_fieldBIdx, _NONE, INT_NONE, INT_NONE, UNIT_PATCHGROUP, NSIDE_00,
                         false, OPT__BC_FLU, BC_POT_NONE, -1.0, -1.0, -1.0, -1.0, false );
   }

} // FUNCITON : DerivedFuncWithName_PatchGroup



#else  // #ifdef LIBYT_USE_PATCH_GROUP
#ifdef MHD
//-------------------------------------------------------------------------------------------------------
// Function    :  MagX/Y/Z_DerivedFunc
// Description :  Derived function for CCMagX/Y/Z
//
// Note        :  1. This function's pointer will be passed into libyt.
//                2. The argument should be declared like this, in order to match the libyt API.
//                3. yt_getGridInfo_Dimensions() gets the grid_dimensions[0][1][2] in [x][y][z] coordinate.
//
// Parameter   :  list_len    : length of list_gid
//                list_gid    : a list of grid id to prepare.
//                data_array  : store data here, will be returned and wrapped by libyt.
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void MagX_DerivedFunc( const int list_len, const long *list_gid, const char *field_name, yt_array *data_array )
{

   for (int lid=0; lid<list_len; lid++)
   {
//    get the dimension of the grid, and the data array pointer of the grid
      int Dimensions[3];
      yt_data DataRaw;
      yt_getGridInfo_Dimensions( list_gid[lid], &Dimensions );
      yt_getGridInfo_FieldData( list_gid[lid], "CCMagX", &DataRaw );

//    cast the DataRaw
      real *Data = (real *) DataRaw.data_ptr;

//    compute converted field via data. [z, y, x] direction.
      for (int k=0; k<Dimensions[2]; k++) {
      for (int j=0; j<Dimensions[1]; j++) {
      for (int i=0; i<Dimensions[0]; i++) {
         int idx_Bx = IDX321_BX(i, j, k, Dimensions[0], Dimensions[1]);
         int idx_cc = i + j * Dimensions[0] + k * Dimensions[0] * Dimensions[1];
         ((real *) data_array[lid].data_ptr)[idx_cc] = 0.5 * ( Data[ idx_Bx ] + Data[ idx_Bx + 1 ] );
      }}} // i, j, k
   } // for(int lid=0; lid<list_len; lid++)

} // FUNCTION : MagX_DerivedFunc



void MagY_DerivedFunc( const int list_len, const long *list_gid, const char *field_name, yt_array *data_array )
{

   for (int lid=0; lid<list_len; lid++)
   {
//    get the dimension of the grid, and the data array pointer of the grid
      int Dimensions[3];
      yt_data DataRaw;
      yt_getGridInfo_Dimensions( list_gid[lid], &Dimensions );
      yt_getGridInfo_FieldData( list_gid[lid], "CCMagY", &DataRaw );

//    cast the DataRaw
      real *Data = (real *) DataRaw.data_ptr;

//    compute converted field via data. [z, y, x] direction.
      for (int k=0; k<Dimensions[2]; k++) {
      for (int j=0; j<Dimensions[1]; j++) {
      for (int i=0; i<Dimensions[0]; i++) {
         int idx_By = IDX321_BY(i, j, k, Dimensions[0], Dimensions[1]);
         int idx_cc = i + j * Dimensions[0] + k * Dimensions[0] * Dimensions[1];
         ((real *) data_array[lid].data_ptr)[idx_cc] = 0.5 * ( Data[ idx_By ] + Data[ idx_By + Dimensions[0] ] );
      }}} // i, j, k
    } // for (int lid=0; lid<list_len; lid++)

} // FUNCTION : MagY_DerivedFunc



void MagZ_DerivedFunc( const int list_len, const long *list_gid, const char *field_name, yt_array *data_array )
{

   for (int lid=0; lid<list_len; lid++)
   {
//    get the dimension of the grid, and the data array pointer of the grid
      int Dimensions[3];
      yt_data DataRaw;
      yt_getGridInfo_Dimensions( list_gid[lid], &Dimensions );
      yt_getGridInfo_FieldData( list_gid[lid], "CCMagZ", &DataRaw );

//    cast the DataRaw
      real *Data = (real *) DataRaw.data_ptr;

//    compute converted field via data. [z, y, x] direction.
      for (int k=0; k<Dimensions[2]; k++) {
      for (int j=0; j<Dimensions[1]; j++) {
      for (int i=0; i<Dimensions[0]; i++) {
         int idx_Bz = IDX321_BZ(i, j, k, Dimensions[0], Dimensions[1]);
         int idx_cc = i + j * Dimensions[0] + k * Dimensions[0] * Dimensions[1];
         ((real *) data_array[lid].data_ptr)[idx_cc] = 0.5 * ( Data[ idx_Bz ] + Data[ idx_Bz + Dimensions[0] * Dimensions[1] ] );
      }}} // i, j, k
   } // for (int lid=0; lid<list_len; lid++)

} // FUNCTION : MagZ_DerivedFunc
#endif // #ifdef MHD



#if ( MODEL == HYDRO )
//-------------------------------------------------------------------------------------------------------
// Function    :  Temperature_DerivedFunc
// Description :  Derived function for fluid temperature invoke by EoS routine
//
// Note        :  1. This function's pointer will be passed into libyt.
//                2. The argument should be declared like this, in order to match the libyt API.
//                3. yt_getGridInfo_Dimensions() gets the grid_dimensions[0][1][2] in [x][y][z] coordinate.
//
// Parameter   :  list_len    : length of list_gid
//                list_gid    : a list of grid id to prepare.
//                data_array  : store data here, will be returned and wrapped by libyt.
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Temperature_DerivedFunc( const int list_len, const long *list_gid, const char *field_name, yt_array *data_array )
{

// universal
   char *CCMagLabel[] = {"CCMagX", "CCMagY", "CCMagZ"};

// loop through list_gid
   for (int lid=0; lid<list_len; lid++)
   {
//    get dim of the grid to be return, and all the other NCOMP_TOTAL fields.
      int Dimensions[3];
      yt_data   DataRaw[NCOMP_TOTAL];
      real     *Data[NCOMP_TOTAL];
      yt_getGridInfo_Dimensions( list_gid[lid], &Dimensions );
      for (int v=0; v<NCOMP_TOTAL; v++)
      {
         yt_getGridInfo_FieldData( list_gid[lid], FieldLabel[v], &(DataRaw[v]));
         Data[v] = (real *) DataRaw[v].data_ptr;
      }

//    preparation for getting Passive and Emag
      real Passive[NCOMP_PASSIVE];
      real Emag;
#     ifdef MHD
      yt_data   MagDataRaw[NCOMP_MAG];
      real     *MagData[NCOMP_MAG];
      for (int v=0; v<NCOMP_MAG; v++)
      {
         yt_getGridInfo_FieldData( list_gid[lid], CCMagLabel[v], &(MagDataRaw[v]));
         MagData[v] = (real *) MagDataRaw[v].data_ptr;
      }
#     endif

//    get temperature cell-by-cell
      for (int k=0; k<Dimensions[2]; k++) {
      for (int j=0; j<Dimensions[1]; j++) {
      for (int i=0; i<Dimensions[0]; i++) {
         int idx = i + j * Dimensions[0] + k * Dimensions[0] * Dimensions[1];

//       get Passive
         for (int v=NCOMP_FLUID; v<NCOMP_TOTAL; v++)
         {
            Passive[v-NCOMP_FLUID] = Data[v][idx];
         }

//       get Emag
#        ifdef MHD
         Emag = MHD_GetCellCenteredBEnergy( MagData[MAGX], MagData[MAGY], MagData[MAGZ],
                                            Dimensions[0], Dimensions[1], Dimensions[2], i, j, k );
#        else
         Emag = NULL_REAL;
#        endif

//       get temperature
         ((real *) data_array[lid].data_ptr)[idx] = Hydro_Con2Temp( Data[DENS][idx], Data[MOMX][idx], Data[MOMY][idx], Data[MOMZ][idx], Data[ENGY][idx],
                                                                    Passive, false, NULL_REAL, Emag,
                                                                    EoS_DensEint2Temp_CPUPtr, EoS_AuxArray_Flt, EoS_AuxArray_Int, h_EoS_Table );
      }}} // i, j, k
    } // for (int lid=0; lid<list_len; lid++)

} // FUNCITON : Temperature_DerivedFunc
#endif // #if ( MODEL == HYDRO )

#endif // #ifdef LIBYT_USE_PATCH_GROUP

#endif // #ifdef SUPPORT_LIBYT
