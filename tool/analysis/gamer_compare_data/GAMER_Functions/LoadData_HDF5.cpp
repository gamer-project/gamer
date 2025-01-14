#ifdef SUPPORT_HDF5

#include "CompareData.h"
#include "hdf5.h"
#include <typeinfo>

#ifdef FLOAT8
#  define H5T_GAMER_REAL H5T_NATIVE_DOUBLE
#else
#  define H5T_GAMER_REAL H5T_NATIVE_FLOAT
#endif

#ifdef FLOAT8_PAR
#  define H5T_GAMER_REAL_PAR H5T_NATIVE_DOUBLE
#else
#  define H5T_GAMER_REAL_PAR H5T_NATIVE_FLOAT
#endif

#ifdef INT8_PAR
#  define H5T_GAMER_LONG_PAR H5T_NATIVE_LONG
#else
#  define H5T_GAMER_LONG_PAR H5T_NATIVE_INT
#endif

#ifdef GAMER_DEBUG
#  define DEBUG_HDF5
#endif

template <typename T>
static herr_t LoadField( const char *FieldLabel, void *FieldPtr, const hid_t H5_SetID_Target,
                         const hid_t H5_TypeID_Target, const bool Fatal_Nonexist,
                         const T *ComprPtr, const int NCompr, const bool Fatal_Compr );
static void LoadOnePatch( AMR_t &amr, const hid_t H5_FileID, const int lv, const int GID,
                          const bool Recursive, const int *SonList, const int (*CrList)[3],
                          const hid_t *H5_SetID_Field, const hid_t  H5_SpaceID_Field, const hid_t  H5_MemID_Field, const int NField,
                          const hid_t *H5_SetID_Mag,   const hid_t *H5_SpaceID_Mag,   const hid_t *H5_MemID_Mag,   const int NMag );




//-------------------------------------------------------------------------------------------------------
// Function    :  LoadData_HDF5
// Description :  Load data from the input HDF5 file
//
// Note        :  1. Only work for format version >= 2440
//                   --> NField and NMag are supported only after 2440
//                2. Unlike LoadData(), this function always loads and allocates both leaf and non-leaf patches
//                3. *Label[MAX_STRING] and ParData[] will be allocated here and must be deallocated manually
//
// Parameter   :  FileName       : Name of the input file
//                amr            : Target AMR_t pointer
//                Format         : 1/2 --> C-binary/HDF5
//                NField         : Number of cell-centered fields stored in the file
//                NMag           : Number of face-centered magnetic components stored in the file
//                NParAttFlt     : Number of particle floating-point attributes stored in the file
//                NParAttInt     : Number of particle integer        attributes stored in the file
//                NPar           : NUmber of particles
//                ParFltData     : Particle floating-point data array (allocated here --> must be dallocated manually later)
//                ParIntData     : Particle integer        data array (allocated here --> must be dallocated manually later)
//                FieldLabel     : Labels of cell-centered fields
//                MagLabel       : Labels of face-centered magnetic components
//                ParAttFltLabel : Labels of particle floating-point attributes
//                ParAttIntLabel : Labels of particle integer        attributes
//
// Return      :  amr, Format, NField, NMag, NParAttFlt, NParAttInt, NPar, ParData
//                FieldLabel, MagLabel, ParAttFltLabel, ParAttIntLabel
//-------------------------------------------------------------------------------------------------------
void LoadData_HDF5( const char *FileName, AMR_t &amr, int &Format, int &NField, int &NMag, int &NParAttFlt,
                    int &NParAttInt, long &NPar, real_par **&ParFltData, long_par **&ParIntData, char (*&FieldLabel)[MAX_STRING],
                    char (*&MagLabel)[MAX_STRING], char (*&ParAttFltLabel)[MAX_STRING], char (*&ParAttIntLabel)[MAX_STRING] )
{

   Aux_Message( stdout, "Loading HDF5 data %s ...\n", FileName );


// check
   if ( !Aux_CheckFileExist(FileName) )
      Aux_Error( ERROR_INFO, "input file \"%s\" does not exist !!\n", FileName );

   if ( !H5Fis_hdf5(FileName) )
      Aux_Error( ERROR_INFO, "input file \"%s\" is not in the HDF5 format !!\n", FileName );


// 1. load the simulation info
   Aux_Message( stdout, "   Loading simulation information ...\n" );

   const bool    Fatal     = true;
   const bool NonFatal     = false;
   const int  PatchSize_RT = PS1;
   const int  NLevel_RT    = NLEVEL;
#  ifdef FLOAT8
   const int  Float8_RT    = 1;
#  else
   const int  Float8_RT    = 0;
#  endif

   int    PatchSize_RS, NLevel_RS, Float8_RS;
   int    FormatVersion, NPatchTotal[NLEVEL], NPatchAllLv, DumpID;
   long   Step;
   double Time[NLEVEL];
   char **FieldLabel_In=NULL, **MagLabel_In=NULL, **ParAttFltLabel_In=NULL, **ParAttIntLabel_In=NULL;
   int   *NullPtr=NULL;
   int    WithPar;

   hid_t  H5_FileID, H5_SetID_KeyInfo, H5_TypeID_KeyInfo, H5_SetID_InputPara, H5_TypeID_InputPara;
   hid_t  H5_SetID_Cr, H5_SetID_Son;
   herr_t H5_Status;


// 1-1. open the HDF5 file
   H5_FileID = H5Fopen( FileName, H5F_ACC_RDONLY, H5P_DEFAULT );
   if ( H5_FileID < 0 )
      Aux_Error( ERROR_INFO, "failed to open the input HDF5 file \"%s\" !!\n", FileName );


// 1-2. load the dataset and datatype of KeyInfo and InputPara
   H5_SetID_KeyInfo  = H5Dopen( H5_FileID, "Info/KeyInfo", H5P_DEFAULT );
   if ( H5_SetID_KeyInfo < 0 )
      Aux_Error( ERROR_INFO, "failed to open the dataset \"%s\" !!\n", "Info/KeyInfo" );

   H5_TypeID_KeyInfo = H5Dget_type( H5_SetID_KeyInfo );
   if ( H5_TypeID_KeyInfo < 0 )
      Aux_Error( ERROR_INFO, "failed to open the datatype of \"%s\" !!\n", "Info/KeyInfo" );

   H5_SetID_InputPara  = H5Dopen( H5_FileID, "Info/InputPara", H5P_DEFAULT );
   if ( H5_SetID_InputPara < 0 )
      Aux_Error( ERROR_INFO, "failed to open the dataset \"%s\" !!\n", "Info/InputPara" );

   H5_TypeID_InputPara = H5Dget_type( H5_SetID_InputPara );
   if ( H5_TypeID_InputPara < 0 )
      Aux_Error( ERROR_INFO, "failed to open the datatype of \"%s\" !!\n", "Info/InputPara" );


// 1-3. load relevant target fields in KeyInfo and InputPara
// load and check the format version of HDF5 output
   LoadField( "FormatVersion",        &FormatVersion,        H5_SetID_KeyInfo,   H5_TypeID_KeyInfo,   Fatal,  NullPtr,        -1, NonFatal );

   Aux_Message( stdout, "      The format version of the HDF5 file = %ld\n", FormatVersion );

   if ( FormatVersion < 2440 )
      Aux_Error( ERROR_INFO, "unsupported data format version (only support version >= 2440) !!\n" );

// general info
   LoadField( "Float8",               &Float8_RS,           H5_SetID_KeyInfo,    H5_TypeID_KeyInfo,   Fatal,  &Float8_RT,        1,    Fatal );
   LoadField( "NLevel",               &NLevel_RS,           H5_SetID_KeyInfo,    H5_TypeID_KeyInfo,   Fatal,  &NLevel_RT,        1,    Fatal );
   LoadField( "PatchSize",            &PatchSize_RS,        H5_SetID_KeyInfo,    H5_TypeID_KeyInfo,   Fatal,  &PatchSize_RT,     1,    Fatal );
   LoadField( "NPatch",                NPatchTotal,         H5_SetID_KeyInfo,    H5_TypeID_KeyInfo,   Fatal,   NullPtr,         -1, NonFatal );
   LoadField( "NX0",                   amr.nx0_tot,         H5_SetID_KeyInfo,    H5_TypeID_KeyInfo,   Fatal,   NullPtr,         -1, NonFatal );
   LoadField( "DumpID",               &DumpID,              H5_SetID_KeyInfo,    H5_TypeID_KeyInfo,   Fatal,   NullPtr,         -1, NonFatal );
   LoadField( "Step",                 &Step,                H5_SetID_KeyInfo,    H5_TypeID_KeyInfo,   Fatal,   NullPtr,         -1, NonFatal );
   LoadField( "Time",                  Time,                H5_SetID_KeyInfo,    H5_TypeID_KeyInfo,   Fatal,   NullPtr,         -1, NonFatal );

// numbers of grid fields, particle attributes, and particles
   LoadField( "NFieldStored",         &NField,              H5_SetID_KeyInfo,    H5_TypeID_KeyInfo,   Fatal,   NullPtr,         -1, NonFatal );
   LoadField( "NMagStored",           &NMag,                H5_SetID_KeyInfo,    H5_TypeID_KeyInfo,   Fatal,   NullPtr,         -1, NonFatal );
   LoadField( "Particle",             &WithPar,             H5_SetID_KeyInfo,    H5_TypeID_KeyInfo,   Fatal,   NullPtr,         -1, NonFatal );
   if ( WithPar ) {
   int Float8_Par_RS;
   int Int8_Par_RS;
#  ifdef FLOAT8_PAR
   const int  Float8_Par_RT    = 1;
#  else
   const int  Float8_Par_RT    = 0;
#  endif
#  ifdef INT8_PAR
   const int  Int8_Par_RT      = 1;
#  else
   const int  Int8_Par_RT      = 0;
#  endif
   int Float8_Par_check_flag;
   LoadField( "Par_NPar",             &NPar,                H5_SetID_KeyInfo,    H5_TypeID_KeyInfo,    Fatal,   NullPtr,         -1, NonFatal );
   Float8_Par_check_flag =
   LoadField( "Float8_Par",           &Float8_Par_RS,       H5_SetID_KeyInfo,    H5_TypeID_KeyInfo, NonFatal,  &Float8_Par_RT,    1,    Fatal );
   if ( FormatVersion >= 2500 )
   {
   LoadField( "Par_NAttFltStored",    &NParAttFlt,          H5_SetID_KeyInfo,    H5_TypeID_KeyInfo,    Fatal,   NullPtr,         -1, NonFatal );
   LoadField( "Par_NAttIntStored",    &NParAttInt,          H5_SetID_KeyInfo,    H5_TypeID_KeyInfo,    Fatal,   NullPtr,         -1, NonFatal );
   LoadField( "Int8_Par",             &Int8_Par_RS,         H5_SetID_KeyInfo,    H5_TypeID_KeyInfo,    Fatal,  &Int8_Par_RT,      1,    Fatal );
   }
   else
   {
   LoadField( "Par_NAttStored",       &NParAttFlt,          H5_SetID_KeyInfo,    H5_TypeID_KeyInfo,    Fatal,   NullPtr,         -1, NonFatal );
   NParAttInt = 0;
   LoadField( "Int8_Par",             &Int8_Par_RS,         H5_SetID_KeyInfo,    H5_TypeID_KeyInfo, NonFatal,  &Int8_Par_RT,      1,    Fatal );
   }
   if ( Float8_Par_check_flag != 0  &&  sizeof(real) != sizeof(real_par) )
      Aux_Error( ERROR_INFO, "Must adopt FLOAT8_PAR=FLOAT8 in Makefile when Float8_Par is not stored in the snapshot !!\n");
   } // if ( WithPar )
   else {
      NParAttFlt = 0;
      NParAttInt = 0;
      NPar       = 0;
   } // if ( WithPar ) ... else ...

// field and particle attribute labels
   FieldLabel_In = new char* [NField];
   for (int v=0; v<NField; v++)
   {
      char Key[MAX_STRING];
      sprintf( Key, "FieldLabel%02d", v );
      LoadField( Key,                &FieldLabel_In[v],     H5_SetID_InputPara,  H5_TypeID_InputPara, Fatal,   NullPtr,         -1, NonFatal );
   }

   MagLabel_In = new char* [NMag];
   for (int v=0; v<NMag; v++)
   {
      char Key[MAX_STRING];
      sprintf( Key, "MagLabel%02d", v );
      LoadField( Key,                &MagLabel_In[v],       H5_SetID_InputPara,  H5_TypeID_InputPara, Fatal,   NullPtr,         -1, NonFatal );
   }

   ParAttFltLabel_In = new char* [NParAttFlt];
   for (int v=0; v<NParAttFlt; v++)
   {
      char Key[MAX_STRING];
      if ( FormatVersion >= 2500 )   sprintf( Key, "ParAttFltLabel%02d", v );
      else                           sprintf( Key,    "ParAttLabel%02d", v );
      LoadField( Key,                &ParAttFltLabel_In[v], H5_SetID_InputPara,  H5_TypeID_InputPara, Fatal,   NullPtr,         -1, NonFatal );
   }

   ParAttIntLabel_In = new char* [NParAttInt];
   for (int v=0; v<NParAttInt; v++)
   {
      char Key[MAX_STRING];
      sprintf( Key, "ParAttIntLabel%02d", v );
      LoadField( Key,                &ParAttIntLabel_In[v], H5_SetID_InputPara,  H5_TypeID_InputPara, Fatal,   NullPtr,         -1, NonFatal );
   }

// total number of patches
   NPatchAllLv = 0;
   for (int lv=0; lv<NLEVEL; lv++)  NPatchAllLv += NPatchTotal[lv];


// 1-4. close all objects
   H5_Status = H5Tclose( H5_TypeID_KeyInfo );
   H5_Status = H5Dclose( H5_SetID_KeyInfo );
   H5_Status = H5Tclose( H5_TypeID_InputPara );
   H5_Status = H5Dclose( H5_SetID_InputPara );


   Aux_Message( stdout, "   Loading simulation information ... done\n" );



// 2. load the tree information (i.e., corner and son) of all patches
// 2-1. corner
   Aux_Message( stdout, "   Loading corner table ...\n" );

// allocate memory
   int (*CrList_AllLv)[3] = new int [ NPatchAllLv ][3];

// load data
   H5_SetID_Cr = H5Dopen( H5_FileID, "Tree/Corner", H5P_DEFAULT );
   if ( H5_SetID_Cr < 0 )     Aux_Error( ERROR_INFO, "failed to open the dataset \"%s\" !!\n", "Tree/Corner" );
   H5_Status = H5Dread( H5_SetID_Cr, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, CrList_AllLv );
   H5_Status = H5Dclose( H5_SetID_Cr );

   Aux_Message( stdout, "   Loading corner table ... done\n" );


// 2-2. son
   Aux_Message( stdout, "   Loading son table ...\n" );

// allocate memory
   int *SonList_AllLv = new int [ NPatchAllLv ];

// load data
   H5_SetID_Son = H5Dopen( H5_FileID, "Tree/Son", H5P_DEFAULT );
   if ( H5_SetID_Son < 0 )    Aux_Error( ERROR_INFO, "failed to open the dataset \"%s\" !!\n", "Tree/Son" );
   H5_Status = H5Dread( H5_SetID_Son, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, SonList_AllLv );
   H5_Status = H5Dclose( H5_SetID_Son );

   Aux_Message( stdout, "   Loading son table ... done\n" );



// 3. load and allocate patches
   Aux_Message( stdout, "   Loading patches ...\n" );

   const bool Recursive_Yes = true;

   FieldLabel = new char [NField][MAX_STRING];
   MagLabel   = new char [NMag  ][MAX_STRING];

   hsize_t H5_SetDims_Field[4], H5_MemDims_Field[4];
   hid_t   H5_SetID_Field[NField], H5_MemID_Field, H5_SpaceID_Field, H5_GroupID_GridData;
   hsize_t H5_SetDims_Mag[4], H5_MemDims_Mag[4];
   hid_t   H5_SetID_Mag[NMag], H5_MemID_Mag[NMag], H5_SpaceID_Mag[NMag];


// 3-1. set the names of all grid variables
   for (int v=0; v<NField; v++)  sprintf( FieldLabel[v], "%s", FieldLabel_In[v] );
   for (int v=0; v<NMag;   v++)  sprintf(   MagLabel[v], "%s",   MagLabel_In[v] );


// 3-2. initialize relevant HDF5 objects
   H5_SetDims_Field[0] = NPatchAllLv;
   H5_SetDims_Field[1] = PS1;
   H5_SetDims_Field[2] = PS1;
   H5_SetDims_Field[3] = PS1;

   H5_SpaceID_Field = H5Screate_simple( 4, H5_SetDims_Field, NULL );
   if ( H5_SpaceID_Field < 0 )   Aux_Error( ERROR_INFO, "failed to create the space \"%s\" !!\n", "H5_SpaceID_Field" );

   H5_MemDims_Field[0] = 1;
   H5_MemDims_Field[1] = PS1;
   H5_MemDims_Field[2] = PS1;
   H5_MemDims_Field[3] = PS1;

   H5_MemID_Field = H5Screate_simple( 4, H5_MemDims_Field, NULL );
   if ( H5_MemID_Field < 0 )  Aux_Error( ERROR_INFO, "failed to create the space \"%s\" !!\n", "H5_MemDims_Field" );

   for (int v=0; v<NMag; v++)
   {
      H5_SetDims_Mag[0] = NPatchAllLv;
      for (int t=1; t<4; t++)
      H5_SetDims_Mag[t] = ( 3-t == v ) ? PS1P1 : PS1;

      H5_SpaceID_Mag[v] = H5Screate_simple( 4, H5_SetDims_Mag, NULL );
      if ( H5_SpaceID_Mag[v] < 0 )
         Aux_Error( ERROR_INFO, "failed to create the space \"%s[%d]\" !!\n", "H5_SpaceID_Mag", v );

      H5_MemDims_Mag[0] = 1;
      for (int t=1; t<4; t++)
      H5_MemDims_Mag[t] = ( 3-t == v ) ? PS1P1 : PS1;

      H5_MemID_Mag[v] = H5Screate_simple( 4, H5_MemDims_Mag, NULL );
      if ( H5_MemID_Mag[v] < 0 )
         Aux_Error( ERROR_INFO, "failed to create the space \"%s[%d]\" !!\n", "H5_MemID_Mag", v );
   }


// 3-3. open the target datasets just once
   H5_GroupID_GridData = H5Gopen( H5_FileID, "GridData", H5P_DEFAULT );
   if ( H5_GroupID_GridData < 0 )   Aux_Error( ERROR_INFO, "failed to open the group \"%s\" !!\n", "GridData" );

   for (int v=0; v<NField; v++)
   {
      H5_SetID_Field[v] = H5Dopen( H5_GroupID_GridData, FieldLabel[v], H5P_DEFAULT );
      if ( H5_SetID_Field[v] < 0 )  Aux_Error( ERROR_INFO, "failed to open the dataset \"%s\" !!\n", FieldLabel[v] );
   }

   for (int v=0; v<NMag; v++)
   {
      H5_SetID_Mag  [v] = H5Dopen( H5_GroupID_GridData, MagLabel[v],   H5P_DEFAULT );
      if ( H5_SetID_Mag[v] < 0 )    Aux_Error( ERROR_INFO, "failed to open the dataset \"%s\" !!\n", MagLabel[v] );
   }


// 3-4. begin to load data
// loop over all root-level patches
   const int TenPercent = MAX( NPatchTotal[0]/10, 1 );

   for (int GID=0; GID<NPatchTotal[0]; GID++)
   {
      if ( GID % TenPercent == 0 )
      Aux_Message( stdout, "      %5.1lf%% completed ...\n", 100.0*GID/NPatchTotal[0] );

//    load the entire patch family recursively (actually it's not necessary here since we load all patches anyway)
      LoadOnePatch( amr, H5_FileID, 0, GID, Recursive_Yes, SonList_AllLv, CrList_AllLv,
                    H5_SetID_Field, H5_SpaceID_Field, H5_MemID_Field, NField,
                    H5_SetID_Mag,   H5_SpaceID_Mag,   H5_MemID_Mag,   NMag );
   }

// free HDF5 objects
   for (int v=0; v<NField; v++)  H5_Status = H5Dclose( H5_SetID_Field[v] );
   for (int v=0; v<NMag;   v++)  H5_Status = H5Dclose( H5_SetID_Mag  [v] );

   H5_Status = H5Gclose( H5_GroupID_GridData );
   H5_Status = H5Sclose( H5_SpaceID_Field );
   H5_Status = H5Sclose( H5_MemID_Field );
   for (int v=0; v<NMag; v++)
   {
      H5_Status = H5Sclose( H5_SpaceID_Mag[v] );
      H5_Status = H5Sclose( H5_MemID_Mag  [v] );
   }


// 3-5. record the total number of loaded patches at each level
   Aux_Message( stdout, "\n" );
   for (int lv=0; lv<NLEVEL; lv++)     Aux_Message( stdout, "      Lv %2d: %10d patches loaded\n", lv, amr.num[lv] );
   Aux_Message( stdout, "\n" );

   for (int lv=0; lv<NLEVEL; lv++)
   {
      if ( amr.num[lv] != NPatchTotal[lv] )
         Aux_Error( ERROR_INFO, "Lv %2d: amr.num (%d) != NPatchTotal (%d) !!\n",
                    lv, amr.num[lv], NPatchTotal[lv] );
   }

   Aux_Message( stdout, "   Loading patches ... done\n" );



// 4. load particles
// 4-1. set the names of all particle attributes
   ParAttFltLabel = new char [NParAttFlt][MAX_STRING];
   ParAttIntLabel = new char [NParAttInt][MAX_STRING];
   for (int v=0; v<NParAttFlt; v++)    sprintf( ParAttFltLabel[v], "%s", ParAttFltLabel_In[v] );
   for (int v=0; v<NParAttInt; v++)    sprintf( ParAttIntLabel[v], "%s", ParAttIntLabel_In[v] );

   if ( NPar > 0 )
   {
//    4-2. allocate the particle data array
      Aux_AllocateArray2D<real_par>( ParFltData, NParAttFlt, NPar );
      Aux_AllocateArray2D<long_par>( ParIntData, NParAttInt, NPar );


//    4-3. initialize HDF5 objects
      hsize_t H5_SetDims_ParData[1], H5_MemDims_ParData[1];
      hid_t   H5_GroupID_ParData, H5_SpaceID_ParData, H5_MemID_ParData, H5_SetID_ParFltData[NParAttFlt], H5_SetID_ParIntData[NParAttInt];

      H5_GroupID_ParData = H5Gopen( H5_FileID, "Particle", H5P_DEFAULT );
      if ( H5_GroupID_ParData < 0 )    Aux_Error( ERROR_INFO, "failed to open the group \"%s\" !!\n", "Particle" );

      H5_SetDims_ParData[0] = NPar;
      H5_SpaceID_ParData    = H5Screate_simple( 1, H5_SetDims_ParData, NULL );
      if ( H5_SpaceID_ParData < 0 )    Aux_Error( ERROR_INFO, "failed to create the space \"%s\" !!\n", "H5_SpaceID_ParData" );

      H5_MemDims_ParData[0] = NPar;
      H5_MemID_ParData      = H5Screate_simple( 1, H5_MemDims_ParData, NULL );
      if ( H5_MemID_ParData < 0 )      Aux_Error( ERROR_INFO, "failed to create the space \"%s\" !!\n", "H5_MemID_ParData" );

      for (int v=0; v<NParAttFlt; v++)
      {
         H5_SetID_ParFltData[v] = H5Dopen( H5_GroupID_ParData, ParAttFltLabel[v], H5P_DEFAULT );
         if ( H5_SetID_ParFltData[v] < 0 )   Aux_Error( ERROR_INFO, "failed to open the dataset \"%s\" !!\n", ParAttFltLabel[v] );
      }

      for (int v=0; v<NParAttInt; v++)
      {
         H5_SetID_ParIntData[v] = H5Dopen( H5_GroupID_ParData, ParAttIntLabel[v], H5P_DEFAULT );
         if ( H5_SetID_ParIntData[v] < 0 )   Aux_Error( ERROR_INFO, "failed to open the dataset \"%s\" !!\n", ParAttIntLabel[v] );
      }


//    4-4. load particle data
      for (int v=0; v<NParAttFlt; v++)
      {
         H5_Status = H5Dread( H5_SetID_ParFltData[v], H5T_GAMER_REAL_PAR, H5_MemID_ParData, H5_SpaceID_ParData, H5P_DEFAULT, ParFltData[v] );
         if ( H5_Status < 0 )    Aux_Error( ERROR_INFO, "failed to load the particle floating-point attribute %d !!\n", v );
      }

      for (int v=0; v<NParAttInt; v++)
      {
         H5_Status = H5Dread( H5_SetID_ParIntData[v], H5T_GAMER_LONG_PAR, H5_MemID_ParData, H5_SpaceID_ParData, H5P_DEFAULT, ParIntData[v] );
         if ( H5_Status < 0 )    Aux_Error( ERROR_INFO, "failed to load the particle integer attribute %d !!\n", v );
      }


//    4-5. close HDF5 objects
      for (int v=0; v<NParAttFlt; v++)    H5_Status = H5Dclose( H5_SetID_ParFltData[v] );
      for (int v=0; v<NParAttInt; v++)    H5_Status = H5Dclose( H5_SetID_ParIntData[v] );
      H5_Status = H5Sclose( H5_MemID_ParData );
      H5_Status = H5Sclose( H5_SpaceID_ParData );
      H5_Status = H5Gclose( H5_GroupID_ParData );
   } // if ( NPar > 0 )



// 5. close all HDF5 objects and free memory
   for (int v=0; v<NField;     v++)    free(     FieldLabel_In[v] );
   for (int v=0; v<NMag;       v++)    free(       MagLabel_In[v] );
   for (int v=0; v<NParAttFlt; v++)    free( ParAttFltLabel_In[v] );
   for (int v=0; v<NParAttInt; v++)    free( ParAttIntLabel_In[v] );
   delete []     FieldLabel_In;
   delete []       MagLabel_In;
   delete [] ParAttFltLabel_In;
   delete [] ParAttIntLabel_In;
   delete [] CrList_AllLv;
   delete [] SonList_AllLv;

   H5_Status = H5Fclose( H5_FileID );



// 6. record parameters
   Format = 2;

   Aux_Message( stdout, "\n" );
   Aux_Message( stdout, "   Data information\n" );
   Aux_Message( stdout, "   ==============================\n" );
   Aux_Message( stdout, "   DumpID      = %d\n",  DumpID      );
   Aux_Message( stdout, "   Step        = %ld\n", Step        );
   Aux_Message( stdout, "   Time        = %lf\n", Time[0]     );
   Aux_Message( stdout, "   Format      = %d\n",  Format      );
   Aux_Message( stdout, "   WithPar     = %d\n",  WithPar     );
   Aux_Message( stdout, "   NField      = %d\n",  NField      );
   Aux_Message( stdout, "   NMag        = %d\n",  NMag        );
   if ( WithPar ) {
   Aux_Message( stdout, "   NParAttFlt  = %d\n",  NParAttFlt  );
   Aux_Message( stdout, "   NParAttInt  = %d\n",  NParAttInt  );
   Aux_Message( stdout, "   NPar        = %ld\n", NPar        ); }
   for (int lv=0; lv<NLEVEL; lv++)
   Aux_Message( stdout, "   NPatch[%2d] = %d\n",  lv, amr.num[lv] );
   Aux_Message( stdout, "   ==============================\n" );
   Aux_Message( stdout, "\n" );


   Aux_Message( stdout, "Loading HDF5 data %s ... done\n", FileName );

} // FUNCTION : LoadData_HDF5



//-------------------------------------------------------------------------------------------------------
// Function    :  LoadField
// Description :  Load a single field from the input compound dataset
//
// Note        :  1. This function works for arbitary datatype (int, float, char, 1D array ...)
//                2. Memory must be allocated for FieldPtr in advance with sufficent size (except for "char *")
//                3. For loading a string, which has (type(FieldPtr) = (char *)), the memory must be freed
//                   manually by calling "free()"
//                4. It can also compare the loaded variables (FieldPtr) with the reference value (ComprPtr)
//                   (perform comparison only if "NCompr > 0")
//                   --> Please make sure that "FieldPtr" and "ComprPtr" point to the same type since we
//                       use the type of "ComprPtr" to typecast "FieldPtr"
//
// Parameter   :  FieldLabel        : Name of the target field
//                FieldPtr          : Pointer to store the retrieved data
//                H5_SetID_Target   : HDF5 dataset  ID of the target compound variable
//                H5_TypeID_Target  : HDF5 datatype ID of the target compound variable
//                Fatal_Nonexist    : Whether or not the nonexistence of the target field is fatal
//                                    --> true  : terminate the program     if the target field cannot be found
//                                        false : display a warning message if the target field cannot be found
//                ComprPtr          : Pointer to store the reference values for comparison
//                NCompr            : Number of elements to be compared
//                Fatal_Compr       : Whether or not the comparison result is fatal
//                                    --> true  : terminate the program     if "FieldPtr[X] != ComprPtr[X]"
//                                        false : display a warning message if "FieldPtr[X] != ComprPtr[X]"
//
// Return      :  Success/fail <-> 0/-1
//-------------------------------------------------------------------------------------------------------
template <typename T>
herr_t LoadField( const char *FieldLabel, void *FieldPtr, const hid_t H5_SetID_Target,
                  const hid_t H5_TypeID_Target, const bool Fatal_Nonexist,
                  const T *ComprPtr, const int NCompr, const bool Fatal_Compr )
{

#  ifdef DEBUG_HDF5
   if ( NCompr > 0  &&  ComprPtr == NULL )
      Aux_Error( ERROR_INFO, "ComprPtr == NULL for NCompr = %d > 0 !!\n", NCompr );
#  endif


   bool   Check_Pass = true;
   int    H5_FieldIdx;
   size_t H5_FieldSize;
   hid_t  H5_TypeID_Field;    // datatype ID of the target field in the compound variable
   hid_t  H5_TypeID_Load;     // datatype ID for loading the target field
   herr_t H5_Status;


// load
   H5_FieldIdx = H5Tget_member_index( H5_TypeID_Target, FieldLabel );

   if ( H5_FieldIdx >= 0 )
   {
      H5_TypeID_Field  = H5Tget_member_type( H5_TypeID_Target, H5_FieldIdx );
      H5_FieldSize     = H5Tget_size( H5_TypeID_Field );

      H5_TypeID_Load   = H5Tcreate( H5T_COMPOUND, H5_FieldSize );
      H5_Status        = H5Tinsert( H5_TypeID_Load, FieldLabel, 0, H5_TypeID_Field );

      H5_Status        = H5Dread( H5_SetID_Target, H5_TypeID_Load, H5S_ALL, H5S_ALL, H5P_DEFAULT, FieldPtr );
      if ( H5_Status < 0 )    Aux_Error( ERROR_INFO, "failed to load the field \"%s\" !!\n", FieldLabel );

      H5_Status        = H5Tclose( H5_TypeID_Field );
      H5_Status        = H5Tclose( H5_TypeID_Load  );
   } // if ( H5_FieldIdx >= 0 )

   else
   {
      if ( Fatal_Nonexist )
         Aux_Error( ERROR_INFO, "target field \"%s\" does not exist in the input file !!\n", FieldLabel );

      else
         Aux_Message( stderr, "WARNING : target field \"%s\" does not exist in the input file !!\n", FieldLabel );

      return -1;
   } // if ( H5_FieldIdx >= 0 ) ... else ...


// comparison
   char ArrayIdx[10];

   for (int t=0; t<NCompr; t++)
   {
      if ( NCompr == 1 )   sprintf( ArrayIdx, "%s", "" );
      else                 sprintf( ArrayIdx, "[%d]", t );

      if ( ((T*)FieldPtr)[t] != ComprPtr[t] )
      {
         if (  typeid(T) == typeid( int)  ||  typeid(T) == typeid( long)  ||
               typeid(T) == typeid(uint)  ||  typeid(T) == typeid(ulong)    )
         {
            if ( Fatal_Compr )
            {
               Aux_Error( ERROR_INFO, "\"%s%s\" : input file (%ld) != runtime (%ld) !!\n",
                          FieldLabel, ArrayIdx, (long)((T*)FieldPtr)[t], (long)ComprPtr[t] );
               return -2;
            }

            else
            {
               Aux_Message( stderr, "WARNING : \"%s%s\" : input file (%ld) != runtime (%ld) !!\n",
                            FieldLabel, ArrayIdx, (long)((T*)FieldPtr)[t], (long)ComprPtr[t] );
               Check_Pass = false;
            }
         }

         else if (  typeid(T) == typeid(float)  ||  typeid(T) == typeid(double)  )
         {
            if ( Fatal_Compr )
            {
               Aux_Error( ERROR_INFO, "\"%s%s\" : input file (%20.14e) != runtime (%20.14e) !!\n",
                          FieldLabel, ArrayIdx,  ((T*)FieldPtr)[t], ComprPtr[t] );
               return -2;
            }

            else
            {
               Aux_Message( stderr, "WARNING : \"%s%s\" : input file (%20.14e) != runtime (%20.14e) !!\n",
                            FieldLabel, ArrayIdx, ((T*)FieldPtr)[t], ComprPtr[t] );
               Check_Pass = false;
            }
         }

         else
         {
            Aux_Message( stderr, "WARNING : \"%s%s\" : unsupported data type !!\n",
                         FieldLabel, ArrayIdx );
            Check_Pass = false;
         }

      } // if ( ((T*)FieldPtr)[t] != ComprPtr[t] )
   } // for (int t=0; t<NCompr; t++)

   return (Check_Pass) ? 0 : -3;

} // FUNCTION : LoadField



//-------------------------------------------------------------------------------------------------------
// Function    :  LoadOnePatch
// Description :  Allocate and load all fields for one patch
//
// Note        :  1. Both leaf and non-leaf patches will store data in the HDF5 output
//                2. If "Recursive == true", this function will be invoked recursively to find all children
//                   (and children's children, ...) patches
//
// Parameter   :  amr              : Targeted AMR_t pointer
//                H5_FileID        : HDF5 file ID of the input file
//                lv               : Target level
//                GID              : Target GID
//                Recursive        : Find all children (and childrens' children, ...) recuresively
//                SonList          : List of son indices
//                                   --> Set only when LOAD_BALANCE is not defined
//                CrList           : List of patch corners
//                H5_SetID_Field   : HDF5 dataset ID           for cell-centered fields
//                H5_SpaceID_Field : HDF5 dataset dataspace ID for cell-centered fields
//                H5_MemID_Field   : HDF5 memory dataspace ID  for cell-centered fields
//                NField           : Number of cell-centered fields stored in the file
//                H5_SetID_Mag     : HDF5 dataset ID           for face-centered magnetic field
//                H5_SpaceID_Mag   : HDF5 dataset dataspace ID for face-centered magnetic field
//                H5_MemID_Mag     : HDF5 memory dataspace ID  for face-centered magnetic field
//                NMag             : Number of face-centered magnetic components stored in the file
//-------------------------------------------------------------------------------------------------------
void LoadOnePatch( AMR_t &amr, const hid_t H5_FileID, const int lv, const int GID,
                   const bool Recursive, const int *SonList, const int (*CrList)[3],
                   const hid_t *H5_SetID_Field, const hid_t  H5_SpaceID_Field, const hid_t  H5_MemID_Field, const int NField,
                   const hid_t *H5_SetID_Mag,   const hid_t *H5_SpaceID_Mag,   const hid_t *H5_MemID_Mag,   const int NMag )
{

   if ( SonList == NULL )  Aux_Error( ERROR_INFO, "SonList == NULL (lv %d, GID %d) !!\n", lv, GID );


   const int  PScale   = PS1*amr.scale[lv];
   const bool WithData = true;   // load data for all patches

   hsize_t H5_Count_Field[4], H5_Offset_Field[4];
   hsize_t H5_Count_Mag[4], H5_Offset_Mag[4];
   herr_t  H5_Status;
   int     SonGID0, PID, CrL[3], CrR[3];


// allocate patch
   for (int d=0; d<3; d++)
   {
      CrL[d] = CrList[GID][d];
      CrR[d] = CrL[d] + PScale;
   }

   amr.pnew( lv, CrL[0], CrL[1], CrL[2], -1, WithData, NField, NMag );


// distinguish leaf and non-leaf patches
// --> but no need to store the real son PID
   PID     = amr.num[lv] - 1;
   SonGID0 = SonList[GID];

   amr.patch[lv][PID]->son = ( SonGID0 == -1 ) ? -1 : 1;


// a. cell-centered fields
// a1. determine the subset of dataspace
   H5_Offset_Field[0] = GID;
   H5_Offset_Field[1] = 0;
   H5_Offset_Field[2] = 0;
   H5_Offset_Field[3] = 0;

   H5_Count_Field [0] = 1;
   H5_Count_Field [1] = PS1;
   H5_Count_Field [2] = PS1;
   H5_Count_Field [3] = PS1;

   H5_Status = H5Sselect_hyperslab( H5_SpaceID_Field, H5S_SELECT_SET, H5_Offset_Field, NULL, H5_Count_Field, NULL );
   if ( H5_Status < 0 )   Aux_Error( ERROR_INFO, "failed to create a hyperslab for the grid data !!\n" );

// a2. load data
   for (int v=0; v<NField; v++)
   {
      H5_Status = H5Dread( H5_SetID_Field[v], H5T_GAMER_REAL, H5_MemID_Field, H5_SpaceID_Field, H5P_DEFAULT,
                           amr.patch[lv][PID]->field[v] );
      if ( H5_Status < 0 )
         Aux_Error( ERROR_INFO, "failed to load a field variable (lv %d, GID %d, v %d) !!\n", lv, GID, v );
   }


// b. face-centered magnetic fields
   for (int v=0; v<NMag; v++)
   {
//    b1. determine the subset of dataspace
      H5_Offset_Mag[0] = GID;
      H5_Offset_Mag[1] = 0;
      H5_Offset_Mag[2] = 0;
      H5_Offset_Mag[3] = 0;

      H5_Count_Mag [0] = 1;
      for (int t=1; t<4; t++)
      H5_Count_Mag [t] = ( 3-t == v ) ? PS1P1 : PS1;

      H5_Status = H5Sselect_hyperslab( H5_SpaceID_Mag[v], H5S_SELECT_SET, H5_Offset_Mag, NULL, H5_Count_Mag, NULL );
      if ( H5_Status < 0 )   Aux_Error( ERROR_INFO, "failed to create a hyperslab for the magnetic field %d !!\n", v );

//    load data
      H5_Status = H5Dread( H5_SetID_Mag[v], H5T_GAMER_REAL, H5_MemID_Mag[v], H5_SpaceID_Mag[v], H5P_DEFAULT,
                           amr.patch[lv][PID]->mag[v] );
      if ( H5_Status < 0 )
         Aux_Error( ERROR_INFO, "failed to load magnetic field (lv %d, GID %d, v %d) !!\n", lv, GID, v );
   } // for (int v=0; v<NMag; v++)


// enter the next level
   if ( Recursive )
   {
      if ( SonGID0 != -1 )
      {
         for (int SonGID=SonGID0; SonGID<SonGID0+8; SonGID++)
            LoadOnePatch( amr, H5_FileID, lv+1, SonGID, Recursive, SonList, CrList,
                          H5_SetID_Field, H5_SpaceID_Field, H5_MemID_Field, NField,
                          H5_SetID_Mag,   H5_SpaceID_Mag,   H5_MemID_Mag,   NMag );
      }
   }

} // FUNCTION : LoadOnePatch



#endif // #ifdef SUPPORT_HDF5
