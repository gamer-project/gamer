#ifdef SUPPORT_HDF5

#include "CompareData.h"
#include "hdf5.h"
#include <typeinfo>

#ifdef FLOAT8
#  define H5T_GAMER_REAL H5T_NATIVE_DOUBLE
#else
#  define H5T_GAMER_REAL H5T_NATIVE_FLOAT
#endif

#ifdef GAMER_DEBUG
#  define DEBUG_HDF5
#endif

template <typename T>
static herr_t LoadField( const char *FieldName, void *FieldPtr, const hid_t H5_SetID_Target,
                         const hid_t H5_TypeID_Target, const bool Fatal_Nonexist,
                         const T *ComprPtr, const int NCompr, const bool Fatal_Compr );
static void LoadOnePatch( AMR_t &amr, const hid_t H5_FileID, const int lv, const int GID, const bool Recursive, const int *SonList,
                          const int (*CrList)[3], const hid_t *H5_SetID_Field, const hid_t H5_SpaceID_Field, const hid_t H5_MemID_Field,
                          const bool WithPot, const int WithParDens );




//-------------------------------------------------------------------------------------------------------
// Function    :  LoadData_HDF5
// Description :  Load data from the input HDF5 file
//
// Note        :  1. Only work for format version >= 2300
//                2. Unlike LoadData(), this function always loads and allocates both leaf and non-leaf patches
//
// Parameter   :  amr         : Targeted AMR_t pointer
//                FileName    : Name of the input file
//                WithPot     : true --> the loaded data contain potential field
//                WithParDens : > 0  --> the loaded data contain particle density on grids
//                WithPar     : true --> the loaded data contain particles
//                NParVarOut  : Number of particle attributes stored in the file
//                NPar        : NUmber of particles
//                ParData     : Particle data array (allocated here --> must be dallocated manually later)
//                WithMagCC   : true --> the loaded data contain cell-centered magnetic field
//                WithMagFC   : true --> the loaded data contain face-centered magnetic field
//
// Return      :  amr, WithPot, WithParDens, WithPar, NParVarOut, NPar, ParData, WithMagCC, WithMagFC
//-------------------------------------------------------------------------------------------------------
void LoadData_HDF5( AMR_t &amr, const char *FileName, bool &WithPot, int &WithParDens, bool &WithPar,
                    int &NParVarOut, long &NPar, real **&ParData, bool &WithMagCC, bool &WithMagFC )
{

   Aux_Message( stdout, "Loading HDF5 data %s ...\n", FileName );


// check
   if ( !Aux_CheckFileExist(FileName) )
      Aux_Error( ERROR_INFO, "input file \"%s\" does not exist !!\n", FileName );

   if ( !H5Fis_hdf5(FileName) )
      Aux_Error( ERROR_INFO, "input file \"%s\" is not in the HDF5 format !!\n", FileName );


// 1. load the simulation info
   Aux_Message( stdout, "   Loading simulation information ...\n" );

   const bool    Fatal        = true;
   const bool NonFatal        = false;
   const int  Model_RT        = MODEL;
   const int  PatchSize_RT    = PS1;
   const int  NLevel_RT       = NLEVEL;
   const int  NCompFluid_RT   = NCOMP_FLUID;
   const int  NCompPassive_RT = NCOMP_PASSIVE;
#  ifdef FLOAT8
   const int  Float8_RT       = 1;
#  else
   const int  Float8_RT       = 0;
#  endif
   const int  MaxString       = 512;

   int    Model_RS, PatchSize_RS, NLevel_RS, NCompFluid_RS, NCompPassive_RS, Float8_RS;
   int    FormatVersion, Gravity, NPatchTotal[NLEVEL], NPatchAllLv;
   char  *FieldName_In[NCOMP_TOTAL];
   int   *NullPtr = NULL;

#  if ( MODEL == HYDRO )
   int    Magnetohydrodynamics;
#  endif

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

   if ( FormatVersion < 2300 )
      Aux_Error( ERROR_INFO, "unsupported data format version (only support version >= 2300) !!\n" );

   LoadField( "Model",                &Model_RS,            H5_SetID_KeyInfo,    H5_TypeID_KeyInfo,   Fatal, &Model_RT,        1,    Fatal );
   LoadField( "Float8",               &Float8_RS,           H5_SetID_KeyInfo,    H5_TypeID_KeyInfo,   Fatal, &Float8_RT,       1,    Fatal );
   LoadField( "NLevel",               &NLevel_RS,           H5_SetID_KeyInfo,    H5_TypeID_KeyInfo,   Fatal, &NLevel_RT,       1,    Fatal );
   LoadField( "PatchSize",            &PatchSize_RS,        H5_SetID_KeyInfo,    H5_TypeID_KeyInfo,   Fatal, &PatchSize_RT,    1,    Fatal );
   LoadField( "NCompFluid",           &NCompFluid_RS,       H5_SetID_KeyInfo,    H5_TypeID_KeyInfo,   Fatal, &NCompFluid_RT,   1,    Fatal );
   LoadField( "NCompPassive",         &NCompPassive_RS,     H5_SetID_KeyInfo,    H5_TypeID_KeyInfo,   Fatal, &NCompPassive_RT, 1,    Fatal );
   LoadField( "Gravity",              &Gravity,             H5_SetID_KeyInfo,    H5_TypeID_KeyInfo,   Fatal,  NullPtr,        -1, NonFatal );
   LoadField( "Particle",             &WithPar,             H5_SetID_KeyInfo,    H5_TypeID_KeyInfo,   Fatal,  NullPtr,        -1, NonFatal );
   LoadField( "Par_NAttStored",       &NParVarOut,          H5_SetID_KeyInfo,    H5_TypeID_KeyInfo,   Fatal,  NullPtr,        -1, NonFatal );
   LoadField( "Par_NPar",             &NPar,                H5_SetID_KeyInfo,    H5_TypeID_KeyInfo,   Fatal,  NullPtr,        -1, NonFatal );
   LoadField( "NPatch",                NPatchTotal,         H5_SetID_KeyInfo,    H5_TypeID_KeyInfo,   Fatal,  NullPtr,        -1, NonFatal );
   LoadField( "NX0",                   amr.nx0_tot,         H5_SetID_KeyInfo,    H5_TypeID_KeyInfo,   Fatal,  NullPtr,        -1, NonFatal );

   if ( Gravity )
   LoadField( "Opt__Output_Pot",      &WithPot,             H5_SetID_InputPara,  H5_TypeID_InputPara, Fatal,  NullPtr,        -1, NonFatal );
   else
   WithPot = false;

   if ( WithPar )
   LoadField( "Opt__Output_ParDens",  &WithParDens,         H5_SetID_InputPara,  H5_TypeID_InputPara, Fatal,  NullPtr,        -1, NonFatal );
   else
   WithParDens = 0;

// check if B field is included
#  if ( MODEL == HYDRO )
   if ( FormatVersion >= 2400 )
   LoadField( "Magnetohydrodynamics", &Magnetohydrodynamics,H5_SetID_KeyInfo,    H5_TypeID_KeyInfo,   Fatal,   NullPtr,       -1, NonFatal );
   else
   Magnetohydrodynamics = 0;
#  endif

// field labels
   for (int v=0; v<NCOMP_TOTAL; v++)
   {
      char Key[MaxString];
      sprintf( Key, "FieldLabel%02d", v );

      LoadField( Key,                 &FieldName_In[v],     H5_SetID_InputPara,  H5_TypeID_InputPara,    Fatal,   NullPtr,    -1, NonFatal );
   }

   NPatchAllLv = 0;
   for (int lv=0; lv<NLEVEL; lv++)  NPatchAllLv += NPatchTotal[lv];


// 1-4. close all objects
   H5_Status = H5Tclose( H5_TypeID_KeyInfo );
   H5_Status = H5Dclose( H5_SetID_KeyInfo );
   H5_Status = H5Tclose( H5_TypeID_InputPara );
   H5_Status = H5Dclose( H5_SetID_InputPara );


   Aux_Message( stdout, "   Loading simulation information ... done\n" );



// 3. load the tree information (i.e., corner and son) of all patches
// 3-1. corner
   Aux_Message( stdout, "   Loading corner table ...\n" );

// allocate memory
   int (*CrList_AllLv)[3] = new int [ NPatchAllLv ][3];

// load data
   H5_SetID_Cr = H5Dopen( H5_FileID, "Tree/Corner", H5P_DEFAULT );
   if ( H5_SetID_Cr < 0 )     Aux_Error( ERROR_INFO, "failed to open the dataset \"%s\" !!\n", "Tree/Corner" );
   H5_Status = H5Dread( H5_SetID_Cr, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, CrList_AllLv );
   H5_Status = H5Dclose( H5_SetID_Cr );

   Aux_Message( stdout, "   Loading corner table ... done\n" );


// 3-2. son
   Aux_Message( stdout, "   Loading son table ...\n" );

// allocate memory
   int *SonList_AllLv = new int [ NPatchAllLv ];

// load data
   H5_SetID_Son = H5Dopen( H5_FileID, "Tree/Son", H5P_DEFAULT );
   if ( H5_SetID_Son < 0 )    Aux_Error( ERROR_INFO, "failed to open the dataset \"%s\" !!\n", "Tree/Son" );
   H5_Status = H5Dread( H5_SetID_Son, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, SonList_AllLv );
   H5_Status = H5Dclose( H5_SetID_Son );

   Aux_Message( stdout, "   Loading son table ... done\n" );



// 4. load and allocate patches
   Aux_Message( stdout, "   Loading patches ...\n" );

   const bool Recursive_Yes = true;
   const int  NCOMP_ADD     = 2;    // maximum number of additional variables
                                    // --> currently 2: potential and particle/total density

   char (*FieldName)[MaxString] = new char [ NCOMP_TOTAL + NCOMP_ADD ][MaxString];

   hsize_t H5_SetDims_Field[4], H5_MemDims_Field[4];
   hid_t   H5_SetID_Field[ NCOMP_TOTAL + NCOMP_ADD ], H5_MemID_Field, H5_SpaceID_Field, H5_GroupID_GridData;


// 4-1. set the names of all grid variables
   for (int v=0; v<NCOMP_TOTAL; v++)
   sprintf( FieldName[v], FieldName_In[v] );

// set the names of potential and particle/total density
   sprintf( FieldName[ NCOMP_TOTAL + 0 ], (WithPot)?"Pote":"None" );
   sprintf( FieldName[ NCOMP_TOTAL + 1 ], (WithParDens==0)?"None":
                                          (WithParDens==1)?"ParDens":
                                          (WithParDens==2)?"TotalDens":
                                                           "Unknown" );


// 4-2. initialize relevant HDF5 objects
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


// 4-3. open the target datasets just once
   H5_GroupID_GridData = H5Gopen( H5_FileID, "GridData", H5P_DEFAULT );
   if ( H5_GroupID_GridData < 0 )   Aux_Error( ERROR_INFO, "failed to open the group \"%s\" !!\n", "GridData" );

   for (int v=0; v<NCOMP_TOTAL; v++)
   {
      H5_SetID_Field[v] = H5Dopen( H5_GroupID_GridData, FieldName[v], H5P_DEFAULT );
      if ( H5_SetID_Field[v] < 0 )  Aux_Error( ERROR_INFO, "failed to open the dataset \"%s\" !!\n", FieldName[v] );
   }

   if ( WithPot )
   {
      const int v = NCOMP_TOTAL + 0;
      H5_SetID_Field[v] = H5Dopen( H5_GroupID_GridData, FieldName[v], H5P_DEFAULT );
      if ( H5_SetID_Field[v] < 0 )  Aux_Error( ERROR_INFO, "failed to open the dataset \"%s\" !!\n", FieldName[v] );
   }

   if ( WithParDens )
   {
      const int v = NCOMP_TOTAL + 1;
      H5_SetID_Field[v] = H5Dopen( H5_GroupID_GridData, FieldName[v], H5P_DEFAULT );
      if ( H5_SetID_Field[v] < 0 )  Aux_Error( ERROR_INFO, "failed to open the dataset \"%s\" !!\n", FieldName[v] );
   }


// 4-4. begin to load data
// loop over all root-level patches
   const int TenPercent = MAX( NPatchTotal[0]/10, 1 );

   for (int GID=0; GID<NPatchTotal[0]; GID++)
   {
      if ( GID % TenPercent == 0 )
      Aux_Message( stdout, "      %5.1lf%% completed ...\n", 100.0*GID/NPatchTotal[0] );

//    load the entire patch family recursively (actually it's not necessary here since we load all patches anyway)
      LoadOnePatch( amr, H5_FileID, 0, GID, Recursive_Yes, SonList_AllLv, CrList_AllLv,
                    H5_SetID_Field, H5_SpaceID_Field, H5_MemID_Field, WithPot, WithParDens );
   }

// free HDF5 objects
   for (int v=0; v<NCOMP_TOTAL; v++)   H5_Status = H5Dclose( H5_SetID_Field[            v] );
   if ( WithPot )                      H5_Status = H5Dclose( H5_SetID_Field[NCOMP_TOTAL+0] );
   if ( WithParDens )                  H5_Status = H5Dclose( H5_SetID_Field[NCOMP_TOTAL+1] );

   H5_Status = H5Gclose( H5_GroupID_GridData );
   H5_Status = H5Fclose( H5_FileID );
   H5_Status = H5Sclose( H5_SpaceID_Field );
   H5_Status = H5Sclose( H5_MemID_Field );


// 4-5. record the total number of loaded patches at each level
   Aux_Message( stdout, "\n" );
   for (int lv=0; lv<NLEVEL; lv++)     Aux_Message( stdout, "      Lv %2d: %10d patches loaded)\n", lv, amr.num[lv] );
   Aux_Message( stdout, "\n" );

   for (int lv=0; lv<NLEVEL; lv++)
   {
      if ( amr.num[lv] != NPatchTotal[lv] )
         Aux_Error( ERROR_INFO, "Lv %2d: amr.num (%d) != NPatchTotal (%d) !!\n",
                    lv, amr.num[lv], NPatchTotal[lv] );
   }

   Aux_Message( stdout, "   Loading patches ... done\n" );



// 5. close all HDF5 objects and free memory
   delete [] FieldName;
   delete [] CrList_AllLv;
   delete [] SonList_AllLv;


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
// Parameter   :  FieldName         : Name of the target field
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
herr_t LoadField( const char *FieldName, void *FieldPtr, const hid_t H5_SetID_Target,
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
   H5_FieldIdx = H5Tget_member_index( H5_TypeID_Target, FieldName );

   if ( H5_FieldIdx >= 0 )
   {
      H5_TypeID_Field  = H5Tget_member_type( H5_TypeID_Target, H5_FieldIdx );
      H5_FieldSize     = H5Tget_size( H5_TypeID_Field );

      H5_TypeID_Load   = H5Tcreate( H5T_COMPOUND, H5_FieldSize );
      H5_Status        = H5Tinsert( H5_TypeID_Load, FieldName, 0, H5_TypeID_Field );

      H5_Status        = H5Dread( H5_SetID_Target, H5_TypeID_Load, H5S_ALL, H5S_ALL, H5P_DEFAULT, FieldPtr );
      if ( H5_Status < 0 )    Aux_Error( ERROR_INFO, "failed to load the field \"%s\" !!\n", FieldName );

      H5_Status        = H5Tclose( H5_TypeID_Field );
      H5_Status        = H5Tclose( H5_TypeID_Load  );
   } // if ( H5_FieldIdx >= 0 )

   else
   {
      if ( Fatal_Nonexist )
         Aux_Error( ERROR_INFO, "target field \"%s\" does not exist in the input file !!\n", FieldName );

      else
         Aux_Message( stderr, "WARNING : target field \"%s\" does not exist in the input file !!\n", FieldName );

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
                          FieldName, ArrayIdx, (long)((T*)FieldPtr)[t], (long)ComprPtr[t] );
               return -2;
            }

            else
            {
               Aux_Message( stderr, "WARNING : \"%s%s\" : input file (%ld) != runtime (%ld) !!\n",
                            FieldName, ArrayIdx, (long)((T*)FieldPtr)[t], (long)ComprPtr[t] );
               Check_Pass = false;
            }
         }

         else if (  typeid(T) == typeid(float)  ||  typeid(T) == typeid(double)  )
         {
            if ( Fatal_Compr )
            {
               Aux_Error( ERROR_INFO, "\"%s%s\" : input file (%20.14e) != runtime (%20.14e) !!\n",
                          FieldName, ArrayIdx,  ((T*)FieldPtr)[t], ComprPtr[t] );
               return -2;
            }

            else
            {
               Aux_Message( stderr, "WARNING : \"%s%s\" : input file (%20.14e) != runtime (%20.14e) !!\n",
                            FieldName, ArrayIdx, ((T*)FieldPtr)[t], ComprPtr[t] );
               Check_Pass = false;
            }
         }

         else
         {
            Aux_Message( stderr, "WARNING : \"%s%s\" : unsupported data type !!\n",
                         FieldName, ArrayIdx );
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
//                H5_SetID_Field   : HDF5 dataset ID for grid data
//                H5_SpaceID_Field : HDF5 dataset dataspace ID for grid data
//                H5_MemID_Field   : HDF5 memory dataspace ID for grid data
//                WithPot          : load potential field
//                WithParDens      : load particle density on grids
//-------------------------------------------------------------------------------------------------------
void LoadOnePatch( AMR_t &amr, const hid_t H5_FileID, const int lv, const int GID, const bool Recursive, const int *SonList,
                   const int (*CrList)[3], const hid_t *H5_SetID_Field, const hid_t H5_SpaceID_Field, const hid_t H5_MemID_Field,
                   const bool WithPot, const int WithParDens )

{

   if ( SonList == NULL )  Aux_Error( ERROR_INFO, "SonList == NULL (lv %d, GID %d) !!\n", lv, GID );


   const int  PScale   = PS1*amr.scale[lv];
   const bool WithData = true;   // load data for all patches

   hsize_t H5_Count_Field[4], H5_Offset_Field[4];
   herr_t  H5_Status;
   int     SonGID0, PID, CrL[3], CrR[3];


// allocate patch
   for (int d=0; d<3; d++)
   {
      CrL[d] = CrList[GID][d];
      CrR[d] = CrL[d] + PScale;
   }

   amr.pnew( lv, CrL[0], CrL[1], CrL[2], -1, WithData );


// distinguish leaf and non-leaf patches
// --> but no need to store the real son PID
   PID     = amr.num[lv] - 1;
   SonGID0 = SonList[GID];

   amr.patch[lv][PID]->son = ( SonGID0 == -1 ) ? -1 : 1;


// load field data from disk
   if ( WithData )
   {
//    determine the subset of dataspace for grid data
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

//    load the fluid data
      for (int v=0; v<NCOMP_TOTAL; v++)
      {
         H5_Status = H5Dread( H5_SetID_Field[v], H5T_GAMER_REAL, H5_MemID_Field, H5_SpaceID_Field, H5P_DEFAULT,
                              amr.patch[lv][PID]->fluid[v] );
         if ( H5_Status < 0 )
            Aux_Error( ERROR_INFO, "failed to load a field variable (lv %d, GID %d, v %d) !!\n", lv, GID, v );
      }

//    load the potential data
      if ( WithPot )
      {
         const int v = NCOMP_TOTAL + 0;
         H5_Status = H5Dread( H5_SetID_Field[v], H5T_GAMER_REAL, H5_MemID_Field, H5_SpaceID_Field, H5P_DEFAULT,
                              amr.patch[lv][PID]->pot );
         if ( H5_Status < 0 )
            Aux_Error( ERROR_INFO, "failed to load the potential data (lv %d, GID %d) !!\n", lv, GID );
      }

//    load the particle (or total) density data
      if ( WithParDens )
      {
         const int v = NCOMP_TOTAL + 1;
         H5_Status = H5Dread( H5_SetID_Field[v], H5T_GAMER_REAL, H5_MemID_Field, H5_SpaceID_Field, H5P_DEFAULT,
                              amr.patch[lv][PID]->par_dens );
         if ( H5_Status < 0 )
            Aux_Error( ERROR_INFO, "failed to load the particle density data (lv %d, GID %d) !!\n", lv, GID );
      }
   }


// enter the next level
   if ( Recursive )
   {
      if ( SonGID0 != -1 )
      {
         for (int SonGID=SonGID0; SonGID<SonGID0+8; SonGID++)
            LoadOnePatch( amr, H5_FileID, lv+1, SonGID, Recursive, SonList, CrList,
                          H5_SetID_Field, H5_SpaceID_Field, H5_MemID_Field, WithPot, WithParDens );
      }
   }

} // FUNCTION : LoadOnePatch



#endif // #ifdef SUPPORT_HDF5
