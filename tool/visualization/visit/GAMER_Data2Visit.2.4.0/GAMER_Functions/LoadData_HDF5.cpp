#include "Data2Visit.h"
#include <typeinfo>

template <typename T>
static herr_t LoadField( const char *FieldName, void *FieldPtr, const hid_t H5_SetID_Target,
                         const hid_t H5_TypeID_Target, const bool Fatal_Nonexist,
                         const T *ComprPtr, const int NCompr, const bool Fatal_Compr );
static void LoadOnePatchGroup( const hid_t H5_FileID, const int lv, const int GID0, const bool LoadPot,
                               const hid_t H5_TypeID_PatchInfo, const hid_t H5_TypeID_PatchData,
                               PatchInfo_t *PatchInfo, PatchData_t *PatchData, const int CanMin[], const int CanMax[],
                               const int MaxDsetPerGroup );




//-------------------------------------------------------------------------------------------------------
// Function    :  LoadData_HDF5
// Description :  Load data from the input HDF5 file
//-------------------------------------------------------------------------------------------------------
void LoadData_HDF5()
{

   cout << "Loading HDF5 data \"" << FileName_In << "\" ..." << endl;


   if ( !Aux_CheckFileExist(FileName_In) )
      Aux_Error( ERROR_INFO, "restart HDF5 file \"%s\" does not exist !!\n", FileName_In );

   if ( !H5Fis_hdf5(FileName_In) )
      Aux_Error( ERROR_INFO, "restart HDF5 file \"%s\" is not in the HDF5 format !!\n", FileName_In );


// 1. load the simulation info
// =================================================================================================
   const bool    Fatal  = true;
   const bool NonFatal  = false;
   const int  Model_RT  = MODEL;
#  ifdef FLOAT8
   const int  Float8_RT = 1;
#  else
   const int  Float8_RT = 0;
#  endif
   const int  NLevel_RT = NLEVEL;

   int    Model_RS, Float8_RS, NLevel_RS, H5_MaxDsetPerGroup;  // RT = RunTime, RS = ReStart
   int    FormatVersion, Gravity, NPatchTotal[NLEVEL], NX0_Tot[3];
   long   Step;
   double BoxSize;
   int    LoadPot = 0;     // must be an integer
   int   *NullPtr = NULL;

#  if ( MODEL == ELBDM )
   double ELBDM_Mass;
   double ELBDM_PlanckConst;
#  endif

   hid_t  H5_FileID;
   hid_t  H5_SetID_KeyInfo, H5_SetID_InputPara, H5_TypeID_KeyInfo, H5_TypeID_InputPara;
   herr_t H5_Status;


   Aux_Message( stdout, "   Loading simulation information ...\n" );


// 1-1. open the HDF5 file
   H5_FileID = H5Fopen( FileName_In, H5F_ACC_RDONLY, H5P_DEFAULT );

   if ( H5_FileID < 0 )
      Aux_Error( ERROR_INFO, "failed to open the restart HDF5 file \"%s\" !!\n", FileName_In );


// 1-2. load the KeyInfo and InputPara datasets and datatypes
   H5_SetID_KeyInfo    = H5Dopen( H5_FileID, "SimuInfo/KeyInfo", H5P_DEFAULT );
   H5_TypeID_KeyInfo   = H5Dget_type( H5_SetID_KeyInfo );

   H5_SetID_InputPara  = H5Dopen( H5_FileID, "SimuInfo/InputPara", H5P_DEFAULT );
   H5_TypeID_InputPara = H5Dget_type( H5_SetID_InputPara );


// 1-3. load all target fields in KeyInfo one-by-one (by all ranks)
   LoadField( "FormatVersion",     &FormatVersion,     H5_SetID_KeyInfo,   H5_TypeID_KeyInfo,   Fatal,  NullPtr,   -1, NonFatal );

// format version for HDF5 output starts from 2000
   Aux_Message( stdout, "      The format version of the HDF5 RESTART file = %ld\n", FormatVersion );

   if ( FormatVersion < 2000 )
      Aux_Error( ERROR_INFO, "unsupported data format version (only support version >= 2000) !!\n" );

   LoadField( "Model",             &Model_RS,          H5_SetID_KeyInfo,   H5_TypeID_KeyInfo,   Fatal, &Model_RT,   1,    Fatal );
   LoadField( "Float8",            &Float8_RS,         H5_SetID_KeyInfo,   H5_TypeID_KeyInfo,   Fatal, &Float8_RT,  1,    Fatal );
   LoadField( "NLevel",            &NLevel_RS,         H5_SetID_KeyInfo,   H5_TypeID_KeyInfo,   Fatal, &NLevel_RT,  1,    Fatal );

   LoadField( "Gravity",           &Gravity,           H5_SetID_KeyInfo,   H5_TypeID_KeyInfo,   Fatal,  NullPtr,   -1, NonFatal );
   if ( Gravity ) {
   LoadField( "OutputPot",         &LoadPot,           H5_SetID_KeyInfo,   H5_TypeID_KeyInfo,   Fatal,  NullPtr,   -1, NonFatal ); }
   LoadField( "DumpID",            &DumpID,            H5_SetID_KeyInfo,   H5_TypeID_KeyInfo,   Fatal,  NullPtr,   -1, NonFatal );
   LoadField( "NX0",                NX0_Tot,           H5_SetID_KeyInfo,   H5_TypeID_KeyInfo,   Fatal,  NullPtr,   -1, NonFatal );
   LoadField( "NPatch",             NPatchTotal,       H5_SetID_KeyInfo,   H5_TypeID_KeyInfo,   Fatal,  NullPtr,   -1, NonFatal );
   LoadField( "Step",              &Step,              H5_SetID_KeyInfo,   H5_TypeID_KeyInfo,   Fatal,  NullPtr,   -1, NonFatal );
   LoadField( "Time",               Time,              H5_SetID_KeyInfo,   H5_TypeID_KeyInfo,   Fatal,  NullPtr,   -1, NonFatal );
   LoadField( "CellSize",           amr.dh,            H5_SetID_KeyInfo,   H5_TypeID_KeyInfo,   Fatal,  NullPtr,   -1, NonFatal );
   LoadField( "BoxScale",           amr.BoxScale,      H5_SetID_KeyInfo,   H5_TypeID_KeyInfo,   Fatal,  NullPtr,   -1, NonFatal );
   LoadField( "BoxSize",            amr.BoxSize,       H5_SetID_KeyInfo,   H5_TypeID_KeyInfo,   Fatal,  NullPtr,   -1, NonFatal );
   LoadField( "H5_MaxDsetPerGroup",&H5_MaxDsetPerGroup,H5_SetID_KeyInfo,   H5_TypeID_KeyInfo,   Fatal,  NullPtr,   -1, NonFatal );

#  if   ( MODEL == HYDRO )
   LoadField( "Gamma",             &GAMMA,             H5_SetID_InputPara, H5_TypeID_InputPara, Fatal,  NullPtr,   -1, NonFatal );

#  elif ( MODEL == ELBDM )
   LoadField( "ELBDM_Mass",        &ELBDM_Mass,        H5_SetID_InputPara, H5_TypeID_InputPara, Fatal,  NullPtr,   -1, NonFatal );
   LoadField( "ELBDM_PlanckConst", &ELBDM_PlanckConst, H5_SetID_InputPara, H5_TypeID_InputPara, Fatal,  NullPtr,   -1, NonFatal );

   ELBDM_ETA = ELBDM_Mass / ELBDM_PlanckConst;
#  endif


// 1-4. Check whether or not to load the potential data
   if ( Gravity )
      OutputPot = (bool)LoadPot;          // always output the potential data if they exist

   else
      OutputPot = false;


// 1-5. close all objects
   H5_Status = H5Tclose( H5_TypeID_KeyInfo );
   H5_Status = H5Dclose( H5_SetID_KeyInfo );
   H5_Status = H5Tclose( H5_TypeID_InputPara );
   H5_Status = H5Dclose( H5_SetID_InputPara );


   Aux_Message( stdout, "      DumpID      = %d\n",      DumpID          );
   Aux_Message( stdout, "      Time        = %13.7e\n",  Time[0]         );
   Aux_Message( stdout, "      Step        = %ld\n",     Step            );
   Aux_Message( stdout, "      NX0_Tot [0] = %8d\n",     NX0_Tot[0]      );
   Aux_Message( stdout, "      NX0_Tot [1] = %8d\n",     NX0_Tot[1]      );
   Aux_Message( stdout, "      NX0_Tot [2] = %8d\n",     NX0_Tot[2]      );
   Aux_Message( stdout, "      BoxSize [0] = %14.5lf\n", amr.BoxSize[0]  );
   Aux_Message( stdout, "      BoxSize [1] = %14.5lf\n", amr.BoxSize[1]  );
   Aux_Message( stdout, "      BoxSize [2] = %14.5lf\n", amr.BoxSize[2]  );
   Aux_Message( stdout, "      BoxScale[0] = %8d\n",     amr.BoxScale[0] );
   Aux_Message( stdout, "      BoxScale[1] = %8d\n",     amr.BoxScale[1] );
   Aux_Message( stdout, "      BoxScale[2] = %8d\n",     amr.BoxScale[2] );


   Aux_Message( stdout, "   Loading simulation information ... done\n" );



// 2. set up and check parameters dedicated in this tool
// =================================================================================================
// set parameters to the default values
   SetDefaultParameter();

// truncate the targeted region
   TruncateBox();

// check parameters
   CheckParameter();

// set up the candidate box, only patches with cells inside the candidate box will be allocated
   int CanMax[3], CanMin[3];
   for (int dim=0; dim<3; dim++)
   {
      CanMin[dim] = CornerScale[dim];
      CanMax[dim] = CornerScale[dim] + SizeScale[dim];
   }



// 3. load the simulation data
// =================================================================================================
   Aux_Message( stdout, "   Loading data by walking tree ...\n" );

// 3-1. get the datatypes of patch info and data
   hid_t H5_AttID_PatchInfo, H5_SetID_Patch, H5_TypeID_PatchInfo, H5_TypeID_PatchData;

   char PatchName[100];
   PatchInfo_t PatchInfo;
   PatchData_t PatchData;

   H5_SetID_Patch = H5Dopen( H5_FileID, "Level_00/Cluster_000000000/Patch_000000000", H5P_DEFAULT );

   if ( H5_SetID_Patch < 0 )
      Aux_Error( ERROR_INFO, "failed to open the patch dataset \"%s\" !!\n", "Level_00/Cluster_000000000/Patch_000000000" );

   H5_AttID_PatchInfo = H5Aopen( H5_SetID_Patch, "Info", H5P_DEFAULT );

   if ( H5_AttID_PatchInfo < 0 )
      Aux_Error( ERROR_INFO, "failed to open the info attribute of \"%s\" !!\n", "Level_00/Cluster_000000000/Patch_000000000" );

   H5_TypeID_PatchData = H5Dget_type( H5_SetID_Patch );
   H5_TypeID_PatchInfo = H5Aget_type( H5_AttID_PatchInfo );

   H5_Status = H5Aclose( H5_AttID_PatchInfo );
   H5_Status = H5Dclose( H5_SetID_Patch );


// 3-2. load data
   const int TenPercent = MAX( NPatchTotal[0]/10-NPatchTotal[0]/10%8, 1 );

   for (int GID0=0; GID0<NPatchTotal[0]; GID0+=8)
   {
      if ( GID0%TenPercent == 0 )
         Aux_Message( stdout, "      %5.1f%% completed ...\n", 100.0*(double)GID0/NPatchTotal[0] );

      LoadOnePatchGroup( H5_FileID, 0, GID0, OutputPot, H5_TypeID_PatchInfo, H5_TypeID_PatchData,
                         &PatchInfo, &PatchData, CanMin, CanMax, H5_MaxDsetPerGroup );
   }

   Aux_Message( stdout, "   Loading data by walking tree ... done\n" );

   Aux_Message( stdout, "\n" );
   for (int lv=0; lv<NLEVEL; lv++)     Aux_Message( stdout, "      Lv %2d: %10d patches loaded\n", lv, amr.num[lv] );
   Aux_Message( stdout, "\n" );


// 3-3. close all HDF5 objects and free memory (especially for the variable-length string)
   H5_Status = H5Fclose( H5_FileID );
   H5_Status = H5Tclose( H5_TypeID_PatchInfo );
   H5_Status = H5Tclose( H5_TypeID_PatchData );


   cout << "Loading HDF5 data \"" << FileName_In << "\" ... done" << endl;

} // FUNCTION : LoadData_HDFt



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
         Aux_Error( ERROR_INFO, "target field \"%s\" does not exist !!\n", FieldName );

      else
         Aux_Message( stderr, "WARNING : target field \"%s\" does not exist !!\n", FieldName );

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
               Aux_Error( ERROR_INFO, "\"%s%s\" : RESTART file (%ld) != runtime (%ld) !!\n",
                          FieldName, ArrayIdx, (long)((T*)FieldPtr)[t], (long)ComprPtr[t] );
               return -2;
            }

            else
            {
               Aux_Message( stderr, "WARNING : \"%s%s\" : RESTART file (%ld) != runtime (%ld) !!\n",
                            FieldName, ArrayIdx, (long)((T*)FieldPtr)[t], (long)ComprPtr[t] );
               Check_Pass = false;
            }
         }

         else if (  typeid(T) == typeid(float)  ||  typeid(T) == typeid(double)  )
         {
            if ( Fatal_Compr )
            {
               Aux_Error( ERROR_INFO, "\"%s%s\" : RESTART file (%20.14e) != runtime (%20.14e) !!\n",
                          FieldName, ArrayIdx,  ((T*)FieldPtr)[t], ComprPtr[t] );
               return -2;
            }

            else
            {
               Aux_Message( stderr, "WARNING : \"%s%s\" : RESTART file (%20.14e) != runtime (%20.14e) !!\n",
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
// Function    :  LoadOnePatchGroup
// Description :  Load one patch group from the GAMER HDF5 output
//
// Note        :  1. In this function one should specify the criteria for loading a patch group
//                2. Both leaf and non-leaf patches will store data in the HDF5 output
//                3. This function will be invoked recursively to load all child patches
//
// Parameter   :  H5_FileID            : HDF5 file ID of the input file
//                lv                   : Target level
//                GID0                 : GID of the target patch group
//                LoadPot              : True --> load the gravitational potential
//                H5_TypeID_PatchInfo  : HDF5 type ID of patch info
//                H5_TypeID_PatchData  : HDF5 type ID of patch data
//                PatchInfo            : Pointer for storing the patch info
//                PatchData            : Pointer for storing the patch data
//                CanMax/Min           : Maximum/Minimum scale of the target region
//                MaxDsetPerGroup      : Maximum number of datasets per group
//-------------------------------------------------------------------------------------------------------
void LoadOnePatchGroup( const hid_t H5_FileID, const int lv, const int GID0, const bool LoadPot,
                        const hid_t H5_TypeID_PatchInfo, const hid_t H5_TypeID_PatchData,
                        PatchInfo_t *PatchInfo, PatchData_t *PatchData, const int CanMin[], const int CanMax[],
                        const int MaxDsetPerGroup )
{

   const int OneVarSize  = CUBE(PS1)*sizeof(real);
   const int PScale      = PATCH_SIZE*amr.scale[lv];

   int  CrL[3], CrR[3], PID, SonGID0, CID;
   char PatchName[100];
   bool Within;

   hid_t  H5_AttID_PatchInfo, H5_SetID_Patch;
   herr_t H5_Status;


// loop over all local patches
   for (int GID=GID0; GID<GID0+8; GID++)
   {
//    0. get the cluster ID
      CID = GID / MaxDsetPerGroup;


//    1. set the names of the target datasets
      sprintf( PatchName, "Level_%d%d/Cluster_%d%d%d%d%d%d%d%d%d/Patch_%d%d%d%d%d%d%d%d%d",
               lv/10, lv%10, CID/100000000,
                            (CID%100000000)/10000000,
                            (CID%10000000 )/1000000,
                            (CID%1000000  )/100000,
                            (CID%100000   )/10000,
                            (CID%10000    )/1000,
                            (CID%1000     )/100,
                            (CID%100      )/10,
                            (CID%10       ),
                             GID/100000000,
                            (GID%100000000)/10000000,
                            (GID%10000000 )/1000000,
                            (GID%1000000  )/100000,
                            (GID%100000   )/10000,
                            (GID%10000    )/1000,
                            (GID%1000     )/100,
                            (GID%100      )/10,
                            (GID%10       )           );


//    2. load the patch info
      H5_SetID_Patch = H5Dopen( H5_FileID, PatchName, H5P_DEFAULT );

      if ( H5_SetID_Patch < 0 )  Aux_Error( ERROR_INFO, "failed to open the patch dataset \"%s\" !!\n", PatchName );

      H5_AttID_PatchInfo = H5Aopen( H5_SetID_Patch, "Info", H5P_DEFAULT );

      if ( H5_AttID_PatchInfo < 0 )    Aux_Error( ERROR_INFO, "failed to open the patch attribute \"%s\" !!\n", PatchName );

      H5_Status = H5Aread( H5_AttID_PatchInfo, H5_TypeID_PatchInfo, PatchInfo );
      H5_Status = H5Aclose( H5_AttID_PatchInfo );
      H5_Status = H5Dclose( H5_SetID_Patch );


//    3. check whether this patch lies within the target range
      for (int d=0; d<3; d++)
      {
         CrL[d] = PatchInfo->Corner[d];
         CrR[d] = CrL[d] + PScale;
      }

      if ( CrL[0] < CanMax[0]  &&  CrR[0] > CanMin[0]  &&
           CrL[1] < CanMax[1]  &&  CrR[1] > CanMin[1]  &&
           CrL[2] < CanMax[2]  &&  CrR[2] > CanMin[2]    )
         Within = true;

      else
         Within = false;


//    4. check whether it's a leaf patch or not
      SonGID0 = PatchInfo->Son;


//    5. allocate patch and load data (only for leaf patches lying within the target domain)
      if ( Within )
      {
         if ( SonGID0 == -1 )
         {
//          allocate patch
            PID = amr.num[lv];
            amr.pnew( lv, CrL[0], CrL[1], CrL[2] );

//          load the patch data (only if it's a leaf patch lying within the target domain)
            H5_SetID_Patch = H5Dopen( H5_FileID, PatchName, H5P_DEFAULT );
            H5_Status      = H5Dread( H5_SetID_Patch, H5_TypeID_PatchData, H5S_ALL, H5S_ALL, H5P_DEFAULT, PatchData );
            H5_Status      = H5Dclose( H5_SetID_Patch );

#           if   ( MODEL == HYDRO )
            memcpy( amr.patch[lv][PID]->fluid[DENS], PatchData->Dens, OneVarSize );
            memcpy( amr.patch[lv][PID]->fluid[MOMX], PatchData->MomX, OneVarSize );
            memcpy( amr.patch[lv][PID]->fluid[MOMY], PatchData->MomY, OneVarSize );
            memcpy( amr.patch[lv][PID]->fluid[MOMZ], PatchData->MomZ, OneVarSize );
            memcpy( amr.patch[lv][PID]->fluid[ENGY], PatchData->Engy, OneVarSize );

#           elif ( MODEL == MHD )
#           warning : WAIT MHD !!!

#           elif ( MODEL == ELBDM )
            memcpy( amr.patch[lv][PID]->fluid[DENS], PatchData->Dens, OneVarSize );
            memcpy( amr.patch[lv][PID]->fluid[REAL], PatchData->Real, OneVarSize );
            memcpy( amr.patch[lv][PID]->fluid[IMAG], PatchData->Imag, OneVarSize );

#           else
#           error : unsupported MODEL !!
#           endif // MODEL

            if ( LoadPot )
            memcpy( amr.patch[lv][PID]->pot,         PatchData->Pote, OneVarSize );
         } // if ( SonGID0 == -1 )

//       find leaf patches recursively
         else
         LoadOnePatchGroup( H5_FileID, lv+1, SonGID0, LoadPot, H5_TypeID_PatchInfo, H5_TypeID_PatchData,
                            PatchInfo, PatchData, CanMin, CanMax, MaxDsetPerGroup );
      } // if ( Within )
   } // for (int GID=GID0; GID<GID0+8; GID++)

} // FUNCTION : LoadOnePatchGroup
