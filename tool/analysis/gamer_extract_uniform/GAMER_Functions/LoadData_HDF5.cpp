#ifdef SUPPORT_HDF5

#include "ExtractUniform.h"
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
static void LoadOnePatch( const hid_t H5_FileID, const int lv, const int GID, const bool Recursive, const int *SonList,
                          const int (*CrList)[3], const hid_t *H5_SetID_Field, const hid_t H5_SpaceID_Field, const hid_t H5_MemID_Field );




//-------------------------------------------------------------------------------------------------------
// Function    :  LoadData_HDF5
// Description :  Load data from the input HDF5 file
//
// Note        :  1. Only work for format version >= 2200
//
// Parameter   :  FileName : Target file name
//-------------------------------------------------------------------------------------------------------
void LoadData_HDF5( const char *FileName )
{

   if ( MyRank == 0 )   Aux_Message( stdout, "%s ...\n", __FUNCTION__ );


// check
   if ( MyRank == 0 )
   {
      if ( !Aux_CheckFileExist(FileName) )
         Aux_Error( ERROR_INFO, "restart HDF5 file \"%s\" does not exist !!\n", FileName );

      if ( !H5Fis_hdf5(FileName) )
         Aux_Error( ERROR_INFO, "restart file \"%s\" is not in the HDF5 format !!\n", FileName );
   }

   MPI_Barrier( MPI_COMM_WORLD );


// 1. load the simulation info
   if ( MyRank == 0 )   Aux_Message( stdout, "   Loading simulation information ...\n" );

   const bool    Fatal        = true;
   const bool NonFatal        = false;
   const int  Model_RT        = MODEL;
   const int  PatchSize_RT    = PATCH_SIZE;
   const int  NLevel_RT       = NLEVEL;
   const int  NCompFluid_RT   = NCOMP_FLUID;
   const int  NCompPassive_RT = NCOMP_PASSIVE;
#  ifdef FLOAT8
   const int  Float8_RT       = 1;
#  else
   const int  Float8_RT       = 0;
#  endif
   const int MaxString        = 512;

   int    Model_RS, PatchSize_RS, NLevel_RS, NCompFluid_RS, NCompPassive_RS, Float8_RS;
   int    FormatVersion, Gravity, Particle, ExtBC_RS[6], NPatchTotal[NLEVEL], NPatchAllLv;
   int    LoadPot = 0;     // must be integer
   char  *PassiveFieldName_Grid[NCOMP_PASSIVE]; // for format version <  2300
   char  *FieldName_In[NCOMP_TOTAL];            // for format version >= 2300
   int   *NullPtr = NULL;

#  if ( MODEL == ELBDM )
   double ELBDM_Mass;
   double ELBDM_PlanckConst;
#  endif

   hid_t  H5_FileID, H5_SetID_KeyInfo, H5_TypeID_KeyInfo, H5_SetID_InputPara, H5_TypeID_InputPara;
   hid_t  H5_SetID_Makefile, H5_TypeID_Makefile;
   hid_t  H5_SetID_Cr, H5_SetID_Son;
   herr_t H5_Status;


// 1-1. open the HDF5 file
   H5_FileID = H5Fopen( FileName, H5F_ACC_RDONLY, H5P_DEFAULT );
   if ( H5_FileID < 0 )
      Aux_Error( ERROR_INFO, "failed to open the restart HDF5 file \"%s\" !!\n", FileName );


// 1-2. load the dataset and datatype of KeyInfo, InputPara, and Makefile
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

   H5_SetID_Makefile  = H5Dopen( H5_FileID, "Info/Makefile", H5P_DEFAULT );
   if ( H5_SetID_Makefile < 0 )
      Aux_Error( ERROR_INFO, "failed to open the dataset \"%s\" !!\n", "Info/Makefile" );

   H5_TypeID_Makefile = H5Dget_type( H5_SetID_Makefile );
   if ( H5_TypeID_Makefile < 0 )
      Aux_Error( ERROR_INFO, "failed to open the datatype of \"%s\" !!\n", "Info/Makefile" );


// 1-3. load all target fields one-by-one (by all ranks)
   LoadField( "FormatVersion",  &FormatVersion,  H5_SetID_KeyInfo, H5_TypeID_KeyInfo, Fatal,  NullPtr,       -1, NonFatal );

// format version for HDF5 output
   if ( MyRank == 0 )
   {
      Aux_Message( stdout, "      The format version of the HDF5 RESTART file = %ld\n", FormatVersion );

      if ( FormatVersion < 2200 )
         Aux_Error( ERROR_INFO, "unsupported data format version (only support version >= 2200) !!\n" );
   }

   MPI_Barrier( MPI_COMM_WORLD );

// initialize variables that may not exist
   int Comoving_int, WithUnit_int;  // boolean variables are stored as integers in the HDF5 output
   WithUnit_int = 0;
   H0 = Unit_L = Unit_M = Unit_T = Unit_V = Unit_D = Unit_E = Unit_P = MU = WRONG;

   LoadField( "Model",               &Model_RS,          H5_SetID_KeyInfo,   H5_TypeID_KeyInfo,      Fatal, &Model_RT,        1,    Fatal );
   LoadField( "Float8",              &Float8_RS,         H5_SetID_KeyInfo,   H5_TypeID_KeyInfo,      Fatal, &Float8_RT,       1,    Fatal );
   LoadField( "NLevel",              &NLevel_RS,         H5_SetID_KeyInfo,   H5_TypeID_KeyInfo,      Fatal, &NLevel_RT,       1,    Fatal );
   LoadField( "PatchSize",           &PatchSize_RS,      H5_SetID_KeyInfo,   H5_TypeID_KeyInfo,      Fatal, &PatchSize_RT,    1,    Fatal );
   LoadField( "NCompFluid",          &NCompFluid_RS,     H5_SetID_KeyInfo,   H5_TypeID_KeyInfo,      Fatal, &NCompFluid_RT,   1,    Fatal );
   LoadField( "NCompPassive",        &NCompPassive_RS,   H5_SetID_KeyInfo,   H5_TypeID_KeyInfo,      Fatal, &NCompPassive_RT, 1,    Fatal );

   LoadField( "Gravity",             &Gravity,           H5_SetID_KeyInfo,   H5_TypeID_KeyInfo,      Fatal,  NullPtr,        -1, NonFatal );
   LoadField( "Particle",            &Particle,          H5_SetID_KeyInfo,   H5_TypeID_KeyInfo,      Fatal,  NullPtr,        -1, NonFatal );
   LoadField( "DumpID",              &DumpID,            H5_SetID_KeyInfo,   H5_TypeID_KeyInfo,      Fatal,  NullPtr,        -1, NonFatal );
   LoadField( "NX0",                  NX0_TOT,           H5_SetID_KeyInfo,   H5_TypeID_KeyInfo,      Fatal,  NullPtr,        -1, NonFatal );
   LoadField( "NPatch",               NPatchTotal,       H5_SetID_KeyInfo,   H5_TypeID_KeyInfo,      Fatal,  NullPtr,        -1, NonFatal );
   LoadField( "Step",                &Step,              H5_SetID_KeyInfo,   H5_TypeID_KeyInfo,      Fatal,  NullPtr,        -1, NonFatal );
   LoadField( "Time",                 Time,              H5_SetID_KeyInfo,   H5_TypeID_KeyInfo,      Fatal,  NullPtr,        -1, NonFatal );
   LoadField( "CellSize",             amr.dh,            H5_SetID_KeyInfo,   H5_TypeID_KeyInfo,      Fatal,  NullPtr,        -1, NonFatal );
   LoadField( "BoxScale",             amr.BoxScale,      H5_SetID_KeyInfo,   H5_TypeID_KeyInfo,      Fatal,  NullPtr,        -1, NonFatal );
   LoadField( "BoxSize",              amr.BoxSize,       H5_SetID_KeyInfo,   H5_TypeID_KeyInfo,      Fatal,  NullPtr,        -1, NonFatal );

   LoadField( "Comoving",            &Comoving_int,      H5_SetID_Makefile,  H5_TypeID_Makefile,     Fatal,  NullPtr,        -1, NonFatal );
   Comoving = (bool)Comoving_int;

   if ( Comoving ) {
   LoadField( "Hubble0",             &H0,                H5_SetID_InputPara, H5_TypeID_InputPara, NonFatal,  NullPtr,        -1, NonFatal ); }

   LoadField( "Opt__Unit",           &WithUnit_int,      H5_SetID_InputPara, H5_TypeID_InputPara, NonFatal,  NullPtr,        -1, NonFatal );
   WithUnit = (bool)WithUnit_int;

   if ( WithUnit ) {
   LoadField( "Unit_L",              &Unit_L,            H5_SetID_InputPara, H5_TypeID_InputPara,    Fatal,  NullPtr,        -1, NonFatal );
   LoadField( "Unit_M",              &Unit_M,            H5_SetID_InputPara, H5_TypeID_InputPara,    Fatal,  NullPtr,        -1, NonFatal );
   LoadField( "Unit_T",              &Unit_T,            H5_SetID_InputPara, H5_TypeID_InputPara,    Fatal,  NullPtr,        -1, NonFatal );
   LoadField( "Unit_V",              &Unit_V,            H5_SetID_InputPara, H5_TypeID_InputPara,    Fatal,  NullPtr,        -1, NonFatal );
   LoadField( "Unit_D",              &Unit_D,            H5_SetID_InputPara, H5_TypeID_InputPara,    Fatal,  NullPtr,        -1, NonFatal );
   LoadField( "Unit_E",              &Unit_E,            H5_SetID_InputPara, H5_TypeID_InputPara,    Fatal,  NullPtr,        -1, NonFatal );
   LoadField( "Unit_P",              &Unit_P,            H5_SetID_InputPara, H5_TypeID_InputPara,    Fatal,  NullPtr,        -1, NonFatal ); }

   LoadField( "Opt__BC_Flu",          ExtBC_RS,          H5_SetID_InputPara, H5_TypeID_InputPara,    Fatal,  NullPtr,        -1, NonFatal );

   if ( Gravity ) {
   LoadField( "Opt__Output_Pot",     &LoadPot,           H5_SetID_InputPara, H5_TypeID_InputPara,    Fatal,  NullPtr,        -1, NonFatal ); }

   if ( Particle ) {
   LoadField( "Opt__Output_ParDens", &OutputParDens,     H5_SetID_InputPara, H5_TypeID_InputPara,    Fatal,  NullPtr,        -1, NonFatal ); }

#  if   ( MODEL == HYDRO )
   LoadField( "Gamma",               &GAMMA,             H5_SetID_InputPara, H5_TypeID_InputPara,    Fatal,  NullPtr,        -1, NonFatal );
   LoadField( "MolecularWeight",     &MU,                H5_SetID_InputPara, H5_TypeID_InputPara, NonFatal,  NullPtr,        -1, NonFatal );
#  elif ( MODEL == ELBDM )
   LoadField( "ELBDM_Mass",          &ELBDM_Mass,        H5_SetID_InputPara, H5_TypeID_InputPara,    Fatal,  NullPtr,        -1, NonFatal );
   LoadField( "ELBDM_PlanckConst",   &ELBDM_PlanckConst, H5_SetID_InputPara, H5_TypeID_InputPara,    Fatal,  NullPtr,        -1, NonFatal );

   ELBDM_ETA = ELBDM_Mass / ELBDM_PlanckConst;
#  endif

// field labels
   if ( FormatVersion >= 2300 )
   {
      for (int v=0; v<NCOMP_TOTAL; v++)
      {
         char Key[MaxString];
         sprintf( Key, "FieldLabel%02d", v );

         LoadField( Key, &FieldName_In[v], H5_SetID_InputPara, H5_TypeID_InputPara, Fatal, NullPtr, -1, NonFatal );
      }
   }

   else
   {
      for (int v=0; v<NCOMP_PASSIVE; v++)
      {
         char Key[MaxString];
         sprintf( Key, "PassiveFieldName_Grid%02d", v );

         LoadField( Key, &PassiveFieldName_Grid[v], H5_SetID_InputPara, H5_TypeID_InputPara, Fatal, NullPtr, -1, NonFatal );
      }
   }


// ExtBC_RS == 1 and ExtBC == 0 are for periodic BC
   if ( ExtBC_RS[0] != 1  &&  ExtBC == 0 )
      Aux_Error( ERROR_INFO, "please set \"-B 1\" since the simulation domain is non-periodic !!\n" );

   if ( ExtBC_RS[0] == 1  &&  ExtBC != 0  &&  MyRank == 0 )
      Aux_Message( stderr, "WARNING : you might want to set \"-B 0\" since the simulation domain is periodic !!\n" );

   for (int d=0; d<3; d++)    NX0[d] = NX0_TOT[d] / NGPU_X[d];

   NPatchAllLv = 0;
   for (int lv=0; lv<NLEVEL; lv++)  NPatchAllLv += NPatchTotal[lv];


// 1-4. Check whether or not to load the potential and particle density data
   if ( Gravity )
   {
      OutputPot = (bool)LoadPot;          // always output the potential data if they exist

      if ( OutputPot )
      {
         NLoad ++;                        // NLoad and NOut are global variables
         NOut  ++;
      }
   }

   else
      OutputPot = false;

   if ( Particle )
   {
      if ( OutputParDens )
      {
         NLoad ++;                        // NLoad and NOut are global variables
         NOut  ++;
      }
   }

   else
      OutputParDens = 0;


// 1-5. close all objects
   H5_Status = H5Tclose( H5_TypeID_KeyInfo );
   H5_Status = H5Dclose( H5_SetID_KeyInfo );
   H5_Status = H5Tclose( H5_TypeID_InputPara );
   H5_Status = H5Dclose( H5_SetID_InputPara );
   H5_Status = H5Tclose( H5_TypeID_Makefile );
   H5_Status = H5Dclose( H5_SetID_Makefile );
   H5_Status = H5Fclose( H5_FileID );


   MPI_Barrier( MPI_COMM_WORLD );

   if ( MyRank == 0 )   Aux_Message( stdout, "   Loading simulation information ... done\n" );



// 2. set up and check parameters dedicated in this tool
// initialize parameters related to the targeted domain
   Init_TargetDomain();

// allocate memory for the array "BaseP"
   Init_MemAllocate();

// set the conversion factor for calculating temperature
// --> recalculate it only when Convert2Temp < 0.0 so that it does NOT overwrite the value of -u
#  if ( MODEL == HYDRO  ||  MODEL == MHD )
   if ( OutputTemp  &&  Convert2Temp < 0.0 )    Init_Convert2Temp();
#  endif

// verify the input parameters
   if ( MyRank == 0 )   CheckParameter();

// set the candidate box
   GetCandidateBox();

// set the range of the targeted sub-domain
   int TRange_Min[3], TRange_Max[3];

   for (int d=0; d<3; d++)
   {
      TRange_Min[d] = MyRank_X[d]*NX0[d]*amr.scale[0];
      TRange_Max[d] = TRange_Min[d] + NX0[d]*amr.scale[0];
   }



// 3. load the tree information (i.e., corner and son) of all patches (by all ranks)
   H5_FileID = H5Fopen( FileName, H5F_ACC_RDONLY, H5P_DEFAULT );
   if ( H5_FileID < 0 )
      Aux_Error( ERROR_INFO, "failed to open the restart HDF5 file \"%s\" !!\n", FileName );

// 3-1. corner
   if ( MyRank == 0 )   Aux_Message( stdout, "   Loading corner table ...\n" );

// allocate memory
   int (*CrList_AllLv)[3] = new int [ NPatchAllLv ][3];

// load data
   H5_SetID_Cr = H5Dopen( H5_FileID, "Tree/Corner", H5P_DEFAULT );
   if ( H5_SetID_Cr < 0 )     Aux_Error( ERROR_INFO, "failed to open the dataset \"%s\" !!\n", "Tree/Corner" );
   H5_Status = H5Dread( H5_SetID_Cr, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, CrList_AllLv );
   H5_Status = H5Dclose( H5_SetID_Cr );

// shift the corner scale (usually used when the target region lies outside the simulation box)
   if ( Shift2Center )
   {
      for (int t=0; t<NPatchAllLv; t++)
      for (int d=0; d<3; d++)
         CrList_AllLv[t][d] = ( CrList_AllLv[t][d] + ShiftScale[d] + amr.BoxScale[d] ) % amr.BoxScale[d];
   }

   if ( MyRank == 0 )   Aux_Message( stdout, "   Loading corner table ... done\n" );


// 3-2. son
   if ( MyRank == 0 )   Aux_Message( stdout, "   Loading son table ...\n" );

// allocate memory
   int *SonList_AllLv = new int [ NPatchAllLv ];

// load data
   H5_SetID_Son = H5Dopen( H5_FileID, "Tree/Son", H5P_DEFAULT );
   if ( H5_SetID_Son < 0 )    Aux_Error( ERROR_INFO, "failed to open the dataset \"%s\" !!\n", "Tree/Son" );
   H5_Status = H5Dread( H5_SetID_Son, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, SonList_AllLv );
   H5_Status = H5Dclose( H5_SetID_Son );

   if ( MyRank == 0 )   Aux_Message( stdout, "   Loading son table ... done\n" );


   H5_Status = H5Fclose( H5_FileID );



// 4. load and allocate patches
   if ( MyRank == 0 )   Aux_Message( stdout, "   Loading patches ...\n" );

   const bool Recursive_Yes = true;
   const int  NCOMP_ADD     = 2;    // maximum number of additional variables
                                    // --> currently 2: potential and particle/total density

   char (*FieldName)[MaxString] = new char [ NCOMP_TOTAL + NCOMP_ADD ][MaxString];

   hsize_t H5_SetDims_Field[4], H5_MemDims_Field[4];
   hid_t   H5_SetID_Field[ NCOMP_TOTAL + NCOMP_ADD ], H5_MemID_Field, H5_SpaceID_Field, H5_GroupID_GridData;


// 4-1. set the names of all grid variables
   if ( FormatVersion >= 2300 )
   {
      for (int v=0; v<NCOMP_TOTAL; v++)
      sprintf( FieldName[v], FieldName_In[v] );
   }

   else
   {
#     if   ( MODEL == HYDRO )
      sprintf( FieldName[DENS], "Dens" );
      sprintf( FieldName[MOMX], "MomX" );
      sprintf( FieldName[MOMY], "MomY" );
      sprintf( FieldName[MOMZ], "MomZ" );
      sprintf( FieldName[ENGY], "Engy" );

#     elif ( MODEL == ELBDM )
      sprintf( FieldName[DENS], "Dens" );
      sprintf( FieldName[REAL], "Real" );
      sprintf( FieldName[IMAG], "Imag" );

#     else
#     error : ERROR : unsupported MODEL !!
#     endif

      for (int v=0; v<NCOMP_PASSIVE; v++)
      sprintf( FieldName[ NCOMP_FLUID + v ], PassiveFieldName_Grid[v] );
   }

// set the names of potential and particle/total density
   sprintf( FieldName[ NCOMP_TOTAL + 0 ], (OutputPot)?"Pote":"None" );
   sprintf( FieldName[ NCOMP_TOTAL + 1 ], (OutputParDens==0)?"None":
                                          (OutputParDens==1)?"ParDens":
                                          (OutputParDens==2)?"TotalDens":
                                                             "Unknown" );


// 4-2. initialize relevant HDF5 objects
   H5_SetDims_Field[0] = NPatchAllLv;
   H5_SetDims_Field[1] = PATCH_SIZE;
   H5_SetDims_Field[2] = PATCH_SIZE;
   H5_SetDims_Field[3] = PATCH_SIZE;

   H5_SpaceID_Field = H5Screate_simple( 4, H5_SetDims_Field, NULL );
   if ( H5_SpaceID_Field < 0 )   Aux_Error( ERROR_INFO, "failed to create the space \"%s\" !!\n", "H5_SpaceID_Field" );

   H5_MemDims_Field[0] = 1;
   H5_MemDims_Field[1] = PATCH_SIZE;
   H5_MemDims_Field[2] = PATCH_SIZE;
   H5_MemDims_Field[3] = PATCH_SIZE;

   H5_MemID_Field = H5Screate_simple( 4, H5_MemDims_Field, NULL );
   if ( H5_MemID_Field < 0 )  Aux_Error( ERROR_INFO, "failed to create the space \"%s\" !!\n", "H5_MemDims_Field" );


// load data in one rank at a time
   for (int TRank=0; TRank<NGPU; TRank++)
   {
      if ( MyRank == TRank )
      {
//       4-3. open the target datasets just once
         H5_FileID = H5Fopen( FileName, H5F_ACC_RDONLY, H5P_DEFAULT );
         if ( H5_FileID < 0 )
            Aux_Error( ERROR_INFO, "failed to open the restart HDF5 file \"%s\" !!\n", FileName );

         H5_GroupID_GridData = H5Gopen( H5_FileID, "GridData", H5P_DEFAULT );
         if ( H5_GroupID_GridData < 0 )   Aux_Error( ERROR_INFO, "failed to open the group \"%s\" !!\n", "GridData" );

         for (int v=0; v<NCOMP_TOTAL; v++)
         {
            H5_SetID_Field[v] = H5Dopen( H5_GroupID_GridData, FieldName[v], H5P_DEFAULT );
            if ( H5_SetID_Field[v] < 0 )  Aux_Error( ERROR_INFO, "failed to open the dataset \"%s\" !!\n", FieldName[v] );
         }

         if ( OutputPot )
         {
            const int Idx = NCOMP_TOTAL + 0;
            H5_SetID_Field[Idx] = H5Dopen( H5_GroupID_GridData, FieldName[Idx], H5P_DEFAULT );
            if ( H5_SetID_Field[Idx] < 0 )   Aux_Error( ERROR_INFO, "failed to open the dataset \"%s\" !!\n", FieldName[Idx] );
         }

         if ( OutputParDens )
         {
            const int Idx = NCOMP_TOTAL + 1;
            H5_SetID_Field[Idx] = H5Dopen( H5_GroupID_GridData, FieldName[Idx], H5P_DEFAULT );
            if ( H5_SetID_Field[Idx] < 0 )   Aux_Error( ERROR_INFO, "failed to open the dataset \"%s\" !!\n", FieldName[Idx] );
         }


//       4-4. begin to load data
//       loop over the corners of all root-level patches
         const int TenPercent = MAX( NPatchTotal[0]/10, 1 );

         for (int GID=0; GID<NPatchTotal[0]; GID++)
         {
            if ( GID % TenPercent == 0 )
            Aux_Message( stdout, "      Rank %2d: %5.1lf%% completed ...\n",
                         MyRank, 100.0*GID/NPatchTotal[0] );

//          load the info and data of the entire patch family recursively if the root patch is within the target range
            if (  CrList_AllLv[GID][0] >= TRange_Min[0]  &&  CrList_AllLv[GID][0] < TRange_Max[0]  &&
                  CrList_AllLv[GID][1] >= TRange_Min[1]  &&  CrList_AllLv[GID][1] < TRange_Max[1]  &&
                  CrList_AllLv[GID][2] >= TRange_Min[2]  &&  CrList_AllLv[GID][2] < TRange_Max[2]     )
               LoadOnePatch( H5_FileID, 0, GID, Recursive_Yes, SonList_AllLv, CrList_AllLv,
                             H5_SetID_Field, H5_SpaceID_Field, H5_MemID_Field );
         }

//       free resource
         for (int v=0; v<NCOMP_TOTAL; v++)   H5_Status = H5Dclose( H5_SetID_Field[            v] );
         if ( OutputPot )                    H5_Status = H5Dclose( H5_SetID_Field[NCOMP_TOTAL+0] );
         if ( OutputParDens )                H5_Status = H5Dclose( H5_SetID_Field[NCOMP_TOTAL+1] );
         H5_Status = H5Gclose( H5_GroupID_GridData );
         H5_Status = H5Fclose( H5_FileID );
      } // if ( MyRank == TargetRank )

      MPI_Barrier( MPI_COMM_WORLD );
   } // for (int TargetRank=0; TargetRank<NGPU; TargetRank++)

// free HDF5 objects
   H5_Status = H5Sclose( H5_SpaceID_Field );
   H5_Status = H5Sclose( H5_MemID_Field );


// 4-5. record the number of real patches
   for (int lv=0; lv<NLEVEL; lv++)
   {
      NPatchComma[lv][0] = 0;

      for (int m=1; m<28; m++)   NPatchComma[lv][m] = amr.num[lv];
   }


// 4-6. collect the total number of loaded patches at each level
   int NLoadPatch[NLEVEL];
   MPI_Reduce( amr.num, NLoadPatch, NLEVEL, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD );

   if ( MyRank == 0 )
   {
      Aux_Message( stdout, "\n" );
      for (int lv=0; lv<NLEVEL; lv++)     Aux_Message( stdout, "      Lv %2d: %10d patches loaded)\n", lv, NLoadPatch[lv] );
      Aux_Message( stdout, "\n" );

#     ifdef DEBUG_HDF5
      for (int lv=0; lv<=TargetLevel; lv++)
      {
         if ( NLoadPatch[lv] != NPatchTotal[lv] )
            Aux_Error( ERROR_INFO, "Lv %2d: NLoadPatch (%d) != NPatchTotal (%d) !!\n",
                       lv, NLoadPatch[lv], NPatchTotal[lv] );
      }
#     endif

      Aux_Message( stdout, "   Loading patches ... done\n" );
   }



// 5. close all HDF5 objects and free memory
   if ( FieldName     != NULL )  delete [] FieldName;
   if ( CrList_AllLv  != NULL )  delete [] CrList_AllLv;
   if ( SonList_AllLv != NULL )  delete [] SonList_AllLv;



// 6. complete the tree structure
   for (int lv=0; lv<NLEVEL; lv++)
   {
//    construct the relation : father <-> son
      if ( lv > 0 )     FindFather( lv );

//    allocate the buffer patches
      Buf_AllocateBufferPatch( lv );

//    set up the BaseP List
      if ( lv == 0 )    Init_RecordBasePatch();

//    set up the BounP_IDMap
      Buf_RecordBoundaryPatch( lv );

//    construct the sibling relation
      SiblingSearch( lv );

//    get the IDs of patches for sending and receiving data between neighbor ranks
      Buf_RecordExchangeDataPatchID( lv );
   }

// fill up the data in the buffer patches
   for (int lv=0; lv<NLEVEL; lv++)
   {
      Buf_GetBufferData( lv, 1, BufSize );

      if ( OutputPot )
      Buf_GetBufferData( lv, 2, BufSize );

      if ( OutputParDens )
      Buf_GetBufferData( lv, 4, BufSize );
   }


   if ( MyRank == 0 )   Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

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
         Aux_Error( ERROR_INFO, "target field \"%s\" does not exist in the restart file !!\n", FieldName );

      else if ( MyRank == 0 )
         Aux_Message( stderr, "WARNING : target field \"%s\" does not exist in the restart file !!\n", FieldName );

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
               if ( MyRank == 0 )
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
               if ( MyRank == 0 )
               Aux_Message( stderr, "WARNING : \"%s%s\" : RESTART file (%20.14e) != runtime (%20.14e) !!\n",
                            FieldName, ArrayIdx, ((T*)FieldPtr)[t], ComprPtr[t] );
               Check_Pass = false;
            }
         }

         else
         {
            if ( MyRank == 0 )
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
// Parameter   :  H5_FileID         : HDF5 file ID of the restart file
//                lv                : Target level
//                GID               : Target GID
//                Recursive         : Find all children (and childrens' children, ...) recuresively
//                SonList           : List of son indices
//                                    --> Set only when LOAD_BALANCE is not defined
//                CrList            : List of patch corners
//                H5_SetID_Field    : HDF5 dataset ID for grid data
//                H5_SpaceID_Field  : HDF5 dataset dataspace ID for grid data
//                H5_MemID_Field    : HDF5 memory dataspace ID for grid data
//-------------------------------------------------------------------------------------------------------
void LoadOnePatch( const hid_t H5_FileID, const int lv, const int GID, const bool Recursive, const int *SonList,
                   const int (*CrList)[3], const hid_t *H5_SetID_Field, const hid_t H5_SpaceID_Field, const hid_t H5_MemID_Field )
{

   const bool WithData = WithinCandidateBox( CrList[GID], PATCH_SIZE*amr.scale[lv], CanBuf );

   hsize_t H5_Count_Field[4], H5_Offset_Field[4];
   herr_t  H5_Status;
   int     SonGID0, PID;

// allocate patch (note that data arrays are allocated only if this patch lies within the target domain)
   amr.pnew( lv, CrList[GID][0], CrList[GID][1], CrList[GID][2], -1, WithData );

   PID = amr.num[lv] - 1;


// load field data from disk (only if this patch lies within the target domain)
   if ( WithData )
   {
//    determine the subset of dataspace for grid data
      H5_Offset_Field[0] = GID;
      H5_Offset_Field[1] = 0;
      H5_Offset_Field[2] = 0;
      H5_Offset_Field[3] = 0;

      H5_Count_Field [0] = 1;
      H5_Count_Field [1] = PATCH_SIZE;
      H5_Count_Field [2] = PATCH_SIZE;
      H5_Count_Field [3] = PATCH_SIZE;

      H5_Status = H5Sselect_hyperslab( H5_SpaceID_Field, H5S_SELECT_SET, H5_Offset_Field, NULL, H5_Count_Field, NULL );
      if ( H5_Status < 0 )   Aux_Error( ERROR_INFO, "failed to create a hyperslab for the grid data !!\n" );

//    load fluid data
      for (int v=0; v<NCOMP_TOTAL; v++)
      {
         H5_Status = H5Dread( H5_SetID_Field[v], H5T_GAMER_REAL, H5_MemID_Field, H5_SpaceID_Field, H5P_DEFAULT,
                              amr.patch[lv][PID]->fluid[v] );
         if ( H5_Status < 0 )
            Aux_Error( ERROR_INFO, "failed to load a field variable (lv %d, GID %d, v %d) !!\n", lv, GID, v );
      }

//    load potential data
      if ( OutputPot )
      {
         const int Idx = NCOMP_TOTAL + 0;
         H5_Status = H5Dread( H5_SetID_Field[Idx], H5T_GAMER_REAL, H5_MemID_Field, H5_SpaceID_Field, H5P_DEFAULT,
                              amr.patch[lv][PID]->pot );
         if ( H5_Status < 0 )
            Aux_Error( ERROR_INFO, "failed to load the potential data (lv %d, GID %d) !!\n", lv, GID );
      }

//    load particle (or total) density data
      if ( OutputParDens )
      {
         const int Idx = NCOMP_TOTAL + 1;
         H5_Status = H5Dread( H5_SetID_Field[Idx], H5T_GAMER_REAL, H5_MemID_Field, H5_SpaceID_Field, H5P_DEFAULT,
                              amr.patch[lv][PID]->par_dens );
         if ( H5_Status < 0 )
            Aux_Error( ERROR_INFO, "failed to load the particle density data (lv %d, GID %d) !!\n", lv, GID );
      }
   }


// enter the next level
// --> even if this patch lies outside the target domain
// --> in other words, we always allocate ALL patches (but some of them may not have data arrays allocated)
   if ( Recursive )
   {
      if ( SonList == NULL )  Aux_Error( ERROR_INFO, "SonList == NULL (lv %d, GID %d) !!\n", lv, GID );

      SonGID0 = SonList[GID];

      if ( SonGID0 != -1 )
      {
         for (int SonGID=SonGID0; SonGID<SonGID0+8; SonGID++)
            LoadOnePatch( H5_FileID, lv+1, SonGID, Recursive, SonList, CrList,
                          H5_SetID_Field, H5_SpaceID_Field, H5_MemID_Field );
      }
   }

} // FUNCTION : LoadOnePatch



#endif // #ifdef SUPPORT_HDF5
