#include "ExtractUniform.h"
#ifdef SUPPORT_HDF5
#include "hdf5.h"
#endif

static void Load_Parameter_Before_2000( FILE *File, const int FormatVersion, bool &DataOrder_xyzv,
                                        bool &LoadPot, int *NX0_Tot, double &BoxSize, double &Gamma, double &ELBDM_Eta,
                                        const long HeaderOffset_Makefile, const long HeaderOffset_Constant,
                                        const long HeaderOffset_Parameter, bool &Comoving );
static void Load_Parameter_After_2000( FILE *File, const int FormatVersion, bool &LoadPot, int *NX0_Tot,
                                       double &BoxSize, double &Gamma, double &ELBDM_Eta,
                                       const long HeaderOffset_Makefile, const long HeaderOffset_Constant,
                                       const long HeaderOffset_Parameter, bool &LoadPar, int &LoadParDens, int &NParVarOut,
                                       double &H0, bool &WithUnit, double &Unit_L, double &Unit_M, double &Unit_T, double &Unit_V,
                                       double &Unit_D, double &Unit_E, double &Unit_P, double &Mu, bool &Comoving );
static void CompareVar( const char *VarName, const bool   RestartVar, const bool   RuntimeVar, const bool Fatal );
static void CompareVar( const char *VarName, const int    RestartVar, const int    RuntimeVar, const bool Fatal );
static void CompareVar( const char *VarName, const long   RestartVar, const long   RuntimeVar, const bool Fatal );
static void CompareVar( const char *VarName, const real   RestartVar, const real   RuntimeVar, const bool Fatal );
static void CompareVar( const char *VarName, const double RestartVar, const double RuntimeVar, const bool Fatal );
static void LoadOnePatch( FILE *File, const int lv, const int LoadPID, const bool LoadPot, const bool LoadParDens );

#ifdef SUPPORT_HDF5
void LoadData_HDF5( const char *FileName );
#endif




//-------------------------------------------------------------------------------------------------------
// Function    :  LoadData
// Description :  Load an output data of GAMER
//
// Parameter   :  None
//-------------------------------------------------------------------------------------------------------
void LoadData()
{

// load the HDF5 data
#  ifdef SUPPORT_HDF5
   if (  Aux_CheckFileExist(FileName_In)  &&  H5Fis_hdf5(FileName_In)  )
   {
      LoadData_HDF5( FileName_In );
      return;
   }
#  endif


   if ( MyRank == 0 )   cout << "LoadData \"" << FileName_In << "\" ... "<< endl;


   FILE *File = fopen( FileName_In, "rb" );

   if ( File == NULL  &&  MyRank == 0 )   Aux_Error( ERROR_INFO, "input file \"%s\" does not exist !!\n", FileName_In );


// initialize the NPatchComma list as 0
   for (int lv=0; lv<NLEVEL; lv++)  for (int m=0; m<28; m++)      NPatchComma[lv][m] = 0;


// record the size of different data types
   const int size_bool   = sizeof( bool   );
   const int size_int    = sizeof( int    );
   const int size_uint   = sizeof( uint   );
   const int size_long   = sizeof( long   );
   const int size_ulong  = sizeof( ulong  );
   const int size_real   = sizeof( real   );
   const int size_double = sizeof( double );

   if ( size_int != size_uint )
   {
      fprintf( stderr, "ERROR : sizeof(int) = %d != sizeof(uint) = %d !!\n", size_int, size_uint );
      MPI_Exit();
   }

   if ( size_long != size_ulong )
   {
      fprintf( stderr, "ERROR : sizeof(long) = %d != sizeof(ulong) = %d !!\n", size_long, size_ulong );
      MPI_Exit();
   }


// a. load the information of data format
// =================================================================================================
   long FormatVersion, CheckCode;
   long HeaderSize_Format, HeaderSize_Makefile, HeaderSize_Constant, HeaderSize_Parameter, HeaderSize_SimuInfo;
   long HeaderOffset_Format, HeaderOffset_Makefile, HeaderOffset_Constant, HeaderOffset_Parameter, HeaderOffset_SimuInfo;
   long HeaderSize_Total;


// verify the input data format version
   fread( &FormatVersion, sizeof(long), 1, File );

   if ( MyRank == 0 )
   Aux_Message( stdout, "   The format version of the input file = %ld\n", FormatVersion );

   if ( FormatVersion < 1200 )
      Aux_Error( ERROR_INFO, "unsupported data format version (only support version >= 1200) !!\n" );

   else if ( FormatVersion >= 2000  &&  UseTree )
      Aux_Error( ERROR_INFO, "FormatVersion >= 2000 does NOT support the \"UseTree\" option (-T) !!\n" );


// determine the header size
   if ( FormatVersion < 2000 )
   {
      fseek( File, sizeof(long), SEEK_CUR );

      HeaderSize_Format    =  256;
      HeaderSize_Makefile  =  256;
      HeaderSize_Constant  =  256;
      HeaderSize_Parameter = 1280;  // since "HeaderSize" = 2048 for FormatVersion < 2000
      HeaderSize_SimuInfo  = 1024;
   }

   else
   {
      fread( &HeaderSize_Format,    sizeof(long), 1, File );
      fread( &HeaderSize_Makefile,  sizeof(long), 1, File );
      fread( &HeaderSize_Constant,  sizeof(long), 1, File );
      fread( &HeaderSize_Parameter, sizeof(long), 1, File );
      fread( &HeaderSize_SimuInfo,  sizeof(long), 1, File );
   }

// determine the offset of each header (in bytes)
   HeaderOffset_Format    = 0;    // it must be zero
   HeaderOffset_Makefile  = HeaderOffset_Format    + HeaderSize_Format;
   HeaderOffset_Constant  = HeaderOffset_Makefile  + HeaderSize_Makefile;
   HeaderOffset_Parameter = HeaderOffset_Constant  + HeaderSize_Constant;
   HeaderOffset_SimuInfo  = HeaderOffset_Parameter + HeaderSize_Parameter;

   HeaderSize_Total       = HeaderOffset_SimuInfo  + HeaderSize_SimuInfo;


// load the check code
   fread( &CheckCode, sizeof(long), 1, File );


// check if the size of different data types are consistent (only for version >= 1200)
   int size_bool_restart, size_int_restart, size_long_restart, size_real_restart, size_double_restart;

   fread( &size_bool_restart,   sizeof(int), 1, File );
   fread( &size_int_restart,    sizeof(int), 1, File );
   fread( &size_long_restart,   sizeof(int), 1, File );
   fread( &size_real_restart,   sizeof(int), 1, File );
   fread( &size_double_restart, sizeof(int), 1, File );

   if ( size_bool_restart != size_bool )
   {
      fprintf( stderr, "ERROR : sizeof(bool) is inconsistent : RESTART file = %d, runtime = %d !!\n",
               size_bool_restart, size_bool );
      MPI_Exit();
   }

   if ( size_int_restart != size_int )
   {
      fprintf( stderr, "ERROR : sizeof(int) is inconsistent : RESTART file = %d, runtime = %d !!\n",
               size_int_restart, size_int );
      MPI_Exit();
   }

   if ( size_long_restart != size_long )
   {
      fprintf( stderr, "ERROR : sizeof(long) is inconsistent : RESTART file = %d, runtime = %d !!\n",
               size_long_restart, size_long );
      MPI_Exit();
   }

   if ( size_real_restart != size_real )
   {
      fprintf( stderr, "ERROR : sizeof(real) is inconsistent : RESTART file = %d, runtime = %d !!\n",
               size_real_restart, size_real );
      MPI_Exit();
   }

   if ( size_double_restart != size_double )
   {
      fprintf( stderr, "ERROR : sizeof(double) is inconsistent : RESTART file = %d, runtime = %d !!\n",
               size_double_restart, size_double );
      MPI_Exit();
   }



// b. load all simulation parameters
// =================================================================================================
   bool   DataOrder_xyzv, OutputPar;
   int    NParVarOut;
   double BoxSize;
#  if ( MODEL != ELBDM )
   double ELBDM_ETA;
#  endif

   if ( FormatVersion < 2000 )
   {
      Load_Parameter_Before_2000( File, FormatVersion, DataOrder_xyzv, OutputPot, NX0_TOT, BoxSize, GAMMA, ELBDM_ETA,
                                  HeaderOffset_Makefile, HeaderOffset_Constant, HeaderOffset_Parameter, Comoving );
      OutputPar     = false;
      OutputParDens = 0;
      NParVarOut    = -1;

      WithUnit      = false;
      H0 = Unit_L = Unit_M = Unit_T = Unit_V = Unit_D = Unit_E = Unit_P = MU = WRONG;
   }

   else
   {
      Load_Parameter_After_2000( File, FormatVersion, OutputPot, NX0_TOT, BoxSize, GAMMA, ELBDM_ETA,
                                 HeaderOffset_Makefile, HeaderOffset_Constant, HeaderOffset_Parameter,
                                 OutputPar, OutputParDens, NParVarOut,
                                 H0, WithUnit, Unit_L, Unit_M, Unit_T, Unit_V, Unit_D, Unit_E, Unit_P, MU, Comoving );
      DataOrder_xyzv = false;
   }



// c. load the simulation information
// =================================================================================================
// verify the check code
   long checkcode;

   fseek( File, HeaderOffset_SimuInfo, SEEK_SET );

   fread( &checkcode, sizeof(long), 1, File );

   if ( checkcode != CheckCode )
      Aux_Error( ERROR_INFO, "incorrect check code in the RESTART file (input %ld <-> expeect %ld) !!\n",
                 checkcode, CheckCode );


// load information necessary for restart
   long AdvanceCounter[NLEVEL], NPar;
   int  NPatchTotal[NLEVEL], NDataPatch_Total[NLEVEL];

   fread( &DumpID,          sizeof(int),         1, File );
   fread( Time,             sizeof(double), NLEVEL, File );
   fread( &Step,            sizeof(long),        1, File );
   fread( NPatchTotal,      sizeof(int),    NLEVEL, File );
   fread( NDataPatch_Total, sizeof(int),    NLEVEL, File );
   fread( AdvanceCounter,   sizeof(long),   NLEVEL, File );
   fseek( File, sizeof(double), SEEK_CUR );

   if ( OutputPar )
   fread( &NPar,            sizeof(long),        1, File );


// verify the size of the RESTART file
   if ( MyRank == 0 )   Aux_Message( stdout, "   Verifying the file size ...\n" );

   long DataSize[NLEVEL], ExpectSize, InputSize, PatchDataSize;

   if ( OutputPot )
   {
      NLoad ++;
      NOut  ++;
   }

   if ( OutputParDens )
   {
      NLoad ++;
      NOut  ++;
   }

   PatchDataSize = PATCH_SIZE*PATCH_SIZE*PATCH_SIZE*NLoad*sizeof(real);
   ExpectSize    = HeaderSize_Total;

   for (int lv=0; lv<NLEVEL; lv++)
   {
      DataSize[lv]  = 0;
      DataSize[lv] += NPatchTotal[lv]*4*sizeof(int);        // 4 = corner(3) + son(1)
      DataSize[lv] += NDataPatch_Total[lv]*PatchDataSize;

      ExpectSize   += DataSize[lv];
   }

   if ( OutputPar )
   {
      for (int lv=0; lv<NLEVEL; lv++)
      {
//       2 = NPar + starting particle index stored in each leaf patch
         const long ParInfoSize = (long)NDataPatch_Total[lv]*2*sizeof(long);

         DataSize[lv] += ParInfoSize;
         ExpectSize   += ParInfoSize;
      }

      ExpectSize += (long)NParVarOut*NPar*sizeof(real);
   }

   fseek( File, 0, SEEK_END );
   InputSize = ftell( File );

   if ( InputSize != ExpectSize )
      Aux_Error( ERROR_INFO, "the size of the file <%s> is incorrect --> input = %ld <-> expect = %ld !!\n",
                 FileName_In, InputSize, ExpectSize );

   fclose( File );

   if ( MyRank == 0 )   Aux_Message( stdout, "   Verifying the file size ... passed\n" );


// set up the simulation box
   int NX0_Max;
   NX0_Max = ( NX0_TOT[0] > NX0_TOT[1] ) ? NX0_TOT[0] : NX0_TOT[1];
   NX0_Max = ( NX0_TOT[2] > NX0_Max    ) ? NX0_TOT[2] : NX0_Max;

   for (int lv=0; lv<NLEVEL; lv++)     amr.dh[lv] = BoxSize / (double)( NX0_Max*(1<<lv) );

   for (int d=0; d<3; d++)
   {
      amr.BoxSize [d] = NX0_TOT[d]*amr.dh   [0];
      amr.BoxScale[d] = NX0_TOT[d]*amr.scale[0];
   }

   NX0[0] = NX0_TOT[0] / NGPU_X[0];
   NX0[1] = NX0_TOT[1] / NGPU_X[1];
   NX0[2] = NX0_TOT[2] / NGPU_X[2];



// d. set up and check parameters dedicated in this tool
// =================================================================================================
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
   int TargetRange_Min[3], TargetRange_Max[3];

   for (int d=0; d<3; d++)
   {
      TargetRange_Min[d] = MyRank_X[d]*NX0[d]*amr.scale[0];
      TargetRange_Max[d] = TargetRange_Min[d] + NX0[d]*amr.scale[0];
   }



// e. load the simulation data
// =================================================================================================
   long Offset = HeaderSize_Total;
   int  LoadCorner[3], LoadSon, PID;
   bool GotYou;

// array for re-ordering the fluid data from "xyzv" to "vxyz"
   real (*InvData_Flu)[PATCH_SIZE][PATCH_SIZE][NCOMP_TOTAL] = NULL;
   if ( DataOrder_xyzv )
   {
      if ( UseTree )    Aux_Error( ERROR_INFO, "UseTree option does not support DataOrder_xyzv !!\n" );

      InvData_Flu = new real [PATCH_SIZE][PATCH_SIZE][PATCH_SIZE][NCOMP_TOTAL];
   }


   if ( UseTree )
   {
      if ( MyRank == 0 )   Aux_Message( stdout, "   Loading data by walking tree ...\n" );

      for (int TargetRank=0; TargetRank<NGPU; TargetRank++)
      {
         if ( MyRank == TargetRank )
         {
            File = fopen( FileName_In, "rb" );

//          loop over all root-level patches
            const int TenPercent = MAX( NPatchTotal[0]/10, 1 );

            for (int LoadPID=0; LoadPID<NPatchTotal[0]; LoadPID++)
            {
               if ( LoadPID%TenPercent == 0 )
               Aux_Message( stdout, "      %5.1f%% completed ...\n", 100.0*(double)LoadPID/NPatchTotal[0] );

               for (int d=0; d<3; d++)    LoadCorner[d]= tree->block[0][LoadPID].corner[d];

//             shift the corner scale (usually used when the target region lies outside the simulation box)
               if ( Shift2Center )
               {
                  for (int d=0; d<3; d++)
                     LoadCorner[d] = ( LoadCorner[d] + ShiftScale[d] + amr.BoxScale[d] ) % amr.BoxScale[d];
               }

//             verify that the loaded patch is within the target range of MyRank
               if (  LoadCorner[0] >= TargetRange_Min[0]  &&  LoadCorner[0] < TargetRange_Max[0]  &&
                     LoadCorner[1] >= TargetRange_Min[1]  &&  LoadCorner[1] < TargetRange_Max[1]  &&
                     LoadCorner[2] >= TargetRange_Min[2]  &&  LoadCorner[2] < TargetRange_Max[2]     )
                  LoadOnePatch( File, 0, LoadPID, OutputPot, OutputParDens );
            }

            fclose( File );
         } // if ( MyRank == TargetRank )

         MPI_Barrier( MPI_COMM_WORLD );
      } // for (int TargetRank=0; TargetRank<NGPU; TargetRank++)

//    collect the total number of loaded blocks at each level
      int NLoadBlock[NLEVEL];
      MPI_Reduce( amr.num, NLoadBlock, NLEVEL, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD );

      if ( MyRank == 0 )
      {
         Aux_Message( stdout, "\n" );
         for (int lv=0; lv<NLEVEL; lv++)     Aux_Message( stdout, "      Lv %2d: %10d blocks loaded)\n", lv, NLoadBlock[lv] );
         Aux_Message( stdout, "\n" );

#        ifdef GAMER_DEBUG
         for (int lv=0; lv<NLEVEL; lv++)
         {
            if ( NLoadBlock[lv] != tree->nblock[lv] )
               Aux_Error( ERROR_INFO, "Lv %2d: NLoadBlock (%d) != tree->nblock (%d) !!\n",
                          lv, NLoadBlock[lv], tree->nblock[lv] );
         }
#        endif

         Aux_Message( stdout, "   Loading data by walking tree ... done\n" );
      }
   } // if ( UseTree )

   else
   for (int lv=0; lv<NLEVEL; lv++)
   {
      for (int TargetRank=0; TargetRank<NGPU; TargetRank++)
      {
         if ( MyRank == 0 )
         {
            fprintf( stdout, "   Loading data: level %2d, MyRank %3d ... ", lv, TargetRank );
            fflush( stdout );
         }

         if ( MyRank == TargetRank )
         {
            File = fopen( FileName_In, "rb" );
            fseek( File, Offset, SEEK_SET );

            for (int LoadPID=0; LoadPID<NPatchTotal[lv]; LoadPID++)
            {
//             e1. load the patch information
               fread(  LoadCorner, sizeof(int), 3, File );
               fread( &LoadSon,    sizeof(int), 1, File );


//             shift the corner scale (usually used when the target region lies outside the simulation box)
               if ( Shift2Center )
               {
                  for (int d=0; d<3; d++)
                     LoadCorner[d] = ( LoadCorner[d] + ShiftScale[d] + amr.BoxScale[d] ) % amr.BoxScale[d];
               }


//             verify that the loaded patch is within the target range of MyRank
               if (  LoadCorner[0] >= TargetRange_Min[0]  &&  LoadCorner[0] < TargetRange_Max[0]  &&
                     LoadCorner[1] >= TargetRange_Min[1]  &&  LoadCorner[1] < TargetRange_Max[1]  &&
                     LoadCorner[2] >= TargetRange_Min[2]  &&  LoadCorner[2] < TargetRange_Max[2]     )
               {

//                verify that the loaded patch is within the candidate box
                  GotYou = WithinCandidateBox( LoadCorner, PATCH_SIZE*amr.scale[lv], CanBuf );

                  amr.pnew( lv, LoadCorner[0], LoadCorner[1], LoadCorner[2], -1, GotYou );

//                e2. load the physical data if it is a leaf patch
                  if ( LoadSon == -1 )
                  {
                     if ( GotYou )
                     {
                        PID = amr.num[lv] - 1;

//                      e2-0. skip particle information
                        if ( OutputPar )  fseek( File, 2*sizeof(long), SEEK_CUR );

//                      e2-1. load the fluid variables
                        if ( DataOrder_xyzv )
                        {
                           fread( InvData_Flu, sizeof(real), PATCH_SIZE*PATCH_SIZE*PATCH_SIZE*NCOMP_TOTAL, File );

                           for (int v=0; v<NCOMP_TOTAL; v++)
                           for (int k=0; k<PATCH_SIZE; k++)
                           for (int j=0; j<PATCH_SIZE; j++)
                           for (int i=0; i<PATCH_SIZE; i++)
                              amr.patch[lv][PID]->fluid[v][k][j][i] = InvData_Flu[k][j][i][v];
                        }

                        else
                           fread( amr.patch[lv][PID]->fluid,    sizeof(real),
                                  PATCH_SIZE*PATCH_SIZE*PATCH_SIZE*NCOMP_TOTAL, File );

//                      e2-2. load the gravitational potential
                        if ( OutputPot )
                           fread( amr.patch[lv][PID]->pot,      sizeof(real), PATCH_SIZE*PATCH_SIZE*PATCH_SIZE, File );

//                      e2-3. load the particle density on grids
                        if ( OutputParDens )
                           fread( amr.patch[lv][PID]->par_dens, sizeof(real), PATCH_SIZE*PATCH_SIZE*PATCH_SIZE, File );
                     } // if ( GotYou )

                     else
                        fseek( File, PatchDataSize + ( (OutputPar)?2*sizeof(long):0 ), SEEK_CUR );
                  } // if ( LoadSon == -1 )
               } // LoadCorner within target range of MyRank

               else if ( LoadSon == -1 )  fseek( File, PatchDataSize + ( (OutputPar)?2*sizeof(long):0 ), SEEK_CUR );

            } // for (int LoadPID=0; LoadPID<NPatchTotal[lv]; LoadPID++)

            fclose( File );

            Offset += DataSize[lv];

         } // if ( MyRank == TargetRank )

         MPI_Barrier( MPI_COMM_WORLD );

         if ( MyRank == 0 )
         {
            fprintf( stdout, "done\n" );
            fflush( stdout );
         }

      } // for (int TargetRank=0; TargetRank<NGPU; TargetRank++)
   } // for (int lv=0; lv<NLEVEL; lv++)


   if ( DataOrder_xyzv )  delete [] InvData_Flu;


// record the number of the real patches
   for (int lv=0; lv<NLEVEL; lv++)
   for (int m=1; m<28; m++)               NPatchComma[lv][m] = amr.num[lv];


// complete all levels
//=====================================================================================
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


// fill up the data for patches that are not leaf patches
   for (int lv=NLEVEL-2; lv>=0; lv--)     Flu_Restrict( lv, OutputPot, OutputParDens );


// fill up the data in the buffer patches
   for (int lv=0; lv<NLEVEL; lv++)
   {
      Buf_GetBufferData( lv, 1, BufSize );

      if ( OutputPot )
      Buf_GetBufferData( lv, 2, BufSize );

      if ( OutputParDens )
      Buf_GetBufferData( lv, 4, BufSize );
   }


   if ( MyRank == 0 )   cout << "LoadData ... done" << endl;

} // FUNCTION : LoadData



//-------------------------------------------------------------------------------------------------------
// Function    :  Load_Parameter_Before_2000
// Description :  Load all simulation parameters from the RESTART file with format version < 2000 (but >= 1200)
//
// Note        :  "OPT__RESTART_HEADER == RESTART_HEADER_INHERIC" can be used in this function
//
// Parameter   :  File           : RESTART file pointer
//                FormatVersion  : Format version of the RESTART file
//                DataOrder_xyzv : Order of data stored in the RESTART file (true/false --> xyzv/vxyz)
//                LoadPot        : Whether or not the RESTART file stores the potential data
//                NX0_Tot        : Total number of base-level cells along each direction
//                BoxSize        : Physical size of the simulation box
//                Gamma          : ratio of specific heat
//                ELBDM_Eta      : Particle mass / Planck constant in ELBDM
//                HeaderOffset_X : Offsets of different headers
//                Comoving       : true --> cosmological simulations in comoving coords.
//
// Return      :  DataOrder_xyzv, LoadPot, NX0_Tot, BoxSize, Gamma, ELBDM_Eta, Comoving
//-------------------------------------------------------------------------------------------------------
void Load_Parameter_Before_2000( FILE *File, const int FormatVersion, bool &DataOrder_xyzv,
                                 bool &LoadPot, int *NX0_Tot, double &BoxSize, double &Gamma, double &ELBDM_Eta,
                                 const long HeaderOffset_Makefile, const long HeaderOffset_Constant,
                                 const long HeaderOffset_Parameter, bool &Comoving )
{

   if ( MyRank == 0 )   Aux_Message( stdout, "   Loading simulation parameters ...\n" );


// a. load the simulation options and parameters defined in the Makefile
// =================================================================================================
   bool gravity, individual_timestep, comoving, gpu, gamer_optimization, gamer_debug, timing, timing_solver;
   bool intel, float8, serial, ooc, overlap_mpi, openmp, fermi;
   int  model, pot_scheme, flu_scheme, lr_scheme, rsolver, load_balance, nlevel, max_patch;

   fseek( File, HeaderOffset_Makefile, SEEK_SET );

   fread( &model,                      sizeof(int),                     1,             File );
   fread( &gravity,                    sizeof(bool),                    1,             File );
   fread( &pot_scheme,                 sizeof(int),                     1,             File );
   fread( &individual_timestep,        sizeof(bool),                    1,             File );
   fread( &comoving,                   sizeof(bool),                    1,             File );
   fread( &flu_scheme,                 sizeof(int),                     1,             File );
   fread( &lr_scheme,                  sizeof(int),                     1,             File );
   fread( &rsolver,                    sizeof(int),                     1,             File );
   fread( &gpu,                        sizeof(bool),                    1,             File );
   fread( &gamer_optimization,         sizeof(bool),                    1,             File );
   fread( &gamer_debug,                sizeof(bool),                    1,             File );
   fread( &timing,                     sizeof(bool),                    1,             File );
   fread( &timing_solver,              sizeof(bool),                    1,             File );
   fread( &intel,                      sizeof(bool),                    1,             File );
   fread( &float8,                     sizeof(bool),                    1,             File );
   fread( &serial,                     sizeof(bool),                    1,             File );
   fread( &ooc,                        sizeof(bool),                    1,             File );
   fread( &load_balance,               sizeof(int),                     1,             File );
   fread( &overlap_mpi,                sizeof(bool),                    1,             File );
   fread( &openmp,                     sizeof(bool),                    1,             File );
   fread( &fermi,                      sizeof(bool),                    1,             File );
   fread( &nlevel,                     sizeof(int),                     1,             File );
   fread( &max_patch,                  sizeof(int),                     1,             File );


// b. load the symbolic constants defined in "Macro.h, CUPOT.h, and CUFLU.h"
// =================================================================================================
   bool enforce_positive, char_reconstruction, hll_no_ref_state, hll_include_all_waves, waf_dissipate;
   bool use_psolver_10to14;
   int  ncomp_fluid, patch_size, flu_ghost_size, pot_ghost_size, gra_ghost_size, check_intermediate;
   int  flu_block_size_x, flu_block_size_y, pot_block_size_x, pot_block_size_z, gra_block_size_z;
   real min_value, max_error;

   fseek( File, HeaderOffset_Constant, SEEK_SET );

   fread( &ncomp_fluid,                sizeof(int),                     1,             File );
   fread( &patch_size,                 sizeof(int),                     1,             File );
   fread( &min_value,                  sizeof(real),                    1,             File );
   fread( &flu_ghost_size,             sizeof(int),                     1,             File );
   fread( &pot_ghost_size,             sizeof(int),                     1,             File );
   fread( &gra_ghost_size,             sizeof(int),                     1,             File );
   fread( &enforce_positive,           sizeof(bool),                    1,             File );
   fread( &char_reconstruction,        sizeof(bool),                    1,             File );
   fread( &check_intermediate,         sizeof(int),                     1,             File );
   fread( &hll_no_ref_state,           sizeof(bool),                    1,             File );
   fread( &hll_include_all_waves,      sizeof(bool),                    1,             File );
   fread( &waf_dissipate,              sizeof(bool),                    1,             File );
   fread( &max_error,                  sizeof(real),                    1,             File );
   fread( &flu_block_size_x,           sizeof(int),                     1,             File );
   fread( &flu_block_size_y,           sizeof(int),                     1,             File );
   fread( &use_psolver_10to14,         sizeof(bool),                    1,             File );
   fread( &pot_block_size_x,           sizeof(int),                     1,             File );
   fread( &pot_block_size_z,           sizeof(int),                     1,             File );
   fread( &gra_block_size_z,           sizeof(int),                     1,             File );


// c. load the simulation parameters recorded in the file "Input__Parameter"
// =================================================================================================
   bool   opt__adaptive_dt, opt__dt_user, opt__flag_rho, opt__flag_rho_gradient, opt__flag_pres_gradient;
   bool   opt__flag_engy_density, opt__flag_user, opt__fixup_flux, opt__fixup_restrict, opt__overlap_mpi;
   bool   opt__gra_p5_gradient, opt__int_time, opt__output_error, opt__output_base, opt__output_pot;
   bool   opt__timing_barrier, opt__int_phase;
   int    nx0_tot[3], gamer_nrank, gamer_nrank_x[3], omp_nthread, ooc_nrank, ooc_nrank_x[3], regrid_count;
   int    flag_buffer_size, max_level, opt__lr_limiter, opt__waf_limiter, flu_gpu_npgroup, gpu_nstream;
   int    sor_max_iter, sor_min_iter, mg_max_iter, mg_npre_smooth, mg_npost_smooth, pot_gpu_npgroup;
   int    opt__flu_int_scheme, opt__pot_int_scheme, opt__rho_int_scheme;
   int    opt__gra_int_scheme, opt__ref_flu_int_scheme, opt__ref_pot_int_scheme;
   int    opt__output_total, opt__output_part, opt__output_mode, output_step;
   long   end_step;
   real   lb_wli_max, gamma, minmod_coeff, ep_coeff, elbdm_mass, planck_const, newton_g, sor_omega;
   real   mg_tolerated_error, output_part_x, output_part_y, output_part_z;
   double box_size, end_t, omega_m0, dt__fluid, dt__gravity, dt__phase, dt__max_delta_a, output_dt;

   fseek( File, HeaderOffset_Parameter, SEEK_SET );

   fread( &box_size,                   sizeof(double),                  1,             File );
   fread(  nx0_tot,                    sizeof(int),                     3,             File );
   fread( &gamer_nrank,                sizeof(int),                     1,             File );
   fread(  gamer_nrank_x,              sizeof(int),                     3,             File );
   fread( &omp_nthread,                sizeof(int),                     1,             File );
   fread( &end_t,                      sizeof(double),                  1,             File );
   fread( &end_step,                   sizeof(long),                    1,             File );
   fread( &ooc_nrank,                  sizeof(int),                     1,             File );
   fread(  ooc_nrank_x,                sizeof(int),                     3,             File );
   fread( &omega_m0,                   sizeof(double),                  1,             File );
   fread( &dt__fluid,                  sizeof(double),                  1,             File );
   fread( &dt__gravity,                sizeof(double),                  1,             File );
   fread( &dt__phase,                  sizeof(double),                  1,             File );
   fread( &dt__max_delta_a,            sizeof(double),                  1,             File );
   fread( &opt__adaptive_dt,           sizeof(bool),                    1,             File );
   fread( &opt__dt_user,               sizeof(bool),                    1,             File );
   fread( &regrid_count,               sizeof(int),                     1,             File );
   fread( &flag_buffer_size,           sizeof(int),                     1,             File );
   fread( &max_level,                  sizeof(int),                     1,             File );
   fread( &opt__flag_rho,              sizeof(bool),                    1,             File );
   fread( &opt__flag_rho_gradient,     sizeof(bool),                    1,             File );
   fread( &opt__flag_pres_gradient,    sizeof(bool),                    1,             File );
   fread( &opt__flag_engy_density,     sizeof(bool),                    1,             File );
   fread( &opt__flag_user,             sizeof(bool),                    1,             File );
   fread( &lb_wli_max,                 sizeof(real),                    1,             File );
   fread( &gamma,                      sizeof(real),                    1,             File );
   fread( &minmod_coeff,               sizeof(real),                    1,             File );
   fread( &ep_coeff,                   sizeof(real),                    1,             File );
   fread( &opt__lr_limiter,            sizeof(int),                     1,             File );
   fread( &opt__waf_limiter,           sizeof(int),                     1,             File );
   fread( &elbdm_mass,                 sizeof(real),                    1,             File );
   fread( &planck_const,               sizeof(real),                    1,             File );
   fread( &flu_gpu_npgroup,            sizeof(int),                     1,             File );
   fread( &gpu_nstream,                sizeof(int),                     1,             File );
   fread( &opt__fixup_flux,            sizeof(bool),                    1,             File );
   fread( &opt__fixup_restrict,        sizeof(bool),                    1,             File );
   fread( &opt__overlap_mpi,           sizeof(bool),                    1,             File );
   fread( &newton_g,                   sizeof(real),                    1,             File );
   fread( &sor_omega,                  sizeof(real),                    1,             File );
   fread( &sor_max_iter,               sizeof(int),                     1,             File );
   fread( &sor_min_iter,               sizeof(int),                     1,             File );
   fread( &mg_max_iter,                sizeof(int),                     1,             File );
   fread( &mg_npre_smooth,             sizeof(int),                     1,             File );
   fread( &mg_npost_smooth,            sizeof(int),                     1,             File );
   fread( &mg_tolerated_error,         sizeof(real),                    1,             File );
   fread( &pot_gpu_npgroup,            sizeof(int),                     1,             File );
   fread( &opt__gra_p5_gradient,       sizeof(bool),                    1,             File );
   fread( &opt__int_time,              sizeof(bool),                    1,             File );
   fread( &opt__int_phase,             sizeof(bool),                    1,             File );
   fread( &opt__flu_int_scheme,        sizeof(int),                     1,             File );
   fread( &opt__pot_int_scheme,        sizeof(int),                     1,             File );
   fread( &opt__rho_int_scheme,        sizeof(int),                     1,             File );
   fread( &opt__gra_int_scheme,        sizeof(int),                     1,             File );
   fread( &opt__ref_flu_int_scheme,    sizeof(int),                     1,             File );
   fread( &opt__ref_pot_int_scheme,    sizeof(int),                     1,             File );
   fread( &opt__output_total,          sizeof(int),                     1,             File );
   fread( &opt__output_part,           sizeof(int),                     1,             File );
   fread( &opt__output_error,          sizeof(bool),                    1,             File );
   fread( &opt__output_base,           sizeof(bool),                    1,             File );
   fread( &opt__output_pot,            sizeof(bool),                    1,             File );
   fread( &opt__output_mode,           sizeof(int),                     1,             File );
   fread( &output_step,                sizeof(int),                     1,             File );
   fread( &output_dt,                  sizeof(double),                  1,             File );
   fread( &output_part_x,              sizeof(real),                    1,             File );
   fread( &output_part_y,              sizeof(real),                    1,             File );
   fread( &output_part_z,              sizeof(real),                    1,             File );
   fread( &opt__timing_barrier,        sizeof(bool),                    1,             File );

   if ( MyRank == 0 )   Aux_Message( stdout, "   Loading simulation parameters ... done\n" );


// d. check parameters (before loading any size-dependent parameters)
// =================================================================================================
   if ( MyRank == 0 )   Aux_Message( stdout, "   Checking loaded parameters ...\n" );

   const bool Fatal    = true;
   const bool NonFatal = false;

#  ifdef FLOAT8
   if ( !float8 )
   {
      fprintf( stderr, "ERROR : %s : RESTART file (%s) != runtime (%s) !!\n", "FLOAT8", "OFF", "ON" );
      MPI_Exit();
   }
#  else
   if (  float8 )
   {
      fprintf( stderr, "ERROR : %s : RESTART file (%s) != runtime (%s) !!\n", "FLOAT8", "ON", "OFF" );
      MPI_Exit();
   }
#  endif

   CompareVar( "MODEL",                   model,                  MODEL,                        Fatal );
   CompareVar( "NLEVEL",                  nlevel,                 NLEVEL,                       Fatal );
   CompareVar( "NCOMP_FLUID",             ncomp_fluid,            NCOMP_FLUID,                  Fatal );
   CompareVar( "PATCH_SIZE",              patch_size,             PATCH_SIZE,                   Fatal );

   if ( MyRank == 0 )   Aux_Message( stdout, "   Checking loaded parameters ... done\n" );


// set the returned variables
   if ( FormatVersion >= 1210 )
   DataOrder_xyzv = false;
   else
   DataOrder_xyzv = ( opt__output_total == 1 ) ? true : false;
   LoadPot        = opt__output_pot;
   BoxSize        = box_size;
   Gamma          = gamma;
   for (int d=0; d<3; d++)
   NX0_Tot[d]     = nx0_tot[d];
   ELBDM_Eta      = elbdm_mass / planck_const;
   Comoving       = comoving;

} // FUNCTION : Load_Parameter_Before_2000



//-------------------------------------------------------------------------------------------------------
// Function    :  Load_Parameter_After_2000
// Description :  Load all simulation parameters from the RESTART file with format version >= 2000
//
// Note        :  All floating-point variables are declared as double after version 2000
//
// Parameter   :  File           : RESTART file pointer
//                FormatVersion  : Format version of the RESTART file
//                LoadPot        : Whether or not the RESTART file stores the potential data
//                NX0_Tot        : Total number of base-level cells along each direction
//                BoxSize        : Physical size of the simulation box
//                Gamma          : ratio of specific heat
//                ELBDM_Eta      : mass/planck_const in ELBDM
//                HeaderOffset_X : Offsets of different headers
//                LoadPar        : Whether or not the RESTART file stores the particle data
//                LoadParDens    : Whether or not the RESTART file stores the particle density on grids
//                NParVarOut     : Number of particles attributes stored (for checking the file size only)
//                H0             : Dimensionless Hubble parameter
//                WithUnit       : true --> restart file stores the code units
//                Unit_*         : Code units
//                Mu             : Mean molecular weight
//                Comoving       : true --> cosmological simulations in comoving coords.
//
// Return      :  DataOrder_xyzv, LoadPot, NX0_Tot, BoxSize, Gamma, ELBDM_Eta, LoadPar, LoadParDens, NParVarOut,
//                H0, WithUnit, Unit_*, Mu, Comoving
//-------------------------------------------------------------------------------------------------------
void Load_Parameter_After_2000( FILE *File, const int FormatVersion, bool &LoadPot, int *NX0_Tot,
                                double &BoxSize, double &Gamma, double &ELBDM_Eta,
                                const long HeaderOffset_Makefile, const long HeaderOffset_Constant,
                                const long HeaderOffset_Parameter, bool &LoadPar, int &LoadParDens, int &NParVarOut,
                                double &H0, bool &WithUnit, double &Unit_L, double &Unit_M, double &Unit_T, double &Unit_V,
                                double &Unit_D, double &Unit_E, double &Unit_P, double &Mu, bool &Comoving )
{

   if ( MyRank == 0 )   Aux_Message( stdout, "   Loading simulation parameters ...\n" );


// a. load the simulation options and parameters defined in the Makefile
// =================================================================================================
   bool gravity, individual_timestep, comoving, gpu, gamer_optimization, gamer_debug, timing, timing_solver;
   bool intel, float8, serial, overlap_mpi, openmp, store_pot_ghost, unsplit_gravity, particle;
   bool conserve_mass, laplacian_4th, self_interaction, laohu, support_hdf5;
   int  model, pot_scheme, flu_scheme, lr_scheme, rsolver, load_balance, nlevel, max_patch, ncomp_passive, gpu_arch;

   fseek( File, HeaderOffset_Makefile, SEEK_SET );

   fread( &model,                      sizeof(int),                     1,             File );
   fread( &gravity,                    sizeof(bool),                    1,             File );
   fread( &pot_scheme,                 sizeof(int),                     1,             File );
   fread( &individual_timestep,        sizeof(bool),                    1,             File );
   fread( &comoving,                   sizeof(bool),                    1,             File );
   fread( &flu_scheme,                 sizeof(int),                     1,             File );
   fread( &lr_scheme,                  sizeof(int),                     1,             File );
   fread( &rsolver,                    sizeof(int),                     1,             File );
   fread( &gpu,                        sizeof(bool),                    1,             File );
   fread( &gamer_optimization,         sizeof(bool),                    1,             File );
   fread( &gamer_debug,                sizeof(bool),                    1,             File );
   fread( &timing,                     sizeof(bool),                    1,             File );
   fread( &timing_solver,              sizeof(bool),                    1,             File );
   fread( &intel,                      sizeof(bool),                    1,             File );
   fread( &float8,                     sizeof(bool),                    1,             File );
   fread( &serial,                     sizeof(bool),                    1,             File );
   fread( &load_balance,               sizeof(int),                     1,             File );
   fread( &overlap_mpi,                sizeof(bool),                    1,             File );
   fread( &openmp,                     sizeof(bool),                    1,             File );
   fread( &gpu_arch,                   sizeof(int),                     1,             File );
   fread( &nlevel,                     sizeof(int),                     1,             File );
   fread( &max_patch,                  sizeof(int),                     1,             File );
   fread( &store_pot_ghost,            sizeof(bool),                    1,             File );
   fread( &unsplit_gravity,            sizeof(bool),                    1,             File );
   fread( &particle,                   sizeof(bool),                    1,             File );
   fread( &ncomp_passive,              sizeof(int),                     1,             File );
   fread( &conserve_mass,              sizeof(bool),                    1,             File );
   fread( &laplacian_4th,              sizeof(bool),                    1,             File );
   fread( &self_interaction,           sizeof(bool),                    1,             File );
   fread( &laohu,                      sizeof(bool),                    1,             File );
   fread( &support_hdf5,               sizeof(bool),                    1,             File );

// reset parameters
   if ( FormatVersion < 2100 )   particle = false;


// b. load the symbolic constants defined in "Macro.h, CUPOT.h, and CUFLU.h"
// =================================================================================================
   bool   enforce_positive, char_reconstruction, hll_no_ref_state, hll_include_all_waves, waf_dissipate;
   bool   use_psolver_10to14;
   int    ncomp_fluid, patch_size, flu_ghost_size, pot_ghost_size, gra_ghost_size, check_intermediate;
   int    flu_block_size_x, flu_block_size_y, pot_block_size_x, pot_block_size_z, gra_block_size_z;
   int    par_nvar, par_npassive;
   double min_value, max_error;

   fseek( File, HeaderOffset_Constant, SEEK_SET );

   fread( &ncomp_fluid,                sizeof(int),                     1,             File );
   fread( &patch_size,                 sizeof(int),                     1,             File );
   fread( &min_value,                  sizeof(double),                  1,             File );
   fread( &flu_ghost_size,             sizeof(int),                     1,             File );
   fread( &pot_ghost_size,             sizeof(int),                     1,             File );
   fread( &gra_ghost_size,             sizeof(int),                     1,             File );
   fread( &enforce_positive,           sizeof(bool),                    1,             File );
   fread( &char_reconstruction,        sizeof(bool),                    1,             File );
   fread( &check_intermediate,         sizeof(int),                     1,             File );
   fread( &hll_no_ref_state,           sizeof(bool),                    1,             File );
   fread( &hll_include_all_waves,      sizeof(bool),                    1,             File );
   fread( &waf_dissipate,              sizeof(bool),                    1,             File );
   fread( &max_error,                  sizeof(double),                  1,             File );
   fread( &flu_block_size_x,           sizeof(int),                     1,             File );
   fread( &flu_block_size_y,           sizeof(int),                     1,             File );
   fread( &use_psolver_10to14,         sizeof(bool),                    1,             File );
   fread( &pot_block_size_x,           sizeof(int),                     1,             File );
   fread( &pot_block_size_z,           sizeof(int),                     1,             File );
   fread( &gra_block_size_z,           sizeof(int),                     1,             File );
   fread( &par_nvar,                   sizeof(int),                     1,             File );
   fread( &par_npassive,               sizeof(int),                     1,             File );

// reset parameters
   if ( !particle )
   {
      par_nvar     = WRONG;
      par_npassive = WRONG;
   }


// c. load the simulation parameters recorded in the file "Input__Parameter"
// =================================================================================================
   bool   opt__adaptive_dt, opt__dt_user, opt__flag_rho, opt__flag_rho_gradient, opt__flag_pres_gradient;
   bool   opt__flag_engy_density, opt__flag_user, opt__fixup_flux, opt__fixup_restrict, opt__overlap_mpi;
   bool   opt__gra_p5_gradient, opt__int_time, opt__output_test_error, opt__output_base, opt__output_pot;
   bool   opt__output_baseps, opt__timing_balance, opt__int_phase, opt__corr_unphy, opt__unit;
   int    nx0_tot[3], mpi_nrank, mpi_nrank_x[3], omp_nthread, regrid_count, opt__output_par_dens;
   int    flag_buffer_size, max_level, opt__lr_limiter, opt__waf_limiter, flu_gpu_npgroup, gpu_nstream;
   int    sor_max_iter, sor_min_iter, mg_max_iter, mg_npre_smooth, mg_npost_smooth, pot_gpu_npgroup;
   int    opt__flu_int_scheme, opt__pot_int_scheme, opt__rho_int_scheme;
   int    opt__gra_int_scheme, opt__ref_flu_int_scheme, opt__ref_pot_int_scheme;
   int    opt__output_total, opt__output_part, opt__output_mode, output_step, opt__corr_unphy_scheme;
   long   end_step;
   double lb_wli_max, gamma, minmod_coeff, ep_coeff, elbdm_mass, elbdm_planck_const, newton_g, sor_omega;
   double mg_tolerated_error, output_part_x, output_part_y, output_part_z, molecular_weight;
   double box_size, end_t, omega_m0, dt__fluid, dt__gravity, dt__phase, dt__max_delta_a, output_dt, hubble0;
   double unit_l, unit_m, unit_t, unit_v, unit_d, unit_e, unit_p;

   fseek( File, HeaderOffset_Parameter, SEEK_SET );

   fread( &box_size,                   sizeof(double),                  1,             File );
   fread(  nx0_tot,                    sizeof(int),                     3,             File );
   fread( &mpi_nrank,                  sizeof(int),                     1,             File );
   fread(  mpi_nrank_x,                sizeof(int),                     3,             File );
   fread( &omp_nthread,                sizeof(int),                     1,             File );
   fread( &end_t,                      sizeof(double),                  1,             File );
   fread( &end_step,                   sizeof(long),                    1,             File );
   fread( &omega_m0,                   sizeof(double),                  1,             File );
   fread( &dt__fluid,                  sizeof(double),                  1,             File );
   fread( &dt__gravity,                sizeof(double),                  1,             File );
   fread( &dt__phase,                  sizeof(double),                  1,             File );
   fread( &dt__max_delta_a,            sizeof(double),                  1,             File );
   fread( &opt__adaptive_dt,           sizeof(bool),                    1,             File );
   fread( &opt__dt_user,               sizeof(bool),                    1,             File );
   fread( &regrid_count,               sizeof(int),                     1,             File );
   fread( &flag_buffer_size,           sizeof(int),                     1,             File );
   fread( &max_level,                  sizeof(int),                     1,             File );
   fread( &opt__flag_rho,              sizeof(bool),                    1,             File );
   fread( &opt__flag_rho_gradient,     sizeof(bool),                    1,             File );
   fread( &opt__flag_pres_gradient,    sizeof(bool),                    1,             File );
   fread( &opt__flag_engy_density,     sizeof(bool),                    1,             File );
   fread( &opt__flag_user,             sizeof(bool),                    1,             File );
   fread( &lb_wli_max,                 sizeof(double),                  1,             File );
   fread( &gamma,                      sizeof(double),                  1,             File );
   fread( &minmod_coeff,               sizeof(double),                  1,             File );
   fread( &ep_coeff,                   sizeof(double),                  1,             File );
   fread( &opt__lr_limiter,            sizeof(int),                     1,             File );
   fread( &opt__waf_limiter,           sizeof(int),                     1,             File );
   fread( &elbdm_mass,                 sizeof(double),                  1,             File );
   fread( &elbdm_planck_const,         sizeof(double),                  1,             File );
   fread( &flu_gpu_npgroup,            sizeof(int),                     1,             File );
   fread( &gpu_nstream,                sizeof(int),                     1,             File );
   fread( &opt__fixup_flux,            sizeof(bool),                    1,             File );
   fread( &opt__fixup_restrict,        sizeof(bool),                    1,             File );
   fread( &opt__overlap_mpi,           sizeof(bool),                    1,             File );
   fread( &newton_g,                   sizeof(double),                  1,             File );
   fread( &sor_omega,                  sizeof(double),                  1,             File );
   fread( &sor_max_iter,               sizeof(int),                     1,             File );
   fread( &sor_min_iter,               sizeof(int),                     1,             File );
   fread( &mg_max_iter,                sizeof(int),                     1,             File );
   fread( &mg_npre_smooth,             sizeof(int),                     1,             File );
   fread( &mg_npost_smooth,            sizeof(int),                     1,             File );
   fread( &mg_tolerated_error,         sizeof(double),                  1,             File );
   fread( &pot_gpu_npgroup,            sizeof(int),                     1,             File );
   fread( &opt__gra_p5_gradient,       sizeof(bool),                    1,             File );
   fread( &opt__int_time,              sizeof(bool),                    1,             File );
   fread( &opt__int_phase,             sizeof(bool),                    1,             File );
   fread( &opt__flu_int_scheme,        sizeof(int),                     1,             File );
   fread( &opt__pot_int_scheme,        sizeof(int),                     1,             File );
   fread( &opt__rho_int_scheme,        sizeof(int),                     1,             File );
   fread( &opt__gra_int_scheme,        sizeof(int),                     1,             File );
   fread( &opt__ref_flu_int_scheme,    sizeof(int),                     1,             File );
   fread( &opt__ref_pot_int_scheme,    sizeof(int),                     1,             File );
   fread( &opt__output_total,          sizeof(int),                     1,             File );
   fread( &opt__output_part,           sizeof(int),                     1,             File );
   fread( &opt__output_test_error,     sizeof(bool),                    1,             File );
   fread( &opt__output_base,           sizeof(bool),                    1,             File );
   fread( &opt__output_pot,            sizeof(bool),                    1,             File );
   fread( &opt__output_mode,           sizeof(int),                     1,             File );
   fread( &output_step,                sizeof(int),                     1,             File );
   fread( &output_dt,                  sizeof(double),                  1,             File );
   fread( &output_part_x,              sizeof(double),                  1,             File );
   fread( &output_part_y,              sizeof(double),                  1,             File );
   fread( &output_part_z,              sizeof(double),                  1,             File );
   fread( &opt__timing_balance,        sizeof(bool),                    1,             File );
   fread( &opt__output_baseps,         sizeof(bool),                    1,             File );
   fread( &opt__corr_unphy,            sizeof(bool),                    1,             File );
   fread( &opt__corr_unphy_scheme,     sizeof(int),                     1,             File );
   fread( &opt__output_par_dens,       sizeof(int),                     1,             File );
   fread( &hubble0,                    sizeof(double),                  1,             File );
   fread( &opt__unit,                  sizeof(bool),                    1,             File );
   fread( &unit_l,                     sizeof(double),                  1,             File );
   fread( &unit_m,                     sizeof(double),                  1,             File );
   fread( &unit_t,                     sizeof(double),                  1,             File );
   fread( &unit_v,                     sizeof(double),                  1,             File );
   fread( &unit_d,                     sizeof(double),                  1,             File );
   fread( &unit_e,                     sizeof(double),                  1,             File );
   fread( &unit_p,                     sizeof(double),                  1,             File );
   fread( &molecular_weight,           sizeof(double),                  1,             File );

// reset parameters
   if ( !particle )  opt__output_par_dens = 0;

   if ( FormatVersion < 2111 )
   {
      hubble0 = unit_l = unit_m = unit_t = unit_v = unit_d = unit_e = unit_p = molecular_weight = WRONG;
      opt__unit = false;
   }

   if ( MyRank == 0 )   Aux_Message( stdout, "   Loading simulation parameters ... done\n" );


// d. check parameters (before loading any size-dependent parameters)
// =================================================================================================
   if ( MyRank == 0 )   Aux_Message( stdout, "   Checking loaded parameters ...\n" );


   const bool Fatal    = true;
   const bool NonFatal = false;

#  ifdef FLOAT8
   if ( !float8 )
   {
      fprintf( stderr, "ERROR : %s : RESTART file (%s) != runtime (%s) !!\n", "FLOAT8", "OFF", "ON" );
      exit( 1 );
   }
#  else
   if (  float8 )
   {
      fprintf( stderr, "ERROR : %s : RESTART file (%s) != runtime (%s) !!\n", "FLOAT8", "ON", "OFF" );
      exit( 1 );
   }
#  endif

   CompareVar( "MODEL",                   model,                  MODEL,                        Fatal );
   CompareVar( "NLEVEL",                  nlevel,                 NLEVEL,                       Fatal );
   CompareVar( "NCOMP_FLUID",             ncomp_fluid,            NCOMP_FLUID,                  Fatal );
   CompareVar( "NCOMP_PASSIVE",           ncomp_passive,          NCOMP_PASSIVE,                Fatal );
   CompareVar( "PATCH_SIZE",              patch_size,             PATCH_SIZE,                   Fatal );

   if ( MyRank == 0 )   Aux_Message( stdout, "   Checking loaded parameters ... done\n" );


// set the returned variables
   LoadPot     = opt__output_pot;
   LoadPar     = particle;
   LoadParDens = opt__output_par_dens;
   BoxSize     = box_size;
   Gamma       = gamma;
   ELBDM_Eta   = elbdm_mass / elbdm_planck_const;
   for (int d=0; d<3; d++)
   NX0_Tot[d]  = nx0_tot[d];
   H0          = hubble0;
   WithUnit    = opt__unit;
   Unit_L      = unit_l;
   Unit_M      = unit_m;
   Unit_T      = unit_t;
   Unit_V      = unit_v;
   Unit_D      = unit_d;
   Unit_E      = unit_e;
   Unit_P      = unit_p;
   Mu          = molecular_weight;
   Comoving    = comoving;

   if ( LoadPar )
   {
      if ( FormatVersion > 2131 )   NParVarOut = par_nvar;           // after version 2131, par_nvar = PAR_NATT_STORED
      else                          NParVarOut = 7 + par_npassive;   // mass, position x/y/z, velocity x/y/z, and passive variables
   }
   else
                                    NParVarOut = -1;

} // FUNCTION : Load_Parameter_After_2000


//-------------------------------------------------------------------------------------------------------
// Function    :  CompareVar
// Description :  Compare the input variables
//
// Note        :  This function is overloaded to work with different data types
//
// Parameter   :  VarName     : Name of the targeted variable
//                RestartVar  : Variable loaded from the RESTART file
//                RuntimeVar  : Variable loaded from the Input__Parameter
//                Fatal       : Whether or not the difference between RestartVar and RuntimeVar is fatal
//                              --> true  : terminate the program if the input variables are different
//                                  false : display warning message if the input variables are different
//-------------------------------------------------------------------------------------------------------
void CompareVar( const char *VarName, const bool RestartVar, const bool RuntimeVar, const bool Fatal )
{

   if ( RestartVar != RuntimeVar )
   {
      if ( Fatal )
      {
         fprintf( stderr, "ERROR : %s : RESTART file (%d) != runtime (%d) !!\n",
                  VarName, RestartVar, RuntimeVar );
         MPI_Exit();
      }
      else
         fprintf( stderr, "WARNING : %s : RESTART file (%d) != runtime (%d) !!\n",
                  VarName, RestartVar, RuntimeVar );
   }

} // FUNCTION : CompareVar (bool)



//-------------------------------------------------------------------------------------------------------
// overloaded function for type "int"
//-------------------------------------------------------------------------------------------------------
void CompareVar( const char *VarName, const int RestartVar, const int RuntimeVar, const bool Fatal )
{

   if ( RestartVar != RuntimeVar )
   {
      if ( Fatal )
      {
         fprintf( stderr, "ERROR : %s : RESTART file (%d) != runtime (%d) !!\n",
                  VarName, RestartVar, RuntimeVar );
         MPI_Exit();
      }
      else
         fprintf( stderr, "WARNING : %s : RESTART file (%d) != runtime (%d) !!\n",
                  VarName, RestartVar, RuntimeVar );
   }

} // FUNCTION : CompareVar (int)



//-------------------------------------------------------------------------------------------------------
// overloaded function for type "long"
//-------------------------------------------------------------------------------------------------------
void CompareVar( const char *VarName, const long RestartVar, const long RuntimeVar, const bool Fatal )
{

   if ( RestartVar != RuntimeVar )
   {
      if ( Fatal )
      {
         fprintf( stderr, "ERROR : %s : RESTART file (%ld) != runtime (%ld) !!\n",
                  VarName, RestartVar, RuntimeVar );
         MPI_Exit();
      }
      else
         fprintf( stderr, "WARNING : %s : RESTART file (%ld) != runtime (%ld) !!\n",
                  VarName, RestartVar, RuntimeVar );
   }

} // FUNCTION : CompareVar (long)



//-------------------------------------------------------------------------------------------------------
// overloaded function for type "float"
//-------------------------------------------------------------------------------------------------------
void CompareVar( const char *VarName, const float RestartVar, const float RuntimeVar, const bool Fatal )
{

   if ( RestartVar != RuntimeVar )
   {
      if ( Fatal )
      {
         fprintf( stderr, "ERROR : %s : RESTART file (%20.14e) != runtime (%20.14e) !!\n",
                  VarName, RestartVar, RuntimeVar );
         MPI_Exit();
      }
      else
         fprintf( stderr, "WARNING : %s : RESTART file (%20.14e) != runtime (%20.14e) !!\n",
                  VarName, RestartVar, RuntimeVar );
   }

} // FUNCTION : CompareVar (float)



//-------------------------------------------------------------------------------------------------------
// overloaded function for type "double"
//-------------------------------------------------------------------------------------------------------
void CompareVar( const char *VarName, const double RestartVar, const double RuntimeVar, const bool Fatal )
{

   if ( RestartVar != RuntimeVar )
   {
      if ( Fatal )
      {
         fprintf( stderr, "ERROR : %s : RESTART file (%20.14e) != runtime (%20.14e) !!\n",
                  VarName, RestartVar, RuntimeVar );
         MPI_Exit();
      }
      else
         fprintf( stderr, "WARNING : %s : RESTART file (%20.14e) != runtime (%20.14e) !!\n",
                  VarName, RestartVar, RuntimeVar );
   }

} // FUNCTION : CompareVar (double)



//-------------------------------------------------------------------------------------------------------
// Function    :  LoadOnePatch
// Description :  Load one patch from the GAMER data
//
// Note        :  1. In this function one should specify the criteria for loading a patch
//                2. Only leaf patches will store data. If the target patch is not a leaf patch, this function
//                   will be invoked recursively to find the leaf amr.
//
// Parameter   :  File        : GAMER data file
//                lv          : Target level
//                LoadPID     : Target patch index to load data
//                LoadPot     : True --> load the gravitational potential
//                LoadParDens : True --> load the particle density on grids
//-------------------------------------------------------------------------------------------------------
void LoadOnePatch( FILE *File, const int lv, const int LoadPID, const bool LoadPot, const bool LoadParDens )
{

   const int LoadSonPID0 = tree->block[lv][LoadPID].son;
   int  Cr[3], PID;
   bool GotYou;

// load corner
   for (int d=0; d<3; d++)    Cr[d] = tree->block[lv][LoadPID].corner[d];

// shift the corner scale (usually used when the target region lies outside the simulation box)
   if ( Shift2Center )
   {
      for (int d=0; d<3; d++)    Cr[d] = ( Cr[d] + ShiftScale[d] + amr.BoxScale[d] ) % amr.BoxScale[d];
   }

// verify that the loaded patch is within the candidate box
   GotYou = WithinCandidateBox( Cr, PATCH_SIZE*amr.scale[lv], CanBuf );

// allocate patch
   amr.pnew( lv, Cr[0], Cr[1], Cr[2], -1, GotYou );

// load data
   if ( GotYou  &&  LoadSonPID0 == -1 )
   {
      PID = amr.num[lv] - 1;

//    load data
      fseek( File, tree->block[lv][LoadPID].data_pos, SEEK_SET );

      fread( amr.patch[lv][PID]->fluid,    sizeof(real), NCOMP_TOTAL*BLOCK_SIZE*BLOCK_SIZE*BLOCK_SIZE, File );

      if ( LoadPot )
      fread( amr.patch[lv][PID]->pot,      sizeof(real),             BLOCK_SIZE*BLOCK_SIZE*BLOCK_SIZE, File );

      if ( LoadParDens )
      fread( amr.patch[lv][PID]->par_dens, sizeof(real),             BLOCK_SIZE*BLOCK_SIZE*BLOCK_SIZE, File );
   } // if ( GotYou  &&  LoadSonPID0 == -1 )

// enter the next level
// --> even if this patch lies outside the target domain
// --> in other words, we always allocate ALL patches (but some of them may not have data arrays allocated)
   if ( LoadSonPID0 != -1 )
   {
      for (int LoadSonPID=LoadSonPID0; LoadSonPID<LoadSonPID0+8; LoadSonPID++)
         LoadOnePatch( File, lv+1, LoadSonPID, LoadPot, LoadParDens );
   }

} // FUNCTION : LoadOnePatch


