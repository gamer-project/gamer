#include "ExtractProfile.h"
#ifdef SUPPORT_HDF5
#include "hdf5.h"
#endif

static void Load_Parameter_Before_2000( FILE *File, const int FormatVersion, bool &DataOrder_xyzv,
                                        bool &LoadPot, int *NX0_Tot, double &BoxSize, double &Gamma, double &ELBDM_Eta,
                                        const long HeaderOffset_Makefile, const long HeaderOffset_Constant,
                                        const long HeaderOffset_Parameter );
void Load_Parameter_After_2000( FILE *File, const int FormatVersion, bool &LoadPot, int *NX0_Tot,
                                double &BoxSize, double &Gamma, double &ELBDM_Eta,
                                const long HeaderOffset_Makefile, const long HeaderOffset_Constant,
                                const long HeaderOffset_Parameter, bool &LoadPar, int &LoadParDens, int &NParVarOut );
static void CompareVar( const char *VarName, const bool   RestartVar, const bool   RuntimeVar, const bool Fatal );
static void CompareVar( const char *VarName, const int    RestartVar, const int    RuntimeVar, const bool Fatal );
static void CompareVar( const char *VarName, const long   RestartVar, const long   RuntimeVar, const bool Fatal );
static void CompareVar( const char *VarName, const real   RestartVar, const real   RuntimeVar, const bool Fatal );
static void CompareVar( const char *VarName, const double RestartVar, const double RuntimeVar, const bool Fatal );
static void LoadOnePatchGroup( FILE *File, const int lv, const int PID0, const bool LoadPot, const int LoadParDens,
                               const int CanMax_x1, const int CanMin_x1, const int CanMax_y1, const int CanMin_y1,
                               const int CanMax_z1, const int CanMin_z1, const int CanMax_x2, const int CanMin_x2,
                               const int CanMax_y2, const int CanMin_y2, const int CanMax_z2, const int CanMin_z2 );
bool CheckWithinTargetRegion( const bool Periodic, const int CrL[3], const int CrR[3],
                              const int CanMax_x1, const int CanMin_x1, const int CanMax_y1, const int CanMin_y1,
                              const int CanMax_z1, const int CanMin_z1, const int CanMax_x2, const int CanMin_x2,
                              const int CanMax_y2, const int CanMin_y2, const int CanMax_z2, const int CanMin_z2 );
#ifdef SUPPORT_HDF5
void LoadData_HDF5( const char *FileName );
#endif




//-------------------------------------------------------------------------------------------------------
// Function    :  LoadData
// Description :  Load data from the input file
//
// Note        :  1. Unlike LoadData_HDF5, this function only loads and allocates leaf patches inside the target domain
//                   if NeedGhost is off
//                2. AMR data structure always assumes periodicity (even if Periodic == false)
//                   --> Non-periodicity is taken into account by the variables "Center and  Center_Map"
//                       --> Periodic     BC: Center_Map == Center mapped by periodicity
//                           Non-periodic BC: Center_Map == Center
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


   cout << "Loading data \"" << FileName_In << "\" ..." << endl;


   FILE *File = fopen( FileName_In, "rb" );

   if ( File == NULL )
   {
      fprintf( stderr, "\nERROR : the file \"%s\" does not exist !!\n", FileName_In );
      exit( 1 );
   }


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
      exit( 1 );
   }

   if ( size_long != size_ulong )
   {
      fprintf( stderr, "ERROR : sizeof(long) = %d != sizeof(ulong) = %d !!\n", size_long, size_ulong );
      exit( 1 );
   }



// a. load the information of data format
// =================================================================================================
   long FormatVersion, CheckCode;
   long HeaderSize_Format, HeaderSize_Makefile, HeaderSize_Constant, HeaderSize_Parameter, HeaderSize_SimuInfo;
   long HeaderOffset_Format, HeaderOffset_Makefile, HeaderOffset_Constant, HeaderOffset_Parameter, HeaderOffset_SimuInfo;
   long HeaderSize_Total;


// verify the input data format version
   fread( &FormatVersion, sizeof(long), 1, File );

   fprintf( stdout, "   The version of the RESTART file's format = %ld\n", FormatVersion );

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


// check if the size of different data types are consistent
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
      exit( 1 );
   }

   if ( size_int_restart != size_int )
   {
      fprintf( stderr, "ERROR : sizeof(int) is inconsistent : RESTART file = %d, runtime = %d !!\n",
               size_int_restart, size_int );
      exit( 1 );
   }

   if ( size_long_restart != size_long )
   {
      fprintf( stderr, "ERROR : sizeof(long) is inconsistent : RESTART file = %d, runtime = %d !!\n",
               size_long_restart, size_long );
      exit( 1 );
   }

   if ( size_real_restart != size_real )
   {
      fprintf( stderr, "ERROR : sizeof(real) is inconsistent : RESTART file = %d, runtime = %d !!\n",
               size_real_restart, size_real );
      exit( 1 );
   }

   if ( size_double_restart != size_double )
   {
      fprintf( stderr, "ERROR : sizeof(double) is inconsistent : RESTART file = %d, runtime = %d !!\n",
               size_double_restart, size_double );
      exit( 1 );
   }



// b. load all simulation parameters
// =================================================================================================
   bool   DataOrder_xyzv, LoadPot, LoadPar;
   int    LoadParDens, NParVarOut;
   double BoxSize;
#  if ( MODEL != ELBDM )
   double ELBDM_ETA = NULL_REAL;
#  endif

   if ( FormatVersion < 2000 )
   {
      Load_Parameter_Before_2000( File, FormatVersion, DataOrder_xyzv, LoadPot, NX0_TOT, BoxSize, GAMMA, ELBDM_ETA,
                                  HeaderOffset_Makefile, HeaderOffset_Constant, HeaderOffset_Parameter );
      LoadPar     = false;
      LoadParDens = 0;
      NParVarOut  = -1;
   }

   else
   {
      Load_Parameter_After_2000( File, FormatVersion, LoadPot, NX0_TOT, BoxSize, GAMMA, ELBDM_ETA,
                                 HeaderOffset_Makefile, HeaderOffset_Constant, HeaderOffset_Parameter,
                                 LoadPar, LoadParDens, NParVarOut );
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

   if ( LoadPar )
   fread( &NPar,            sizeof(long),        1, File );


// verify the size of the input file
   long DataSize[NLEVEL], ExpectSize, InputSize, PatchDataSize;

   Aux_Message( stdout, "   Verifying the file size ...\n" );

   NIn = NCOMP_TOTAL;            // NIn is a global variable
   if ( LoadPot     )   NIn ++;
   if ( LoadParDens )   NIn ++;

   PatchDataSize = PATCH_SIZE*PATCH_SIZE*PATCH_SIZE*NIn*sizeof(real);
   ExpectSize    = HeaderSize_Total;

   for (int lv=0; lv<NLEVEL; lv++)
   {
      DataSize[lv]  = 0;
      DataSize[lv] += NPatchTotal[lv]*4*sizeof(int);        // 4 = corner(3) + son(1)
      DataSize[lv] += NDataPatch_Total[lv]*PatchDataSize;

      ExpectSize   += DataSize[lv];
   }

   if ( LoadPar )
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
   {
      fprintf( stderr, "ERROR : the size of the file <%s> is incorrect --> input = %ld <-> expect = %ld !!\n",
               FileName_In, InputSize, ExpectSize );
      exit( 1 );
   }

   Aux_Message( stdout, "   Verifying the file size ... passed\n" );


// set up the simulation box
   int NX0_Max;
   NX0_Max = ( NX0_TOT[0] > NX0_TOT[1] ) ? NX0_TOT[0] : NX0_TOT[1];
   NX0_Max = ( NX0_TOT[2] > NX0_Max    ) ? NX0_TOT[2] : NX0_Max;

   for (int lv=0; lv<NLEVEL; lv++)  amr.dh[lv] = BoxSize / (double)( NX0_Max*(1<<lv) );

   for (int d=0; d<3; d++)
   {
      amr.BoxSize [d] = NX0_TOT[d]*amr.dh   [0];
      amr.BoxScale[d] = NX0_TOT[d]*amr.scale[0];
   }



// d. set up and check parameters dedicated in this tool
// =================================================================================================
   static bool FirstTime = true;

   if ( FirstTime )
   {

//    number of output variables
      NOut = NCOMP_TOTAL;

      if ( LoadPot )       NOut ++;

#     if   ( MODEL == HYDRO )
#     elif ( MODEL == MHD )
#     elif ( MODEL == ELBDM )
      if ( ELBDM_GetVir )  NOut += 8;
#     else
#     error : ERROR : unsupported MODEL !!
#     endif // MODEL

      if ( LoadParDens )   NOut ++;


//    convert physical coordinates to cell scales
      if ( !InputScale )
      {
         const double _dh_min = 1.0 / amr.dh[NLEVEL-1];

         for (int d=0; d<3; d++)
         {
            if ( Center[d] != WRONG )  Center[d] *= _dh_min;
         }

         if ( MaxRadius      != WRONG )   MaxRadius      *= _dh_min;
         if ( UseMaxRhoPos_R != WRONG )   UseMaxRhoPos_R *= _dh_min;
      }


//    initialize the sphere information if they are not provided by users
      if ( NShell == WRONG )
      {
         NShell = ( NX0_TOT[0] < NX0_TOT[1] ) ? NX0_TOT[0] : NX0_TOT[1];
         NShell = ( NX0_TOT[2] < NShell     ) ? NX0_TOT[2] : NShell;
         NShell = NShell*amr.scale[0]/2;
      }

      if ( Center[0] == WRONG )  Center[0] = 0.5*amr.BoxScale[0];
      if ( Center[1] == WRONG )  Center[1] = 0.5*amr.BoxScale[1];
      if ( Center[2] == WRONG )  Center[2] = 0.5*amr.BoxScale[2];

      if ( MaxRadius == WRONG )
      {
         MaxRadius = fmin( 0.5*amr.BoxScale[0], 0.5*amr.BoxScale[1] );
         MaxRadius = fmin( 0.5*amr.BoxScale[2], MaxRadius );
      }

      if ( UseMaxRhoPos_R <= 0.0 )  UseMaxRhoPos_R = MaxRadius;

      FirstTime = false;
   } // if ( FirstTime )


// set up the "mapped" sphere center for the periodic B.C.
   const double HalfBox[3] = { 0.5*amr.BoxScale[0], 0.5*amr.BoxScale[1], 0.5*amr.BoxScale[2] };

   if ( Periodic )
      for (int d=0; d<3; d++)
         Center_Map[d] = ( Center[d] > HalfBox[d] ) ? Center[d]-2.0*HalfBox[d] : Center[d]+2.0*HalfBox[d];

   else
      for (int d=0; d<3; d++)
         Center_Map[d] = Center[d];


// set up the candidate box --> only patches with cells inside the candidate box will be allocated
   const int PScale0   = PATCH_SIZE*amr.scale[0];

   const int CanMax_x1 = Center    [0] + MaxRadius + PScale0;  // extend one more base-level patch to ensure the data of all
   const int CanMin_x1 = Center    [0] - MaxRadius - PScale0;  // sibling patches are also loaded (useful when stencil is
   const int CanMax_y1 = Center    [1] + MaxRadius + PScale0;  // required)
   const int CanMin_y1 = Center    [1] - MaxRadius - PScale0;
   const int CanMax_z1 = Center    [2] + MaxRadius + PScale0;
   const int CanMin_z1 = Center    [2] - MaxRadius - PScale0;

   const int CanMax_x2 = Center_Map[0] + MaxRadius + PScale0;
   const int CanMin_x2 = Center_Map[0] - MaxRadius - PScale0;
   const int CanMax_y2 = Center_Map[1] + MaxRadius + PScale0;
   const int CanMin_y2 = Center_Map[1] - MaxRadius - PScale0;
   const int CanMax_z2 = Center_Map[2] + MaxRadius + PScale0;
   const int CanMin_z2 = Center_Map[2] - MaxRadius - PScale0;



// e. load the simulation data
// =================================================================================================
   int  Cr0[3], Cr[3], CrL[3], CrR[3], PScale, LoadSon, PID;
   bool GotYou, PatchGroupAllocated;

// array for re-ordering the fluid data from "xyzv" to "vxyz"
   real (*InvData_Flu)[PATCH_SIZE][PATCH_SIZE][NCOMP_TOTAL] = NULL;
   if ( DataOrder_xyzv )
   {
      if ( UseTree )    Aux_Error( ERROR_INFO, "UseTree option does not support DataOrder_xyzv !!\n" );

      InvData_Flu = new real [PATCH_SIZE][PATCH_SIZE][PATCH_SIZE][NCOMP_TOTAL];
   }

   fseek( File, HeaderSize_Total, SEEK_SET );

   if ( UseTree )
   {
      Aux_Message( stdout, "   Loading data by walking tree ...\n" );

//    loop over all root-level patch groups
      const int TenPercent = MAX( NPatchTotal[0]/10-NPatchTotal[0]/10%8, 1 );

      for (int LoadPID0=0; LoadPID0<NPatchTotal[0]; LoadPID0+=8)
      {
         if ( LoadPID0%TenPercent == 0 )
         Aux_Message( stdout, "      %5.1f%% completed ...\n", 100.0*(double)LoadPID0/NPatchTotal[0] );

         LoadOnePatchGroup( File, 0, LoadPID0, LoadPot, LoadParDens,
                            CanMax_x1, CanMin_x1, CanMax_y1, CanMin_y1, CanMax_z1, CanMin_z1,
                            CanMax_x2, CanMin_x2, CanMax_y2, CanMin_y2, CanMax_z2, CanMin_z2 );
      }

      Aux_Message( stdout, "   Loading data by walking tree ... done\n" );

      Aux_Message( stdout, "\n" );
      for (int lv=0; lv<NLEVEL; lv++)     Aux_Message( stdout, "      Lv %2d: %10d patches loaded\n", lv, amr.num[lv] );
      Aux_Message( stdout, "\n" );
   }

   else
   for (int lv=0; lv<NLEVEL; lv++)
   {
      cout << "   Loading level " << lv << " ... ";

      PScale = PATCH_SIZE*amr.scale[lv];

      for (int LoadPID0=0; LoadPID0<NPatchTotal[lv]; LoadPID0+=8)
      {
         fread( Cr0, sizeof(int), 3, File );       // load the corner of the patch group
         fseek( File, -3*sizeof(int), SEEK_CUR );  // go back to the original position for the next load

         PatchGroupAllocated = false;

         for (int LoadPID=LoadPID0, t=0; LoadPID<LoadPID0+8; LoadPID++, t++)
         {
//          e1. load the patch information
            fread( Cr,       sizeof(int), 3, File );
            fread( &LoadSon, sizeof(int), 1, File );


//          e2. check whether or not the target patch is within the candidate box
            GotYou = false;

            for (int dim=0; dim<3; dim++)
            {
               CrL[dim] = Cr[dim];
               CrR[dim] = Cr[dim] + PScale;
            }

            GotYou = CheckWithinTargetRegion( Periodic, CrL, CrR, CanMax_x1, CanMin_x1, CanMax_y1, CanMin_y1,
                                              CanMax_z1, CanMin_z1, CanMax_x2, CanMin_x2,
                                              CanMax_y2, CanMin_y2, CanMax_z2, CanMin_z2 );


//          e3. when ghost zone are required:
//          (1) always allocate patch (even when it lies outside the target domain)
//          (2) always allocate data array if it lies within the target domain (even for non-leaf patches)
            if ( NeedGhost )    amr.pnew( lv, CrL[0], CrL[1], CrL[2], -1, GotYou );


//          e4. load the physical data if it is a leaf patch
            if ( LoadSon == -1 )
            {
//             skip particle info
               if ( LoadPar )    fseek( File, 2*sizeof(long), SEEK_CUR );

               if ( GotYou )
               {
//                e4-0. set the target PID (and allocate patches when ghost zones are not required)
                  if ( NeedGhost )
                     PID = amr.num[lv] - 1;

                  else
                  {
                     if ( !PatchGroupAllocated )
                     {
//                      allocate all patches in the same patch group to ensure the functionality of "Prepare_PatchData"
                        amr.pnew( lv, Cr0[0]+0*PScale, Cr0[1]+0*PScale, Cr0[2]+0*PScale, -1, false );
                        amr.pnew( lv, Cr0[0]+1*PScale, Cr0[1]+0*PScale, Cr0[2]+0*PScale, -1, false );
                        amr.pnew( lv, Cr0[0]+0*PScale, Cr0[1]+1*PScale, Cr0[2]+0*PScale, -1, false );
                        amr.pnew( lv, Cr0[0]+0*PScale, Cr0[1]+0*PScale, Cr0[2]+1*PScale, -1, false );
                        amr.pnew( lv, Cr0[0]+1*PScale, Cr0[1]+1*PScale, Cr0[2]+0*PScale, -1, false );
                        amr.pnew( lv, Cr0[0]+0*PScale, Cr0[1]+1*PScale, Cr0[2]+1*PScale, -1, false );
                        amr.pnew( lv, Cr0[0]+1*PScale, Cr0[1]+0*PScale, Cr0[2]+1*PScale, -1, false );
                        amr.pnew( lv, Cr0[0]+1*PScale, Cr0[1]+1*PScale, Cr0[2]+1*PScale, -1, false );

                        PatchGroupAllocated = true;
                     }

                     PID = amr.num[lv] - 8 + t;
                     amr.patch[lv][PID]->dnew();
                  }

//                e4-1. load the fluid variables
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
                     fread( amr.patch[lv][PID]->fluid,    sizeof(real), PATCH_SIZE*PATCH_SIZE*PATCH_SIZE*NCOMP_TOTAL, File );

//                e4-2. load the gravitational potential
                  if ( LoadPot )
                     fread( amr.patch[lv][PID]->pot,      sizeof(real), PATCH_SIZE*PATCH_SIZE*PATCH_SIZE,             File );

//                e4-3. load the particle density on grids
                  if ( LoadParDens )
                     fread( amr.patch[lv][PID]->par_dens, sizeof(real), PATCH_SIZE*PATCH_SIZE*PATCH_SIZE,             File );

               } // if ( GotYou )

               else
                  fseek( File, PatchDataSize, SEEK_CUR );

            } // if ( LoadSon == -1 )
         } // for (int LoadPID=LoadPID0; LoadPID<LoadPID0+8; LoadPID++)
      } // for (int LoadPID0=0; LoadPID0<NPatchTotal[lv]; LoadPID0+=8)

      Aux_Message( stdout, "done (%10d patches loaded)\n", amr.num[lv] );

#     ifdef GAMER_DEBUG
      if ( amr.num[lv]%8 != 0 )   Aux_Error( ERROR_INFO, "amr.num[%d] % 8 = (%d) != 0 !!\n", lv, amr.num[lv]%8 );
#     endif

   } // for (int lv=0; lv<NLEVEL; lv++)


   if ( InvData_Flu != NULL )    delete [] InvData_Flu;

   fclose( File );



// f. construct the octree structure when ghost zones are required
//=====================================================================================
   if ( NeedGhost )
   {
      for (int lv=0; lv<NLEVEL; lv++)
      {
//       set up the BaseP List
         if ( lv == 0 )    Init_RecordBasePatch();

//       construct the relation : father <-> son
         else              FindFather( lv );

//       construct the sibling relation
         SiblingSearch( lv );
      }

//    fill up data for non-leaf patches
      for (int lv=NLEVEL-2; lv>=0; lv--)  Flu_Restrict( lv, LoadPot, LoadParDens );
   }



// g. set the shell width and other parameters
// =================================================================================================
   if ( LogBin > 1.0 )
   {
      ShellWidth = GetMinShellWidth( Center, Center_Map );

      if ( GetNShell > 0.0 )  ShellWidth *= GetNShell;   // minimum shell width in the log bin

      NShell = int( log(MaxRadius/ShellWidth)/log(LogBin) ) + 2;
   }

   else  // linear bin
   {
      if ( GetNShell > 0.0 )
      {
         ShellWidth = GetNShell * GetMinShellWidth( Center, Center_Map );
         NShell     = (int)ceil( MaxRadius / ShellWidth );
      }

      else
         ShellWidth = MaxRadius / (double)NShell;
   }

// always output additional grid data if they exist
   OutputPot     = LoadPot;
   OutputParDens = LoadParDens;


   cout << "Loading data \"" << FileName_In << "\" ... done" << endl;

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
//                ELBDM_Eta      : mass/planck_const in ELBDM
//                HeaderOffset_X : Offsets of different headers
//
// Return      :  DataOrder_xyzv, LoadPot, NX0_Tot, BoxSize, Gamma, ELBDM_Eta
//-------------------------------------------------------------------------------------------------------
void Load_Parameter_Before_2000( FILE *File, const int FormatVersion, bool &DataOrder_xyzv,
                                 bool &LoadPot, int *NX0_Tot, double &BoxSize, double &Gamma, double &ELBDM_Eta,
                                 const long HeaderOffset_Makefile, const long HeaderOffset_Constant,
                                 const long HeaderOffset_Parameter )
{

   fprintf( stdout, "   Loading simulation parameters ...\n" );


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


   fprintf( stdout, "   Loading simulation parameters ... done\n" );


// d. check parameters (before loading any size-dependent parameters)
// =================================================================================================
   fprintf( stdout, "   Checking loaded parameters ...\n" );


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
   CompareVar( "PATCH_SIZE",              patch_size,             PATCH_SIZE,                   Fatal );


   fprintf( stdout, "   Checking loaded parameters ... done\n" );


// set the returned variables
   if ( FormatVersion >= 1210 )
   DataOrder_xyzv = false;
   else
   DataOrder_xyzv = ( opt__output_total == 1 ) ? true : false;
   LoadPot        = opt__output_pot;
   BoxSize        = box_size;
   Gamma          = gamma;
   ELBDM_Eta      = elbdm_mass / planck_const;
   for (int d=0; d<3; d++)
   NX0_Tot[d]     = nx0_tot[d];

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
//                Gamma          : Ratio of specific heat
//                ELBDM_Eta      : Mass/Planck_const in ELBDM
//                HeaderOffset_X : Offsets of different headers
//                LoadPar        : Whether or not the RESTART file stores the particle data
//                LoadParDens    : Whether or not the RESTART file stores the particle density data on grids
//                NParVarOut     : Number of particle variables stored in the file (for check file size only)
//
// Return      :  DataOrder_xyzv, LoadPot, NX0_Tot, BoxSize, Gamma, ELBDM_Eta, LoadPar, LoadParDens, NParVarOut
//-------------------------------------------------------------------------------------------------------
void Load_Parameter_After_2000( FILE *File, const int FormatVersion, bool &LoadPot, int *NX0_Tot,
                                double &BoxSize, double &Gamma, double &ELBDM_Eta,
                                const long HeaderOffset_Makefile, const long HeaderOffset_Constant,
                                const long HeaderOffset_Parameter, bool &LoadPar, int &LoadParDens,
                                int &NParVarOut )
{

   Aux_Message( stdout, "   Loading simulation parameters ...\n" );


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

   if ( FormatVersion >= 2100 )
   fread( &particle,                   sizeof(bool),                    1,             File );
   else
   particle = false;

   fread( &ncomp_passive,              sizeof(int),                     1,             File );
   fread( &conserve_mass,              sizeof(bool),                    1,             File );
   fread( &laplacian_4th,              sizeof(bool),                    1,             File );
   fread( &self_interaction,           sizeof(bool),                    1,             File );
   fread( &laohu,                      sizeof(bool),                    1,             File );
   fread( &support_hdf5,               sizeof(bool),                    1,             File );


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

   if ( particle ) {
   fread( &par_nvar,                   sizeof(int),                     1,             File );
   fread( &par_npassive,               sizeof(int),                     1,             File ); }


// c. load the simulation parameters recorded in the file "Input__Parameter"
// =================================================================================================
   bool   opt__adaptive_dt, opt__dt_user, opt__flag_rho, opt__flag_rho_gradient, opt__flag_pres_gradient;
   bool   opt__flag_engy_density, opt__flag_user, opt__fixup_flux, opt__fixup_restrict, opt__overlap_mpi;
   bool   opt__gra_p5_gradient, opt__int_time, opt__output_test_error, opt__output_base, opt__output_pot;
   bool   opt__output_baseps, opt__timing_balance, opt__int_phase, opt__corr_unphy;
   int    nx0_tot[3], mpi_nrank, mpi_nrank_x[3], omp_nthread, regrid_count, opt__output_par_dens;
   int    flag_buffer_size, max_level, opt__lr_limiter, opt__waf_limiter, flu_gpu_npgroup, gpu_nstream;
   int    sor_max_iter, sor_min_iter, mg_max_iter, mg_npre_smooth, mg_npost_smooth, pot_gpu_npgroup;
   int    opt__flu_int_scheme, opt__pot_int_scheme, opt__rho_int_scheme;
   int    opt__gra_int_scheme, opt__ref_flu_int_scheme, opt__ref_pot_int_scheme;
   int    opt__output_total, opt__output_part, opt__output_mode, output_step, opt__corr_unphy_scheme;
   long   end_step;
   double lb_wli_max, gamma, minmod_coeff, ep_coeff, elbdm_mass, elbdm_planck_const, newton_g, sor_omega;
   double mg_tolerated_error, output_part_x, output_part_y, output_part_z;
   double box_size, end_t, omega_m0, dt__fluid, dt__gravity, dt__phase, dt__max_delta_a, output_dt;

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

   if ( particle )
   fread( &opt__output_par_dens,       sizeof(int),                     1,             File );
   else
   opt__output_par_dens = 0;


// d. check parameters (before loading any size-dependent parameters)
// =================================================================================================
   Aux_Message( stdout, "   Checking loaded parameters ...\n" );


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


   Aux_Message( stdout, "   Checking loaded parameters ... done\n" );


// set the returned variables
   LoadPot     = opt__output_pot;
   LoadPar     = particle;
   LoadParDens = opt__output_par_dens;
   BoxSize     = box_size;
   Gamma       = gamma;
   ELBDM_Eta   = elbdm_mass / elbdm_planck_const;
   for (int d=0; d<3; d++)
   NX0_Tot[d]  = nx0_tot[d];

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
         exit( 1 );
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
         exit( 1 );
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
         exit( 1 );
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
         exit( 1 );
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
         exit( 1 );
      }
      else
         fprintf( stderr, "WARNING : %s : RESTART file (%20.14e) != runtime (%20.14e) !!\n",
                  VarName, RestartVar, RuntimeVar );
   }

} // FUNCTION : CompareVar (double)



//-------------------------------------------------------------------------------------------------------
// Function    :  LoadOnePatchGroup
// Description :  Load one patch group from the GAMER data
//
// Note        :  1. In this function one should specify the criteria for loading a patch group
//                2. Only leaf patches will store data. If the target patch is not a leaf patch, this function
//                   will be invoked recursively to find the leaf patch.
//
// Parameter   :  File        : GAMER data file
//                lv          : Target level
//                LoadPID0    : Target patch group index to load data
//                LoadPot     : True --> load the gravitational potential
//                LoadParDens : > 0 --> load the particle density on grids
//                CanMax/Min1 : Maximum/Minimum scale of the target region
//                CanMax/Min2 : Maximum/Minimum scale of the target region (for peroidic BC. only)
//-------------------------------------------------------------------------------------------------------
void LoadOnePatchGroup( FILE *File, const int lv, const int LoadPID0, const bool LoadPot, const int LoadParDens,
                        const int CanMax_x1, const int CanMin_x1, const int CanMax_y1, const int CanMin_y1,
                        const int CanMax_z1, const int CanMin_z1, const int CanMax_x2, const int CanMin_x2,
                        const int CanMax_y2, const int CanMin_y2, const int CanMax_z2, const int CanMin_z2 )
{

   const int PScale = PATCH_SIZE*amr.scale[lv];
   const int *Cr0   = tree->block[lv][LoadPID0].corner;

   bool PatchGroupAllocated = false;
   int *CrL, CrR[3], PID, LoadSonPID0;
   bool GotYou;

   for (int LoadPID=LoadPID0, t=0; LoadPID<LoadPID0+8; LoadPID++, t++)
   {
      CrL = tree->block[lv][LoadPID].corner;
      for (int d=0; d<3; d++)    CrR[d] = CrL[d] + PScale;

      GotYou = CheckWithinTargetRegion( Periodic, CrL, CrR, CanMax_x1, CanMin_x1, CanMax_y1, CanMin_y1,
                                        CanMax_z1, CanMin_z1, CanMax_x2, CanMin_x2,
                                        CanMax_y2, CanMin_y2, CanMax_z2, CanMin_z2 );

      if ( GotYou )
      {
         LoadSonPID0 = tree->block[lv][LoadPID].son;

//       load data only if it's a leaf block
         if ( LoadSonPID0 == -1 )
         {
//          allocate the entire patch group to ensure the functionality of "Prepare_PatchData"
            if ( !PatchGroupAllocated )
            {
               amr.pnew( lv, Cr0[0]+0*PScale, Cr0[1]+0*PScale, Cr0[2]+0*PScale, -1, false );
               amr.pnew( lv, Cr0[0]+1*PScale, Cr0[1]+0*PScale, Cr0[2]+0*PScale, -1, false );
               amr.pnew( lv, Cr0[0]+0*PScale, Cr0[1]+1*PScale, Cr0[2]+0*PScale, -1, false );
               amr.pnew( lv, Cr0[0]+0*PScale, Cr0[1]+0*PScale, Cr0[2]+1*PScale, -1, false );
               amr.pnew( lv, Cr0[0]+1*PScale, Cr0[1]+1*PScale, Cr0[2]+0*PScale, -1, false );
               amr.pnew( lv, Cr0[0]+0*PScale, Cr0[1]+1*PScale, Cr0[2]+1*PScale, -1, false );
               amr.pnew( lv, Cr0[0]+1*PScale, Cr0[1]+0*PScale, Cr0[2]+1*PScale, -1, false );
               amr.pnew( lv, Cr0[0]+1*PScale, Cr0[1]+1*PScale, Cr0[2]+1*PScale, -1, false );

               PatchGroupAllocated = true;
            }

            PID = amr.num[lv] - 8 + t;
            amr.patch[lv][PID]->dnew();

//          load data
            fseek( File, tree->block[lv][LoadPID].data_pos, SEEK_SET );
            fread( amr.patch[lv][PID]->fluid,    sizeof(real), NCOMP_TOTAL*BLOCK_SIZE*BLOCK_SIZE*BLOCK_SIZE, File );

            if ( LoadPot )
            fread( amr.patch[lv][PID]->pot,      sizeof(real),             BLOCK_SIZE*BLOCK_SIZE*BLOCK_SIZE, File );

            if ( LoadParDens )
            fread( amr.patch[lv][PID]->par_dens, sizeof(real),             BLOCK_SIZE*BLOCK_SIZE*BLOCK_SIZE, File );
         } // if ( LoadSonPID0 == -1 )

//       not a leaf block --> go to the next level
         else
         {
            LoadOnePatchGroup( File, lv+1, LoadSonPID0, LoadPot, LoadParDens,
                               CanMax_x1, CanMin_x1, CanMax_y1, CanMin_y1, CanMax_z1, CanMin_z1,
                               CanMax_x2, CanMin_x2, CanMax_y2, CanMin_y2, CanMax_z2, CanMin_z2 );
         } // if ( LoadSonPID0 == -1 ) ... else ...
      } // if ( GotYou )
   } // for (int LoadPID=LoadPID0, t=0; LoadPID<LoadPID0+8; LoadPID++, t++)

} // FUNCTION : LoadOnePatchGroup



//-------------------------------------------------------------------------------------------------------
// Function    :  CheckWithinTargetRegion
// Description :  Check whether the input corners are within the target region
//
// Note        :  Work with both periodic and non-periodic BC.
//
// Parameter   :  Periodic    : True/False <--> Periodic/Non-periodic
//                CrL         : Left corner of the input range
//                CrR         : Right corner of the input range
//                CanMax/Min1 : Maximum/Minimum scale of the target region
//                CanMax/Min2 : Maximum/Minimum scale of the target region (for peroidic BC. only)
//
// Return      : True --> the input corners are WITHIN the target region
//-------------------------------------------------------------------------------------------------------
bool CheckWithinTargetRegion( const bool Periodic, const int CrL[3], const int CrR[3],
                              const int CanMax_x1, const int CanMin_x1, const int CanMax_y1, const int CanMin_y1,
                              const int CanMax_z1, const int CanMin_z1, const int CanMax_x2, const int CanMin_x2,
                              const int CanMax_y2, const int CanMin_y2, const int CanMax_z2, const int CanMin_z2 )
{

   bool Within = false;

   if ( Periodic )
   {
      if (  ( (CrL[0] <= CanMax_x1 && CrR[0] >= CanMin_x1) || (CrL[0] <= CanMax_x2 && CrR[0] >= CanMin_x2) )  &&
            ( (CrL[1] <= CanMax_y1 && CrR[1] >= CanMin_y1) || (CrL[1] <= CanMax_y2 && CrR[1] >= CanMin_y2) )  &&
            ( (CrL[2] <= CanMax_z1 && CrR[2] >= CanMin_z1) || (CrL[2] <= CanMax_z2 && CrR[2] >= CanMin_z2) )    )
         Within = true;
   }

   else
   {
      if (   CrL[0] <= CanMax_x1  &&  CrR[0] >= CanMin_x1  &&
             CrL[1] <= CanMax_y1  &&  CrR[1] >= CanMin_y1  &&
             CrL[2] <= CanMax_z1  &&  CrR[2] >= CanMin_z1    )
        Within = true;
   }

   return Within;

} // FUNCTION : CheckWithinTargetRegion
