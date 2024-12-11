#include "ExtractUniform.h"
#ifdef SUPPORT_HDF5
#include "HDF5_Typedef.h"
#endif

int         OutputXYZ         = WRONG;                      // output mode (x,y,z,x-proj,y-proj,z-proj,3D)
char       *FileName_In       = NULL;                       // name of the input file
char       *FileName_Tree     = NULL;                       // name of the tree file
char       *Suffix            = NULL;                       // suffix attached to the output file name
int         OutputFormat      = 2;                          // output data format: 1=text, 2=HDF5, 3=C-binary
bool        OutputPot         = false;                      // true --> output gravitational potential
int         OutputParDens     = 0;                          // 0=off, 1=output particle density, 2=output total density
bool        OutputPres        = false;                      // true --> output pressure
bool        OutputTemp        = false;                      // true --> output temperature
bool        OutputDivVel      = false;                      // true --> output divergence( velocity )
bool        OutputCurlVel     = false;                      // true --> output curl( velocity ) = vorticity
bool        OutputELBDM_Vel   = false;                      // true --> output velocity field in ELBDM
bool        ConvertCom2PhyV   = false;                      // convert peculiar velocity from comoving to physical coords (in km/s)
bool        InputScale        = false;                      // true --> input cell scales instead of coordinates
bool        Shift2Center      = false;                      // true --> shift target region to the box center (assume periodic)
int         ShiftScale[3]     = { 0, 0, 0 };                // scale shifted by Shift2Center
int         ExtBC             = 0;                          // external boundary condition (0/1:periodic/outflow)
int         OMP_NThread       = -1;                         // number of OpenMP threads
int         TargetLevel       = 0;                          // targeted level
int         Scale_Start   [3] = { WRONG, WRONG, WRONG };    // targeted x, y, z starting scale
int         Scale_Size    [3] = { WRONG, WRONG, WRONG };    // targeted x, y, z scale size
int         Idx_Start     [3] = { WRONG, WRONG, WRONG };    // targeted x, y, z array indices
int         Idx_Size      [3] = { WRONG, WRONG, WRONG };    // targeted x, y, z array size
int         Idx_MyStart   [3] = { WRONG, WRONG, WRONG };    // starting x, y, z array indices of this process
int         Idx_MySize    [3] = { WRONG, WRONG, WRONG };    // array size of this process
double      PhyCoord_Start[3] = { WRONG, WRONG, WRONG };    // starting physical coordinates (may be affected by Shift2Center)
double      PhyCoord_Size [3] = { WRONG, WRONG, WRONG };    // targeted size in physical coordinates
double      OutCoord_Start[3] = { WRONG, WRONG, WRONG };    // real starting physical coordinates independent of Shift2Center
int         NGPU_X[3]         = { 1, 1, 1 };                // number of MPI ranks in each direction
int         CanBuf            = WRONG;                      // buffer size for the candidate box
int         NLoad             = NCOMP_TOTAL;                // number of variables loaded from the input file
int         NOut              = NCOMP_TOTAL;                // number of variables to be outputted
real       *OutputArray       = NULL;                       // array storing the output data
int         BufSize           = WRONG;                      // buffer size of the prepapred patch data
double      Int_MonoCoeff     = 2.0;                        // coefficient for the interpolation monotonicity (1<=coeff<=4)
double      Convert2Temp      = WRONG;                      // coefficient for converting pressure/density to temperature
                                                            // --> velocity_code_unit^2*hydrogen_mass*mean_atomic_weight/Boltzmann_constant
IntScheme_t IntScheme         = INT_DEFAULT;                // interpolation scheme
bool        UseTree           = false;                      // true --> use the tree file to improve the I/O performance
tree_t     *tree              = NULL;                       // tree structure when UseTree is on
int         CanMin[3][3], CanMax[3][3];                     // range of the candidate box for loading patch data

AMR_t     amr;
ParaVar_t ParaVar;
double    GAMMA, MU, H0;
int       NX0_TOT[3], NX0[3], MyRank, MyRank_X[3], SibRank[26], NGPU, DumpID;
int      *BaseP, *BounP_IDMap[NLEVEL][26], SendP_NList[NLEVEL][26], *SendP_IDList[NLEVEL][26];
int       RecvP_NList[NLEVEL][26], *RecvP_IDList[NLEVEL][26], NPatchComma[NLEVEL][28];
double    Time[NLEVEL];
long      Step;
bool      WithUnit, Comoving;
double    Unit_L, Unit_M, Unit_T, Unit_V, Unit_D, Unit_E, Unit_P;

#if ( MODEL == ELBDM )
double    ELBDM_ETA;
double    ELBDM_MASS;
bool      ELBDM_IntPhase = true;
#endif




//-------------------------------------------------------------------------------------------------------
// Function    :  ReadOption
// Description :  Read the command-line options
//-------------------------------------------------------------------------------------------------------
void ReadOption( int argc, char **argv )
{

   if ( MyRank == 0 )   Aux_Message( stdout, "%s ...\n", __FUNCTION__ );


   double Temp_Start[3] = { WRONG, WRONG, WRONG };
   double Temp_Size [3] = { WRONG, WRONG, WRONG };
   int c;

   while ( (c = getopt(argc, argv, "hPUdvVTscwki:o:l:n:x:y:z:X:Y:Z:p:q:r:I:m:t:e:B:u:f:")) != -1 )
   {
      switch ( c )
      {
         case 'i': FileName_In     = optarg;
                   break;
         case 'e': FileName_Tree   = optarg;
                   break;
         case 'o': Suffix          = optarg;
                   break;
         case 'l': TargetLevel     = atoi(optarg);
                   break;
         case 'B': ExtBC           = atoi(optarg);
                   break;
         case 'n': OutputXYZ       = atoi(optarg);
                   break;
         case 'x': Temp_Start[0]   = atof(optarg);
                   break;
         case 'y': Temp_Start[1]   = atof(optarg);
                   break;
         case 'z': Temp_Start[2]   = atof(optarg);
                   break;
         case 'X': Temp_Size[0]    = atof(optarg);
                   break;
         case 'Y': Temp_Size[1]    = atof(optarg);
                   break;
         case 'Z': Temp_Size[2]    = atof(optarg);
                   break;
#        ifndef SERIAL
         case 'p': NGPU_X[0]       = atoi(optarg);
                   break;
         case 'q': NGPU_X[1]       = atoi(optarg);
                   break;
         case 'r': NGPU_X[2]       = atoi(optarg);
                   break;
#        endif
         case 'I': IntScheme       = (IntScheme_t)atoi(optarg);
                   break;
         case 'm': Int_MonoCoeff   = atof(optarg);
                   break;
         case 't': OMP_NThread     = atoi(optarg);
                   break;
         case 'f': OutputFormat    = atoi(optarg);
                   break;
         case 'P': OutputPres      = true;
                   break;
         case 'U': OutputTemp      = true;
                   break;
         case 'd': OutputDivVel    = true;
                   break;
         case 'v': OutputCurlVel   = true;
                   break;
         case 's': InputScale      = true;
                   break;
         case 'c': Shift2Center    = true;
                   break;
         case 'T': UseTree         = true;
                   break;
         case 'k': ConvertCom2PhyV = true;
                   break;
         case 'u': Convert2Temp    = atof(optarg);
                   break;
#        if ( MODEL == ELBDM )
         case 'V': OutputELBDM_Vel = true;
                   break;
         case 'w': ELBDM_IntPhase  = false;
                   break;
#        endif
         case 'h':
         case '?': cerr << endl << "usage: " << argv[0]
                        << " [-h (for help)] [-i input fileName] [-o suffix to the output file [none]]"
                        << endl << "                             "
                        << " [-n output mode (1~7 : X-slice, Y-slice, Z-slice, X-proj, Y-proj, Z-proj, 3D)"
                        << endl << "                             "
                        << " [-x/y/z starting coordinate in x/y/z [0]] [-X/Y/Z target size in x/y/z [BoxSize]]"
#                       ifndef SERIAL
                        << endl << "                             "
                        << " [-p/q/r # of ranks in the x/y/z directions [1,1,1]]"
#                       endif
#                       if   ( MODEL == HYDRO  ||  MODEL == MHD )
                        << endl << "                             "
                        << " [-d (output div(vel)) [off]] [-v (output curl(vel)) [off]] [-P (output pressure) [off]]"
                        << endl << "                             "
                        << " [-U (output temperature) [off]] "
                        << endl << "                             "
                        << " [-u coefficient for converting \"pressure over density\" to temperature"
                        << endl << "                             "
                        << "     = velocity_code_unit^2*hydrogen_mass*mean_atomic_weight/Boltzmann_constant]"
#                       elif ( MODEL == ELBDM )
                        << endl << "                             "
                        << " [-V (output ELBDM vel) [off]] [-w (do not interpolate on phase) [interpolate on phase]]"
#                       endif
                        << endl << "                             "
                        << " [-B external boundary condition (0/1: periodic/outflow) [0]]"
                        << endl << "                             "
                        << " [-k (convert peculiar velocity from comoving to physical coords in km/s) [false]]"
                        << endl << "                             "
                        << " [-f output data format (1=text, 2=HDF5, 3=C-binary) [2]] [-l targeted level [0]]"
                        << endl << "                             "
                        << " [-s [input cell scales instead of physical coordinates to specify the range] [off]]"
                        << endl << "                             "
                        << " [-t number of OpenMP threads [omp_get_max_threads]]"
                        << endl << "                             "
                        << " [-c shift the target region to the box center assuming periodicity [off]]"
                        << endl << "                             "
#                       if ( MODEL == ELBDM )
                        << " [-I interpolation scheme (1,2,3,4,5,6,7->MinMod-3D,MinMod-1D,vanLeer,CQuad,Quad,CQuar,Quar) [6]]"
#                       else
                        << " [-I interpolation scheme (1,2,3,4,5,6,7->MinMod-3D,MinMod-1D,vanLeer,CQuad,Quad,CQuar,Quar) [4]]"
#                       endif
                        << endl << "                             "
                        << " [-m coefficient for monotonic interpolation (1->4) [2.0])"
                        << endl << "                             "
                        << " [-T (load the tree file) [off]] [-e tree file]"
                        << endl << endl;
                   exit( 1 );

      } // switch ( c ) ...
   } // while ...


// reset the useless options
#  if ( MODEL != HYDRO  &&  MODEL != MHD )
   if ( OutputPres )
   {
      OutputPres = false;

      if ( MyRank == 0 )
      fprintf( stderr, "WARNING : option -P (OutputPres) is useless in this MODEL (%d) !!\n", MODEL );
   }

   if ( OutputTemp )
   {
      OutputTemp = false;

      if ( MyRank == 0 )
      fprintf( stderr, "WARNING : option -U (OutputTemp) is useless in this MODEL (%d) !!\n", MODEL );
   }

   if ( OutputDivVel )
   {
      OutputDivVel = false;

      if ( MyRank == 0 )
      fprintf( stderr, "WARNING : option -d (OutputDivVel) is useless in this MODEL (%d) !!\n", MODEL );
   }

   if ( OutputCurlVel )
   {
      OutputCurlVel = false;

      if ( MyRank == 0 )
      fprintf( stderr, "WARNING : option -v (OutputCurlVel) is useless in this MODEL (%d) !!\n", MODEL );
   }
#  endif // #if ( MODEL != HYDRO  &&  MODEL != MHD )

#  if   ( MODEL == ELBDM )
   if ( !OutputELBDM_Vel  &&  ConvertCom2PhyV )
   {
      ConvertCom2PhyV = false;

      if ( MyRank == 0 )
      fprintf( stderr, "WARNING : option -k (ConvertCom2PhyV) is useless when option -V is off !!\n" );
   }

#  elif ( MODEL == HYDRO )
   if ( !OutputDivVel  &&  !OutputCurlVel  &&  ConvertCom2PhyV )
   {
      ConvertCom2PhyV = false;

      if ( MyRank == 0 )
      fprintf( stderr, "WARNING : option -k (ConvertCom2PhyV) is useless when both options -d and -v are off !!\n" );
   }

#  else
   if ( ConvertCom2PhyV )
   {
      ConvertCom2PhyV = false;

      if ( MyRank == 0 )
      fprintf( stderr, "WARNING : option -k (ConvertCom2PhyV) is useless in this MODEL (%d) !!\n", MODEL );
   }

#  endif // MODEL

#  if ( MODEL != ELBDM )
   if ( OutputELBDM_Vel )
   {
      OutputELBDM_Vel = false;

      if ( MyRank == 0 )
      fprintf( stderr, "WARNING : option -V (OutputELBDM_Vel) is useless in this MODEL (%d) !!\n", MODEL );
   }
#  endif

#  ifdef SERIAL
   NGPU_X[0] = 1;
   NGPU_X[1] = 1;
   NGPU_X[2] = 1;
#  endif


// default interpolation scheme
   if ( IntScheme == INT_DEFAULT )
   {
#     if ( MODEL == ELBDM )
      IntScheme = INT_CQUAR;
#     else
      IntScheme = INT_CQUAD;
#     endif
   }


// buffer size
   int NSide, NGhost;
   Int_Table( IntScheme, NSide, NGhost );
   BufSize = NGhost*2;
// BufSize = 8;


// set up OpenMP
#  ifdef OPENMP
   const int OMP_Max_NThread = omp_get_max_threads();

   if ( OMP_NThread <= 0 )
   {
      OMP_NThread = OMP_Max_NThread;

      if ( MyRank == 0 )  fprintf( stdout, "NOTE : parameter \"%s\" is set to the default value = %d\n",
                                   "OMP_NThread", OMP_NThread );
   }

   else if ( OMP_NThread > OMP_Max_NThread   &&  MyRank == 0 )
      fprintf( stderr, "WARNING : OMP_NThread (%d) > omp_get_max_threads (%d) !!\n", OMP_NThread, OMP_Max_NThread );

   omp_set_num_threads( OMP_NThread );
   omp_set_nested( false );

#  else
   OMP_NThread = 1;
#  endif // #ifdef OPENMP ... else ...


// set target range
   if ( InputScale )
   {
      for (int d=0; d<3; d++)
      {
         Scale_Start[d] = (int)Temp_Start[d];
         Scale_Size [d] = (int)Temp_Size [d];
      }
   }
   else
   {
      for (int d=0; d<3; d++)
      {
         PhyCoord_Start[d] = (double)Temp_Start[d];
         PhyCoord_Size [d] = (double)Temp_Size [d];
      }
   }

   NGPU = NGPU_X[0]*NGPU_X[1]*NGPU_X[2];

   if ( OutputPres      )  NOut ++;
   if ( OutputTemp      )  NOut ++;
   if ( OutputDivVel    )  NOut ++;
   if ( OutputCurlVel   )  NOut +=3;
   if ( OutputELBDM_Vel )  NOut +=3;


// check the input command-line options
   if ( FileName_In == NULL )
   {
      cerr << "ERROR : please provide the name of the input file (-i FileName) !!" << endl;
      exit( 1 );
   }

   if ( UseTree && fopen( FileName_Tree, "rb" ) == NULL )
   {
      fprintf( stderr, "ERROR : the input tree file \"%s\" does not exist (-e tree file) !!\n", FileName_Tree );
      exit( 1 );
   }

   if ( OutputXYZ == WRONG )
   {
      cerr << "ERROR : please provide the targeted output data (-n output mode) !!" << endl;
      exit( 1 );
   }

   if ( TargetLevel >= NLEVEL  ||  TargetLevel < 0 )
   {
      cerr << "ERROR : incorrect TargetLevel input (-l TargetLevel) !!" << endl;
      exit( 1 );
   }

   if ( OutputXYZ < 1  ||  OutputXYZ > 7 )
   {
      cerr << "ERROR : incorrect OutputXYZ input (-n output mode) !!" << endl;
      exit( 1 );
   }

   if ( OutputFormat < 1  ||  OutputFormat > 3 )
   {
      cerr << "ERROR : incorrect output format (-f 1/2/3) !!" << endl;
      exit( 1 );
   }

#  ifndef SUPPORT_HDF5
   if ( OutputFormat == 2 )
   {
      cerr << "ERROR : must enable SUPPORT_HDF5 in the Makefile for the HDF5 output format (i.e., -f 2) !!" << endl;
      exit( 1 );
   }
#  endif

   if (  ( OutputFormat == 2 || OutputFormat == 3 )  &&  ( NGPU_X[0] != 1 || NGPU_X[1] != 1 )  )
   {
      cerr << "ERROR : only slab decomposition is allowed (-p 1 -q 1 -r NP) for the HDF5 and C-binary output formats !!" << endl;
      exit( 1 );
   }

   if ( ExtBC < 0  ||  ExtBC > 1 )
   {
      cerr << "ERROR : incorrect ExtBC input (-B 0/1) !!" << endl;
      exit( 1 );
   }

   if ( ExtBC != 0  &&  Shift2Center )
   {
      cerr << "ERROR : option -c only works with the periodic BC. (-B 0) !!" << endl;
      exit( 1 );
   }

   if ( Int_MonoCoeff < 1.0  ||  Int_MonoCoeff > 4.0 )
      Aux_Error( ERROR_INFO, "Int_MonoCoeff (%14.7e) is not within the correct range (1<=Int_MonoCoeff<=4) !!\n",
                 Int_MonoCoeff );

   if ( IntScheme != INT_MINMOD3D  &&  IntScheme != INT_MINMOD1D  &&
        IntScheme != INT_VANLEER   &&  IntScheme != INT_CQUAD     &&
        IntScheme != INT_QUAD      &&  IntScheme != INT_CQUAR     &&
        IntScheme != INT_QUAR )
      Aux_Error( ERROR_INFO, "unsupported interpolation scheme (%d) !!\n", IntScheme );

   if ( BufSize <= 0  ||  BufSize > PS1  ||  BufSize%2 != 0 )
      Aux_Error( ERROR_INFO, "unsupported buffer size (%d) --> must be a positive even number <= %d !!\n", BufSize, PS1 );

#  if ( MODEL == ELBDM )
   if ( ELBDM_IntPhase  &&  IntScheme == INT_MINMOD1D )
      Aux_Error( ERROR_INFO, "MinMod1D does not support interpolation on phase !!\n" );
#  endif


   if ( MyRank == 0 )   Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : ReadOption



//-------------------------------------------------------------------------------------------------------
// Function    :  TakeNote
// Description :  Output the simulation parameters
//
// Parameter   :  argc, argv  : command-line arguments
//-------------------------------------------------------------------------------------------------------
void TakeNote( int argc, char *argv[] )
{

   cout << "TakeNote ... " << endl;

   cout << "==========================================================" << endl;

   Aux_Message( stdout, "Command-line arguments :\n" );
   for (int v=0; v<argc; v++)    Aux_Message( stdout, " %s", argv[v] );
   Aux_Message( stdout, "\n" );
   cout << "----------------------------------------------------------"   << endl;

   cout << "DumpID            = " << DumpID                               << endl;
   cout << "Time              = " << Time[0]                              << endl;
   cout << "Step              = " << Step                                 << endl;
   cout << "NX0_TOT[0]        = " << NX0_TOT[0]                           << endl;
   cout << "NX0_TOT[1]        = " << NX0_TOT[1]                           << endl;
   cout << "NX0_TOT[2]        = " << NX0_TOT[2]                           << endl;
   cout << "BoxScale[0]       = " << amr.BoxScale[0]                      << endl;
   cout << "BoxScale[1]       = " << amr.BoxScale[1]                      << endl;
   cout << "BoxScale[2]       = " << amr.BoxScale[2]                      << endl;
   cout << "BoxSize[0]        = " << amr.BoxSize[0]                       << endl;
   cout << "BoxSize[1]        = " << amr.BoxSize[1]                       << endl;
   cout << "BoxSize[2]        = " << amr.BoxSize[2]                       << endl;
   cout << "----------------------------------------------------------"   << endl;

#  if   ( MODEL == HYDRO )
   cout << "MODEL             = " << "HYDRO"                              << endl;
#  elif ( MODEL == MHD )
   cout << "MODEL             = " << "MHD"                                << endl;
#  elif ( MODEL == ELBDM )
   cout << "MODEL             = " << "ELBDM"                              << endl;
#  else
#  error : ERROR : unsupported MODEL !!
#  endif // MODEL
#  ifdef FLOAT8
   cout << "FLOAT8            = " << "ON"                                 << endl;
#  else
   cout << "FLOAT8            = " << "OFF"                                << endl;
#  endif
#  ifdef OPENMP
   cout << "OPENMP            = " << "ON"                                 << endl;
   cout << "Number of threads = " << OMP_NThread                          << endl;
#  else
   cout << "OPENMP            = " << "OFF"                                << endl;
#  endif
#  ifdef GAMER_DEBUG
   cout << "GAMER_DEBUG       = " << "ON"                                 << endl;
#  else
   cout << "GAMER_DEBUG       = " << "OFF"                                << endl;
#  endif
#  ifdef SUPPORT_HDF5
   cout << "SUPPORT_HDF5      = " << "ON"                                 << endl;
#  else
   cout << "SUPPORT_HDF5      = " << "OFF"                                << endl;
#  endif
#  ifdef SERIAL
   cout << "SERIAL            = " << "ON"                                 << endl;
#  else
   cout << "SERIAL            = " << "OFF"                                << endl;
#  endif
   cout << "BufSize           = " << BufSize                              << endl;
   cout << "NLoad             = " << NLoad                                << endl;
   cout << "NOut              = " << NOut                                 << endl;
   cout << "OutputFormat      = " << OutputFormat                         << endl;
   cout << "OutputPot         = " << ( (OutputPot      ) ? "YES" : "NO" ) << endl;
   cout << "OutputParDens     = " << OutputParDens                        << endl;
   cout << "OutputPres        = " << ( (OutputPres     ) ? "YES" : "NO" ) << endl;
   cout << "OutputTemp        = " << ( (OutputTemp     ) ? "YES" : "NO" ) << endl;
   cout << "OutputDivVel      = " << ( (OutputDivVel   ) ? "YES" : "NO" ) << endl;
   cout << "OutputCurlVel     = " << ( (OutputCurlVel  ) ? "YES" : "NO" ) << endl;
   cout << "OutputELBDM_Vel   = " << ( (OutputELBDM_Vel) ? "YES" : "NO" ) << endl;
   cout << "ConvertCom2PhyV   = " << ( (ConvertCom2PhyV) ? "YES" : "NO" ) << endl;
   if ( OutputTemp )
   cout << "Convert2Temp      = " << Convert2Temp                         << endl;
   cout << "Output option     = " << OutputXYZ                            << endl;
   cout << "NGPU_X[0]         = " << NGPU_X[0]                            << endl;
   cout << "NGPU_X[1]         = " << NGPU_X[1]                            << endl;
   cout << "NGPU_X[2]         = " << NGPU_X[2]                            << endl;

   cout << "----------------------------------------------------------"   << endl;
   cout << "TargetLv          = " << TargetLevel                          << endl;
   cout << "dh[TargetLv]      = " << amr.dh[TargetLevel]                  << endl;
   cout << "dh[FinestLv]      = " << amr.dh[NLEVEL-1]                     << endl;
   cout << "scale[TargetLv]   = " << amr.scale[TargetLevel]               << endl;
   cout << "InputScale        = " << ( (InputScale) ? "YES" : "NO" )      << endl;
   cout << "Shift2Center      = " << ( (Shift2Center) ? "YES" : "NO" )    << endl;
   if ( Shift2Center ) {
   cout << "ShiftScale[0]     = " << ShiftScale[0]                        << endl;
   cout << "ShiftScale[1]     = " << ShiftScale[1]                        << endl;
   cout << "ShiftScale[2]     = " << ShiftScale[2]                        << endl; }
   cout << "ExtBC             = " << ExtBC                                << endl;
   cout << "Scale_Start[0]    = " << Scale_Start[0]                       << endl;
   cout << "Scale_Start[1]    = " << Scale_Start[1]                       << endl;
   cout << "Scale_Start[2]    = " << Scale_Start[2]                       << endl;
   cout << "PhyCoord_Start[0] = " << PhyCoord_Start[0]                    << endl;
   cout << "PhyCoord_Start[1] = " << PhyCoord_Start[1]                    << endl;
   cout << "PhyCoord_Start[2] = " << PhyCoord_Start[2]                    << endl;
   cout << "Idx_Start[0]      = " << Idx_Start[0]                         << endl;
   cout << "Idx_Start[1]      = " << Idx_Start[1]                         << endl;
   cout << "Idx_Start[2]      = " << Idx_Start[2]                         << endl;
   cout << "Scale_Size[0]     = " << Scale_Size[0]                        << endl;
   cout << "Scale_Size[1]     = " << Scale_Size[1]                        << endl;
   cout << "Scale_Size[2]     = " << Scale_Size[2]                        << endl;
   cout << "PhyCoord_Size[0]  = " << PhyCoord_Size[0]                     << endl;
   cout << "PhyCoord_Size[1]  = " << PhyCoord_Size[1]                     << endl;
   cout << "PhyCoord_Size[2]  = " << PhyCoord_Size[2]                     << endl;
   cout << "Idx_Size[0]       = " << Idx_Size[0]                          << endl;
   cout << "Idx_Size[1]       = " << Idx_Size[1]                          << endl;
   cout << "Idx_Size[2]       = " << Idx_Size[2]                          << endl;
   printf( "IntScheme         = %s\n", ( IntScheme == INT_MINMOD3D ) ? "MINMOD3D" :
                                       ( IntScheme == INT_MINMOD1D ) ? "MINMOD1D" :
                                       ( IntScheme == INT_VANLEER  ) ? "VANLEER"  :
                                       ( IntScheme == INT_CQUAD    ) ? "CQUAD"    :
                                       ( IntScheme == INT_QUAD     ) ? "QUAD"     :
                                       ( IntScheme == INT_CQUAR    ) ? "CQUAR"    :
                                       ( IntScheme == INT_QUAR     ) ? "QUAR"     :
                                       "UNKNOWN" );
   printf( "Int_MonoCoeff     = %13.7e\n", Int_MonoCoeff              );
   printf( "UseTree           = %s\n",     (UseTree)?"YES":"NO"       );
   printf( "FileName_Tree     = %s\n",     FileName_Tree              );
   printf( "Comoving          = %d\n",     Comoving                   );
   if ( Comoving )
   printf( "H0                = %13.7e\n", H0                         );
#  if   ( MODEL == HYDRO )
   printf( "GAMMA             = %13.7e\n", GAMMA                      );
   printf( "MU                = %13.7e\n", MU                         );
#  elif ( MODEL == ELBDM )
   printf( "ELBDM_ETA         = %13.7e\n", ELBDM_ETA );
   printf( "ELBDM_MASS        = %13.7e\n", ELBDM_MASS );
   printf( "ELBDM_IntPhase    = %s\n",    (ELBDM_IntPhase)?"YES":"NO" );
#  endif
   printf( "WithUnit          = %d\n",             WithUnit           );
   if ( WithUnit ) {
   printf( "Unit_L            = %13.7e cm\n",         Unit_L          );
   printf( "Unit_M            = %13.7e g\n",          Unit_M          );
   printf( "Unit_T            = %13.7e s\n",          Unit_T          );
   printf( "Unit_V            = %13.7e cm/s\n",       Unit_V          );
   printf( "Unit_D            = %13.7e g/cm^3\n",     Unit_D          );
   printf( "Unit_E            = %13.7e g*cm^2/s^2\n", Unit_E          );
   printf( "Unit_P            = %13.7e g/cm/s^2\n",   Unit_P          ); }

   cout << endl << "Summary :" << endl;
   cout << "------------------------------------------------------------------------------------------------------"   << endl;
   printf( " Output region: (x,y,z) = (%13.7e, %13.7e, %13.7e) --> (%13.7e, %13.7e, %13.7e)\n",
            OutCoord_Start[0], OutCoord_Start[1], OutCoord_Start[2],
            OutCoord_Start[0]+PhyCoord_Size[0], OutCoord_Start[1]+PhyCoord_Size[1], OutCoord_Start[2]+PhyCoord_Size[2] );
   printf( " Array size   : %5d * %5d * %5d\n", (OutputXYZ==4)?1:Idx_Size[0],
                                                (OutputXYZ==5)?1:Idx_Size[1],
                                                (OutputXYZ==6)?1:Idx_Size[2] );
   cout << "------------------------------------------------------------------------------------------------------"   << endl;
   cout << endl;

   cout << "=========================================================="   << endl;

} // FUNCTION : TakeNote



//-------------------------------------------------------------------------------------------------------
// Function    :  Add1Field
// Description :  Add the name and unit of a single field
//
// Parameter   :  Idx       : Current field index (call-by-reference)
//                FieldName : Field name array
//                FieldUnit : Field unit array
//                Name      : Name of the target field
//                Unit      : Unit of the target field
//
// Return      :  Next field index (Idx+1)
//-------------------------------------------------------------------------------------------------------
void Add1Field( int &Idx, char FieldName[][MAX_STRING], char FieldUnit[][MAX_STRING],
                const char Name[], const char Unit[] )
{

   strcpy( FieldName[Idx], Name );
   strcpy( FieldUnit[Idx], Unit );

   Idx ++;

} // FUNCTION : Add1Field



//-------------------------------------------------------------------------------------------------------
// Function    :  Output
// Description :  Output data
//-------------------------------------------------------------------------------------------------------
void Output()
{

   if ( MyRank == 0 )   cout << "Output ... " << endl;


   const int    scale   = amr.scale[TargetLevel];
   const double scale_2 = 0.5*scale;
   const double dh_min  = amr.dh[NLEVEL-1];
   const long   Size1v  = (long)Idx_MySize[0]*Idx_MySize[1]*Idx_MySize[2];   // total array size of one component


// 1. set the output file name(s)
   char FileName_Out[MAX_STRING], FileName_Out_Binary[NOut][MAX_STRING];
   char FieldName[NOut][MAX_STRING], FieldUnit[NOut][MAX_STRING], PassiveName[MAX_STRING];
   char *FieldUnitPtr[1]={NULL};
   int  NextIdx;

   switch ( OutputXYZ )
   {
      case 1:  sprintf( FileName_Out, "%s", "SliceX" );  break;
      case 2:  sprintf( FileName_Out, "%s", "SliceY" );  break;
      case 3:  sprintf( FileName_Out, "%s", "SliceZ" );  break;
      case 4:  sprintf( FileName_Out, "%s", "ProjX"  );  break;
      case 5:  sprintf( FileName_Out, "%s", "ProjY"  );  break;
      case 6:  sprintf( FileName_Out, "%s", "ProjZ"  );  break;
      case 7:  sprintf( FileName_Out, "%s", "Cube"   );  break;
   }


// determine the field names and units
   NextIdx = 0;

#  if   ( MODEL == HYDRO )
   Add1Field( NextIdx, FieldName, FieldUnit, "Dens", "code_mass/code_length**3" );
   Add1Field( NextIdx, FieldName, FieldUnit, "MomX", "code_mass/(code_length**2*code_time)" );
   Add1Field( NextIdx, FieldName, FieldUnit, "MomY", "code_mass/(code_length**2*code_time)" );
   Add1Field( NextIdx, FieldName, FieldUnit, "MomZ", "code_mass/(code_length**2*code_time)" );
   Add1Field( NextIdx, FieldName, FieldUnit, "Engy", "code_mass/(code_length*code_time**2)" );

#  elif ( MODEL == MHD )
#  warning : WAIT MHD !!!

#  elif ( MODEL == ELBDM )
   Add1Field( NextIdx, FieldName, FieldUnit, "Dens", "code_mass/code_length**3" );
   Add1Field( NextIdx, FieldName, FieldUnit, "Real", "code_mass**0.5/code_length**1.5" );
   Add1Field( NextIdx, FieldName, FieldUnit, "Imag", "code_mass**0.5/code_length**1.5" );

#  else
#  error : ERROR : unsupported MODEL !!
#  endif // MODEL

// here we have temporarily assumed that all passive scalars have the unit of mass density
// --> wrong for entropy
   for (int v=0; v<NCOMP_PASSIVE; v++) {
      sprintf( PassiveName, "Passive%02d", v );
      Add1Field( NextIdx, FieldName, FieldUnit, PassiveName, "code_mass/code_length**3" );
   }

   if ( OutputPot )
      Add1Field( NextIdx, FieldName, FieldUnit, "Pote", "code_length**2/code_time**2" );

   if ( OutputParDens )
      Add1Field( NextIdx, FieldName, FieldUnit, (OutputParDens==1)?"ParDens":"TotalDens", "code_mass/code_length**3" );

#  if   ( MODEL == HYDRO )
   if ( OutputPres )
      Add1Field( NextIdx, FieldName, FieldUnit, "Pres", "code_mass/(code_length*code_time**2)" );

   if ( OutputTemp )
      Add1Field( NextIdx, FieldName, FieldUnit, "Temp", "K" );

   if ( OutputDivVel )
      Add1Field( NextIdx, FieldName, FieldUnit, "DivV", "1/code_time" );

   if ( OutputCurlVel ) {
      Add1Field( NextIdx, FieldName, FieldUnit, "CurlVx", "1/code_time" );
      Add1Field( NextIdx, FieldName, FieldUnit, "CurlVy", "1/code_time" );
      Add1Field( NextIdx, FieldName, FieldUnit, "CurlVz", "1/code_time" );
   }

#  elif ( MODEL == ELBDM )
   if ( OutputELBDM_Vel ) {
      Add1Field( NextIdx, FieldName, FieldUnit, "VelX", "code_length/code_time" );
      Add1Field( NextIdx, FieldName, FieldUnit, "VelY", "code_length/code_time" );
      Add1Field( NextIdx, FieldName, FieldUnit, "VelZ", "code_length/code_time" );
   }
#  endif


// add field names
   if ( OutputFormat == 3 )
   {
      for (int v=0; v<NOut; v++)
         sprintf( FileName_Out_Binary[v], "%s_%s", FileName_Out, FieldName[v] );
   }


// add suffix
   if ( Suffix != NULL )
   {
      if ( OutputFormat == 3 )
         for (int v=0; v<NOut; v++)   strcat( FileName_Out_Binary[v], Suffix );
      else
         strcat( FileName_Out, Suffix );
   }


// add extension
   switch ( OutputFormat )
   {
      case 1:                                strcat( FileName_Out, ".txt" );              break;
      case 2:                                strcat( FileName_Out, ".hdf5" );             break;
      case 3:  for (int v=0; v<NOut; v++)    strcat( FileName_Out_Binary[v], ".cbin" );   break;
      default: Aux_Error( ERROR_INFO, "unsupported output format !!\n" );                 break;
   }



// 2. remove existing files
   if ( MyRank == 0 )
   {
      if ( OutputFormat == 3 )
         for (int v=0; v<NOut; v++)    remove( FileName_Out_Binary[v] );
      else
         remove( FileName_Out );
   }
   MPI_Barrier( MPI_COMM_WORLD );



// 3. initialize the HDF5 output
#  ifdef SUPPORT_HDF5
   hid_t   H5_FileID, H5_GroupID_Data, H5_SpaceID_Data;
   herr_t  H5_Status;
   int     Idx_AllRankSizeZ[NGPU];

   if ( OutputFormat == 2 )
   {
      hid_t   H5_DataCreatePropList, H5_AttID_Unit, H5_SpaceID_Scalar, H5_TypeID_VarStr;
      hsize_t H5_SetDims_Data[3];

//    create the data space
      H5_SetDims_Data[0] = ( OutputXYZ == 6 ) ? 1 : Idx_Size[2];
      H5_SetDims_Data[1] = ( OutputXYZ == 5 ) ? 1 : Idx_Size[1];
      H5_SetDims_Data[2] = ( OutputXYZ == 4 ) ? 1 : Idx_Size[0];

      H5_SpaceID_Data = H5Screate_simple( 3, H5_SetDims_Data, NULL );
      if ( H5_SpaceID_Data < 0 )   Aux_Error( ERROR_INFO, "failed to create the space \"%s\" !!\n", "H5_SpaceID_Data" );

      if ( MyRank == 0 )
      {
//       create file
         H5_FileID = H5Fcreate( FileName_Out, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT );
         if ( H5_FileID < 0 )    Aux_Error( ERROR_INFO, "failed to create the HDF5 file \"%s\" !!\n", FileName_Out );

//       create and set the "Info" compound datatype
         SetHDF5Info( H5_FileID );

//       create the "Data" group
         H5_GroupID_Data = H5Gcreate( H5_FileID, "Data", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
         if ( H5_GroupID_Data < 0 )    Aux_Error( ERROR_INFO, "failed to create the group \"%s\" !!\n", "Data" );

//       create the dataset of each field
//       --> do NOT write fill values to any dataset for higher I/O performance
         H5_DataCreatePropList = H5Pcreate( H5P_DATASET_CREATE );
         H5_Status             = H5Pset_fill_time( H5_DataCreatePropList, H5D_FILL_TIME_NEVER );
         H5_SpaceID_Scalar     = H5Screate( H5S_SCALAR );
         H5_TypeID_VarStr      = H5Tcopy( H5T_C_S1 );
         H5_Status             = H5Tset_size( H5_TypeID_VarStr, H5T_VARIABLE );

         for (int v=0; v<NOut; v++)
         {
            hid_t H5_SetID_Data = H5Dcreate( H5_GroupID_Data, FieldName[v], H5T_GAMER_REAL, H5_SpaceID_Data,
                                             H5P_DEFAULT, H5_DataCreatePropList, H5P_DEFAULT );
            if ( H5_SetID_Data < 0 )  Aux_Error( ERROR_INFO, "failed to create the dataset \"%s\" !!\n", FieldName[v] );

//          add the field unit as an attribute
            H5_AttID_Unit = H5Acreate( H5_SetID_Data, "Unit", H5_TypeID_VarStr, H5_SpaceID_Scalar,
                                       H5P_DEFAULT, H5P_DEFAULT );
            if ( H5_AttID_Unit < 0 )   Aux_Error( ERROR_INFO, "failed to create the attribute \"%s\" !!\n", "Unit" );
            FieldUnitPtr[0] = FieldUnit[v];  // somehow using FieldUnit[v] below doesn't work ...
            H5_Status = H5Awrite( H5_AttID_Unit, H5_TypeID_VarStr, FieldUnitPtr );

            H5_Status = H5Aclose( H5_AttID_Unit );
            H5_Status = H5Dclose( H5_SetID_Data );
         }

//       close the file and group
         H5_Status = H5Pclose( H5_DataCreatePropList );
         H5_Status = H5Sclose( H5_SpaceID_Scalar );
         H5_Status = H5Tclose( H5_TypeID_VarStr );
         H5_Status = H5Gclose( H5_GroupID_Data );
         H5_Status = H5Fclose( H5_FileID );
      } // if ( MyRank == 0 )

//    get the array size of all ranks (assuming slab decomposition)
      MPI_Allgather( &Idx_MySize[2], 1, MPI_INT, Idx_AllRankSizeZ, 1, MPI_INT, MPI_COMM_WORLD );
   } // if ( OutputFormat == 2 )
#  endif // #ifdef SUPPORT_HDF5



// 3. output data
   int    ii, jj, kk;
   long   ID;
   double x, y, z;
   real   u[NOut];

   for (int TargetRank=0; TargetRank<NGPU; TargetRank++)
   {
//    for the output-slice operation, the useless rank will NOT output any data because one of the Idx_MySize[x]
//    will be equal to zero
      if (   MyRank == TargetRank  &&  (      OutputXYZ  < 4
                                         ||   OutputXYZ == 7
                                         || ( OutputXYZ == 4 && MyRank_X[0] == 0 )
                                         || ( OutputXYZ == 5 && MyRank_X[1] == 0 )
                                         || ( OutputXYZ == 6 && MyRank_X[2] == 0 )  )   )
      {
//       3-1. text file --> all components will be outputted to the same file
         if ( OutputFormat == 1 )
         {
            FILE *File = fopen( FileName_Out, "a" );

            if ( TargetRank == 0 )
            {
               NextIdx = 0;

//             general fields
               fprintf( File, "#%9s %10s %10s %20s %20s %20s", "i", "j", "k", "x", "y", "z" );
               for (int v=0; v<NCOMP_FLUID; v++)      fprintf( File, " %13s", FieldName[NextIdx++] );
               for (int v=0; v<NCOMP_PASSIVE; v++)    fprintf( File, " %13s", FieldName[NextIdx++] );
               if ( OutputPot )                       fprintf( File, " %13s", FieldName[NextIdx++] );
               if ( OutputParDens )                   fprintf( File, " %13s", FieldName[NextIdx++] );

//             model-specific fields
#              if   ( MODEL == HYDRO )
               if ( OutputPres )                      fprintf( File, " %13s", FieldName[NextIdx++] );
               if ( OutputTemp )                      fprintf( File, " %13s", FieldName[NextIdx++] );
               if ( OutputDivVel )                    fprintf( File, " %13s", FieldName[NextIdx++] );
               if ( OutputCurlVel ) {                 fprintf( File, " %13s", FieldName[NextIdx++] );
                                                      fprintf( File, " %13s", FieldName[NextIdx++] );
                                                      fprintf( File, " %13s", FieldName[NextIdx++] ); }
#              elif ( MODEL == MHD )
#              warning : WAIT MHD !!!

#              elif ( MODEL == ELBDM )
               if ( OutputELBDM_Vel ) {               fprintf( File, " %13s", FieldName[NextIdx++] );
                                                      fprintf( File, " %13s", FieldName[NextIdx++] );
                                                      fprintf( File, " %13s", FieldName[NextIdx++] ); }

#              else
#              error : ERROR : unsupported MODEL !!
#              endif // MODEL

               fprintf( File, "\n" );
            } // if ( TargetRank == 0 )


            for (int k=0; k<Idx_MySize[2]; k++)  {  kk = ( k + Idx_MyStart[2] )*scale;    z = (kk+scale_2)*dh_min;
            for (int j=0; j<Idx_MySize[1]; j++)  {  jj = ( j + Idx_MyStart[1] )*scale;    y = (jj+scale_2)*dh_min;
            for (int i=0; i<Idx_MySize[0]; i++)  {  ii = ( i + Idx_MyStart[0] )*scale;    x = (ii+scale_2)*dh_min;

               NextIdx = NCOMP_TOTAL;
               ID      = ( (long)k*Idx_MySize[1] + j )*Idx_MySize[0] + i;

               for (int v=0; v<NOut; v++)    u[v] = OutputArray[ ID + (long)v*Size1v ];

               fprintf( File, "%10d %10d %10d %20.14e %20.14e %20.14e", ii, jj, kk, x, y, z );

               for (int v=0; v<NCOMP_TOTAL; v++)   fprintf( File, " %13.6e", u[v] );

               if ( OutputPot )                    fprintf( File, " %13.6e", u[NextIdx++] );
               if ( OutputParDens )                fprintf( File, " %13.6e", u[NextIdx++] );

#              if   ( MODEL == HYDRO )
               if ( OutputPres )                   fprintf( File, " %13.6e", u[NextIdx++] );
               if ( OutputTemp )                   fprintf( File, " %13.6e", u[NextIdx++] );
               if ( OutputDivVel )                 fprintf( File, " %13.6e", u[NextIdx++] );
               if ( OutputCurlVel ) {              fprintf( File, " %13.6e", u[NextIdx++] );
                                                   fprintf( File, " %13.6e", u[NextIdx++] );
                                                   fprintf( File, " %13.6e", u[NextIdx++] ); }

#              elif ( MODEL == MHD )
#              warning : WAIT MHD !!!

#              elif ( MODEL == ELBDM )
               if ( OutputELBDM_Vel ) {            fprintf( File, " %13.6e", u[NextIdx++] );
                                                   fprintf( File, " %13.6e", u[NextIdx++] );
                                                   fprintf( File, " %13.6e", u[NextIdx++] ); }

#              else
#              error : ERROR : unsupported MODEL !!
#              endif // MODEL

               fprintf( File, "\n" );

            }}} // i,j,k

            fclose( File );
         } // if ( OutputFormat == 1 )


//       3-2. HDF5 file --> different components will be outputted to different datasets in the same file
#        ifdef SUPPORT_HDF5
         else if ( OutputFormat == 2 )
         {
            hid_t   H5_MemID_Data;
            hsize_t H5_MemDims_Data[3], H5_Count_Data[3], H5_Offset_Data[3];

//          HDF5 file must be synchronized before being written by the next rank
            SyncHDF5File( FileName_Out );

//          reopen the file and group
            H5_FileID = H5Fopen( FileName_Out, H5F_ACC_RDWR, H5P_DEFAULT );
            if ( H5_FileID < 0 )    Aux_Error( ERROR_INFO, "failed to open the HDF5 file \"%s\" !!\n", FileName_Out );

            H5_GroupID_Data = H5Gopen( H5_FileID, "Data", H5P_DEFAULT );
            if ( H5_GroupID_Data < 0 )   Aux_Error( ERROR_INFO, "failed to open the group \"%s\" !!\n", "Data" );

//          set the memory space
            H5_MemDims_Data[0] = Idx_MySize[2];
            H5_MemDims_Data[1] = Idx_MySize[1];
            H5_MemDims_Data[2] = Idx_MySize[0];

            H5_MemID_Data = H5Screate_simple( 3, H5_MemDims_Data, NULL );
            if ( H5_MemID_Data < 0 )  Aux_Error( ERROR_INFO, "failed to create the space \"%s\" !!\n", "H5_MemID_Data" );

//          set the subset of the dataspace
            H5_Offset_Data[0] = 0;
            H5_Offset_Data[1] = 0;
            H5_Offset_Data[2] = 0;
            for (int r=0; r<MyRank; r++)  H5_Offset_Data[0] += Idx_AllRankSizeZ[r];

            H5_Count_Data [0] = Idx_MySize[2];
            H5_Count_Data [1] = Idx_MySize[1];
            H5_Count_Data [2] = Idx_MySize[0];

            H5_Status = H5Sselect_hyperslab( H5_SpaceID_Data, H5S_SELECT_SET, H5_Offset_Data, NULL, H5_Count_Data, NULL );
            if ( H5_Status < 0 )   Aux_Error( ERROR_INFO, "failed to create a hyperslab !!\n" );

//          write data to disk (one field at a time)
            for (int v=0; v<NOut; v++)
            {
               hid_t H5_SetID_Data = H5Dopen( H5_GroupID_Data, FieldName[v], H5P_DEFAULT );
               H5_Status           = H5Dwrite( H5_SetID_Data, H5T_GAMER_REAL, H5_MemID_Data, H5_SpaceID_Data, H5P_DEFAULT,
                                               OutputArray+(long)v*Size1v );
               if ( H5_Status < 0 )   Aux_Error( ERROR_INFO, "failed to output the field \"%s\" !!\n", FieldName );
               H5_Status = H5Dclose( H5_SetID_Data );
            }

            H5_Status = H5Sclose( H5_MemID_Data );
            H5_Status = H5Gclose( H5_GroupID_Data );
            H5_Status = H5Fclose( H5_FileID );
         } // else if ( OutputFormat == 2 )
#        endif // #ifdef SUPPORT_HDF5


//       3-3. C-binary file --> different components will be outputted to different files
         else if ( OutputFormat == 3 )
         {
            for (int v=0; v<NOut; v++)
            {
               FILE *File = fopen( FileName_Out_Binary[v], "ab" );

               fwrite( OutputArray+(long)v*Size1v, sizeof(real), Size1v, File );

               fclose( File );
            }
         } // else if ( OutputFormat == 3 )


         else
            Aux_Error( ERROR_INFO, "unsupported output format !!\n" );

      } // ( MyRank == TargetRank  &&  ... )

      MPI_Barrier( MPI_COMM_WORLD );

   } // for (int TargetRank=0; TargetRank<NGPU; TargetRank++)


// close all HDF5 objects
#  ifdef SUPPORT_HDF5
   if ( OutputFormat == 2 )
   {
      H5_Status = H5Sclose( H5_SpaceID_Data );
   }
#  endif


   if ( MyRank == 0 )   cout << "Output ... done" << endl;

} // FUNCTION : Output



#ifdef SUPPORT_HDF5
//-------------------------------------------------------------------------------------------------------
// Function    :  SetHDF5Info
// Description :  Create and set the HDF5 compound datatype "Info"
//
// Note        :  1. Data structure "Info_t" is defined in "HDF5_Typedef.h"
//
// Parameter   :  H5_FileID : HDF5 file ID (call-by-reference)
//-------------------------------------------------------------------------------------------------------
void SetHDF5Info( hid_t &H5_FileID )
{

// 1. create the HDF5 compound datatype "Info"
   hid_t  H5_ComID_Info;
   herr_t H5_Status;

// create the array type
   const hsize_t H5_ArrDims_3Var  = 3;  // array size of [3]
   const hid_t   H5_ArrID_3Double = H5Tarray_create( H5T_NATIVE_DOUBLE, 1, &H5_ArrDims_3Var );
   const hid_t   H5_ArrID_3Int    = H5Tarray_create( H5T_NATIVE_INT,    1, &H5_ArrDims_3Var );

// create the "variable-length string" datatype
   /*
   hid_t  H5_ArrID_VarStr;
   herr_t H5_Status;
   H5_ArrID_VarStr = H5Tcopy( H5T_C_S1 );
   H5_Status       = H5Tset_size( H5_ArrID_VarStr, H5T_VARIABLE );
   */

// get the compound type
   H5_ComID_Info = H5Tcreate( H5T_COMPOUND, sizeof(Info_t) );
   H5Tinsert( H5_ComID_Info, "DumpID",             HOFFSET(Info_t,DumpID           ),  H5T_NATIVE_INT    );
   H5Tinsert( H5_ComID_Info, "GridDimension",      HOFFSET(Info_t,GridDimension    ),  H5_ArrID_3Int     );
   H5Tinsert( H5_ComID_Info, "WithUnit",           HOFFSET(Info_t,WithUnit         ),  H5T_NATIVE_INT    );
   H5Tinsert( H5_ComID_Info, "OutputMode",         HOFFSET(Info_t,OutputMode       ),  H5T_NATIVE_INT    );
   H5Tinsert( H5_ComID_Info, "Time",               HOFFSET(Info_t,Time             ),  H5T_NATIVE_DOUBLE );
   H5Tinsert( H5_ComID_Info, "CellWidth",          HOFFSET(Info_t,CellWidth        ),  H5T_NATIVE_DOUBLE );
   H5Tinsert( H5_ComID_Info, "SubdomainSize",      HOFFSET(Info_t,SubdomainSize    ),  H5_ArrID_3Double  );
   H5Tinsert( H5_ComID_Info, "SubdomainLeftEdge",  HOFFSET(Info_t,SubdomainLeftEdge),  H5_ArrID_3Double  );
   H5Tinsert( H5_ComID_Info, "Unit_L",             HOFFSET(Info_t,Unit_L           ),  H5T_NATIVE_DOUBLE );
   H5Tinsert( H5_ComID_Info, "Unit_M",             HOFFSET(Info_t,Unit_M           ),  H5T_NATIVE_DOUBLE );
   H5Tinsert( H5_ComID_Info, "Unit_T",             HOFFSET(Info_t,Unit_T           ),  H5T_NATIVE_DOUBLE );
   H5Tinsert( H5_ComID_Info, "Unit_V",             HOFFSET(Info_t,Unit_V           ),  H5T_NATIVE_DOUBLE );
   H5Tinsert( H5_ComID_Info, "Unit_D",             HOFFSET(Info_t,Unit_D           ),  H5T_NATIVE_DOUBLE );
   H5Tinsert( H5_ComID_Info, "Unit_E",             HOFFSET(Info_t,Unit_E           ),  H5T_NATIVE_DOUBLE );
   H5Tinsert( H5_ComID_Info, "Unit_P",             HOFFSET(Info_t,Unit_P           ),  H5T_NATIVE_DOUBLE );
#  if ( MODEL == ELBDM )
   H5Tinsert( H5_ComID_Info, "ELBDM_Mass",         HOFFSET(Info_t,ELBDM_Mass       ),  H5T_NATIVE_DOUBLE );
#  endif

// free memory
   H5_Status = H5Tclose( H5_ArrID_3Int    );
   H5_Status = H5Tclose( H5_ArrID_3Double );
// H5_Status = H5Tclose( H5_ArrID_VarStr  );


// 2. fill in the structure "Info"
   Info_t Info;

   Info.DumpID     = DumpID;
   Info.WithUnit   = WithUnit;
   Info.Time       = Time[0];
   Info.OutputMode = OutputXYZ;
   Info.CellWidth  = amr.dh[TargetLevel];
   Info.Unit_L     = (WithUnit) ? Unit_L : 1.0;
   Info.Unit_M     = (WithUnit) ? Unit_M : 1.0;
   Info.Unit_T     = (WithUnit) ? Unit_T : 1.0;
   Info.Unit_V     = (WithUnit) ? Unit_V : 1.0;
   Info.Unit_D     = (WithUnit) ? Unit_D : 1.0;
   Info.Unit_E     = (WithUnit) ? Unit_E : 1.0;
   Info.Unit_P     = (WithUnit) ? Unit_P : 1.0;
#  if ( MODEL == ELBDM )
   const double Const_eV = 1.6021766208e-12;
   const double Const_c  = 2.99792458e10;
   Info.ELBDM_Mass = (WithUnit) ? ELBDM_MASS*Unit_M/( Const_eV/(Const_c*Const_c) ) : ELBDM_MASS;
#  endif

   for (int d=0; d<3; d++)
   {
//    for projection, we record the coordinates of the original 3D region to be projected
      Info.GridDimension    [d] = ( OutputXYZ == 4+d ) ? 1 : Idx_Size[d];
      Info.SubdomainSize    [d] = PhyCoord_Size [d];
      Info.SubdomainLeftEdge[d] = OutCoord_Start[d];
   }


// 3. write to disk
   hid_t H5_SpaceID_Scalar, H5_SetID_Info;

   H5_SpaceID_Scalar = H5Screate( H5S_SCALAR );
   H5_SetID_Info     = H5Dcreate( H5_FileID, "Info", H5_ComID_Info, H5_SpaceID_Scalar, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
   if ( H5_SetID_Info < 0 )   Aux_Error( ERROR_INFO, "failed to create the dataset \"%s\" !!\n", "Info" );

   H5_Status = H5Dwrite( H5_SetID_Info, H5_ComID_Info, H5S_ALL, H5S_ALL, H5P_DEFAULT, &Info );

   H5_Status = H5Dclose( H5_SetID_Info );
   H5_Status = H5Sclose( H5_SpaceID_Scalar );
   H5_Status = H5Tclose( H5_ComID_Info );

} // FUNCTION : SetHDF5Info
#endif // #ifdef SUPPORT_HDF5



//-------------------------------------------------------------------------------------------------------
// Function    :  CheckParameter
// Description :  Verify the input parameters
//-------------------------------------------------------------------------------------------------------
void CheckParameter()
{

   if ( MyRank == 0 )   cout << "   CheckParameter ... " << endl;


#  if ( MODEL != HYDRO  &&  MODEL != ELBDM )
#     error : ERROR : unsupported MODEL !!
#  endif

#  if ( defined OPENMP  &&  !defined _OPENMP )
#     error : ERROR : something is wrong in OpenMP, the macro "_OPENMP" is NOT defined !!
#  endif

   int NRank;
#  ifdef SERIAL
   NRank = 1;
#  else
   MPI_Comm_size( MPI_COMM_WORLD, &NRank );
   if ( NRank != NGPU )
   {
      fprintf( stderr, "ERROR : number of ranks (%d) != number of GPUs (%d) !!\n", NRank, NGPU );
      MPI_Exit();
   }
#  endif

   for (int d=0; d<3; d++)
   {
      if ( PhyCoord_Size[d] < amr.dh[NLEVEL-1] )
      {
         fprintf( stderr, "ERROR : incorrect PhyCoord_Size[%d] (%14.7e) --> must be >= finest cell (%14.7e)!!\n",
                  d, PhyCoord_Size[d], amr.dh[NLEVEL-1] );
         MPI_Exit();
      }

      if ( Scale_Size[d] <= 0 )
      {
         fprintf( stderr, "ERROR : incorrect Scale_Size[%d] (%d) --> must be positive !!\n", d, Scale_Size[d] );
         MPI_Exit();
      }

      if ( Idx_Size[d] <= 0 )
      {
         fprintf( stderr, "ERROR : incorrect Idx_Size[%d] (%d) --> must be positive !!\n", d, Idx_Size[d] );
         MPI_Exit();
      }

      if ( PhyCoord_Start[d] < 0.0  ||  PhyCoord_Start[d] >= amr.BoxSize[d] )
      {
         fprintf( stderr, "ERROR : incorrect PhyCoord_Start[%d] (%14.7e) --> out of range !!\n",
                  d, PhyCoord_Start[d] );
         MPI_Exit();
      }

      if ( Scale_Start[d] < 0  ||  Scale_Start[d] >= amr.BoxScale[d] )
      {
         fprintf( stderr, "ERROR : incorrect Scale_Start[%d] (%d) --> out of range !!\n", d, Scale_Start[d] );
         MPI_Exit();
      }

      if ( Idx_Start[d] < 0  ||  Idx_Start[d] >= NX0_TOT[d]*(1<<TargetLevel) )
      {
         fprintf( stderr, "ERROR : incorrect Idx_Start[%d] (%d) --> out of range !!\n", d, Idx_Start[d] );
         MPI_Exit();
      }

      if ( PhyCoord_Start[d]+PhyCoord_Size[d] > amr.BoxSize[d] )
      {
         fprintf( stderr, "ERROR : incorrect PhyCoord_Start[%d]+PhyCoord_Size[%d] (%14.7e) --> out of range !!\n",
                  d, d, PhyCoord_Start[d]+PhyCoord_Size[d] );
         MPI_Exit();
      }

      if ( Scale_Start[d]+Scale_Size[d] > amr.BoxScale[d] )
      {
         fprintf( stderr, "ERROR : incorrect Scale_Start[%d]+Scale_Size[%d] (%d) --> out of range !!\n",
                  d, d, Scale_Start[d]+Scale_Size[d] );
         MPI_Exit();
      }

      if ( Idx_Start[d]+Idx_Size[d] > NX0_TOT[d]*(1<<TargetLevel) )
      {
         fprintf( stderr, "ERROR : incorrect Idx_Start[%d]+Idx_Size[%d] (%d) --> out of range !!\n",
                  d, d, Idx_Start[d]+Idx_Size[d] );
         MPI_Exit();
      }
   } // for (int d=0; d<3; d++)

   if ( NX0_TOT[0]%(PS2*NGPU_X[0]) != 0  ||  NX0_TOT[1]%(PS2*NGPU_X[1]) != 0  ||  NX0_TOT[2]%(PS2*NGPU_X[2]) != 0 )
   {
      fprintf( stderr, "ERROR : number of base-level patches in each direction in one rank must be \"%s\" !!\n",
                 "a multiple of TWO" );
      MPI_Exit();
   }

   if ( !Comoving  &&  ConvertCom2PhyV )
   {
      ConvertCom2PhyV = false;

      if ( MyRank == 0 )
      fprintf( stderr, "WARNING : option -k (ConvertCom2PhyV) is useless for non-cosmological simulations !!\n" );
   }

   if ( OutputPres )
   {
      if ( MyRank == 0 )
      fprintf( stderr, "WARNING : option -P (OutputPres) assumes constant-gamma EoS !!\n" );
   }

#  if ( MODEL == HYDRO  ||  MODEL == MHD )
   if ( OutputTemp  &&  Convert2Temp < 0.0 )
      Aux_Error( ERROR_INFO, "Convert2Temp (%14.7e) < 0.0 for calculating temperature (use -u) !!\n",
                 Convert2Temp );
#  endif


   if ( MyRank == 0 )   cout << "   CheckParameter ... done" << endl;

} // FUNCTION : CheckParameter



//-------------------------------------------------------------------------------------------------------
// Function    :  WithinCandidateBox
// Description :  Check whether or not the input patch is within the candidate box
//
// Parameter   :  Corner   : Pointer to the three corner coordinates
//                Size     : Size of the targeted patch
//                Buf      : Buffer size attached to each side of the candidate box
//-------------------------------------------------------------------------------------------------------
bool WithinCandidateBox( const int *Corner, const int Size, const int Buf )
{

   bool Inside[3];
   int cr1[3], cr2[3];

   for (int d=0; d<3; d++)
   {
      Inside[d] = false;
      cr1   [d] = Corner[d];
      cr2   [d] = Corner[d] + Size;
   }


   for (int d=0; d<3; d++)
   {
      for (int s=0; s<3; s++)
      {
         if (  cr1[d] < CanMax[d][s]+Buf  &&  cr2[d] > CanMin[d][s]-Buf  )
         {
            Inside[d] = true;

            break;
         }
      }
   }

   return ( Inside[0] && Inside[1] && Inside[2] );

} // FUNCTION : WithinCandidateBox



//-------------------------------------------------------------------------------------------------------
// Function    :  GetCandidateBox
// Description :  Evaluate the range of the candidate box (considering the periodic boundary condition)
//-------------------------------------------------------------------------------------------------------
void GetCandidateBox()
{

   if ( MyRank == 0 )   cout << "   GetCandidateBox ... " << flush;


   const int scale = amr.scale[TargetLevel];

   for (int d=0; d<3; d++)
   for (int s=0; s<3; s++)
   {
      CanMin[d][s] = (Idx_Start[d]            )*scale + amr.BoxScale[d]*(s-1);
      CanMax[d][s] = (Idx_Start[d]+Idx_Size[d])*scale + amr.BoxScale[d]*(s-1);
   }


// set the buffer size to the size of two base-level patches
   CanBuf = 2*PATCH_SIZE*amr.scale[0];


   if ( MyRank == 0 )   cout << "done" << endl;

} // FUNCTION : GetCandidateBox



//-------------------------------------------------------------------------------------------------------
// Function    :  Init_TargetDomain
// Description :  Initialize the parameters related to parallelization and the targeted domain
//-------------------------------------------------------------------------------------------------------
void Init_TargetDomain()
{

   if ( MyRank == 0 )   cout << "   Init_TargetDomain ... " << endl;


// set parameters for parallelization
   MyRank_X[0] = MyRank % NGPU_X[0];
   MyRank_X[1] = (MyRank/NGPU_X[0]) % NGPU_X[1];
   MyRank_X[2] = (MyRank/NGPU_X[0]) / NGPU_X[1];

// sibling mpi ranks
   if ( ExtBC == 0 )    // periodic BC.
   {
      SibRank[ 0] = ( MyRank_X[2]             )          *NGPU_X[0]*NGPU_X[1] +
                    ( MyRank_X[1]             )          *NGPU_X[0]           +
                    ( MyRank_X[0]+NGPU_X[0]-1 )%NGPU_X[0];

      SibRank[ 1] = ( MyRank_X[2]             )          *NGPU_X[0]*NGPU_X[1] +
                    ( MyRank_X[1]             )          *NGPU_X[0]           +
                    ( MyRank_X[0]          +1 )%NGPU_X[0];

      SibRank[ 2] = ( MyRank_X[2]             )          *NGPU_X[0]*NGPU_X[1] +
                    ( MyRank_X[1]+NGPU_X[1]-1 )%NGPU_X[1]*NGPU_X[0]           +
                    ( MyRank_X[0]             );

      SibRank[ 3] = ( MyRank_X[2]             )          *NGPU_X[0]*NGPU_X[1] +
                    ( MyRank_X[1]          +1 )%NGPU_X[1]*NGPU_X[0]           +
                    ( MyRank_X[0]             );

      SibRank[ 4] = ( MyRank_X[2]+NGPU_X[2]-1 )%NGPU_X[2]*NGPU_X[0]*NGPU_X[1] +
                    ( MyRank_X[1]             )          *NGPU_X[0]           +
                    ( MyRank_X[0]             );

      SibRank[ 5] = ( MyRank_X[2]          +1 )%NGPU_X[2]*NGPU_X[0]*NGPU_X[1] +
                    ( MyRank_X[1]             )          *NGPU_X[0]           +
                    ( MyRank_X[0]             );

      SibRank[ 6] = ( MyRank_X[2]             )          *NGPU_X[0]*NGPU_X[1] +
                    ( MyRank_X[1]+NGPU_X[1]-1 )%NGPU_X[1]*NGPU_X[0]           +
                    ( MyRank_X[0]+NGPU_X[0]-1 )%NGPU_X[0];

      SibRank[ 7] = ( MyRank_X[2]             )          *NGPU_X[0]*NGPU_X[1] +
                    ( MyRank_X[1]+NGPU_X[1]-1 )%NGPU_X[1]*NGPU_X[0]           +
                    ( MyRank_X[0]          +1 )%NGPU_X[0];

      SibRank[ 8] = ( MyRank_X[2]             )          *NGPU_X[0]*NGPU_X[1] +
                    ( MyRank_X[1]          +1 )%NGPU_X[1]*NGPU_X[0]           +
                    ( MyRank_X[0]+NGPU_X[0]-1 )%NGPU_X[0];

      SibRank[ 9] = ( MyRank_X[2]             )          *NGPU_X[0]*NGPU_X[1] +
                    ( MyRank_X[1]          +1 )%NGPU_X[1]*NGPU_X[0]           +
                    ( MyRank_X[0]          +1 )%NGPU_X[0];

      SibRank[10] = ( MyRank_X[2]+NGPU_X[2]-1 )%NGPU_X[2]*NGPU_X[0]*NGPU_X[1] +
                    ( MyRank_X[1]+NGPU_X[1]-1 )%NGPU_X[1]*NGPU_X[0]           +
                    ( MyRank_X[0]             );

      SibRank[11] = ( MyRank_X[2]+NGPU_X[2]-1 )%NGPU_X[2]*NGPU_X[0]*NGPU_X[1] +
                    ( MyRank_X[1]          +1 )%NGPU_X[1]*NGPU_X[0]           +
                    ( MyRank_X[0]             );

      SibRank[12] = ( MyRank_X[2]          +1 )%NGPU_X[2]*NGPU_X[0]*NGPU_X[1] +
                    ( MyRank_X[1]+NGPU_X[1]-1 )%NGPU_X[1]*NGPU_X[0]           +
                    ( MyRank_X[0]             );

      SibRank[13] = ( MyRank_X[2]          +1 )%NGPU_X[2]*NGPU_X[0]*NGPU_X[1] +
                    ( MyRank_X[1]          +1 )%NGPU_X[1]*NGPU_X[0]           +
                    ( MyRank_X[0]             );

      SibRank[14] = ( MyRank_X[2]+NGPU_X[2]-1 )%NGPU_X[2]*NGPU_X[0]*NGPU_X[1] +
                    ( MyRank_X[1]             )          *NGPU_X[0]           +
                    ( MyRank_X[0]+NGPU_X[0]-1 )%NGPU_X[0];

      SibRank[15] = ( MyRank_X[2]          +1 )%NGPU_X[2]*NGPU_X[0]*NGPU_X[1] +
                    ( MyRank_X[1]             )*NGPU_X[0]                     +
                    ( MyRank_X[0]+NGPU_X[0]-1 )%NGPU_X[0];

      SibRank[16] = ( MyRank_X[2]+NGPU_X[2]-1 )%NGPU_X[2]*NGPU_X[0]*NGPU_X[1] +
                    ( MyRank_X[1]             )          *NGPU_X[0]           +
                    ( MyRank_X[0]          +1 )%NGPU_X[0];

      SibRank[17] = ( MyRank_X[2]          +1 )%NGPU_X[2]*NGPU_X[0]*NGPU_X[1] +
                    ( MyRank_X[1]             )          *NGPU_X[0]           +
                    ( MyRank_X[0]          +1 )%NGPU_X[0];

      SibRank[18] = ( MyRank_X[2]+NGPU_X[2]-1 )%NGPU_X[2]*NGPU_X[0]*NGPU_X[1] +
                    ( MyRank_X[1]+NGPU_X[1]-1 )%NGPU_X[1]*NGPU_X[0]           +
                    ( MyRank_X[0]+NGPU_X[0]-1 )%NGPU_X[0];

      SibRank[19] = ( MyRank_X[2]+NGPU_X[2]-1 )%NGPU_X[2]*NGPU_X[0]*NGPU_X[1] +
                    ( MyRank_X[1]+NGPU_X[1]-1 )%NGPU_X[1]*NGPU_X[0]           +
                    ( MyRank_X[0]          +1 )%NGPU_X[0];

      SibRank[20] = ( MyRank_X[2]+NGPU_X[2]-1 )%NGPU_X[2]*NGPU_X[0]*NGPU_X[1] +
                    ( MyRank_X[1]          +1 )%NGPU_X[1]*NGPU_X[0]           +
                    ( MyRank_X[0]+NGPU_X[0]-1 )%NGPU_X[0];

      SibRank[21] = ( MyRank_X[2]+NGPU_X[2]-1 )%NGPU_X[2]*NGPU_X[0]*NGPU_X[1] +
                    ( MyRank_X[1]          +1 )%NGPU_X[1]*NGPU_X[0]           +
                    ( MyRank_X[0]          +1 )%NGPU_X[0];

      SibRank[22] = ( MyRank_X[2]          +1 )%NGPU_X[2]*NGPU_X[0]*NGPU_X[1] +
                    ( MyRank_X[1]+NGPU_X[1]-1 )%NGPU_X[1]*NGPU_X[0]           +
                    ( MyRank_X[0]+NGPU_X[0]-1 )%NGPU_X[0];

      SibRank[23] = ( MyRank_X[2]          +1 )%NGPU_X[2]*NGPU_X[0]*NGPU_X[1] +
                    ( MyRank_X[1]+NGPU_X[1]-1 )%NGPU_X[1]*NGPU_X[0]           +
                    ( MyRank_X[0]          +1 )%NGPU_X[0];

      SibRank[24] = ( MyRank_X[2]          +1 )%NGPU_X[2]*NGPU_X[0]*NGPU_X[1] +
                    ( MyRank_X[1]          +1 )%NGPU_X[1]*NGPU_X[0]           +
                    ( MyRank_X[0]+NGPU_X[0]-1 )%NGPU_X[0];

      SibRank[25] = ( MyRank_X[2]          +1 )%NGPU_X[2]*NGPU_X[0]*NGPU_X[1] +
                    ( MyRank_X[1]          +1 )%NGPU_X[1]*NGPU_X[0]           +
                    ( MyRank_X[0]          +1 )%NGPU_X[0];
   } // if ( ExtBC == 0 )

   else  // non-periodic BC.
   {
      const int Buf     = 1;
      const int Size[3] = { NGPU_X[0]+2*Buf, NGPU_X[1]+2*Buf, NGPU_X[2]+2*Buf };

      int ii, jj, kk, Width[3], Disp[3], ID1, ID2;
      int *RankMap = new int [ Size[0]*Size[1]*Size[2] ];

//    interior ranks
      for (int k=Buf; k<Buf+NGPU_X[2]; k++)  {  kk = k - Buf;
      for (int j=Buf; j<Buf+NGPU_X[1]; j++)  {  jj = j - Buf;
      for (int i=Buf; i<Buf+NGPU_X[0]; i++)  {  ii = i - Buf;

         ID1          = ( k *       Size[1] + j  )*       Size[0] + i;
         ID2          = ( kk*NGPU_X[1] + jj )*NGPU_X[0] + ii;
         RankMap[ID1] = ID2;

      }}}

//    exterior ranks (set to "SIB_OFFSET_NONPERIODI-sib" --> similar to sibling patches lying outside the simulation domain)
      for (int s=0; s<26; s++)
      {
         for (int d=0; d<3; d++)
         {
            Width[d] = TABLE_01( s, 'x'+d, Buf, NGPU_X[d], Buf                );
            Disp [d] = TABLE_01( s, 'x'+d,   0,            Buf, Buf+NGPU_X[d] );
         }

         for (int k=Disp[2]; k<Disp[2]+Width[2]; k++)
         for (int j=Disp[1]; j<Disp[1]+Width[1]; j++)
         for (int i=Disp[0]; i<Disp[0]+Width[0]; i++)
         {
            ID1          = ( k*Size[1] + j )*Size[0] + i;
            RankMap[ID1] = SIB_OFFSET_NONPERIODIC - s;   // set to "SIB_OFFSET_NONPERIODIC-s" for the non-periodic B.C.
         }
      }

//    record the sibling ranks
      for (int s=0; s<26; s++)
      {
         for (int d=0; d<3; d++)    Disp[d] = TABLE_01( s, 'x'+d, -1, 0, +1 );

         ID1        = ( (Buf+MyRank_X[2]+Disp[2])*Size[1] + (Buf+MyRank_X[1]+Disp[1]) )*Size[0] + (Buf+MyRank_X[0]+Disp[0]);
         SibRank[s] = RankMap[ID1];
      }

      delete [] RankMap;

   } // if ( ExtBC == 0 ) ... else ...


// set the default (start/size) of (cell scales/physical coordinates)
   if ( InputScale )
   {
      for (int d=0; d<3; d++)
      {
         if ( Scale_Start[d] == WRONG )   Scale_Start[d] = 0;
         if ( Scale_Size [d] == WRONG )   Scale_Size [d] = amr.BoxScale[d] - Scale_Start[d];

         PhyCoord_Start[d] = amr.BoxSize[d]*( (double)Scale_Start[d]/(double)amr.BoxScale[d] );
         PhyCoord_Size [d] = amr.BoxSize[d]*( (double)Scale_Size [d]/(double)amr.BoxScale[d] );
      }
   }

   else
   {
      for (int d=0; d<3; d++)
      {
         if ( PhyCoord_Start[d] == WRONG )   PhyCoord_Start[d] = 0.0;
         if ( PhyCoord_Size [d] == WRONG )   PhyCoord_Size [d] = amr.BoxSize[d] - PhyCoord_Start[d];

         Scale_Start[d] = int(amr.BoxScale[d]*PhyCoord_Start[d]/amr.BoxSize[d]);
         Scale_Size [d] = int(amr.BoxScale[d]*PhyCoord_Size [d]/amr.BoxSize[d]);
      }
   }

   switch ( OutputXYZ )
   {
      case 1 : Scale_Size[0] = 1;   PhyCoord_Size[0] = amr.dh[NLEVEL-1];   break;
      case 2 : Scale_Size[1] = 1;   PhyCoord_Size[1] = amr.dh[NLEVEL-1];   break;
      case 3 : Scale_Size[2] = 1;   PhyCoord_Size[2] = amr.dh[NLEVEL-1];   break;
   }


   if ( ExtBC == 0 )
   {
//    reset the starting coords to be inside the box
      for (int d=0; d<3; d++)
      {
         const int    Scale_Start0    = Scale_Start   [d];
         const double PhyCoord_Start0 = PhyCoord_Start[d];

         if      ( Scale_Start[d] < 0 )                  Scale_Start[d] += amr.BoxScale[d];
         else if ( Scale_Start[d] >= amr.BoxScale[d] )   Scale_Start[d] -= amr.BoxScale[d];

         if ( Scale_Start[d] != Scale_Start0 )
         {
            PhyCoord_Start[d] = amr.BoxSize[d]*( (double)Scale_Start[d]/(double)amr.BoxScale[d] );

            if ( MyRank == 0 )
            {
               Aux_Message( stderr, "WARNING : starting Scale/PhyCoord[%d] (%d/%13.7e) lies outside the simulation box\n",
                            d, Scale_Start0, PhyCoord_Start0 );
               Aux_Message( stderr, "          --> reset starting Scale/PhyCoord to (%d/%13.7e)\n",
                            Scale_Start[d], PhyCoord_Start[d] );
            }
         }
      }


//    turn on "Shift2Center" if the target region lying outside the simulation box and periodic BC. is adopted
      for (int d=0; d<3; d++)
      {
         if ( Scale_Start[d]+Scale_Size[d] > amr.BoxScale[d] )
         {
            Shift2Center = true;

            if ( MyRank == 0 )
            {
               Aux_Message( stderr, "WARNING : target region along [%d] lies outside the simulation box !!\n", d );
               Aux_Message( stderr, "          Scale   : %13d + %13d = %13d > %13d\n",
                            Scale_Start[d], Scale_Size[d], Scale_Start[d]+Scale_Size[d], amr.BoxScale[d] );
               Aux_Message( stderr, "          PhyCoord: %13.7e + %13.7e = %13.7e > %13.7e\n",
                            PhyCoord_Start[d], PhyCoord_Size[d], PhyCoord_Start[d]+PhyCoord_Size[d], amr.BoxSize[d] );
               Aux_Message( stderr, "          --> Shift2Center is turned ON automatically\n" );
            }

            break;
         }
      }
   }


// shift the target region to be close to the box center
   if ( Shift2Center )
   {
      const int ShiftUnit = PS2*amr.scale[0]; // ShiftScale must be a multiple of the base-level patch group

      for (int d=0; d<3; d++)
      {
         ShiftScale[d] = ShiftUnit*(  (int)round( (0.5*amr.BoxScale[d]-Scale_Start[d]-0.5*Scale_Size[d])/ShiftUnit )  );

         Scale_Start   [d] = ( Scale_Start[d] + ShiftScale[d] + amr.BoxScale[d] ) % amr.BoxScale[d];
         PhyCoord_Start[d] = amr.BoxSize[d]*( (double)Scale_Start[d]/(double)amr.BoxScale[d] );

//       check
         if ( Scale_Size[d] > amr.BoxScale[d]-ShiftUnit  &&  MyRank == 0 )
         {
            Aux_Message( stderr, "WARNING : Scale_Size[%d] (%d) > Box.Scale[%d]-ShiftUnit (%d)\n",
                         d, Scale_Size[d], d, amr.BoxScale[d]-ShiftUnit );
            Aux_Message( stderr, "          --> It may fail unless the target region properly aligns with grids\n" );
         }

         if ( PhyCoord_Start[d] < 0.0  ||  PhyCoord_Start[d] >= amr.BoxSize[d] )
         {
            fprintf( stderr, "ERROR : incorrect PhyCoord_Start[%d] (%14.7e) --> out of range !!\n",
                     d, PhyCoord_Start[d] );
            MPI_Exit();
         }

         if ( Scale_Start[d] < 0  ||  Scale_Start[d] >= amr.BoxScale[d] )
         {
            fprintf( stderr, "ERROR : incorrect Scale_Start[%d] (%d) --> out of range !!\n", d, Scale_Start[d] );
            MPI_Exit();
         }

         if ( PhyCoord_Start[d]+PhyCoord_Size[d] > amr.BoxSize[d] )
         {
            fprintf( stderr, "ERROR : incorrect PhyCoord_Start[%d]+PhyCoord_Size[%d] (%14.7e) --> out of range !!\n",
                     d, d, PhyCoord_Start[d]+PhyCoord_Size[d] );
            MPI_Exit();
         }

         if ( Scale_Start[d]+Scale_Size[d] > amr.BoxScale[d] )
         {
            fprintf( stderr, "ERROR : incorrect Scale_Start[%d]+Scale_Size[%d] (%d) --> out of range !!\n",
                     d, d, Scale_Start[d]+Scale_Size[d] );
            MPI_Exit();
         }
      } // for (int d=0; d<3; d++)
   } // if ( Shift2Center )


// set the targeted array indices and ranges
   for (int d=0; d<3; d++)
   {
      Idx_Start[d] = Scale_Start[d]   /amr.scale[TargetLevel];
      Idx_Size [d] = (Scale_Size[d]-1)/amr.scale[TargetLevel] + 1;
   }


// recalculate scale and physical coords according to Idx
   for (int d=0; d<3; d++)
   {
      if ( d+1 == OutputXYZ )    continue;

      Scale_Start   [d] = Idx_Start[d]*amr.scale[TargetLevel];
      Scale_Size    [d] = Idx_Size [d]*amr.scale[TargetLevel];

      PhyCoord_Start[d] = Idx_Start[d]*amr.dh   [TargetLevel];
      PhyCoord_Size [d] = Idx_Size [d]*amr.dh   [TargetLevel];
   }


// set the local array size and the starting indices
   int Min, Max, NGrid, X_End;

   for (int d=0; d<3; d++)
   {
      NGrid = NX0[d]*(1<<TargetLevel);
      Min   = (MyRank_X[d]  )*NGrid;
      Max   = (MyRank_X[d]+1)*NGrid-1;

//    check whether or not the sub-domain is within the targeted range
      if ( Min < Idx_Start[d]+Idx_Size[d]  &&  Max >= Idx_Start[d] )
      {
         Idx_MyStart[d] = ( Idx_Start[d]             > Min ) ? Idx_Start[d]-Min : 0;
         X_End          = ( Idx_Start[d]+Idx_Size[d] > Max ) ? NGrid-1 : NGrid-1-(Max-Idx_Start[d]-Idx_Size[d]+1);
         Idx_MySize [d] = X_End - Idx_MyStart[d] + 1;
         Idx_MyStart[d] += MyRank_X[d]*NGrid;

         if (   (  d == 0  &&  ( OutputXYZ == 1 || OutputXYZ == 4 )  )  ||
                (  d == 1  &&  ( OutputXYZ == 2 || OutputXYZ == 5 )  )  ||
                (  d == 2  &&  ( OutputXYZ == 3 || OutputXYZ == 6 )  )     )     Idx_MySize[d] = 1;


//       verify the variables "Idx_MyStart and Idx_MySize"
         if ( Idx_MyStart[d] < Idx_Start[d]  ||  Idx_MyStart[d] >= Idx_Start[d]+Idx_Size[d] )
         {
            fprintf( stderr, "ERROR : incorrect Idx_MyStart[%d] (%d) --> out of range !!\n", d, Idx_MyStart[d] );
            MPI_Exit();
         }

         if ( Idx_MySize[d] < 0  ||  Idx_MySize[d] > NX0[d]*(1<<TargetLevel) )
         {
            fprintf( stderr, "ERROR : incorrect Idx_MySize[%d] (%d) --> out of range !!\n", d, Idx_MySize[d] );
            MPI_Exit();
         }

         if ( Idx_MyStart[d]+Idx_MySize[d] > Idx_Start[d]+Idx_Size[d] )
         {
            fprintf( stderr, "ERROR : incorrect Idx_MyStart[%d]+Idx_MySize[%d] (%d) --> out of range (%d) !!\n",
                     d, d, Idx_MyStart[d]+Idx_MySize[d], Idx_Start[d]+Idx_Size[d] );
            MPI_Exit();
         }

      } // if ( Min < Idx_Start[d]+Idx_Size[d]  &&  Max >= Idx_Start[d] )

      else
      {
         Idx_MyStart[d] = 0;
         Idx_MySize [d] = ( OutputXYZ == d+4 ) ? 1 : 0;  // set Idx_MySize[X] = 1 for the function "SumOverRanks"
      }
   } // for (int d=0; d<3; d++)


// record the target subdomain
   const double PhyCoord_HalfSize[3] = { 0.5*PhyCoord_Size[0], 0.5*PhyCoord_Size[1], 0.5*PhyCoord_Size[2] };
   for (int d=0; d<3; d++)
   {
      OutCoord_Start[d] = PhyCoord_Start[d] - ShiftScale[d]*amr.dh[NLEVEL-1];

      while ( OutCoord_Start[d]+PhyCoord_HalfSize[d] < 0.0  ||  OutCoord_Start[d]+PhyCoord_HalfSize[d] >= amr.BoxSize[d] )
      {
         if ( OutCoord_Start[d]+PhyCoord_HalfSize[d] <  0.0            )   OutCoord_Start[d] += amr.BoxSize[d];
         if ( OutCoord_Start[d]+PhyCoord_HalfSize[d] >= amr.BoxSize[d] )   OutCoord_Start[d] -= amr.BoxSize[d];
      }
   }


   if ( MyRank == 0 )   cout << "   Init_TargetDomain ... done" << endl;

} // FUNCTION : Init_TargetDomain



//-------------------------------------------------------------------------------------------------------
// Function    :  Init_Convert2Temp
// Description :  Set the conversion factor for calculating temperature
//
// Note        :  Must work with data providing both units and molecular weight (MU)
//-------------------------------------------------------------------------------------------------------
void Init_Convert2Temp()
{

   const double kB  = 1.38064852e-16;     // Boltzmann constant in erg/K
   const double amu = 1.660539040e-24;    // atomic mass unit in g

   if ( WithUnit  &&  MU > 0.0 )
   {
      if ( Unit_V <= 0.0 )
         Aux_Error( ERROR_INFO, "Unit_V = %14.7e <= 0.0 !!\n", Unit_V );

      Convert2Temp = Unit_V*Unit_V*amu*MU/kB;
   }

   else
      Aux_Error( ERROR_INFO, "cannot determine the conversion factor for calculating temperature (WithUnit %d, MU %14.7e) !!\n",
                 WithUnit, MU );

} // FUCNTION : Init_Convert2Temp



#ifndef SERIAL
//-------------------------------------------------------------------------------------------------------
// Function    :  Init_MPI
// Description :  Initialize MPI and parameters for parallelization
//-------------------------------------------------------------------------------------------------------
void Init_MPI( int *argc, char ***argv )
{

   if ( MPI_Init( argc, argv ) != MPI_SUCCESS )
   {
      cerr << "MPI_Init failed !!" << endl;
      exit( 1 );
   }

   if ( MPI_Comm_rank( MPI_COMM_WORLD, &MyRank ) != MPI_SUCCESS )
   {
      cerr << "MPI_Comm_rank failed !!" << endl;
      exit( 1 );
   }

   if ( MyRank == 0 )   cout << "Init_MPI ... done" << endl;

} // FUNCTION : Init_MPI
#endif



//-------------------------------------------------------------------------------------------------------
// Function    :  PreparePatch
// Description :  Prepare the central and buffer data of a single patch
//
// Note        :  If any of the sibling patch of PID does NOT exist, it will interpolate on the CData to fill up
//                the buffer data in FData
//
// Parameter   :  lv       : Refinement level of the targeted patch
//                PID      : Patch index of the targeted patch
//                Buffer   : Size of buffer (ghost zone)
//                FData    : Array to store the prepared data ( width = PATCH_SIZE + 2*Buffer )
//                CData    : Array for interpolation ( width = PATCH_SIZE + 2*Buffer )
//                BC_Face  : Priority of the BC. along different boundary faces (z>y>x)
//-------------------------------------------------------------------------------------------------------
void PreparePatch( const int lv, const int PID, const int Buffer, real FData[], real CData[], const int BC_Face[] )
{

   const int  Size     = PATCH_SIZE + 2*Buffer;
   const int  dii      = 1;
   const int  djj      = Size;
   const long dkk      = Size*Size;
   const long dvv      = Size*Size*Size;
   const int  Size3[3] = { Size, Size, Size };

   int  i0, j0, k0, ii0, jj0, kk0, i_loop, j_loop, k_loop, ID, CStart[3], CRange[3], FStart[3];
   int  SibPID, LocalID, BC_Sibling, BC_Idx_Start[3], BC_Idx_End[3], PotIdx, ParDensIdx, NextIdx;


// determine the array indices for outputting potential and particle densities
   NextIdx    = NCOMP_TOTAL;
   PotIdx     = -1;
   ParDensIdx = -1;

   if ( OutputPot     )    PotIdx     = NextIdx ++;
   if ( OutputParDens )    ParDensIdx = NextIdx ++;


// determine which variable needs to be monotonic
   const bool PhaseUnwrapping_Yes    = true;
   const bool PhaseUnwrapping_No     = false;
   const bool EnsureMonotonicity_Yes = true;
   const bool EnsureMonotonicity_No  = false;

   bool Monotonicity[NLoad];

   for (int v=0; v<NLoad; v++)
   {
//    we now apply monotonic interpolation to ALL fluid variables (which helps alleviate the issue of negative density/pressure)
#     if ( MODEL == HYDRO )
      /*
      if ( v == DENS  ||  v == ENGY  ||  v == ParDensIdx )  Monotonicity[v] = EnsureMonotonicity_Yes;
      else                                                  Monotonicity[v] = EnsureMonotonicity_No;
      */
                                                            Monotonicity[v] = EnsureMonotonicity_Yes;

#     elif ( MODEL == MHD )
#     warning : WAIT MHD !!!

#     elif ( MODEL == ELBDM )
      if ( v == DENS  ||  v == ParDensIdx  ||  ( v >= NCOMP_FLUID && v < NCOMP_TOTAL )  )
                                                            Monotonicity[v] = EnsureMonotonicity_Yes;
      else
                                                            Monotonicity[v] = EnsureMonotonicity_No;

#     else
                                                            Monotonicity[v] = EnsureMonotonicity_No;
#     warning : WARNING : DO YOU WANT TO ENSURE THE POSITIVITY OF INTERPOLATION ??
#     endif // MODEL
   }


// array for ELBDM_IntPhase
#  if ( MODEL == ELBDM )
   real *CData_Phase = ( ELBDM_IntPhase ) ? (new real [CUBE(Size)]) : NULL;

   bool GotPhaseAlready = false;
#  endif


// prepare the central data
#  ifdef GAMER_DEBUG
   if ( amr.patch[lv][PID]->fluid == NULL )
      Aux_Error( ERROR_INFO, "lv %d, PID %d, fluid array is NOT allocated !!\n", lv, PID );

   if ( OutputPot  &&  amr.patch[lv][PID]->pot == NULL )
      Aux_Error( ERROR_INFO, "lv %d, PID %d, potential array is NOT allocated !!\n", lv, PID );

   if ( OutputParDens  &&  amr.patch[lv][PID]->par_dens == NULL )
      Aux_Error( ERROR_INFO, "lv %d, PID %d, particle density array is NOT allocated !!\n", lv, PID );
#  endif

   for (int v=0; v<NCOMP_TOTAL; v++)
   for (int k=0, kk=Buffer; k<PATCH_SIZE; k++, kk++)
   for (int j=0, jj=Buffer; j<PATCH_SIZE; j++, jj++)
   for (int i=0, ii=Buffer; i<PATCH_SIZE; i++, ii++)
   {
      ID        = (long)v*dvv + kk*dkk + jj*djj + ii;
      FData[ID] = amr.patch[lv][PID]->fluid[v][k][j][i];
   }

   if ( OutputPot )
   for (int k=0, kk=Buffer; k<PATCH_SIZE; k++, kk++)
   for (int j=0, jj=Buffer; j<PATCH_SIZE; j++, jj++)
   for (int i=0, ii=Buffer; i<PATCH_SIZE; i++, ii++)
   {
      ID        = (long)PotIdx*dvv + kk*dkk + jj*djj + ii;
      FData[ID] = amr.patch[lv][PID]->pot[k][j][i];
   }

   if ( OutputParDens )
   for (int k=0, kk=Buffer; k<PATCH_SIZE; k++, kk++)
   for (int j=0, jj=Buffer; j<PATCH_SIZE; j++, jj++)
   for (int i=0, ii=Buffer; i<PATCH_SIZE; i++, ii++)
   {
      ID        = (long)ParDensIdx*dvv + kk*dkk + jj*djj + ii;
      FData[ID] = amr.patch[lv][PID]->par_dens[k][j][i];
   }


// prepare the buffer data
   for (int sib=0; sib<26; sib++)
   {
      SibPID = amr.patch[lv][PID]->sibling[sib];

//    (1) if the targeted sibling patch exists --> just copy data from the nearby patches at the same level
      if ( SibPID >= 0 )
      {
#        ifdef GAMER_DEBUG
         if ( amr.patch[lv][SibPID]->fluid == NULL )
            Aux_Error( ERROR_INFO, "lv %d, SibPID %d, fluid array is NOT allocated !!\n", lv, SibPID );

         if ( OutputPot  &&  amr.patch[lv][SibPID]->pot == NULL )
            Aux_Error( ERROR_INFO, "lv %d, SibPID %d, potential array is NOT allocated !!\n", lv, SibPID );

         if ( OutputParDens &&  amr.patch[lv][SibPID]->par_dens == NULL )
            Aux_Error( ERROR_INFO, "lv %d, SibPID %d, particle density array is NOT allocated !!\n", lv, SibPID );
#        endif

         i0     = TABLE_01( sib, 'x', PATCH_SIZE-Buffer, 0, 0 );
         j0     = TABLE_01( sib, 'y', PATCH_SIZE-Buffer, 0, 0 );
         k0     = TABLE_01( sib, 'z', PATCH_SIZE-Buffer, 0, 0 );
         ii0    = TABLE_01( sib, 'x', 0, Buffer, PATCH_SIZE+Buffer );
         jj0    = TABLE_01( sib, 'y', 0, Buffer, PATCH_SIZE+Buffer );
         kk0    = TABLE_01( sib, 'z', 0, Buffer, PATCH_SIZE+Buffer );
         i_loop = TABLE_01( sib, 'x', Buffer, PATCH_SIZE, Buffer );
         j_loop = TABLE_01( sib, 'y', Buffer, PATCH_SIZE, Buffer );
         k_loop = TABLE_01( sib, 'z', Buffer, PATCH_SIZE, Buffer );

         for (int v=0; v<NCOMP_TOTAL; v++)
         for (int k=k0, kk=kk0; k<k0+k_loop; k++, kk++)
         for (int j=j0, jj=jj0; j<j0+j_loop; j++, jj++)
         for (int i=i0, ii=ii0; i<i0+i_loop; i++, ii++)
         {
            ID        = (long)v*dvv + kk*dkk + jj*djj + ii;
            FData[ID] = amr.patch[lv][SibPID]->fluid[v][k][j][i];
         }

         if ( OutputPot )
         for (int k=k0, kk=kk0; k<k0+k_loop; k++, kk++)
         for (int j=j0, jj=jj0; j<j0+j_loop; j++, jj++)
         for (int i=i0, ii=ii0; i<i0+i_loop; i++, ii++)
         {
            ID        = (long)PotIdx*dvv + kk*dkk + jj*djj + ii;
            FData[ID] = amr.patch[lv][SibPID]->pot[k][j][i];
         }

         if ( OutputParDens )
         for (int k=k0, kk=kk0; k<k0+k_loop; k++, kk++)
         for (int j=j0, jj=jj0; j<j0+j_loop; j++, jj++)
         for (int i=i0, ii=ii0; i<i0+i_loop; i++, ii++)
         {
            ID        = (long)ParDensIdx*dvv + kk*dkk + jj*djj + ii;
            FData[ID] = amr.patch[lv][SibPID]->par_dens[k][j][i];
         }
      } // if ( SibPID >= 0 )



//    (2) if the targeted sibling patch does not exist --> interpolate from the coarse-grid array "CData"
      else if ( SibPID == -1 )
      {
//       the CData array must be properly prepared if any of the sibling patch does NOT exist
#        ifdef GAMER_DEBUG
         if ( CData == NULL )
            Aux_Error( ERROR_INFO, "CData == NULL, Rank = %d, lv = %d, PID = %d !!\n", MyRank, lv, PID );
#        endif


         LocalID = PID%8;

         for (int d=0; d<3; d++)
         {
            CStart[d] = TABLE_01( sib, 'x'+d, Buffer/2, Buffer, Buffer+PS1/2 ) + TABLE_02( LocalID, 'x'+d, 0, PS1/2 );
            CRange[d] = TABLE_01( sib, 'x'+d, Buffer/2, PS1/2, Buffer/2 );
            FStart[d] = TABLE_01( sib, 'x'+d, 0, Buffer, Buffer+PS1 );
         }


#        if ( MODEL == ELBDM )
         if ( ELBDM_IntPhase )
         {
            real *const CData_Dens = CData + DENS*CUBE(Size);
            real *const CData_Real = CData + REAL*CUBE(Size);
            real *const CData_Imag = CData + IMAG*CUBE(Size);
            real *const FData_Dens = FData + DENS*CUBE(Size);
            real *const FData_Real = FData + REAL*CUBE(Size);
            real *const FData_Imag = FData + IMAG*CUBE(Size);

//          get the wrapped phase
            if ( !GotPhaseAlready )
            {
               for (int t=0; t<CUBE(Size); t++)    CData_Phase[t] = ATAN2( CData_Imag[t], CData_Real[t] );

               GotPhaseAlready = true;
            }

//          interpolate density
            Interpolate( CData_Dens,  Size3, CStart, CRange, FData_Dens, Size3, FStart, 1, IntScheme,
                         PhaseUnwrapping_No, EnsureMonotonicity_Yes, Int_MonoCoeff );

//          interpolate phase (store in the REAL component)
            Interpolate( CData_Phase, Size3, CStart, CRange, FData_Real, Size3, FStart, 1, IntScheme,
                         PhaseUnwrapping_Yes, EnsureMonotonicity_No, Int_MonoCoeff );

//          retrieve real and imaginary parts
            real Amp, Phase, Rho;

            for (int k=FStart[2]; k<FStart[2]+2*CRange[2]; k++)
            for (int j=FStart[1]; j<FStart[1]+2*CRange[1]; j++)
            for (int i=FStart[0]; i<FStart[0]+2*CRange[0]; i++)
            {
               ID  = (long)k*dkk + j*djj + i;

               Phase = FData_Real[ID];
               Rho   = FData_Dens[ID];

//             be careful about the negative density introduced from the round-off errors
               if ( Rho < (real)0.0 )
               {
                  FData_Dens[ID] = (real)0.0;
                  Rho            = (real)0.0;
               }

               Amp            = SQRT( Rho );
               FData_Real[ID] = Amp*COS( Phase );
               FData_Imag[ID] = Amp*SIN( Phase );
            } // i,j,k


//          interpolate passive scalars
            for (int v=NCOMP_FLUID; v<NCOMP_TOTAL; v++)
            Interpolate( CData+v*CUBE(Size), Size3, CStart, CRange, FData+v*CUBE(Size), Size3, FStart, 1,
                         IntScheme, PhaseUnwrapping_No, Monotonicity[v], Int_MonoCoeff );


//          interpolate potential
            if ( OutputPot )
            Interpolate( CData+PotIdx*CUBE(Size), Size3, CStart, CRange, FData+PotIdx*CUBE(Size), Size3, FStart, 1,
                         IntScheme, PhaseUnwrapping_No, EnsureMonotonicity_No, Int_MonoCoeff );


//          interpolate particle density
            if ( OutputParDens  )
            Interpolate( CData+ParDensIdx*CUBE(Size), Size3, CStart, CRange, FData+ParDensIdx*CUBE(Size), Size3, FStart, 1,
                         IntScheme, PhaseUnwrapping_No, EnsureMonotonicity_Yes, Int_MonoCoeff );

         } // if ( ELBDM_IntPhase )

         else // if ( ELBDM_IntPhase ) ... else ...
#        endif // #if ( MODEL == ELBDM )
         {
            for (int v=0; v<NLoad; v++)
               Interpolate( CData+v*CUBE(Size), Size3, CStart, CRange, FData+v*CUBE(Size), Size3, FStart, 1,
                            IntScheme, PhaseUnwrapping_No, Monotonicity[v], Int_MonoCoeff );
         }
      } // else if ( SibPID == -1 ) ...



//    (3) if the targeted sibling patch lies outside the simulation domain --> apply the specified B.C.
      else if ( SibPID <= SIB_OFFSET_NONPERIODIC )
      {
         for (int d=0; d<3; d++)
         {
            BC_Idx_Start[d] = TABLE_01( sib, 'x'+d, 0, Buffer, Buffer+PATCH_SIZE );
            BC_Idx_End  [d] = TABLE_01( sib, 'x'+d, Buffer, PATCH_SIZE, Buffer ) + BC_Idx_Start[d] - 1;
         }

         BC_Sibling = SIB_OFFSET_NONPERIODIC - SibPID;

         switch ( ExtBC )
         {
            case 1:
               Hydro_BoundaryCondition_Outflow( FData, BC_Face[BC_Sibling], NLoad, Buffer,
                                                Size, Size, Size, BC_Idx_Start, BC_Idx_End );
            break;

            default:
               Aux_Error( ERROR_INFO, "unsupported ExtBC (%d) !!\n", ExtBC );

         } // switch ( FluBC[ BC_Face[BC_Sibling] ] )
      } // else if ( SibPID <= SIB_OFFSET_NONPERIODIC )

      else
         Aux_Error( ERROR_INFO, "SibPID == %d (PID %d, Sib %d) !!\n", SibPID, PID, sib );

   } // for (int sib=0; sib<26; sib++)


// free memory
#  if ( MODEL == ELBDM )
   if ( CData_Phase != NULL )    delete [] CData_Phase;
#  endif

} // FUNCTION : PreparePatch




//-------------------------------------------------------------------------------------------------------
// Function    :  StoreData
// Description :  Store data in the array "Out"
//
// Note        :  Projection is also performed in this function
//
// Parameter   :  lv       : Refinement level of the input patch
//                PID      : Pathc index of the input patch
//                FData    : Array storing the data to be outputted
//                Buffer   : Size of buffer
//                Out      : Array to store the results
//-------------------------------------------------------------------------------------------------------
void StoreData( const int lv, const int PID, real FData[], const int Buffer, real *Out )
{

   const long Size1v    = (long)Idx_MySize[0]*Idx_MySize[1]*Idx_MySize[2];    // total array size of one component
   const int  PatchSize = PATCH_SIZE*( 1<<(TargetLevel-lv) );
   const int  FSize     = PatchSize + 2*Buffer;
   const int  Corner[3] = { amr.patch[lv][PID]->corner[0]/amr.scale[TargetLevel],
                            amr.patch[lv][PID]->corner[1]/amr.scale[TargetLevel],
                            amr.patch[lv][PID]->corner[2]/amr.scale[TargetLevel]  };

   int    ijk_min[3], ijk_max[3], NSum_local, ProjDir, ii, jj, kk, Stride[3]={0,0,0};
   long   Jump, ID1, ID2;
   double SumData, NAve;


// calculate the index range
   for (int d=0; d<3; d++)
   {
      if ( Corner[d] < Idx_Start[d]+Idx_Size[d]  &&  Corner[d]+PatchSize > Idx_Start[d] )
      {
         ijk_min[d] = ( Idx_Start[d] > Corner[d] ) ? Idx_Start[d]-Corner[d]+Buffer : Buffer;
         ijk_max[d] = ( Idx_Start[d]+Idx_Size[d] >= Corner[d]+PatchSize ) ?
            Buffer+PatchSize-1 : Buffer+PatchSize-1-(Corner[d]+PatchSize-Idx_Start[d]-Idx_Size[d]);
      }

      else
         return;
   }


// projection
   if ( OutputXYZ == 4  ||  OutputXYZ == 5  ||  OutputXYZ == 6 )
   {
      ProjDir          = OutputXYZ - 4;
      NSum_local       = ijk_max[ProjDir] - ijk_min[ProjDir] + 1;
      ijk_max[ProjDir] = ijk_min[ProjDir];
      NAve             = (double)Idx_Size[ProjDir];
      Jump             = 1;

      for (int t=0; t<ProjDir; t++)    Jump *= FSize;

      for (int v=0; v<NOut; v++)
      for (int k=ijk_min[2]; k<=ijk_max[2]; k++)
      for (int j=ijk_min[1]; j<=ijk_max[1]; j++)
      for (int i=ijk_min[0]; i<=ijk_max[0]; i++)
      {
         SumData = 0.0;
         ID1     = (((long)v*FSize + k)*FSize + j)*FSize + i;

         for (int t=0; t<NSum_local; t++)    SumData += FData[ ID1 + (long)t*Jump ];

         FData[ID1] = SumData / NAve;
      }
   }


// store data in the OutputArray array
   switch ( OutputXYZ )
   {
      case 1 : case 4 :    Stride[0] = 0;    Stride[1] = 1;             Stride[2] = Idx_MySize[1];
                           break;

      case 2 : case 5 :    Stride[0] = 1;    Stride[1] = 0;             Stride[2] = Idx_MySize[0];
                           break;

      case 3 : case 6 :    Stride[0] = 1;    Stride[1] = Idx_MySize[0];  Stride[2] = 0;
                           break;

      case 7 :             Stride[0] = 1;    Stride[1] = Idx_MySize[0];  Stride[2] = Idx_MySize[1]*Idx_MySize[0];
                           break;
   }

   for (int v=0; v<NOut; v++)                   {
   for (int k=ijk_min[2]; k<=ijk_max[2]; k++)   {  kk  = k - Buffer + Corner[2] - Idx_MyStart[2];
   for (int j=ijk_min[1]; j<=ijk_max[1]; j++)   {  jj  = j - Buffer + Corner[1] - Idx_MyStart[1];
   for (int i=ijk_min[0]; i<=ijk_max[0]; i++)   {  ii  = i - Buffer + Corner[0] - Idx_MyStart[0];
                                                   ID1 = (((long)v*FSize + k)*FSize + j)*FSize + i;
                                                   ID2 = (long)v*Size1v + (long)kk*Stride[2] + (long)jj*Stride[1] + (long)ii*Stride[0];

      Out[ID2] += FData[ID1];

   }}}}

} // FUNCTION : StoreData



#if ( MODEL == HYDRO )
//-------------------------------------------------------------------------------------------------------
// Function    :  GetPres
// Description :  Evaluate pressure
//
// Parameter   :  FData    : Array to store the output data
//                FSize    : Width of the FData array
//                Buffer   : Size of buffer
//                NextIdx  : Index to store the evaluated results
//-------------------------------------------------------------------------------------------------------
void GetPres( real FData[], const int FSize, const int Buffer, const int NextIdx )
{

   const int  di       = 1;
   const int  dj       = FSize;
   const long dk       = FSize*FSize;
   const long dv       = FSize*FSize*FSize;
   const real Gamma_m1 = GAMMA - 1.0;

   const real *Dens = FData +       0*dv;
   const real *MomX = FData +       1*dv;
   const real *MomY = FData +       2*dv;
   const real *MomZ = FData +       3*dv;
   const real *Engy = FData +       4*dv;
         real *Out  = FData + NextIdx*dv;

   long Idx;


//###OPTIMIZATION : perform calculations only for cells lying on the targeted slice
   for (int k=Buffer; k<FSize-Buffer; k++)
   for (int j=Buffer; j<FSize-Buffer; j++)
   for (int i=Buffer; i<FSize-Buffer; i++)
   {
      Idx      = (long)k*dk + j*dj + i*di;
      Out[Idx] = (  Engy[Idx] - 0.5*( MomX[Idx]*MomX[Idx] + MomY[Idx]*MomY[Idx] + MomZ[Idx]*MomZ[Idx] )/Dens[Idx] )*Gamma_m1;
   }

} // FUNCTION : GetPres



//-------------------------------------------------------------------------------------------------------
// Function    :  GetTemp
// Description :  Evaluate temperature
//
// Parameter   :  FData    : Array to store the output data
//                FSize    : Width of the FData array
//                Buffer   : Size of buffer
//                NextIdx  : Index to store the evaluated results
//                Con2T    : Factor for converting pressure/density to temperature
//                           --> Temperature is defined as "Con2T*Pressure/Density"
//-------------------------------------------------------------------------------------------------------
void GetTemp( real FData[], const int FSize, const int Buffer, const int NextIdx, const real Con2T )
{

   const int  di       = 1;
   const int  dj       = FSize;
   const long dk       = FSize*FSize;
   const long dv       = FSize*FSize*FSize;
   const real Gamma_m1 = GAMMA - 1.0;

   const real *Dens = FData +       0*dv;
   const real *MomX = FData +       1*dv;
   const real *MomY = FData +       2*dv;
   const real *MomZ = FData +       3*dv;
   const real *Engy = FData +       4*dv;
         real *Out  = FData + NextIdx*dv;

   real Pres, Temp;
   long Idx;


//###OPTIMIZATION : perform calculations only for cells lying on the targeted slice
   for (int k=Buffer; k<FSize-Buffer; k++)
   for (int j=Buffer; j<FSize-Buffer; j++)
   for (int i=Buffer; i<FSize-Buffer; i++)
   {
      Idx      = (long)k*dk + j*dj + i*di;
      Pres     = (  Engy[Idx] - 0.5*( MomX[Idx]*MomX[Idx] + MomY[Idx]*MomY[Idx] + MomZ[Idx]*MomZ[Idx] )/Dens[Idx] )*Gamma_m1;
      Temp     = Con2T*Pres/Dens[Idx];
      Out[Idx] = Temp;
   }

} // FUNCTION : GetTemp



//-------------------------------------------------------------------------------------------------------
// Function    :  GetCurlVel
// Description :  Evaluate the curl( velocity ) == vorticity
//
// Parameter   :  FData    : Array to store the output data
//                FSize    : Width of the FData array
//                Buffer   : Size of buffer
//                NextIdx  : Index to store the evaluated results
//                dh       : Cell size
//                Com2Phy  : Convert peculiar velocity from comoving to physical coords. (in km/s)
//                           --> by multiplying 100/ScaleA
//-------------------------------------------------------------------------------------------------------
void GetCurlVel( real FData[], const int FSize, const int Buffer, const int NextIdx, const real dh,
                 const bool Com2Phy )
{

   const int  di    = 1;
   const int  dj    = FSize;
   const long dk    = FSize*FSize;
   const long dv    = FSize*FSize*FSize;
   const real Coeff = (Com2Phy) ? 0.5/dh*100.0/Time[0] : 0.5/dh;

   const real *Rho  = FData +           0*dv;
   const real *MomX = FData +           1*dv;
   const real *MomY = FData +           2*dv;
   const real *MomZ = FData +           3*dv;
         real *Wx   = FData + (NextIdx+0)*dv;
         real *Wy   = FData + (NextIdx+1)*dv;
         real *Wz   = FData + (NextIdx+2)*dv;

   real *Vx = new real [dv];
   real *Vy = new real [dv];
   real *Vz = new real [dv];

   real _Rho, PxVy, PxVz, PyVx, PyVz, PzVx, PzVy;
   long Idx;


// get velocity
   for (int t=0; t<dv; t++)
   {
      _Rho  = 1.0/Rho[t];
      Vx[t] = _Rho*MomX[t];
      Vy[t] = _Rho*MomY[t];
      Vz[t] = _Rho*MomZ[t];
   }


//###OPTIMIZATION : perform calculations only for cells lying on the targeted slice
   for (int k=Buffer; k<FSize-Buffer; k++)
   for (int j=Buffer; j<FSize-Buffer; j++)
   for (int i=Buffer; i<FSize-Buffer; i++)
   {
      Idx  = (long)k*dk + j*dj + i*di;
      PxVy = Vy[Idx+di] - Vy[Idx-di];
      PxVz = Vz[Idx+di] - Vz[Idx-di];
      PyVx = Vx[Idx+dj] - Vx[Idx-dj];
      PyVz = Vz[Idx+dj] - Vz[Idx-dj];
      PzVx = Vx[Idx+dk] - Vx[Idx-dk];
      PzVy = Vy[Idx+dk] - Vy[Idx-dk];

      Wx[Idx] = Coeff*( PyVz - PzVy );
      Wy[Idx] = Coeff*( PzVx - PxVz );
      Wz[Idx] = Coeff*( PxVy - PyVx );
   }


   delete [] Vx;
   delete [] Vy;
   delete [] Vz;

} // FUNCTION : GetCurlVel



//-------------------------------------------------------------------------------------------------------
// Function    :  GetDivVel
// Description :  Evaluate the divergence( velocity )
//
// Parameter   :  FData    : Array to store the output data
//                FSize    : Width of the FData array
//                Buffer   : Size of buffer
//                NextIdx  : Index to store the evaluated results
//                dh       : Cell size
//                Com2Phy  : Convert peculiar velocity from comoving to physical coords. (in km/s)
//                           --> by multiplying 100/ScaleA
//-------------------------------------------------------------------------------------------------------
void GetDivVel( real FData[], const int FSize, const int Buffer, const int NextIdx, const real dh,
                const bool Com2Phy )
{

   const int  di    = 1;
   const int  dj    = FSize;
   const long dk    = FSize*FSize;
   const long dv    = FSize*FSize*FSize;
   const real Coeff = (Com2Phy) ? 0.5/dh*100.0/Time[0] : 0.5/dh;

   const real *Rho  = FData +       0*dv;
   const real *MomX = FData +       1*dv;
   const real *MomY = FData +       2*dv;
   const real *MomZ = FData +       3*dv;
         real *Out  = FData + NextIdx*dv;

   real *Vx = new real [dv];
   real *Vy = new real [dv];
   real *Vz = new real [dv];

   real _Rho;
   long Idx;


// get velocity
   for (int t=0; t<dv; t++)
   {
      _Rho  = 1.0/Rho[t];
      Vx[t] = _Rho*MomX[t];
      Vy[t] = _Rho*MomY[t];
      Vz[t] = _Rho*MomZ[t];
   }


//###OPTIMIZATION : perform calculations only for cells lying on the targeted slice
   for (int k=Buffer; k<FSize-Buffer; k++)
   for (int j=Buffer; j<FSize-Buffer; j++)
   for (int i=Buffer; i<FSize-Buffer; i++)
   {
      Idx      = (long)k*dk + j*dj + i*di;
      Out[Idx] = Coeff*( ( Vz[Idx+dk] - Vz[Idx-dk] ) + ( Vy[Idx+dj] - Vy[Idx-dj] ) + ( Vx[Idx+di] - Vx[Idx-di] ) );
   }


   delete [] Vx;
   delete [] Vy;
   delete [] Vz;

} // FUNCTION : GetDivVel
#endif // #if ( MODEL == HYDRO )



#if ( MODEL == ELBDM )
//-------------------------------------------------------------------------------------------------------
// Function    :  GetELBDM_Vel
// Description :  Evaluate the velocity field in ELBDM
//
// Parameter   :  FData    : Array to store the output data
//                FSize    : Width of the FData array
//                Buffer   : Size of buffer
//                NextIdx  : Index to store the evaluated results
//                dh       : Cell size
//                Com2Phy  : Convert peculiar velocity from comoving to physical coords. (in km/s)
//                           --> by multiplying 100/ScaleA
//-------------------------------------------------------------------------------------------------------
void GetELBDM_Vel( real FData[], const int FSize, const int Buffer, const int NextIdx, const real dh,
                   const bool Com2Phy )
{

   const int   di   = 1;
   const int   dj   = FSize;
   const long  dk   = FSize*FSize;
   const long  dv   = FSize*FSize*FSize;
   const real _dh2  = 0.5/dh;
   const real _Eta  = 1.0/ELBDM_ETA;
   const real Coeff = (Com2Phy) ? _Eta*100.0/Time[0] : _Eta;

   const real *Dens = FData +        DENS*dv;
   const real *Real = FData +        REAL*dv;
   const real *Imag = FData +        IMAG*dv;
         real *Vx   = FData + (NextIdx+0)*dv;
         real *Vy   = FData + (NextIdx+1)*dv;
         real *Vz   = FData + (NextIdx+2)*dv;

   real GradR[3], GradI[3], _Dens;
   long Idx;


// evaluate the ELBDM velocity
//###OPTIMIZATION : perform calculations only for cells lying on the targeted slice
   for (int k=Buffer; k<FSize-Buffer; k++)
   for (int j=Buffer; j<FSize-Buffer; j++)
   for (int i=Buffer; i<FSize-Buffer; i++)
   {
      Idx  = (long)k*dk + j*dj + i*di;

      GradR[0] = _dh2*( Real[Idx+di] - Real[Idx-di] );
      GradI[0] = _dh2*( Imag[Idx+di] - Imag[Idx-di] );
      GradR[1] = _dh2*( Real[Idx+dj] - Real[Idx-dj] );
      GradI[1] = _dh2*( Imag[Idx+dj] - Imag[Idx-dj] );
      GradR[2] = _dh2*( Real[Idx+dk] - Real[Idx-dk] );
      GradI[2] = _dh2*( Imag[Idx+dk] - Imag[Idx-dk] );

      _Dens    = 1.0 / Dens[Idx];
      Vx[Idx] = Coeff*_Dens*( Real[Idx]*GradI[0] - Imag[Idx]*GradR[0] );
      Vy[Idx] = Coeff*_Dens*( Real[Idx]*GradI[1] - Imag[Idx]*GradR[1] );
      Vz[Idx] = Coeff*_Dens*( Real[Idx]*GradI[2] - Imag[Idx]*GradR[2] );
   }

} // FUNCTION : GetELBDM_Vel
#endif // #if ( MODEL == ELBDM )



//-------------------------------------------------------------------------------------------------------
// Function    :  Refine2TargetLevel
// Description :  Refine to the level "TargetLevel" patch by patch
//-------------------------------------------------------------------------------------------------------
void Refine2TargetLevel()
{

   if ( MyRank == 0 )   cout << "Refine2TargetLevel ... " << endl;


   const int  Buffer         = BufSize;
   const int  CStart[3]      = { Buffer/2, Buffer/2, Buffer/2 };
   const int  FStart[3]      = { 0, 0, 0 };
#  if ( MODEL != ELBDM )
   const bool ELBDM_IntPhase = false;
#  endif

   int   FaPID, NextIdx, CSize3[3], CRange[3], FSize3[3];
   long  CSize, FSize;
   real *CData = NULL, *FData = NULL;

#  ifdef OPENMP
   const long OutSize = (long)NOut*Idx_MySize[0]*Idx_MySize[1]*Idx_MySize[2];
   real **OutputArray_OMP = NULL;
   int TID;

   if ( OutputXYZ == 4  ||  OutputXYZ == 5  ||  OutputXYZ == 6 )
   {
      OutputArray_OMP = new real* [OMP_NThread];

      for (int t=0; t<OMP_NThread; t++)
      {
         OutputArray_OMP[t] = new real [OutSize];

         for (long i=0; i<OutSize; i++)   OutputArray_OMP[t][i] = 0.0;
      }
   }
#  endif // #ifdef OPENMP


// determine the array indices for outputting potential and particle densities
   int PotIdx, ParDensIdx;

   NextIdx    = NCOMP_TOTAL;
   PotIdx     = -1;
   ParDensIdx = -1;

   if ( OutputPot     )    PotIdx     = NextIdx ++;
   if ( OutputParDens )    ParDensIdx = NextIdx ++;


// determine the priority of different boundary faces (z>y>x) to set the corner cells properly for the non-periodic B.C.
   int BC_Face[26], BC_Face_tmp[3];

   for (int s=0; s<26; s++)
   {
      BC_Face_tmp[0] = TABLE_01( s, 'x', 0, -1, 1 );
      BC_Face_tmp[1] = TABLE_01( s, 'y', 2, -1, 3 );
      BC_Face_tmp[2] = TABLE_01( s, 'z', 4, -1, 5 );

//    z > y > x
      if      ( BC_Face_tmp[2] != -1 )   BC_Face[s] = BC_Face_tmp[2];
      else if ( BC_Face_tmp[1] != -1 )   BC_Face[s] = BC_Face_tmp[1];
      else if ( BC_Face_tmp[0] != -1 )   BC_Face[s] = BC_Face_tmp[0];
   }


// determine which variable needs to be monotonic
   const bool PhaseUnwrapping_Yes    = true;
   const bool PhaseUnwrapping_No     = false;
   const bool EnsureMonotonicity_Yes = true;
   const bool EnsureMonotonicity_No  = false;

   bool Monotonicity[NLoad];

   for (int v=0; v<NLoad; v++)
   {
//    we now apply monotonic interpolation to ALL fluid variables (which helps alleviate the issue of negative density/pressure)
#     if ( MODEL == HYDRO )
      /*
      if ( v == DENS  ||  v == ENGY  ||  v == ParDensIdx )  Monotonicity[v] = EnsureMonotonicity_Yes;
      else                                                  Monotonicity[v] = EnsureMonotonicity_No;
      */
                                                            Monotonicity[v] = EnsureMonotonicity_Yes;

#     elif ( MODEL == MHD )
#     warning : WAIT MHD !!!

#     elif ( MODEL == ELBDM )
      if ( v == DENS  ||  v == ParDensIdx  ||  ( v >= NCOMP_FLUID && v < NCOMP_TOTAL )  )
                                                            Monotonicity[v] = EnsureMonotonicity_Yes;
      else
                                                            Monotonicity[v] = EnsureMonotonicity_No;

#     else
                                                            Monotonicity[v] = EnsureMonotonicity_No;
#     warning : WARNING : DO YOU WANT TO ENSURE THE POSITIVITY OF INTERPOLATION ??
#     endif // MODEL
   }


// loop over all levels
   for (int lv=0; lv<=TargetLevel; lv++)
   {
      if ( MyRank == 0 )   cout << "  Level = " << lv << " ... " << flush;

#     pragma omp parallel private( FaPID, FSize, CSize, NextIdx, CData, FData, CSize3, CRange, FSize3, TID )
      {
#        ifdef OPENMP
         TID = omp_get_thread_num();
#        endif

#        pragma omp for
         for (int PID=0; PID<NPatchComma[lv][1]; PID++)
         {
//          find the target patches
            if (  ( amr.patch[lv][PID]->son == -1 || lv == TargetLevel )  &&
                  WithinCandidateBox( amr.patch[lv][PID]->corner, PATCH_SIZE*amr.scale[lv], 0 )  )
            {
//             allocate memory
               CSize = PATCH_SIZE + 2*Buffer;
               FSize = PATCH_SIZE + 2*Buffer;
               CData = new real [ (long)NOut*CUBE(CSize) ];
               FData = new real [ (long)NOut*CUBE(FSize) ];


//###OPTIMIZATION : prepare the coarse-grid data "only once" for all son patches
//             prepare the coarse-grid data (level = lv-1)
               if ( lv != 0 )
               {
                  FaPID = amr.patch[lv][PID]->father;

                  PreparePatch( lv-1, FaPID, Buffer, CData, NULL, BC_Face );
               }


//             prepare the fine-grid data (level = lv)
               PreparePatch( lv, PID, Buffer, FData, CData, BC_Face );


//             perform interpolation to reach the TagetLevel resolution
               for (int FineLv=lv+1; FineLv<=TargetLevel; FineLv++)
               {
//                reset and reallocate CData and FData pointers
                  delete [] CData;
                  CSize = FSize;
                  CData = FData;
                  FSize = ( FSize - 2*Buffer )*2 + 2*Buffer;
                  FData = new real [ (long)NOut*CUBE(FSize) ];

//                interpolate from CData (level = FineLv-1) to FData (level = FineLv)
                  for (int d=0; d<3; d++)
                  {
                     CSize3[d] = CSize;
                     CRange[d] = CSize - Buffer;
                     FSize3[d] = FSize;
                  }

#                 if ( MODEL == ELBDM )
                  if ( ELBDM_IntPhase )
                  {
                     real *const CData_Dens = CData + DENS*CUBE(CSize);
                     real *const CData_Real = CData + REAL*CUBE(CSize);
                     real *const CData_Imag = CData + IMAG*CUBE(CSize);
                     real *const FData_Dens = FData + DENS*CUBE(FSize);
                     real *const FData_Real = FData + REAL*CUBE(FSize);
                     real *const FData_Imag = FData + IMAG*CUBE(FSize);

//                   get the wrapped phase (store in the REAL component)
                     for (long t=0; t<CUBE(CSize); t++)  CData_Real[t] = ATAN2( CData_Imag[t], CData_Real[t] );

//                   interpolate density
                     Interpolate( CData_Dens, CSize3, CStart, CRange, FData_Dens, FSize3, FStart, 1, IntScheme,
                                  PhaseUnwrapping_No, EnsureMonotonicity_Yes, Int_MonoCoeff );

//                   interpolate phase
                     Interpolate( CData_Real, CSize3, CStart, CRange, FData_Real, FSize3, FStart, 1, IntScheme,
                                  PhaseUnwrapping_Yes, EnsureMonotonicity_No, Int_MonoCoeff );

//                   retrieve real and imaginary parts
                     real Amp, Phase, Rho;

                     for (long t=0; t<CUBE(FSize); t++)
                     {
                        Phase = FData_Real[t];
                        Rho   = FData_Dens[t];

//                      be careful about the negative density introduced from the round-off errors
                        if ( Rho < (real)0.0 )
                        {
                           FData_Dens[t] = (real)0.0;
                           Rho           = (real)0.0;
                        }

                        Amp           = SQRT( Rho );
                        FData_Real[t] = Amp*COS( Phase );
                        FData_Imag[t] = Amp*SIN( Phase );
                     } // for (int t=0; t<CUBE(FSize); t++)

//                   interpolate passive scalars
                     for (int v=NCOMP_FLUID; v<NCOMP_TOTAL; v++)
                     Interpolate( CData+v*CUBE(CSize), CSize3, CStart, CRange, FData+v*CUBE(FSize),
                                  FSize3, FStart, 1, IntScheme, PhaseUnwrapping_No, Monotonicity[v], Int_MonoCoeff );
//                   interpolate potential
                     if ( OutputPot )
                     Interpolate( CData+PotIdx*CUBE(CSize), CSize3, CStart, CRange, FData+PotIdx*CUBE(FSize),
                                  FSize3, FStart, 1, IntScheme, PhaseUnwrapping_No, EnsureMonotonicity_No, Int_MonoCoeff );

//                   interpolate potential
                     if ( OutputParDens )
                     Interpolate( CData+ParDensIdx*CUBE(CSize), CSize3, CStart, CRange, FData+ParDensIdx*CUBE(FSize),
                                  FSize3, FStart, 1, IntScheme, PhaseUnwrapping_No, EnsureMonotonicity_Yes, Int_MonoCoeff );

                  } // if ( ELBDM_IntPhase )

                  else // if ( ELBDM_IntPhase ) ... else ...
#                 endif // #if ( MODEL == ELBDM )
                  {
                     for (int v=0; v<NLoad; v++)
                        Interpolate( CData+v*CUBE(CSize), CSize3, CStart, CRange, FData+v*CUBE(FSize), FSize3, FStart, 1,
                                     IntScheme, PhaseUnwrapping_No, Monotonicity[v], Int_MonoCoeff );
                  }
               } // for (int FineLv=lv+1; FineLv<=TargetLevel; FineLv++)


               NextIdx = NLoad;

#              if ( MODEL == HYDRO )
//             calculate pressure
               if ( OutputPres )
               {
                  GetPres( FData, FSize, Buffer, NextIdx );
                  NextIdx ++;
               }


//             calculate temperature
               if ( OutputTemp )
               {
                  GetTemp( FData, FSize, Buffer, NextIdx, Convert2Temp );
                  NextIdx ++;
               }


//             calculate divergence(vel)
               if ( OutputDivVel )
               {
                  GetDivVel( FData, FSize, Buffer, NextIdx, amr.dh[TargetLevel], ConvertCom2PhyV );
                  NextIdx ++;
               }


//             calculate curl(vel)
               if ( OutputCurlVel )
               {
                  GetCurlVel( FData, FSize, Buffer, NextIdx, amr.dh[TargetLevel], ConvertCom2PhyV );
                  NextIdx += 3;
               }
#              endif


#              if ( MODEL == ELBDM )
//             calculate ELBDM velocity field
               if ( OutputELBDM_Vel )
               {
                  GetELBDM_Vel( FData, FSize, Buffer, NextIdx, amr.dh[TargetLevel], ConvertCom2PhyV );
                  NextIdx += 3;
               }
#              endif


//             store the final result
#              ifdef OPENMP
               if ( OutputXYZ == 4  ||  OutputXYZ == 5  ||  OutputXYZ == 6 )
               StoreData( lv, PID, FData, Buffer, OutputArray_OMP[TID] );
               else
               StoreData( lv, PID, FData, Buffer, OutputArray );

#              else
               StoreData( lv, PID, FData, Buffer, OutputArray );
#              endif


//             free memory
               delete [] CData;
               delete [] FData;
               CData = NULL;
               FData = NULL;

            } // if ( amr.patch[lv][PID]->son == -1 ) ...
         } // for (int PID=0; PID<NPatchComma[lv][1]; PID++)
      } // OpenMP parallel region

      MPI_Barrier( MPI_COMM_WORLD );
      if ( MyRank == 0 )   cout << "done" << endl;

   } // for (int lv=0; lv<=TargetLevel; lv++)


#  ifdef OPENMP
   if ( OutputXYZ == 4  ||  OutputXYZ == 5  ||  OutputXYZ == 6 )
   {
      for (int t=0; t<OMP_NThread; t++)
      {
         for (long i=0; i<OutSize; i++)   OutputArray[i] += OutputArray_OMP[t][i];

         delete [] OutputArray_OMP[t];
      }

      delete [] OutputArray_OMP;
   }
#  endif // #ifdef OPENMP


   if ( MyRank == 0 )   cout << "Refine2TargetLevel ... done" << endl;

} // FUNCTION : Refine2TargetLevel



//-------------------------------------------------------------------------------------------------------
// Function    :  AllocateOutputArray
// Description :  Allocate memory for the array "OutputArray"
//-------------------------------------------------------------------------------------------------------
void AllocateOutputArray()
{

   const long Size = (long)NOut*Idx_MySize[0]*Idx_MySize[1]*Idx_MySize[2];

   OutputArray = new real [Size];

// initialize it as zero (necessary for the projection operation)
   for (long t=0; t<Size; t++)    OutputArray[t] = 0.0;

} // FUNCTION : AllocateOutputArray



#ifndef SERIAL
//-------------------------------------------------------------------------------------------------------
// Function    :  SumOverRanks
// Description :  When doing projection, this function will collect the projection results from different ranks
//-------------------------------------------------------------------------------------------------------
void SumOverRanks()
{

   if ( MyRank == 0 )   cout << "SumOverRanks ... " << flush;


   const long Size = (long)NOut*Idx_MySize[0]*Idx_MySize[1]*Idx_MySize[2];

   int  RecvRank=0, NSend=0;
   int  *SendRank=NULL;
   real *RecvBuffer=NULL;

   switch ( OutputXYZ )
   {
      case 4:     NSend    = NGPU_X[0]-1;
                  RecvRank = MyRank_X[2]*NGPU_X[1]*NGPU_X[0] + MyRank_X[1]*NGPU_X[0];
                  SendRank = new int [NSend];

                  for (int ID=0; ID<NSend; ID++)
                     SendRank[ID] = MyRank_X[2]*NGPU_X[1]*NGPU_X[0] + MyRank_X[1]*NGPU_X[0] + (ID+1);

                  break;

      case 5:     NSend    = NGPU_X[1]-1;
                  RecvRank = MyRank_X[2]*NGPU_X[1]*NGPU_X[0] + MyRank_X[0];
                  SendRank = new int [NSend];

                  for (int ID=0; ID<NSend; ID++)
                     SendRank[ID] = MyRank_X[2]*NGPU_X[1]*NGPU_X[0] + (ID+1)*NGPU_X[0] + MyRank_X[0];

                  break;

      case 6:     NSend    = NGPU_X[2]-1;
                  RecvRank = MyRank_X[1]*NGPU_X[0] + MyRank_X[0];
                  SendRank = new int [NSend];

                  for (int ID=0; ID<NSend; ID++)
                     SendRank[ID] = (ID+1)*NGPU_X[1]*NGPU_X[0] + MyRank_X[1]*NGPU_X[0] + MyRank_X[0];

                  break;

      default :
                  cerr << "ERROR : incorrect OutputXYZ in the function \"SumOverRanks\" !!" << endl;
                  MPI_Exit();
   }


   RecvBuffer = new real [Size];


   for (int SendID=0; SendID<NSend; SendID++)
   {
      if ( MyRank == RecvRank )
      {
#        ifdef FLOAT8
         MPI_Recv( RecvBuffer, Size, MPI_DOUBLE, SendRank[SendID], 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE );
#        else
         MPI_Recv( RecvBuffer, Size, MPI_FLOAT,  SendRank[SendID], 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE );
#        endif

         for (long t=0; t<Size; t++)  OutputArray[t] += RecvBuffer[t];
      }

      else if ( MyRank == SendRank[SendID] )
      {
#        ifdef FLOAT8
         MPI_Send( OutputArray, Size, MPI_DOUBLE, RecvRank, 0, MPI_COMM_WORLD );
#        else
         MPI_Send( OutputArray, Size, MPI_FLOAT,  RecvRank, 0, MPI_COMM_WORLD );
#        endif
      }

      MPI_Barrier( MPI_COMM_WORLD );
   }


   delete [] SendRank;
   delete [] RecvBuffer;


   if ( MyRank == 0 )   cout << "done" << endl;

} // FUNCTION : SumOverRanks
#endif // #ifndef SERIAL



//-------------------------------------------------------------------------------------------------------
// Function    :  LoadTree
// Description :  Load the octree structure and verify it is consistent with the data file
//-------------------------------------------------------------------------------------------------------
void LoadTree()
{

   if ( MyRank == 0 )   Aux_Message( stdout, "%s ...\n", __FUNCTION__ );


// 1. load the tree
   for (int TargetRank=0; TargetRank<NGPU; TargetRank++)
   {
      if ( MyRank == 0 )   Aux_Message( stdout, "   MPI_Rank %3d ... ", TargetRank );

      if ( MyRank == TargetRank )
      {
         int     nlv, data_id, nblock_tot;
         long    step;
         double  time, box_size[3];
         int    *nblock = NULL;
         double *dh     = NULL;

         FILE *FileTree = fopen( FileName_Tree, "rb" );

         if ( FileTree == NULL )    Aux_Error( ERROR_INFO, "the input tree file \"%s\" does not exist !!\n", FileName_Tree );

         fread( &nlv,        sizeof(int),    1,   FileTree );
         fread( &data_id,    sizeof(int),    1,   FileTree );
         fread( &nblock_tot, sizeof(int),    1,   FileTree );
         fread( &step,       sizeof(long),   1,   FileTree );
         fread( &time,       sizeof(double), 1,   FileTree );
         fread(  box_size,   sizeof(double), 3,   FileTree );

         nblock = new int    [nlv];
         dh     = new double [nlv];
         fread(  nblock,     sizeof(int),    nlv, FileTree );
         fread(  dh,         sizeof(double), nlv, FileTree );

         tree = new tree_t( nlv, data_id, time, step, nblock, dh, box_size );

         fread( tree->block[0], sizeof(block_t), tree->nblock_tot, FileTree );

         fclose( FileTree );

         delete [] nblock;
         delete [] dh;
      } // if ( MyRank == TargetRank )

      MPI_Barrier( MPI_COMM_WORLD );

      if ( MyRank == 0 )   Aux_Message( stdout, "done\n" );
   } // for (int TargetRank=0; TargetRank<NGPU; TargetRank++)


// 2. verify the loaded tree is consistent with the data file
   if ( MyRank == 0 )
   {
      FILE *FileGAMER = fopen( FileName_In, "rb" );

      if ( FileGAMER == NULL )   Aux_Error( ERROR_INFO, "the input data file \"%s\" does not exist !!\n", FileName_In );

      const int Offset_Float8 = 256 + 5*sizeof(int) + 9*sizeof(bool);
      const int Offset_NLv    = Offset_Float8 + 1*sizeof(int) + 6*sizeof(bool);
      long    FormatVersion, HeaderSize, CheckCode_Ref, Step, CheckCode;
      int     NLv, DataID;
      bool    Float8;
      double *Time   = NULL;
      int    *NBlock = NULL;

//    load
      fread( &FormatVersion, sizeof(long),     1, FileGAMER );
      fread( &HeaderSize,    sizeof(long),     1, FileGAMER );
      fread( &CheckCode_Ref, sizeof(long),     1, FileGAMER );

      fseek( FileGAMER, Offset_Float8, SEEK_SET );
      fread( &Float8,        sizeof(bool),     1, FileGAMER );

      fseek( FileGAMER, Offset_NLv,    SEEK_SET );
      fread( &NLv,           sizeof(int),      1, FileGAMER );

      fseek( FileGAMER, HeaderSize,    SEEK_SET );
      fread( &CheckCode,     sizeof(long),     1, FileGAMER );

      Time   = new double [NLv];
      NBlock = new int    [NLv];

      fread( &DataID,        sizeof(int),      1, FileGAMER );
      fread( Time,           sizeof(double), NLv, FileGAMER );
      fread( &Step,          sizeof(long),     1, FileGAMER );
      fread( NBlock,         sizeof(int),    NLv, FileGAMER );

//    check
      fprintf( stdout, "   File format = %ld\n", FormatVersion );

      if ( FormatVersion < 1200 )
      {
         fprintf( stderr, "ERROR : unsupported data format version (only support version >= 1200) !!\n" );
         exit( 1 ) ;
      }

      if ( CheckCode != CheckCode_Ref )
      {
         fprintf( stderr, "ERROR : incorrect check code in the data file (input %ld <-> expeect %ld) !!\n",
                  CheckCode, CheckCode_Ref );
         exit( 1 );
      }

#     ifdef FLOAT8
      if ( !Float8 )
      {
         fprintf( stderr, "ERROR : %s : data file (%s) != runtime (%s) !!\n", "FLOAT8", "OFF", "ON" );
         exit( 1 ) ;
      }
#     else
      if (  Float8 )
      {
         fprintf( stderr, "ERROR : %s : data file (%s) != runtime (%s) !!\n", "FLOAT8", "ON", "OFF" );
         exit( 1 ) ;
      }
#     endif

      if ( NLv != tree->nlv )
      {
         fprintf( stderr, "ERROR : inconsistent \"%s\" (data: %d <-> tree: %d) !!\n",
                  "NLv", NLv, tree->nlv );
         exit( 1 );
      }

      if ( DataID != tree->data_id )
      {
         fprintf( stderr, "ERROR : inconsistent \"%s\" (data: %d <-> tree: %d) !!\n",
                  "DataID", DataID, tree->data_id );
         exit( 1 );
      }

      if ( Time[0] != tree->time )
      {
         fprintf( stderr, "ERROR : inconsistent \"%s\" (data: %13.7e <-> tree: %13.7e) !!\n",
                  "Time", Time[0], tree->time );
         exit( 1 );
      }

      if ( Step != tree->step )
      {
         fprintf( stderr, "ERROR : inconsistent \"%s\" (data: %ld <-> tree: %ld) !!\n",
                  "Step", Step, tree->step );
         exit( 1 );
      }

      for (int lv=0; lv<NLv; lv++)
      {
         if ( NBlock[lv] != tree->nblock[lv] )
         {
            fprintf( stderr, "ERROR : inconsistent \"%s\" at level %d (data: %d <-> tree: %d) !!\n",
                     "NBlock", lv, NBlock[lv], tree->nblock[lv] );
            exit( 1 );
         }
      }

      fclose( FileGAMER );

      delete [] Time;
      delete [] NBlock;

      Aux_Message( stdout, "   Check ... passed\n" );
   } // if ( MyRank == 0 )

   MPI_Barrier( MPI_COMM_WORLD );


   if ( MyRank == 0 )   Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : LoadTree



//-------------------------------------------------------------------------------------------------------
// Function    :  main
// Description :
//-------------------------------------------------------------------------------------------------------
int main( int argc, char ** argv )
{

   Init_MPI( &argc, &argv );

   ReadOption( argc, argv );

   if ( UseTree )    LoadTree();

   LoadData();

   if ( MyRank == 0 )   TakeNote( argc, argv );

   AllocateOutputArray();

   Refine2TargetLevel();

#  ifndef SERIAL
   if ( OutputXYZ == 4  ||  OutputXYZ == 5  ||  OutputXYZ == 6 )     SumOverRanks();
#  endif

   End_MemFree();

   Output();

   delete [] OutputArray;

   MPI_Finalize();

   if ( MyRank == 0 )   cout << "Program terminated successfully" << endl << endl;

   return 0;

} // FUNCTION : main


