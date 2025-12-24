#include "GAMER.h"
#include "TestProb.h"



// problem-specific global variables
// =======================================================================================
static int        *TargetCols = new int [6];      // Index of columns read from the turbulence table 
static int        ColIdx_VelX;                    // Column index of x direction velocity in the turbulence table 
static int        ColIdx_VelY;                    // Column index of y direction velocity in the turbulence table 
static int        ColIdx_VelZ;                    // Column index of z direction velocity in the turbulence table 
static int        tur_table_NBin;                 // number of row in turbulence table obtained by Aux_LoadTable
static int        size;                           // turbulence cell number                  

static double     *Table_Rescaled_VelX;           // Table recording the rescaled velocity in x direction
static double     *Table_Rescaled_VelY;           // Table recording the rescaled velocity in y direction
static double     *Table_Rescaled_VelZ;           // Table recording the rescaled velocity in z direction

static double     SinkParTest_R0;                 // The size of the cloud
static double     SinkParTest_Omega0;             // The angular velocity of the cloud
static double     SinkParTest_Core_Mass;          // The mass of the cloud
static double     SinkParTest_Delta_Dens;         // The density of the cloud
static double     SinkParTest_Dens_Contrast;      // The density contrast between the cloud and background
static double     SinkParTest_B0;                 // The magnetic field
static double     SinkParTest_theta_B;            // The angle between the magnetic field and the z-axis
static double     SinkParTest_Mach_num;           // The Mach number of the turbulence
static char       Tur_Table[MAX_STRING];          // The turbulence table
static double     SinkParTest_Cs;                 // sound spped
static double     SinkParTest_Rho0;               // The density of the cloud

double            SinkParTest_rho_AD;             // The adiabatic density threshold



// =======================================================================================

#  if ( EOS == EOS_USER )
void EoS_Init_Barotropic_SinkParTest();
#  endif


//-------------------------------------------------------------------------------------------------------
// Function    :  Validate
// Description :  Validate the compilation flags and runtime parameters for this test problem
//
// Note        :  None
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Validate()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Validating test problem %d ...\n", TESTPROB_ID );


#  if ( MODEL != HYDRO )
   Aux_Error( ERROR_INFO, "MODEL != HYDRO !!\n" );
#  endif

#  ifndef GRAVITY
   Aux_Error( ERROR_INFO, "GRAVITY must be enabled !!\n" );
#  endif

#  ifndef STAR_FORMATION
   Aux_Error( ERROR_INFO, "STAR_FORMATION must be enabled !!\n" );
#  endif

#  ifndef FEEDBACK
   Aux_Error( ERROR_INFO, "FEEDBACK must be enabled !!\n" );
#  endif

#  ifndef PARTICLE
   Aux_Error( ERROR_INFO, "PARTICLE must be enabled !!\n" );
#  endif

#  ifndef MHD
   Aux_Error( ERROR_INFO, "MHD must be enabled !!\n" );
#  endif

#  if ( EOS != EOS_USER )
   Aux_Error( ERROR_INFO, "EOS != EOS_USER !!\n" );
#  endif

#  ifdef PARTICLE
   if ( OPT__INIT == INIT_BY_FUNCTION  &&  amr->Par->Init != PAR_INIT_BY_FUNCTION )
      Aux_Error( ERROR_INFO, "please set PAR_INIT = 1 (by FUNCTION) !!\n" );
#  endif

#  ifdef GRAVITY
   if ( !OPT__SELF_GRAVITY )
   Aux_Error( ERROR_INFO, "must enable OPT__SELF_GRAVITY !!\n" );

   if ( OPT__EXT_ACC )
   Aux_Error( ERROR_INFO, "must disable OPT__EXT_ACC !!\n" );

   if ( OPT__EXT_POT )
   Aux_Error( ERROR_INFO, "must disable OPT__EXT_POT !!\n" );
#  endif

   if ( MPI_Rank == 0 )
   {
#     ifndef DUAL_ENERGY
         Aux_Message( stderr, "WARNING : it's recommended to enable DUAL_ENERGY for this test\n" );
#     endif

      for (int f=0; f<6; f++)
      if ( amr->BoxSize[0] != amr->BoxSize[1]  ||  amr->BoxSize[0] != amr->BoxSize[2] )
         Aux_Error( ERROR_INFO, "simulation domain is not cubic !!\n" );
   }


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Validating test problem %d ... done\n", TESTPROB_ID );

} // FUNCTION : Validate



#if ( MODEL == HYDRO )
//-------------------------------------------------------------------------------------------------------
// Function    :  LoadInputTestProb
// Description :  Read problem-specific runtime parameters from Input__TestProb and store them in HDF5 snapshots (Data_*)
//
// Note        :  1. Invoked by SetParameter() to read parameters
//                2. Invoked by Output_DumpData_Total_HDF5() using the function pointer Output_HDF5_InputTest_Ptr to store parameters
//                3. If there is no problem-specific runtime parameter to load, add at least one parameter
//                   to prevent an empty structure in HDF5_Output_t
//                   --> Example:
//                       LOAD_PARA( load_mode, "TestProb_ID", &TESTPROB_ID, TESTPROB_ID, TESTPROB_ID, TESTPROB_ID );
//
// Parameter   :  load_mode      : Mode for loading parameters
//                                 --> LOAD_READPARA    : Read parameters from Input__TestProb
//                                     LOAD_HDF5_OUTPUT : Store parameters in HDF5 snapshots
//                ReadPara       : Data structure for reading parameters (used with LOAD_READPARA)
//                HDF5_InputTest : Data structure for storing parameters in HDF5 snapshots (used with LOAD_HDF5_OUTPUT)
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void LoadInputTestProb( const LoadParaMode_t load_mode, ReadPara_t *ReadPara, HDF5_Output_t *HDF5_InputTest )
{

#  ifndef SUPPORT_HDF5
   if ( load_mode == LOAD_HDF5_OUTPUT )   Aux_Error( ERROR_INFO, "please turn on SUPPORT_HDF5 in the Makefile for load_mode == LOAD_HDF5_OUTPUT !!\n" );
#  endif

   if ( load_mode == LOAD_READPARA     &&  ReadPara       == NULL )   Aux_Error( ERROR_INFO, "load_mode == LOAD_READPARA and ReadPara == NULL !!\n" );
   if ( load_mode == LOAD_HDF5_OUTPUT  &&  HDF5_InputTest == NULL )   Aux_Error( ERROR_INFO, "load_mode == LOAD_HDF5_OUTPUT and HDF5_InputTest == NULL !!\n" );

// add parameters in the following format:
// --> note that VARIABLE, DEFAULT, MIN, and MAX must have the same data type
// --> some handy constants (e.g., NoMin_int, Eps_float, ...) are defined in "include/ReadPara.h"
// --> LOAD_PARA() is defined in "include/TestProb.h"
// ***************************************************************************************************************************
// LOAD_PARA( load_mode, "KEY_IN_THE_FILE",      &VARIABLE,                  DEFAULT,       MIN,              MAX               );
// ***************************************************************************************************************************
   LOAD_PARA( load_mode, "R0",                   &SinkParTest_R0,            3.0e16,        Eps_double,       NoMax_double      );
   LOAD_PARA( load_mode, "Omega0",               &SinkParTest_Omega0,        0.0,           NoMin_double,     NoMax_double      );
   LOAD_PARA( load_mode, "Core_Mass",            &SinkParTest_Core_Mass,     1.0,           0.0,              NoMax_double      );
   LOAD_PARA( load_mode, "Delta_Dens",           &SinkParTest_Delta_Dens,    0.1,           0.0,              NoMax_double      );
   LOAD_PARA( load_mode, "Dens_Contrast",        &SinkParTest_Dens_Contrast, 1.0e2,         0.0,              NoMax_double      );
   LOAD_PARA( load_mode, "B0"     ,              &SinkParTest_B0,            0.0,           NoMin_double,     NoMax_double      );
   LOAD_PARA( load_mode, "theta_B",              &SinkParTest_theta_B,       0.0,           NoMin_double,     NoMax_double      );
   LOAD_PARA( load_mode, "Mach_num",             &SinkParTest_Mach_num,      0.0,           0.0,              NoMax_double      );
   LOAD_PARA( load_mode, "rho_AD",               &SinkParTest_rho_AD,        1e-14,         0.0,              NoMax_double      );
   LOAD_PARA( load_mode, "Tur_Table",            Tur_Table,                  NoDef_str,     Useless_str,      Useless_str       );

} // FUNCITON : LoadInputTestProb

//-------------------------------------------------------------------------------------------------------
// Function    :  SetParameter
// Description :  Load and set the problem-specific runtime parameters
//
// Note        :  1. Filename is set to "Input__TestProb" by default
//                2. Major tasks in this function:
//                   (1) load the problem-specific runtime parameters
//                   (2) set the problem-specific derived parameters
//                   (3) reset other general-purpose parameters if necessary
//                   (4) make a note of the problem-specific parameters
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void SetParameter()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Setting runtime parameters ...\n" );

// (1) load the problem-specific runtime parameters
// (1-1) read parameters from Input__TestProb
   const char FileName[] = "Input__TestProb";
   ReadPara_t *ReadPara  = new ReadPara_t;

   LoadInputTestProb( LOAD_READPARA, ReadPara, NULL );

   ReadPara->Read( FileName );

   delete ReadPara;

// (1-2) set the default values
   TargetCols[0] =  0;
   TargetCols[1] =  1;
   TargetCols[2] =  2;
   TargetCols[3] =  3;
   TargetCols[4] =  4;
   TargetCols[5] =  5;
   ColIdx_VelX   =  3;
   ColIdx_VelY   =  4;
   ColIdx_VelZ   =  5;
   size = 129;


// (2) set the problem-specific derived parameters
   SinkParTest_Core_Mass *= Const_Msun/UNIT_M;
   SinkParTest_Cs         = SQRT( ( Const_kB*ISO_TEMP/UNIT_E ) / ( MOLECULAR_WEIGHT*Const_amu/UNIT_M ));
   SinkParTest_R0         /= UNIT_L;
   SinkParTest_Rho0       = 3.0 * SinkParTest_Core_Mass / (4.0 * M_PI * CUBE(SinkParTest_R0));
   SinkParTest_Omega0     /= 1/UNIT_T;
   SinkParTest_rho_AD     /= UNIT_D;
   SinkParTest_theta_B    = SinkParTest_theta_B*M_PI/180; // degree to radian
   SinkParTest_B0         /= UNIT_B;

// (3) reset other general-purpose parameters
//     --> a helper macro PRINT_WARNING is defined in TestProb.h
   const long   End_Step_Default = __INT_MAX__;
   const double End_T_Default    = 1.3065165e+12;

   if ( END_STEP < 0 ) {
      END_STEP = End_Step_Default;
      PRINT_RESET_PARA( END_STEP, FORMAT_LONG, "" );
   }

   if ( END_T < 0.0 ) {
      END_T = End_T_Default;
      PRINT_RESET_PARA( END_T, FORMAT_REAL, "" );
   }


// (4) make a note
   if ( MPI_Rank == 0 )
   {
      Aux_Message( stdout, "=============================================================================\n" );
      Aux_Message( stdout, "  test problem ID       = %d\n",             TESTPROB_ID                         );
      Aux_Message( stdout, "  Sound speed           = %13.7e km/s\n",   SinkParTest_Cs*UNIT_V/Const_km       );
      Aux_Message( stdout, "  R0                    = %13.7e cm\n",     SinkParTest_R0*UNIT_L                );
      Aux_Message( stdout, "  Rho0                  = %13.7e g/cm3\n",  SinkParTest_Rho0*UNIT_D              );
      Aux_Message( stdout, "  Omega0                = %13.7e /s\n",     SinkParTest_Omega0*UNIT_T            );
      Aux_Message( stdout, "  Core Mass             = %13.7e M_sun\n",  SinkParTest_Core_Mass/Const_Msun     );
      Aux_Message( stdout, "  B field               = %13.7e G\n",      SinkParTest_B0                       );
      Aux_Message( stdout, "  B angle               = %13.7e radian\n", SinkParTest_theta_B                  );
      Aux_Message( stdout, "  Mach number           = %13.7e \n",       SinkParTest_Mach_num                 );
      Aux_Message( stdout, "  Turbulence table      = %s\n",            Tur_Table                            );
      Aux_Message( stdout, "  rho_AD                = %13.7e g/cm3\n",  SinkParTest_rho_AD                   );
      Aux_Message( stdout, "=============================================================================\n" );
   }


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Setting runtime parameters ... done\n" );

} // FUNCTION : SetParameter

//-------------------------------------------------------------------------------------------------------
// Function    :  Load_Turbulence_SinkParTest
// Description : load turbulence and calculate Vrms_Scale
// Note        :
// Parameter   :
// Return      :
//-------------------------------------------------------------------------------------------------------
void Load_Turbulence_SinkParTest()
{
   const bool RowMajor_No  = false;           // load data into the column major
   const bool AllocMem_Yes = true;            // allocate memory

   double *Table_VelX, *Table_VelY, *Table_VelZ;     // used to store the readed data
   double *tur_table = NULL;                         // used to store turbulence (1D)

   double Total_VelX = 0.0; // used to calculate Vrms
   double Total_VelY = 0.0;
   double Total_VelZ = 0.0;
   double Total_VelX_SQR = 0.0;
   double Total_VelY_SQR = 0.0;
   double Total_VelZ_SQR = 0.0;
   double Vrms, Vrms_Scale;                     // used to rescale velocity

   tur_table_NBin = Aux_LoadTable( tur_table, Tur_Table, 6, TargetCols, RowMajor_No, AllocMem_Yes );

   Table_VelX  = tur_table + ColIdx_VelX * tur_table_NBin;
   Table_VelY  = tur_table + ColIdx_VelY * tur_table_NBin;
   Table_VelZ  = tur_table + ColIdx_VelZ * tur_table_NBin;

   Table_Rescaled_VelX = tur_table + ColIdx_VelX * tur_table_NBin;
   Table_Rescaled_VelY = tur_table + ColIdx_VelY * tur_table_NBin;
   Table_Rescaled_VelZ = tur_table + ColIdx_VelZ * tur_table_NBin;

   for ( int i = 0; i < tur_table_NBin; i++ )
   {
      Total_VelX += Table_VelX[i];
      Total_VelY += Table_VelY[i];
      Total_VelZ += Table_VelZ[i];

      Total_VelX_SQR += SQR(Table_VelX[i]);
      Total_VelY_SQR += SQR(Table_VelY[i]);
      Total_VelZ_SQR += SQR(Table_VelZ[i]);
   }

   // Vrms = SQRT( ( Vx^2 + Vy^2 + Vz^2 ) / N + ( Vx + Vy + Vz / N) ^ 2 )
   Vrms = SQRT( (Total_VelX_SQR + Total_VelY_SQR + Total_VelZ_SQR) / tur_table_NBin - 
                SQR( (Total_VelX + Total_VelY + Total_VelZ) / tur_table_NBin ) );
   Vrms_Scale = SinkParTest_Mach_num * SinkParTest_Cs / Vrms;

   // Rescale velocity
   for ( int i = 0; i < tur_table_NBin; i++ )
   {
      Table_Rescaled_VelX[i] = Vrms_Scale * ( Table_VelX[i] - Total_VelX / tur_table_NBin );
      Table_Rescaled_VelY[i] = Vrms_Scale * ( Table_VelY[i] - Total_VelY / tur_table_NBin );
      Table_Rescaled_VelZ[i] = Vrms_Scale * ( Table_VelZ[i] - Total_VelZ / tur_table_NBin );
   }

   delete [] tur_table;
} // Function : Load_Turbulence_SinkParTest

//-------------------------------------------------------------------------------------------------------
// Function    :  SetGridIC
// Description :  Set the problem-specific initial condition on grids
//
// Note        :  1. This function may also be used to estimate the numerical errors when OPT__OUTPUT_USER is enabled
//                   --> In this case, it should provide the analytical solution at the given "Time"
//                2. This function will be invoked by multiple OpenMP threads when OPENMP is enabled
//                   --> Please ensure that everything here is thread-safe
//                3. Even when DUAL_ENERGY is adopted for HYDRO, one does NOT need to set the dual-energy variable here
//                   --> It will be calculated automatically
//
// Parameter   :  fluid    : Fluid field to be initialized
//                x/y/z    : Physical coordinates
//                Time     : Physical time
//                lv       : Target refinement level
//                AuxArray : Auxiliary array
//
// Return      :  fluid
//-------------------------------------------------------------------------------------------------------
void SetGridIC( real fluid[], const double x, const double y, const double z, const double Time,
                const int lv, double AuxArray[] )
{
   const double dx = x - amr->BoxCenter[0];
   const double dy = y - amr->BoxCenter[1];
   const double dz = z - amr->BoxCenter[2];
   const double Rs = SQRT( SQR(dx) + SQR(dy) + SQR(dz) );

   real Dens, MomX, MomY, MomZ, Pres, Eint, Etot;
   real VelX, VelY, VelZ;

   // Load the turbulence field if not restarting
   if ( OPT__INIT != INIT_BY_RESTART ) Load_Turbulence_SinkParTest();

   int i = (int) ( ( x / amr->BoxSize[0] ) * size );    // turbulence box index (cude)
   int j = (int) ( ( y / amr->BoxSize[1] ) * size );
   int k = (int) ( ( z / amr->BoxSize[2] ) * size );
   int index = i * SQR(size) + j * size + k;
   if ( i < 0 || i > size   ) Aux_Error( ERROR_INFO, "index is out of bound,  i = %d", i  );
   if ( j < 0 || j > size   ) Aux_Error( ERROR_INFO, "index is out of bound,  j = %d", j  );
   if ( k < 0 || k > size   ) Aux_Error( ERROR_INFO, "index is out of bound,  k = %d", k );
   if ( index < 0 || index > tur_table_NBin ) Aux_Error( ERROR_INFO, "index is out of bound, index = %d", index );

   VelX = Table_Rescaled_VelX[ index ];
   VelY = Table_Rescaled_VelY[ index ];
   VelZ = Table_Rescaled_VelZ[ index ];

   if ( Rs < SinkParTest_R0 )
   {
      Dens = SinkParTest_Rho0 * (1 + SinkParTest_Delta_Dens * COS(2 * ATAN(dy/dx)));
      VelX -= SinkParTest_Omega0 * dy;
      VelY += SinkParTest_Omega0 * dx;
   }
   else
   {
      Dens = SinkParTest_Rho0 / SinkParTest_Dens_Contrast;
   }

   Pres = EoS_DensTemp2Pres_CPUPtr( Dens, ISO_TEMP, NULL, EoS_AuxArray_Flt,
                                    EoS_AuxArray_Int, h_EoS_Table );
   Eint = EoS_DensPres2Eint_CPUPtr( Dens, Pres, NULL, EoS_AuxArray_Flt,
                                    EoS_AuxArray_Int, h_EoS_Table );

   MomX = Dens * VelX;
   MomY = Dens * VelY;
   MomZ = Dens * VelZ;
   Etot = Hydro_ConEint2Etot( Dens, MomX, MomY, MomZ, Eint, 0.0 );     // do NOT include magnetic energy here

   fluid[DENS] = Dens;
   fluid[MOMX] = MomX;
   fluid[MOMY] = MomY;
   fluid[MOMZ] = MomZ;
   fluid[ENGY] = Etot;

} // FUNCTION : SetGridIC

#ifdef MHD
//-------------------------------------------------------------------------------------------------------
// Function    :  SetBFieldIC
// Description :  Set the problem-specific initial condition of magnetic field
//
// Note        :  1. This function will be invoked by multiple OpenMP threads when OPENMP is enabled
//                   (unless OPT__INIT_GRID_WITH_OMP is disabled)
//                   --> Please ensure that everything here is thread-safe
//
// Parameter   :  magnetic : Array to store the output magnetic field
//                x/y/z    : Target physical coordinates
//                Time     : Target physical time
//                lv       : Target refinement level
//                AuxArray : Auxiliary array
//
// Return      :  magnetic
//-------------------------------------------------------------------------------------------------------
void SetBFieldIC( real magnetic[], const double x, const double y, const double z, const double Time,
                  const int lv, double AuxArray[] )
{
   double MagX = 0.0, MagY = 0.0, MagZ = 0.0;
   MagZ = SinkParTest_B0*COS(SinkParTest_theta_B);
   MagY = SinkParTest_B0*SIN(SinkParTest_theta_B);

   magnetic[MAGX] = MagX;
   magnetic[MAGY] = MagY;
   magnetic[MAGZ] = MagZ;

} // FUNCTION : SetBFieldIC
#endif // #ifdef MHD
#endif // #if ( MODEL == HYDRO )

#  ifdef PARTICLE
//-------------------------------------------------------------------------------------------------------
// Function    :  AddNewParticleAttribute_SinkParTest
// Description :  Add the problem-specific particle attributes
//
// Note        :  1. Ref: https://github.com/gamer-project/gamer/wiki/Adding-New-Simulations#v-add-problem-specific-grid-fields-and-particle-attributes
//                2. Invoke AddParticleField() for each of the problem-specific particle attribute:
//                   --> Attribute label sent to AddParticleField() will be used as the output name of the attribute
//                   --> Attribute index returned by AddParticleField() can be used to access the particle attribute data
//                3. Pre-declared attribute indices are put in Field.h
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void AddNewParticleAttribute_SinkParTest()
{

// "Idx_ParMetalFrac" has been predefined in Field.h
   if ( Idx_ParID == Idx_Undefined )
      Idx_ParID = AddParticleAttributeInt( "PAR_ID" );

} // FUNCTION : AddNewParticleAttribute_SinkParTest


void Par_Init_ByFunction_SinkParTest( const long NPar_ThisRank, const long NPar_AllRank,
                                   real_par *ParMass, real_par *ParPosX, real_par *ParPosY, real_par *ParPosZ,
                                   real_par *ParVelX, real_par *ParVelY, real_par *ParVelZ, real_par *ParTime,
                                   long_par *ParType, real_par *AllAttributeFlt[PAR_NATT_FLT_TOTAL],
                                   long_par *AllAttributeInt[PAR_NATT_INT_TOTAL] )
{
   // we only use star particles, so keep here empty.
}  // FUNCTION : Par_Init_ByFunction
#  endif



//-------------------------------------------------------------------------------------------------------
// Function    :  Init_TestProb_Hydro_SinkParTest
// Description :  Test problem initializer
//
// Note        :  None
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Init_TestProb_Hydro_SinkParTest()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );


// validate the compilation flags and runtime parameters
   Validate();


#  if ( MODEL == HYDRO )
// set the problem-specific runtime parameters
   SetParameter();

// set the function pointers of various problem-specific routines
   Init_Function_User_Ptr            = SetGridIC;
#  ifdef MHD
   Init_Function_BField_User_Ptr     = SetBFieldIC;
#  endif
#  endif // #if ( MODEL == HYDRO )

#  ifdef PARTICLE
   Par_Init_ByFunction_Ptr           = Par_Init_ByFunction_SinkParTest; // option: PAR_INIT=1;              example: Particle/Par_Init_ByFunction.cpp
   Par_Init_Attribute_User_Ptr       = AddNewParticleAttribute_SinkParTest; // set PAR_NATT_USER;               example: TestProblem/Hydro/AGORA_IsolatedGalaxy/Init_TestProb_Hydro_AGORA_IsolatedGalaxy.cpp --> AddNewParticleAttribute()
#  endif
#  if ( EOS == EOS_USER )
   EoS_Init_Ptr                      = EoS_Init_Barotropic_SinkParTest; // option: EOS in the Makefile;     example: EoS/User_Template/CPU_EoS_User_Template.cpp
   EoS_End_Ptr                       = NULL;
#  endif
#  ifdef SUPPORT_HDF5
   Output_HDF5_InputTest_Ptr         = LoadInputTestProb;
#  endif


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : Init_TestProb_Hydro_SinkParTest