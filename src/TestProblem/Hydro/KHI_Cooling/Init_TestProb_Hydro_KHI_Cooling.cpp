#include "GAMER.h"

static double RandomNumber(RandomNumber_t *RNG, const double Min, const double Max );



// problem-specific global variables
// =======================================================================================
static int    KHC_RSeed;             // random seed
static double KHC_RAmp;              // random number amplitude in both vx, vy, and vz
static double KHC_Vx1;               // background velocity x in the upper region [cm/s]
static double KHC_Vx2;               // background velocity x in the lower region [cm/s]
static double KHC_Vz1;               // background velocity z in the upper region [cm/s]
static double KHC_Vz2;               // background velocity z in the lower region [cm/s]
static double KHC_Temp1;             // background temperature in the upper region [K]
static double KHC_Temp2;             // background temperature in the lower region [K]
static double KHC_n1;                // background total number density in the upper region [cm^-3]
static double KHC_n2;                // background total number density in the lower region [cm^-3]
       int    KHC_RefineShearMaxLv;  // refine the regions around the shear plane to this level
       int    KHC_PeriodicYFactor;   // decompose the simulation domain further into N periodic subdomains along the y axis

static RandomNumber_t *RNG = NULL;
// =======================================================================================

// problem-specific function prototypes
bool Flag_KHI_Cooling( const int i, const int j, const int k, const int lv, const int PID, const double *Threshold );




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

#  ifdef GRAVITY
   Aux_Error( ERROR_INFO, "GRAVITY must be disabled !!\n" );
#  endif

#  ifdef COMOVING
   Aux_Error( ERROR_INFO, "COMOVING must be disabled !!\n" );
#  endif

#  ifdef PARTICLE
   Aux_Error( ERROR_INFO, "PARTICLE must be disabled !!\n" );
#  endif


   if ( OPT__INIT_GRID_WITH_OMP )
      Aux_Error( ERROR_INFO, "please disable \"OPT__INIT_GRID_WITH_OMP\" !!\n" );


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Validating test problem %d ... done\n", TESTPROB_ID );

//   if ( KH_PeriodicYFactor != 1 )
//   Aux_Error( ERROR_INFO,
//              "KH_PeriodicYFactor must be 1 for this setup !!\n" );

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
// **************************************************************************************************************************
// LOAD_PARA( load_mode, "KEY_IN_THE_FILE",     &VARIABLE,               DEFAULT,      MIN,              MAX               );
// **************************************************************************************************************************
   LOAD_PARA( load_mode, "KHC_RSeed",            &KHC_RSeed,              123,             0,                NoMax_int         );
   LOAD_PARA( load_mode, "KHC_RAmp",             &KHC_RAmp,               2.0e4,           Eps_double,       NoMax_double      );
   LOAD_PARA( load_mode, "KHC_Vx1",              &KHC_Vx1,               -1.0e6,           NoMin_double,     NoMax_double      );
   LOAD_PARA( load_mode, "KHC_Vx2",              &KHC_Vx2,                1.0e6,           NoMin_double,     NoMax_double      );
   LOAD_PARA( load_mode, "KHC_Vz1",              &KHC_Vz1,                0.0,             NoMin_double,     NoMax_double      );
   LOAD_PARA( load_mode, "KHC_Vz2",              &KHC_Vz2,                0.0,             NoMin_double,     NoMax_double      );
   LOAD_PARA( load_mode, "KHC_Temp1",            &KHC_Temp1,              3.74e4,          Eps_double,       NoMax_double      );
   LOAD_PARA( load_mode, "KHC_Temp2",            &KHC_Temp2,              7.48e4,          Eps_double,       NoMax_double      );
   LOAD_PARA( load_mode, "KHC_n1",               &KHC_n1,                 2.0e-2,          Eps_double,       NoMax_double      );
   LOAD_PARA( load_mode, "KHC_n2",               &KHC_n2,                 1.0e-2,          Eps_double,       NoMax_double      );
   LOAD_PARA( load_mode, "KHC_RefineShearMaxLv", &KHC_RefineShearMaxLv,   2,               0,                TOP_LEVEL         );
   LOAD_PARA( load_mode, "KHC_PeriodicYFactor",  &KHC_PeriodicYFactor,    1,               1,                NoMax_int         );

} // FUNCTION : LoadInputTestProb



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

// get the number of OpenMP threads
   int NT;
#  ifdef OPENMP
#  pragma omp parallel
#  pragma omp master
   {  NT = omp_get_num_threads();  }
#  else
   {  NT = 1;                      }
#  endif

// allocate RNG
   RNG = new RandomNumber_t( NT );

// set the random seed of each MPI rank
   for (int t=0; t<NT; t++) {
      RNG->SetSeed(t, KHC_RSeed+MPI_Rank*1000+t);
   }

// (2) set the problem-specific derived parameters


// (3) reset other general-purpose parameters
//     --> a helper macro PRINT_RESET_PARA is defined in Macro.h
   const double End_T_Default    = 3.0e5;
   const long   End_Step_Default = __INT_MAX__;

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
      Aux_Message( stdout, "  test problem ID     = %d\n",      TESTPROB_ID );
      Aux_Message( stdout, "  KHC_RSeed            = %d\n",      KHC_RSeed );
      Aux_Message( stdout, "  KHC_RAmp             = % 14.7e\n", KHC_RAmp );
      Aux_Message( stdout, "  KHC_Vx1              = % 14.7e\n", KHC_Vx1 );
      Aux_Message( stdout, "  KHC_Vx2              = % 14.7e\n", KHC_Vx2 );
      Aux_Message( stdout, "  KHC_Vz1              = % 14.7e\n", KHC_Vz1 );
      Aux_Message( stdout, "  KHC_Vz2              = % 14.7e\n", KHC_Vz2 );
      Aux_Message( stdout, "  KHC_Temp1            = % 14.7e\n", KHC_Temp1 );
      Aux_Message( stdout, "  KHC_Temp2            = % 14.7e\n", KHC_Temp2 );
      Aux_Message( stdout, "  KHC_n1               = % 14.7e\n", KHC_n1 );
      Aux_Message( stdout, "  KHC_n2               = % 14.7e\n", KHC_n2 );
      Aux_Message( stdout, "  KHC_RefineShearMaxLv = %d\n",      KHC_RefineShearMaxLv );
      Aux_Message( stdout, "  KHC_PeriodicYFactor  = %d\n",      KHC_PeriodicYFactor );
      Aux_Message( stdout, "=============================================================================\n" );
   }


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Setting runtime parameters ... done\n" );

} // FUNCTION : SetParameter




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

   const double dy_periodic = amr->BoxSize[1] / KHC_PeriodicYFactor;
   const double y_shear     = 0.5*dy_periodic;
   const double y_periodic  = fmod( y, dy_periodic );   

   double Dens, MomX, MomY, MomZ, Pres, Eint, Etot;

// gas settings
   const double X_H = 0.7;                                                        // hydrogen mass fraction
   const double Z   = 0.018;                                                      // metal mass fraction
   const double mu  = 1.0 / ( 2.0 * X_H + 0.75 * ( 1.0 - X_H - Z ) + 0.5 * Z );   // mean molecular weight

// region 1 : cold gas
   if ( y_periodic >= y_shear )
   {
      Dens = ( KHC_n1 * mu * MU_NORM ) / UNIT_D;
      MomX = Dens*( KHC_Vx1 + RandomNumber(RNG,-KHC_RAmp,KHC_RAmp) ) / ( UNIT_L / UNIT_T );
      MomY = Dens*(     0.0 + RandomNumber(RNG,-KHC_RAmp,KHC_RAmp) ) / ( UNIT_L / UNIT_T );
      MomZ = Dens*( KHC_Vz1 + RandomNumber(RNG,-KHC_RAmp,KHC_RAmp) ) / ( UNIT_L / UNIT_T );
      Pres = EoS_DensTemp2Pres_CPUPtr( Dens, KHC_Temp1, NULL, EoS_AuxArray_Flt, EoS_AuxArray_Int, h_EoS_Table );
   }

// region 2 : hot gas
   else
   {
      Dens = ( KHC_n2 * mu * MU_NORM ) / UNIT_D;
      MomX = Dens*( KHC_Vx2 + RandomNumber(RNG,-KHC_RAmp,KHC_RAmp) ) / ( UNIT_L / UNIT_T );
      MomY = Dens*(     0.0 + RandomNumber(RNG,-KHC_RAmp,KHC_RAmp) ) / ( UNIT_L / UNIT_T );
      MomZ = Dens*( KHC_Vz2 + RandomNumber(RNG,-KHC_RAmp,KHC_RAmp) ) / ( UNIT_L / UNIT_T );
      Pres = EoS_DensTemp2Pres_CPUPtr( Dens, KHC_Temp2, NULL, EoS_AuxArray_Flt, EoS_AuxArray_Int, h_EoS_Table );
   }

// compute the total gas energy
   Eint = EoS_DensPres2Eint_CPUPtr( Dens, Pres, NULL, EoS_AuxArray_Flt,
                                    EoS_AuxArray_Int, h_EoS_Table );   // assuming EoS requires no passive scalars
   Etot = Hydro_ConEint2Etot( Dens, MomX, MomY, MomZ, Eint, 0.0 );     // do NOT include magnetic energy here

// set the output array
   fluid[DENS] = Dens;
   fluid[MOMX] = MomX;
   fluid[MOMY] = MomY;
   fluid[MOMZ] = MomZ;
   fluid[ENGY] = Etot;

} // FUNCTION : SetGridIC
#endif // #if ( MODEL == HYDRO )



//-------------------------------------------------------------------------------------------------------
// Function    :  RandomNumber
// Description :  Generate a single random number
//
// Parameter   :  RNG : Random number generator
//                Min : Lower limit of the random number
//                Max : Upper limit of the random number
//
// Return      :  Random number
//-------------------------------------------------------------------------------------------------------
double RandomNumber(RandomNumber_t *RNG, const double Min, const double Max )
{

// thread-private variables
#  ifdef OPENMP
   const int TID = omp_get_thread_num();
#  else
   const int TID = 0;
#  endif

   return RNG->GetValue( TID, Min, Max );

} // FUNCTION : RandomNumber




//-------------------------------------------------------------------------------------------------------
// Function    :  Init_TestProb_Hydro_KHI_Cooling
// Description :  Test problem initializer
//
// Note        :  None
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Init_TestProb_Hydro_KHI_Cooling()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );


// validate the compilation flags and runtime parameters
   Validate();


#  if ( MODEL == HYDRO )
// set the problem-specific runtime parameters
   SetParameter();


// set the function pointers of various problem-specific routines
   Init_Function_User_Ptr    = SetGridIC;
   Flag_User_Ptr             = Flag_KHI_Cooling;
#  ifdef SUPPORT_HDF5
   Output_HDF5_InputTest_Ptr = LoadInputTestProb;
#  endif
#  endif // #if ( MODEL == HYDRO )


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : Init_TestProb_Hydro_KHI_Cooling

