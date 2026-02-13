#include "GAMER.h"



// problem-specific global variables
// =======================================================================================
static double VorPairLin_BgAmp;     // psi_vorpair(x,y) = BgAmp + WaveAmp*cos(ky*y)*exp( i*(kx*x-Omega*t+Phase0) )
static double VorPairLin_WaveAmp;
static double VorPairLin_Phase0;

static double VorPairLin_kx;
static double VorPairLin_ky;
static double VorPairLin_Omega;

// optional:
static double VorPairLin_ZWaveAmp;  // psi(x, y, z) = psi_vorpair(x,y) + ZWaveAmp * exp( i*(kz*z-ZWaveOmega*t) )
static double VorPairLin_kz;
static double VorPairLin_ZWaveOmega;
// =======================================================================================




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


// errors
#  if ( MODEL != ELBDM )
   Aux_Error( ERROR_INFO, "MODEL != ELBDM !!\n" );
#  endif

#  ifdef GRAVITY
   Aux_Error( ERROR_INFO, "GRAVITY must be disabled !!\n" );
#  endif

#  ifdef PARTICLE
   Aux_Error( ERROR_INFO, "PARTICLE must be disabled !!\n" );
#  endif

// -> when switching from wave scheme back to fluid scheme
//    real and imaginary part need to be converted back to phase
// -> the two vortices are connected with a two pi phase jump
//    while the wave scheme does not see this jump, the fluid scheme does
// -> the wave physics is invariant under deformations of this jump contour.
//    as a result, reconstruction of the phase with and without ELBDM_MATCH_PHASE
//    leads to different jump line contours
// -> when ELBDM_MATCH_PHASE is enabled the jump line contour leads to overrefinement
//    which is why we should disable this option for the vortex pair tests
#  if ( ELBDM_SCHEME == ELBDM_HYBRID )
   if ( ELBDM_MATCH_PHASE )
      Aux_Message( stderr, "WARNING: ELBDM_MATCH_PHASE should be disabled in vortex pair tests !!\n" );
#  endif


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Validating test problem %d ... done\n", TESTPROB_ID );

} // FUNCTION : Validate



#if ( MODEL == ELBDM )
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
// LOAD_PARA( load_mode, "KEY_IN_THE_FILE",      &VARIABLE,              DEFAULT,       MIN,              MAX               );
// ***************************************************************************************************************************
   LOAD_PARA( load_mode, "VorPairLin_BgAmp",     &VorPairLin_BgAmp,      1.0,           0.0,              NoMax_double      );
   LOAD_PARA( load_mode, "VorPairLin_WaveAmp",   &VorPairLin_WaveAmp,   -1.0,           Eps_double,       NoMax_double      );
   LOAD_PARA( load_mode, "VorPairLin_ZWaveAmp",  &VorPairLin_ZWaveAmp,   0.0,           0.0,              NoMax_double      );
   LOAD_PARA( load_mode, "VorPairLin_Phase0",    &VorPairLin_Phase0,     0.0,           NoMin_double,     NoMax_double      );
   LOAD_PARA( load_mode, "VorPairLin_kx",        &VorPairLin_kx,         1.0,           1.0,              NoMax_double      );
   LOAD_PARA( load_mode, "VorPairLin_ky",        &VorPairLin_ky,         1.0,           1.0,              NoMax_double      );
   LOAD_PARA( load_mode, "VorPairLin_kz",        &VorPairLin_kz,         1.0,           1.0,              NoMax_double      );

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

// (1-3) check the runtime parameters


// (2) set the problem-specific derived parameters
   VorPairLin_kx        *= 2.0*M_PI/amr->BoxSize[0];  // by default we set wavelength equal to multiples of the box size
   VorPairLin_ky        *= 2.0*M_PI/amr->BoxSize[1];
   VorPairLin_kz        *= 2.0*M_PI/amr->BoxSize[2];
   VorPairLin_Omega      = 0.5/ELBDM_ETA*( SQR(VorPairLin_kx) + SQR(VorPairLin_ky) );
   VorPairLin_ZWaveOmega = 0.5/ELBDM_ETA*  SQR(VorPairLin_kz);


// (3) reset other general-purpose parameters
//     --> a helper macro PRINT_RESET_PARA is defined in Macro.h
   const long End_Step_Default = __INT_MAX__;
   double End_T_Default        = 1.0*2.0*M_PI/VorPairLin_Omega;   // 1 period for 2D test

// in 3D test choose 1 period of the z-wave or x-y-wave depending on which is longer
   if ( VorPairLin_ZWaveAmp > 0 ) {
      End_T_Default = FMAX( End_T_Default, 1.0*2.0*M_PI/VorPairLin_ZWaveOmega );
   }

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
      Aux_Message( stdout, "  test problem ID       = %d\n",     TESTPROB_ID           );
      Aux_Message( stdout, "  VorPairLin_BgAmp      = %13.7e\n", VorPairLin_BgAmp      );
      Aux_Message( stdout, "  VorPairLin_WaveAmp    = %13.7e\n", VorPairLin_WaveAmp    );
      Aux_Message( stdout, "  VorPairLin_Phase0     = %13.7e\n", VorPairLin_Phase0     );
      Aux_Message( stdout, "  VorPairLin_ZWaveAmp   = %13.7e\n", VorPairLin_ZWaveAmp   );
      Aux_Message( stdout, "  VorPairLin_kx         = %13.7e\n", VorPairLin_kx         );
      Aux_Message( stdout, "  VorPairLin_ky         = %13.7e\n", VorPairLin_ky         );
      Aux_Message( stdout, "  VorPairLin_kz         = %13.7e\n", VorPairLin_kz         );
      Aux_Message( stdout, "  VorPairLin_Omega      = %13.7e\n", VorPairLin_Omega      );
      Aux_Message( stdout, "  VorPairLin_ZWaveOmega = %13.7e\n", VorPairLin_ZWaveOmega );
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

   const double phase  = VorPairLin_kx*x - VorPairLin_Omega*Time + VorPairLin_Phase0;
   const double amp    = VorPairLin_WaveAmp*cos( VorPairLin_ky*y );
   const double zphase = VorPairLin_kz*z - VorPairLin_ZWaveOmega*Time;
   const double zamp   = VorPairLin_ZWaveAmp;
   const double Re     = VorPairLin_BgAmp + amp*cos( phase ) + zamp*cos( zphase );
   const double Im     =                  + amp*sin( phase ) + zamp*sin( zphase );

   fluid[DENS] = SQR( Re ) + SQR( Im );

#  if ( ELBDM_SCHEME == ELBDM_HYBRID )
   if ( amr->use_wave_flag[lv] ) {
#  endif
   fluid[REAL] = Re;
   fluid[IMAG] = Im;
#  if ( ELBDM_SCHEME == ELBDM_HYBRID )
   } else {
   fluid[PHAS] = SATAN2( Im, Re );
   fluid[STUB] = 0.0;
   }
#  endif

} // FUNCTION : SetGridIC



//-------------------------------------------------------------------------------------------------------
// Function    :  OutputVortexPairLinearError
// Description :  Output the L1 error
//
// Note        :  1. Invoke Output_L1Error()
//                2. Use SetGridIC() to provide the analytical solution at any given time
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void OutputVortexPairLinearError()
{

   const char Prefix[100]     = "VortexPairLinear";
   const OptOutputPart_t Part = OUTPUT_X;

   Output_L1Error( SetGridIC, NULL, Prefix, Part, OUTPUT_PART_X, OUTPUT_PART_Y, OUTPUT_PART_Z );

} // FUNCTION : OutputVortexPairLinearError
#endif // #if ( MODEL == ELBDM )



//-------------------------------------------------------------------------------------------------------
// Function    :  Init_TestProb_ELBDM_VortexPairLinear
// Description :  Test problem initializer
//
// Note        :  None
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Init_TestProb_ELBDM_VortexPairLinear()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );


// validate the compilation flags and runtime parameters
   Validate();


#  if ( MODEL == ELBDM )
// set the problem-specific runtime parameters
   SetParameter();


   Init_Function_User_Ptr    = SetGridIC;
   Output_User_Ptr           = OutputVortexPairLinearError;
#  ifdef SUPPORT_HDF5
   Output_HDF5_InputTest_Ptr = LoadInputTestProb;
#  endif
#  endif // #if ( MODEL == ELBDM )


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : Init_TestProb_ELBDM_VortexPairLinear
