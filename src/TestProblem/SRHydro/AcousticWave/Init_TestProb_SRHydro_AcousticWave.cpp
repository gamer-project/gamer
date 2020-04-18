#include "GAMER.h"
#include "TestProb.h"



static void OutputError();

// problem-specific global variables
// =======================================================================================
static double Acoustic_RhoAmp;      // amplitude of the proper number density perturbation (assuming background density = 1.0)
static double Acoustic_Phase0;      // initial phase shift
static double Acoustic_Temp0;
static double Acoustic_Rho0;
static double Acoustic_WaveLength;  // wavelength
static double SQRT_3;
static double Cs_Sqr;
static double Cs;
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

#  if ( MODEL != SR_HYDRO )
   Aux_Error( ERROR_INFO, "MODEL != SR_HYDRO !!\n" );
#  endif

#  ifdef GRAVITY
   Aux_Error( ERROR_INFO, "GRAVITY must be disabled !!\n" );
#  endif

#  ifndef FLOAT8
   Aux_Error( ERROR_INFO, "FLOAT8 must be enabled !!\n" );
#  endif

#  ifdef COMOVING
   Aux_Error( ERROR_INFO, "COMOVING must be disabled !!\n" );
#  endif

#  ifdef PARTICLE
   Aux_Error( ERROR_INFO, "PARTICLE must be disabled !!\n" );
#  endif

   if ( amr->BoxSize[0] != amr->BoxSize[1]  ||  amr->BoxSize[0] != amr->BoxSize[2] )
      Aux_Error( ERROR_INFO, "simulation domain must be cubic !!\n" );

   for (int f=0; f<6; f++)
   if ( OPT__BC_FLU[f] != BC_FLU_PERIODIC )
      Aux_Error( ERROR_INFO, "please set \"OPT__BC_FLU_* = 1\" (i.e., periodic BC) !!\n" );


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Validating test problem %d ... done\n", TESTPROB_ID );

} // FUNCTION : Validate



#if ( MODEL == SR_HYDRO )
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
   const char FileName[] = "Input__TestProb";
   ReadPara_t *ReadPara  = new ReadPara_t;

// add parameters in the following format:
// --> note that VARIABLE, DEFAULT, MIN, and MAX must have the same data type
// --> some handy constants (e.g., NoMin_int, Eps_float, ...) are defined in "include/ReadPara.h"
// ********************************************************************************************************************************
// ReadPara->Add( "KEY_IN_THE_FILE",   &VARIABLE,              DEFAULT,       MIN,              MAX               );
// ********************************************************************************************************************************
   ReadPara->Add( "Acoustic_RhoAmp",   &Acoustic_RhoAmp,       -1.0,          Eps_double,       NoMax_double      );
   ReadPara->Add( "Acoustic_Phase0",   &Acoustic_Phase0,        0.0,          NoMin_double,     NoMax_double      );
   ReadPara->Add( "Acoustic_Temp0",    &Acoustic_Temp0,         1.0,          Eps_double,       NoMax_double      );
   ReadPara->Add( "Acoustic_Rho0",     &Acoustic_Rho0,          1.0,          Eps_double,       NoMax_double      );

   ReadPara->Read( FileName );

   delete ReadPara;

   SQRT_3 = sqrt(3.0);

   Acoustic_WaveLength = amr->BoxSize[0] / SQRT_3;

  Cs_Sqr = SoundSpeedSquare( Acoustic_Temp0, (double)GAMMA );
  Cs = sqrt(Cs_Sqr);
// (3) reset other general-purpose parameters
//     --> a helper macro PRINT_WARNING is defined in TestProb.h
   const double End_T_Default    = Acoustic_WaveLength * Cs;
   const long   End_Step_Default = __INT_MAX__;

   if ( END_STEP < 0 ) {
      END_STEP = End_Step_Default;
      PRINT_WARNING( "END_STEP", END_STEP, FORMAT_LONG );
   }

   if ( END_T < 0.0 ) {
      END_T = End_T_Default;
      PRINT_WARNING( "END_T", END_T, FORMAT_REAL );
   }


// (4) make a note
   if ( MPI_Rank == 0 )
   {
      Aux_Message( stdout, "=============================================================================\n" );
      Aux_Message( stdout, "  test problem ID     = %d\n",      TESTPROB_ID );
      Aux_Message( stdout, "  density amplitude   = % 14.7e\n", Acoustic_RhoAmp );
      Aux_Message( stdout, "  ambient temperature = % 14.7e\n", Acoustic_Temp0 );
      Aux_Message( stdout, "  ambient density     = % 14.7e\n", Acoustic_Rho0 );
      Aux_Message( stdout, "  initial phase shift = % 14.7e\n", Acoustic_Phase0 );
      Aux_Message( stdout, "  wave length         = % 14.7e\n", Acoustic_WaveLength );
      Aux_Message( stdout, "  sound speed         = % 14.7e\n", Cs );
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
//                3. Even when DUAL_ENERGY is adopted for SR_HYDRO, one does NOT need to set the dual-energy variable here
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
  double Phase;
  double Pri[NCOMP_FLUID], Con[NCOMP_FLUID];
  double K, LorentzFactor;

  LorentzFactor = 1.0 / sqrt( 1.0 - Cs_Sqr );

  K = (double)2.0 * (double)M_PI  / Acoustic_WaveLength;

  Phase = K * ( x + y + z ) / SQRT_3 - K * Cs * Time;

  Pri[0] = Acoustic_Rho0 * ( (double)1.0 + Acoustic_RhoAmp * sin( Phase + Acoustic_Phase0 ) );
  Pri[1] = LorentzFactor * Cs / SQRT_3;
  Pri[2] = LorentzFactor * Cs / SQRT_3;
  Pri[3] = LorentzFactor * Cs / SQRT_3;
  Pri[4] = Acoustic_Temp0 * Acoustic_Rho0;

  SRHydro_Pri2Con ( Pri, Con, (double)GAMMA );

  fluid[0] = (real)Con[0];
  fluid[1] = (real)Con[1];
  fluid[2] = (real)Con[2];
  fluid[3] = (real)Con[3];
  fluid[4] = (real)Con[4];

} // FUNCTION : SetGridIC



//-------------------------------------------------------------------------------------------------------
// Function    :  OutputError
// Description :  Output the L1 error
//
// Note        :  1. Invoke Output_L1Error()
//                2. Use SetGridIC() to provide the analytical solution at any given time
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void OutputError()
{

   const char Prefix[100]     = "AcousticWave";
   const OptOutputPart_t Part = OUTPUT_DIAG;

   Output_L1Error( SetGridIC, Prefix, Part, NULL_REAL, NULL_REAL, NULL_REAL );

} // FUNCTION : OutputError
#endif // #if ( MODEL == SR_HYDRO )



//-------------------------------------------------------------------------------------------------------
// Function    :  Init_TestProb_Hydro_AcousticWave
// Description :  Test problem initializer
//
// Note        :  None
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Init_TestProb_SRHydro_AcousticWave()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );


// validate the compilation flags and runtime parameters
   Validate();


#  if ( MODEL == SR_HYDRO )
// set the problem-specific runtime parameters
   SetParameter();


// set the function pointers of various problem-specific routines
   Init_Function_User_Ptr   = SetGridIC;
   Output_User_Ptr          = OutputError;
   Flag_User_Ptr            = NULL;
   Mis_GetTimeStep_User_Ptr = NULL;
   Aux_Record_User_Ptr      = NULL;
   BC_User_Ptr              = NULL;
   Flu_ResetByUser_Func_Ptr = NULL;
   End_User_Ptr             = NULL;
#  endif // #if ( MODEL == SR_HYDRO )


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : Init_TestProb_Hydro_AcousticWave
