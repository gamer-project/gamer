#include "GAMER.h"
#include "TestProb.h"



// problem-specific global variables
// =======================================================================================
static double VorPairLin_BgAmp;     // psi(x,y) = BgAmp + WaveAmp*cos(ky*y)*exp( i*(kx*x-Omega*t+Phase0) )
static double VorPairLin_WaveAmp;
static double VorPairLin_Phase0;

static double VorPairLin_kx;
static double VorPairLin_ky;
static double VorPairLin_Omega;
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


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Validating test problem %d ... done\n", TESTPROB_ID );

} // FUNCTION : Validate



#if ( MODEL == ELBDM )
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

// (1-1) add parameters in the following format:
// --> note that VARIABLE, DEFAULT, MIN, and MAX must have the same data type
// --> some handy constants (e.g., Useless_bool, Eps_double, NoMin_int, ...) are defined in "include/ReadPara.h"
// ********************************************************************************************************************************
// ReadPara->Add( "KEY_IN_THE_FILE",    &VARIABLE,              DEFAULT,       MIN,              MAX               );
// ********************************************************************************************************************************
   ReadPara->Add( "VorPairLin_BgAmp",   &VorPairLin_BgAmp,     -1.0,           Eps_double,       NoMax_double      );
   ReadPara->Add( "VorPairLin_WaveAmp", &VorPairLin_WaveAmp,   -1.0,           Eps_double,       NoMax_double      );
   ReadPara->Add( "VorPairLin_Phase0",  &VorPairLin_Phase0,     0.0,           NoMin_double,     NoMax_double      );

   ReadPara->Read( FileName );

   delete ReadPara;

// (1-2) set the default values

// (1-3) check the runtime parameters


// (2) set the problem-specific derived parameters
   VorPairLin_kx    = 2.0*M_PI/amr->BoxSize[0];   // by default we set wavelength equal to the box size
   VorPairLin_ky    = 2.0*M_PI/amr->BoxSize[1];
   VorPairLin_Omega = 0.5/ELBDM_ETA*( SQR(VorPairLin_kx) + SQR(VorPairLin_ky) );


// (3) reset other general-purpose parameters
//     --> a helper macro PRINT_WARNING is defined in TestProb.h
   const long   End_Step_Default = __INT_MAX__;
   const double End_T_Default    = 1.0*2.0*M_PI/VorPairLin_Omega;    // 1 period

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
      Aux_Message( stdout, "  test problem ID    = %d\n",     TESTPROB_ID        );
      Aux_Message( stdout, "  VorPairLin_BgAmp   = %13.7e\n", VorPairLin_BgAmp   );
      Aux_Message( stdout, "  VorPairLin_WaveAmp = %13.7e\n", VorPairLin_WaveAmp );
      Aux_Message( stdout, "  VorPairLin_Phase0  = %13.7e\n", VorPairLin_Phase0  );
      Aux_Message( stdout, "  VorPairLin_kx      = %13.7e\n", VorPairLin_kx      );
      Aux_Message( stdout, "  VorPairLin_ky      = %13.7e\n", VorPairLin_ky      );
      Aux_Message( stdout, "  VorPairLin_Omega   = %13.7e\n", VorPairLin_Omega   );
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

   const double phase = VorPairLin_kx*x - VorPairLin_Omega*Time + VorPairLin_Phase0;
   const double amp   = VorPairLin_WaveAmp*cos( VorPairLin_ky*y );

   fluid[REAL] = VorPairLin_BgAmp + amp*cos( phase );
   fluid[IMAG] =                  + amp*sin( phase );
   fluid[DENS] = SQR( fluid[REAL] ) + SQR( fluid[IMAG] );

} // FUNCTION : SetGridIC
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


   Init_Function_User_Ptr      = SetGridIC;
   Init_Field_User_Ptr         = NULL;
   Flag_User_Ptr               = NULL;
   Mis_GetTimeStep_User_Ptr    = NULL;
   BC_User_Ptr                 = NULL;
   Flu_ResetByUser_Func_Ptr    = NULL;
   Output_User_Ptr             = NULL;
   Aux_Record_User_Ptr         = NULL;
   End_User_Ptr                = NULL;
#  endif // #if ( MODEL == ELBDM )


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : Init_TestProb_ELBDM_VortexPairLinear
