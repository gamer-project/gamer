#include "GAMER.h"
#include "TestProb.h"



// problem-specific global variables
// =======================================================================================
static double VorPairRot_BgAmp;     // psi(R,phi) = BgAmp - J1Amp*J1( sqrt(2*Eta*Omega)*R )*exp( i*(phi-Omega*t+Phase0) )
static double VorPairRot_J1Amp;
static double VorPairRot_Omega;
static double VorPairRot_Phase0;
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

#  if ( ELBDM_SCHEME == ELBDM_HYBRID )
   if ( OPT__HYBRID_MATCH_PHASE )
      Aux_Error( ERROR_INFO, "OPT__HYBRID_MATCH_PHASE must be disabled !!\n" );
#  endif // #  if ( ELBDM_SCHEME == ELBDM_HYBRID )

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
// ReadPara->Add( "KEY_IN_THE_FILE",   &VARIABLE,              DEFAULT,       MIN,              MAX               );
// ********************************************************************************************************************************
   ReadPara->Add( "VorPairRot_BgAmp",  &VorPairRot_BgAmp,     -1.0,           Eps_double,       NoMax_double      );
   ReadPara->Add( "VorPairRot_J1Amp",  &VorPairRot_J1Amp,     -1.0,           Eps_double,       NoMax_double      );
   ReadPara->Add( "VorPairRot_Omega",  &VorPairRot_Omega,     -1.0,           Eps_double,       NoMax_double      );
   ReadPara->Add( "VorPairRot_Phase0", &VorPairRot_Phase0,     0.0,           NoMin_double,     NoMax_double      );

   ReadPara->Read( FileName );

   delete ReadPara;

// (1-2) set the default values

// (1-3) check the runtime parameters


// (2) set the problem-specific derived parameters


// (3) reset other general-purpose parameters
//     --> a helper macro PRINT_WARNING is defined in TestProb.h
   const long   End_Step_Default = __INT_MAX__;
   const double End_T_Default    = 1.0*2.0*M_PI/VorPairRot_Omega;    // 1 period

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
      Aux_Message( stdout, "  test problem ID   = %d\n",     TESTPROB_ID       );
      Aux_Message( stdout, "  VorPairRot_BgAmp  = %13.7e\n", VorPairRot_BgAmp  );
      Aux_Message( stdout, "  VorPairRot_J1Amp  = %13.7e\n", VorPairRot_J1Amp  );
      Aux_Message( stdout, "  VorPairRot_Omega  = %13.7e\n", VorPairRot_Omega  );
      Aux_Message( stdout, "  VorPairRot_Phase0 = %13.7e\n", VorPairRot_Phase0 );
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

   const double dx    = x - amr->BoxCenter[0];
   const double dy    = y - amr->BoxCenter[1];
   const double phase = atan2( dy, dx ) - VorPairRot_Omega*Time + VorPairRot_Phase0;
   const double R     = sqrt( SQR(dx) + SQR(dy) );
   const double J1    = VorPairRot_J1Amp*j1( sqrt(2.0*ELBDM_ETA*VorPairRot_Omega)*R );
   const double Re    = VorPairRot_BgAmp - J1*cos( phase );
   const double Im    =                  - J1*sin( phase );

   fluid[DENS] = SQR( Re ) + SQR( Im );

#  if ( ELBDM_SCHEME == ELBDM_HYBRID )
   if ( amr->use_wave_flag[lv] ) {
#  endif
   fluid[REAL] = Re;
   fluid[IMAG] = Im;
#  if ( ELBDM_SCHEME == ELBDM_HYBRID )
   } else { // if ( amr->use_wave_flag[lv] )
   fluid[PHAS] = SATAN2(Im, Re);
   //if (fluid[PHAS] > M_PI / 2)
   //   printf("Coords %f %f %f Phase %f", x, y, z, fluid[PHAS]);
   fluid[STUB] = 0.0;
   } // if ( amr->use_wave_flag[lv] ) ... else
#  endif

} // FUNCTION : SetGridIC


//-------------------------------------------------------------------------------------------------------
// Function    :  OutputVortexPairRotatingError
// Description :  Output the L1 error
//
// Note        :  1. Invoke Output_L1Error()
//                2. Use SetGridIC() to provide the analytical solution at any given time
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void OutputVortexPairRotatingError()
{

   const char Prefix[100]     = "VortexPairRotating";
   const OptOutputPart_t Part = OUTPUT_X;

   Output_L1Error( SetGridIC, NULL, Prefix, Part, OUTPUT_PART_X, OUTPUT_PART_Y, OUTPUT_PART_Z );

} // FUNCTION : OutputVortexPairRotatingError

#endif // #if ( MODEL == ELBDM )

//-------------------------------------------------------------------------------------------------------
// Function    :  Init_TestProb_ELBDM_VortexPairRotating
// Description :  Test problem initializer
//
// Note        :  None
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Init_TestProb_ELBDM_VortexPairRotating()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );


// validate the compilation flags and runtime parameters
   Validate();


#  if ( MODEL == ELBDM )
// set the problem-specific runtime parameters
   SetParameter();


   Init_Function_User_Ptr = SetGridIC;
   BC_User_Ptr            = SetGridIC;
   Output_User_Ptr        = OutputVortexPairRotatingError;
#  endif // #if ( MODEL == ELBDM )


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : Init_TestProb_ELBDM_VortexPairRotating
