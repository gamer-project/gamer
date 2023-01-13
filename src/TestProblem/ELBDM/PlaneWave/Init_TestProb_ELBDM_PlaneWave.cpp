#include "GAMER.h"
#include "TestProb.h"



// problem-specific global variables
// =======================================================================================
static double PWave_Lambda;        // plane wave wavelength
static double PWave_Amp;           // plane wave amplitude 
static double PWave_Phase0;        // plane wave phase constant
static int    PWave_XYZ;           // plane wave direction (0/1/2/3 --> x/y/z/diagonal)
static int    PWave_LSR;           // plane wave direction (<0/0/>0 --> negative direction/standing wave/positive direction)

static double PWave_Period;        // plane wave period
static double PWave_WaveK;         // plane wave wavenumber
static double PWave_WaveW;         // plane wave angular frequency 
static double PWave_WaveV;         // plane wave velocity 
// =======================================================================================

static void OutputError();




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

#  ifndef FLOAT8
   Aux_Error( ERROR_INFO, "FLOAT8 must be enabled !!\n" );
#  endif


// warnings
   if ( MPI_Rank == 0 )
   {
      if ( !OPT__OUTPUT_USER )
         Aux_Message( stderr, "WARNING : it's recommended to enable OPT__OUTPUT_USER !!\n" );
   }


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
   ReadPara->Add( "PWave_Lambda",      &PWave_Lambda,          1.0,           Eps_double,       NoMax_double      );
   ReadPara->Add( "PWave_Amp",         &PWave_Amp,             1.0,           0.0,              NoMax_double      );
   ReadPara->Add( "PWave_Phase0",      &PWave_Phase0,          0.0,           NoMin_double,     NoMax_double      );
   ReadPara->Add( "PWave_XYZ",         &PWave_XYZ,             0,             0,                3                 );
   ReadPara->Add( "PWave_LSR",         &PWave_LSR,             1,             NoMin_int,        NoMax_int         );

   ReadPara->Read( FileName );

   delete ReadPara;

// (1-2) set the default values

// (1-3) check the runtime parameters
   if (  PWave_XYZ == 3  &&  ( amr->BoxSize[0] != amr->BoxSize[1] || amr->BoxSize[0] != amr->BoxSize[2] )  )
      Aux_Error( ERROR_INFO, "simulation domain must be CUBIC in the %s test if PWave_XYZ == 3 !!\n", "ELBDM plane wave" );


// (2) set the problem-specific derived parameters
   PWave_WaveK  = 2.0*M_PI/PWave_Lambda;
   PWave_WaveW  = 0.5*SQR( PWave_WaveK )/ELBDM_ETA;
   PWave_Period = 2.0*M_PI/PWave_WaveW;
   PWave_WaveV  = PWave_WaveK/ELBDM_ETA;



// (3) reset other general-purpose parameters
//     --> a helper macro PRINT_WARNING is defined in TestProb.h
   const long   End_Step_Default = __INT_MAX__;
   const double End_T_Default    = 10.0*PWave_Period;

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
      Aux_Message( stdout, "  test problem ID                  = %d\n",     TESTPROB_ID                );
      Aux_Message( stdout, "  plane wave wavelength            = %13.7e\n", PWave_Lambda               );
      Aux_Message( stdout, "  plane wave amplitude             = %13.7e\n", PWave_Amp                  );
      Aux_Message( stdout, "  plane wave phase constant        = %13.7e\n", PWave_Phase0               );
      Aux_Message( stdout, "  plane wave wavenumber            = %13.7e\n", PWave_WaveK                );
      Aux_Message( stdout, "  plane wave angular frequency     = %13.7e\n", PWave_WaveW                );
      Aux_Message( stdout, "  plane wave period                = %13.7e\n", PWave_Period               );
      if ( PWave_LSR == 0){
      Aux_Message( stdout, "  standing wave                    = true\n",                              );
      }
      else{
      Aux_Message( stdout, "  standing wave                    = false\n",                             );
      Aux_Message( stdout, "  plane wave propagation direction = %s%s\n",  ( PWave_LSR > 0 )  ? "+" : "-", 
                                                                           ( PWave_XYZ == 0 ) ? "x" :
                                                                           ( PWave_XYZ == 1 ) ? "y" :
                                                                           ( PWave_XYZ == 2 ) ? "z" : "diagonal" );
      Aux_Message( stdout, "  plane wave propagation velocity  = %13.7e\n", PWave_WaveV                );
      }
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

   double r, PhaseR, PhaseL;
   switch ( PWave_XYZ )
   {
      case 0 : r = x;                              break;
      case 1 : r = y;                              break;
      case 2 : r = z;                              break;
      case 3 : r = 1.0/sqrt(3.0)*( x + y + z );    break;
      default : Aux_Error( ERROR_INFO, "incorrect parameter \"%s = %d [0/1/2/3]\" !!\n", "PWave_XYZ", PWave_XYZ );
      break;
   }

   PhaseR  =  1.0*PWave_WaveK*r - PWave_WaveW*Time + PWave_Phase0;
   PhaseL  = -1.0*PWave_WaveK*r - PWave_WaveW*Time + PWave_Phase0;

// set the real and imaginary parts
   if ( PWave_LSR > 0 ){      // Right-moving wave
      fluid[REAL] = 0.5*PWave_Amp*cos( PhaseR ) + 0.5*PWave_Amp*cos( PhaseR );
      fluid[IMAG] = 0.5*PWave_Amp*sin( PhaseR ) + 0.5*PWave_Amp*sin( PhaseR );
   }
   else if ( PWave_LSR < 0 ){ // Left-moving wave
      fluid[REAL] = 0.5*PWave_Amp*cos( PhaseL ) + 0.5*PWave_Amp*cos( PhaseL );
      fluid[IMAG] = 0.5*PWave_Amp*sin( PhaseL ) + 0.5*PWave_Amp*sin( PhaseL );
   }
   else{ //( PWave_LSR == 0 ) // Standing wave
      fluid[REAL] = 0.5*PWave_Amp*cos( PhaseR ) + 0.5*PWave_Amp*cos( PhaseL );
      fluid[IMAG] = 0.5*PWave_Amp*sin( PhaseR ) + 0.5*PWave_Amp*sin( PhaseL );
   }

// set the density
   fluid[DENS] = SQR( fluid[REAL] ) + SQR( fluid[IMAG] );

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

   const char Prefix[100]     = "PlaneWave";
   const OptOutputPart_t Part = ( PWave_XYZ == 3) ? OUTPUT_DIAG : (OUTPUT_X + PWave_XYZ);

   Output_L1Error( SetGridIC, NULL, Prefix, Part, 0.0, 0.0, 0.0 );

} // FUNCTION : OutputError
#endif // #if ( MODEL == ELBDM )



//-------------------------------------------------------------------------------------------------------
// Function    :  Init_TestProb_ELBDM_PlaneWave
// Description :  Test problem initializer
//
// Note        :  None
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Init_TestProb_ELBDM_PlaneWave()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );


// validate the compilation flags and runtime parameters
   Validate();


#  if ( MODEL == ELBDM )
// set the problem-specific runtime parameters
   SetParameter();


   Init_Function_User_Ptr = SetGridIC;
   BC_User_Ptr            = SetGridIC;
   Output_User_Ptr        = OutputError;
#  endif // #if ( MODEL == ELBDM )


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : Init_TestProb_ELBDM_PlaneWave
