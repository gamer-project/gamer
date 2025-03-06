#include "GAMER.h"



// problem-specific global variables
// =======================================================================================
static double Gau_v0;       // mean velocity
static double Gau_Width;    // Gaussian width
static double Gau_Center;   // Gaussian center
static int    Gau_XYZ;      // wave propagation direction (0/1/2 --> x/y/z)
static int    Gau_PeriodicN;// periodic boundary condition
                            // (0 = non-periodic, >0 = number of periodic images each side)
// =======================================================================================

// defined TargetVariable for AnalyticalSolution_GaussianWavePacket()
// =======================================================================================
const  int    GAUSSIAN_TARGET_VAR_REAL = 1;
const  int    GAUSSIAN_TARGET_VAR_IMAG = 2;
const  int    GAUSSIAN_TARGET_VAR_PHAS = 3;
// =======================================================================================

static double AnalyticalSolution_GaussianWavePacket( const double r, const double Time,
                                                     const int TargetVariable );
static void   ComplexValuesAmpPhaseAdder( const double A_1, const double S_1,
                                          const double A_2, const double S_2,
                                          double* A_sum, double* S_sum );
static void   OutputError();
static void   BC( real Array[], const int ArraySize[], real fluid[], const int NVar_Flu,
                  const int GhostSize, const int idx[], const double pos[], const double Time,
                  const int lv, const int TFluVarIdxList[], double AuxArray[] );




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
// ************************************************************************************************************************
// LOAD_PARA( load_mode, "KEY_IN_THE_FILE",   &VARIABLE,              DEFAULT,       MIN,              MAX               );
// ************************************************************************************************************************
   LOAD_PARA( load_mode, "Gau_v0",            &Gau_v0,                1.0,           NoMin_double,     NoMax_double      );
   LOAD_PARA( load_mode, "Gau_Width",         &Gau_Width,             0.1,           Eps_double,       NoMax_double      );
   LOAD_PARA( load_mode, "Gau_Center",        &Gau_Center,            NoDef_double,  NoMin_double,     NoMax_double      );
   LOAD_PARA( load_mode, "Gau_XYZ",           &Gau_XYZ,               0,             0,                2                 );
   LOAD_PARA( load_mode, "Gau_PeriodicN",     &Gau_PeriodicN,         0,             0,                NoMax_int         );

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
   if ( Gau_Center == NoDef_double )   Gau_Center = amr->BoxCenter[Gau_XYZ];

// (1-3) check the runtime parameters
   if ( OPT__BC_FLU[2*Gau_XYZ] == BC_FLU_PERIODIC  &&  Gau_PeriodicN == 0 )
      Aux_Error( ERROR_INFO, "Gau_PeriodicN should be >0 when adopting periodic BC for fluid !!\n" );

#  if ( ELBDM_SCHEME == ELBDM_HYBRID )
   if ( MPI_Rank == 0 )
   {
      if ( Gau_PeriodicN > 0  &&  END_T > 6.0*Gau_Width*ELBDM_ETA*amr->BoxSize[Gau_XYZ]/(2*M_PI) )
         Aux_Message( stderr, "WARNING : END_T = %13.6e, but the analytical solution of unwrapped phase may be wrong after t > %13.6e !!\n",
                              END_T, 6.0*Gau_Width*ELBDM_ETA*amr->BoxSize[Gau_XYZ]/(2*M_PI) );

      if ( Gau_PeriodicN > 0  &&  Gau_PeriodicN < 3 )
         Aux_Message( stderr, "WARNING : Gau_PeriodicN = %d may be too few for a smooth analytical solution of phase in the hybrid scheme !!\n", Gau_PeriodicN );
   }

   if ( ELBDM_FIRST_WAVE_LEVEL > MAX_LEVEL )
   {
//    the periodic boundary condition is incompatible with the fluid scheme in general
      if ( OPT__BC_FLU[2*Gau_XYZ] == BC_FLU_PERIODIC )
      {
         if ( Gau_v0 != 0.0 )
            Aux_Error( ERROR_INFO, "OPT__BC_FLU[%d] must be %d for travelling Gaussian wave packet in the fluid scheme !!\n", 2*Gau_XYZ, BC_FLU_USER );

         else if ( Gau_PeriodicN == 0  &&  ! Mis_CompareRealValue( Gau_Center, amr->BoxCenter[Gau_XYZ], NULL, false ) )
            Aux_Error( ERROR_INFO, "OPT__BC_FLU[%d] must be %d for non-periodic Gaussian wave packet in the fluid scheme !!\n", 2*Gau_XYZ, BC_FLU_USER );
      }
   }
   else // if ( ELBDM_FIRST_WAVE_LEVEL > MAX_LEVEL )
   {
//    the phase on the wave level is wrapped
      if ( !ELBDM_MATCH_PHASE )
         Aux_Error( ERROR_INFO, "ELBDM_MATCH_PHASE should be enabled to make the phase on fluid levels continuous !!\n" );
   } // if ( ELBDM_FIRST_WAVE_LEVEL > MAX_LEVEL ) ... else
#  endif // # if ( ELBDM_SCHEME == ELBDM_HYBRID )


// (2) set the problem-specific derived parameters


// (3) reset other general-purpose parameters
//     --> a helper macro PRINT_RESET_PARA is defined in Macro.h
   const long   End_Step_Default = __INT_MAX__;
   const double End_T_Default    = 0.50*amr->BoxSize[Gau_XYZ]/fabs( Gau_v0 );

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
      Aux_Message( stdout, "  test problem ID                    = %d\n",     TESTPROB_ID   );
      Aux_Message( stdout, "  mean velocity                      = %14.7e\n", Gau_v0        );
      Aux_Message( stdout, "  Gaussian width                     = %14.7e\n", Gau_Width     );
      Aux_Message( stdout, "  Gaussian center                    = %14.7e\n", Gau_Center    );
      Aux_Message( stdout, "  propagation direction              = %d\n",     Gau_XYZ       );
      Aux_Message( stdout, "  number of periodic image each side = %d\n",     Gau_PeriodicN );
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

   double r;
   switch ( Gau_XYZ )
   {
      case 0: r = x;    break;
      case 1: r = y;    break;
      case 2: r = z;    break;

      default : Aux_Error( ERROR_INFO, "incorrect Gau_XYZ (%d) !!\n", Gau_XYZ );
      break;
   }

   const double Re = AnalyticalSolution_GaussianWavePacket( r, Time, GAUSSIAN_TARGET_VAR_REAL );
   const double Im = AnalyticalSolution_GaussianWavePacket( r, Time, GAUSSIAN_TARGET_VAR_IMAG );

   fluid[DENS] = SQR(Re) + SQR(Im);

#  if ( ELBDM_SCHEME == ELBDM_HYBRID )
   if ( amr->use_wave_flag[lv] ) {
#  endif
   fluid[REAL] = Re;
   fluid[IMAG] = Im;
#  if  ( ELBDM_SCHEME == ELBDM_HYBRID )
   } else {
   fluid[PHAS] = AnalyticalSolution_GaussianWavePacket( r, Time, GAUSSIAN_TARGET_VAR_PHAS );
   fluid[STUB] = 0.0;
   }
#  endif

} // FUNCTION : SetGridIC



//-------------------------------------------------------------------------------------------------------
// Function    :  AnalyticalSolution_GaussianWavePacket
// Description :  Analytical solution of Gaussian wave packet as a function of position and time
//
// Note        :  1. TargetVariable = GAUSSIAN_TARGET_VAR_REAL --> the real part of the wave function
//                   TargetVariable = GAUSSIAN_TARGET_VAR_IMAG --> the imaginary part of the wave function
//                   TargetVariable = GAUSSIAN_TARGET_VAR_PHAS --> the unwrapped phase of the wave function
//                2. After Time > ( 6.0*Gau_Width/BoxSize )*( BoxSize/(2*M_PI/(ELBDM_ETA*fabs(Gau_v0))) )*( BoxSize/fabs(Gau_v0) ),
//                   the analytical solution of unwrapped phase would be wrong due to the complicated interference between images
//                   -> Phase unwrapping with the previous time will be applied to avoid this effect, assuming the time evolution of phase is continuous,
//                      but it may take a longer time to trace the phase recursively and still may be wrong
//
// Parameter   :  r                : Physical coordinates
//                Time             : Physical time
//                TargetVariable   : Target variable to return
//
// Return      :  Re or Im or Ph
//-------------------------------------------------------------------------------------------------------
double AnalyticalSolution_GaussianWavePacket( const double r, const double Time, const int TargetVariable )
{

   const double Gau_Const1 = 1.0 + pow(  Time / ( ELBDM_ETA*SQR(Gau_Width) ), 2.0  );
   const double Gau_Theta1 = -0.5*acos(  pow( Gau_Const1, -0.5 )  );
   double Re=0.0, Im=0.0, Am=0.0, Ph=0.0, Am_main=0.0, Ph_main=0.0;

// n=0, m=0: original wave packet
// n>0, m=0/1: images for periodic BC on the plus(+)/minus(-) direction
   for (int n=0; n<Gau_PeriodicN+1; n++) {
   for (int m=0; m<((n==0)?1:2); m++) {
      const double Center     = Gau_Center + n*(1-2*m)*amr->BoxSize[Gau_XYZ];
      const double dr1        = r -     Gau_v0*Time - Center;
      const double dr2        = r - 0.5*Gau_v0*Time - Center;
      const double Gau_Const2 = pow( SQR(Gau_Width)*M_PI*Gau_Const1, -0.25 )
                                *exp(  -0.5*pow( dr1/Gau_Width, 2.0 )/Gau_Const1  );
      const double Gau_Theta2 = 0.5*pow( dr1, 2.0 )*ELBDM_ETA*Time/(  pow( ELBDM_ETA*SQR(Gau_Width), 2.0) + SQR(Time)  )
                                + Gau_v0*ELBDM_ETA*dr2;

      if      ( TargetVariable == GAUSSIAN_TARGET_VAR_REAL )   Re += Gau_Const2*cos( Gau_Theta1 + Gau_Theta2 );
      else if ( TargetVariable == GAUSSIAN_TARGET_VAR_IMAG )   Im += Gau_Const2*sin( Gau_Theta1 + Gau_Theta2 );
      else if ( TargetVariable == GAUSSIAN_TARGET_VAR_PHAS )
      {
//       shift the phase to make sure the phase is continuous at the transition of different Gaussians
         const double PhaseShift = 2*M_PI*round( Gau_v0*ELBDM_ETA*amr->BoxSize[Gau_XYZ]/(2*M_PI) )*n*(1-2*m);
         const double Am_this    = Gau_Const2;
         const double Ph_this    = Gau_Theta1 + Gau_Theta2 + PhaseShift;

//       check if this Gaussian is at the two ends
         const bool isRightMost  = ( n == Gau_PeriodicN  &&  m == 0                        );
         const bool isLeftMost   = ( n == Gau_PeriodicN  &&  m == ((Gau_PeriodicN==0)?0:1) );

//       the Gaussian with the largest amplitude locally is the main Gaussian
         if ( (                  round( (r-(Gau_Center+Gau_v0*Time))/amr->BoxSize[Gau_XYZ] ) == n*(1-2*m) ) ||
              ( isRightMost  &&  round( (r-(Gau_Center+Gau_v0*Time))/amr->BoxSize[Gau_XYZ] ) >  n*(1-2*m) ) ||
              ( isLeftMost   &&  round( (r-(Gau_Center+Gau_v0*Time))/amr->BoxSize[Gau_XYZ] ) <  n*(1-2*m) ) )
         {
//          keep the amplitude and phase of the main Gaussian to be added lastly
            Am_main = Am_this;
            Ph_main = Ph_this;
         }
         else
         {
//          add the phase if it is not the main Gaussian
            ComplexValuesAmpPhaseAdder( Am_this, Ph_this, Am, Ph, &Am, &Ph );
         }
      }
   }}

   if ( TargetVariable == GAUSSIAN_TARGET_VAR_PHAS )
   {
//    add the main phase lastly to make sure the final phase is close to it
      ComplexValuesAmpPhaseAdder( Am_main, Ph_main, Am, Ph,  &Am,  &Ph );

//    shift the phase to make sure its time evolution is continuous
      if ( Time > 6.0*Gau_Width*ELBDM_ETA*amr->BoxSize[Gau_XYZ]/(2*M_PI) )
      {
//        the delta_t is chosen to make the phase changes less than pi each step in the phase equation
          const double delta_t     = MIN( 2.0*M_PI/( ELBDM_ETA*SQR(Gau_v0) ), 2.0*M_PI*ELBDM*SQR(Gau_Width)*SQR(2*Gau_Width/amr->BoxSize[Gau_XYZ]) );
          const double Ph_previous = AnalyticalSolution_GaussianWavePacket( r, Time-delta_t, GAUSSIAN_TARGET_VAR_PHAS );
          Ph -= 2*M_PI*floor( ((Ph - Ph_previous)+M_PI)/(2*M_PI) );
      }
   }

// return the result
   if      ( TargetVariable == GAUSSIAN_TARGET_VAR_REAL )   return Re;
   else if ( TargetVariable == GAUSSIAN_TARGET_VAR_IMAG )   return Im;
   else if ( TargetVariable == GAUSSIAN_TARGET_VAR_PHAS )   return Ph;

   return 0.0;

} // FUNCTION : AnalyticalSolution_GaussianWavePacket



//-------------------------------------------------------------------------------------------------------
// Function    :  ComplexValuesAmpPhaseAdder
// Description :  Add two complex variables by their amplitudes and phases
//
// Note        :  1. The phase will not be wrapped to -pi to pi after addition.
//                2. The order of two inputs matters, S_sum will be close to S_1, not S_2.
//
// Parameter   :  A_1      : Amplitude of the first complex variable
//                S_1      : Phase of the first complex variable
//                A_2      : Phase of the second complex variable
//                S_2      : Phase of the second complex variable
//
// Return      :  A_sun, S_sum
//-------------------------------------------------------------------------------------------------------
void ComplexValuesAmpPhaseAdder( const double A_1, const double S_1, const double A_2, const double S_2, double* A_sum, double* S_sum )
{

// the summed unwrapped phase follows the phase of the first one, S_a
   const double A_a = A_1;
   const double A_b = A_2;
   const double S_a = S_1;
   const double S_b = S_2;

// find the phase difference
// and subtract two pi to make the relative in the range from -pi to pi
   const double S_difference = S_b - S_a;
   const int    N_winding    = floor( (S_difference+M_PI)/(2*M_PI) );
   const double S_relative   = S_difference - 2*M_PI*N_winding;

// convert to real part and imaginary part
   const double R_a = A_a*cos(        0.0 );
   const double I_a = A_a*sin(        0.0 );
   const double R_b = A_b*cos( S_relative );
   const double I_b = A_b*sin( S_relative );

// add the real parts and imaginary parts
   const double R_c = R_a + R_b;
   const double I_c = I_a + I_b;

// convert back to phase (from -pi to pi )
   const double S_c = atan2( I_c, R_c );

// add the phase back to make it close to S_a
   *A_sum = sqrt( R_c*R_c + I_c*I_c );
   *S_sum = S_a + S_c;

} // FUNCTION : ComplexValuesAmpPhaseAdder



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

   const char Prefix[100]     = "Gaussian";
   const OptOutputPart_t Part = OUTPUT_X + Gau_XYZ;

   Output_L1Error( SetGridIC, NULL, Prefix, Part, OUTPUT_PART_X, OUTPUT_PART_Y, OUTPUT_PART_Z );

} // FUNCTION : OutputError



//-------------------------------------------------------------------------------------------------------
// Function    :  BC
// Description :  Set the extenral boundary condition to the analytical solution
//
// Note        :  1. Linked to the function pointer "BC_User_Ptr"
//
// Parameter   :  Array          : Array to store the prepared data including ghost zones
//                ArraySize      : Size of Array including the ghost zones on each side
//                fluid          : Fluid fields to be set
//                NVar_Flu       : Number of fluid variables to be prepared
//                GhostSize      : Number of ghost zones
//                idx            : Array indices
//                pos            : Physical coordinates
//                Time           : Physical time
//                lv             : Refinement level
//                TFluVarIdxList : List recording the target fluid variable indices ( = [0 ... NCOMP_TOTAL-1] )
//                AuxArray       : Auxiliary array
//
// Return      :  fluid
//-------------------------------------------------------------------------------------------------------
void BC( real Array[], const int ArraySize[], real fluid[], const int NVar_Flu,
         const int GhostSize, const int idx[], const double pos[], const double Time,
         const int lv, const int TFluVarIdxList[], double AuxArray[] )
{

// simply call the IC function
   SetGridIC( fluid, pos[0], pos[1], pos[2], Time, lv, AuxArray );

} // FUNCTION : BC
#endif // #if ( MODEL == ELBDM )



//-------------------------------------------------------------------------------------------------------
// Function    :  Init_TestProb_ELBDM_GaussianWavePacket
// Description :  Test problem initializer
//
// Note        :  None
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Init_TestProb_ELBDM_GaussianWavePacket()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );


// validate the compilation flags and runtime parameters
   Validate();


#  if ( MODEL == ELBDM )
// set the problem-specific runtime parameters
   SetParameter();


   Init_Function_User_Ptr    = SetGridIC;
   BC_User_Ptr               = BC;
   Output_User_Ptr           = OutputError;
#  ifdef SUPPORT_HDF5
   Output_HDF5_InputTest_Ptr = LoadInputTestProb;
#  endif
#  endif // #if ( MODEL == ELBDM )


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : Init_TestProb_ELBDM_GaussianWavePacket
