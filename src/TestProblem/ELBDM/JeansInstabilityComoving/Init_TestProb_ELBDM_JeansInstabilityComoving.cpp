#include "GAMER.h"



// problem-specific global variables
// =======================================================================================
static double Jeans_RealAmp0;       // proportional coefficient of the real part of the wave function
static double Jeans_Phase0;         // initial phase shift

static double Jeans_ImagAmp0;       // proportional coefficient of the imaginary part of the wave function
static double Jeans_Wavelength;     // wavelength
static double Jeans_WaveK;          // wavenumber
static double Jeans_WaveKj;         // critical wavenumber
static bool   Jeans_Stable;         // true/false --> Jeans stable/unstable
// =======================================================================================



// inline functions to calculate the amplitudes of the real and imaginary parts
// =======================================================================================
inline double Jeans_RealAmp( const double RealAmp0, const double y )
{
   const double y2 = y*y;

   return RealAmp0 * ( 3.0*cos(y) + 3.0*y*sin(y) - y2*cos(y) ) / y2;
}

inline double Jeans_ImagAmp( const double ImagAmp0, const double y )
{
   const double y2 = y*y;
   const double y3 = y*y2;
   const double y4 = y2*y2;

   return -ImagAmp0 * ( -6.0*y*cos(y) - 6.0*y2*sin(y) + 3.0*y3*cos(y) + y4*sin(y) ) / y4;
}
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

#  ifndef GRAVITY
   Aux_Error( ERROR_INFO, "GRAVITY must be enabled !!\n" );
#  endif

#  ifndef COMOVING
   Aux_Error( ERROR_INFO, "COMOVING must be enabled !!\n" );
#  endif

#  ifdef PARTICLE
   Aux_Error( ERROR_INFO, "PARTICLE must be disabled !!\n" );
#  endif

   for (int f=0; f<6; f++)
   if ( OPT__BC_FLU[f] != BC_FLU_PERIODIC )
      Aux_Error( ERROR_INFO, "must adopt periodic BC for fluid --> reset OPT__BC_FLU* !!\n" );

#  ifdef GRAVITY
   if ( OPT__BC_POT != BC_POT_PERIODIC )
      Aux_Error( ERROR_INFO, "must adopt periodic BC for gravity --> reset OPT__BC_POT !!\n" );
#  endif

#  ifdef COMOVING
   if ( OMEGA_M0 != 1.0 )
      Aux_Error( ERROR_INFO, "must adopt \"OMEGA_M0 = 1.0\" !!\n" );
#  endif

   if ( amr->BoxSize[0] != amr->BoxSize[1]  ||  amr->BoxSize[0] != amr->BoxSize[2] )
      Aux_Error( ERROR_INFO, "simulation domain must be cubic --> reset BOX_SIZE !!\n" );


// warnings
   if ( MPI_Rank == 0 )
   {
#     ifndef FLOAT8
      Aux_Message( stderr, "WARNING : it's recommended to enable FLOAT8 for this test !!\n" );
#     endif

      if ( !OPT__OUTPUT_USER )
         Aux_Message( stderr, "WARNING : it's recommended to enable OPT__OUTPUT_USER !!\n" );
   }


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Validating test problem %d ... done\n", TESTPROB_ID );

} // FUNCTION : Validate



#if ( MODEL == ELBDM  &&  defined GRAVITY  &&  defined COMOVING )
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
   LOAD_PARA( load_mode, "Jeans_RealAmp0",    &Jeans_RealAmp0,       -1.0,           Eps_double,       NoMax_double      );
   LOAD_PARA( load_mode, "Jeans_Phase0",      &Jeans_Phase0,          0.0,           NoMin_double,     NoMax_double      );

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
   Jeans_ImagAmp0   = Jeans_RealAmp0;
   Jeans_Wavelength = amr->BoxSize[0]/sqrt(3.0);   // assuming cubic simulation domain
   Jeans_WaveK      = 2.0*M_PI/Jeans_Wavelength;
   Jeans_WaveKj     = POW( 6.0*Time[0]*SQR(ELBDM_ETA), 0.25 );
   Jeans_Stable     = ( Jeans_WaveK > Jeans_WaveKj ) ? true : false;


// (3) reset other general-purpose parameters
//     --> a helper macro PRINT_RESET_PARA is defined in Macro.h
   const long   End_Step_Default = __INT_MAX__;
// End_T : (stable/unstable) --> (1 period in the high-k limit / grow by a factor of 50 in the low-k limit)
   const double End_T_Default    = ( Jeans_Stable) ?
                                    A_INIT + pow( 0.5*Jeans_WaveK*Jeans_WaveK/M_PI/ELBDM_ETA, 2.0 ) :
                                    A_INIT*50;
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
      const double y = SQR(Jeans_WaveK)/ELBDM_ETA*pow( Time[0], -0.5 );

      Aux_Message( stdout, "=============================================================================\n" );
      Aux_Message( stdout, "  test problem ID       = %d\n",     TESTPROB_ID                      );
      Aux_Message( stdout, "  real part coefficient = %13.7e\n", Jeans_RealAmp0                   );
      Aux_Message( stdout, "  imag part coefficient = %13.7e\n", Jeans_ImagAmp0                   );
      Aux_Message( stdout, "  initial phase shift   = %13.7e\n", Jeans_Phase0                     );
      Aux_Message( stdout, "  real part amplitude   = %13.7e\n", Jeans_RealAmp(Jeans_RealAmp0, y) );
      Aux_Message( stdout, "  imag part amplitude   = %13.7e\n", Jeans_ImagAmp(Jeans_ImagAmp0, y) );
      Aux_Message( stdout, "  wavelength            = %13.7e\n", Jeans_Wavelength                 );
      Aux_Message( stdout, "  wavenumber            = %13.7e\n", Jeans_WaveK                      );
      Aux_Message( stdout, "  critical wavenumber   = %13.7e\n", Jeans_WaveKj                     );
      Aux_Message( stdout, "  Jeans stable          = %s\n",     (Jeans_Stable)?"YES":"NO"        );
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

   const double r       = 1.0/sqrt(3.0)*( x + y + z );
   const double Jeans_y = SQR(Jeans_WaveK)/ELBDM_ETA*pow( Time, -0.5 );
   const double Phase   = Jeans_WaveK*r + Jeans_Phase0;
   const double Re = 1.0 + Jeans_RealAmp( Jeans_RealAmp0, Jeans_y )*cos( Phase );
   const double Im =       Jeans_ImagAmp( Jeans_ImagAmp0, Jeans_y )*cos( Phase );

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

   const char Prefix[100]     = "JeansInstabilityComoving";
   const OptOutputPart_t Part = OUTPUT_DIAG;

   Output_L1Error( SetGridIC, NULL, Prefix, Part, NULL_REAL, NULL_REAL, NULL_REAL );

} // FUNCTION : OutputError
#endif // #if ( MODEL == ELBDM  &&  defined GRAVITY  &&  defined COMOVING )



//-------------------------------------------------------------------------------------------------------
// Function    :  Init_TestProb_ELBDM_JeansInstabilityComoving
// Description :  Test problem initializer
//
// Note        :  None
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Init_TestProb_ELBDM_JeansInstabilityComoving()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );


// validate the compilation flags and runtime parameters
   Validate();


#  if ( MODEL == ELBDM  &&  defined GRAVITY  &&  defined COMOVING )
// set the problem-specific runtime parameters
   SetParameter();


   Init_Function_User_Ptr    = SetGridIC;
   Output_User_Ptr           = OutputError;
#  ifdef SUPPORT_HDF5
   Output_HDF5_InputTest_Ptr = LoadInputTestProb;
#  endif
#  endif // #if ( MODEL == ELBDM  &&  defined GRAVITY  &&  defined COMOVING )


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : Init_TestProb_ELBDM_JeansInstabilityComoving
