#include "GAMER.h"

#ifdef PARTICLE



extern void (*Init_Function_Ptr)( real fluid[], const double x, const double y, const double z, const double Time );
extern void (*Output_TestProbErr_Ptr)( const bool BaseOnly );

static void LoadTestProbParameter();
static void Par_TestProbSol_TwoParOrbit( real fluid[], const double x, const double y, const double z, const double Time );


// global variables in the two particles orbit test
// =======================================================================================
real TwoParOrbit_M;          // particle mass
real TwoParOrbit_v;             // particle velocity
real TwoParOrbit_R;             // radius of the orbit
real TwoParOrbit_T;             // period of the orbit
// =======================================================================================




//-------------------------------------------------------------------------------------------------------
// Function    :  Init_TestProb
// Description :  Initialize parameters for the two particles orbit test
//
// Note        :  1. Please copy this file to "GAMER/src/Init/Init_TestProb.cpp"
//                2. Test problem parameters can be set in the input file "Input__TestProb"
//
// Parameter   :  None
//-------------------------------------------------------------------------------------------------------
void Init_TestProb()
{

   const char *TestProb = "Two particles orbit";

// check
# ifndef PARTICLE
# error : ERROR : "PARTICLE is NOT defined" in the two particles orbit test !!
# endif

# ifndef GRAVITY
# error : ERROR : "GRAVITY must be ON" in the two particles orbit test !!
# endif

# ifdef COMOVING
# error : ERROR : "COMOVING must be OFF" in the two particles orbit test !!
# endif

   if ( OPT__BC_POT != BC_POT_ISOLATED )
      Aux_Error( ERROR_INFO, "pleaset set parameter %s to %d in the %s test !!\n",
                 "BC_POT_ISOLATED", BC_POT_ISOLATED, TestProb );

   if ( amr->Par->NPar_Active_AllRank != 2 )
      Aux_Error( ERROR_INFO, "please set parameter %s to %d in the %s test !!\n",
                 "amr->Par->NPar_Active_AllRank", amr->Par->NPar_Active_AllRank, TestProb );

   if ( amr->Par->Init != PAR_INIT_BY_FUNCTION )
      Aux_Error( ERROR_INFO, "please set parameter %s to %d in the %s test !!\n",
                 "amr->Par->Init", amr->Par->Init, TestProb );

   if ( OPT__INIT != INIT_STARTOVER )
      Aux_Error( ERROR_INFO, "please set parameter %s to %d in the %s test !!\n",
                 "OPT__INIT", OPT__INIT, TestProb );


// set the initialization and output functions
   Init_Function_Ptr      = Par_TestProbSol_TwoParOrbit;
   Output_TestProbErr_Ptr = NULL;


// load the test problem parameters
   LoadTestProbParameter();


// set the test problem parameters
   TwoParOrbit_v = SQRT( (NEWTON_G*TwoParOrbit_M) / (4.0*TwoParOrbit_R) );
   TwoParOrbit_T = 2.0*M_PI*TwoParOrbit_R / TwoParOrbit_v;


// record the test problem parameters
   if ( MPI_Rank == 0 )
   {
      Aux_Message( stdout, "\n" );
      Aux_Message( stdout, "%s test :\n", TestProb );
      Aux_Message( stdout, "=============================================================================\n" );
      Aux_Message( stdout, " Note: particle mass                                = %13.7e\n", TwoParOrbit_M );
      Aux_Message( stdout, "       particle velocity                            = %13.7e\n", TwoParOrbit_v );
      Aux_Message( stdout, "       radius of the orbit                          = %13.7e\n", TwoParOrbit_R );
      Aux_Message( stdout, "       period of the orbit                          = %13.7e\n", TwoParOrbit_T );
      Aux_Message( stdout, "=============================================================================\n" );
      Aux_Message( stdout, "\n" );
   }


// set some default parameters
// End_T : 1 period
   const double End_T_Default    = 1.0*TwoParOrbit_T;
   const long   End_Step_Default = __INT_MAX__;

   if ( END_STEP < 0 )
   {
      END_STEP = End_Step_Default;

      if ( MPI_Rank == 0 )
         Aux_Message( stdout, "NOTE : parameter %s is set to %ld in the %s test !!\n", "END_STEP", END_STEP, TestProb );
   }

   if ( END_T < 0.0 )
   {
      END_T = End_T_Default;

      if ( MPI_Rank == 0 )
         Aux_Message( stdout, "NOTE : parameter %s is set to %13.7e in the %s test !!\n", "END_T", END_T, TestProb );
   }

} // FUNCTION : Init_TestProb



//-------------------------------------------------------------------------------------------------------
// Function    :  Par_TestProbSol_TwoParOrbit
// Description :  Initialize the background density field as zero for the two particles orbit test
//
// Note        :  1. This test works for both ELBDM and HYDRO models
//                   ELBDM : test external potential
//                   HYDRO : test external acceleration
//                2. Invoked by "ELBDM/Hydro_Init_StartOver_AssignData"
//
// Parameter   :  fluid : Fluid field to be initialized
//                x/y/z : Target physical coordinates
//                Time  : Target physical time
//
// Return      :  fluid
//-------------------------------------------------------------------------------------------------------
void Par_TestProbSol_TwoParOrbit( real *fluid, const double x, const double y, const double z, const double Time )
{

#  if ( MODEL == HYDRO )
// set density to negligibly small
   fluid[DENS] = 1.0e-20;
   fluid[MOMX] = 0.0;
   fluid[MOMY] = 0.0;
   fluid[MOMZ] = 0.0;
   fluid[ENGY] = 1.0e-22;

#  elif ( MODEL == ELBDM )
// set wave function as zero everywhere
   fluid[REAL] = 0.0;
   fluid[IMAG] = 0.0;
   fluid[DENS] = 0.0;

#  else
#  error : ERROR : unsupported model !!
#  endif

} // FUNCTION : Par_TestProbSol_TwoParOrbit



//-------------------------------------------------------------------------------------------------------
// Function    :  LoadTestProbParameter
// Description :  Load parameters for the test problem
//
// Note        :  This function is invoked by "Init_TestProb"
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void LoadTestProbParameter()
{

   const char FileName[] = "Input__TestProb";

   FILE *File = fopen( FileName, "r" );

   if ( File == NULL )  Aux_Error( ERROR_INFO, "the file \"%s\" does not exist !!\n", FileName );

   char  *input_line = NULL;
   char   string[100];
   size_t len = 0;
   int    temp_int;


   getline( &input_line, &len, File );

#  ifdef FLOAT8
   getline( &input_line, &len, File );
   sscanf( input_line, "%lf%s",  &TwoParOrbit_M,              string );

   getline( &input_line, &len, File );
   sscanf( input_line, "%lf%s",  &TwoParOrbit_R,              string );
#  else
   getline( &input_line, &len, File );
   sscanf( input_line, "%f%s",   &TwoParOrbit_M,              string );

   getline( &input_line, &len, File );
   sscanf( input_line, "%f%s",   &TwoParOrbit_R,              string );
#  endif // #ifdef FLOAT8 ... else ...

   fclose( File );
   if ( input_line != NULL )     free( input_line );


// check
   if ( TwoParOrbit_M < 0.0 )    Aux_Error( ERROR_INFO, "TwoParOrbit_M (%14.7e) < 0.0 !!\n", TwoParOrbit_M );
   if ( TwoParOrbit_R < 0.0 )    Aux_Error( ERROR_INFO, "TwoParOrbit_R (%14.7e) < 0.0 !!\n", TwoParOrbit_R );

} // FUNCTION : LoadTestProbParameter



#endif // #ifdef PARTICLE
