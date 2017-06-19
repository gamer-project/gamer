#include "GAMER.h"

#if ( MODEL == HYDRO )



extern void (*Init_Function_Ptr)( real fluid[], const double x, const double y, const double z, const double Time );
extern void (*Output_TestProbErr_Ptr)( const bool BaseOnly );

static void HYDRO_TestProbSol_Caustic( real fluid[], const double x, const double y, const double z,
                                                         const double Time );
static void LoadTestProbParameter();


// global variables in the HYDRO caustic test
// =======================================================================================
static double Caustic_VelPeak;      // peak velocity
static double Caustic_Dens;         // background density
static double Caustic_Pres;         // background pressure
static int    Caustic_Dir;          // spatial direction: (0/1) --> (x/diagonal)
// =======================================================================================




//-------------------------------------------------------------------------------------------------------
// Function    :  Init_TestProb
// Description :  Initialize parameters for the HYDRO caustic test
//
// Note        :  1. Please copy this file to "GAMER/src/Init/Init_TestProb.cpp"
//                2. Global variables declared here will also be used in the function
//                   "HYDRO_TestProbSol_Caustic"
//
// Parameter   :  None
//-------------------------------------------------------------------------------------------------------
void Init_TestProb()
{

   const char *TestProb = "HYDRO caustic";

// check
#  if ( MODEL != HYDRO )
#  error : ERROR : "MODEL != HYDRO" in the HYDRO caustic test !!
#  endif

#  ifdef GRAVITY
#  error : ERROR : "GRAVITY must be OFF" in the HYDRO caustic test !!
#  endif

#  ifdef COMOVING
#  error : ERROR : "COMOVING must be OFF" in the HYDRO caustic test !!
#  endif

#  ifdef PARTICLE
#  error : ERROR : "PARTICLE must be OFF" in the HYDRO caustic test !!
#  endif

   if ( OPT__BC_FLU[0] != BC_FLU_PERIODIC )
      Aux_Error( ERROR_INFO, "Please set \"OPT__BC_FLU = 1\" for the %s test!!\n", TestProb );


// set the initialization and output functions
   Init_Function_Ptr      = HYDRO_TestProbSol_Caustic;
   Output_TestProbErr_Ptr = NULL;


// load the test problem parameters
   LoadTestProbParameter();


// record the test problem parameters
   if ( MPI_Rank == 0 )
   {
      Aux_Message( stdout, "\n" );
      Aux_Message( stdout, "%s test :\n", TestProb );
      Aux_Message( stdout, "=============================================================================\n" );
      Aux_Message( stdout, "NOTE : peak velocity         = % 14.7e\n", Caustic_VelPeak );
      Aux_Message( stdout, "       background density    = % 14.7e\n", Caustic_Dens    );
      Aux_Message( stdout, "       background pressure   = % 14.7e\n", Caustic_Pres    );
      Aux_Message( stdout, "       direction             = %s\n",      ( Caustic_Dir == 0 ) ? "X" :
                                                                       ( Caustic_Dir == 1 ) ? "Diagonal" :
                                                                                              "Unknown" );
      Aux_Message( stdout, "=============================================================================\n" );
      Aux_Message( stdout, "\n" );
   }


// set some default parameters
   const double End_T_Default    = ( Caustic_Dir == 0 ) ? 3.0 : 2.0;
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
// Function    :  HYDRO_TestProbSol_Caustic
// Description :  Calculate the initial condition in the HYDRO caustic test
//
// Note        :  1. This function is invoked by "HYDRO_Init_StartOver_AssignData"
//
// Parameter   :  fluid : Array to store the analytical solution to be returned
//                x/y/z : Target physical coordinates
//                Time  : Target physical time
//
// Return      :  fluid
//-------------------------------------------------------------------------------------------------------
void HYDRO_TestProbSol_Caustic( real fluid[], const double x, const double y, const double z, const double Time )
{

// x direction
   if ( Caustic_Dir == 0 )
   {
      const double WaveK = 2.0*M_PI/amr->BoxSize[0];

      fluid[MOMX] = Caustic_Dens*Caustic_VelPeak*sin( WaveK*x );
      fluid[MOMY] = 0.0;
      fluid[MOMZ] = 0.0;
   }

// diagonal direction
   else if ( Caustic_Dir == 1 )
   {
      const double WaveK = 2.0*M_PI/amr->BoxSize[0]*sqrt(3.0);
      const double r     = 1.0/sqrt(3.0)*( x + y + z );

      fluid[MOMX] = Caustic_Dens*Caustic_VelPeak*sin( WaveK*r ) / sqrt(3.0);
      fluid[MOMY] = fluid[MOMX];
      fluid[MOMZ] = fluid[MOMX];
   }

   else
      Aux_Error( ERROR_INFO, "Caustic_Dir = %d is NOT supported [0/1] !!\n", Caustic_Dir );

   fluid[DENS] = Caustic_Dens;
   fluid[ENGY] = Caustic_Pres/(GAMMA-1.0) + 0.5*( SQR(fluid[MOMX]) + SQR(fluid[MOMY]) + SQR(fluid[MOMZ]) ) / fluid[DENS];

} // FUNCTION : HYDRO_TestProbSol_Caustic



//-------------------------------------------------------------------------------------------------------
// Function    :  LoadTestProbParameter
// Description :  Load parameters for the test problem
//
// Note        :  This function is invoked by "Init_TestProb"
//
// Parameter   :  None
//-------------------------------------------------------------------------------------------------------
void LoadTestProbParameter()
{

   const char FileName[] = "Input__TestProb";

   FILE *File = fopen( FileName, "r" );

   if ( File == NULL )  Aux_Error( ERROR_INFO, "the file \"%s\" does not exist !!\n", FileName );

   int    tmp_int;
   char  *input_line = NULL;
   char   string[100];
   size_t len = 0;

   getline( &input_line, &len, File );

   getline( &input_line, &len, File );
   sscanf( input_line, "%lf%s",  &Caustic_VelPeak,     string );

   getline( &input_line, &len, File );
   sscanf( input_line, "%lf%s",  &Caustic_Dens,        string );

   getline( &input_line, &len, File );
   sscanf( input_line, "%lf%s",  &Caustic_Pres,        string );

   getline( &input_line, &len, File );
   sscanf( input_line, "%d%s",   &Caustic_Dir,         string );

   fclose( File );
   if ( input_line != NULL )     free( input_line );


// check
   if ( Caustic_VelPeak < 0.0 )  Aux_Error( ERROR_INFO, "%s = %14.7e < 0.0 !!\n", "Caustic_VelPeak", Caustic_VelPeak );
   if ( Caustic_Dens    < 0.0 )  Aux_Error( ERROR_INFO, "%s = %14.7e < 0.0 !!\n", "Caustic_Dens",    Caustic_Dens );
   if ( Caustic_Pres    < 0.0 )  Aux_Error( ERROR_INFO, "%s = %14.7e < 0.0 !!\n", "Caustic_Pres",    Caustic_Pres );

   if ( Caustic_Dir < 0  ||  Caustic_Dir > 1 )
      Aux_Error( ERROR_INFO, "%s = %d is NOT supported [0/1] !!\n", "Caustic_Dir", Caustic_Dir );

   if (  Caustic_Dir == 1  &&  ( amr->BoxSize[0] != amr->BoxSize[1] || amr->BoxSize[0] != amr->BoxSize[2] )  )
         Aux_Error( ERROR_INFO, "simulation domain must be CUBIC for %s = %d !!\n", "Caustic_Dir", Caustic_Dir );

} // FUNCTION : LoadTestProbParameter



#endif // #if ( MODEL == HYDRO )
