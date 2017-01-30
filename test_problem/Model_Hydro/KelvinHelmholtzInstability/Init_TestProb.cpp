#include "GAMER.h"

#if ( MODEL == HYDRO )



extern void (*Init_Function_Ptr)( real fluid[], const double x, const double y, const double z, const double Time );
extern void (*Output_TestProbErr_Ptr)( const bool BaseOnly );

       void HYDRO_TestProbSol_KH( real fluid[], const double x, const double y, const double z, const double Time );
static void LoadTestProbParameter();
static double RandomNumber( struct drand48_data *Buf, const double Min, const double Max );


// global variables in the HYDRO Kelvin-Helmholtz instability test
// =======================================================================================
int    KH_RSeed;        // random seed
double KH_RAmp;         // random number amplitude in both vx, vy, and vz
double KH_Pres;         // background pressure
double KH_Vx1;          // velocity x in the upper region
double KH_Vx2;          // velocity x in the lower region
double KH_Vy1;          // velocity y in the upper region
double KH_Vy2;          // velocity y in the lower region
double KH_Rho1;         // density in the upper region
double KH_Rho2;         // density in the lower region
bool   KH_AllRankSame;  // all MPI ranks assign the same initial condition --> suitable for the weak scaling test

static struct drand48_data drand_buf;
// =======================================================================================




//-------------------------------------------------------------------------------------------------------
// Function    :  Init_TestProb
// Description :  Initialize parameters for the HYDRO Kelvin-Helmholtz instability test
//
// Note        :  1. Please copy this file to "GAMER/src/Init/Init_TestProb.cpp"
//                2. Global variables declared here will also be used in the function
//                   "HYDRO_TestProbSol_KH"
//
// Parameter   :  None
//-------------------------------------------------------------------------------------------------------
void Init_TestProb()
{

   const char *TestProb = "HYDRO Kelvin-Helmholtz instability";

// check
#  if ( MODEL != HYDRO )
#  error : ERROR : "MODEL != HYDRO" in the HYDRO Kelvin-Helmholtz instability test !!
#  endif

#  ifdef GRAVITY
#  error : ERROR : "GRAVITY must be OFF" in the HYDRO Kelvin-Helmholtz instability test !!
#  endif

#  ifdef COMOVING
#  error : ERROR : "COMOVING must be OFF" in the HYDRO Kelvin-Helmholtz instability test !!
#  endif

#  ifdef PARTICLE
#  error : ERROR : "PARTICLE must be OFF" in the HYDRO Kelvin-Helmholtz instability test !!
#  endif


// set the initialization and output functions
   Init_Function_Ptr      = HYDRO_TestProbSol_KH;
   Output_TestProbErr_Ptr = NULL;


// load the test problem parameters
   LoadTestProbParameter();


// set the global variables
   if ( KH_AllRankSame )   srand48_r( KH_RSeed,          &drand_buf );
   else                    srand48_r( KH_RSeed+MPI_Rank, &drand_buf );


// record the test problem parameters
   if ( MPI_Rank == 0 )
   {
      Aux_Message( stdout, "\n" );
      Aux_Message( stdout, "%s test :\n", TestProb );
      Aux_Message( stdout, "=============================================================================\n" );
      Aux_Message( stdout, "  KH_RSeed               = %d\n",      KH_RSeed );
      Aux_Message( stdout, "  KH_RAmp                = % 14.7e\n", KH_RAmp );
      Aux_Message( stdout, "  KH_Pres                = % 14.7e\n", KH_Pres );
      Aux_Message( stdout, "  KH_Vx1                 = % 14.7e\n", KH_Vx1 );
      Aux_Message( stdout, "  KH_Vx2                 = % 14.7e\n", KH_Vx2 );
      Aux_Message( stdout, "  KH_Vy1                 = % 14.7e\n", KH_Vy1 );
      Aux_Message( stdout, "  KH_Vy2                 = % 14.7e\n", KH_Vy2 );
      Aux_Message( stdout, "  KH_Rho1                = % 14.7e\n", KH_Rho1 );
      Aux_Message( stdout, "  KH_Rho2                = % 14.7e\n", KH_Rho2 );
      Aux_Message( stdout, "  KH_AllRankSame         = %d\n",      KH_AllRankSame );
      Aux_Message( stdout, "=============================================================================\n" );
      Aux_Message( stdout, "\n" );
   } // if ( MPI_Rank == 0 )


// set some default parameters
   const double End_T_Default    = 1.0;
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

   if ( OPT__BC_FLU[0] != BC_FLU_PERIODIC )
      Aux_Error( ERROR_INFO, "Please set \"OPT__BC_FLU = 1\" for the %s test!!\n", TestProb );

} // FUNCTION : Init_TestProb



//-------------------------------------------------------------------------------------------------------
// Function    :  HYDRO_TestProbSol_KH
// Description :  Generate initial condition for the HYDRO Kelvin-Helmholtz instability test
//
// Note        :  1. This function is invoked by "HYDRO_Init_StartOver_AssignData"
//
// Parameter   :  fluid : Array to store the analytical solution to be returned
//                x/y/z : Target physical coordinates
//                Time  : Target physical time
//
// Return      :  fluid
//-------------------------------------------------------------------------------------------------------
void HYDRO_TestProbSol_KH( real fluid[], const double x, const double y, const double z, const double Time )
{

   const double dz_periodic = ( KH_AllRankSame ) ? amr->BoxSize[2] / MPI_NRank_X[2] : amr->BoxSize[2];
   const double z_shear     = 0.5*dz_periodic;
   const double z_periodic  = fmod( z, dz_periodic );

// region 1
   if ( z_periodic >= z_shear )
   {
      fluid[DENS] = KH_Rho1;
      fluid[MOMX] = KH_Rho1*( KH_Vx1 + RandomNumber(&drand_buf,-KH_RAmp,KH_RAmp) );
      fluid[MOMY] = KH_Rho1*( KH_Vy1 + RandomNumber(&drand_buf,-KH_RAmp,KH_RAmp) );
      fluid[MOMZ] = KH_Rho1*(    0.0 + RandomNumber(&drand_buf,-KH_RAmp,KH_RAmp) );
   }

// region 2
   else
   {
      fluid[DENS] = KH_Rho2;
      fluid[MOMX] = KH_Rho2*( KH_Vx2 + RandomNumber(&drand_buf,-KH_RAmp,KH_RAmp) );
      fluid[MOMY] = KH_Rho2*( KH_Vy2 + RandomNumber(&drand_buf,-KH_RAmp,KH_RAmp) );
      fluid[MOMZ] = KH_Rho2*(    0.0 + RandomNumber(&drand_buf,-KH_RAmp,KH_RAmp) );
   }

   fluid[ENGY] = KH_Pres/(GAMMA-1.0) + 0.5*( SQR(fluid[MOMX]) + SQR(fluid[MOMY]) + SQR(fluid[MOMZ]) ) / fluid[DENS];

} // FUNCTION : HYDRO_TestProbSol_KH



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
   sscanf( input_line, "%d%s",   &KH_RSeed,            string );

   getline( &input_line, &len, File );
   sscanf( input_line, "%lf%s",  &KH_RAmp,             string );

   getline( &input_line, &len, File );
   sscanf( input_line, "%lf%s",  &KH_Pres,             string );

   getline( &input_line, &len, File );
   sscanf( input_line, "%lf%s",  &KH_Vx1,              string );

   getline( &input_line, &len, File );
   sscanf( input_line, "%lf%s",  &KH_Vx2,              string );

   getline( &input_line, &len, File );
   sscanf( input_line, "%lf%s",  &KH_Vy1,              string );

   getline( &input_line, &len, File );
   sscanf( input_line, "%lf%s",  &KH_Vy2,              string );

   getline( &input_line, &len, File );
   sscanf( input_line, "%lf%s",  &KH_Rho1,             string );

   getline( &input_line, &len, File );
   sscanf( input_line, "%lf%s",  &KH_Rho2,             string );

   getline( &input_line, &len, File );
   sscanf( input_line, "%d%s",  &tmp_int,              string );
   KH_AllRankSame = (bool)tmp_int;


   fclose( File );
   if ( input_line != NULL )     free( input_line );

} // FUNCTION : LoadTestProbParameter



//-------------------------------------------------------------------------------------------------------
// Function    :  RandomNumber
// Description :  Generate a single random number
//
// Note        :  1. Use drand48_r() instead of rand(). The latter (i) is not thread-safe and (ii) doesn't
//                   seem to work well with MPI even with a single thread (the generated random numbers are
//                   found to be irreproducible with MPI + AMR)
//                2. Must call srand48_r() in advance to set the random seed
//
// Parameter   :  Buf : Buffer used by drand48_r() for generating a random number
//                Min : Lower limit of the random number
//                Max : Upper limit of the random number
//
// Return      :  Random number
//-------------------------------------------------------------------------------------------------------
double RandomNumber( struct drand48_data *Buf, const double Min, const double Max )
{

   double RVal;

   drand48_r( Buf, &RVal );

// drand48_r returns [0.0, 1.0)
   return RVal*(Max-Min) + Min;

} // FUNCTION : RandomNumber



#endif // #if ( MODEL == HYDRO )
