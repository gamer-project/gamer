#include "GAMER.h"

#if ( MODEL == HYDRO )



extern void (*Init_Function_Ptr)( real fluid[], const real x, const real y, const real z, const double Time );
extern void (*Output_TestProbErr_Ptr)( const bool BaseOnly );

static void LoadTestProbParameter();
static void Hydro_TestProbSol_SphericalCollapse( real fluid[], const real x, const real y, const real z, const double Time );


// global variables in the HYDRO spherical collapse test
// =======================================================================================
real SCollapse_Dens_Bg;       // background mass density
real SCollapse_Dens_Delta;    // top-hat mass density --> total density = Dens_Bg*( 1 + Dens_Delta )
real SCollapse_Engy_Bg;       // background energy density
real SCollapse_Radius;        // top-hat radius
real SCollapse_Center[3];     // top-hat center
// =======================================================================================




//-------------------------------------------------------------------------------------------------------
// Function    :  Init_TestProb
// Description :  Initialize parameters for the HYDRO spherical collapse test
//
// Note        :  1. Please copy this file to "GAMER/src/Init/Init_TestProb.cpp"
//                2. Test problem parameters can be set in the input file "Input__TestProb"
//
// Parameter   :  None 
//-------------------------------------------------------------------------------------------------------
void Init_TestProb()
{  

   const char *TestProb = "HYDRO spherical collapse";

// check
#  if ( MODEL != HYDRO )
#  error : ERROR : "MODEL != HYDRO" in the HYDRO spherical collapse test !!
#  endif

#  ifndef GRAVITY
#  error : ERROR : "GRAVITY must be ON" in the HYDRO spherical collapse test !!
#  endif

#  ifndef COMOVING
#  error : ERROR : "COMOVING must be ON" in the HYDRO spherical collapse test !!
#  endif


// set the initialization and output functions
   Init_Function_Ptr      = Hydro_TestProbSol_SphericalCollapse;
   Output_TestProbErr_Ptr = NULL;


// load the test problem parameters
   LoadTestProbParameter();


// record the test problem parameters
   if ( MPI_Rank == 0 )
   {
      Aux_Message( stdout, "\n" );
      Aux_Message( stdout, "%s test :\n", TestProb );
      Aux_Message( stdout, "=============================================================================\n" );
      Aux_Message( stdout, "NOTE : background mass density   = %13.7e\n", SCollapse_Dens_Bg );
      Aux_Message( stdout, "     : top-hat over density      = %13.7e\n", SCollapse_Dens_Delta );
      Aux_Message( stdout, "       background energy density = %13.7e\n", SCollapse_Engy_Bg);
      Aux_Message( stdout, "       top-hat radius            = %13.7e\n", SCollapse_Radius );
      Aux_Message( stdout, "       top-hat center            = (%13.7e, %13.7e, %13.7e)\n", SCollapse_Center[0], 
                                                                                            SCollapse_Center[1],
                                                                                            SCollapse_Center[2] );
      Aux_Message( stdout, "=============================================================================\n" );
      Aux_Message( stdout, "\n" );
   }


// set some default parameters
   const double End_T_Default    = 5.0e-2;
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
         Aux_Message( stdout, "NOTE : parameter %s is reset to %13.7e in the %s test !!\n",
                      "END_T", END_T, TestProb );
   }

   if ( !OPT__INIT_RESTRICT )
   {
      OPT__INIT_RESTRICT = true;

      if ( MPI_Rank == 0 )
         Aux_Message( stdout, "NOTE : parameter %s is reset to %d in the %s test !!\n",
                      "OPT__INIT_RESTRICT", OPT__INIT_RESTRICT, TestProb );
   }

   if ( OPT__OUTPUT_TEST_ERROR  &&  Output_TestProbErr_Ptr == NULL )
   {
      OPT__OUTPUT_TEST_ERROR = false;

      if ( MPI_Rank == 0 )
         Aux_Message( stdout, "NOTE : option %s is not supported for the %s test and has been disabled !!\n",
                      "OPT__OUTPUT_TEST_ERROR", TestProb );
   }
   else if ( !OPT__OUTPUT_TEST_ERROR  &&  Output_TestProbErr_Ptr != NULL )
   {
      OPT__OUTPUT_TEST_ERROR = true;

      if ( MPI_Rank == 0 )
         Aux_Message( stdout, "NOTE : option %s is enabled for the %s test !!\n", 
                      "OPT__OUTPUT_TEST_ERROR", TestProb );
   }

} // FUNCTION : Init_TestProb



//-------------------------------------------------------------------------------------------------------
// Function    :  Hydro_TestProbSol_SphericalCollapse
// Description :  Initialize the HYDRO spherical collapse test  
//
// Note        :  Invoked by "Hydro_Init_StartOver_AssignData"
//
// Parameter   :  fluid : Fluid field to be initialized
//                x/y/z : Target physical coordinates
//                Time  : Target physical time
//
// Return      :  fluid
//-------------------------------------------------------------------------------------------------------
void Hydro_TestProbSol_SphericalCollapse( real *fluid, const real x, const real y, const real z, const double Time )
{

   const real r = SQRT( SQR(x-SCollapse_Center[0]) + SQR(y-SCollapse_Center[1]) + SQR(z-SCollapse_Center[2]) );

   fluid[DENS] = SCollapse_Dens_Bg;
   fluid[MOMX] = 0.0;
   fluid[MOMY] = 0.0;
   fluid[MOMZ] = 0.0;
   fluid[ENGY] = SCollapse_Engy_Bg;

   if ( r <= SCollapse_Radius )   fluid[DENS] *= ( 1.0 + SCollapse_Dens_Delta );

} // FUNCTION : Hydro_TestProbSol_SphericalCollapse



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

   int    temp_int;
   char  *input_line = NULL;
   char   string[100];
   size_t len = 0;

   getline( &input_line, &len, File );

#  ifdef FLOAT8
   getline( &input_line, &len, File );
   sscanf( input_line, "%lf%s",  &SCollapse_Dens_Bg,        string );

   getline( &input_line, &len, File );
   sscanf( input_line, "%lf%s",  &SCollapse_Dens_Delta,     string );

   getline( &input_line, &len, File );
   sscanf( input_line, "%lf%s",  &SCollapse_Engy_Bg,        string );

   getline( &input_line, &len, File );
   sscanf( input_line, "%lf%s",  &SCollapse_Radius,         string );

   getline( &input_line, &len, File );
   sscanf( input_line, "%lf%s",  &SCollapse_Center[0],      string );

   getline( &input_line, &len, File );
   sscanf( input_line, "%lf%s",  &SCollapse_Center[1],      string );

   getline( &input_line, &len, File );
   sscanf( input_line, "%lf%s",  &SCollapse_Center[2],      string );

#  else

   getline( &input_line, &len, File );
   sscanf( input_line, "%f%s",   &SCollapse_Dens_Bg,        string );

   getline( &input_line, &len, File );
   sscanf( input_line, "%f%s",   &SCollapse_Dens_Delta,     string );

   getline( &input_line, &len, File );
   sscanf( input_line, "%f%s",   &SCollapse_Engy_Bg,        string );

   getline( &input_line, &len, File );
   sscanf( input_line, "%f%s",   &SCollapse_Radius,         string );

   getline( &input_line, &len, File );
   sscanf( input_line, "%f%s",   &SCollapse_Center[0],      string );

   getline( &input_line, &len, File );
   sscanf( input_line, "%f%s",   &SCollapse_Center[1],      string );

   getline( &input_line, &len, File );
   sscanf( input_line, "%f%s",   &SCollapse_Center[2],      string );
#  endif // #ifdef FLOAT8 ... else ...

   fclose( File );
   if ( input_line != NULL )     free( input_line );


// set the default explosion center
   for (int d=0; d<3; d++)    
      if ( SCollapse_Center[d] < 0.0 )    SCollapse_Center[d] = 0.5*amr->BoxSize[d];

} // FUNCTION : LoadTestProbParameter



#endif // #if ( MODEL == HYDRO )
