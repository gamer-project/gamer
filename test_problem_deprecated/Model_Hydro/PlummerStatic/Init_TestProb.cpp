#include "GAMER.h"

#if ( MODEL == HYDRO )



extern void (*Init_Function_Ptr)( real fluid[], const double x, const double y, const double z, const double Time );
extern void (*Output_TestProbErr_Ptr)( const bool BaseOnly );

       void HYDRO_TestProbSol_PlummerStatic( real fluid[], const double x, const double y, const double z, const double Time );
static void LoadTestProbParameter();


// global variables in the HYDRO Plummer hydrostatic test
// =======================================================================================
static real Plummer_r0;       // scale radius
static real Plummer_M;        // total mass
static real Plummer_rho0;     // peak density
static real Plummer_FreeT;    // free-fall time at Plummer_r0
// =======================================================================================




//-------------------------------------------------------------------------------------------------------
// Function    :  Init_TestProb
// Description :  Initialize parameters for the HYDRO Plummer hydrostatic test
//
// Note        :  1. Please copy this file to "GAMER/src/Init/Init_TestProb.cpp"
//                2. Global variables declared here will also be used in the function
//                   "HYDRO_TestProbSol_PlummerStatic"
//
// Parameter   :  None 
//-------------------------------------------------------------------------------------------------------
void Init_TestProb()
{  

   const char *TestProb = "HYDRO Plummer hydrostatic";

// check
#  if ( MODEL != HYDRO )
#  error : ERROR : "MODEL != HYDRO" in the HYDRO Plummer hydrostatic test !!
#  endif

#  ifndef GRAVITY
#  error : ERROR : "GRAVITY must be ON" in the HYDRO Plummer hydrostatic test !!
#  endif

#  ifdef COMOVING
#  error : ERROR : "COMOVING must be OFF" in the HYDRO Plummer hydrostatic test !!
#  endif

   if ( amr->BoxSize[0] != amr->BoxSize[1]  ||  amr->BoxSize[0] != amr->BoxSize[2] )
      Aux_Error( ERROR_INFO, "simulation domain must be CUBIC in the %s test !!\n", TestProb );


// set the initialization and output functions
   Init_Function_Ptr      = HYDRO_TestProbSol_PlummerStatic;
   Output_TestProbErr_Ptr = NULL;


// load the test problem parameters
   LoadTestProbParameter();


// set global variables
   Plummer_rho0  = 3.0*Plummer_M / ( 4.0*M_PI*CUBE(Plummer_r0) );
   Plummer_FreeT = sqrt( (3.0*M_PI*pow(2.0,1.5)) / (32.0*NEWTON_G*Plummer_rho0) );


// record the test problem parameters
   if ( MPI_Rank == 0 )
   {
      Aux_Message( stdout, "\n" );
      Aux_Message( stdout, "%s test :\n", TestProb );
      Aux_Message( stdout, "=============================================================================\n" );
      Aux_Message( stdout, "NOTE : scale radius           = %13.7e\n", Plummer_r0    );
      Aux_Message( stdout, "       total mass             = %13.7e\n", Plummer_M     );
      Aux_Message( stdout, "       peak density           = %13.7e\n", Plummer_rho0  );
      Aux_Message( stdout, "       free fall time         = %13.7e\n", Plummer_FreeT );
      Aux_Message( stdout, "=============================================================================\n" );
      Aux_Message( stdout, "\n" );
   }


// set some default parameters
// End_T : 5 free-fall time at the scale radius
   const double End_T_Default    = 5.0*Plummer_FreeT;
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

   if ( OPT__BC_POT != BC_POT_ISOLATED )
      Aux_Error( ERROR_INFO, "Please set \"OPT__BC_POT = 1\" for the %s test!!\n", TestProb );

} // FUNCTION : Init_TestProb



//-------------------------------------------------------------------------------------------------------
// Function    :  HYDRO_TestProbSol_PlummerStatic
// Description :  Calculate the analytical solution in the HYDRO Plummer hydrostatic test  
//
// Note        :  1. Wave vector is along the diagonal direction
//                2. Background density is assumed to be ONE
//                3. This function is invoked by "HYDRO_Init_StartOver_AssignData" and "Output_TestProbErr" 
//
// Parameter   :  fluid : Array to store the analytical solution to be returned
//                x/y/z : Target physical coordinates 
//                Time  : Target physical time
//
// Return      :  fluid
//-------------------------------------------------------------------------------------------------------
void HYDRO_TestProbSol_PlummerStatic( real fluid[], const double x, const double y, const double z, const double Time )
{

   const double Cen[3]    = { 0.5*amr->BoxSize[0],
                            0.5*amr->BoxSize[1],
                            0.5*amr->BoxSize[2] };
   const double r         = sqrt( SQR(x-Cen[0]) + SQR(y-Cen[1]) + SQR(z-Cen[2]) );
   const double a         = r / Plummer_r0;
   const double a2        = a*a;
   const double _Gamma_m1 = 1.0/(GAMMA-1.0);

   fluid[DENS] = Plummer_rho0 * pow( 1.0 + a2, -2.5 );
   fluid[MOMX] = 0.0;
   fluid[MOMY] = 0.0;
   fluid[MOMZ] = 0.0;
   fluid[ENGY] = (  NEWTON_G*Plummer_rho0*Plummer_M / ( 6.0*Plummer_r0*CUBE(1.0+a2) )  )*_Gamma_m1;

} // FUNCTION : HYDRO_TestProbSol_PlummerStatic



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
   sscanf( input_line, "%lf%s",  &Plummer_r0,               string );

   getline( &input_line, &len, File );
   sscanf( input_line, "%lf%s",  &Plummer_M,                string );

#  else

   getline( &input_line, &len, File );
   sscanf( input_line, "%f%s",  &Plummer_r0,                string );

   getline( &input_line, &len, File );
   sscanf( input_line, "%f%s",  &Plummer_M,                 string );
#  endif // #ifdef FLOAT8 ... else ...

   fclose( File );
   if ( input_line != NULL )     free( input_line );

} // FUNCTION : LoadTestProbParameter



#endif // #if ( MODEL == HYDRO )
