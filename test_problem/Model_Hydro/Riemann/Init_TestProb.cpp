#include "GAMER.h"

#if ( MODEL == HYDRO )



extern void (*Init_Function_Ptr)( real fluid[], const real x, const real y, const real z, const double Time );
extern void (*Output_TestProbErr_Ptr)( const bool BaseOnly );

static void LoadTestProbParameter();
static void Hydro_TestProbSol_Riemann( real fluid[], const real x, const real y, const real z, const double Time );


// global variables in the HYDRO Riemann problem test
// =======================================================================================
enum Riemann_t { SOD_SHOCK_TUBE=0, STRONG_SHOCK=1, TWO_SHOCKS=2, EINFELDT_1203=3, EINFELDT_1125=4, SONIC_RARE=5 };

static Riemann_t Riemann_Prob;   // target Riemann problem
static char Riemann_Name[100];   // name of the target Riemann problem
static real Riemann_RhoL;        // left-state density
static real Riemann_VelL;        // left-state velocity
static real Riemann_VelL_T;      // left-state transverse velocity
static real Riemann_PreL;        // left-state pressure
static real Riemann_RhoR;        // right-state density
static real Riemann_VelR;        // right-state velocity
static real Riemann_VelR_T;      // right-state transverse velocity
static real Riemann_PreR;        // right-state pressure
static real Riemann_EndT;        // end physical time
static int  Riemann_LR;          // wave propagation direction (>0/<0 --> positive/negative direction)
static int  Riemann_XYZ;         // wave propagation direction (0/1/2 --> x/y/z)
// =======================================================================================




//-------------------------------------------------------------------------------------------------------
// Function    :  Init_TestProb
// Description :  Initialize parameters for the HYDRO Riemann problem test
//
// Note        :  1. Please copy this file to "GAMER/src/Init/Init_TestProb.cpp"
//                2. Test problem parameters cab be set in the input file "Input__TestProb"
//             
// Parameter   :  None 
//-------------------------------------------------------------------------------------------------------
void Init_TestProb()
{  

   const char *TestProb = "HYDRO Riemann problem";

// check
#  if ( MODEL != HYDRO )
#  error : ERROR : "MODEL != HYDRO" in the HYDRO Riemann problem test !!
#  endif

#  ifdef GRAVITY
#  error : ERROR : "GRAVITY must be OFF" in the HYDRO Riemann problem test !!
#  endif

#  ifdef COMOVING
#  error : ERROR : "COMOVING must be OFF" in the HYDRO Riemann problem test !!
#  endif


// set the initialization and output functions
   Init_Function_Ptr      = Hydro_TestProbSol_Riemann;
   Output_TestProbErr_Ptr = NULL;


// load the test problem parameters
   LoadTestProbParameter();


// set the test problem parameters
   switch ( Riemann_Prob )
   {
      case SOD_SHOCK_TUBE : Riemann_RhoL = 1.0;    Riemann_VelL = 0.0;  Riemann_PreL = 1.0;  Riemann_VelL_T = 0.0;
                            Riemann_RhoR = 0.125;  Riemann_VelR = 0.0;  Riemann_PreR = 0.1;  Riemann_VelR_T = 0.0;
                            Riemann_EndT = 0.1;
                            sprintf( Riemann_Name, "Sod's shock tube" );
                            break;

      case STRONG_SHOCK   : Riemann_RhoL = 1250.0;  Riemann_VelL = 0.0;  Riemann_PreL = 500.0;  Riemann_VelL_T = 0.0;
                            Riemann_RhoR =  125.0;  Riemann_VelR = 0.0;  Riemann_PreR =   5.0;  Riemann_VelR_T = 0.0;
                            Riemann_EndT = 0.4;
                            sprintf( Riemann_Name, "strong shock" );
                            break;

      case TWO_SHOCKS     : Riemann_RhoL = 1.0;  Riemann_VelL = 3.0;  Riemann_PreL = 1.0;  Riemann_VelL_T = 0.0;
                            Riemann_RhoR = 2.0;  Riemann_VelR = 1.0;  Riemann_PreR = 1.0;  Riemann_VelR_T = 0.0;
                            Riemann_EndT = 0.1;
                            sprintf( Riemann_Name, "two shocks" );
                            break;

      case EINFELDT_1203  : Riemann_RhoL = 1.0;  Riemann_VelL = -2.0;  Riemann_PreL = GAMMA-1.0;  Riemann_VelL_T = 0.0;
                            Riemann_RhoR = 1.0;  Riemann_VelR = +2.0;  Riemann_PreR = GAMMA-1.0;  Riemann_VelR_T = 0.0;
                            Riemann_EndT = 0.1;
                            sprintf( Riemann_Name, "Einfeldt's 1-2-0-3" );
                            break;

      case EINFELDT_1125  : Riemann_RhoL = 1.0;  Riemann_VelL = -1.0;  Riemann_PreL = 2.5*(GAMMA-1.0);  Riemann_VelL_T = -2.0;
                            Riemann_RhoR = 1.0;  Riemann_VelR = +1.0;  Riemann_PreR = 2.5*(GAMMA-1.0);  Riemann_VelR_T = +2.0;
                            Riemann_EndT = 0.1;
                            sprintf( Riemann_Name, "Einfeldt's 1-1-2-5" );
                            break;

      case SONIC_RARE     : Riemann_RhoL = 1.0;    Riemann_VelL = 0.75;  Riemann_PreL = 1.0;  Riemann_VelL_T = 0.0;
                            Riemann_RhoR = 0.125;  Riemann_VelR = 0.0;   Riemann_PreR = 0.1;  Riemann_VelR_T = 0.0;
                            Riemann_EndT = 0.1;
                            sprintf( Riemann_Name, "sonic rarefaction wave" );
                            break;

      default : Aux_Error( ERROR_INFO, "unsupported Riemann problem (%d) !!\n", Riemann_Prob );
   } // switch ( Riemann_Prob )


// record the test problem parameters
   if ( MPI_Rank == 0 )
   {
      Aux_Message( stdout, "\n" );
      Aux_Message( stdout, "%s test :\n", TestProb );
      Aux_Message( stdout, "=============================================================================\n" );
      Aux_Message( stdout, "NOTE : target problem                   = %s\n",     Riemann_Name );
      Aux_Message( stdout, "       left-state density               = %13.7e\n", Riemann_RhoL );
      Aux_Message( stdout, "       left-state velocity              = %13.7e\n", Riemann_VelL );
      Aux_Message( stdout, "       left-state transverse velocity   = %13.7e\n", Riemann_VelL_T );
      Aux_Message( stdout, "       left-state pressure              = %13.7e\n", Riemann_PreL );
      Aux_Message( stdout, "       right-state density              = %13.7e\n", Riemann_RhoR );
      Aux_Message( stdout, "       right-state velocity             = %13.7e\n", Riemann_VelR );
      Aux_Message( stdout, "       right-state transverse velocity  = %13.7e\n", Riemann_VelR_T );
      Aux_Message( stdout, "       right-state pressure             = %13.7e\n", Riemann_PreR );
      Aux_Message( stdout, "       propagation direction            = %s%s\n",   ( Riemann_LR > 0 ) ? "+" : "-", 
                                                                                 ( Riemann_XYZ == 0 ) ? "x" :
                                                                                 ( Riemann_XYZ == 1 ) ? "y" : "z" );
      Aux_Message( stdout, "=============================================================================\n" );
      Aux_Message( stdout, "\n" );
   }


// set some default parameters
   const long End_Step_Default = __INT_MAX__;

   if ( END_STEP < 0 )
   {
      END_STEP = End_Step_Default;

      if ( MPI_Rank == 0 )
         Aux_Message( stdout, "NOTE : parameter %s is set to %ld in the %s test !!\n", "END_STEP", END_STEP, TestProb );
   }

   if ( END_T < 0.0 )
   {
      END_T = Riemann_EndT;

      if ( MPI_Rank == 0 )
         Aux_Message( stdout, "NOTE : parameter %s is set to %13.7e in the %s test !!\n", "END_T", END_T, TestProb );
   }

   if (  ( Riemann_XYZ == 0 && OPT__OUTPUT_PART != OUTPUT_X )  || 
         ( Riemann_XYZ == 1 && OPT__OUTPUT_PART != OUTPUT_Y )  || 
         ( Riemann_XYZ == 2 && OPT__OUTPUT_PART != OUTPUT_Z )    )
   {
      OPT__OUTPUT_PART = ( Riemann_XYZ == 0 ) ? OUTPUT_X : ( Riemann_XYZ == 1 ) ? OUTPUT_Y : OUTPUT_Z;

      if ( MPI_Rank == 0 )
         Aux_Message( stdout, "NOTE : parameter %s is reset to %d in the %s test !!\n", 
                      "OPT__OUTPUT_PART", OPT__OUTPUT_PART, TestProb );
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
// Function    :  Hydro_TestProbSol_Riemann
// Description :  Assign the initial condition for the HYDRO Riemann problem test 
//
// Note        :  This function is invoked by "Hydro_Init_StartOver_AssignData"
//
// Parameter   :  fluid : Array to store the initial fluid attributes to be returned
//                x/y/z : Target physical coordinates 
//                Time  : Target physical time
//
// Return      :  fluid
//-------------------------------------------------------------------------------------------------------
void Hydro_TestProbSol_Riemann( real fluid[], const real x, const real y, const real z, const double Time )
{

   const real _Gamma_m1 = 1.0/(GAMMA-1.0);
   real r, BoxCen;
   int  TVar[NCOMP];

   switch ( Riemann_XYZ )
   {
      case 0 : r=x; TVar[0]=DENS; TVar[1]=MOMX; TVar[2]=MOMY; TVar[3]=MOMZ; TVar[4]=ENGY; BoxCen=0.5*amr->BoxSize[0]; break;
      case 1 : r=y; TVar[0]=DENS; TVar[1]=MOMY; TVar[2]=MOMZ; TVar[3]=MOMX; TVar[4]=ENGY; BoxCen=0.5*amr->BoxSize[1]; break;
      case 2 : r=z; TVar[0]=DENS; TVar[1]=MOMZ; TVar[2]=MOMX; TVar[3]=MOMY; TVar[4]=ENGY; BoxCen=0.5*amr->BoxSize[2]; break;
      default : Aux_Error( ERROR_INFO, "incorrect parameter %s = %d !!\n", "Riemann_XYZ", Riemann_XYZ );
   }

   if (  ( Riemann_LR > 0 && r < BoxCen )  ||  ( Riemann_LR < 0 && r > BoxCen )  )
   {
      fluid[ TVar[0] ] = Riemann_RhoL;
      fluid[ TVar[1] ] = Riemann_RhoL*Riemann_VelL;
      fluid[ TVar[2] ] = Riemann_RhoL*Riemann_VelL_T;
      fluid[ TVar[3] ] = 0.0;
      fluid[ TVar[4] ] = 0.5*( SQR(fluid[MOMX]) + SQR(fluid[MOMY]) + SQR(fluid[MOMZ]) )/fluid[DENS] + Riemann_PreL*_Gamma_m1;
   }

   else
   {
      fluid[ TVar[0] ] = Riemann_RhoR;
      fluid[ TVar[1] ] = Riemann_RhoR*Riemann_VelR;
      fluid[ TVar[2] ] = Riemann_RhoR*Riemann_VelR_T;
      fluid[ TVar[3] ] = 0.0;
      fluid[ TVar[4] ] = 0.5*( SQR(fluid[MOMX]) + SQR(fluid[MOMY]) + SQR(fluid[MOMZ]) )/fluid[DENS] + Riemann_PreR*_Gamma_m1;
   }

   if ( Riemann_LR < 0 )  
   {
      fluid[ TVar[1] ] = -fluid[ TVar[1] ];
      fluid[ TVar[2] ] = -fluid[ TVar[2] ];
   }

} // FUNCTION : Hydro_TestProbSol_Riemann



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

   getline( &input_line, &len, File );
   sscanf( input_line, "%d%s",   &temp_int,                 string );
   Riemann_Prob = (Riemann_t)temp_int;

   getline( &input_line, &len, File );
   getline( &input_line, &len, File );
   getline( &input_line, &len, File );
   getline( &input_line, &len, File );
   getline( &input_line, &len, File );

   getline( &input_line, &len, File );
   sscanf( input_line, "%d%s",   &Riemann_LR,               string );

   getline( &input_line, &len, File );
   sscanf( input_line, "%d%s",   &Riemann_XYZ,              string );

   fclose( File );
   if ( input_line != NULL )     free( input_line );


// check
   if ( Riemann_Prob < 0  ||  Riemann_Prob > 5 )
      Aux_Error( ERROR_INFO, "incorrect parameter \"%s = %d\" !!\n", "Riemann_Prob", Riemann_Prob );

   if ( Riemann_XYZ < 0  ||  Riemann_XYZ > 2 )
      Aux_Error( ERROR_INFO, "incorrect parameter \"%s = %d\" !!\n", "Riemann_XYZ", Riemann_XYZ );

} // FUNCTION : LoadTestProbParameter



#endif // #if ( MODEL == HYDRO )
