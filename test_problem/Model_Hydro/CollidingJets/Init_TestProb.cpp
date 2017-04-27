#include "GAMER.h"

#if ( MODEL == HYDRO )



extern void (*Init_Function_Ptr)( real fluid[], const double x, const double y, const double z, const double Time );
extern void (*Output_TestProbErr_Ptr)( const bool BaseOnly );

static void HYDRO_TestProbSol_Jet( real fluid[], const double x, const double y, const double z, const double Time );
static void LoadTestProbParameter();


// global variables in the HYDRO colliding jets test
// =======================================================================================
double   Jet_BgDens;                   // ambient density
double   Jet_BgTemp;                   // ambient temperature
int      Jet_NJet;                     // number of jets (1/2)
double  *Jet_Radius        = NULL;     // radius of the cylinder-shape jet source
double  *Jet_HalfHeight    = NULL;     // half height of the cylinder-shape jet source
double  *Jet_SrcVel        = NULL;     // jet velocity
double  *Jet_SrcDens       = NULL;     // jet density
double  *Jet_SrcTemp       = NULL;     // jet temperature
double (*Jet_Vec)[3]       = NULL;     // jet orientation vector (x,y,z) (NOT necessary to be a unit vector)
double (*Jet_CenOffset)[3] = NULL;     // jet central coordinates offset

double   Jet_BgEint;                   // ambient internal energy
double  *Jet_SrcEint       = NULL;     // jet internal energy
double (*Jet_Cen)[3]       = NULL;     // jet central coordinates
double  *Jet_WaveK         = NULL;     // jet wavenumber used in the sin() function to have smooth bidirectional jets
// =======================================================================================




//-------------------------------------------------------------------------------------------------------
// Function    :  Init_TestProb
// Description :  Initialize parameters for the HYDRO colliding jets test
//
// Note        :  1. Please copy this file to "GAMER/src/Init/Init_TestProb.cpp"
//                2. Global variables declared here will also be used in the function
//                   "HYDRO_TestProbSol_Jet"
//
// Parameter   :  None
//-------------------------------------------------------------------------------------------------------
void Init_TestProb()
{

   const char *TestProb = "HYDRO colliding jets";

// check
#  if ( MODEL != HYDRO )
#  error : ERROR : "MODEL != HYDRO" in the HYDRO colliding jets test !!
#  endif

#  ifdef GRAVITY
#  error : ERROR : "GRAVITY must be OFF" in the HYDRO colliding jets test !!
#  endif

#  ifdef COMOVING
#  error : ERROR : "COMOVING must be OFF" in the HYDRO colliding jets test !!
#  endif

#  ifdef PARTICLE
#  error : ERROR : "PARTICLE must be OFF" in the HYDRO colliding jets test !!
#  endif


// set the initialization and output functions
   Init_Function_Ptr      = HYDRO_TestProbSol_Jet;
   Output_TestProbErr_Ptr = NULL;


// load the test problem parameters
   LoadTestProbParameter();


// set the derived jet parameters
   Jet_BgEint = Jet_BgDens*Jet_BgTemp/(MOLECULAR_WEIGHT*Const_amu/UNIT_M) / ( GAMMA - (real)1.0 );

   for (int n=0; n<Jet_NJet; n++)
   {
      Jet_SrcEint[n] = Jet_SrcDens[n]*Jet_SrcTemp[n]/(MOLECULAR_WEIGHT*Const_amu/UNIT_M) / ( GAMMA - (real)1.0 );

      for (int d=0; d<3; d++)    Jet_Cen[n][d] = 0.5*amr->BoxSize[d] + Jet_CenOffset[n][d];

      Jet_WaveK[n] = 0.5*M_PI/Jet_HalfHeight[n];
   }


// record the test problem parameters
   if ( MPI_Rank == 0 )
   {
      Aux_Message( stdout, "\n" );
      Aux_Message( stdout, "%s test :\n", TestProb );
      Aux_Message( stdout, "=============================================================================\n" );
      Aux_Message( stdout, "Jet_BgDens          = % 14.7e g/cm^3\n", Jet_BgDens*UNIT_D                     );
      Aux_Message( stdout, "Jet_BgTemp          = % 14.7e keV\n",    Jet_BgTemp*UNIT_E/Const_keV           );
      Aux_Message( stdout, "Jet_NJet            = %d\n",             Jet_NJet                              );
      for (int n=0; n<Jet_NJet; n++) {
      Aux_Message( stdout, "Jet # %d\n", n );
      Aux_Message( stdout, "   Jet_Radius       = % 14.7e kpc\n",    Jet_Radius    [n]*UNIT_L/Const_kpc    );
      Aux_Message( stdout, "   Jet_HalfHeight   = % 14.7e kpc\n",    Jet_HalfHeight[n]*UNIT_L/Const_kpc    );
      Aux_Message( stdout, "   Jet_SrcVel       = % 14.7e cm/s\n",   Jet_SrcVel    [n]*UNIT_V              );
      Aux_Message( stdout, "   Jet_SrcDens      = % 14.7e g/cm^3\n", Jet_SrcDens   [n]*UNIT_D              );
      Aux_Message( stdout, "   Jet_SrcTemp      = % 14.7e keV\n",    Jet_SrcTemp   [n]*UNIT_E/Const_keV    );
      Aux_Message( stdout, "   Jet_Vec[x]       = % 14.7e\n",        Jet_Vec       [n][0]                  );
      Aux_Message( stdout, "   Jet_Vec[y]       = % 14.7e\n",        Jet_Vec       [n][1]                  );
      Aux_Message( stdout, "   Jet_Vec[z]       = % 14.7e\n",        Jet_Vec       [n][2]                  );
      Aux_Message( stdout, "   Jet_CenOffset[x] = % 14.7e kpc\n",    Jet_CenOffset [n][0]*UNIT_L/Const_kpc );
      Aux_Message( stdout, "   Jet_CenOffset[y] = % 14.7e kpc\n",    Jet_CenOffset [n][1]*UNIT_L/Const_kpc );
      Aux_Message( stdout, "   Jet_CenOffset[z] = % 14.7e kpc\n",    Jet_CenOffset [n][2]*UNIT_L/Const_kpc );
      }
      Aux_Message( stdout, "=============================================================================\n" );
      Aux_Message( stdout, "\n" );
   }


// set some default parameters
   const double End_T_Default    = 50.0*Const_Myr/UNIT_T;
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
// Function    :  HYDRO_TestProbSol_Jet
// Description :  Calculate the initial condition in the HYDRO colliding jets test
//
// Note        :  1. This function is invoked by "HYDRO_Init_StartOver_AssignData"
//                2. Here we only set the initial condition for the **ambient gas**
//                   --> The source of jet(s) is set by Flu_ResetByUser()
//
// Parameter   :  fluid : Array to store the analytical solution to be returned
//                x/y/z : Target physical coordinates
//                Time  : Target physical time
//
// Return      :  fluid
//-------------------------------------------------------------------------------------------------------
void HYDRO_TestProbSol_Jet( real fluid[], const double x, const double y, const double z, const double Time )
{

   fluid[DENS] = Jet_BgDens;
   fluid[MOMX] = 0.0;
   fluid[MOMY] = 0.0;
   fluid[MOMZ] = 0.0;
   fluid[ENGY] = Jet_BgEint + 0.5*( SQR(fluid[MOMX]) + SQR(fluid[MOMY]) + SQR(fluid[MOMZ]) ) / fluid[DENS];

} // FUNCTION : HYDRO_TestProbSol_Jet



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

   if ( File == NULL )  Aux_Error( ERROR_INFO, "file \"%s\" does not exist !!\n", FileName );

   int    tmp_int;
   char  *input_line = NULL;
   char   string[100];
   size_t len = 0;

   getline( &input_line, &len, File );

// load the input parameters independent of the number of jets
   getline( &input_line, &len, File );
   sscanf( input_line, "%lf%s",  &Jet_BgDens,            string );

   getline( &input_line, &len, File );
   sscanf( input_line, "%lf%s",  &Jet_BgTemp,            string );

   getline( &input_line, &len, File );
   sscanf( input_line, "%d%s",   &Jet_NJet,              string );


// check the number of jets
   if ( Jet_NJet < 1 )  Aux_Error( ERROR_INFO, "%s = %d < 1 !!\n", "Jet_NJet", Jet_NJet );


// allocate memory for variables whose size depends on the number of jets
   Jet_Radius     = new double [Jet_NJet];
   Jet_HalfHeight = new double [Jet_NJet];
   Jet_SrcVel     = new double [Jet_NJet];
   Jet_SrcDens    = new double [Jet_NJet];
   Jet_SrcTemp    = new double [Jet_NJet];
   Jet_Vec        = new double [Jet_NJet][3];
   Jet_CenOffset  = new double [Jet_NJet][3];
   Jet_SrcEint    = new double [Jet_NJet];
   Jet_Cen        = new double [Jet_NJet][3];
   Jet_WaveK      = new double [Jet_NJet];


// load the input parameters of each jet
   for (int n=0; n<Jet_NJet; n++) {

   getline( &input_line, &len, File );
   sscanf( input_line, "%lf%s",  &Jet_Radius[n],         string );

   getline( &input_line, &len, File );
   sscanf( input_line, "%lf%s",  &Jet_HalfHeight[n],     string );

   getline( &input_line, &len, File );
   sscanf( input_line, "%lf%s",  &Jet_SrcVel[n],         string );

   getline( &input_line, &len, File );
   sscanf( input_line, "%lf%s",  &Jet_SrcDens[n],        string );

   getline( &input_line, &len, File );
   sscanf( input_line, "%lf%s",  &Jet_SrcTemp[n],        string );

   for (int d=0; d<3; d++) {
   getline( &input_line, &len, File );
   sscanf( input_line, "%lf%s",  &Jet_Vec[n][d],         string ); }

   for (int d=0; d<3; d++) {
   getline( &input_line, &len, File );
   sscanf( input_line, "%lf%s",  &Jet_CenOffset[n][d],   string ); }

   } // for (int n=0; n<Jet_NJet; n++)

   fclose( File );
   if ( input_line != NULL )     free( input_line );


// check
   if ( Jet_BgDens <= 0.0 )   Aux_Error( ERROR_INFO, "%s = %14.7e <= 0.0 !!\n", "Jet_BgDens", Jet_BgDens );
   if ( Jet_BgTemp <= 0.0 )   Aux_Error( ERROR_INFO, "%s = %14.7e <= 0.0 !!\n", "Jet_BgTemp", Jet_BgTemp );

   for (int n=0; n<Jet_NJet; n++) {

   if ( Jet_Radius    [n] <= 0.0 )  Aux_Error( ERROR_INFO, "jet %d: %s = %14.7e <= 0.0 !!\n", n, "Jet_Radius",     Jet_Radius    [n] );
   if ( Jet_HalfHeight[n] <= 0.0 )  Aux_Error( ERROR_INFO, "jet %d: %s = %14.7e <= 0.0 !!\n", n, "Jet_HalfHeight", Jet_HalfHeight[n] );
   if ( Jet_SrcVel    [n] <= 0.0 )  Aux_Error( ERROR_INFO, "jet %d: %s = %14.7e <= 0.0 !!\n", n, "Jet_SrcVel",     Jet_SrcVel    [n] );
   if ( Jet_SrcDens   [n] <= 0.0 )  Aux_Error( ERROR_INFO, "jet %d: %s = %14.7e <= 0.0 !!\n", n, "Jet_SrcDens",    Jet_SrcDens   [n] );
   if ( Jet_SrcTemp   [n] <= 0.0 )  Aux_Error( ERROR_INFO, "jet %d: %s = %14.7e <= 0.0 !!\n", n, "Jet_SrcTemp",    Jet_SrcTemp   [n] );

   }


// convert to code units
   Jet_BgDens           *= 1.0       / UNIT_D;
   Jet_BgTemp           *= Const_keV / UNIT_E;

   for (int n=0; n<Jet_NJet; n++) {

   Jet_Radius    [n]    *= Const_kpc / UNIT_L;
   Jet_HalfHeight[n]    *= Const_kpc / UNIT_L;
   Jet_SrcVel    [n]    *= 1.0       / UNIT_V;
   Jet_SrcDens   [n]    *= 1.0       / UNIT_D;
   Jet_SrcTemp   [n]    *= Const_keV / UNIT_E;
   for (int d=0; d<3; d++) {
   Jet_CenOffset [n][d] *= Const_kpc / UNIT_L; }
   }

} // FUNCTION : LoadTestProbParameter



#endif // #if ( MODEL == HYDRO )
