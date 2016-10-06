#include "GAMER.h"

#ifdef PARTICLE



extern void (*Init_Function_Ptr)( real fluid[], const double x, const double y, const double z, const double Time );
extern void (*Output_TestProbErr_Ptr)( const bool BaseOnly );

static void LoadTestProbParameter();
static void Par_TestProbSol_ClusterMerger( real fluid[], const double x, const double y, const double z, const double Time );

//extern double MassProf_ClusterMerger( const double r );


// global variables in the Cluster Merger test
// =======================================================================================
int    ClusterMerger_RanSeed;          // random seed for setting particle position and velocity
double ClusterMerger_Mvir;             // total: virial mass
double ClusterMerger_Rvir;             // total: virial radius (and also the cut-off radius)
double ClusterMerger_DM_C;             // dark matter: NFW concentration parameter
double ClusterMerger_DM_Rzero;         // dark matter: radius at which the dark matter velocity dispersion is set to zero
double ClusterMerger_Gas_Rcore;        // gas: core radius in the beta model
double ClusterMerger_Gas_TvirM14;      // gas: temperature at the virial radius for a 1e14 Msun cluster
double ClusterMerger_Gas_MTScaling;    // gas: M-T scaling index for temperature normalization --> Tvir ~ Mvir^MT_Scaling
double ClusterMerger_Gas_MFrac;        // gas: mass fraction
int    ClusterMerger_NBin_MassProf;    // number of radial bins for the mass profile table
bool   ClusterMerger_Coll;             // (true/false) --> test (cluster merger / single cluster)

/*
double ClusterMerger_Coll_D;           // distance between two ClusterMerger clouds for the ClusterMerger collision test
double ClusterMerger_Coll_ImpactPara;  // impact parameter
double ClusterMerger_Coll_BulkVel[3];  // bulk velocity
*/

double ClusterMerger_Gas_Tvir;         // gas: temperature at the virial radius for the target cluster mass
/*
double ClusterMerger_Rho0;                   // peak density
double ClusterMerger_R0;                     // scale radius
*/
// =======================================================================================




//-------------------------------------------------------------------------------------------------------
// Function    :  Init_TestProb
// Description :  Initialize parameters for the cluster merger test
//
// Note        :  1. Please copy this file to "GAMER/src/Init/Init_TestProb.cpp"
//                2. Test problem parameters can be set in the input file "Input__TestProb"
//
// Parameter   :  None
//-------------------------------------------------------------------------------------------------------
void Init_TestProb()
{

   const char *TestProb = "cluster merger";

// check
# if ( MODEL != HYDRO )
# error : ERROR : only support the HYDRO model !!
# endif

# ifndef PARTICLE
# error : ERROR : "PARTICLE must be ON" in the cluster merger test !!
# endif

# ifndef GRAVITY
# error : ERROR : "GRAVITY must be ON" in the cluster merger test !!
# endif

# ifdef COMOVING
# error : ERROR : "COMOVING must be OFF" in the cluster merger test !!
# endif

   if ( !OPT__UNIT )
      Aux_Error( ERROR_INFO, "please turn on \"OPT__UNIT\" and set units properly for this test problem !!\n" );


// set the initialization and output functions
   Init_Function_Ptr      = Par_TestProbSol_ClusterMerger;
   Output_TestProbErr_Ptr = NULL;


// load the test problem parameters
   LoadTestProbParameter();


// convert parameters to code units
   ClusterMerger_Mvir        *= Const_Msun / UNIT_M;
   ClusterMerger_Rvir        *= Const_Mpc  / UNIT_L;
   ClusterMerger_DM_Rzero    *= Const_Mpc  / UNIT_L;
   ClusterMerger_Gas_Rcore   *= Const_Mpc  / UNIT_L;
   ClusterMerger_Gas_TvirM14 *= Const_keV  / UNIT_E;


// set the test problem parameters
// (1) temperature normalization for the target cluster mass
   ClusterMerger_Gas_Tvir = ClusterMerger_Gas_TvirM14
                            *pow(  ( ClusterMerger_Mvir / (1.0e14*Const_Msun/UNIT_M) ), ClusterMerger_Gas_MTScaling  );


// record the test problem parameters
   if ( MPI_Rank == 0 )
   {
      Aux_Message( stdout, "\n" );
      Aux_Message( stdout, "%s test :\n", TestProb );
      Aux_Message( stdout, "=============================================================================\n" );
      Aux_Message( stdout, " Note: random seed for initializing particles             = %d\n",          ClusterMerger_RanSeed );
      Aux_Message( stdout, "       [total]       virial mass                          = %13.7e Msun\n", ClusterMerger_Mvir*UNIT_M/Const_Msun );
      Aux_Message( stdout, "       [total]       virial radius                        = %13.7e Mpc\n",  ClusterMerger_Rvir*UNIT_L/Const_Mpc );
      Aux_Message( stdout, "       [dark matter] concentration parameter              = %13.7e\n",      ClusterMerger_DM_C );
      Aux_Message( stdout, "       [dark matter] radius with zero velocity dispersion = %13.7e Mpc\n",  ClusterMerger_DM_Rzero*UNIT_L/Const_Mpc );
      Aux_Message( stdout, "       [gas]         core radius in beta model            = %13.7e Mpc\n",  ClusterMerger_Gas_Rcore*UNIT_L/Const_Mpc );
      Aux_Message( stdout, "       [gas]         temperature normalization at Rvir    = %13.7e keV\n",  ClusterMerger_Gas_Tvir*UNIT_E/Const_keV );
      Aux_Message( stdout, "                                                          = %13.7e K\n",    ClusterMerger_Gas_Tvir*UNIT_E/Const_kB );
      Aux_Message( stdout, "       [gas]         M-T scaling index                    = %13.7e\n",      ClusterMerger_Gas_MTScaling );
      Aux_Message( stdout, "       [gas]         mass fraction                        = %13.7e\n",      ClusterMerger_Gas_MFrac );
      Aux_Message( stdout, "       number of bins for interpolating mass profile      = %d\n",          ClusterMerger_NBin_MassProf );
      Aux_Message( stdout, "       test mode                                          = %s\n",          (ClusterMerger_Coll)?
                                                                                                        "cluster merger":"single cluster" );
      /*
      if ( ClusterMerger_Coll )
      Aux_Message( stdout, "       initial distance between two clouds          = %13.7e\n",  ClusterMerger_Coll_D );
      for (int d=0; d<3; d++)
      Aux_Message( stdout, "       bulk velocity [%d]                            = %14.7e\n", d, ClusterMerger_BulkVel[d] );
      */
      Aux_Message( stdout, "=============================================================================\n" );
      Aux_Message( stdout, "\n" );
   }


// set some default parameters
   const double End_T_Default    = 10.0*Const_Gyr/UNIT_T;
   const long   End_Step_Default = __INT_MAX__;

   if ( END_STEP < 0 )
   {
      END_STEP = End_Step_Default;

      if ( MPI_Rank == 0 )
         Aux_Message( stdout, "NOTE : parameter %s is set to %ld in the %s test\n", "END_STEP", END_STEP, TestProb );
   }

   if ( END_T < 0.0 )
   {
      END_T = End_T_Default;

      if ( MPI_Rank == 0 )
         Aux_Message( stdout, "NOTE : parameter %s is set to %13.7e Gyr in the %s test\n", "END_T", END_T*UNIT_T/Const_Gyr, TestProb );
   }

   if ( OPT__OUTPUT_TEST_ERROR )
   {
      OPT__OUTPUT_TEST_ERROR = false;

      if ( MPI_Rank == 0 )
         Aux_Message( stderr, "WARNING : parameter %s is reset to %d in the %s test !!\n",
                      "OPT__OUTPUT_TEST_ERROR", OPT__OUTPUT_TEST_ERROR, TestProb );
   }

   if ( OPT__INIT == INIT_STARTOVER  &&  amr->Par->Init != PAR_INIT_BY_FUNCTION )
      Aux_Error( ERROR_INFO, "please set PAR_INIT = PAR_INIT_BY_FUNCTION !!\n" );

} // FUNCTION : Init_TestProb



//-------------------------------------------------------------------------------------------------------
// Function    :  Par_TestProbSol_ClusterMerger
// Description :  Initialize the background density field as zero for the cluster merger test
//
// Note        :  1. Isothermal beta model for gas
//                   --> Currently beta is fixed to 1
//
// Parameter   :  fluid : Fluid field to be initialized
//                x/y/z : Target physical coordinates
//                Time  : Target physical time
//
// Return      :  fluid
//-------------------------------------------------------------------------------------------------------
void Par_TestProbSol_ClusterMerger( real *fluid, const double x, const double y, const double z, const double Time )
{

   /*
#  if   ( MODEL == HYDRO )
// gas share the same density profile as particles (except for different total masses)
   const double TotM    = 4.0/3.0*M_PI*CUBE(Plummer_R0)*Plummer_Rho0;
   const double GasRho0 = Plummer_Rho0*Plummer_GasMFrac;
   const double PresBg  = 0.0;   // background pressure
   double r2, a2;

   if ( Plummer_Coll )
   {
      const double Coll_Offset = 0.5*Plummer_Coll_D/sqrt(3.0);
      double Center[3];

      fluid[DENS] = 0.0;
      fluid[ENGY] = 0.0;

      for (int t=-1; t<=1; t+=2)
      {
         for (int d=0; d<3; d++)    Center[d] = Plummer_Center[d] + Coll_Offset*(double)t;

         r2 = SQR(x-Center[0])+ SQR(y-Center[1]) + SQR(z-Center[2]);
         a2 = r2 / SQR(Plummer_R0);

         fluid[DENS] += GasRho0 * pow( 1.0 + a2, -2.5 );
         fluid[ENGY] += (  NEWTON_G*TotM*GasRho0 / ( 6.0*Plummer_R0*CUBE(1.0 + a2) ) + PresBg  ) / ( GAMMA - 1.0 );
      }

      fluid[MOMX]  = fluid[DENS]*Plummer_BulkVel[0];
      fluid[MOMY]  = fluid[DENS]*Plummer_BulkVel[1];
      fluid[MOMZ]  = fluid[DENS]*Plummer_BulkVel[2];
      fluid[ENGY] += 0.5*( SQR(fluid[MOMX]) + SQR(fluid[MOMY]) + SQR(fluid[MOMZ]) ) / fluid[DENS];
   }

   else
   {
      r2 = SQR(x-Plummer_Center[0])+ SQR(y-Plummer_Center[1]) + SQR(z-Plummer_Center[2]);
      a2 = r2 / SQR(Plummer_R0);

      fluid[DENS] = GasRho0 * pow( 1.0 + a2, -2.5 );
      fluid[MOMX] = fluid[DENS]*Plummer_BulkVel[0];
      fluid[MOMY] = fluid[DENS]*Plummer_BulkVel[1];
      fluid[MOMZ] = fluid[DENS]*Plummer_BulkVel[2];
      fluid[ENGY] = (  NEWTON_G*TotM*GasRho0 / ( 6.0*Plummer_R0*CUBE(1.0 + a2) ) + PresBg  ) / ( GAMMA - 1.0 )
                    + 0.5*( SQR(fluid[MOMX]) + SQR(fluid[MOMY]) + SQR(fluid[MOMZ]) ) / fluid[DENS];
   } // if ( Plummer_Coll ) ... else ...

#  elif ( MODEL == ELBDM )
// set wave function as zero everywhere
   fluid[REAL] = 0.0;
   fluid[IMAG] = 0.0;
   fluid[DENS] = 0.0;

#  else
#  error : unsupported MODEL !!
#  endif
   */

} // FUNCTION : Par_TestProbSol_Plummer



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


// skip the header
   getline( &input_line, &len, File );

   getline( &input_line, &len, File );
   sscanf( input_line, "%d%s",   &ClusterMerger_RanSeed,       string );

   getline( &input_line, &len, File );
   sscanf( input_line, "%lf%s",  &ClusterMerger_Mvir,          string );

   getline( &input_line, &len, File );
   sscanf( input_line, "%lf%s",  &ClusterMerger_Rvir,          string );

   getline( &input_line, &len, File );
   sscanf( input_line, "%lf%s",  &ClusterMerger_DM_C,          string );

   getline( &input_line, &len, File );
   sscanf( input_line, "%lf%s",  &ClusterMerger_DM_Rzero,      string );

   getline( &input_line, &len, File );
   sscanf( input_line, "%lf%s",  &ClusterMerger_Gas_Rcore,     string );

   getline( &input_line, &len, File );
   sscanf( input_line, "%lf%s",  &ClusterMerger_Gas_TvirM14,   string );

   getline( &input_line, &len, File );
   sscanf( input_line, "%lf%s",  &ClusterMerger_Gas_MTScaling, string );

   getline( &input_line, &len, File );
   sscanf( input_line, "%lf%s",  &ClusterMerger_Gas_MFrac,     string );

   getline( &input_line, &len, File );
   sscanf( input_line, "%d%s",   &ClusterMerger_NBin_MassProf, string );

   getline( &input_line, &len, File );
   sscanf( input_line, "%d%s",   &temp_int,                    string );
   ClusterMerger_Coll = (bool)temp_int;

   fclose( File );
   if ( input_line != NULL )     free( input_line );


// set the default parameters
   if ( ClusterMerger_DM_Rzero <= 0.0 )
   {
      ClusterMerger_DM_Rzero = 1.0e4*ClusterMerger_Rvir;

      if ( MPI_Rank == 0 )  Aux_Message( stdout, "NOTE : parameter \"%s\" is set to the default value = %13.7e\n",
                                         "ClusterMerger_DM_Rzero", ClusterMerger_DM_Rzero );
   }

   if ( ClusterMerger_Gas_TvirM14 <= 0.0 )
   {
      ClusterMerger_Gas_TvirM14 = 2.5;    // in keV

      if ( MPI_Rank == 0 )  Aux_Message( stdout, "NOTE : parameter \"%s\" is set to the default value = %13.7e keV\n",
                                         "ClusterMerger_Gas_TvirM14", ClusterMerger_Gas_TvirM14 );
   }

   if ( ClusterMerger_Gas_MTScaling <= 0.0 )
   {
      ClusterMerger_Gas_MTScaling = 2.0/3.0;

      if ( MPI_Rank == 0 )  Aux_Message( stdout, "NOTE : parameter \"%s\" is set to the default value = %13.7e\n",
                                         "ClusterMerger_Gas_MTScaling", ClusterMerger_Gas_MTScaling );
   }

   if ( ClusterMerger_NBin_MassProf <= 1 )
   {
      ClusterMerger_NBin_MassProf = 10000;

      if ( MPI_Rank == 0 )  Aux_Message( stdout, "NOTE : parameter \"%s\" is set to the default value = %d\n",
                                         "ClusterMerger_NBin_MassProf", ClusterMerger_NBin_MassProf );
   }


// check
   if ( ClusterMerger_RanSeed < 0 )
      Aux_Error( ERROR_INFO, "ClusterMerger_RanSeed (%d) < 0 !!\n", ClusterMerger_RanSeed );

   if ( ClusterMerger_Mvir <= 0.0 ) 
      Aux_Error( ERROR_INFO, "ClusterMerger_Mvir (%14.7e) <= 0.0 !!\n", ClusterMerger_Mvir );

   if ( ClusterMerger_Rvir <= 0.0 )
      Aux_Error( ERROR_INFO, "ClusterMerger_Rvir (%14.7e) <= 0.0 !!\n", ClusterMerger_Rvir );

   if ( ClusterMerger_DM_C <= 0.0 )
      Aux_Error( ERROR_INFO, "ClusterMerger_DM_C (%14.7e) <= 0.0 !!\n", ClusterMerger_DM_C );

   if ( ClusterMerger_DM_Rzero < ClusterMerger_Rvir )
      Aux_Error( ERROR_INFO, "ClusterMerger_DM_Rzero (%14.7e) < ClusterMerger_Rvir (%14.7e) !!\n",
                 ClusterMerger_DM_Rzero, ClusterMerger_Rvir );

   if ( ClusterMerger_Gas_Rcore >= ClusterMerger_Rvir )
      Aux_Error( ERROR_INFO, "ClusterMerger_Gas_Rcore (%14.7e) >= ClusterMerger_Rvir (%14.7e) !!\n",
                 ClusterMerger_Gas_Rcore, ClusterMerger_Rvir );

   if ( ClusterMerger_Gas_TvirM14 <= 0.0 )
      Aux_Error( ERROR_INFO, "ClusterMerger_Gas_TvirM14 (%14.7e) <= 0.0 !!\n", ClusterMerger_Gas_TvirM14 );

   if ( ClusterMerger_Gas_MTScaling <= 0.0 )
      Aux_Error( ERROR_INFO, "ClusterMerger_Gas_MTScaling (%14.7e) <= 0.0 !!\n", ClusterMerger_Gas_MTScaling );

   if ( ClusterMerger_Gas_MFrac <= 0.0  ||  ClusterMerger_Gas_MFrac > 1.0 )
      Aux_Error( ERROR_INFO, "ClusterMerger_Gas_MFrac (%14.7e) is not within the correct range [0.0 < x <= 1.0] !!\n",
                 ClusterMerger_Gas_MFrac );

   if ( ClusterMerger_NBin_MassProf <= 1   )
      Aux_Error( ERROR_INFO, "ClusterMerger_NBin_MassProf (%d) <= 1 !!\n", ClusterMerger_NBin_MassProf );

} // FUNCTION : LoadTestProbParameter



#endif // #ifdef PARTICLE
