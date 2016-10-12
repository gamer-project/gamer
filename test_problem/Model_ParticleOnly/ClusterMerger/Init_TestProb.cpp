#include "GAMER.h"
#ifdef SUPPORT_GSL
#include <gsl/gsl_integration.h>
#else
#error : ERROR : please turn on SUPPORT_GSL for the cluster merger test problem !!
#endif

#if ( MODEL == HYDRO  &&  defined PARTICLE )



extern void (*Init_Function_Ptr)( real fluid[], const double x, const double y, const double z, const double Time );
extern void (*Output_TestProbErr_Ptr)( const bool BaseOnly );

static void LoadTestProbParameter();
static void Par_TestProbSol_ClusterMerger( real fluid[], const double x, const double y, const double z, const double Time );

static void SetTable_Gas_PresProf( const int NBin, double *r, double *Pres );
static void SetTable_DM_MassProf ( const int NBin, double *r, double *Mass );
static void SetTable_DM_SigmaProf( const int NBin, double *r, double *Sigma );

static double GSL_IntFunc_Gas_PresProf( double r, void *IntPara );
static double GSL_IntFunc_DM_SigmaProf( double r, void *IntPara );

static double MassProf_Total( const double r );
       double MassProf_DM( const double r );
static double MassProf_Gas( const double r );
static double DensProf_DM( const double r );
static double DensProf_Gas( const double r );


// global variables in the Cluster Merger test
// =======================================================================================
int    ClusterMerger_RanSeed;          // random seed for setting particle position and velocity
double ClusterMerger_Mvir;             // total: virial mass
double ClusterMerger_Rvir;             // total: virial radius for both gas and dark matter (i.e., total mass within Rvir = Mvir)
double ClusterMerger_Rcut;             // total: cut-off radius for both gas and dark matter
                                       // --> outside Rcut, gas density and pressure are set to the values at Rcut and dark matter
                                       //     density is set to zero (i.e., no dark matter particles locate outside Rcut)
                                       // --> note that Rcut is different from Rvir, where the latter is the radius enclosing Mvir
double ClusterMerger_DM_C;             // dark matter: NFW concentration parameter
double ClusterMerger_DM_Rzero;         // dark matter: radius at which the dark matter velocity dispersion is set to zero
double ClusterMerger_Gas_Rcore;        // gas: core radius in the beta model
bool   ClusterMerger_Gas_NormTvir;     // gas: (true/false) => normalize temperature at Rvir by using (TvirM14/pressure=0 at Rzero)
double ClusterMerger_Gas_TvirM14;      // gas: temperature at the virial radius for a 1e14 Msun cluster
double ClusterMerger_Gas_MTScaling;    // gas: M-T scaling index for temperature normalization
                                       // --> Tvir ~ Mvir^MT_Scaling, where Mvir is the **total** cluster mass
double ClusterMerger_Gas_MFrac;        // gas: mass fraction
double ClusterMerger_Gas_DensBgRatio;  // gas: background density  normalized to the density  at Rcut
double ClusterMerger_Gas_PresBgRatio;  // gas: background pressure normalized to the pressure at Rcut
int    ClusterMerger_NBin_PresProf;    // number of radial bins for the table of gas pressure profile
int    ClusterMerger_NBin_MassProf;    // number of radial bins for the table of dark matter mass profile
int    ClusterMerger_NBin_SigmaProf;   // number of radial bins for the table of dark matter velocity dispersion profile
double ClusterMerger_IntRmin;          // minimum radius in the interpolation table
bool   ClusterMerger_Coll;             // (true/false) --> test (cluster merger / single cluster)

/*
double ClusterMerger_Coll_D;           // distance between two ClusterMerger clouds for the ClusterMerger collision test
double ClusterMerger_Coll_ImpactPara;  // impact parameter
double ClusterMerger_Coll_BulkVel[3];  // bulk velocity
*/

double ClusterMerger_DM_M;             // dark matter: total dark matter mass
double ClusterMerger_DM_Rs;            // dark matter: NFW scale factor
double ClusterMerger_DM_Rho0;          // dark matter: NFW density parameter
double ClusterMerger_Gas_M;            // gas: total gas mass
double ClusterMerger_Gas_Rho0;         // gas: density parameter in the isothermal beta model
double ClusterMerger_Gas_Tvir;         // gas: temperature at the virial radius for the target cluster mass

double *ClusterMerger_Gas_PresProf   = NULL;    // table of gas pressure for interpolation
double *ClusterMerger_DM_MassProf    = NULL;    // table of dark matter mass profile for interpolation
double *ClusterMerger_DM_SigmaProf   = NULL;    // table of dark matter velocity dispersion profile for interpolation

double *ClusterMerger_Gas_PresProf_R = NULL;    // radii of ClusterMerger_Gas_PresProf table
double *ClusterMerger_DM_MassProf_R  = NULL;    // radii of ClusterMerger_DM_MassProf table
double *ClusterMerger_DM_SigmaProf_R = NULL;    // radii of ClusterMerger_DM_SigmaProf table
// =======================================================================================

// GSL parameters
#ifdef SUPPORT_GSL
const int    GSL_IntRule  = GSL_INTEG_GAUSS41;     // Gauss-Kronrod integration rules (options: 15,21,31,41,51,61)
const size_t GSL_WorkSize = 1000000;               // work size used by GSL

gsl_integration_workspace *GSL_WorkSpace = NULL;
#endif




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
#  if ( MODEL != HYDRO )
#  error : ERROR : only support the HYDRO model !!
#  endif

#  ifndef PARTICLE
#  error : ERROR : "PARTICLE must be ON" in the cluster merger test !!
#  endif

#  ifndef GRAVITY
#  error : ERROR : "GRAVITY must be ON" in the cluster merger test !!
#  endif

#  ifdef COMOVING
#  error : ERROR : "COMOVING must be OFF" in the cluster merger test !!
#  endif

#  ifndef SUPPORT_GSL
#  error : ERROR : "SUPPORT_GSL must be ON" in the cluster merger test !!
#  endif

   if ( !OPT__UNIT )
      Aux_Error( ERROR_INFO, "please turn on \"OPT__UNIT\" and set units properly for this test problem !!\n" );

   if ( OPT__BC_FLU[0] == BC_FLU_PERIODIC  ||  OPT__BC_POT == BC_POT_PERIODIC )
      Aux_Error( ERROR_INFO, "one should not adopt periodic boundary condition for this test problem !!\n" );


// set the initialization and output functions
   Init_Function_Ptr      = Par_TestProbSol_ClusterMerger;
   Output_TestProbErr_Ptr = NULL;


// load the test problem parameters
   LoadTestProbParameter();


// convert parameters to code units
   ClusterMerger_Mvir        *= Const_Msun / UNIT_M;
   ClusterMerger_Rvir        *= Const_Mpc  / UNIT_L;
   ClusterMerger_Rcut        *= Const_Mpc  / UNIT_L;
   ClusterMerger_DM_Rzero    *= Const_Mpc  / UNIT_L;
   ClusterMerger_Gas_Rcore   *= Const_Mpc  / UNIT_L;
   ClusterMerger_Gas_TvirM14 *= Const_keV  / UNIT_E;
   ClusterMerger_IntRmin     *= Const_Mpc  / UNIT_L;


// set the test problem parameters
// (1) temperature normalization for the target cluster mass
   ClusterMerger_Gas_Tvir = ClusterMerger_Gas_TvirM14
                            *pow(  ( ClusterMerger_Mvir / (1.0e14*Const_Msun/UNIT_M) ), ClusterMerger_Gas_MTScaling  );

// (2) NFW
   const double C = ClusterMerger_DM_C;
   ClusterMerger_DM_M      = ClusterMerger_Mvir * ( 1.0 - ClusterMerger_Gas_MFrac );
   ClusterMerger_DM_Rs     = ClusterMerger_Rvir / C;
   ClusterMerger_DM_Rho0   = ClusterMerger_DM_M / (  4.0*M_PI*CUBE(ClusterMerger_DM_Rs)*( -C/(1.0+C) + log(1.0+C) )  );

// (3) isothermal beta model
   const double x = ClusterMerger_Rvir / ClusterMerger_Gas_Rcore;
   ClusterMerger_Gas_M     = ClusterMerger_Mvir * ClusterMerger_Gas_MFrac;
   ClusterMerger_Gas_Rho0  = ClusterMerger_Gas_M / (  4.0*M_PI*CUBE(ClusterMerger_Gas_Rcore)*( asinh(x) - x/sqrt(1.0+x*x) )  );


// initialize GSL
   GSL_WorkSpace = gsl_integration_workspace_alloc( GSL_WorkSize );


// set up interpolation tables
// (1) allocate memory
   ClusterMerger_Gas_PresProf   = new double [ ClusterMerger_NBin_PresProf  ];
   ClusterMerger_Gas_PresProf_R = new double [ ClusterMerger_NBin_PresProf  ];
   ClusterMerger_DM_MassProf    = new double [ ClusterMerger_NBin_MassProf  ];
   ClusterMerger_DM_MassProf_R  = new double [ ClusterMerger_NBin_MassProf  ];
   ClusterMerger_DM_SigmaProf   = new double [ ClusterMerger_NBin_SigmaProf ];
   ClusterMerger_DM_SigmaProf_R = new double [ ClusterMerger_NBin_SigmaProf ];

// (2) gas pressure
   SetTable_Gas_PresProf( ClusterMerger_NBin_PresProf, ClusterMerger_Gas_PresProf_R, ClusterMerger_Gas_PresProf );

// (3) dark matter mass profile
   SetTable_DM_MassProf( ClusterMerger_NBin_MassProf, ClusterMerger_DM_MassProf_R, ClusterMerger_DM_MassProf );

// (4) dark matter velocity dispersion
   SetTable_DM_SigmaProf( ClusterMerger_NBin_SigmaProf, ClusterMerger_DM_SigmaProf_R, ClusterMerger_DM_SigmaProf );


// record the test problem parameters
   if ( MPI_Rank == 0 )
   {
      Aux_Message( stdout, "\n" );
      Aux_Message( stdout, "%s test :\n", TestProb );
      Aux_Message( stdout, "=============================================================================\n" );
      Aux_Message( stdout, "   random seed for initializing particles             = %d\n",            ClusterMerger_RanSeed );
      Aux_Message( stdout, "   [total]       virial mass                          = %13.7e Msun\n",   ClusterMerger_Mvir*UNIT_M/Const_Msun );
      Aux_Message( stdout, "   [total]       virial radius (Rvir)                 = %13.7e Mpc\n",    ClusterMerger_Rvir*UNIT_L/Const_Mpc );
      Aux_Message( stdout, "   [total]       cut-off radius (Rcut)                = %13.7e Mpc\n",    ClusterMerger_Rcut*UNIT_L/Const_Mpc );
      Aux_Message( stdout, "   [total]       enclosed mass within Rcut            = %13.7e Msun\n",   MassProf_Total(ClusterMerger_Rcut)*UNIT_M/Const_Msun );
      Aux_Message( stdout, "   [dark matter] enclosed mass within Rvir            = %13.7e Msun\n",   ClusterMerger_DM_M*UNIT_M/Const_Msun );
      Aux_Message( stdout, "   [dark matter] enclosed mass within Rcut            = %13.7e Msun\n",   MassProf_DM(ClusterMerger_Rcut)*UNIT_M/Const_Msun );
      Aux_Message( stdout, "   [dark matter] NFW concentration parameter          = %13.7e\n",        ClusterMerger_DM_C );
      Aux_Message( stdout, "   [dark matter] NFW scale radius                     = %13.7e Mpc\n",    ClusterMerger_DM_Rs*UNIT_L/Const_Mpc );
      Aux_Message( stdout, "   [dark matter] NFW density parameter                = %13.7e g/cm^3\n", ClusterMerger_DM_Rho0*UNIT_D );
      Aux_Message( stdout, "   [dark matter] radius with zero velocity dispersion = %13.7e Mpc\n",    ClusterMerger_DM_Rzero*UNIT_L/Const_Mpc );
      Aux_Message( stdout, "   [gas]         enclosed mass within Rvir            = %13.7e Msun\n",   ClusterMerger_Gas_M*UNIT_M/Const_Msun );
      Aux_Message( stdout, "   [gas]         enclosed mass within Rcut            = %13.7e Msun\n",   MassProf_Gas(ClusterMerger_Rcut)*UNIT_M/Const_Msun );
      Aux_Message( stdout, "   [gas]         core radius in beta model            = %13.7e Mpc\n",    ClusterMerger_Gas_Rcore*UNIT_L/Const_Mpc );
      Aux_Message( stdout, "   [gas]         density parameter in beta model      = %13.7e g/cm^3\n", ClusterMerger_Gas_Rho0*UNIT_D );
      Aux_Message( stdout, "   [gas]         temperature normalization            = %d\n",            ClusterMerger_Gas_NormTvir );
      if ( ClusterMerger_Gas_NormTvir ) {
      Aux_Message( stdout, "   [gas]         temperature normalization at Rvir    = %13.7e keV\n",    ClusterMerger_Gas_Tvir*UNIT_E/Const_keV );
      Aux_Message( stdout, "                                                      = %13.7e K\n",      ClusterMerger_Gas_Tvir*UNIT_E/Const_kB );
      Aux_Message( stdout, "   [gas]         M-T scaling index                    = %13.7e\n",        ClusterMerger_Gas_MTScaling ); }
      Aux_Message( stdout, "   [gas]         mass fraction                        = %13.7e\n",        ClusterMerger_Gas_MFrac );
      Aux_Message( stdout, "   [gas]         density  ratio: background / Rcut    = %13.7e\n",        ClusterMerger_Gas_DensBgRatio );
      Aux_Message( stdout, "   [gas]         pressure ratio: background / Rcut    = %13.7e\n",        ClusterMerger_Gas_PresBgRatio );
      Aux_Message( stdout, "   # of bins in the gas pressure profile              = %d\n",            ClusterMerger_NBin_PresProf );
      Aux_Message( stdout, "   # of bins in the dark matter mass profile          = %d\n",            ClusterMerger_NBin_MassProf );
      Aux_Message( stdout, "   # of bins in the dark matter sigma profile         = %d\n",            ClusterMerger_NBin_SigmaProf );
      Aux_Message( stdout, "   minimum radius in the interpolation table          = %13.7e Mpc\n",    ClusterMerger_IntRmin*UNIT_L/Const_Mpc );
      Aux_Message( stdout, "   test mode                                          = %s\n",            (ClusterMerger_Coll)?
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


// free GSL resource
   gsl_integration_workspace_free( GSL_WorkSpace );

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

   const double BoxCenter[3] = { 0.5*amr->BoxSize[0], 0.5*amr->BoxSize[1], 0.5*amr->BoxSize[2] };
   const double DensBg       = ClusterMerger_Gas_DensBgRatio*DensProf_Gas( ClusterMerger_Rcut );
   const double PresBg       = ClusterMerger_Gas_PresBgRatio*ClusterMerger_Gas_PresProf[ ClusterMerger_NBin_PresProf - 1 ];

   if ( ClusterMerger_Coll )
   {
      Aux_Error( ERROR_INFO, "NOT SUPPORETD YET !!\n" );

      /*
      const double Coll_Offset = 0.5*ClusterMerger_Coll_D/sqrt(3.0);
      double Center[3];

      fluid[DENS] = 0.0;
      fluid[ENGY] = 0.0;

      for (int t=-1; t<=1; t+=2)
      {
         for (int d=0; d<3; d++)    Center[d] = ClusterMerger_Center[d] + Coll_Offset*(double)t;

         r2 = SQR(x-Center[0])+ SQR(y-Center[1]) + SQR(z-Center[2]);
         a2 = r2 / SQR(ClusterMerger_R0);

         fluid[DENS] += GasRho0 * pow( 1.0 + a2, -2.5 );
         fluid[ENGY] += (  NEWTON_G*TotM*GasRho0 / ( 6.0*ClusterMerger_R0*CUBE(1.0 + a2) ) + PresBg  ) / ( GAMMA - 1.0 );
      }

      fluid[MOMX]  = fluid[DENS]*ClusterMerger_BulkVel[0];
      fluid[MOMY]  = fluid[DENS]*ClusterMerger_BulkVel[1];
      fluid[MOMZ]  = fluid[DENS]*ClusterMerger_BulkVel[2];
      fluid[ENGY] += 0.5*( SQR(fluid[MOMX]) + SQR(fluid[MOMY]) + SQR(fluid[MOMZ]) ) / fluid[DENS];
      */
   }

   else
   {
      double r, Dens, Pres;

      r = sqrt( SQR(x-BoxCenter[0]) + SQR(y-BoxCenter[1]) + SQR(z-BoxCenter[2]) );

      if ( r >= ClusterMerger_Rcut )
      {
         Dens = DensBg;
         Pres = PresBg;
      }

      else if ( r < ClusterMerger_IntRmin )
         Aux_Error( ERROR_INFO, "r (%20.14e) < minimum radius in the interpolation table (%20.14e) !!\n",
                    r, ClusterMerger_IntRmin );

      else
      {
         Dens = DensProf_Gas( r );
         Pres = Mis_InterpolateFromTable( ClusterMerger_NBin_PresProf, ClusterMerger_Gas_PresProf_R,
                                          ClusterMerger_Gas_PresProf, r );

         if ( Pres == NULL_REAL )   Aux_Error( ERROR_INFO, "interpolation on pressure failed !!\n" );
      }

      fluid[DENS] = Dens;
      fluid[MOMX] = 0.0;
      fluid[MOMY] = 0.0;
      fluid[MOMZ] = 0.0;
      fluid[ENGY] = Pres / ( GAMMA - 1.0 )
                    + 0.5*( SQR(fluid[MOMX]) + SQR(fluid[MOMY]) + SQR(fluid[MOMZ]) ) / fluid[DENS];

   } // if ( ClusterMerger_Coll ) ... else ...

} // FUNCTION : Par_TestProbSol_ClusterMerger



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
   sscanf( input_line, "%d%s",   &ClusterMerger_RanSeed,         string );

   getline( &input_line, &len, File );
   sscanf( input_line, "%lf%s",  &ClusterMerger_Mvir,            string );

   getline( &input_line, &len, File );
   sscanf( input_line, "%lf%s",  &ClusterMerger_Rvir,            string );

   getline( &input_line, &len, File );
   sscanf( input_line, "%lf%s",  &ClusterMerger_Rcut,            string );

   getline( &input_line, &len, File );
   sscanf( input_line, "%lf%s",  &ClusterMerger_DM_C,            string );

   getline( &input_line, &len, File );
   sscanf( input_line, "%lf%s",  &ClusterMerger_DM_Rzero,        string );

   getline( &input_line, &len, File );
   sscanf( input_line, "%lf%s",  &ClusterMerger_Gas_Rcore,       string );

   getline( &input_line, &len, File );
   sscanf( input_line, "%d%s",   &temp_int,                      string );
   ClusterMerger_Gas_NormTvir = (bool)temp_int;

   getline( &input_line, &len, File );
   sscanf( input_line, "%lf%s",  &ClusterMerger_Gas_TvirM14,     string );

   getline( &input_line, &len, File );
   sscanf( input_line, "%lf%s",  &ClusterMerger_Gas_MTScaling,   string );

   getline( &input_line, &len, File );
   sscanf( input_line, "%lf%s",  &ClusterMerger_Gas_MFrac,       string );

   getline( &input_line, &len, File );
   sscanf( input_line, "%lf%s",  &ClusterMerger_Gas_DensBgRatio, string );

   getline( &input_line, &len, File );
   sscanf( input_line, "%lf%s",  &ClusterMerger_Gas_PresBgRatio, string );

   getline( &input_line, &len, File );
   sscanf( input_line, "%d%s",   &ClusterMerger_NBin_PresProf,   string );

   getline( &input_line, &len, File );
   sscanf( input_line, "%d%s",   &ClusterMerger_NBin_MassProf,   string );

   getline( &input_line, &len, File );
   sscanf( input_line, "%d%s",   &ClusterMerger_NBin_SigmaProf,  string );

   getline( &input_line, &len, File );
   sscanf( input_line, "%lf%s",  &ClusterMerger_IntRmin,         string );

   getline( &input_line, &len, File );
   sscanf( input_line, "%d%s",   &temp_int,                      string );
   ClusterMerger_Coll = (bool)temp_int;

   fclose( File );
   if ( input_line != NULL )     free( input_line );


// set the default parameters
   if ( ClusterMerger_Rcut <= 0.0 )
   {
      for (int d=0; d<3; d++)    ClusterMerger_Rcut = fmax( ClusterMerger_Rcut, 0.5*amr->BoxSize[d]*sqrt(3.0) );

      if ( MPI_Rank == 0 )  Aux_Message( stdout, "NOTE : parameter \"%s\" is set to the default value = %13.7e\n",
                                         "ClusterMerger_Rcut", ClusterMerger_Rcut );
   }

   if ( ClusterMerger_DM_Rzero <= 0.0 )
   {
      ClusterMerger_DM_Rzero = 1.0e4*ClusterMerger_Rcut;

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

   if ( ClusterMerger_Gas_DensBgRatio <= 0.0 )
   {
      ClusterMerger_Gas_DensBgRatio = 1.0e-2;

      if ( MPI_Rank == 0 )  Aux_Message( stdout, "NOTE : parameter \"%s\" is set to the default value = %13.7e\n",
                                         "ClusterMerger_Gas_DensBgRatio", ClusterMerger_Gas_DensBgRatio );
   }

   if ( ClusterMerger_Gas_PresBgRatio <= 0.0 )
   {
      ClusterMerger_Gas_PresBgRatio = 1.0e-2;

      if ( MPI_Rank == 0 )  Aux_Message( stdout, "NOTE : parameter \"%s\" is set to the default value = %13.7e\n",
                                         "ClusterMerger_Gas_PresBgRatio", ClusterMerger_Gas_PresBgRatio );
   }

   if ( ClusterMerger_NBin_PresProf <= 1 )
   {
      ClusterMerger_NBin_PresProf = 10000;

      if ( MPI_Rank == 0 )  Aux_Message( stdout, "NOTE : parameter \"%s\" is set to the default value = %d\n",
                                         "ClusterMerger_NBin_PresProf", ClusterMerger_NBin_PresProf );
   }

   if ( ClusterMerger_NBin_MassProf <= 1 )
   {
      ClusterMerger_NBin_MassProf = 10000;

      if ( MPI_Rank == 0 )  Aux_Message( stdout, "NOTE : parameter \"%s\" is set to the default value = %d\n",
                                         "ClusterMerger_NBin_MassProf", ClusterMerger_NBin_MassProf );
   }

   if ( ClusterMerger_NBin_SigmaProf <= 1 )
   {
      ClusterMerger_NBin_SigmaProf = 10000;

      if ( MPI_Rank == 0 )  Aux_Message( stdout, "NOTE : parameter \"%s\" is set to the default value = %d\n",
                                         "ClusterMerger_NBin_SigmaProf", ClusterMerger_NBin_SigmaProf );
   }

   if ( ClusterMerger_IntRmin <= 0.0 )
   {
      ClusterMerger_IntRmin = 1.0e-5;     // in Mpc

      if ( MPI_Rank == 0 )  Aux_Message( stdout, "NOTE : parameter \"%s\" is set to the default value = %13.7e Mpc\n",
                                         "ClusterMerger_IntRmin", ClusterMerger_IntRmin );
   }


// check
   if ( ClusterMerger_RanSeed < 0 )
      Aux_Error( ERROR_INFO, "ClusterMerger_RanSeed (%d) < 0 !!\n", ClusterMerger_RanSeed );

   if ( ClusterMerger_Mvir <= 0.0 )
      Aux_Error( ERROR_INFO, "ClusterMerger_Mvir (%14.7e) <= 0.0 !!\n", ClusterMerger_Mvir );

   if ( ClusterMerger_Rvir <= 0.0  ||  ClusterMerger_Rvir > ClusterMerger_Rcut )
      Aux_Error( ERROR_INFO, "incorrect ClusterMerger_Rvir = %13.7e (ClusterMerger_Rcut = %13.7e) !!\n",
                 ClusterMerger_Rvir, ClusterMerger_Rcut );

   if ( ClusterMerger_Rcut <= 0.0 )
      Aux_Error( ERROR_INFO, "ClusterMerger_Rcut (%14.7e) <= 0.0 !!\n", ClusterMerger_Rcut );

   if ( ClusterMerger_DM_C <= 0.0 )
      Aux_Error( ERROR_INFO, "ClusterMerger_DM_C (%14.7e) <= 0.0 !!\n", ClusterMerger_DM_C );

   if ( ClusterMerger_DM_Rzero < ClusterMerger_Rcut )
      Aux_Error( ERROR_INFO, "ClusterMerger_DM_Rzero (%13.7e) < ClusterMerger_Rcut (%13.7e) !!\n",
                 ClusterMerger_DM_Rzero, ClusterMerger_Rcut );

   if ( ClusterMerger_Gas_Rcore <= 0.0  ||  ClusterMerger_Gas_Rcore > ClusterMerger_Rcut )
      Aux_Error( ERROR_INFO, "incorrect ClusterMerger_Gas_Rcore = %13.7e (ClusterMerger_Rcut = %13.7e) !!\n",
                 ClusterMerger_Gas_Rcore, ClusterMerger_Rcut );

   if ( ClusterMerger_Gas_TvirM14 <= 0.0 )
      Aux_Error( ERROR_INFO, "ClusterMerger_Gas_TvirM14 (%14.7e) <= 0.0 !!\n", ClusterMerger_Gas_TvirM14 );

   if ( ClusterMerger_Gas_MTScaling <= 0.0 )
      Aux_Error( ERROR_INFO, "ClusterMerger_Gas_MTScaling (%14.7e) <= 0.0 !!\n", ClusterMerger_Gas_MTScaling );

   if ( ClusterMerger_Gas_MFrac <= 0.0  ||  ClusterMerger_Gas_MFrac > 1.0 )
      Aux_Error( ERROR_INFO, "ClusterMerger_Gas_MFrac (%13.7e) is not within the correct range [0.0 < x <= 1.0] !!\n",
                 ClusterMerger_Gas_MFrac );

   if ( ClusterMerger_Gas_DensBgRatio <= 0.0 )
      Aux_Error( ERROR_INFO, "ClusterMerger_Gas_DensBgRatio (%14.7e) <= 0.0 !!\n", ClusterMerger_Gas_DensBgRatio );

   if ( ClusterMerger_Gas_PresBgRatio <= 0.0 )
      Aux_Error( ERROR_INFO, "ClusterMerger_Gas_PresBgRatio (%14.7e) <= 0.0 !!\n", ClusterMerger_Gas_PresBgRatio );

   if ( ClusterMerger_NBin_PresProf <= 1   )
      Aux_Error( ERROR_INFO, "ClusterMerger_NBin_PresProf (%d) <= 1 !!\n", ClusterMerger_NBin_PresProf );

   if ( ClusterMerger_NBin_MassProf <= 1   )
      Aux_Error( ERROR_INFO, "ClusterMerger_NBin_MassProf (%d) <= 1 !!\n", ClusterMerger_NBin_MassProf );

   if ( ClusterMerger_NBin_SigmaProf <= 1   )
      Aux_Error( ERROR_INFO, "ClusterMerger_NBin_SigmaProf (%d) <= 1 !!\n", ClusterMerger_NBin_SigmaProf );

   if ( ClusterMerger_IntRmin <= 0.0 )
      Aux_Error( ERROR_INFO, "ClusterMerger_IntRmin (%14.7e) <= 0.0 !!\n", ClusterMerger_IntRmin );

} // FUNCTION : LoadTestProbParameter



//-------------------------------------------------------------------------------------------------------
// Function    :  MassProf_Total
// Description :  Total mass profile (dark matter + gas)
//
// Note        :  Calculate the enclosed mass for a given radius
//
// Parameter   :  r  : Input radius
//
// Return      :  Total enclosed mass
//-------------------------------------------------------------------------------------------------------
double MassProf_Total( const double r )
{

   return MassProf_DM( r ) + MassProf_Gas( r );

} // FUNCTION : MassProf_Total



//-------------------------------------------------------------------------------------------------------
// Function    :  MassProf_DM
// Description :  Mass profile of the NFW model for dark matter
//
// Note        :  Calculate the enclosed mass for a given radius
//
// Parameter   :  r  : Input radius
//
// Return      :  Dark matter enclosed mass
//-------------------------------------------------------------------------------------------------------
double MassProf_DM( const double r )
{

   static double Norm = 4.0*M_PI*CUBE(ClusterMerger_DM_Rs)*ClusterMerger_DM_Rho0;

   const double x = r / ClusterMerger_DM_Rs;

   return Norm*( -x/(1.0+x) + log(1.0+x) );

} // FUNCTION : MassProf_DM



//-------------------------------------------------------------------------------------------------------
// Function    :  MassProf_Gas
// Description :  Mass profile of the isothermal beta model for gas
//
// Note        :  1. Calculate the enclosed mass for a given radius
//                2. Currently only work with beta = 1.0
//
// Parameter   :  r  : Input radius
//
// Return      :  Gas enclosed mass
//-------------------------------------------------------------------------------------------------------
double MassProf_Gas( const double r )
{

   static double Norm = 4.0*M_PI*CUBE(ClusterMerger_Gas_Rcore)*ClusterMerger_Gas_Rho0;

   const double x = r / ClusterMerger_Gas_Rcore;

   return Norm*( asinh(x) - x/sqrt(1.0+x*x) );

} // FUNCTION : MassProf_Gas



//-------------------------------------------------------------------------------------------------------
// Function    :  DensProf_DM
// Description :  Density profile of the NFW model for dark matter
//
// Note        :  Calculate the mass density for a given radius
//
// Parameter   :  r  : Input radius
//
// Return      :  Dark matter mass density
//-------------------------------------------------------------------------------------------------------
double DensProf_DM( const double r )
{

   const double x = r / ClusterMerger_DM_Rs;

   return ClusterMerger_DM_Rho0 / ( x*SQR(1.0+x) );

} // FUNCTION : DensProf_DM



//-------------------------------------------------------------------------------------------------------
// Function    :  DensProf_Gas
// Description :  Density profile of the isothermal beta model for gas
//
// Note        :  1. Calculate the mass density for a given radius
//                2. Currently only work with beta = 1.0
//
// Parameter   :  r  : Input radius
//
// Return      :  Gas mass density
//-------------------------------------------------------------------------------------------------------
double DensProf_Gas( const double r )
{

   const double x    = r / ClusterMerger_Gas_Rcore;
   const double beta = 1.0;

   return ClusterMerger_Gas_Rho0 * pow( 1.0+x*x, -1.5*beta );

} // FUNCTION : DensProf_Gas



//-------------------------------------------------------------------------------------------------------
// Function    :  SetTable_Gas_PresProf
// Description :  Calculate the interpolation table of gas pressure profile
//
// Note        :  1. Calculate pressure from hydrostatic equilibrium
//                   --> grad( P(r) ) = -rho(r)*grad( phi(r) ) = -G*rho(r)*M_tot(r)/r^2
//                   --> P(r1) = P(r2) + G*Integrate( rho(r)*M_tot(r)/r^2, [r, r1, r2] )
//                   --> Note that rho is gas density and M_tot is **total** mass profile
//                2. P(Rvir) is fixed by T(Rvir)/rho(Rvir), where T(Rvir) is ClusterMerger_Gas_Tvir,
//                   if ClusterMerger_Gas_NormTvir is on
//                   --> Otherwise it's normalized to P(Rzero) = 0
//                3. Evenly sample in log space
//                   --> Table_R[     0] = IntRmin
//                       Table_R[NBin-1] = Rcut
//                4. r and Pres arrays must be preallocated with the size NBin
//
// Parameter   :  NBin  : Number of radial bins in the table
//                r     : Radius at each bin
//                Pres  : Gas pressure at each bin
//
// Return      :  r, Pres
//-------------------------------------------------------------------------------------------------------
void SetTable_Gas_PresProf( const int NBin, double *r, double *Pres )
{

// check
   if ( r    == NULL )  Aux_Error( ERROR_INFO, "r == NULL !!\n" );
   if ( Pres == NULL )  Aux_Error( ERROR_INFO, "Pres == NULL !!\n" );


// set up GSL
   const double GSL_MaxAbsErr = 0.0;      // maximum allowed absolute error (0 --> disable)
   const double GSL_MaxRelErr = 1.0e-6;   // maximum allowed relative error (0 --> disable)
   gsl_function GSL_Func;

   GSL_Func.function = &GSL_IntFunc_Gas_PresProf;
   GSL_Func.params   = NULL;


// set up radius
// --> r_min should be much smaller than dh_min in case that cell center and cluster center are very close
//     (possible for the cluster merger case where cluster centers can be arbitrarily close to a cell center)
   const double r_min = ClusterMerger_IntRmin;
   const double r_max = ClusterMerger_Rcut;
   const double dr    = pow( r_max/r_min, 1.0/(NBin-1.0) );

   for (int b=0; b<NBin; b++)    r[b] = r_min*pow( dr, (double)b );


// calculate pressure at Rcut by having P(Rzero) = 0
   double GSL_Result, GSL_AbsErr;
   gsl_integration_qag( &GSL_Func, r[NBin-1], ClusterMerger_DM_Rzero, GSL_MaxAbsErr, GSL_MaxRelErr, GSL_WorkSize,
                        GSL_IntRule, GSL_WorkSpace, &GSL_Result, &GSL_AbsErr );

// gravitational constant is not included in the GSL integrand for better performance
   Pres[NBin-1] = NEWTON_G*GSL_Result;


// calculate pressure at r < Rcut
   for (int b=NBin-2; b>=0; b--)
   {
      gsl_integration_qag( &GSL_Func, r[b], r[b+1], GSL_MaxAbsErr, GSL_MaxRelErr, GSL_WorkSize,
                           GSL_IntRule, GSL_WorkSpace, &GSL_Result, &GSL_AbsErr );

      Pres[b] = Pres[b+1] + NEWTON_G*GSL_Result;
   }


// renormalize pressure so that T(Rvir) == ClusterMerger_Gas_Tvir
// --> not recommended ...
   if ( ClusterMerger_Gas_NormTvir )
   {
      const double PresNorm = DensProf_Gas( ClusterMerger_Rvir ) * ClusterMerger_Gas_Tvir / ( Const_amu/UNIT_M*MOLECULAR_WEIGHT );
      const int    Idx      = Mis_BinarySearch_Real( r, 0, NBin-1, ClusterMerger_Rvir );

      if ( Idx < 0 )    Aux_Error( ERROR_INFO, "cannot determine the array index with r ~ Rvir !!\n" );

//    don't bother applying interpolation since it doesn't need to be very precise here
      const double PresDiff = PresNorm - Pres[Idx];

      for (int b=0; b<NBin; b++)    Pres[b] += PresDiff;
   }

} // FUNCTION : SetTable_Gas_PresProf



//-------------------------------------------------------------------------------------------------------
// Function    :  GSL_IntFunc_Gas_PresProf
// Description :  GSL integrand for calculating the gas pressure profile
//
// Note        :  1. Return rho_gas(r)*M_total(r)/r^2
//                2. Gavitational constant is not included in this integrand and must be multiplied later
//
// Parameter   :  r        : Radius
//                IntPara  : Integration parameters
//
// Return      :  See note above
//-------------------------------------------------------------------------------------------------------
double GSL_IntFunc_Gas_PresProf( double r, void *IntPara )
{

   return DensProf_Gas( r ) * MassProf_Total( r ) / SQR(r);

} // FUNCTION : GSL_IntFunc_Gas_PresProf



//-------------------------------------------------------------------------------------------------------
// Function    :  SetTable_DM_MassProf
// Description :  Calculate the interpolation table of dark matter mass profile
//
// Note        :  1. NFW has analytical mass profile
//                   --> This table is for estimating r(M) instead of M(r)
//                3. Evenly sample in log space
//                   --> Table_R[     0] = IntRmin
//                       Table_R[NBin-1] = Rcut
//                4. r and Mass arrays must be preallocated with the size NBin
//
// Parameter   :  NBin  : Number of radial bins in the table
//                r     : Radius at each bin
//                Mass  : Dark matter mass profile at each bin
//
// Return      :  r, Mass
//-------------------------------------------------------------------------------------------------------
void SetTable_DM_MassProf( const int NBin, double *r, double *Mass )
{

// check
   if ( r    == NULL )  Aux_Error( ERROR_INFO, "r == NULL !!\n" );
   if ( Mass == NULL )  Aux_Error( ERROR_INFO, "Mass == NULL !!\n" );


// set up radius
// --> r_min should be much smaller than dh_min since the particle position can be arbitrarily close to the cluster center
   const double r_min = ClusterMerger_IntRmin;
   const double r_max = ClusterMerger_Rcut;
   const double dr    = pow( r_max/r_min, 1.0/(NBin-1.0) );

   for (int b=0; b<NBin; b++)    r[b] = r_min*pow( dr, (double)b );


// get mass profile at each radial bin
   for (int b=0; b<NBin; b++)    Mass[b] = MassProf_DM( r[b] );

} // FUNCTION : SetTable_DM_MassProf



//-------------------------------------------------------------------------------------------------------
// Function    :  SetTable_DM_SigmaProf
// Description :  Calculate the interpolation table of dark matter velocity dispersion profile
//
// Note        :  1. Calculate velocity dispersion from Jeans equation
//                   --> grad( rho(r)*sigma(r)^2 ) = -rho(r)*grad( phi(r) ) = -G*rho(r)*M_tot(r)/r^2
//                   --> rho(r1)*sigma(r1)^2 = rho(r2)*sigma(r2)^2 + G*Integrate( rho(r)*M_tot(r)/r^2, [r, r1, r2] )
//                   --> Note that rho is dark matter density and M_tot is **total** mass profile
//                   --> It's essentially the same as calculating gas pressure
//                       --> Just replace pressure_gas(r) by rho_dm(r)*sigma(r)^2
//                2. Normalized to rho(Rzero)*sigma(Rzero)^2 = 0
//                3. Evenly sample in log space
//                   --> Table_R[     0] = IntRmin
//                       Table_R[NBin-1] = Rcut
//                4. r and Sigma arrays must be preallocated with the size NBin
//
// Parameter   :  NBin  : Number of radial bins in the table
//                r     : Radius at each bin
//                Sigma : Dark matter velocity dispersion at each bin
//
// Return      :  r, Sigma
//-------------------------------------------------------------------------------------------------------
void SetTable_DM_SigmaProf( const int NBin, double *r, double *Sigma )
{

// check
   if ( r     == NULL )    Aux_Error( ERROR_INFO, "r == NULL !!\n" );
   if ( Sigma == NULL )    Aux_Error( ERROR_INFO, "Sigma == NULL !!\n" );


// set up GSL
   const double GSL_MaxAbsErr = 0.0;      // maximum allowed absolute error (0 --> disable)
   const double GSL_MaxRelErr = 1.0e-6;   // maximum allowed relative error (0 --> disable)
   gsl_function GSL_Func;

   GSL_Func.function = &GSL_IntFunc_DM_SigmaProf;
   GSL_Func.params   = NULL;


// set up radius
// --> r_min should be much smaller than dh_min since the particle position can be arbitrarily close to the cluster center
   const double r_min = ClusterMerger_IntRmin;
   const double r_max = ClusterMerger_Rcut;
   const double dr    = pow( r_max/r_min, 1.0/(NBin-1.0) );

   for (int b=0; b<NBin; b++)    r[b] = r_min*pow( dr, (double)b );


// calculate rho*sigma^2 at Rcut by having rho(Rzero)*sigma(Rzero)^2 = 0
   double GSL_Result, GSL_AbsErr;
   gsl_integration_qag( &GSL_Func, r[NBin-1], ClusterMerger_DM_Rzero, GSL_MaxAbsErr, GSL_MaxRelErr, GSL_WorkSize,
                        GSL_IntRule, GSL_WorkSpace, &GSL_Result, &GSL_AbsErr );

// gravitational constant is not included in the GSL integrand for better performance
   Sigma[NBin-1] = NEWTON_G*GSL_Result;


// calculate rho*sigma^2 at r < Rcut
   for (int b=NBin-2; b>=0; b--)
   {
      gsl_integration_qag( &GSL_Func, r[b], r[b+1], GSL_MaxAbsErr, GSL_MaxRelErr, GSL_WorkSize,
                           GSL_IntRule, GSL_WorkSpace, &GSL_Result, &GSL_AbsErr );

      Sigma[b] = Sigma[b+1] + NEWTON_G*GSL_Result;
   }


// convert rho*sigma^2 to sigma
   for (int b=0; b<NBin; b++)    Sigma[b] = sqrt(  Sigma[b] / DensProf_DM( r[b] )  );

} // FUNCTION : SetTable_DM_SigmaProf



//-------------------------------------------------------------------------------------------------------
// Function    :  GSL_IntFunc_DM_SigmaProf
// Description :  GSL integrand for calculating the dark matter velocity dispersion profile
//
// Note        :  1. Return rho_dm(r)*M_total(r)/r^2
//                2. Gavitational constant is not included in this integrand and must be multiplied later
//
// Parameter   :  r        : Radius
//                IntPara  : Integration parameters
//
// Return      :  See note above
//-------------------------------------------------------------------------------------------------------
double GSL_IntFunc_DM_SigmaProf( double r, void *IntPara )
{

   return DensProf_DM( r ) * MassProf_Total( r ) / SQR(r);

} // FUNCTION : GSL_IntFunc_DM_SigmaProf



#endif // #if ( MODEL == HYDRO  &&  defined PARTICLE )
