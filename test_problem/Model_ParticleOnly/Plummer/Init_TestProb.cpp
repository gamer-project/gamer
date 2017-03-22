#include "GAMER.h"

#ifdef PARTICLE



extern void (*Init_Function_Ptr)( real fluid[], const double x, const double y, const double z, const double Time );
extern void (*Output_TestProbErr_Ptr)( const bool BaseOnly );

static void LoadTestProbParameter();
static void Par_TestProbSol_Plummer( real fluid[], const double x, const double y, const double z, const double Time );

extern double MassProf_Plummer( const double r );


// global variables in the Plummer model test
// =======================================================================================
int  Plummer_RSeed;        // random seed for setting particle position and velocity
real Plummer_MaxR;         // maximum distance from the origin
real Plummer_Rho0;         // peak density
real Plummer_R0;           // scale radius
int  Plummer_NBinR;        // number of radial bins in the mass profile table
bool Plummer_Collision;    // true/false --> test Plummer collision/single Plummer
real Plummer_Collision_D;  // distance between two Plummer clouds for the Plummer collision test
real Plummer_Center[3];    // central coordinates
real Plummer_BulkVel[3];   // bulk velocity
#if ( MODEL == HYDRO )
real Plummer_GasMFrac;     // gas mass fraction (for HYDRO only)
#endif

real Plummer_FreeT;        // free-fall time at Plummer_R0
// =======================================================================================




//-------------------------------------------------------------------------------------------------------
// Function    :  Init_TestProb
// Description :  Initialize parameters for the Plummer model test
//
// Note        :  1. Please copy this file to "GAMER/src/Init/Init_TestProb.cpp"
//                2. Test problem parameters can be set in the input file "Input__TestProb"
//
// Parameter   :  None
//-------------------------------------------------------------------------------------------------------
void Init_TestProb()
{

   const char *TestProb = "Plummer model";

// check
# if ( MODEL != ELBDM  &&  MODEL != HYDRO )
# error : ERROR : only support ELBDM and HYDRO models !!
# endif

# ifndef PARTICLE
# error : ERROR : "PARTICLE is NOT defined" in the Plummer model test !!
# endif

# ifndef GRAVITY
# error : ERROR : "GRAVITY must be ON" in the Plummer model test !!
# endif

# ifdef COMOVING
# error : ERROR : "COMOVING must be OFF" in the Plummer model test !!
# endif


// set the initialization and output functions
   Init_Function_Ptr      = Par_TestProbSol_Plummer;
   Output_TestProbErr_Ptr = NULL;


// load the test problem parameters
   LoadTestProbParameter();


// set the test problem parameters
   Plummer_FreeT = sqrt( (3.0*M_PI*pow(2.0,1.5)) / (32.0*NEWTON_G*Plummer_Rho0) );


// record the test problem parameters
   if ( MPI_Rank == 0 )
   {
      Aux_Message( stdout, "\n" );
      Aux_Message( stdout, "%s test :\n", TestProb );
      Aux_Message( stdout, "=============================================================================\n" );
      Aux_Message( stdout, " Note: random seed for setting particle position    = %d\n",      Plummer_RSeed );
      Aux_Message( stdout, "       maximum radius                               = %13.7e\n",  Plummer_MaxR );
      Aux_Message( stdout, "       density parameter                            = %13.7e\n",  Plummer_Rho0 );
      Aux_Message( stdout, "       scale radius                                 = %13.7e\n",  Plummer_R0 );
      Aux_Message( stdout, "       number of radial bins in the mass profile    = %d\n",      Plummer_NBinR );
      Aux_Message( stdout, "       free fall time at scale radius               = %13.7e\n",  Plummer_FreeT );
      Aux_Message( stdout, "       test mode                                    = %s\n",     (Plummer_Collision)?
                                                                                             "collision":"single cloud" );
      for (int d=0; d<3; d++)
      Aux_Message( stdout, "       central coordinate [%d]                       = %14.7e\n", d, Plummer_Center[d] );
      if ( Plummer_Collision )
      Aux_Message( stdout, "       initial distance between two clouds          = %13.7e\n",  Plummer_Collision_D );
      for (int d=0; d<3; d++)
      Aux_Message( stdout, "       bulk velocity [%d]                            = %14.7e\n", d, Plummer_BulkVel[d] );
#     if ( MODEL == HYDRO )
      Aux_Message( stdout, "       gas mass fraction                            = %13.7e\n", Plummer_GasMFrac );
#     endif
      Aux_Message( stdout, "=============================================================================\n" );
      Aux_Message( stdout, "\n" );
   }


// set some default parameters
// End_T : 20 free-fall time at the scale radius
   const double End_T_Default    = (Plummer_Collision) ? 50.0 : 20.0*Plummer_FreeT;
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

   if ( OPT__OUTPUT_TEST_ERROR )
   {
      OPT__OUTPUT_TEST_ERROR = false;

      if ( MPI_Rank == 0 )
         Aux_Message( stdout, "NOTE : parameter %s is reset to %d in the %s test !!\n",
                      "OPT__OUTPUT_TEST_ERROR", OPT__OUTPUT_TEST_ERROR, TestProb );
   }

   if ( OPT__INIT == INIT_STARTOVER  &&  amr->Par->Init != PAR_INIT_BY_FUNCTION )
      Aux_Error( ERROR_INFO, "please set PAR_INIT = PAR_INIT_BY_FUNCTION !!\n" );

   if ( OPT__EXTERNAL_POT )
   {
      OPT__EXTERNAL_POT = false;

      if ( MPI_Rank == 0 )
         Aux_Message( stdout, "NOTE : option %s is not supported for the %s test and has been disabled !!\n",
                      "OPT__EXTERNAL_POT", TestProb );
   }

} // FUNCTION : Init_TestProb



//-------------------------------------------------------------------------------------------------------
// Function    :  Par_TestProbSol_Plummer
// Description :  Initialize the background density field as zero for the Plummer model test
//
// Note        :  1. This test supports two modes:
//                   (1) MODEL == ELBDM: grid density is set to zero. So essentially there are particles only
//                   (2) MODEL == HYDRO: gas and particles share the same density profile, and the mass ratio is
//                                       determined by the input parameter "Plummer_GasMFrac"
//
// Parameter   :  fluid : Fluid field to be initialized
//                x/y/z : Target physical coordinates
//                Time  : Target physical time
//
// Return      :  fluid
//-------------------------------------------------------------------------------------------------------
void Par_TestProbSol_Plummer( real *fluid, const double x, const double y, const double z, const double Time )
{

#  if   ( MODEL == HYDRO )
// gas share the same density profile as particles (except for different total masses)
   const double TotM    = 4.0/3.0*M_PI*CUBE(Plummer_R0)*Plummer_Rho0;
   const double GasRho0 = Plummer_Rho0*Plummer_GasMFrac;
   const double PresBg  = 0.0;   // background pressure
   double r2, a2, Dens;

   if ( Plummer_Collision )
   {
      const double Coll_Offset = 0.5*Plummer_Collision_D/sqrt(3.0);
      double Center[3];

      fluid[DENS] = 0.0;
      fluid[ENGY] = 0.0;
      for (int v=NCOMP_FLUID; v<NCOMP_TOTAL; v++)  fluid[v] = 0.0;

      for (int t=-1; t<=1; t+=2)
      {
         for (int d=0; d<3; d++)    Center[d] = Plummer_Center[d] + Coll_Offset*(double)t;

         r2   = SQR(x-Center[0]) + SQR(y-Center[1]) + SQR(z-Center[2]);
         a2   = r2 / SQR(Plummer_R0);
         Dens = GasRho0 * pow( 1.0 + a2, -2.5 );

         fluid[DENS] += Dens;
         fluid[ENGY] += (  NEWTON_G*TotM*GasRho0 / ( 6.0*Plummer_R0*CUBE(1.0 + a2) ) + PresBg  ) / ( GAMMA - 1.0 );

#        if ( NCOMP_PASSIVE > 0 )
#        if ( NCOMP_PASSIVE == 2 )
         fluid[ (t==-1)?CLOUD0:CLOUD1 ] = Dens;
#        else
#        error : ERROR : please specify how to set the passive scalars here !!
#        endif
#        endif // #if ( NCOMP_PASSIVE > 0 )
      }

      fluid[MOMX]  = fluid[DENS]*Plummer_BulkVel[0];
      fluid[MOMY]  = fluid[DENS]*Plummer_BulkVel[1];
      fluid[MOMZ]  = fluid[DENS]*Plummer_BulkVel[2];
      fluid[ENGY] += 0.5*( SQR(fluid[MOMX]) + SQR(fluid[MOMY]) + SQR(fluid[MOMZ]) ) / fluid[DENS];
   }

   else
   {
      r2   = SQR(x-Plummer_Center[0]) + SQR(y-Plummer_Center[1]) + SQR(z-Plummer_Center[2]);
      a2   = r2 / SQR(Plummer_R0);
      Dens = GasRho0 * pow( 1.0 + a2, -2.5 );

      fluid[DENS] = Dens;
      fluid[MOMX] = fluid[DENS]*Plummer_BulkVel[0];
      fluid[MOMY] = fluid[DENS]*Plummer_BulkVel[1];
      fluid[MOMZ] = fluid[DENS]*Plummer_BulkVel[2];
      fluid[ENGY] = (  NEWTON_G*TotM*GasRho0 / ( 6.0*Plummer_R0*CUBE(1.0 + a2) ) + PresBg  ) / ( GAMMA - 1.0 )
                    + 0.5*( SQR(fluid[MOMX]) + SQR(fluid[MOMY]) + SQR(fluid[MOMZ]) ) / fluid[DENS];

//    just set all passive scalars as zero
      for (int v=NCOMP_FLUID; v<NCOMP_TOTAL; v++)  fluid[v] = 0.0;
   } // if ( Plummer_Collision ) ... else ...

#  elif ( MODEL == ELBDM )
// set wave function as zero everywhere
   fluid[REAL] = 0.0;
   fluid[IMAG] = 0.0;
   fluid[DENS] = 0.0;

#  else
#  error : unsupported MODEL !!
#  endif

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


   getline( &input_line, &len, File );

   getline( &input_line, &len, File );
   sscanf( input_line, "%d%s",   &Plummer_RSeed,          string );

#  ifdef FLOAT8
   getline( &input_line, &len, File );
   sscanf( input_line, "%lf%s",  &Plummer_MaxR,           string );

   getline( &input_line, &len, File );
   sscanf( input_line, "%lf%s",  &Plummer_Rho0,           string );

   getline( &input_line, &len, File );
   sscanf( input_line, "%lf%s",  &Plummer_R0,             string );
#  else
   getline( &input_line, &len, File );
   sscanf( input_line, "%f%s",   &Plummer_MaxR,           string );

   getline( &input_line, &len, File );
   sscanf( input_line, "%f%s",  &Plummer_Rho0,            string );

   getline( &input_line, &len, File );
   sscanf( input_line, "%f%s",  &Plummer_R0,              string );
#  endif // #ifdef FLOAT8 ... else ...

   getline( &input_line, &len, File );
   sscanf( input_line, "%d%s",   &Plummer_NBinR,          string );

   getline( &input_line, &len, File );
   sscanf( input_line, "%d%s",   &temp_int,               string );
   Plummer_Collision = (bool)temp_int;

#  ifdef FLOAT8
   getline( &input_line, &len, File );
   sscanf( input_line, "%lf%s",  &Plummer_Collision_D,    string );

   getline( &input_line, &len, File );
   sscanf( input_line, "%lf%s",  &Plummer_Center[0],      string );

   getline( &input_line, &len, File );
   sscanf( input_line, "%lf%s",  &Plummer_Center[1],      string );

   getline( &input_line, &len, File );
   sscanf( input_line, "%lf%s",  &Plummer_Center[2],      string );

   getline( &input_line, &len, File );
   sscanf( input_line, "%lf%s",  &Plummer_BulkVel[0],     string );

   getline( &input_line, &len, File );
   sscanf( input_line, "%lf%s",  &Plummer_BulkVel[1],     string );

   getline( &input_line, &len, File );
   sscanf( input_line, "%lf%s",  &Plummer_BulkVel[2],     string );

   getline( &input_line, &len, File );
#  if ( MODEL == HYDRO )
   sscanf( input_line, "%lf%s",  &Plummer_GasMFrac,       string );
#  endif

#  else // #ifdef FLOAT8

   getline( &input_line, &len, File );
   sscanf( input_line, "%f%s",   &Plummer_Collision_D,    string );

   getline( &input_line, &len, File );
   sscanf( input_line, "%f%s",   &Plummer_Center[0],      string );

   getline( &input_line, &len, File );
   sscanf( input_line, "%f%s",   &Plummer_Center[1],      string );

   getline( &input_line, &len, File );
   sscanf( input_line, "%f%s",   &Plummer_Center[2],      string );

   getline( &input_line, &len, File );
   sscanf( input_line, "%f%s",   &Plummer_BulkVel[0],     string );

   getline( &input_line, &len, File );
   sscanf( input_line, "%f%s",   &Plummer_BulkVel[1],     string );

   getline( &input_line, &len, File );
   sscanf( input_line, "%f%s",   &Plummer_BulkVel[2],     string );

   getline( &input_line, &len, File );
#  if ( MODEL == HYDRO )
   sscanf( input_line, "%f%s",  &Plummer_GasMFrac,        string );
#  endif

#  endif // #ifdef FLOAT8 ... else ...

   fclose( File );
   if ( input_line != NULL )     free( input_line );


// set the default test problem parameters
   if ( Plummer_MaxR <= 0.0 )
   {
//    since the outmost base patches are not allowed to be refined for the isolated gravity solver,
//    we set MaxR = 0.5*L-2*PatchSize to ensure that all particles can be refined
//    Plummer_MaxR = 0.5*amr->BoxSize[0] - 2.0*PS1*amr->dh[0];
      Plummer_MaxR = 0.375;

      if ( MPI_Rank == 0 )  Aux_Message( stdout, "NOTE : parameter \"%s\" is set to the default value = %13.7e\n",
                                         "Plummer_MaxR", Plummer_MaxR );
   }

   if ( Plummer_Rho0 <= 0.0 )
   {
      Plummer_Rho0 = 1.0;

      if ( MPI_Rank == 0 )  Aux_Message( stdout, "NOTE : parameter \"%s\" is set to the default value = %13.7e\n",
                                         "Plummer_Rho0", Plummer_Rho0 );
   }

   if ( Plummer_R0 <= 0.0 )
   {
      Plummer_R0 = 0.1;

      if ( MPI_Rank == 0 )  Aux_Message( stdout, "NOTE : parameter \"%s\" is set to the default value = %13.7e\n",
                                         "Plummer_R0", Plummer_R0 );
   }

   if ( Plummer_NBinR <= 1 )
   {
      Plummer_NBinR = 1000;

      if ( MPI_Rank == 0 )  Aux_Message( stdout, "NOTE : parameter \"%s\" is set to the default value = %d\n",
                                         "Plummer_NBinR", Plummer_NBinR );
   }

   if ( Plummer_Collision  &&  Plummer_Collision_D < 0.0 )
   {
      Plummer_Collision_D = 1.5;

      if ( MPI_Rank == 0 )  Aux_Message( stdout, "NOTE : parameter \"%s\" is set to the default value = %13.7e\n",
                                         "Plummer_Collision_D", Plummer_Collision_D );
   }

   for (int d=0; d<3; d++)
   if ( Plummer_Center[d] < 0.0 )
   {
      Plummer_Center[d] = 0.5*amr->BoxSize[d];

      if ( MPI_Rank == 0 )  Aux_Message( stdout, "NOTE : parameter \"%s[%d]\" is set to the default value = %13.7e\n",
                                         "Plummer_Center", d, Plummer_Center[d] );
   }

   if ( Plummer_RSeed < 0   )    Aux_Error( ERROR_INFO, "Plummer_RSeed (%d) < 0 !!\n",         Plummer_RSeed );
   if ( Plummer_MaxR  < 0.0 )    Aux_Error( ERROR_INFO, "Plummer_MaxR (%14.7e) < 0.0 !!\n",    Plummer_MaxR );
   if ( Plummer_Rho0  < 0.0 )    Aux_Error( ERROR_INFO, "Plummer_Rho0 (%14.7e) < 0.0 !!\n",    Plummer_Rho0 );
   if ( Plummer_R0    < 0.0 )    Aux_Error( ERROR_INFO, "Plummer_R0 (%14.7e) < 0.0 !!\n",      Plummer_R0   );
   if ( Plummer_NBinR <= 1  )    Aux_Error( ERROR_INFO, "Plummer_NBinR (%d) <= 1 !!\n",        Plummer_NBinR );
#  if ( MODEL == HYDRO )
   if ( Plummer_GasMFrac <= 0.0  ||  Plummer_GasMFrac > 1.0 )
      Aux_Error( ERROR_INFO, "Plummer_GasMFrac (%14.7e) is not within the correct range [0.0 < x <= 1.0] !!\n", Plummer_GasMFrac );
#  endif

} // FUNCTION : LoadTestProbParameter



#endif // #ifdef PARTICLE
