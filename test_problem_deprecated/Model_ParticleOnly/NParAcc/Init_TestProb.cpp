#include "GAMER.h"

#ifdef PARTICLE



extern void (*Init_Function_Ptr)( real fluid[], const double x, const double y, const double z, const double Time );
extern void (*Output_TestProbErr_Ptr)( const bool BaseOnly );

static void LoadTestProbParameter();
static void Par_TestProbSol_NParAcc( real fluid[], const double x, const double y, const double z, const double Time );
static void Par_OutputError_NParAcc( const bool Useless );

extern double MassProf_Plummer( const double r );
extern double MassProf_Burkert( const double r );
extern double MassProf_NFW    ( const double r );


// global variables in the N particles force test
// =======================================================================================
int  NParAcc_Mode;         // 1/2/3 : Plummer/Burkert/NFW
int  NParAcc_RSeed;        // random seed for setting particle position
real NParAcc_MaxR;         // maximum distance from the origin
real NParAcc_Rho0;         // density parameter in the assumed density profile
real NParAcc_R0;           // scale radius in the assumed density profile
int  NParAcc_NBinR;        // number of radial bins in the mass profile table
bool NParAcc_ComprDirN;    // calculate the direct N-body acceleration for comparison
int  NParAcc_AddSoft;      // 0/1/2: no soften / (e^2+r^2)^(3/2) / (e^5+r^5)^(5/3)
real NParAcc_Soft;         // soften length
// =======================================================================================




//-------------------------------------------------------------------------------------------------------
// Function    :  Init_TestProb
// Description :  Initialize parameters for the N particles force test
//
// Note        :  1. Please copy this file to "GAMER/src/Init/Init_TestProb.cpp"
//                2. Test problem parameters can be set in the input file "Input__TestProb"
//
// Parameter   :  None
//-------------------------------------------------------------------------------------------------------
void Init_TestProb()
{

   const char *TestProb = "N particles force";

// check
# if ( MODEL != ELBDM )
# error : ERROR : "MODEL != ELBDM" in the N particles force test !!
# endif

# ifndef PARTICLE
# error : ERROR : "PARTICLE is NOT defined" in the N particles force test !!
# endif

# ifndef STORE_PAR_ACC
# error : ERROR : "STORE_PAR_ACC must be ON" in the two particles force test !!
# endif

# ifndef GRAVITY
# error : ERROR : "GRAVITY must be ON" in the N particles force test !!
# endif

# ifdef COMOVING
# error : ERROR : "COMOVING must be OFF" in the N particles force test !!
# endif

# ifndef SERIAL
# error : ERROR : "Currently only support SERIAL" in the N particles force test !!
# endif

# ifdef TIMING
# error : ERROR : "TIMING must be OFF" in the N particles force test !!
# endif

   if ( OPT__BC_POT != BC_POT_ISOLATED )
      Aux_Error( ERROR_INFO, "pleaset set parameter %s to %d in the %s test !!\n",
                 "BC_POT_ISOLATED", BC_POT_ISOLATED, TestProb );

   if ( amr->Par->Init != PAR_INIT_BY_FUNCTION )
      Aux_Error( ERROR_INFO, "please set parameter %s to %d in the %s test !!\n",
                 "amr->Par->Init", amr->Par->Init, TestProb );

   if ( OPT__INIT != INIT_STARTOVER )
      Aux_Error( ERROR_INFO, "please set parameter %s to %d in the %s test !!\n",
                 "OPT__INIT", OPT__INIT, TestProb );

   if ( OPT__EXTERNAL_POT )
      Aux_Error( ERROR_INFO, "please disable the parameter %s in the %s test !!\n",
                 "OPT__EXTERNAL_POT", TestProb );

   if ( !OPT__OUTPUT_TEST_ERROR )
      Aux_Error( ERROR_INFO, "please turn on the option %s in the %s test !!\n",
                 "OPT__OUTPUT_TEST_ERROR", OPT__OUTPUT_TEST_ERROR, TestProb );


// set the initialization and output functions
   Init_Function_Ptr      = Par_TestProbSol_NParAcc;
   Output_TestProbErr_Ptr = Par_OutputError_NParAcc;


// load the test problem parameters
   LoadTestProbParameter();


// set the test problem parameters


// record the test problem parameters
   if ( MPI_Rank == 0 )
   {
      Aux_Message( stdout, "\n" );
      Aux_Message( stdout, "%s test :\n", TestProb );
      Aux_Message( stdout, "=============================================================================\n"   );
      Aux_Message( stdout, " Note: test mode                                    = %d\n",     NParAcc_Mode      );
      Aux_Message( stdout, "       random seed for setting particle position    = %d\n",     NParAcc_RSeed     );
      Aux_Message( stdout, "       maximum radius                               = %13.7e\n", NParAcc_MaxR      );
      Aux_Message( stdout, "       density parameter                            = %13.7e\n", NParAcc_Rho0      );
      Aux_Message( stdout, "       scale radius                                 = %13.7e\n", NParAcc_R0        );
      Aux_Message( stdout, "       number of radial bins in the mass profile    = %d\n",     NParAcc_NBinR     );
      Aux_Message( stdout, "       compute direct N-body force for comparison   = %d\n",     NParAcc_ComprDirN );
      if ( NParAcc_ComprDirN ) {
      Aux_Message( stdout, "       soften length mode                           = %d\n",     NParAcc_AddSoft   );
      if ( NParAcc_AddSoft )
      Aux_Message( stdout, "       soften length                                = %13.7e\n", NParAcc_Soft      ); }
      Aux_Message( stdout, "=============================================================================\n"   );
      Aux_Message( stdout, "\n" );
   }


// set some default parameters
   const double End_T_Default    = 0.0;
   const long   End_Step_Default = 0;

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
// Function    :  Par_TestProbSol_NParAcc
// Description :  Initialize the background density field as zero for the N particles force test
//
// Note        :  1. Currently particle test must work with the ELBDM model
//                2. Invoked by "ELBDM_Init_StartOver_AssignData"
//
// Parameter   :  fluid : Fluid field to be initialized
//                x/y/z : Target physical coordinates
//                Time  : Target physical time
//
// Return      :  fluid
//-------------------------------------------------------------------------------------------------------
void Par_TestProbSol_NParAcc( real *fluid, const double x, const double y, const double z, const double Time )
{

// set wave function as zero everywhere
   fluid[REAL] = 0.0;
   fluid[IMAG] = 0.0;
   fluid[DENS] = 0.0;

} // FUNCTION : Par_TestProbSol_NParAcc



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
   sscanf( input_line, "%d%s",   &NParAcc_Mode,           string );

   getline( &input_line, &len, File );
   sscanf( input_line, "%d%s",   &NParAcc_RSeed,          string );

#  ifdef FLOAT8
   getline( &input_line, &len, File );
   sscanf( input_line, "%lf%s",  &NParAcc_MaxR,           string );

   getline( &input_line, &len, File );
   sscanf( input_line, "%lf%s",  &NParAcc_Rho0,           string );

   getline( &input_line, &len, File );
   sscanf( input_line, "%lf%s",  &NParAcc_R0,             string );
#  else
   getline( &input_line, &len, File );
   sscanf( input_line, "%f%s",   &NParAcc_MaxR,           string );

   getline( &input_line, &len, File );
   sscanf( input_line, "%f%s",  &NParAcc_Rho0,            string );

   getline( &input_line, &len, File );
   sscanf( input_line, "%f%s",  &NParAcc_R0,              string );
#  endif // #ifdef FLOAT8 ... else ...

   getline( &input_line, &len, File );
   sscanf( input_line, "%d%s",   &NParAcc_NBinR,          string );

   getline( &input_line, &len, File );
   sscanf( input_line, "%d%s",   &temp_int,                 string );
   NParAcc_ComprDirN = (bool)temp_int;

   getline( &input_line, &len, File );
   sscanf( input_line, "%d%s",   &NParAcc_AddSoft,        string );

#  ifdef FLOAT8
   getline( &input_line, &len, File );
   sscanf( input_line, "%lf%s",  &NParAcc_Soft,           string );
#  else
   getline( &input_line, &len, File );
   sscanf( input_line, "%f%s",   &NParAcc_Soft,           string );
#  endif

   fclose( File );
   if ( input_line != NULL )     free( input_line );


// set the default test problem parameters
   if ( NParAcc_MaxR <= 0.0 )
   {
//    since the outmost base patches are not allowed to be refined for the isolated gravity solver,
//    we set MaxR = 0.5*L-2*PatchSize to ensure that all particles can be refined
      NParAcc_MaxR = 0.5*amr->BoxSize[0] - 2.0*PS1*amr->dh[0];

      if ( MPI_Rank == 0 )  Aux_Message( stdout, "NOTE : parameter \"%s\" is set to the default value = %13.7e\n",
                                         "NParAcc_MaxR", NParAcc_MaxR );
   }

   if ( NParAcc_Rho0 <= 0.0 )
   {
      NParAcc_Rho0 = 1.0;

      if ( MPI_Rank == 0 )  Aux_Message( stdout, "NOTE : parameter \"%s\" is set to the default value = %13.7e\n",
                                         "NParAcc_Rho0", NParAcc_Rho0 );
   }

   if ( NParAcc_R0 <= 0.0 )
   {
      NParAcc_R0 = PS1*amr->dh[0];

      if ( MPI_Rank == 0 )  Aux_Message( stdout, "NOTE : parameter \"%s\" is set to the default value = %13.7e\n",
                                         "NParAcc_R0", NParAcc_R0 );
   }

   if ( NParAcc_NBinR <= 1 )
   {
      NParAcc_NBinR = 1000;

      if ( MPI_Rank == 0 )  Aux_Message( stdout, "NOTE : parameter \"%s\" is set to the default value = %d\n",
                                         "NParAcc_NBinR", NParAcc_NBinR );
   }

   if ( NParAcc_ComprDirN  &&  NParAcc_AddSoft != 0  &&  NParAcc_Soft < 0.0 )
   {
      if ( MAX_LEVEL == 0 )   NParAcc_Soft = 1.20*amr->dh[MAX_LEVEL];
      else                    NParAcc_Soft = 1.40*amr->dh[MAX_LEVEL];

      if ( MPI_Rank == 0 )  Aux_Message( stdout, "NOTE : parameter \"%s\" is set to the default value = %13.7e\n",
                                         "NParAcc_Soft", NParAcc_Soft );
   }

// check
   if ( NParAcc_Mode < 1  ||  NParAcc_Mode > 3 )
      Aux_Error( ERROR_INFO, "NParAcc_Mode (%d) != 1/2/3 !!\n", NParAcc_Mode );

   if ( NParAcc_RSeed < 0   )    Aux_Error( ERROR_INFO, "NParAcc_RSeed (%d) < 0 !!\n",         NParAcc_RSeed );
   if ( NParAcc_MaxR  < 0.0 )    Aux_Error( ERROR_INFO, "NParAcc_MaxR (%14.7e) < 0.0 !!\n",    NParAcc_MaxR );
   if ( NParAcc_Rho0  < 0.0 )    Aux_Error( ERROR_INFO, "NParAcc_Rho0 (%14.7e) < 0.0 !!\n",    NParAcc_Rho0 );
   if ( NParAcc_R0    < 0.0 )    Aux_Error( ERROR_INFO, "NParAcc_R0 (%14.7e) < 0.0 !!\n",      NParAcc_R0   );
   if ( NParAcc_NBinR <= 1  )    Aux_Error( ERROR_INFO, "NParAcc_NBinR (%d) <= 1 !!\n",        NParAcc_NBinR );
   if (  NParAcc_ComprDirN  &&  (NParAcc_AddSoft < 0 || NParAcc_AddSoft > 2 )  )
                                 Aux_Error( ERROR_INFO, "NParAcc_AddSoft (%d) != 0/1/2  !!\n", NParAcc_AddSoft );
   if (  NParAcc_ComprDirN  &&  NParAcc_AddSoft != 0  &&  NParAcc_Soft < 0.0  )
                                 Aux_Error( ERROR_INFO, "NParAcc_Soft (%14.7e) < 0.0 !!\n",    NParAcc_Soft);

} // FUNCTION : LoadTestProbParameter



//------------------------------------------------------------------------------------------------------
// Function    :  Par_OutputError_NParAcc
// Description :  Calculate and output the particle acceleration
//
// Parameter   :  Useless  : BaseOnly parameter in "Output_TestProbErr_Ptr" --> useless here
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Par_OutputError_NParAcc( const bool Useless )
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Calculating particle acceleration with AMR ...\n" );


   const bool StoreAcc_Yes    = true;
   const bool UseStoredAcc_No = false;

   for (int lv=0; lv<=MAX_LEVEL; lv++)
      Par_UpdateParticle( lv, amr->PotSgTime[lv][ amr->PotSg[lv] ], NULL_REAL, PAR_UPSTEP_ACC_ONLY,
                          StoreAcc_Yes, UseStoredAcc_No );


// calculate the analytical (direct N-body) force
   const double Cen[3] = { 0.5*amr->BoxSize[0],
                           0.5*amr->BoxSize[1],
                           0.5*amr->BoxSize[2] };
   const double Eps2   = SQR( NParAcc_Soft );
   const double Eps5   = Eps2*CUBE( NParAcc_Soft );

   real *Mass           =   amr->Par->Mass;
   real *Pos[3]         = { amr->Par->PosX, amr->Par->PosY, amr->Par->PosZ };
   real *Acc[3]         = { amr->Par->AccX, amr->Par->AccY, amr->Par->AccZ };
   double *AccAna[3]    = { NULL, NULL, NULL };
   double dr[3], r, _r, r2, r3, Fac, AccNumR, AccNumT, AccNumA, AccAnaR, AccAnaT, AccAnaA, AccProR;

   for (int d=0; d<3; d++)    AccAna[d] = new double [amr->Par->NPar_Active_AllRank];

   for (int d=0; d<3; d++)
   for (int p=0; p<amr->Par->NPar_Active_AllRank; p++)   AccAna[d][p] = 0.0;

   if ( NParAcc_ComprDirN )
   {
      if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Calculating direct N-body acceleration...\n" );

      for (int p1=0;    p1<amr->Par->NPar_Active_AllRank; p1++)
      for (int p2=p1+1; p2<amr->Par->NPar_Active_AllRank; p2++)
      {
         for (int d=0; d<3; d++)    dr[d] = Pos[d][p2] - Pos[d][p1];

         r   = sqrt( SQR(dr[0]) + SQR(dr[1]) + SQR(dr[2]) );
         r2  = r*r;
         r3  = r*r2;

         if      ( NParAcc_AddSoft == 0 )    Fac = NEWTON_G/r3;
         else if ( NParAcc_AddSoft == 1 )    Fac = NEWTON_G*pow( r2   +Eps2, -3.0/2.0 );
         else if ( NParAcc_AddSoft == 2 )    Fac = NEWTON_G*pow( r2*r3+Eps5, -3.0/5.0 );

         for (int d=0; d<3; d++)
         {
            AccAna[d][p1] += Fac*dr[d]*Mass[p2];
            AccAna[d][p2] -= Fac*dr[d]*Mass[p1];
         }
      }
   } // if ( NParAcc_ComprDirN )


// set the mass profile of the target model
   double (*MassProf)( const double r ) = NULL;

   switch ( NParAcc_Mode )
   {
      case 1:  MassProf = MassProf_Plummer;   break;
      case 2:  MassProf = MassProf_Burkert;   break;
      case 3:  MassProf = MassProf_NFW;       break;

      default: Aux_Error( ERROR_INFO, "unsupported mode (%d) !!\n", NParAcc_Mode );
   }


// output results
   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Outputting results ...\n" );

   const char FileName[] = "Record__NParAcc";

   if ( Aux_CheckFileExist(FileName) )
      Aux_Message( stderr, "WARNING : file \"%s\" already exists and will be overwritten !!\n", FileName );

   FILE *File = fopen( FileName, "w" );
   fprintf( File, "#ID    : Particle ID\n" );
   fprintf( File, "#Mass  : Particle mass\n" );
   fprintf( File, "#X/Y/Z : Particle positions relative to box center or acceleration directions\n" );
   fprintf( File, "#R/T/A : Radial, tangential, absolute components\n" );
   fprintf( File, "#r     : Distance from the origin\n" );
   fprintf( File, "#Acc   : Gravitational acceleration on particles\n" );
   fprintf( File, "#Num   : Numerical result\n" );
   fprintf( File, "#Ana   : Analytical result based on direct N-body calculation\n" );
   fprintf( File, "#Pro   : Analytical result based on spherically symmetric density profile\n" );
   fprintf( File, "#------------------------------------------------------------------------------\n" );
   fprintf( File, "#%9s  %13s  %13s  %13s  %13s  %13s  %13s  %13s  %13s",
            "[1]ID", "[2]Mass", "[3]X", "[4]Y", "[5]Z", "[6]r", "[7]AccNumX", "[8]AccNumY", "[9]AccNumZ" );
   fprintf( File, "  %13s  %13s  %13s  %13s  %13s  %13s  %13s",
            "[10]AccNumR", "[11]AccAnaR", "[12]AccProR", "[13]AccNumT", "[14]AccAnaT", "[15]AccNumA", "[16]AccAnaA" );
   fprintf( File, "\n" );

   for (long p=0; p<amr->Par->NPar_Active_AllRank; p++)
   {
      for (int d=0; d<3; d++)    dr[d] = Pos[d][p] - Cen[d];

      r  = sqrt( SQR(dr[0]) + SQR(dr[1]) + SQR(dr[2]) );
      _r = 1.0/r;

      AccNumA = sqrt( SQR(Acc   [0][p]) + SQR(Acc   [1][p]) + SQR(Acc   [2][p]) );
      AccAnaA = sqrt( SQR(AccAna[0][p]) + SQR(AccAna[1][p]) + SQR(AccAna[2][p]) );

      AccNumR = ( Acc   [0][p]*dr[0] + Acc   [1][p]*dr[1] + Acc   [2][p]*dr[2] )*_r;
      AccAnaR = ( AccAna[0][p]*dr[0] + AccAna[1][p]*dr[1] + AccAna[2][p]*dr[2] )*_r;
      AccProR = -NEWTON_G*MassProf(r)*SQR(_r);

      AccNumT = sqrt( SQR(AccNumA) - SQR(AccNumR) );
      AccAnaT = sqrt( SQR(AccAnaA) - SQR(AccAnaR) );

      fprintf( File, "%10ld", p );
      fprintf( File, "  %13.6e  %13.6e  %13.6e  %13.6e  %13.6e  %13.6e  %13.6e  %13.6e",
               Mass[p], dr[0], dr[1], dr[2], r, Acc[0][p], Acc[1][p], Acc[2][p] );
      fprintf( File, "  %13.6e  %13.6e  %13.6e  %13.6e  %13.6e  %13.6e  %13.6e",
               AccNumR, AccAnaR, AccProR, AccNumT, AccAnaT, AccNumA, AccAnaA );

      fprintf( File, "\n" );
   } // for (long p=0; p<amr->Par->NPar_Active_AllRank; p++)


   fclose( File );

   for (int d=0; d<3; d++)    delete [] AccAna[d];


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : Par_OutputError_NParAcc



#endif // #ifdef PARTICLE
