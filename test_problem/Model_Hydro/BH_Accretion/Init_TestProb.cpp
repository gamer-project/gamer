#include "GAMER.h"

#if ( MODEL == HYDRO )



extern void (*Init_Function_Ptr)( real fluid[], const double x, const double y, const double z, const double Time );
extern void (*Output_TestProbErr_Ptr)( const bool BaseOnly );

       void HYDRO_TestProbSol_BHAccretion( real fluid[], const double x, const double y, const double z, const double Time );
static void LoadTestProbParameter();


// global variables in the HYDRO BH accretion test
// =======================================================================================
double BH_MassBH;          // black hole mass (in internal units)
double BH_Rho0;            // background density (in internal units)
double BH_T0;              // background temperature (in internal units)
double BH_Mu;              // mean atomic weight
double BH_RefineRadius0;   // refinement radius at the base level (in internal units)
                           // NOTE: refinement radius at Lv is set to BH_RefineRadius0*2^(-Lv)
                           // --> all refinement shells have roughly the same number of cells at its level
                           //     (except Lv=MAX_LEVEL, which will have twice the number of cells along the radius
                           //      unless BH_HalfMaxLvRefR is on)
bool   BH_HalfMaxLvRefR;   // halve the refinement radius at the maximum level
double BH_InBC_Rho;        // density     inside the void region (in internal units)
double BH_InBC_T;          // temperature inside the void region (in internal units)
double BH_InBC_NCell;      // number of finest cells for the radius of the void region
double BH_Soften_NCell;    // number of finest cells for the soften length (<0.0 ==> disable)

double BH_InBC_R;          // radius of the void region (=BH_InBC_NCell*dh[MAX_LEVEL])
double BH_InBC_E;          // energy inside the void region (
double BH_Soften_R;        // soften length (=BH_Soften_NCell*dh[MAX_LEVEL])
double BH_Cs;              // background sound speed
double BH_RS;              // Schwarzschild radius (in kpc)
double BH_RB;              // Bondi radius (in kpc)
double BH_TimeB;           // Bondi time (in Myr)

double BH_SinkMass;        // total mass             in the void region removed in one global time-step
double BH_SinkMomX;        // total x-momentum       ... 
double BH_SinkMomY;        // total y-momentum       ... 
double BH_SinkMomZ;        // total z-momentum       ... 
double BH_SinkMomXAbs;     // total |x-momentum|     ... 
double BH_SinkMomYAbs;     // total |y-momentum|     ... 
double BH_SinkMomZAbs;     // total |z-momentum|     ... 
double BH_SinkEk;          // total kinematic energy ...
double BH_SinkEt;          // total thermal   energy ...
int    BH_SinkNCell;       // total number of finest cells within the void region
// =======================================================================================

// some constants and units (in CGS)
// =======================================================================================
// constants
double keV        = 1.60217646e-9;
double kpc        = 3.08568025e21;
double km         = 1.0e5;
double Myr        = 3.15576e13;
double Mole       = 6.02214078e23;
double Msun       = 1.9885e33;
double HydroMass  = 1.0;               // mass of neutral atomic hydrogen per mole
double SpeedC     = 2.99792458e10;     // speed of light
double NewtonG    = 6.67408e-8;        // gravitational constant

// internal units
double UnitI_Leng = kpc;
double UnitI_Dens = 5.0e-25;
double UnitI_Velo = 1.0e8;                                           // 10^3 km/s
double UnitI_Mass = UnitI_Dens*CUBE(UnitI_Leng);
double UnitI_Engy = UnitI_Mass*SQR(UnitI_Velo);
double UnitI_Time = UnitI_Leng/UnitI_Velo;                           // ~1.0 Myr

// external units
double UnitE_Leng = kpc;
double UnitE_Dens = 1.0;
double UnitE_Mass = Msun;
double UnitE_Engy = keV;
// =======================================================================================




//-------------------------------------------------------------------------------------------------------
// Function    :  Init_TestProb
// Description :  Initialize parameters for the HYDRO BH accretion test
//
// Note        :  1. Please copy this file to "GAMER/src/Init/Init_TestProb.cpp"
//                2. Global variables declared here will also be used in the function
//                   "HYDRO_TestProbSol_BHAccretion"
//
// Parameter   :  None 
//-------------------------------------------------------------------------------------------------------
void Init_TestProb()
{  

   const char *TestProb = "HYDRO BH accretion";

// check
#  if ( MODEL != HYDRO )
#  error : ERROR : "MODEL != HYDRO" in the HYDRO BH accretion test !!
#  endif

#  ifndef GRAVITY
#  error : ERROR : "GRAVITY must be ON" in the HYDRO BH accretion test !!
#  endif

#  ifdef COMOVING
#  error : ERROR : "COMOVING must be OFF" in the HYDRO BH accretion test !!
#  endif


// set the initialization and output functions
   Init_Function_Ptr      = HYDRO_TestProbSol_BHAccretion;
   Output_TestProbErr_Ptr = NULL;


// load the test problem parameters
   LoadTestProbParameter();


// set global variables
   if ( NEWTON_G <= 0.0 )
   NEWTON_G    = NewtonG/CUBE(UnitI_Leng)*UnitI_Mass*SQR(UnitI_Time);   // set G in internal units

   BH_InBC_R   = BH_InBC_NCell*amr->dh[MAX_LEVEL];
   BH_InBC_E   = BH_InBC_Rho*BH_InBC_T/(BH_Mu*HydroMass/Mole/UnitI_Mass)/(GAMMA-1.0);
   BH_Soften_R = BH_Soften_NCell*amr->dh[MAX_LEVEL];
   BH_Cs       = sqrt( GAMMA*BH_T0/(BH_Mu*HydroMass/Mole/UnitI_Mass) );
   BH_RS       = 2.0*NEWTON_G*BH_MassBH/SQR(SpeedC/UnitI_Velo);
   BH_RB       =     NEWTON_G*BH_MassBH/SQR(BH_Cs);
   BH_TimeB    = BH_RB/BH_Cs;


// record the test problem parameters
   if ( MPI_Rank == 0 )
   {
      Aux_Message( stdout, "\n" );
      Aux_Message( stdout, "%s test :\n", TestProb );
      Aux_Message( stdout, "=============================================================================\n" );
      Aux_Message( stdout, "Physical constants :\n" );
      Aux_Message( stdout, "  BH_MassBH              = %13.7e (%13.7e Msun)\n",   BH_MassBH, BH_MassBH*UnitI_Mass/Msun              );
      Aux_Message( stdout, "  BH_Rho0                = %13.7e (%13.7e g/cm^3)\n", BH_Rho0, BH_Rho0*UnitI_Dens                       );
      Aux_Message( stdout, "  BH_T0                  = %13.7e (%13.7e keV)\n",    BH_T0, BH_T0*UnitI_Engy/keV                       );
      Aux_Message( stdout, "  BH_Mu                  = %13.7e\n",                 BH_Mu                                             );
      Aux_Message( stdout, "  BH_RefineRadius0       = %13.7e (%13.7e kpc)\n",    BH_RefineRadius0, BH_RefineRadius0*UnitI_Leng/kpc );
      Aux_Message( stdout, "  BH_HalfMaxLvRefR       = %s\n",                     (BH_HalfMaxLvRefR)?"YES":"NO"                     );
      Aux_Message( stdout, "  BH_InBC_Rho            = %13.7e (%13.7e g/cm^3)\n", BH_InBC_Rho, BH_InBC_Rho*UnitI_Dens               );
      Aux_Message( stdout, "  BH_InBC_T              = %13.7e (%13.7e keV)\n",    BH_InBC_T, BH_InBC_T*UnitI_Engy/keV               );
      Aux_Message( stdout, "  BH_InBC_NCell          = %13.7e\n",                 BH_InBC_NCell                                     );
      Aux_Message( stdout, "  BH_InBC_R              = %13.7e (%13.7e kpc)\n",    BH_InBC_R, BH_InBC_R*UnitI_Leng/kpc               );
      Aux_Message( stdout, "  BH_InBC_E              = %13.7e\n",                 BH_InBC_E                                         );
      Aux_Message( stdout, "  BH_Soften_NCell        = %13.7e\n",                 BH_Soften_NCell                                   );
      Aux_Message( stdout, "  BH_Soften_R            = %13.7e (%13.7e kpc)\n",    BH_Soften_R, BH_Soften_R*UnitI_Leng/kpc           );
      Aux_Message( stdout, "  BH_Cs                  = %13.7e (%13.7e km/s)\n",   BH_Cs, BH_Cs*UnitI_Velo/km                        );
      Aux_Message( stdout, "  Schwarzschild radius   = %13.7e (%13.7e kpc)\n",    BH_RS, BH_RS*UnitI_Leng/kpc                       );
      Aux_Message( stdout, "  Bondi         radius   = %13.7e (%13.7e kpc)\n",    BH_RB, BH_RB*UnitI_Leng/kpc                       );
      Aux_Message( stdout, "  Bondi         time     = %13.7e (%13.7e Myr)\n",    BH_TimeB, BH_TimeB*UnitI_Time/Myr                 );

      Aux_Message( stdout, "\n" );
      Aux_Message( stdout, "Internal units in CGS :\n" );
      Aux_Message( stdout, "  Length                 = %13.7e\n", UnitI_Leng );
      Aux_Message( stdout, "  Density                = %13.7e\n", UnitI_Dens );
      Aux_Message( stdout, "  Velocity               = %13.7e\n", UnitI_Velo );
      Aux_Message( stdout, "  Mass                   = %13.7e\n", UnitI_Mass );
      Aux_Message( stdout, "  Energy                 = %13.7e\n", UnitI_Engy );
      Aux_Message( stdout, "  Time                   = %13.7e\n", UnitI_Time );
      Aux_Message( stdout, "=============================================================================\n" );
      Aux_Message( stdout, "\n" );
   } // if ( MPI_Rank == 0 )


// set some default parameters
// End_T : 2 Bondi time
   const double End_T_Default    = 2.0*BH_TimeB;
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

   if ( MPI_Rank == 0 )
   {
      if ( OPT__GRAVITY_TYPE != GRAVITY_EXTERNAL )
         Aux_Message( stderr, "WARNING : OPT__GRAVITY_TYPE != GRAVITY_EXTERNAL (are you sure?) !!\n" );
   }

} // FUNCTION : Init_TestProb



//-------------------------------------------------------------------------------------------------------
// Function    :  HYDRO_TestProbSol_BHAccretion
// Description :  Calculate the analytical solution in the HYDRO BH accretion test  
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
void HYDRO_TestProbSol_BHAccretion( real fluid[], const double x, const double y, const double z, const double Time )
{

   fluid[DENS] = BH_Rho0;
   fluid[MOMX] = 0.0;
   fluid[MOMY] = 0.0;
   fluid[MOMZ] = 0.0;
   fluid[ENGY] = SQR(BH_Cs)*BH_Rho0/( GAMMA*(GAMMA-1.0) );

} // FUNCTION : HYDRO_TestProbSol_BHAccretion



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

// all floating-point variables are declared as double
   getline( &input_line, &len, File );
   sscanf( input_line, "%lf%s",  &BH_MassBH,           string );

   getline( &input_line, &len, File );
   sscanf( input_line, "%lf%s",  &BH_Rho0,             string );

   getline( &input_line, &len, File );
   sscanf( input_line, "%lf%s",  &BH_T0,               string );

   getline( &input_line, &len, File );
   sscanf( input_line, "%lf%s",  &BH_Mu,               string );

   getline( &input_line, &len, File );
   sscanf( input_line, "%lf%s",  &BH_RefineRadius0,    string );

   getline( &input_line, &len, File );
   sscanf( input_line, "%d%s",   &temp_int,            string );
   BH_HalfMaxLvRefR = (bool)temp_int;

   getline( &input_line, &len, File );
   sscanf( input_line, "%lf%s",  &BH_InBC_Rho,         string );

   getline( &input_line, &len, File );
   sscanf( input_line, "%lf%s",  &BH_InBC_T,           string );

   getline( &input_line, &len, File );
   sscanf( input_line, "%lf%s",  &BH_InBC_NCell,       string );

   getline( &input_line, &len, File );
   sscanf( input_line, "%lf%s",  &BH_Soften_NCell,     string );

   fclose( File );
   if ( input_line != NULL )     free( input_line );


// convert from external to internal units
   BH_MassBH        *= UnitE_Mass/UnitI_Mass;
   BH_Rho0          *= UnitE_Dens/UnitI_Dens;
   BH_T0            *= UnitE_Engy/UnitI_Engy;
   BH_RefineRadius0 *= UnitE_Leng/UnitI_Leng;
   BH_InBC_Rho      *= UnitE_Dens/UnitI_Dens;
   BH_InBC_T        *= UnitE_Engy/UnitI_Engy;

} // FUNCTION : LoadTestProbParameter



#endif // #if ( MODEL == HYDRO )
