#include <random>
#include <limits.h>
#include <math.h>
#include "GAMER.h"
#include "TestProb.h"



void ***calloc_3d_array( size_t nt, size_t nr, size_t nc, size_t size );
void free_3d_array( void ***array );
void Mis_Cartesian2Spherical( const double Cartesian[], double Spherical[] );
void CartesianRotate( double x[], double theta, double phi, bool inverse );
void Interpolation_UM_IC( real x, real y, real z, real ****Pri_input, real **XYZ, real *Pri_output, bool disk );
real TrilinearInterpolation( real *FieldAtVertices, real *xyz000, real *dxyz, real *xyz );
static void SetArrayDisk();
static void SetArrayHVC();
static real *randXYZ;
static int  numClouds   = 0;
static real radiusCloud = 1.0;
void randCloud( real **randXYZ, int numClouds );


// problem-specific global variables
// =======================================================================================

// options
       int      Jet_Ambient;             // [0/1/9]: uniform/Milky-Way/load-from-file
static int      Jet_Fire;                // [0/1/2/3]: no jet/upper jet/lower jet/bipolar jet
static double   Jet_Duration;            // a duration of jet injection from the start of simulation

// general parameters
static double   ParticleMass;            // atomic mass unit in jet source

// uniform background parameters
static double   Amb_UniformDens;         // uniform ambient density
static double   Amb_UniformVel[3];       // uniform ambient 4-velocity
static double   Amb_UniformTemp;         // uniform ambient temperature

// Milky Way parameters
       double   IsothermalSlab_Center[3];

// jet fluid parameters
static double   Jet_SrcVel;              // jet 4-velocity
static double   Jet_SrcGamma;            // jet lorentz factor
static double   Jet_SrcDens;             // jet density
static double   Jet_SrcTemp;             // jet temperature
static bool     Jet_SmoothVel;           // smooth radial component of 4-velocity on cross section

static bool     Jet_SphericalSrc;

#ifdef COSMIC_RAY
static double   Jet_Src_CR_Engy;
static double   Amb_CR_Engy;
#endif

// sound speed
static double   CharacteristicSpeed;     // the characteristic speed of the simulation problem
                                         // the default end-time (END_T) will be estimated from
                                         // `CharacteristicSpeed` and `BOX_SIZE`
static real *buffer;
static real *Header_disk;
static real ***Rhoo_disk;
static real ***VelX_disk;
static real ***VelY_disk;
static real ***VelZ_disk;
static real ***Pres_disk;
static real *Header_hvc;
static real ***Rhoo_hvc;
static real ***VelX_hvc;
static real ***VelY_hvc;
static real ***VelZ_hvc;
static real ***Pres_hvc;
static real *X_disk;
static real *Y_disk;
static real *Z_disk;
static real *X_hvc;
static real *Y_hvc;
static real *Z_hvc;

// fermi bubbles
       real   IsothermalSlab_VelocityDispersion;
       real   IsothermalSlab_PeakDens;
       real   IsothermalSlab_Truncation;
static real   ambientTemperature;
static real   gasDiskTemperature;
static real   gasDiskPeakDens;
       real   interfaceHeight;
static double criticalTemp;
static double gasDisk_highResRadius;
static int    gasDisk_lowRes_LEVEL;
static int    jetSrc_lowRes_LEVEL;

// Dark logarithmic halo potential
       real  v_halo;
       real  distance_h;

void Init_ExtPot_IsothermalSlab();

// =======================================================================================
/*        G       A       C              */
/*          ____________                 */
/*          \     |     /                */
/*           \   E|    /        z        */
/*            \  /|\  /         ^        */
/*             \/_|_\/          |        */
/*             /\O| /\B         |        */
/*            /  \|/  \                  */
/*           /   D|    \                 */
/*          /_____|_____\                */
/*                F                      */
// =======================================================================================
// jet geometry parameters
static double   Jet_Radius;              // length of OB
static double   Jet_HalfHeight;          // length of OA
static double   Jet_HalfOpeningAngle;    // half-opening angle (i.e. ∠ADC)
static double   Jet_CenOffset[3];        // jet central coordinates offset
static double   Jet_Center[3];           // jet central coordinates
static double   Jet_MaxDis;              // maximum distance between the jet source and the jet center


// precession parameters
static double   Jet_AngularVelocity;     // precession angular velocity (degree per code_time)
static double   Jet_PrecessionAngle;     // precession angle in degree
static double   Jet_PrecessionAxis[3];   // cone orientation vector (x,y,z). i.e. vector OA


#if ( NCOMP_PASSIVE_USER > 0 )
static FieldIdx_t Passive_0000 = 5;  // disk
static FieldIdx_t Passive_0001 = 6;  // src, the ejected material is Jet_Dens
static FieldIdx_t Passive_0002 = 7;  // src, the ejected material is CRay
#endif
// =======================================================================================



//-------------------------------------------------------------------------------------------------------
// Function    :  Validate
// Description :  Validate the compilation flags and runtime parameters for this test problem
//
// Note        :  None
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Validate()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Validating test problem %d ...\n", TESTPROB_ID );


// errors
#  if ( MODEL != HYDRO )
   Aux_Error( ERROR_INFO, "MODEL != HYDRO !!\n" );
#  endif

#  ifndef SRHD
   Aux_Error( ERROR_INFO, "SRHD must be enabled !!\n" );
#  endif

#  ifdef PARTICLE
   Aux_Error( ERROR_INFO, "PARTICLE must be disabled !!\n" );
#  endif

#  ifndef GRAVITY
   if ( Jet_Ambient == 1  ||  Jet_Ambient == 2 )
       Aux_Error( ERROR_INFO, "GRAVITY must be enabled !!\n" );
#  endif


// warnings
   if ( MPI_Rank == 0 )
   {
      for (int s=0; s<6; s++)
         if ( OPT__BC_FLU[s] != BC_FLU_OUTFLOW )
            Aux_Message( stderr, "WARNING : it's recommended to use the outflow BC (currently OPT__BC_FLU[%d] = %d != 2)\n",
                         s, OPT__BC_FLU[s] );
   }


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Validating test problem %d ... done\n", TESTPROB_ID );

} // FUNCTION : Validate



#if ( MODEL == HYDRO )
//-------------------------------------------------------------------------------------------------------
// Function    :  LoadInputTestProb
// Description :  Read problem-specific runtime parameters from Input__TestProb and store them in HDF5 snapshots (Data_*)
//
// Note        :  1. Invoked by SetParameter() to read parameters
//                2. Invoked by Output_DumpData_Total_HDF5() using the function pointer Output_HDF5_InputTest_Ptr to store parameters
//                3. If there is no problem-specific runtime parameter to load, add at least one parameter
//                   to prevent an empty structure in HDF5_Output_t
//                   --> Example:
//                       LOAD_PARA( load_mode, "TestProb_ID", &TESTPROB_ID, TESTPROB_ID, TESTPROB_ID, TESTPROB_ID );
//
// Parameter   :  load_mode      : Mode for loading parameters
//                                 --> LOAD_READPARA    : Read parameters from Input__TestProb
//                                     LOAD_HDF5_OUTPUT : Store parameters in HDF5 snapshots
//                ReadPara       : Data structure for reading parameters (used with LOAD_READPARA)
//                HDF5_InputTest : Data structure for storing parameters in HDF5 snapshots (used with LOAD_HDF5_OUTPUT)
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void LoadInputTestProb( const LoadParaMode_t load_mode, ReadPara_t *ReadPara, HDF5_Output_t *HDF5_InputTest )
{

#  ifndef SUPPORT_HDF5
   if ( load_mode == LOAD_HDF5_OUTPUT )   Aux_Error( ERROR_INFO, "please turn on SUPPORT_HDF5 in the Makefile for load_mode == LOAD_HDF5_OUTPUT !!\n" );
#  endif

   if ( load_mode == LOAD_READPARA     &&  ReadPara       == NULL )   Aux_Error( ERROR_INFO, "load_mode == LOAD_READPARA and ReadPara == NULL !!\n" );
   if ( load_mode == LOAD_HDF5_OUTPUT  &&  HDF5_InputTest == NULL )   Aux_Error( ERROR_INFO, "load_mode == LOAD_HDF5_OUTPUT and HDF5_InputTest == NULL !!\n" );

// add parameters in the following format:
// --> note that VARIABLE, DEFAULT, MIN, and MAX must have the same data type
// --> some handy constants (e.g., NoMin_int, Eps_float, ...) are defined in "include/ReadPara.h"
// --> LOAD_PARA() is defined in "include/TestProb.h"
// ****************************************************************************************************************************
// LOAD_PARA( load_mode, "KEY_IN_THE_FILE",         &VARIABLE,                  DEFAULT,       MIN,           MAX            );
// ****************************************************************************************************************************
   LOAD_PARA( load_mode, "Jet_Ambient",             &Jet_Ambient,               1,             0,             9              );
   LOAD_PARA( load_mode, "Jet_Fire",                &Jet_Fire,                  3,             0,             3              );
   LOAD_PARA( load_mode, "Jet_SphericalSrc",        &Jet_SphericalSrc,          false,         Useless_bool,  Useless_bool   );
   LOAD_PARA( load_mode, "Jet_Duration",            &Jet_Duration,              NoMax_double,  0.0,           NoMax_double   );

// load jet fluid parameters
   LOAD_PARA( load_mode, "Jet_SrcVel",              &Jet_SrcVel,               -1.0,           NoMin_double,  NoMax_double   );
   LOAD_PARA( load_mode, "Jet_SmoothVel",           &Jet_SmoothVel,             false,         Useless_bool,  Useless_bool   );
   LOAD_PARA( load_mode, "Jet_SrcDens",             &Jet_SrcDens,              -1.0,           Eps_double,    NoMax_double   );
   LOAD_PARA( load_mode, "Jet_SrcTemp",             &Jet_SrcTemp,              -1.0,           Eps_double,    NoMax_double   );
   LOAD_PARA( load_mode, "gasDisk_highResRadius",   &gasDisk_highResRadius,    -1.0,           NoMin_double,  NoMax_double   );
   LOAD_PARA( load_mode, "gasDisk_lowRes_LEVEL",    &gasDisk_lowRes_LEVEL,     -1,             0,             NoMax_int      );
   LOAD_PARA( load_mode, "jetSrc_lowRes_LEVEL",     &jetSrc_lowRes_LEVEL,      -1,             0,             NoMax_int      );
#  ifdef COSMIC_RAY
   LOAD_PARA( load_mode, "Jet_Src_CR_Engy",         &Jet_Src_CR_Engy,          -1.0,           0.0,           NoMax_double   );
   LOAD_PARA( load_mode, "Amb_CR_Engy",             &Amb_CR_Engy,              -1.0,           0.0,           NoMax_double   );
#  endif

// load source geometry parameters
   LOAD_PARA( load_mode, "Jet_Radius",              &Jet_Radius,               -1.0,           Eps_double,    NoMax_double   );
   LOAD_PARA( load_mode, "Jet_HalfHeight",          &Jet_HalfHeight,           -1.0,           Eps_double,    NoMax_double   );
   LOAD_PARA( load_mode, "Jet_HalfOpeningAngle",    &Jet_HalfOpeningAngle,     -1.0,           0.0,           90.0           );
   LOAD_PARA( load_mode, "Jet_CenOffset_x",         &Jet_CenOffset[0],          NoDef_double,  NoMin_double,  NoMax_double   );
   LOAD_PARA( load_mode, "Jet_CenOffset_y",         &Jet_CenOffset[1],          NoDef_double,  NoMin_double,  NoMax_double   );
   LOAD_PARA( load_mode, "Jet_CenOffset_z",         &Jet_CenOffset[2],          NoDef_double,  NoMin_double,  NoMax_double   );

// load precession parameters
   LOAD_PARA( load_mode, "Jet_AngularVelocity",     &Jet_AngularVelocity,       NoDef_double,  0.0,           NoMax_double   );
   LOAD_PARA( load_mode, "Jet_PrecessionAngle",     &Jet_PrecessionAngle,       NoDef_double,  NoMin_double,  90.0           );
   LOAD_PARA( load_mode, "Jet_PrecessionAxis_x",    &Jet_PrecessionAxis[0],     NoDef_double,  NoMin_double,  NoMax_double   );
   LOAD_PARA( load_mode, "Jet_PrecessionAxis_y",    &Jet_PrecessionAxis[1],     NoDef_double,  NoMin_double,  NoMax_double   );
   LOAD_PARA( load_mode, "Jet_PrecessionAxis_z",    &Jet_PrecessionAxis[2],     NoDef_double,  NoMin_double,  NoMax_double   );

// load uniform background parameters
   LOAD_PARA( load_mode, "Amb_UniformDens",         &Amb_UniformDens,          -1.0,           Eps_double,    NoMax_double   );
   LOAD_PARA( load_mode, "Amb_UniformVel_x",        &Amb_UniformVel[0],         0.0,           NoMin_double,  NoMax_double   );
   LOAD_PARA( load_mode, "Amb_UniformVel_y",        &Amb_UniformVel[1],         0.0,           NoMin_double,  NoMax_double   );
   LOAD_PARA( load_mode, "Amb_UniformVel_z",        &Amb_UniformVel[2],         0.0,           NoMin_double,  NoMax_double   );
   LOAD_PARA( load_mode, "Amb_UniformTemp",         &Amb_UniformTemp,          -1.0,           Eps_double,    NoMax_double   );


   LOAD_PARA( load_mode, "CharacteristicSpeed",     &CharacteristicSpeed,      -1.0,           NoMin_double,  NoMax_double   );
   LOAD_PARA( load_mode, "criticalTemp",            &criticalTemp,             -1.0,           NoMin_double,  NoMax_double   );

// load Milky Way parameters
   LOAD_PARA( load_mode, "IsothermalSlab_Center_x", &IsothermalSlab_Center[0], -1.0,           NoMin_double,  NoMax_double   );
   LOAD_PARA( load_mode, "IsothermalSlab_Center_y", &IsothermalSlab_Center[1], -1.0,           NoMin_double,  NoMax_double   );
   LOAD_PARA( load_mode, "IsothermalSlab_Center_z", &IsothermalSlab_Center[2], -1.0,           NoMin_double,  NoMax_double   );

} // FUNCITON : LoadInputTestProb



//-------------------------------------------------------------------------------------------------------
// Function    :  SetParameter
// Description :  Load and set the problem-specific runtime parameters
//
// Note        :  1. Filename is set to "Input__TestProb" by default
//                2. Major tasks in this function:
//                   (1) load the problem-specific runtime parameters
//                   (2) set the problem-specific derived parameters
//                   (3) reset other general-purpose parameters if necessary
//                   (4) make a note of the problem-specific parameters
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void SetParameter()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Setting runtime parameters ...\n" );


// (1) load the problem-specific runtime parameters
// (1-1) read parameters from Input__TestProb
   const char FileName[] = "Input__TestProb";
   ReadPara_t *ReadPara  = new ReadPara_t;

   LoadInputTestProb( LOAD_READPARA, ReadPara, NULL );

   ReadPara->Read( FileName );

   delete ReadPara;

// Read header for the Fermi bubbles
   SetArrayDisk();
// SetArrayHVC();
   randCloud( &randXYZ, numClouds );

// replace useless parameters with NaN
   if ( Jet_Ambient != 0 )
   {
      Amb_UniformDens   = NAN;
      Amb_UniformVel[0] = NAN;
      Amb_UniformVel[1] = NAN;
      Amb_UniformVel[2] = NAN;
   }

// (1-2) check runtime parameters

// check time-dependent source
   if ( IsothermalSlab_Center[0] == -1.0 )
      IsothermalSlab_Center[0] = 0.5*amr->BoxSize[0];

   if ( IsothermalSlab_Center[1] == -1.0 )
      IsothermalSlab_Center[1] = 0.5*amr->BoxSize[1];

   if ( IsothermalSlab_Center[2] == -1.0 )
      IsothermalSlab_Center[2] = 0.5*amr->BoxSize[2];

   if ( Jet_Ambient == 9  &&  OPT__INIT != 3 )
      Aux_Error( ERROR_INFO, "OPT__INIT must be 3 !!\n" );


// check UNIT_L is in reasonable range
   if ( ( UNIT_L <= 0.5*Const_kpc  ||  2.0*Const_kpc <= UNIT_L )  &&  OPT__UNIT )
      Aux_Error( ERROR_INFO, "UNIT_L=%e is far from %e !!\n", UNIT_L, Const_kpc );

   const double Const_Erg2eV = 6.2415e11;

// (1-2) convert to code unit
   Jet_SrcVel  *= Const_c  / UNIT_V;
   Jet_SrcTemp *= Const_kB / (ParticleMass*Const_c*Const_c);
   Jet_SrcDens *= 1.0      / UNIT_D;

#  ifdef COSMIC_RAY
   Jet_Src_CR_Engy *= 1.0 / UNIT_P;
   Amb_CR_Engy     *= 1.0 / UNIT_P;
#  endif

   Jet_Radius            *= Const_kpc / UNIT_L;
   Jet_HalfHeight        *= Const_kpc / UNIT_L;
   Jet_HalfOpeningAngle  *= M_PI      / 180.0;
   Jet_PrecessionAngle   *= M_PI      / 180.0;

   Jet_CenOffset[0]      *= Const_kpc / UNIT_L;
   Jet_CenOffset[1]      *= Const_kpc / UNIT_L;
   Jet_CenOffset[2]      *= Const_kpc / UNIT_L;

   Jet_Duration          *= Const_Myr / UNIT_T;

   gasDisk_highResRadius *= Const_kpc / UNIT_L;

   if ( Jet_Ambient == 0 )
   {
      Amb_UniformDens   *= 1.0     / UNIT_D;
      Amb_UniformVel[0] *= Const_c / UNIT_V;
      Amb_UniformVel[1] *= Const_c / UNIT_V;
      Amb_UniformVel[2] *= Const_c / UNIT_V;
   }
   else if ( Jet_Ambient == 2  ||  Jet_Ambient == 3  ||  Jet_Ambient == 4 )
   {
     IsothermalSlab_VelocityDispersion  = Header_disk[15];
     IsothermalSlab_PeakDens            = Header_disk[16];
     //IsothermalSlab_Truncation          = Header_disk[21];
     IsothermalSlab_Truncation          = 0.95*0.5*amr->BoxSize[2];
     gasDiskPeakDens                    = Header_disk[19];

     IsothermalSlab_VelocityDispersion *= 1e5       / UNIT_V; // km/s --> 1/c
     IsothermalSlab_PeakDens           *= 1.0       / UNIT_D;
     IsothermalSlab_Truncation         *= Const_kpc / UNIT_L;
     IsothermalSlab_Center[0]          *= Const_kpc / UNIT_L;
     IsothermalSlab_Center[1]          *= Const_kpc / UNIT_L;
     IsothermalSlab_Center[2]          *= Const_kpc / UNIT_L;

     distance_h                         = Header_disk[29];
     distance_h                        *= Const_kpc / UNIT_L;

     v_halo                             = Header_disk[30];
     v_halo                            *= 1.0       / UNIT_V;

     gasDiskPeakDens                   /= UNIT_D;

     ambientTemperature                 = Header_disk[20];
     ambientTemperature                *= Const_kB / (ParticleMass*UNIT_V*UNIT_V);

     gasDiskTemperature                 = Header_disk[18];
     gasDiskTemperature                *= Const_kB / (ParticleMass*UNIT_V*UNIT_V);

     criticalTemp                      *= Const_kB / (ParticleMass*Const_c*Const_c);
   } // if ( Jet_Ambient == 0 ) ... else if ...

   Amb_UniformTemp     *= Const_kB / (ParticleMass*Const_c*Const_c);
   Jet_AngularVelocity *= 1.0;    // the unit of Jet_AngularVelocity is UNIT_T


// (2) set the problem-specific derived parameters
   const double SecAngle = 1.0 / cos( 0.5*Jet_HalfOpeningAngle );
   const double TanAngle = sin( 0.5*Jet_HalfOpeningAngle ) * SecAngle;

   Jet_MaxDis  = sqrt( SQR( Jet_Radius ) + SQR( Jet_HalfHeight * SecAngle ) + 2.0 * Jet_Radius * Jet_HalfHeight * TanAngle );
   Jet_SrcGamma = sqrt(1.0 + SQR(Jet_SrcVel));

   for (int d=0; d<3; d++)    Jet_Center[d] = 0.5*amr->BoxSize[d] + Jet_CenOffset[d];


// (4) reset other general-purpose parameters
//     --> a helper macro PRINT_RESET_PARA is defined in Macro.h
   const long   End_Step_Default = __INT_MAX__;
   const double End_T_Default    = 0.5 * BOX_SIZE * UNIT_L / CharacteristicSpeed / UNIT_V / UNIT_T;

   if ( END_STEP < 0 )
   {
      END_STEP = End_Step_Default;
      PRINT_RESET_PARA( END_STEP, FORMAT_LONG, "" );
   }

   if ( END_T < 0.0 )
   {
      if ( CharacteristicSpeed == -1.0 )   Aux_Error( ERROR_INFO, "CharacteristicSpeed must be provided !!\n" );
      else                                 END_T = End_T_Default;

      PRINT_RESET_PARA( END_T, FORMAT_REAL, "" );
   }

   if ( OUTPUT_DT < 0.0 )
   {
      OUTPUT_DT = END_T / 30.0;
      PRINT_RESET_PARA( OUTPUT_DT, FORMAT_REAL, "" );
   }


// (4) make a note
   if ( MPI_Rank == 0 )
   {
      Aux_Message( stdout, "=============================================================================\n"                   );
      Aux_Message( stdout, "  test problem ID          = %d\n",                TESTPROB_ID                                     );
      Aux_Message( stdout, "  Jet_Ambient              = %d\n",                Jet_Ambient                                     );
      Aux_Message( stdout, "  Jet_Fire                 = %d\n",                Jet_Fire                                        );
      Aux_Message( stdout, "  Jet_SmoothVel            = %d\n",                Jet_SmoothVel                                   );
      Aux_Message( stdout, "  Jet_SphericalSrc         = %d\n",                Jet_SphericalSrc                                );
      Aux_Message( stdout, "  Jet_Duration             = %14.7e Myr \n",       Jet_Duration*UNIT_T/Const_Myr                   );
      Aux_Message( stdout, "  ParticleMass             = %14.7e g\n",          ParticleMass                                    );
      Aux_Message( stdout, "  Jet_SrcVel               = %14.7e c\n",          Jet_SrcVel                                      );
      Aux_Message( stdout, "  Jet_SrcDens              = %14.7e g/cm^3\n",     Jet_SrcDens*UNIT_D                              );
      Aux_Message( stdout, "  Jet_SrcTemp              = %14.7e kT/mc**2\n",   Jet_SrcTemp                                     );
#     ifdef COSMIC_RAY
      Aux_Message( stdout, "  Jet_Src_CR_Engy          = %14.7e \n",           Jet_Src_CR_Engy*UNIT_P                          );
      Aux_Message( stdout, "  Amb_CR_Engy              = %14.7e \n",           Amb_CR_Engy*UNIT_P                              );
#     endif
      Aux_Message( stdout, "  Jet_NumDensSrc           = %14.7e per cc\n",     Jet_SrcDens*UNIT_D/ParticleMass                 );
      Aux_Message( stdout, "  Jet_CenOffset[x]         = %14.7e kpc\n",        Jet_CenOffset[0]*UNIT_L/Const_kpc               );
      Aux_Message( stdout, "  Jet_CenOffset[y]         = %14.7e kpc\n",        Jet_CenOffset[1]*UNIT_L/Const_kpc               );
      Aux_Message( stdout, "  Jet_CenOffset[z]         = %14.7e kpc\n",        Jet_CenOffset[2]*UNIT_L/Const_kpc               );
      Aux_Message( stdout, "  Jet_AngularVelocity      = %14.7e degree/kyr\n", Jet_AngularVelocity                             );
      Aux_Message( stdout, "  Jet_PrecessionAngle      = %14.7e degree\n",     Jet_PrecessionAngle*180.0/M_PI                  );
      Aux_Message( stdout, "  Jet_HalfOpeningAngle     = %14.7e degree\n",     Jet_HalfOpeningAngle*180.0/M_PI                 );
      Aux_Message( stdout, "  Jet_Radius               = %14.7e kpc\n",        Jet_Radius*UNIT_L/Const_kpc                     );
      Aux_Message( stdout, "  Jet_HalfHeight           = %14.7e kpc\n",        Jet_HalfHeight*UNIT_L/Const_kpc                 );
      Aux_Message( stdout, "  Jet_MaxDis               = %14.7e kpc\n",        Jet_MaxDis*UNIT_L/Const_kpc                     );
      Aux_Message( stdout, "  gasDisk_highResRadius    = %14.7e kpc\n",        gasDisk_highResRadius*UNIT_L/Const_kpc          );
      Aux_Message( stdout, "  gasDisk_lowRes_LEVEL     = %d\n",                gasDisk_lowRes_LEVEL                            );
      Aux_Message( stdout, "  jetSrc_lowRes_LEVEL      = %d\n",                jetSrc_lowRes_LEVEL                             );

      if ( Jet_Ambient == 0 )
      {
      Aux_Message( stdout, "  Amb_UniformDens          = %14.7e g/cm^3\n",     Amb_UniformDens*UNIT_D                          );
      Aux_Message( stdout, "  Amb_UniformTemp          = %14.7e kT/mc**2\n",   Amb_UniformTemp                                 );
      Aux_Message( stdout, "  Amb_UniformVel[x]        = %14.7e c\n",          Amb_UniformVel[0]                               );
      Aux_Message( stdout, "  Amb_UniformVel[y]        = %14.7e c\n",          Amb_UniformVel[1]                               );
      Aux_Message( stdout, "  Amb_UniformVel[z]        = %14.7e c\n",          Amb_UniformVel[2]                               );
      Aux_Message( stdout, "  Jet_UniformNumDens       = %14.7e per cc\n",     Amb_UniformDens*UNIT_D/ParticleMass             );
      } else if ( Jet_Ambient == 2 )
      {
      Aux_Message( stdout, "  IsothermalSlab_Center[0] = %14.7e kpc\n",        IsothermalSlab_Center[0]*UNIT_L/Const_kpc       );
      Aux_Message( stdout, "  IsothermalSlab_Center[1] = %14.7e kpc\n",        IsothermalSlab_Center[1]*UNIT_L/Const_kpc       );
      Aux_Message( stdout, "  IsothermalSlab_Center[2] = %14.7e kpc\n",        IsothermalSlab_Center[2]*UNIT_L/Const_kpc       );
      Aux_Message( stdout, "  criticalTemp             = %14.7e K\n",          criticalTemp / ( Const_kB / (ParticleMass*Const_c*Const_c) ) );
      } // if ( Jet_Ambient == 0 ) ... else if ...

      Aux_Message( stdout, "  CharacteristicSpeed      = %14.7e c\n",          CharacteristicSpeed / UNIT_V                    );
      Aux_Message( stdout, "  Jet_PrecessionAxis[x]    = %14.7e\n",            Jet_PrecessionAxis[0]                           );
      Aux_Message( stdout, "  Jet_PrecessionAxis[y]    = %14.7e\n",            Jet_PrecessionAxis[1]                           );
      Aux_Message( stdout, "  Jet_PrecessionAxis[z]    = %14.7e\n",            Jet_PrecessionAxis[2]                           );

      Aux_Message( stdout, "=============================================================================\n"                   );
   } // if ( MPI_Rank == 0 )

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Setting runtime parameters ... done\n" );

} // FUNCTION : SetParameter



void ReadBinFile( char *FileName, real **buffer )
{

   FILE   *pFile;
   long    lSize;
   size_t  result;

   pFile = fopen( FileName, "rb" );
   if ( pFile == NULL )   Aux_Error( ERROR_INFO, "File error !!\n" );

   // obtain file size
   fseek( pFile, 0, SEEK_END );
   lSize = ftell( pFile );
   rewind( pFile );

   // allocate memory to contain the whole file
   *buffer = (real*)calloc( lSize, sizeof(double) );
   if ( *buffer == NULL )   Aux_Error( ERROR_INFO, "Memory error !!\n" );

   // copy the file into the *buffer
   result = fread( *buffer, 1, lSize, pFile );
   if ( result != lSize )   Aux_Error( ERROR_INFO, "Reading error !!\n" );

   fclose( pFile );

} // FUNCTION : ReadBinFile



void randCloud( real **randXYZ, int numClouds )
{

   *randXYZ = (real*)malloc( 3*numClouds*sizeof(real) );

   for (int i=0; i<3*numClouds; i+=3)
   {
      (*randXYZ)[i+0] = (real)rand() / (real)RAND_MAX * amr->BoxSize[0];
      (*randXYZ)[i+1] = (real)rand() / (real)RAND_MAX * amr->BoxSize[1];
      (*randXYZ)[i+2] = (real)rand() / (real)RAND_MAX * amr->BoxSize[2];

      if ( (*randXYZ)[i+2] < amr->BoxSize[2]*0.5+3.0 )   i-=3;
   } // for (int i=0; i<3*numClouds; i+=3)

} // FUNCTION : randCloud



bool checkInsideClouds( real *randXYZ, int numClouds, real x, real y, real z, real *cloudCenter )
{

   for (int i=0; i<3*numClouds; i+=3)
   {
      const real distance = SQRT( SQR( x-randXYZ[i+0] ) + SQR( y-randXYZ[i+1] ) + SQR( z-randXYZ[i+2] ) );

      if ( distance < radiusCloud )
      {
         cloudCenter[0] = randXYZ[i+0];
         cloudCenter[1] = randXYZ[i+1];
         cloudCenter[2] = randXYZ[i+2];
         return true;
      } // if ( distance < radiusCloud )
   } // for (int i=0; i<3*numClouds; i+=3)

   return false;

} // FUNCTION : checkInsideClouds



// Reading table for interpolations in SetGridIC()
void SetArrayDisk()
{

   char TableFileName[] = "FermiBubble_IC";
   ReadBinFile( TableFileName, &buffer );

   int headerSize = (int)buffer[0];

   Header_disk = (real*)malloc( (size_t)headerSize * sizeof(real) );

   memcpy( Header_disk, buffer, (size_t)headerSize * sizeof(real) );

   ParticleMass = Header_disk[8] * Header_disk[9];

   interfaceHeight  = Header_disk[17];
   interfaceHeight *= Const_kpc / UNIT_L;

   const int Nx = (int)Header_disk[22];
   const int Ny = (int)Header_disk[23];
   const int Nz = (int)Header_disk[24];

   if ( 5*Nx*Ny*Nz > INT_MAX )   Aux_Error( ERROR_INFO, "integer overflow !!\n" );

   real Lx = Header_disk[12];
   real Ly = Header_disk[13];
   real Lz = Header_disk[14];
   const int numGhost = (int)Header_disk[25];

   const int NX = Nx + 2*numGhost;
   const int NY = Ny + 2*numGhost;
   const int NZ = Nz + 2*numGhost;

   if ( OPT__INIT == 1 )
   {
      Rhoo_disk = (real***)calloc_3d_array( (size_t)NX, (size_t)NY, (size_t)NZ, sizeof(real) );
      VelX_disk = (real***)calloc_3d_array( (size_t)NX, (size_t)NY, (size_t)NZ, sizeof(real) );
      VelY_disk = (real***)calloc_3d_array( (size_t)NX, (size_t)NY, (size_t)NZ, sizeof(real) );
      VelZ_disk = (real***)calloc_3d_array( (size_t)NX, (size_t)NY, (size_t)NZ, sizeof(real) );
      Pres_disk = (real***)calloc_3d_array( (size_t)NX, (size_t)NY, (size_t)NZ, sizeof(real) );

      X_disk = (real*)calloc( (size_t)NX, sizeof(real) );
      Y_disk = (real*)calloc( (size_t)NY, sizeof(real) );
      Z_disk = (real*)calloc( (size_t)NZ, sizeof(real) );

      if ( X_disk == NULL )   Aux_Error( ERROR_INFO, "X_disk is NULL at %d, %s !!\n", __LINE__, __FUNCTION__ );
      if ( Y_disk == NULL )   Aux_Error( ERROR_INFO, "Y_disk is NULL at %d, %s !!\n", __LINE__, __FUNCTION__ );
      if ( Z_disk == NULL )   Aux_Error( ERROR_INFO, "Z_disk is NULL at %d, %s !!\n", __LINE__, __FUNCTION__ );

      real *Ptr;
      Ptr = buffer + headerSize;

      for (int c=0; c<5*NX*NY*NZ; c++)
      {
         const int cc = c%(NX*NY*NZ);
         const int i  = (cc - cc%(NY*NZ)) / (NY*NZ);
         const int j  = ((cc - cc%NZ) / NZ) % NY;
         const int k  = cc%NZ;

         if ( 0          <= c && c <   NX*NY*NZ )   Rhoo_disk[i][j][k] = Ptr[c];
         if (   NX*NY*NZ <= c && c < 2*NX*NY*NZ )   VelX_disk[i][j][k] = Ptr[c];
         if ( 2*NX*NY*NZ <= c && c < 3*NX*NY*NZ )   VelY_disk[i][j][k] = Ptr[c];
         if ( 3*NX*NY*NZ <= c && c < 4*NX*NY*NZ )   VelZ_disk[i][j][k] = Ptr[c];
         if ( 4*NX*NY*NZ <= c && c < 5*NX*NY*NZ )   Pres_disk[i][j][k] = Ptr[c];
      } // for (int c=0; c<5*NX*NY*NZ; c++)

      Ptr += 5*NX*NY*NZ;
      for (int c=0; c<NX; c++)   X_disk[c] = Ptr[c];

      Ptr += NX;
      for (int c=0; c<NY; c++)   Y_disk[c] = Ptr[c];

      Ptr += NY;
      for (int c=0; c<NZ; c++)   Z_disk[c] = Ptr[c];
   } // if ( Step == 0 )

   free(buffer);

} // FUNCTION : SetArrayDisk



// Reading table for interpolations in SetGridIC()
void SetArrayHVC()
{

   char TableFileName[] = "FermiBubble_IC_HVC";
   ReadBinFile( TableFileName, &buffer );

   int headerSize = (int)buffer[0];

   Header_hvc = (real*)malloc( (size_t)headerSize*sizeof(real) );

   memcpy( Header_hvc, buffer, (size_t)headerSize*sizeof(real) );

   ParticleMass = Header_hvc[8]*Header_hvc[9];

   interfaceHeight  = Header_hvc[17];
   interfaceHeight *= Const_kpc/UNIT_L;

   const int Nx = (int)Header_hvc[22];
   const int Ny = (int)Header_hvc[23];
   const int Nz = (int)Header_hvc[24];

   if ( 5*Nx*Ny*Nz > INT_MAX )   Aux_Error( ERROR_INFO, "integer overflow !!\n" );

   real Lx = Header_hvc[12];
   real Ly = Header_hvc[13];
   real Lz = Header_hvc[14];
   int numGhost = (int)Header_hvc[25];

   const int NX = Nx+2*numGhost;
   const int NY = Ny+2*numGhost;
   const int NZ = Nz+2*numGhost;

   if ( OPT__INIT == 1 )
   {
      Rhoo_hvc = (real***)calloc_3d_array( (size_t)NX, (size_t)NY, (size_t)NZ, sizeof(real) );
      VelX_hvc = (real***)calloc_3d_array( (size_t)NX, (size_t)NY, (size_t)NZ, sizeof(real) );
      VelY_hvc = (real***)calloc_3d_array( (size_t)NX, (size_t)NY, (size_t)NZ, sizeof(real) );
      VelZ_hvc = (real***)calloc_3d_array( (size_t)NX, (size_t)NY, (size_t)NZ, sizeof(real) );
      Pres_hvc = (real***)calloc_3d_array( (size_t)NX, (size_t)NY, (size_t)NZ, sizeof(real) );

      X_hvc = (real*)calloc( (size_t)NX, sizeof(real) );
      Y_hvc = (real*)calloc( (size_t)NY, sizeof(real) );
      Z_hvc = (real*)calloc( (size_t)NZ, sizeof(real) );

      if ( X_hvc == NULL )   Aux_Error( ERROR_INFO, "X_hvc is NULL at %d, %s !!\n", __LINE__, __FUNCTION__ );
      if ( Y_hvc == NULL )   Aux_Error( ERROR_INFO, "Y_hvc is NULL at %d, %s !!\n", __LINE__, __FUNCTION__ );
      if ( Z_hvc == NULL )   Aux_Error( ERROR_INFO, "Z_hvc is NULL at %d, %s !!\n", __LINE__, __FUNCTION__ );

      real *Ptr;
      Ptr = buffer + headerSize;

      for (int c=0; c<5*NX*NY*NZ; c++)
      {
        const int cc = c%(NX*NY*NZ);
        const int i  = (cc - cc%(NY*NZ)) / (NY*NZ);
        const int j  = ((cc - cc%NZ) / NZ) % NY;
        const int k  = cc%NZ;

        if ( 0          <= c && c <   NX*NY*NZ )   Rhoo_hvc[i][j][k] = Ptr[c];
        if (   NX*NY*NZ <= c && c < 2*NX*NY*NZ )   VelX_hvc[i][j][k] = Ptr[c];
        if ( 2*NX*NY*NZ <= c && c < 3*NX*NY*NZ )   VelY_hvc[i][j][k] = Ptr[c];
        if ( 3*NX*NY*NZ <= c && c < 4*NX*NY*NZ )   VelZ_hvc[i][j][k] = Ptr[c];
        if ( 4*NX*NY*NZ <= c && c < 5*NX*NY*NZ )   Pres_hvc[i][j][k] = Ptr[c];
      } // for (int c=0; c<5*NX*NY*NZ; c++)

      Ptr += 5*NX*NY*NZ;
      for (int c=0; c<NX; c++)   X_hvc[c] = Ptr[c];

      Ptr += NX;
      for (int c=0; c<NY; c++)   Y_hvc[c] = Ptr[c];

      Ptr += NY;
      for (int c=0; c<NZ; c++)   Z_hvc[c] = Ptr[c];
   } // if ( Step == 0 )

   free(buffer);

} // FUNCTION : SetArrayHVC



void Interpolation_UM_IC( real x, real y, real z, real ****Pri_input, real **XYZ, real *Pri_output, bool disk )
{

   real *Header = (disk) ? Header_disk : Header_hvc;

   real xyz[3] = {x, y, z};
   const int Nx = (int)Header[22];
   const int Ny = (int)Header[23];
   const int Nz = (int)Header[24];

   real dx = Header[26];
   real dy = Header[27];
   real dz = Header[28];

   real dxyz[3] = {dx, dy, dz};

   const int numGhost = (int)Header[25];

   const int NX = Nx + 2*numGhost;
   const int NY = Ny + 2*numGhost;
   const int NZ = Nz + 2*numGhost;

   const int Idx = Mis_BinarySearch_Real( XYZ[0], 0, NX-1, x );
   const int Jdx = Mis_BinarySearch_Real( XYZ[1], 0, NY-1, y );
   const int Kdx = Mis_BinarySearch_Real( XYZ[2], 0, NZ-1, z );

   if ( Idx < 0  ||  Idx > NX-2 )   Aux_Error( ERROR_INFO, "x=%e is out of range! XYZ[0][0]=%e, XYZ[0][%d]=%e\n", x, XYZ[0][0], NX-1, XYZ[0][NX-1] );
   if ( Jdx < 0  ||  Jdx > NY-2 )   Aux_Error( ERROR_INFO, "y=%e is out of range! XYZ[1][1]=%e, XYZ[1][%d]=%e\n", y, XYZ[1][0], NY-1, XYZ[1][NY-1] );
   if ( Kdx < 0  ||  Kdx > NZ-2 )   Aux_Error( ERROR_INFO, "z=%e is out of range! XYZ[2][2]=%e, XYZ[2][%d]=%e\n", z, XYZ[2][0], NZ-1, XYZ[2][NZ-1] );

   // TODO: NCOMP_TOTAL to NCOMP_TOTAL_PLUS_MAG
   real Vertex000[NCOMP_TOTAL]={0.0}, Vertex001[NCOMP_TOTAL]={0.0}, Vertex010[NCOMP_TOTAL]={0.0}, Vertex100[NCOMP_TOTAL]={0.0},
        Vertex011[NCOMP_TOTAL]={0.0}, Vertex101[NCOMP_TOTAL]={0.0}, Vertex110[NCOMP_TOTAL]={0.0}, Vertex111[NCOMP_TOTAL]={0.0};

   for (int v=0; v<5; v++)
   {
      Vertex000[v] = Pri_input[v][Idx  ][Jdx  ][Kdx  ];
      Vertex001[v] = Pri_input[v][Idx  ][Jdx  ][Kdx+1];
      Vertex010[v] = Pri_input[v][Idx  ][Jdx+1][Kdx  ];
      Vertex100[v] = Pri_input[v][Idx+1][Jdx  ][Kdx  ];
      Vertex011[v] = Pri_input[v][Idx  ][Jdx+1][Kdx+1];
      Vertex101[v] = Pri_input[v][Idx+1][Jdx  ][Kdx+1];
      Vertex110[v] = Pri_input[v][Idx+1][Jdx+1][Kdx  ];
      Vertex111[v] = Pri_input[v][Idx+1][Jdx+1][Kdx+1];
   }

   bool Unphy = false;
   Unphy |= Hydro_IsUnphysical( UNPHY_MODE_PRIM, Vertex000, NULL_REAL, EoS_DensEint2Pres_CPUPtr, EoS_GuessHTilde_CPUPtr, EoS_HTilde2Temp_CPUPtr, EoS_AuxArray_Flt, EoS_AuxArray_Int, h_EoS_Table, __FILE__,  __LINE__, __FUNCTION__, UNPHY_VERBOSE );
   Unphy |= Hydro_IsUnphysical( UNPHY_MODE_PRIM, Vertex001, NULL_REAL, EoS_DensEint2Pres_CPUPtr, EoS_GuessHTilde_CPUPtr, EoS_HTilde2Temp_CPUPtr, EoS_AuxArray_Flt, EoS_AuxArray_Int, h_EoS_Table, __FILE__,  __LINE__, __FUNCTION__, UNPHY_VERBOSE );
   Unphy |= Hydro_IsUnphysical( UNPHY_MODE_PRIM, Vertex010, NULL_REAL, EoS_DensEint2Pres_CPUPtr, EoS_GuessHTilde_CPUPtr, EoS_HTilde2Temp_CPUPtr, EoS_AuxArray_Flt, EoS_AuxArray_Int, h_EoS_Table, __FILE__,  __LINE__, __FUNCTION__, UNPHY_VERBOSE );
   Unphy |= Hydro_IsUnphysical( UNPHY_MODE_PRIM, Vertex100, NULL_REAL, EoS_DensEint2Pres_CPUPtr, EoS_GuessHTilde_CPUPtr, EoS_HTilde2Temp_CPUPtr, EoS_AuxArray_Flt, EoS_AuxArray_Int, h_EoS_Table, __FILE__,  __LINE__, __FUNCTION__, UNPHY_VERBOSE );
   Unphy |= Hydro_IsUnphysical( UNPHY_MODE_PRIM, Vertex011, NULL_REAL, EoS_DensEint2Pres_CPUPtr, EoS_GuessHTilde_CPUPtr, EoS_HTilde2Temp_CPUPtr, EoS_AuxArray_Flt, EoS_AuxArray_Int, h_EoS_Table, __FILE__,  __LINE__, __FUNCTION__, UNPHY_VERBOSE );
   Unphy |= Hydro_IsUnphysical( UNPHY_MODE_PRIM, Vertex110, NULL_REAL, EoS_DensEint2Pres_CPUPtr, EoS_GuessHTilde_CPUPtr, EoS_HTilde2Temp_CPUPtr, EoS_AuxArray_Flt, EoS_AuxArray_Int, h_EoS_Table, __FILE__,  __LINE__, __FUNCTION__, UNPHY_VERBOSE );
   Unphy |= Hydro_IsUnphysical( UNPHY_MODE_PRIM, Vertex101, NULL_REAL, EoS_DensEint2Pres_CPUPtr, EoS_GuessHTilde_CPUPtr, EoS_HTilde2Temp_CPUPtr, EoS_AuxArray_Flt, EoS_AuxArray_Int, h_EoS_Table, __FILE__,  __LINE__, __FUNCTION__, UNPHY_VERBOSE );
   Unphy |= Hydro_IsUnphysical( UNPHY_MODE_PRIM, Vertex111, NULL_REAL, EoS_DensEint2Pres_CPUPtr, EoS_GuessHTilde_CPUPtr, EoS_HTilde2Temp_CPUPtr, EoS_AuxArray_Flt, EoS_AuxArray_Int, h_EoS_Table, __FILE__,  __LINE__, __FUNCTION__, UNPHY_VERBOSE );

   if ( Unphy )   Aux_Error( ERROR_INFO, "Idx=%d, Jdx=%d, Kdx=%d, x=%e, y=%e, z=%e is unphysical !!\n", Idx, Jdx, Kdx, x, y, z );

   real xyz000[3] = { XYZ[0][Idx], XYZ[1][Jdx], XYZ[2][Kdx] };

   for (int v=0; v<5; v++)
   {
     real FieldAtVertices[8] = { Vertex000[v], Vertex001[v], Vertex010[v], Vertex100[v], Vertex011[v], Vertex101[v], Vertex110[v], Vertex111[v] };

     Pri_output[v] = TrilinearInterpolation( FieldAtVertices, xyz000, dxyz, xyz );
   }

} // FUNCTION : Interpolation_UM_IC



#ifdef GRAVITY
//-------------------------------------------------------------------------------------------------------
// Function    :  IsothermalSlab_Pot
// Description :
//
// Note        :
//
// Parameter   :
//
// Return      :
//-------------------------------------------------------------------------------------------------------
real IsothermalSlab_Pot( const real z )
{

   real Pot, Log;

// 1. isothermal slab
   Pot  = 2.0 * M_PI * NEWTON_G * IsothermalSlab_PeakDens;
   Pot /= SQR( IsothermalSlab_VelocityDispersion );
   Pot  = log( cosh( z*sqrt(Pot) ) );
   Pot *= 2.0 * SQR( IsothermalSlab_VelocityDispersion );

// 2. log potential
   Log  = SQR(v_halo) * log( z*z + SQR(distance_h) );
   Log -= SQR(v_halo) * log( SQR(distance_h) );

   return Pot + Log;

} // FUNCTION : IsothermalSlab_Pot
#endif



//-------------------------------------------------------------------------------------------------------
// Function    :  SetGridIC
// Description :  Set the problem-specific initial condition on grids
//
// Note        :  1. This function may also be used to estimate the numerical errors when OPT__OUTPUT_USER is enabled
//                   --> In this case, it should provide the analytical solution at the given "Time"
//                2. This function will be invoked by multiple OpenMP threads when OPENMP is enabled
//                   --> Please ensure that everything here is thread-safe
//
// Parameter   :  fluid    : Fluid field to be initialized
//                x/y/z    : Physical coordinates
//                Time     : Physical time
//                lv       : Target refinement level
//                AuxArray : Auxiliary array
//
// Return      :  fluid
//-------------------------------------------------------------------------------------------------------
void SetGridIC( real fluid[], const double x, const double y, const double z, const double Time,
                const int lv, double AuxArray[] )
{

// variables for jet
   const real xc = x - IsothermalSlab_Center[0];
   const real yc = y - IsothermalSlab_Center[1];
   const real zc = z - IsothermalSlab_Center[2];
   real Pri[NCOMP_TOTAL] = {0.0};

   real ***Pri_disk_input[5] = { Rhoo_disk, VelX_disk, VelY_disk, VelZ_disk, Pres_disk };
   real *** Pri_hvc_input[5] = { Rhoo_hvc,   VelX_hvc,  VelY_hvc,  VelZ_hvc,  Pres_hvc };

   if ( Jet_Ambient == 0 ) // uniform ambient
   {
      Pri[0] = (real)Amb_UniformDens;
      Pri[1] = (real)Amb_UniformVel[0];
      Pri[2] = (real)Amb_UniformVel[1];
      Pri[3] = (real)Amb_UniformVel[2];
      Pri[4] = (real)Amb_UniformTemp * Amb_UniformDens;

      Hydro_Pri2Con( Pri, fluid, false, PassiveNorm_NVar, PassiveNorm_VarIdx, EoS_DensPres2Eint_CPUPtr,
                     EoS_Temp2HTilde_CPUPtr, EoS_HTilde2Temp_CPUPtr, EoS_AuxArray_Flt, EoS_AuxArray_Int,
                     h_EoS_Table, NULL );

      if ( Hydro_IsUnphysical( UNPHY_MODE_PRIM, Pri, NULL_REAL, EoS_DensEint2Pres_CPUPtr, EoS_GuessHTilde_CPUPtr,
                               EoS_HTilde2Temp_CPUPtr, EoS_AuxArray_Flt, EoS_AuxArray_Int, h_EoS_Table,
                               __FILE__,  __LINE__, __FUNCTION__, UNPHY_VERBOSE ) )
         Aux_Error( ERROR_INFO, "Unphysical cell at (%e, %e, %e)\n", x, y, z );

#     if ( NCOMP_PASSIVE_USER > 0 )
      fluid[Passive_0000] = fluid[DENS];
      fluid[Passive_0001] = 0.0;
      fluid[Passive_0002] = 0.0;
#     endif

#     ifdef COSMIC_RAY
      fluid[CRAY] = Amb_CR_Engy * Jet_SrcGamma;
#     endif
   }
   else if ( Jet_Ambient == 2 ) // cold disk in stratified ambient
   {
#     ifdef GRAVITY
      if ( fabs(zc) < interfaceHeight )
      {
         real *XYZ[3] = { X_disk, Y_disk, Z_disk };

         Interpolation_UM_IC( xc, yc, zc, Pri_disk_input, XYZ, Pri, true );

         if ( Hydro_IsUnphysical( UNPHY_MODE_PRIM, Pri, NULL_REAL, EoS_DensEint2Pres_CPUPtr, EoS_GuessHTilde_CPUPtr,
                                  EoS_HTilde2Temp_CPUPtr, EoS_AuxArray_Flt, EoS_AuxArray_Int, h_EoS_Table,
                                  __FILE__,  __LINE__, __FUNCTION__, UNPHY_VERBOSE ) )
            Aux_Error( ERROR_INFO, "Unphysical cell at (%e, %e %e)\n", x, y, z );

         if ( Pri[4]/Pri[0] > criticalTemp )   Pri[0] = Pri[4] / ambientTemperature;

         Hydro_Pri2Con( Pri, fluid, false, PassiveNorm_NVar, PassiveNorm_VarIdx, EoS_DensPres2Eint_CPUPtr,
                        EoS_Temp2HTilde_CPUPtr, EoS_HTilde2Temp_CPUPtr, EoS_AuxArray_Flt, EoS_AuxArray_Int,
                        h_EoS_Table, NULL );

#        if ( NCOMP_PASSIVE_USER > 0 )
         fluid[Passive_0000] = fluid[DENS];
         fluid[Passive_0001] = 0.0;
         fluid[Passive_0002] = 0.0;
#        endif

#        ifdef COSMIC_RAY
         fluid[CRAY] = Amb_CR_Engy * Jet_SrcGamma;
#        endif
      }
      else // if ( fabs(zc) < interfaceHeight )
      {
         real Dens_gDisk_ambient, PotAtZ0, ambientDens;

         PotAtZ0 = IsothermalSlab_Pot( interfaceHeight );

         Dens_gDisk_ambient  = ambientTemperature / gasDiskTemperature;
         Dens_gDisk_ambient *= exp( PotAtZ0 * (ambientTemperature-gasDiskTemperature) / (ambientTemperature*gasDiskTemperature) );

         if ( Dens_gDisk_ambient > HUGE_NUMBER  ||  Dens_gDisk_ambient < -HUGE_NUMBER )
            Aux_Error( ERROR_INFO, "(Dens_gDisk_ambient = %e) not in [-HUGE_NUMBER, HUGE_NUMBER] !! %s: %d\n", Dens_gDisk_ambient, __FUNCTION__, __LINE__ );

         real ambientPeakDens  = gasDiskPeakDens / Dens_gDisk_ambient;

         if ( fabs(zc) > IsothermalSlab_Truncation )
            ambientDens = -IsothermalSlab_Pot(IsothermalSlab_Truncation) / ambientTemperature;
         else
            ambientDens = -IsothermalSlab_Pot(zc) / ambientTemperature;

         ambientDens  = exp(ambientDens);
         ambientDens *= ambientPeakDens;

         Pri[0] = ambientDens;
         Pri[1] = 0.0;
         Pri[2] = 0.0;
         Pri[3] = 0.0;
         Pri[4] = ambientDens * ambientTemperature;

         real cloudCenter[3];

         if ( checkInsideClouds( randXYZ, numClouds, x, y, z, cloudCenter ) )
         {
            real Pri_hvc_output[5];

            real *XYZ[3] =  { X_hvc, Y_hvc, Z_hvc };

            Interpolation_UM_IC( x-cloudCenter[0], y-cloudCenter[1], z-cloudCenter[2],
                                 Pri_hvc_input, XYZ, Pri_hvc_output, false );

            if ( Hydro_IsUnphysical( UNPHY_MODE_PRIM, Pri_hvc_output, NULL_REAL, EoS_DensEint2Pres_CPUPtr,
                                     EoS_GuessHTilde_CPUPtr, EoS_HTilde2Temp_CPUPtr, EoS_AuxArray_Flt,
                                     EoS_AuxArray_Int, h_EoS_Table, __FILE__,  __LINE__, __FUNCTION__, UNPHY_VERBOSE ) )
               Aux_Error( ERROR_INFO, "Unphysical cell at (%e, %e %e)\n", x, y, z );

            Pri[0] = Pri_hvc_output[0]*0.05*Const_mp/UNIT_D;
            Pri[1] = 0.0;
            Pri[2] = 0.0;
            Pri[3] = 0.0;
            Pri[4] = Pri[0];
         }

         Hydro_Pri2Con( Pri, fluid, false, PassiveNorm_NVar, PassiveNorm_VarIdx, EoS_DensPres2Eint_CPUPtr,
                        EoS_Temp2HTilde_CPUPtr, EoS_HTilde2Temp_CPUPtr, EoS_AuxArray_Flt, EoS_AuxArray_Int,
                        h_EoS_Table, NULL );

         if ( Hydro_IsUnphysical( UNPHY_MODE_PRIM, Pri, NULL_REAL, EoS_DensEint2Pres_CPUPtr, EoS_GuessHTilde_CPUPtr,
                                  EoS_HTilde2Temp_CPUPtr, EoS_AuxArray_Flt, EoS_AuxArray_Int, h_EoS_Table,
                                  __FILE__,  __LINE__, __FUNCTION__, UNPHY_VERBOSE ) )
            Aux_Error( ERROR_INFO, "Unphysical cell at (%e, %e %e)\n", x, y, z );

#        if ( NCOMP_PASSIVE_USER > 0 )
         fluid[Passive_0000] = 0.0;
         fluid[Passive_0001] = 0.0;
         fluid[Passive_0002] = 0.0;
#        endif

#        ifdef COSMIC_RAY
         fluid[CRAY] = Amb_CR_Engy * Jet_SrcGamma;
#        endif
      } // if ( fabs(zc) < interfaceHeight ) ... else ...
#  endif
   }
   else if ( Jet_Ambient == 3 ) // cold disk in uniform ambient
   {
#     ifdef GRAVITY
      if ( fabs(zc) < interfaceHeight )
      {
         real *XYZ[3] =  { X_disk, Y_disk, Z_disk };

         Interpolation_UM_IC( xc, yc, zc, Pri_disk_input, XYZ, Pri, true );

         if ( Hydro_IsUnphysical( UNPHY_MODE_PRIM, Pri, NULL_REAL, EoS_DensEint2Pres_CPUPtr, EoS_GuessHTilde_CPUPtr,
                                  EoS_HTilde2Temp_CPUPtr, EoS_AuxArray_Flt, EoS_AuxArray_Int, h_EoS_Table,
                                  __FILE__,  __LINE__, __FUNCTION__, UNPHY_VERBOSE ) )
            Aux_Error( ERROR_INFO, "Unphysical cell at (%e, %e, %e)\n", x, y, z );

         Hydro_Pri2Con( Pri, fluid, false, PassiveNorm_NVar, PassiveNorm_VarIdx, EoS_DensPres2Eint_CPUPtr,
                        EoS_Temp2HTilde_CPUPtr, EoS_HTilde2Temp_CPUPtr, EoS_AuxArray_Flt, EoS_AuxArray_Int,
                        h_EoS_Table, NULL );

#        if ( NCOMP_PASSIVE_USER > 0 )
         fluid[Passive_0000] = 0.0;
         fluid[Passive_0001] = 0.0;
         fluid[Passive_0002] = 0.0;
#        endif

#        ifdef COSMIC_RAY
         fluid[CRAY] = Amb_CR_Engy * Jet_SrcGamma;
#        endif
      }
      else // if ( fabs(zc) < interfaceHeight )
      {
         real Dens_gDisk_ambient, PotAtZ0, ambientDens;

         PotAtZ0 = IsothermalSlab_Pot(interfaceHeight);

         Dens_gDisk_ambient  = ambientTemperature / gasDiskTemperature;
         Dens_gDisk_ambient *= exp( PotAtZ0*(ambientTemperature-gasDiskTemperature)/(ambientTemperature*gasDiskTemperature) );

         if ( Dens_gDisk_ambient > HUGE_NUMBER  ||  Dens_gDisk_ambient < -HUGE_NUMBER )
            Aux_Error( ERROR_INFO, "(Dens_gDisk_ambient = %e) not in [-HUGE_NUMBER, HUGE_NUMBER] !! %s: %d\n", Dens_gDisk_ambient, __FUNCTION__, __LINE__ );

         real ambientPeakDens = gasDiskPeakDens / Dens_gDisk_ambient;

         ambientDens  = -IsothermalSlab_Pot(interfaceHeight) / ambientTemperature;
         ambientDens  = exp(ambientDens);
         ambientDens *= ambientPeakDens;

         Pri[0] = ambientDens;
         Pri[1] = 0.0;
         Pri[2] = 0.0;
         Pri[3] = 0.0;
         Pri[4] = ambientDens * ambientTemperature;

         if ( Hydro_IsUnphysical( UNPHY_MODE_PRIM, Pri, NULL_REAL, EoS_DensEint2Pres_CPUPtr, EoS_GuessHTilde_CPUPtr,
                                  EoS_HTilde2Temp_CPUPtr, EoS_AuxArray_Flt, EoS_AuxArray_Int, h_EoS_Table,
                                  __FILE__,  __LINE__, __FUNCTION__, UNPHY_VERBOSE ) )
            Aux_Error( ERROR_INFO, "Unphysical cell at (%e, %e, %e)\n", x, y, z );

         Hydro_Pri2Con( Pri, fluid, false, PassiveNorm_NVar, PassiveNorm_VarIdx, EoS_DensPres2Eint_CPUPtr,
                        EoS_Temp2HTilde_CPUPtr, EoS_HTilde2Temp_CPUPtr, EoS_AuxArray_Flt, EoS_AuxArray_Int,
                        h_EoS_Table, NULL );

#        if ( NCOMP_PASSIVE_USER > 0 )
         fluid[Passive_0000] = 0.0;
         fluid[Passive_0001] = 0.0;
         fluid[Passive_0002] = 0.0;
#        endif

#        ifdef COSMIC_RAY
         fluid[CRAY] = Amb_CR_Engy * Jet_SrcGamma;
#        endif

      } // if ( fabs(zc) < interfaceHeight ) ... else ...
   } // else if ( Jet_Ambient == 3 )
   else if ( Jet_Ambient == 4 )
   {
      real ambientPeakDens = (real)2.842783e-27 / UNIT_D;
      real ambientDens;

      if ( fabs(zc) > IsothermalSlab_Truncation )
         ambientDens = -IsothermalSlab_Pot(IsothermalSlab_Truncation)/ambientTemperature;
      else
         ambientDens = -IsothermalSlab_Pot(zc)/ambientTemperature;

      ambientDens = ambientPeakDens * exp(ambientDens);

      Pri[0] = ambientDens;
      Pri[1] = 0.0;
      Pri[2] = 0.0;
      Pri[3] = 0.0;
      Pri[4] = ambientDens*ambientTemperature;

      Hydro_Pri2Con( Pri, fluid, false, PassiveNorm_NVar, PassiveNorm_VarIdx, EoS_DensPres2Eint_CPUPtr,
                     EoS_Temp2HTilde_CPUPtr, EoS_HTilde2Temp_CPUPtr, EoS_AuxArray_Flt, EoS_AuxArray_Int,
                     h_EoS_Table, NULL );

      if ( Hydro_IsUnphysical( UNPHY_MODE_PRIM, Pri, NULL_REAL, EoS_DensEint2Pres_CPUPtr, EoS_GuessHTilde_CPUPtr,
                               EoS_HTilde2Temp_CPUPtr, EoS_AuxArray_Flt, EoS_AuxArray_Int, h_EoS_Table,
                               __FILE__,  __LINE__, __FUNCTION__, UNPHY_VERBOSE ) )
         Aux_Error( ERROR_INFO, "Unphysical cell at (%e, %e, %e)\n", x, y, z );

#     if ( NCOMP_PASSIVE_USER > 0 )
      fluid[Passive_0000] = 0.0;
      fluid[Passive_0001] = 0.0;
      fluid[Passive_0002] = 0.0;
#     endif

#     ifdef COSMIC_RAY
      fluid[CRAY] = Amb_CR_Engy * Jet_SrcGamma;
#     endif

#     endif // #ifdef GRAVITY
   } // if ( Jet_Ambient == 0 ) ... else if

} // FUNCTION : SetGridIC



//-------------------------------------------------------------------------------------------------------
// Function    :  Flu_ResetByUser_FermiBubble
// Description :  Function to reset the fluid field
//
// Note        :  1. Invoked by "Flu_ResetByUser_API()" and "Model_Init_ByFunction_AssignData()" using the
//                   function pointer "Flu_ResetByUser_Func_Ptr"
//                2. This function will be invoked when constructing the initial condition
//                    (by calling "Model_Init_ByFunction_AssignData()") and after each update
//                    (by calling "Flu_ResetByUser_API()")
//                3. Input "fluid" array stores the original values
//                4. Even when DUAL_ENERGY is adopted, one does NOT need to set the dual-energy variable here
//                   --> It will be set automatically in "Flu_ResetByUser_API()" and "Model_Init_ByFunction_AssignData()"
//                5. Enabled by the runtime option "OPT__RESET_FLUID"
//
// Parameter   :  fluid    : Fluid array storing both the input (origial) and reset values
//                           --> Including both active and passive variables
//                x/y/z    : Target physical coordinates
//                Time     : Target physical time
//                lv       : Target refinement level
//                AuxArray : Auxiliary array
//
// Return      :  true  : This cell has been reset
//                false : This cell has not been reset
//
// =======================================================================================
/*        G       A       C              */
/*          ____________                 */
/*          \     |     /                */
/*           \   E|    /        z        */
/*            \  /|\  /         ^        */
/*             \/_|_\/          |        */
/*             /\O| /\B         |        */
/*            /  \|/  \                  */
/*           /   D|    \                 */
/*          /_____|_____\                */
/*                F                      */
// =======================================================================================
int Flu_ResetByUser_FermiBubble( real fluid[], const double Emag, const double x, const double y, const double z,
                                 const double Time, const double dt, const int lv, double AuxArray[] )
{

   if ( Jet_Fire == 0 )   return false;

   if ( Jet_Duration < Time )   return false;

   if ( !Jet_SphericalSrc )
   {
      double xp[3], rp[3];
      double Prim[NCOMP_TOTAL] = {0.0}, Cons[NCOMP_TOTAL] = {0.0}, Vel[3];
      real PriReal[NCOMP_TOTAL] = {0.0};
      double PrecessionAxis_Spherical[3], Omega_t;
      bool InsideUpperCone, InsideLowerCone;
      double Jet_SrcVelSmooth;

      Omega_t = Jet_AngularVelocity * Time * M_PI / 180.0;

//    shift the coordinate origin to the source center (the point O)
      xp[0] = x - Jet_Center[0];
      xp[1] = y - Jet_Center[1];
      xp[2] = z - Jet_Center[2];

      if ( Jet_PrecessionAxis[0] != 0.0  ||  Jet_PrecessionAxis[1] != 0.0  ||  Jet_PrecessionAxis[2] == 0.0 )
      {
//       get theta, phi for the first rotation
         Mis_Cartesian2Spherical( Jet_PrecessionAxis, PrecessionAxis_Spherical );
//       rotate coordinate to align z-axis with fixed precession axis
         CartesianRotate( xp, PrecessionAxis_Spherical[1], PrecessionAxis_Spherical[2], false );
      }

//    rotate coordinate to align z-axis with rotating symmetric axis
      CartesianRotate( xp, Jet_PrecessionAngle, Omega_t, false );

//    determine whether or not the point is inside of source
      InsideUpperCone  = SQR(xp[0]) + SQR(xp[1]) <= SQR( +tan(Jet_HalfOpeningAngle)*xp[2] + Jet_Radius );
      InsideUpperCone &= 0.0 <= xp[2] && xp[2] <= Jet_HalfHeight;

      InsideLowerCone  = SQR(xp[0]) + SQR(xp[1]) <= SQR( -tan(Jet_HalfOpeningAngle)*xp[2] + Jet_Radius );
      InsideLowerCone &= -Jet_HalfHeight <= xp[2] && xp[2] <= 0.0;

      if ( Jet_HalfOpeningAngle != 0.0 )
      {
         InsideUpperCone &= SQR(xp[0]) + SQR(xp[1]) + SQR(xp[2] + Jet_Radius/tan(Jet_HalfOpeningAngle))
                         <= SQR(Jet_HalfHeight+Jet_Radius/tan(Jet_HalfOpeningAngle));

         InsideUpperCone &= 0.0 <= xp[2];

         InsideLowerCone &= SQR(xp[0]) + SQR(xp[1]) + SQR(xp[2] - Jet_Radius/tan(Jet_HalfOpeningAngle))
                         <= SQR(Jet_HalfHeight+Jet_Radius/tan(Jet_HalfOpeningAngle));

         InsideLowerCone &= xp[2] <= 0.0;
      }
      else
      {
         InsideUpperCone &= 0.0 <= xp[2]  &&  xp[2] <= Jet_HalfHeight;
         InsideLowerCone &= -Jet_HalfHeight <= xp[2]  &&  xp[2] <= 0.0;
      } // if ( Jet_HalfOpeningAngle != 0.0 ) ... else ...

//    set fluid variable inside source
      if ( ( InsideUpperCone  &&  ( Jet_Fire == 1  ||  Jet_Fire == 3 ) )  ||
           ( InsideLowerCone  &&  ( Jet_Fire == 2  ||  Jet_Fire == 3 ) ) )
      {
         if ( Jet_HalfOpeningAngle == 0.0 )
         {
            Vel[0] = 0.0;
            Vel[1] = 0.0;
            if ( InsideUpperCone == true )   Vel[2] = +Jet_SrcVel;
            else                             Vel[2] = -Jet_SrcVel;

            CartesianRotate( Vel, Jet_PrecessionAngle, Omega_t, true );

            if ( Jet_PrecessionAxis[0] != 0.0  ||  Jet_PrecessionAxis[1] != 0.0  ||  Jet_PrecessionAxis[2] == 0.0 )
               CartesianRotate( Vel, PrecessionAxis_Spherical[1], PrecessionAxis_Spherical[2], true );

            Prim[0] = Jet_SrcDens;
            Prim[1] = Vel[0];
            Prim[2] = Vel[1];
            Prim[3] = Vel[2];
            Prim[4] = Jet_SrcTemp*Jet_SrcDens;
         }
         else // if ( Jet_HalfOpeningAngle == 0.0 )
         {
//          shift origin to the point D/E
            if ( InsideUpperCone == true )   xp[2] += Jet_Radius/tan(Jet_HalfOpeningAngle);
            else                             xp[2] -= Jet_Radius/tan(Jet_HalfOpeningAngle);

            CartesianRotate( xp, Jet_PrecessionAngle, Omega_t, true );

            Mis_Cartesian2Spherical( xp, rp );

            if ( InsideLowerCone == true )   rp[1] -= M_PI;

//          smooth velocity on cross section
            if ( Jet_SmoothVel )   Jet_SrcVelSmooth = Jet_SrcVel*SQR(cos( 0.5 * M_PI * rp[1] / Jet_HalfOpeningAngle ));
            else                   Jet_SrcVelSmooth = Jet_SrcVel;

            if ( Jet_PrecessionAxis[0] != 0.0  ||  Jet_PrecessionAxis[1] != 0.0  ||  Jet_PrecessionAxis[2] == 0.0 )
               CartesianRotate( xp, PrecessionAxis_Spherical[1], PrecessionAxis_Spherical[2], true );

            Mis_Cartesian2Spherical( xp, rp );

            Prim[0] = Jet_SrcDens;
            Prim[1] = Jet_SrcVelSmooth*sin(rp[1])*cos(rp[2]);
            Prim[2] = Jet_SrcVelSmooth*sin(rp[1])*sin(rp[2]);
            Prim[3] = Jet_SrcVelSmooth*cos(rp[1]);
            Prim[4] = Jet_SrcTemp*Jet_SrcDens;
         } // if ( Jet_HalfOpeningAngle == 0.0 ) ... else ...

         PriReal[0] = (real)Prim[0];
         PriReal[1] = (real)Prim[1];
         PriReal[2] = (real)Prim[2];
         PriReal[3] = (real)Prim[3];
         PriReal[4] = (real)Prim[4];

         Hydro_Pri2Con( PriReal, fluid, NULL_BOOL, NULL_INT, NULL, EoS_DensPres2Eint_CPUPtr,
                        EoS_Temp2HTilde_CPUPtr, EoS_HTilde2Temp_CPUPtr,
                        EoS_AuxArray_Flt, EoS_AuxArray_Int, h_EoS_Table, NULL );

#        if ( NCOMP_PASSIVE_USER > 0 )
         fluid[Passive_0000] = 0.0;
         fluid[Passive_0001] = fluid[DENS];
         fluid[Passive_0002] = Jet_Src_CR_Engy * Jet_SrcGamma;
#        endif

#        ifdef COSMIC_RAY
         fluid[CRAY] = Jet_Src_CR_Engy * Jet_SrcGamma;
#        endif

         return true;
      } // if ( ( InsideUpperCone  &&  ( Jet_Fire == 1  ||  Jet_Fire == 3 ) )  ||
//              ( InsideLowerCone  &&  ( Jet_Fire == 2  ||  Jet_Fire == 3 ) ) )
   }
   else // if ( !Jet_SphericalSrc )
   {
      double xp[3], rp[3];
      double Prim[NCOMP_TOTAL] = {0.0}, Cons[NCOMP_TOTAL] = {0.0}, Vel[3];
      real PriReal[NCOMP_TOTAL] = {0.0};

//    shift the coordinate origin to the source center (the point O)
      xp[0] = x - Jet_Center[0];
      xp[1] = y - Jet_Center[1];
      xp[2] = z - Jet_Center[2];

      //Mis_Cartesian2Spherical(xp, rp);

      const double R = SQRT( SQR(xp[0]) + SQR(xp[1]) + SQR(xp[2]) );

      if ( R < Jet_HalfHeight )
      {
         Prim[0] = Jet_SrcDens;
         Prim[1] = Jet_SrcVel*xp[0]/R;
         Prim[2] = Jet_SrcVel*xp[1]/R;
         Prim[3] = Jet_SrcVel*xp[2]/R;
         Prim[4] = Jet_SrcTemp*Jet_SrcDens;

         PriReal[0] = (real)Prim[0];
         PriReal[1] = (real)Prim[1];
         PriReal[2] = (real)Prim[2];
         PriReal[3] = (real)Prim[3];
         PriReal[4] = (real)Prim[4];

         Hydro_Pri2Con( PriReal, fluid, NULL_BOOL, NULL_INT, NULL, EoS_DensPres2Eint_CPUPtr,
                        EoS_Temp2HTilde_CPUPtr, EoS_HTilde2Temp_CPUPtr, EoS_AuxArray_Flt,
                        EoS_AuxArray_Int, h_EoS_Table, NULL );

#        if ( NCOMP_PASSIVE_USER > 0 )
         fluid[Passive_0000] = 0.0;
         fluid[Passive_0001] = fluid[DENS];
         fluid[Passive_0002] = Jet_Src_CR_Engy * Jet_SrcGamma;
#        endif

#        ifdef COSMIC_RAY
         fluid[CRAY] = Jet_Src_CR_Engy * Jet_SrcGamma;
#        endif

         return true;
      } // if ( R < Jet_HalfHeight )
   } // if ( !Jet_SphericalSrc ) ... else ...

   return false;

} // FUNCTION : Flu_ResetByUser_FermiBubble



//-------------------------------------------------------------------------------------------------------
// Function    :  Flag_Region_FermiBubble
// Description :
//
// Note        :  1. Invoked by Flag_Check() using the function pointer "Flag_Region_Ptr",
//                   which must be set by a test problem initializer
//                2. Enabled by the runtime option "OPT__FLAG_REGION"
//
// Parameter   :  i,j,k       : Indices of the target element in the patch ptr[0][lv][PID]
//                lv          : Refinement level of the target patch
//                PID         : ID of the target patch
//
// Return      :  "true/false"  if the input cell "is/is not" within the region allowed for refinement
//-------------------------------------------------------------------------------------------------------
static bool Flag_Region_FermiBubble( const int i, const int j, const int k, const int lv, const int PID )
{

   if ( Step <= 0 )   return true;
   const double dh        = amr->dh[lv]; // grid size
   const double Pos[3]    = { amr->patch[0][lv][PID]->EdgeL[0] + (i+0.5)*dh,  // x,y,z position
                              amr->patch[0][lv][PID]->EdgeL[1] + (j+0.5)*dh,
                              amr->patch[0][lv][PID]->EdgeL[2] + (k+0.5)*dh  };
   const double Center[3] = { 0.5*amr->BoxSize[0],
                              0.5*amr->BoxSize[1],
                              0.5*amr->BoxSize[2] };
   const double dr[3]     = { Pos[0]-Center[0]-Jet_CenOffset[0],
                              Pos[1]-Center[1]-Jet_CenOffset[1],
                              Pos[2]-Center[2]-Jet_CenOffset[2] };
   const double R         = sqrt( SQR(dr[0]) + SQR(dr[1]) );

   const bool   Flag      = R > gasDisk_highResRadius  &&  lv > gasDisk_lowRes_LEVEL  &&  fabs(dr[2]) < 2.0*interfaceHeight;

   if ( Flag )   return false;

   return true;

} // FUNCTION : Flag_Region_FermiBubble



//-------------------------------------------------------------------------------------------------------
// Function    :  Flag_FermiBubble
// Description :  Refine the area near the jet source or the disk
//
// Note        :  1. Invoked by Flag_Check() using the function pointer "Flag_User_Ptr",
//                   which must be set by a test problem initializer
//                2. Enabled by the runtime option "OPT__FLAG_USER"
//
// Parameter   :  i,j,k     : Indices of the target element in the patch ptr[ amr->FluSg[lv] ][lv][PID]
//                lv        : Refinement level of the target patch
//                PID       : ID of the target patch
//                Threshold : User-provided threshold for the flag operation, which is loaded from the
//                            file "Input__Flag_User"
//
// Return      :  "true"  if the flag criteria are satisfied
//                "false" if the flag criteria are not satisfied
//-------------------------------------------------------------------------------------------------------
bool Flag_FermiBubble( const int i, const int j, const int k, const int lv, const int PID, const double *Threshold )
{

   const double dh        = amr->dh[lv];                                                  // grid size
   const double Pos[3]    = { amr->patch[0][lv][PID]->EdgeL[0] + (i+0.5)*dh,              // x,y,z position
                              amr->patch[0][lv][PID]->EdgeL[1] + (j+0.5)*dh,
                              amr->patch[0][lv][PID]->EdgeL[2] + (k+0.5)*dh  };

   const double Center[3] = { Jet_Center[0], Jet_Center[1], Jet_Center[2] };

   const double dR[3]     = { Pos[0]-Center[0], Pos[1]-Center[1], Pos[2]-Center[2] };
   const double R         = sqrt( SQR(dR[0]) + SQR(dR[1]) + SQR(dR[2]) );

   bool Flag, Src = R <= dh*1.8;
   if ( Jet_Ambient != 4 )
   {
      bool Disk = fabs(dR[2]) <= dh*1.8;
      if ( lv >= jetSrc_lowRes_LEVEL )   Disk = false;
      Flag = Src || Disk;
   }
   else
   {
      Flag = Src;
   }
   return Flag;

} // FUNCTION : Flag_FermiBubble



//-------------------------------------------------------------------------------------------------------
// Function    :  CartesianRotate
// Description :
//
// Note        :
//
// Parameter   :
//
// Return      :
//-------------------------------------------------------------------------------------------------------
void CartesianRotate( double x[], double theta, double phi, bool inverse )
{

   double xp[3];

   if ( inverse )
   {
      xp[0] = -            sin(phi)*x[0] - cos(theta)*cos(phi)*x[1] + sin(theta)*cos(phi)*x[2];
      xp[1] = +            cos(phi)*x[0] - cos(theta)*sin(phi)*x[1] + sin(theta)*sin(phi)*x[2];
      xp[2] =                            + sin(theta)*         x[1] + cos(theta)*         x[2];
   } else
   {
      xp[0] = -            sin(phi)*x[0] +            cos(phi)*x[1];
      xp[1] = - cos(theta)*cos(phi)*x[0] - cos(theta)*sin(phi)*x[1] + sin(theta)*         x[2];
      xp[2] = + sin(theta)*cos(phi)*x[0] + sin(theta)*sin(phi)*x[1] + cos(theta)*         x[2];
   }

   for (int i=0; i<3; i++)   x[i] = xp[i];

} // FUNCTION : CartesianRotate



//-------------------------------------------------------------------------------------------------------
// Function    :  AddNewField_FermiBubble
// Description :  Add the problem-specific fields
//
// Note        :  1. Ref: https://github.com/gamer-project/gamer/wiki/Adding-New-Simulations#v-add-problem-specific-grid-fields-and-particle-attributes
//                2. Invoke AddField() for each of the problem-specific field:
//                   --> Field label sent to AddField() will be used as the output name of the field
//                   --> Field index returned by AddField() can be used to access the field data
//                3. Pre-declared field indices are put in Field.h
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void AddNewField_FermiBubble()
{

#  if ( NCOMP_PASSIVE_USER > 0 )
   if ( Passive_0000 == 5 )   Passive_0000 = AddField( "Passive_0000", FIXUP_FLUX_YES, FIXUP_REST_YES, NORMALIZE_NO, INTERP_FRAC_NO );
   if ( Passive_0001 == 6 )   Passive_0001 = AddField( "Passive_0001", FIXUP_FLUX_YES, FIXUP_REST_YES, NORMALIZE_NO, INTERP_FRAC_NO );
   if ( Passive_0002 == 7 )   Passive_0002 = AddField( "Passive_0002", FIXUP_FLUX_YES, FIXUP_REST_YES, NORMALIZE_NO, INTERP_FRAC_NO );
#  endif

} // FUNCTION : AddNewField_Jet
#endif // #if ( MODEL == HYDRO )



//-------------------------------------------------------------------------------------------------------
// Function    :  Init_TestProb_Hydro_FermiBubble
// Description :  Test problem initializer
//
// Note        :  None
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Init_TestProb_Hydro_FermiBubble()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );


// validate the compilation flags and runtime parameters
   Validate();


#  if ( MODEL == HYDRO )
// set the problem-specific runtime parameters
   SetParameter();


   Init_Function_User_Ptr    = SetGridIC;
   Flag_User_Ptr             = Flag_FermiBubble;
   Flag_Region_Ptr           = Flag_Region_FermiBubble;
   Flu_ResetByUser_Func_Ptr  = Flu_ResetByUser_FermiBubble;
#  ifdef GRAVITY
   Init_ExtPot_Ptr           = Init_ExtPot_IsothermalSlab;
#  endif

#  if ( NCOMP_PASSIVE_USER > 0 )
   Init_Field_User_Ptr       = AddNewField_FermiBubble;
#  endif
#  ifdef SUPPORT_HDF5
   Output_HDF5_InputTest_Ptr = LoadInputTestProb;
#  endif
#  endif // #if ( MODEL == HYDRO )


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : Init_TestProb_Hydro_FermiBubble
