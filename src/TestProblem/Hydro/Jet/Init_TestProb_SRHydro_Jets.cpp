#include <random>
#include <limits.h>
#include <math.h>
#include "GAMER.h"
#include "TestProb.h"

void ***calloc_3d_array (size_t nt, size_t nr, size_t nc, size_t size);
void free_3d_array(void ***array);
void Mis_Cartesian2Spherical( const double Cartesian[], double Spherical[] );
void CartesianRotate( double x[], const double theta, const double phi, const bool inverse );
void Interpolation_UM_IC( real x, real y, real z, real ****Pri_input, real **XYZ, real *Pri_output, bool disk );
real TrilinearInterpolation(real *FieldAtVertices, real *xyz000, real *dxyz, real *xyz);
static void SetArrayDisk();
static void SetArrayHVC();
static real *randXYZ;
static int  N_CLOUDS = 0;
static real R_CLOUD  = 1.0;
void randCloud( real **randXYZ, const int numClouds );


#if ( MODEL == HYDRO )

// problem-specific global variables
// =======================================================================================

// options
       int      Jet_Ambient;             // [0/1/9]: uniform/Milky-Way/load-from-file
static bool     Jet_Precession;          // flag: precessing jet source
static bool     Jet_TimeDependentSrc;    // flag: time-dependent fluid variables in source
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
static double   Jet_SrcDens;             // jet density
static double   Jet_SrcTemp;             // jet temperature
static bool     Jet_SmoothVel;           // smooth radial component of 4-velocity on cross section

static bool     Jet_SphericalSrc;

#ifdef COSMIC_RAY
static double Jet_Src_CR_Engy;
static double Amb_CR_Engy;
#endif

// sound speed
static double   CharacteristicSpeed;     // the characteristic speed of the simulation problem
                                         // the default end-time (END_T) will be estimated from
                                         // `CharacteristicSpeed` and `BOX_SIZE`
static real *BUFFER;
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
static double jetSrc_highResRadius;
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
//
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
                                         // --> NOT necessary to be a unit vector

// time-depent source
static double   Jet_BurstStartTime;      // start burst time in jet source
static double   Jet_BurstEndTime;        // end burst time in jet source
static double   Jet_Burst4VelRatio;      // increase 4-velocity     by a factor of `Jet_Burst4VelRatio` during `Jet_BurstStartTime` and `Jet_BurstEndTime`
static double   Jet_BurstDensRatio;      // increase proper density by a factor of `Jet_BurstDensRatio` during `Jet_BurstStartTime` and `Jet_BurstEndTime`
static double   Jet_BurstTempRatio;      // increase temperature    by a factor of `Jet_BurstTempRatio` during `Jet_BurstStartTime` and `Jet_BurstEndTime`
static bool     Flag_Burst4Vel;          // flag: burst 4-velocity
static bool     Flag_BurstDens;          // flag: burst proper density
static bool     Flag_BurstTemp;          // flag: burst temperature

static double   Amb_FluSphereRadius;     //

#if (NCOMP_PASSIVE_USER > 0)
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
#  ifdef PARTICLE
   Aux_Error( ERROR_INFO, "PARTICLE must be disabled !!\n" );
#  endif

#  ifndef GRAVITY
// TODO: 3 and 4 should also be included
   if ( Jet_Ambient == 1 || Jet_Ambient == 2 )
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
   const char FileName[] = "Input__TestProb";
   ReadPara_t *ReadPara  = new ReadPara_t;

// (1-1) add parameters in the following format:
// --> note that VARIABLE, DEFAULT, MIN, and MAX must have the same data type
// --> some handy constants (e.g., Useless_bool, Eps_double, NoMin_int, ...) are defined in "include/ReadPara.h"
// ************************************************************************************************************************
// ReadPara->Add( "KEY_IN_THE_FILE",         &VARIABLE,             DEFAULT,                   MIN,            MAX    );
// ************************************************************************************************************************

// load options
   ReadPara->Add( "Jet_Ambient",             &Jet_Ambient,              1,                       0,              9    );
   ReadPara->Add( "Jet_Fire",                &Jet_Fire,                 3,                       0,              3    );
   ReadPara->Add( "Jet_Precession",          &Jet_Precession,           false,        Useless_bool,   Useless_bool    );
   ReadPara->Add( "Jet_SphericalSrc",        &Jet_SphericalSrc,         false,        Useless_bool,   Useless_bool    );
   ReadPara->Add( "Jet_TimeDependentSrc",    &Jet_TimeDependentSrc,     false,        Useless_bool,   Useless_bool    );
   ReadPara->Add( "Jet_Duration",            &Jet_Duration        ,     NoMax_double,          0.0,   NoMax_double    );

// load jet fluid parameters
   ReadPara->Add( "Jet_SrcVel",              &Jet_SrcVel    ,          -1.0,          NoMin_double,   NoMax_double    );
   ReadPara->Add( "Jet_SmoothVel",           &Jet_SmoothVel ,           false,        Useless_bool,   Useless_bool    );
   ReadPara->Add( "Jet_SrcDens",             &Jet_SrcDens   ,          -1.0,          Eps_double,     NoMax_double    );
   ReadPara->Add( "Jet_SrcTemp",             &Jet_SrcTemp   ,          -1.0,          Eps_double,     NoMax_double    );
   ReadPara->Add( "gasDisk_highResRadius",   &gasDisk_highResRadius,   -1.0,          NoMin_double,   NoMax_double    );
   ReadPara->Add( "gasDisk_lowRes_LEVEL",    &gasDisk_lowRes_LEVEL,    -1,                       0,      NoMax_int    );
   ReadPara->Add( "jetSrc_highResRadius",    &jetSrc_highResRadius,    -1.0,          NoMin_double,   NoMax_double    );
   ReadPara->Add( "jetSrc_lowRes_LEVEL",     &jetSrc_lowRes_LEVEL,     -1,                       0,      NoMax_int    );
#  ifdef COSMIC_RAY
   ReadPara->Add( "Jet_Src_CR_Engy",         &Jet_Src_CR_Engy   ,      -1.0,                   0.0,   NoMax_double    );
   ReadPara->Add( "Amb_CR_Engy",             &Amb_CR_Engy   ,          -1.0,                   0.0,   NoMax_double    );
#  endif


// load source geometry parameters
   ReadPara->Add( "Jet_Radius",              &Jet_Radius,              -1.0,          Eps_double,     NoMax_double    );
   ReadPara->Add( "Jet_HalfHeight",          &Jet_HalfHeight,          -1.0,          Eps_double,     NoMax_double    );
   ReadPara->Add( "Jet_HalfOpeningAngle",    &Jet_HalfOpeningAngle,    -1.0,                   0.0,           90.0    );
   ReadPara->Add( "Jet_CenOffset_x",         &Jet_CenOffset [0],        NoDef_double, NoMin_double,   NoMax_double    );
   ReadPara->Add( "Jet_CenOffset_y",         &Jet_CenOffset [1],        NoDef_double, NoMin_double,   NoMax_double    );
   ReadPara->Add( "Jet_CenOffset_z",         &Jet_CenOffset [2],        NoDef_double, NoMin_double,   NoMax_double    );

// load precission parameters
   ReadPara->Add( "Jet_AngularVelocity",     &Jet_AngularVelocity,      NoDef_double,          0.0,   NoMax_double    );
   ReadPara->Add( "Jet_PrecessionAngle",     &Jet_PrecessionAngle,      NoDef_double, NoMin_double,           90.0    );
   ReadPara->Add( "Jet_PrecessionAxis_x",    &Jet_PrecessionAxis[0],    NoDef_double, NoMin_double,   NoMax_double    );
   ReadPara->Add( "Jet_PrecessionAxis_y",    &Jet_PrecessionAxis[1],    NoDef_double, NoMin_double,   NoMax_double    );
   ReadPara->Add( "Jet_PrecessionAxis_z",    &Jet_PrecessionAxis[2],    NoDef_double, NoMin_double,   NoMax_double    );

// load uniform background parameters
   ReadPara->Add( "Amb_UniformDens",         &Amb_UniformDens,         -1.0,          Eps_double,     NoMax_double    );
   ReadPara->Add( "Amb_UniformVel_x",        &Amb_UniformVel[0],        0.0,          NoMin_double,   NoMax_double    );
   ReadPara->Add( "Amb_UniformVel_y",        &Amb_UniformVel[1],        0.0,          NoMin_double,   NoMax_double    );
   ReadPara->Add( "Amb_UniformVel_z",        &Amb_UniformVel[2],        0.0,          NoMin_double,   NoMax_double    );
   ReadPara->Add( "Amb_UniformTemp",         &Amb_UniformTemp,         -1.0,          Eps_double,     NoMax_double    );


   ReadPara->Add( "Amb_FluSphereRadius",     &Amb_FluSphereRadius,     -1.0,          NoMin_double,   NoMax_double    );
   ReadPara->Add( "CharacteristicSpeed",     &CharacteristicSpeed,     -1.0,          NoMin_double,   NoMax_double    );
   ReadPara->Add( "criticalTemp",            &criticalTemp,            -1.0,          NoMin_double,   NoMax_double    );

// load Milky Way parameters
   ReadPara->Add( "IsothermalSlab_Center_x", &IsothermalSlab_Center[0],          -1.0,          NoMin_double,   NoMax_double    );
   ReadPara->Add( "IsothermalSlab_Center_y", &IsothermalSlab_Center[1],          -1.0,          NoMin_double,   NoMax_double    );
   ReadPara->Add( "IsothermalSlab_Center_z", &IsothermalSlab_Center[2],          -1.0,          NoMin_double,   NoMax_double    );

// load time-dependent source varibles
   ReadPara->Add( "Jet_BurstStartTime",      &Jet_BurstStartTime,      -1.0,          NoMin_double,   NoMax_double    );
   ReadPara->Add( "Jet_BurstEndTime",        &Jet_BurstEndTime,        -1.0,          NoMin_double,   NoMax_double    );

   ReadPara->Read( FileName );

   delete ReadPara;

// Read header for the fermi bubbles
   SetArrayDisk();
// SetArrayHVC();
   randCloud( &randXYZ, N_CLOUDS );

// replace useless parameters with NaN
   if ( Jet_Ambient != 0 )
   {
     Amb_UniformDens      = NAN;
     Amb_UniformVel[0]    = NAN;
     Amb_UniformVel[1]    = NAN;
     Amb_UniformVel[2]    = NAN;
   }

   if ( Amb_FluSphereRadius < 0.0 )
   {
      Amb_FluSphereRadius = NAN;
   }

   if ( !Jet_TimeDependentSrc )
   {
     Jet_BurstDensRatio  = NAN;
     Jet_Burst4VelRatio  = NAN;
     Jet_BurstTempRatio  = NAN;
     Jet_BurstStartTime  = NAN;
     Jet_BurstEndTime    = NAN;
   }

// (1-2) check runtime parameters

// check time-dependent source
   if ( Jet_TimeDependentSrc )
   {
     if ( !Flag_Burst4Vel && !Flag_BurstDens && !Flag_BurstTemp )
       Aux_Error( ERROR_INFO, "One of Flag_Burst4Vel, Flag_BurstDens or Flag_BurstTemp must be enabled !!\n" );

     if ( Jet_BurstEndTime <= Jet_BurstStartTime )
       Aux_Error( ERROR_INFO, "Jet_BurstEndTime <= Jet_BurstStartTime !!\n" );

     if ( Jet_BurstEndTime >= END_T )
       Aux_Error( ERROR_INFO, "Jet_BurstEndTime >= END_T !!\n" );

     if ( Flag_Burst4Vel && Jet_Burst4VelRatio <= Eps_double )
       Aux_Error( ERROR_INFO, "Jet_Burst4VelRatio <= Eps_double !!\n" );

     if ( Flag_BurstDens && Jet_BurstDensRatio <= Eps_double )
       Aux_Error( ERROR_INFO, "Jet_BurstDensRatio <= Eps_double !!\n" );

     if ( Flag_BurstTemp && Jet_BurstTempRatio <= Eps_double )
       Aux_Error( ERROR_INFO, "Jet_BurstTempRatio <= Eps_double !!\n" );
   }

   if ( IsothermalSlab_Center[0] == -1.0 )
        IsothermalSlab_Center[0] = 0.5*amr->BoxSize[0];

   if ( IsothermalSlab_Center[1] == -1.0 )
        IsothermalSlab_Center[1] = 0.5*amr->BoxSize[1];

   if ( IsothermalSlab_Center[2] == -1.0 )
        IsothermalSlab_Center[2] = 0.5*amr->BoxSize[2];

   if ( Jet_Ambient == 9 && OPT__INIT != 3 )
   {
      Aux_Error( ERROR_INFO, "OPT__INIT must be 3 !!\n" );
   }


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
   jetSrc_highResRadius  *= Const_kpc / UNIT_L;

   if ( Jet_Ambient == 0 )
   {
      Amb_UniformDens   *= 1.0     / UNIT_D;
      Amb_UniformVel[0] *= Const_c / UNIT_V;
      Amb_UniformVel[1] *= Const_c / UNIT_V;
      Amb_UniformVel[2] *= Const_c / UNIT_V;
   } else if ( Jet_Ambient == 2  ||  Jet_Ambient == 3  ||  Jet_Ambient == 4 )
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
     ambientTemperature                *= Const_kB/(ParticleMass*UNIT_V*UNIT_V);

     gasDiskTemperature                 = Header_disk[18];
     gasDiskTemperature                *= Const_kB / (ParticleMass*UNIT_V*UNIT_V);

     criticalTemp                      *= Const_kB / (ParticleMass*Const_c*Const_c);
   } // if ( Jet_Ambient == 0 ) ... else if ...

   Amb_UniformTemp     *= Const_kB / (ParticleMass*Const_c*Const_c);
   Jet_AngularVelocity *= 1.0;    // the unit of Jet_AngularVelocity is UNIT_T


   if ( Amb_FluSphereRadius > 0.0 )
   {
      Amb_FluSphereRadius *= Const_kpc / UNIT_L;
   }


   if ( Jet_TimeDependentSrc )
   {
     Jet_BurstStartTime *= 1e3*Const_yr / UNIT_T;
     Jet_BurstEndTime   *= 1e3*Const_yr / UNIT_T;
     Jet_Burst4VelRatio *=     Const_c  / UNIT_V;
     Jet_BurstDensRatio *= 1.0          / UNIT_D;
   }

// (2) set the problem-specific derived parameters
   const double SecAngle = 1.0 / cos(0.5*Jet_HalfOpeningAngle);
   const double TanAngle = sin(0.5*Jet_HalfOpeningAngle) * SecAngle;

   Jet_MaxDis  = sqrt( SQR( Jet_Radius ) + SQR( Jet_HalfHeight * SecAngle ) + 2.0 * Jet_Radius * Jet_HalfHeight * TanAngle );

   for (int d=0; d<3; d++)    Jet_Center[d] = 0.5*amr->BoxSize[d] + Jet_CenOffset[d];


// (4) reset other general-purpose parameters
//     --> a helper macro PRINT_RESET_PARA is defined in Macro.h
   const long   End_Step_Default = __INT_MAX__;
   const double End_T_Default    = 0.5*BOX_SIZE * UNIT_L / (CharacteristicSpeed *UNIT_V) / UNIT_T;

   if ( END_STEP < 0 ) {
      END_STEP = End_Step_Default;
      PRINT_RESET_PARA( END_STEP, FORMAT_LONG, "" );
   }

   if ( END_T < 0.0 ) {
      if ( CharacteristicSpeed == -1.0 ) Aux_Error( ERROR_INFO, "CharacteristicSpeed must be provided !!\n" );
      else                               END_T = End_T_Default;

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
      Aux_Message( stdout, "  Jet_Precession           = %d\n",                Jet_Precession                                  );
      Aux_Message( stdout, "  Jet_SphericalSrc         = %d\n",                Jet_SphericalSrc                                );
      Aux_Message( stdout, "  Jet_TimeDependentSrc     = %d\n",                Jet_TimeDependentSrc                            );
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
      Aux_Message( stdout, "  Jet_CenOffset[x]         = %14.7e kpc\n",        Jet_CenOffset [0]*UNIT_L/Const_kpc              );
      Aux_Message( stdout, "  Jet_CenOffset[y]         = %14.7e kpc\n",        Jet_CenOffset [1]*UNIT_L/Const_kpc              );
      Aux_Message( stdout, "  Jet_CenOffset[z]         = %14.7e kpc\n",        Jet_CenOffset [2]*UNIT_L/Const_kpc              );
      Aux_Message( stdout, "  Jet_AngularVelocity      = %14.7e degree/kyr\n", Jet_AngularVelocity                             );
      Aux_Message( stdout, "  Jet_PrecessionAngle      = %14.7e degree\n",     Jet_PrecessionAngle*180.0/M_PI                  );
      Aux_Message( stdout, "  Jet_HalfOpeningAngle     = %14.7e degree\n",     Jet_HalfOpeningAngle*180.0/M_PI                 );
      Aux_Message( stdout, "  Jet_Radius               = %14.7e kpc\n",        Jet_Radius*UNIT_L/Const_kpc                     );
      Aux_Message( stdout, "  Jet_HalfHeight           = %14.7e kpc\n",        Jet_HalfHeight*UNIT_L/Const_kpc                 );
      Aux_Message( stdout, "  Jet_MaxDis               = %14.7e kpc\n",        Jet_MaxDis*UNIT_L/Const_kpc                     );
      Aux_Message( stdout, "  gasDisk_highResRadius    = %14.7e kpc\n",        gasDisk_highResRadius*UNIT_L/Const_kpc          );
      Aux_Message( stdout, "  gasDisk_lowRes_LEVEL     = %d\n",                gasDisk_lowRes_LEVEL                            );
      Aux_Message( stdout, "  jetSrc_highResRadius     = %14.7e kpc\n",        jetSrc_highResRadius*UNIT_L/Const_kpc           );
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
      Aux_Message( stdout, "  criticalTemp             = %14.7e K\n",          criticalTemp / ( Const_kB /( ParticleMass * Const_c * Const_c ) ) );
      } // if ( Jet_Ambient == 0 ) ... else if ...

      Aux_Message( stdout, "  CharacteristicSpeed      = %14.7e c\n",          CharacteristicSpeed / UNIT_V                    );
      Aux_Message( stdout, "  Jet_PrecessionAxis[x]    = %14.7e\n",            Jet_PrecessionAxis[0]                           );
      Aux_Message( stdout, "  Jet_PrecessionAxis[y]    = %14.7e\n",            Jet_PrecessionAxis[1]                           );
      Aux_Message( stdout, "  Jet_PrecessionAxis[z]    = %14.7e\n",            Jet_PrecessionAxis[2]                           );

      if ( Amb_FluSphereRadius > 0.0 )
      {
      Aux_Message( stdout, "  Amb_FluSphereRadius      = %14.7e kpc\n",        Amb_FluSphereRadius*UNIT_L/Const_kpc            );
      } // if ( Amb_FluSphereRadius > 0.0 )

      if ( Jet_TimeDependentSrc )
      {
      Aux_Message( stdout, "  Jet_BurstStartTime       = %14.7e kyr \n",       Jet_BurstStartTime*UNIT_T/(1e3*Const_yr)        );
      Aux_Message( stdout, "  Jet_BurstEndTime         = %14.7e kyr \n",       Jet_BurstEndTime*UNIT_T/(1e3*Const_yr)          );
      Aux_Message( stdout, "  Jet_Burst4VelRatio       = %14.7e c \n",         Jet_Burst4VelRatio                              );
      Aux_Message( stdout, "  Jet_BurstDensRatio       = %14.7e g/cm^3\n",     Jet_BurstDensRatio*UNIT_D                       );
      Aux_Message( stdout, "  Jet_BurstTempRatio       = %14.7e\n",            Jet_BurstTempRatio                              );
      Aux_Message( stdout, "  Flag_Burst4Vel           = %d\n",                Flag_Burst4Vel                                  );
      Aux_Message( stdout, "  Flag_BurstDens           = %d\n",                Flag_BurstDens                                  );
      Aux_Message( stdout, "  Flag_BurstTemp           = %d\n",                Flag_BurstTemp                                  );
      } // if ( Jet_TimeDependentSrc )

      Aux_Message( stdout, "=============================================================================\n"                   );
   } // if ( MPI_Rank == 0 )

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Setting runtime parameters ... done\n"                                     );

} // FUNCTION : SetParameter



void ReadBinFile( char *FileName, real **buffer )
{
  FILE *pFile;
  long lSize;
  size_t result;

  pFile = fopen( FileName, "rb" );
  if ( pFile == NULL ) Aux_Error( ERROR_INFO, "File error !!\n" );

  // obtain file size
  fseek( pFile, 0, SEEK_END );
  lSize = ftell( pFile );
  rewind( pFile );

  // allocate memory to contain the whole file
  *buffer = (real*) calloc( lSize, sizeof(double) );
  if ( *buffer == NULL ) Aux_Error( ERROR_INFO, "Memory error !!\n" );

  // copy the file into the *buffer
  result = fread( *buffer, 1, lSize, pFile );
  if ( result != lSize ) Aux_Error( ERROR_INFO, "Reading error !!\n" );

  fclose (pFile);

} // FUNCTION : ReadBinFile



void randCloud( real **randXYZ, const int numClouds )
{
   const real Lx = amr->BoxSize[0];
   const real Ly = amr->BoxSize[1];
   const real Lz = amr->BoxSize[2];

   // TODO: should this be a const from input? how to decided this parameter?
   const real z_max = Lz*0.5 + 3.0;

   *randXYZ = (real*)malloc( 3*numClouds*sizeof(real) );

   for (int i=0; i<3*numClouds; i+=3)
   {
      (*randXYZ)[i+0] = (real)rand() / (real)RAND_MAX * Lx;
      (*randXYZ)[i+1] = (real)rand() / (real)RAND_MAX * Ly;
      (*randXYZ)[i+2] = (real)rand() / (real)RAND_MAX * Lz;

      if ( (*randXYZ)[i+2] < z_max )   i-=3;
   } // for (int i=0; i<3*numClouds; i+=3)
} // FUNCTION : randCloud



bool checkInsideClouds( const real *randXYZ, const int numClouds, const real x, const real y, const real z, real *cloudCenter )
{
   real distance;

   for (int i=0; i<3*numClouds; i+=3)
   {
      distance = SQRT( SQR( x-randXYZ[i+0] ) + SQR( y-randXYZ[i+1] ) + SQR( z-randXYZ[i+2] ) );

      if ( distance >= R_CLOUD )   continue;

      for (int d=0; d<3; d++)   cloudCenter[d] = randXYZ[i+d];
      return true;
   }
   return false;
} // FUNCTION : checkInsideClouds



void SetArrayDisk()
{
// Reading table for interpolations in SetGridIC()
   char TableFileName[] = "UM_IC";
   ReadBinFile( TableFileName, &BUFFER );

   const int headerSize = (int)BUFFER[0];

   Header_disk = (real*)malloc( (size_t)headerSize * sizeof(real) );

   memcpy( Header_disk, BUFFER, (size_t)headerSize * sizeof(real) );

   ParticleMass = Header_disk[8] * Header_disk[9];

   real Lx = Header_disk[12];
   real Ly = Header_disk[13];
   real Lz = Header_disk[14];

   interfaceHeight  = Header_disk[17];
   interfaceHeight *= Const_kpc / UNIT_L;

   int Nx = (int)Header_disk[22];
   int Ny = (int)Header_disk[23];
   int Nz = (int)Header_disk[24];

   const int numGhost = (int)Header_disk[25];

   real dx = Lx / (real)Nx;
   real dy = Lx / (real)Nx;
   real dz = Lx / (real)Nx;

   const int NX    = Nx + 2*numGhost;
   const int NY    = Ny + 2*numGhost;
   const int NZ    = Nz + 2*numGhost;
   const int N_TOT = NX * NY * NZ;

   if ( 5*N_TOT > INT_MAX )   Aux_Error( ERROR_INFO, "Integer overflow (5*N_TOT > INT_MAX), N_TOT=NX*NY*NZ, NX=%d, NY=%d, NZ=%d, INT_MAX=%d!! \n", NX, NY, NZ, INT_MAX );

   if ( Step == 0 )
   {
      Rhoo_disk = (real***)calloc_3d_array( (size_t)NX, (size_t)NY, (size_t)NZ, sizeof(real) );
      VelX_disk = (real***)calloc_3d_array( (size_t)NX, (size_t)NY, (size_t)NZ, sizeof(real) );
      VelY_disk = (real***)calloc_3d_array( (size_t)NX, (size_t)NY, (size_t)NZ, sizeof(real) );
      VelZ_disk = (real***)calloc_3d_array( (size_t)NX, (size_t)NY, (size_t)NZ, sizeof(real) );
      Pres_disk = (real***)calloc_3d_array( (size_t)NX, (size_t)NY, (size_t)NZ, sizeof(real) );

      X_disk = (real*)calloc( (size_t)NX, sizeof(real) );
      Y_disk = (real*)calloc( (size_t)NY, sizeof(real) );
      Z_disk = (real*)calloc( (size_t)NZ, sizeof(real) );

      if ( X_disk == NULL ) Aux_Error( ERROR_INFO, "X_disk is NULL at %d, %s !!\n", __LINE__, __FUNCTION__ );
      if ( Y_disk == NULL ) Aux_Error( ERROR_INFO, "Y_disk is NULL at %d, %s !!\n", __LINE__, __FUNCTION__ );
      if ( Z_disk == NULL ) Aux_Error( ERROR_INFO, "Z_disk is NULL at %d, %s !!\n", __LINE__, __FUNCTION__ );

      real *Ptr;
      Ptr = BUFFER + headerSize;

      for (int c=0; c<5*N_TOT; c++) {
         const int cc = c%(N_TOT);
         const int i  = (cc - cc%(NY*NZ)) / (NY*NZ);
         const int j  = ((cc - cc%NZ) / NZ) % NY;
         const int k  = cc%NZ;

         if ( 0       <= c && c <   N_TOT )   Rhoo_disk[i][j][k] = Ptr[c];
         if (   N_TOT <= c && c < 2*N_TOT )   VelX_disk[i][j][k] = Ptr[c];
         if ( 2*N_TOT <= c && c < 3*N_TOT )   VelY_disk[i][j][k] = Ptr[c];
         if ( 3*N_TOT <= c && c < 4*N_TOT )   VelZ_disk[i][j][k] = Ptr[c];
         if ( 4*N_TOT <= c && c < 5*N_TOT )   Pres_disk[i][j][k] = Ptr[c];
      }

      Ptr += 5*N_TOT;
      for (int c=0; c<NX; c++) X_disk[c] = Ptr[c];

      Ptr += NX;
      for (int c=0; c<NY; c++) Y_disk[c] = Ptr[c];

      Ptr += NY;
      for (int c=0; c<NZ; c++) Z_disk[c] = Ptr[c];
   } // if ( Step == 0 )

   free(BUFFER);

} // FUNCTION : SetArrayDisk



void SetArrayHVC()
{
// Reading table for interpolations in SetGridIC()
   char TableFileName[] = "UM_IC_HVC";
   ReadBinFile( TableFileName, &BUFFER );

   const int headerSize = (int)BUFFER[0];

   Header_hvc = (real*)malloc( (size_t)headerSize*sizeof(real) );

   memcpy( Header_hvc, BUFFER, (size_t)headerSize*sizeof(real) );

   ParticleMass = Header_hvc[8]*Header_hvc[9];

   real Lx = Header_hvc[12];
   real Ly = Header_hvc[13];
   real Lz = Header_hvc[14];

   interfaceHeight  = Header_hvc[17];
   interfaceHeight *= Const_kpc/UNIT_L;

   int Nx = (int)Header_hvc[22];
   int Ny = (int)Header_hvc[23];
   int Nz = (int)Header_hvc[24];

   const int numGhost = (int)Header_hvc[25];

   real dx = Lx / (real)Nx;
   real dy = Lx / (real)Nx;
   real dz = Lx / (real)Nx;

   const int NX = Nx + 2*numGhost;
   const int NY = Ny + 2*numGhost;
   const int NZ = Nz + 2*numGhost;
   const int N_TOT = NX * NY * NZ;

   if ( 5*N_TOT > INT_MAX )   Aux_Error( ERROR_INFO, "Integer overflow (5*N_TOT > INT_MAX), N_TOT=NX*NY*NZ, NX=%d, NY=%d, NZ=%d, INT_MAX=%d!! \n", NX, NY, NZ, INT_MAX );

   if ( Step == 0 )
   {
      Rhoo_hvc = (real***)calloc_3d_array( (size_t)NX, (size_t)NY, (size_t)NZ, sizeof(real) );
      VelX_hvc = (real***)calloc_3d_array( (size_t)NX, (size_t)NY, (size_t)NZ, sizeof(real) );
      VelY_hvc = (real***)calloc_3d_array( (size_t)NX, (size_t)NY, (size_t)NZ, sizeof(real) );
      VelZ_hvc = (real***)calloc_3d_array( (size_t)NX, (size_t)NY, (size_t)NZ, sizeof(real) );
      Pres_hvc = (real***)calloc_3d_array( (size_t)NX, (size_t)NY, (size_t)NZ, sizeof(real) );

      X_hvc = (real*)calloc( (size_t)NX,sizeof(real) );
      Y_hvc = (real*)calloc( (size_t)NY,sizeof(real) );
      Z_hvc = (real*)calloc( (size_t)NZ,sizeof(real) );

      if ( X_hvc == NULL ) Aux_Error( ERROR_INFO, "X_hvc is NULL at %d, %s !!\n", __LINE__, __FUNCTION__ );
      if ( Y_hvc == NULL ) Aux_Error( ERROR_INFO, "Y_hvc is NULL at %d, %s !!\n", __LINE__, __FUNCTION__ );
      if ( Z_hvc == NULL ) Aux_Error( ERROR_INFO, "Z_hvc is NULL at %d, %s !!\n", __LINE__, __FUNCTION__ );

      real *Ptr;
      Ptr = BUFFER + headerSize;

      for ( int c=0; c<5*N_TOT; c++ ) {
        const int cc = c%(N_TOT);
        const int i = (cc - cc%(NY*NZ)) / (NY*NZ);
        const int j = ((cc - cc%NZ) / NZ) % NY;
        const int k = cc%NZ;

        if ( 0       <= c && c <   N_TOT )   Rhoo_hvc[i][j][k] = Ptr[c];
        if (   N_TOT <= c && c < 2*N_TOT )   VelX_hvc[i][j][k] = Ptr[c];
        if ( 2*N_TOT <= c && c < 3*N_TOT )   VelY_hvc[i][j][k] = Ptr[c];
        if ( 3*N_TOT <= c && c < 4*N_TOT )   VelZ_hvc[i][j][k] = Ptr[c];
        if ( 4*N_TOT <= c && c < 5*N_TOT )   Pres_hvc[i][j][k] = Ptr[c];
      }

      Ptr += 5*N_TOT;
      for (int c=0; c<NX; c++) X_hvc[c] = Ptr[c];

      Ptr += NX;
      for (int c=0; c<NY; c++) Y_hvc[c] = Ptr[c];

      Ptr += NY;
      for (int c=0; c<NZ; c++) Z_hvc[c] = Ptr[c];
   } // if ( Step == 0 )

   free(BUFFER);
} // FUNCTION : SetArrayHVC



void Interpolation_UM_IC( real x, real y, real z, real ****Pri_input, real **XYZ, real *Pri_output, bool disk )
{
   real *Header = (disk) ? Header_disk : Header_hvc;

   real xyz[3] = {x, y, z};

   const int numGhost = (int)Header[25];

// N{X,Y,Z} = N{x,y,z} + 2*N_ghost
   const int NX = (int)Header[22] + 2*numGhost;
   const int NY = (int)Header[23] + 2*numGhost;
   const int NZ = (int)Header[24] + 2*numGhost;

   real dxyz[3] = { Header[26], Header[27], Header[28] };  // {dx, dy, dz}

   const int Idx = Mis_BinarySearch_Real( XYZ[0], 0, NX-1, x );
   const int Jdx = Mis_BinarySearch_Real( XYZ[1], 0, NY-1, y );
   const int Kdx = Mis_BinarySearch_Real( XYZ[2], 0, NZ-1, z );

   if ( Idx < 0  ||  Idx > NX-2 ) Aux_Error( ERROR_INFO, "x=%e is out of range! XYZ[0][0]=%e, XYZ[0][%d]=%e\n", x, XYZ[0][0], NX-1, XYZ[0][NX-1] );
   if ( Jdx < 0  ||  Jdx > NY-2 ) Aux_Error( ERROR_INFO, "y=%e is out of range! XYZ[1][1]=%e, XYZ[1][%d]=%e\n", y, XYZ[1][0], NY-1, XYZ[1][NY-1] );
   if ( Kdx < 0  ||  Kdx > NZ-2 ) Aux_Error( ERROR_INFO, "z=%e is out of range! XYZ[2][2]=%e, XYZ[2][%d]=%e\n", z, XYZ[2][0], NZ-1, XYZ[2][NZ-1] );

   real Vertices[8][5];
   for (int v=0; v<5; v++)
   {
      Vertices[0][v] = Pri_input[v][Idx  ][Jdx  ][Kdx  ];
      Vertices[1][v] = Pri_input[v][Idx  ][Jdx  ][Kdx+1];
      Vertices[2][v] = Pri_input[v][Idx  ][Jdx+1][Kdx  ];
      Vertices[3][v] = Pri_input[v][Idx+1][Jdx  ][Kdx  ];
      Vertices[4][v] = Pri_input[v][Idx  ][Jdx+1][Kdx+1];
      Vertices[5][v] = Pri_input[v][Idx+1][Jdx  ][Kdx+1];
      Vertices[6][v] = Pri_input[v][Idx+1][Jdx+1][Kdx  ];
      Vertices[7][v] = Pri_input[v][Idx+1][Jdx+1][Kdx+1];
   }

// TODO: One the above is confirmed to be correct, remove the original code
   // real Vertex000[5] = { Pri_input[0][Idx  ][Jdx  ][Kdx  ], Pri_input[1][Idx  ][Jdx  ][Kdx  ], Pri_input[2][Idx  ][Jdx  ][Kdx  ], Pri_input[3][Idx  ][Jdx  ][Kdx  ], Pri_input[4][Idx  ][Jdx  ][Kdx  ] };
   // real Vertex001[5] = { Pri_input[0][Idx  ][Jdx  ][Kdx+1], Pri_input[1][Idx  ][Jdx  ][Kdx+1], Pri_input[2][Idx  ][Jdx  ][Kdx+1], Pri_input[3][Idx  ][Jdx  ][Kdx+1], Pri_input[4][Idx  ][Jdx  ][Kdx+1] };
   // real Vertex010[5] = { Pri_input[0][Idx  ][Jdx+1][Kdx  ], Pri_input[1][Idx  ][Jdx+1][Kdx  ], Pri_input[2][Idx  ][Jdx+1][Kdx  ], Pri_input[3][Idx  ][Jdx+1][Kdx  ], Pri_input[4][Idx  ][Jdx+1][Kdx  ] };
   // real Vertex100[5] = { Pri_input[0][Idx+1][Jdx  ][Kdx  ], Pri_input[1][Idx+1][Jdx  ][Kdx  ], Pri_input[2][Idx+1][Jdx  ][Kdx  ], Pri_input[3][Idx+1][Jdx  ][Kdx  ], Pri_input[4][Idx+1][Jdx  ][Kdx  ] };
   // real Vertex011[5] = { Pri_input[0][Idx  ][Jdx+1][Kdx+1], Pri_input[1][Idx  ][Jdx+1][Kdx+1], Pri_input[2][Idx  ][Jdx+1][Kdx+1], Pri_input[3][Idx  ][Jdx+1][Kdx+1], Pri_input[4][Idx  ][Jdx+1][Kdx+1] };
   // real Vertex101[5] = { Pri_input[0][Idx+1][Jdx  ][Kdx+1], Pri_input[1][Idx+1][Jdx  ][Kdx+1], Pri_input[2][Idx+1][Jdx  ][Kdx+1], Pri_input[3][Idx+1][Jdx  ][Kdx+1], Pri_input[4][Idx+1][Jdx  ][Kdx+1] };
   // real Vertex110[5] = { Pri_input[0][Idx+1][Jdx+1][Kdx  ], Pri_input[1][Idx+1][Jdx+1][Kdx  ], Pri_input[2][Idx+1][Jdx+1][Kdx  ], Pri_input[3][Idx+1][Jdx+1][Kdx  ], Pri_input[4][Idx+1][Jdx+1][Kdx  ] };
   // real Vertex111[5] = { Pri_input[0][Idx+1][Jdx+1][Kdx+1], Pri_input[1][Idx+1][Jdx+1][Kdx+1], Pri_input[2][Idx+1][Jdx+1][Kdx+1], Pri_input[3][Idx+1][Jdx+1][Kdx+1], Pri_input[4][Idx+1][Jdx+1][Kdx+1] };

   bool Unphy = false;
   for (int i=0; i<8; i++)
   {
      for (int v=0; v<5; v++)   Unphy |= ( Vertices[i][v] != Vertices[i][v] );
      Unphy |= ( (real) TINY_NUMBER >= Vertices[i][0] );
      Unphy |= ( (real)-HUGE_NUMBER >= Vertices[i][1] );
      Unphy |= ( (real)-HUGE_NUMBER >= Vertices[i][2] );
      Unphy |= ( (real)-HUGE_NUMBER >= Vertices[i][3] );
      Unphy |= ( (real) TINY_NUMBER >= Vertices[i][4] );
      for (int v=0; v<5; v++)   Unphy |= ( Vertices[i][v] >= (real)HUGE_NUMBER );

      // TODO : Can not use since the function also check for the passive scalar
      // Unphy |= Hydro_IsUnphysical( UNPHY_MODE_PRIM, Vertices[i], NULL, NULL, NULL, NULL, EoS_DensEint2Pres_CPUPtr, EoS_GuessHTilde_CPUPtr, EoS_HTilde2Temp_CPUPtr, EoS_AuxArray_Flt, EoS_AuxArray_Int, h_EoS_Table, __FILE__,  __LINE__, __FUNCTION__, UNPHY_VERBOSE );
   } // for (int i=0; i<8; i++)

   // Unphy |= Hydro_IsUnphysical( UNPHY_MODE_PRIM, Vertex000, NULL, NULL, NULL, NULL, EoS_DensEint2Pres_CPUPtr, EoS_GuessHTilde_CPUPtr, EoS_HTilde2Temp_CPUPtr, EoS_AuxArray_Flt, EoS_AuxArray_Int, h_EoS_Table, __FILE__,  __LINE__, __FUNCTION__, true );
   // Unphy |= Hydro_IsUnphysical( UNPHY_MODE_PRIM, Vertex001, NULL, NULL, NULL, NULL, EoS_DensEint2Pres_CPUPtr, EoS_GuessHTilde_CPUPtr, EoS_HTilde2Temp_CPUPtr, EoS_AuxArray_Flt, EoS_AuxArray_Int, h_EoS_Table, __FILE__,  __LINE__, __FUNCTION__, true );
   // Unphy |= Hydro_IsUnphysical( UNPHY_MODE_PRIM, Vertex010, NULL, NULL, NULL, NULL, EoS_DensEint2Pres_CPUPtr, EoS_GuessHTilde_CPUPtr, EoS_HTilde2Temp_CPUPtr, EoS_AuxArray_Flt, EoS_AuxArray_Int, h_EoS_Table, __FILE__,  __LINE__, __FUNCTION__, true );
   // Unphy |= Hydro_IsUnphysical( UNPHY_MODE_PRIM, Vertex100, NULL, NULL, NULL, NULL, EoS_DensEint2Pres_CPUPtr, EoS_GuessHTilde_CPUPtr, EoS_HTilde2Temp_CPUPtr, EoS_AuxArray_Flt, EoS_AuxArray_Int, h_EoS_Table, __FILE__,  __LINE__, __FUNCTION__, true );
   // Unphy |= Hydro_IsUnphysical( UNPHY_MODE_PRIM, Vertex011, NULL, NULL, NULL, NULL, EoS_DensEint2Pres_CPUPtr, EoS_GuessHTilde_CPUPtr, EoS_HTilde2Temp_CPUPtr, EoS_AuxArray_Flt, EoS_AuxArray_Int, h_EoS_Table, __FILE__,  __LINE__, __FUNCTION__, true );
   // Unphy |= Hydro_IsUnphysical( UNPHY_MODE_PRIM, Vertex110, NULL, NULL, NULL, NULL, EoS_DensEint2Pres_CPUPtr, EoS_GuessHTilde_CPUPtr, EoS_HTilde2Temp_CPUPtr, EoS_AuxArray_Flt, EoS_AuxArray_Int, h_EoS_Table, __FILE__,  __LINE__, __FUNCTION__, true );
   // Unphy |= Hydro_IsUnphysical( UNPHY_MODE_PRIM, Vertex101, NULL, NULL, NULL, NULL, EoS_DensEint2Pres_CPUPtr, EoS_GuessHTilde_CPUPtr, EoS_HTilde2Temp_CPUPtr, EoS_AuxArray_Flt, EoS_AuxArray_Int, h_EoS_Table, __FILE__,  __LINE__, __FUNCTION__, true );
   // Unphy |= Hydro_IsUnphysical( UNPHY_MODE_PRIM, Vertex111, NULL, NULL, NULL, NULL, EoS_DensEint2Pres_CPUPtr, EoS_GuessHTilde_CPUPtr, EoS_HTilde2Temp_CPUPtr, EoS_AuxArray_Flt, EoS_AuxArray_Int, h_EoS_Table, __FILE__,  __LINE__, __FUNCTION__, true );

   // TODO: should we raise an error? Original exit(0)
   if ( Unphy ) Aux_Error( ERROR_INFO, "Unphysical near : Idx=%d, Jdx=%d, Kdx=%d, x=%e, y=%e, z=%e !!\n", Idx, Jdx, Kdx, x, y, z );

   real xyz000[3] = {XYZ[0][Idx], XYZ[1][Jdx], XYZ[2][Kdx] };
   real FieldAtVertices[8];

   for ( int v=0; v<5; v++ ) {
      for (int i=0; i<8; i++)   FieldAtVertices[i] = Vertices[i][v];
      // real FieldAtVertices[8] = { Vertex000[v], Vertex001[v], Vertex010[v], Vertex100[v], Vertex011[v], Vertex101[v], Vertex110[v], Vertex111[v] };
      Pri_output[v] = TrilinearInterpolation( FieldAtVertices, xyz000, dxyz, xyz );
   }

   //free_3d_array((void***)Rhoo_disk);
   //free_3d_array((void***)VelX_disk);
   //free_3d_array((void***)VelY_disk);
   //free_3d_array((void***)VelZ_disk);
   //free_3d_array((void***)Pres_disk);
   //free(X_disk);
   //free(Y_disk);
   //free(Z_disk);
   //free(BUFFER);
} // FUNCTION : Interpolation_UM_IC



#ifdef GRAVITY
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
   real Pri[NCOMP_FLUID] = {0.0};
   real xc = x - IsothermalSlab_Center[0];
   real yc = y - IsothermalSlab_Center[1];
   real zc = z - IsothermalSlab_Center[2];

   real ***Pri_disk_input[5] = { Rhoo_disk, VelX_disk, VelY_disk, VelZ_disk, Pres_disk };
   real *** Pri_hvc_input[5] = { Rhoo_hvc,   VelX_hvc,  VelY_hvc,  VelZ_hvc,  Pres_hvc };

// TODO: Jet_Ambient == 1 is missing
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


      if ( Hydro_IsUnphysical( UNPHY_MODE_PRIM, Pri, NULL, NULL, NULL, NULL, EoS_DensEint2Pres_CPUPtr,
                               EoS_GuessHTilde_CPUPtr, EoS_HTilde2Temp_CPUPtr, EoS_AuxArray_Flt,
                               EoS_AuxArray_Int, h_EoS_Table, __FILE__,  __LINE__, __FUNCTION__, UNPHY_VERBOSE ) )
         Aux_Error( ERROR_INFO, "Unphysical cell at (%e, %e %e)\n", x, y, z);

#     if (NCOMP_PASSIVE_USER > 0)
      fluid[Passive_0000] = fluid[DENS];
      fluid[Passive_0001] = 0.0;
      fluid[Passive_0002] = 0.0;
#     endif
#     ifdef COSMIC_RAY
      fluid[CRAY] = Amb_CR_Engy * SQRT((real)1.0+SQR(Jet_SrcVel));
#     endif
   }
   else if ( Jet_Ambient == 2 ) // cold disk in stratified ambient
   {
#     ifdef GRAVITY
      if ( fabs(zc) < interfaceHeight )
      {
         real *XYZ[3] = { X_disk, Y_disk, Z_disk };

         Interpolation_UM_IC( xc, yc, zc, Pri_disk_input, XYZ, Pri, true );

         if ( Hydro_IsUnphysical( UNPHY_MODE_PRIM, Pri, NULL, NULL, NULL, NULL, EoS_DensEint2Pres_CPUPtr,
                                  EoS_GuessHTilde_CPUPtr, EoS_HTilde2Temp_CPUPtr, EoS_AuxArray_Flt,
                                  EoS_AuxArray_Int, h_EoS_Table, __FILE__,  __LINE__, __FUNCTION__, UNPHY_VERBOSE ) )
            Aux_Error( ERROR_INFO, "Unphysical cell at (%e, %e %e)\n", x, y, z);

         if ( Pri[4]/Pri[0] > criticalTemp ) Pri[0] = Pri[4] / ambientTemperature;


         Hydro_Pri2Con( Pri, fluid, false, PassiveNorm_NVar, PassiveNorm_VarIdx, EoS_DensPres2Eint_CPUPtr,
                        EoS_Temp2HTilde_CPUPtr, EoS_HTilde2Temp_CPUPtr, EoS_AuxArray_Flt, EoS_AuxArray_Int,
                        h_EoS_Table, NULL );

#        if (NCOMP_PASSIVE_USER > 0)
         fluid[Passive_0000] = fluid[DENS];
         fluid[Passive_0001] = 0.0;
         fluid[Passive_0002] = 0.0;
#        endif

#        ifdef COSMIC_RAY
         fluid[CRAY] = Amb_CR_Engy * SQRT((real)1.0+SQR(Jet_SrcVel));
#        endif
      } else // if ( fabs(zc) < interfaceHeight )
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

         if ( checkInsideClouds( randXYZ, N_CLOUDS, x, y, z, cloudCenter ) )
         {
            real Pri_hvc_output[5];

            real *XYZ[3] =  { X_hvc, Y_hvc, Z_hvc };

            Interpolation_UM_IC( x-cloudCenter[0], y-cloudCenter[1], z-cloudCenter[2],
                                 Pri_hvc_input, XYZ, Pri_hvc_output, false );

            if ( Hydro_IsUnphysical( UNPHY_MODE_PRIM, Pri_hvc_output, NULL, NULL, NULL, NULL, EoS_DensEint2Pres_CPUPtr,
                                     EoS_GuessHTilde_CPUPtr, EoS_HTilde2Temp_CPUPtr, EoS_AuxArray_Flt,
                                     EoS_AuxArray_Int, h_EoS_Table, __FILE__,  __LINE__, __FUNCTION__, UNPHY_VERBOSE ) )
               Aux_Error( ERROR_INFO, "Unphysical cell at (%e, %e %e)\n", x, y, z);

            Pri[0] = Pri_hvc_output[0]*0.05*Const_mp/UNIT_D;
            Pri[1] = 0.0;
            Pri[2] = 0.0;
            Pri[3] = 0.0;
            Pri[4] = Pri[0];
         }

         Hydro_Pri2Con( Pri, fluid, false, PassiveNorm_NVar, PassiveNorm_VarIdx, EoS_DensPres2Eint_CPUPtr,
                        EoS_Temp2HTilde_CPUPtr, EoS_HTilde2Temp_CPUPtr, EoS_AuxArray_Flt, EoS_AuxArray_Int,
                        h_EoS_Table, NULL );

         if ( Hydro_IsUnphysical( UNPHY_MODE_PRIM, Pri, NULL, NULL, NULL, NULL, EoS_DensEint2Pres_CPUPtr,
                                  EoS_GuessHTilde_CPUPtr, EoS_HTilde2Temp_CPUPtr, EoS_AuxArray_Flt,
                                  EoS_AuxArray_Int, h_EoS_Table, __FILE__,  __LINE__, __FUNCTION__, UNPHY_VERBOSE ) )
            Aux_Error( ERROR_INFO, "Unphysical cell at (%e, %e %e)\n", x, y, z);

#        if (NCOMP_PASSIVE_USER > 0)
         fluid[Passive_0000] = 0.0;
         fluid[Passive_0001] = 0.0;
         fluid[Passive_0002] = 0.0;
#        endif

#        ifdef COSMIC_RAY
         fluid[CRAY] = Amb_CR_Engy * SQRT((real)1.0+SQR(Jet_SrcVel));
#        endif
      } // if ( fabs(zc) < interfaceHeight ) ... else ...
#  endif
   } else if ( Jet_Ambient == 3 ) // cold disk in uniform ambient
   {
#     ifdef GRAVITY
      if ( fabs(zc) < interfaceHeight )
      {
         real *XYZ[3] =  { X_disk, Y_disk, Z_disk };

         Interpolation_UM_IC( xc, yc, zc, Pri_disk_input, XYZ, Pri, true );

         if ( Hydro_IsUnphysical( UNPHY_MODE_PRIM, Pri, NULL, NULL, NULL, NULL, EoS_DensEint2Pres_CPUPtr,
                                  EoS_GuessHTilde_CPUPtr, EoS_HTilde2Temp_CPUPtr, EoS_AuxArray_Flt,
                                  EoS_AuxArray_Int, h_EoS_Table, __FILE__,  __LINE__, __FUNCTION__, UNPHY_VERBOSE ) )
            Aux_Error( ERROR_INFO, "Unphysical cell at (%e, %e %e)\n", x, y, z);

         Hydro_Pri2Con( Pri, fluid, false, PassiveNorm_NVar, PassiveNorm_VarIdx, EoS_DensPres2Eint_CPUPtr,
                        EoS_Temp2HTilde_CPUPtr, EoS_HTilde2Temp_CPUPtr, EoS_AuxArray_Flt, EoS_AuxArray_Int,
                        h_EoS_Table, NULL );

#        if (NCOMP_PASSIVE_USER > 0)
         fluid[Passive_0000] = 0.0;
         fluid[Passive_0001] = 0.0;
         fluid[Passive_0002] = 0.0;
#        endif

#        ifdef COSMIC_RAY
         fluid[CRAY] = Amb_CR_Engy * SQRT((real)1.0+SQR(Jet_SrcVel));
#        endif
      } else // if ( fabs(zc) < interfaceHeight )
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

         if ( Hydro_IsUnphysical( UNPHY_MODE_PRIM, Pri, NULL, NULL, NULL, NULL, EoS_DensEint2Pres_CPUPtr,
                                  EoS_GuessHTilde_CPUPtr, EoS_HTilde2Temp_CPUPtr, EoS_AuxArray_Flt,
                                  EoS_AuxArray_Int, h_EoS_Table, __FILE__,  __LINE__, __FUNCTION__, UNPHY_VERBOSE ) )
            Aux_Error( ERROR_INFO, "Unphysical cell at (%e, %e %e)\n", x, y, z);

         Hydro_Pri2Con( Pri, fluid, false, PassiveNorm_NVar, PassiveNorm_VarIdx, EoS_DensPres2Eint_CPUPtr,
                        EoS_Temp2HTilde_CPUPtr, EoS_HTilde2Temp_CPUPtr, EoS_AuxArray_Flt, EoS_AuxArray_Int,
                        h_EoS_Table, NULL );


#        if (NCOMP_PASSIVE_USER > 0)
         fluid[Passive_0000] = 0.0;
         fluid[Passive_0001] = 0.0;
         fluid[Passive_0002] = 0.0;
#        endif

#        ifdef COSMIC_RAY
         fluid[CRAY] = Amb_CR_Engy * SQRT((real)1.0+SQR(Jet_SrcVel));
#        endif

      } // if ( fabs(zc) < interfaceHeight ) ... else ...
   } // else if ( Jet_Ambient == 3 )
   else if ( Jet_Ambient == 4 )
   {
      real ambientDens;
//
//    PotAtZ0 = IsothermalSlab_Pot(interfaceHeight);
//
//    Dens_gDisk_ambient  = ambientTemperature / gasDiskTemperature;
//    Dens_gDisk_ambient *= exp( PotAtZ0*(ambientTemperature-gasDiskTemperature)/(ambientTemperature*gasDiskTemperature) );
//
//    if (Dens_gDisk_ambient > HUGE_NUMBER || Dens_gDisk_ambient < -HUGE_NUMBER){
//      printf("Dens_gDisk_ambient=%e! %s: %d\n", Dens_gDisk_ambient, __FUNCTION__, __LINE__);
//      exit(0);
//    }

      real ambientPeakDens  = (real)2.842783e-27/UNIT_D;


      if ( fabs(zc) > IsothermalSlab_Truncation )
         ambientDens  = -IsothermalSlab_Pot(IsothermalSlab_Truncation)/ambientTemperature;
      else
         ambientDens  = -IsothermalSlab_Pot(zc)/ambientTemperature;

      //printf("IsothermalSlab_Pot(zc)=%e\n", IsothermalSlab_Pot(zc));
      //printf("ambientTemperature=%e\n", ambientTemperature);

      ambientDens  = exp(ambientDens);
      ambientDens *= ambientPeakDens;

      Pri[0] = ambientDens;
      Pri[1] = 0.0;
      Pri[2] = 0.0;
      Pri[3] = 0.0;
      Pri[4] = ambientDens*ambientTemperature;

      Hydro_Pri2Con( Pri, fluid, false, PassiveNorm_NVar, PassiveNorm_VarIdx, EoS_DensPres2Eint_CPUPtr,
                     EoS_Temp2HTilde_CPUPtr, EoS_HTilde2Temp_CPUPtr, EoS_AuxArray_Flt, EoS_AuxArray_Int,
                     h_EoS_Table, NULL );


      if ( Hydro_IsUnphysical( UNPHY_MODE_PRIM, Pri, NULL, NULL, NULL, NULL, EoS_DensEint2Pres_CPUPtr,
                               EoS_GuessHTilde_CPUPtr, EoS_HTilde2Temp_CPUPtr, EoS_AuxArray_Flt,
                               EoS_AuxArray_Int, h_EoS_Table, __FILE__,  __LINE__, __FUNCTION__, UNPHY_VERBOSE ) )
            Aux_Error( ERROR_INFO, "Unphysical cell at (%e, %e %e)\n", x, y, z);

#     if (NCOMP_PASSIVE_USER > 0)
      fluid[Passive_0000] = 0.0;
      fluid[Passive_0001] = 0.0;
      fluid[Passive_0002] = 0.0;
#     endif

#     ifdef COSMIC_RAY
      fluid[CRAY] = Amb_CR_Engy * SQRT((real)1.0+SQR(Jet_SrcVel));
#     endif

#     endif // #ifdef GRAVITY
   } else
   {
      // TODO : Add an error for wrong Jet_Ambient
   } // if ( Jet_Ambient == 0 ) ... else if ... else ...
} // FUNCTION : SetGridIC



//-------------------------------------------------------------------------------------------------------
// Function    :  Flu_ResetByUser_Jets
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
//
int Flu_ResetByUser_Jets( real fluid[], const double Emag, const double x, const double y, const double z,
                           const double Time, const double dt, const int lv, double AuxArray[] )
{
  if ( Jet_Fire == 0 ) return false;

  if ( Jet_Duration < Time ) return false;

  if ( !Jet_SphericalSrc ) {
     double xp[3], rp[3];
     double Prim[NCOMP_FLUID] = {0.0}, Cons[NCOMP_FLUID] = {0.0}, Vel[3];
     real PriReal[NCOMP_FLUID] = {0.0};
     double PrecessionAxis_Spherical[3], Omega_t;
     bool InsideUpperCone, InsideLowerCone;
     double Jet_SrcVelSmooth;

     Omega_t = Jet_AngularVelocity * Time * M_PI / 180.0;

//   shift the coordinate origin to the source center (the point O)
     xp[0] = x - Jet_Center[0];
     xp[1] = y - Jet_Center[1];
     xp[2] = z - Jet_Center[2];

     if ( Jet_PrecessionAxis[0] != 0.0  ||  Jet_PrecessionAxis[1] != 0.0  ||  Jet_PrecessionAxis[2] == 0.0 )
     {
//      get theta, phi for the first rotation
        Mis_Cartesian2Spherical( Jet_PrecessionAxis, PrecessionAxis_Spherical );
//      rotate coordinate to align z-axis with fixed precession axis
        CartesianRotate( xp, PrecessionAxis_Spherical[1], PrecessionAxis_Spherical[2], false );
     }

//   rotate coordinate to align z-axis with rotating symmetric axis
     CartesianRotate( xp, Jet_PrecessionAngle, Omega_t, false );

//   determine whether or not the point is inside of source
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
     } else
     {
        InsideUpperCone &= 0.0 <= xp[2] && xp[2] <= Jet_HalfHeight;
        InsideLowerCone &= -Jet_HalfHeight <= xp[2] && xp[2] <= 0.0;
     } // if ( Jet_HalfOpeningAngle != 0.0 ) ... else ...

//   set fluid variable inside source
//   TODO: try to define Jet_upper/lower to avoid the ugly code
     if ( ( InsideUpperCone  &&  ( Jet_Fire == 1  ||  Jet_Fire == 3 ) )  ||
          ( InsideLowerCone  &&  ( Jet_Fire == 2  ||  Jet_Fire == 3 ) ) )
     {
        if ( Jet_HalfOpeningAngle == 0.0 )
        {
           Vel[0] = 0.0;
           Vel[1] = 0.0;
           if ( InsideUpperCone == true ) Vel[2] = +Jet_SrcVel;
           else                           Vel[2] = -Jet_SrcVel;

           CartesianRotate( Vel, Jet_PrecessionAngle, Omega_t, true );

           if ( Jet_PrecessionAxis[0] != 0.0  ||  Jet_PrecessionAxis[1] != 0.0  ||  Jet_PrecessionAxis[2] == 0.0 )
              CartesianRotate( Vel, PrecessionAxis_Spherical[1], PrecessionAxis_Spherical[2], true );

           Prim[0] = Jet_SrcDens;
           Prim[1] = Vel[0];
           Prim[2] = Vel[1];
           Prim[3] = Vel[2];
           Prim[4] = Jet_SrcTemp*Jet_SrcDens;
        } else // if ( Jet_HalfOpeningAngle == 0.0 )
        {
//         shift origin to the point D/E
           if ( InsideUpperCone == true )   xp[2] += Jet_Radius/tan(Jet_HalfOpeningAngle);
           else                             xp[2] -= Jet_Radius/tan(Jet_HalfOpeningAngle);

           CartesianRotate( xp, Jet_PrecessionAngle, Omega_t, true );

           Mis_Cartesian2Spherical( xp, rp );

           if ( InsideLowerCone == true ) rp[1] -= M_PI;

//         smooth velocity on cross section
           if ( Jet_SmoothVel ) Jet_SrcVelSmooth = Jet_SrcVel*SQR(cos( 0.5 * M_PI * rp[1] / Jet_HalfOpeningAngle ));
           else                 Jet_SrcVelSmooth = Jet_SrcVel;

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

#       if (NCOMP_PASSIVE_USER > 0)
        fluid[Passive_0000] = 0.0;
        fluid[Passive_0001] = fluid[DENS];
        fluid[Passive_0002] = Jet_Src_CR_Engy * SQRT((real)1.0+SQR(Jet_SrcVel));
#       endif

#       ifdef COSMIC_RAY
        fluid[CRAY] = Jet_Src_CR_Engy * SQRT((real)1.0+SQR(Jet_SrcVel));
#       endif

        return true;
     } // if ( ( InsideUpperCone  &&  ( Jet_Fire == 1  ||  Jet_Fire == 3 ) )  ||
//             ( InsideLowerCone  &&  ( Jet_Fire == 2  ||  Jet_Fire == 3 ) ) )
  } else // if ( !Jet_SphericalSrc )
  {
     double xp[3], rp[3];
     double Prim[NCOMP_FLUID] = {0.0}, Cons[NCOMP_FLUID] = {0.0}, Vel[3];
     real PriReal[NCOMP_FLUID] = {0.0};

//   shift the coordinate origin to the source center (the point O)
     xp[0] = x - Jet_Center[0];
     xp[1] = y - Jet_Center[1];
     xp[2] = z - Jet_Center[2];

     //Mis_Cartesian2Spherical(xp, rp);

     double R = SQRT( SQR(xp[0]) + SQR(xp[1]) + SQR(xp[2]) );

//   For quick test, we temporarily use Jet_HalfHeight to represent the radius of spherical source
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

#       if (NCOMP_PASSIVE_USER > 0)
        fluid[Passive_0000] = 0.0;
        fluid[Passive_0001] = fluid[DENS];
        fluid[Passive_0002] = Jet_Src_CR_Engy * SQRT((real)1.0+SQR(Jet_SrcVel));
#       endif

#       ifdef COSMIC_RAY
        fluid[CRAY] = Jet_Src_CR_Engy * SQRT((real)1.0+SQR(Jet_SrcVel));
#       endif

        return true;
     } // if ( R < Jet_HalfHeight )
  } // if ( !Jet_SphericalSrc ) ... else ...

  return false;

} // FUNCTION : Flu_ResetByUser_Jets



// (true/false): if the target cell (is/is not) within the region to be refined
static bool Flag_Region( const int i, const int j, const int k, const int lv, const int PID )
{
   if ( Step == 0 )   return true; // TODO: figure why refine at first, or why the code is like this

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

// TODO the last condition seems weird
   bool Flag = R > gasDisk_highResRadius  ||  lv > gasDisk_lowRes_LEVEL  ||  fabs(dr[2]) < 2.0*interfaceHeight;

   return Flag;
//   original flag. TODO: once the question up is solved, remove it
//   bool Flag = R > gasDisk_highResRadius && lv > gasDisk_lowRes_LEVEL && fabs(dr[2]) < 2.0*interfaceHeight;
//   if ( Flag ) return false;
//   else        return true;
} // FUNCTION : Flag_Region



bool Flag_User( const int i, const int j, const int k, const int lv, const int PID, const double *Threshold )
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
      if (lv >= jetSrc_lowRes_LEVEL) Disk = false;
      Flag = Src || Disk;
   }
   else
   {
      Flag = Src;
   }
   return Flag;
} // FUNCTION : Flag_User



void CartesianRotate( double x[], const double theta, const double phi, const bool inverse )
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
// Function    :  Mis_GetTimeStep_User_Template
// Description :  Template of user-defined criteria to estimate the evolution time-step
//
// Note        :  1. This function should be applied to both physical and comoving coordinates and always
//                   return the evolution time-step (dt) actually used in various solvers
//                   --> Physical coordinates : dt = physical time interval
//                       Comoving coordinates : dt = delta(scale_factor) / ( Hubble_parameter*scale_factor^3 )
//                   --> We convert dt back to the physical time interval, which equals "delta(scale_factor)"
//                       in the comoving coordinates, in Mis_GetTimeStep()
//                2. Invoked by Mis_GetTimeStep() using the function pointer "Mis_GetTimeStep_User_Ptr",
//                   which must be set by a test problem initializer
//                3. Enabled by the runtime option "OPT__DT_USER"
//
// Parameter   :  lv       : Target refinement level
//                dTime_dt : dTime/dt (== 1.0 if COMOVING is off)
//
// Return      :  dt
//-------------------------------------------------------------------------------------------------------
double Mis_GetTimeStep_User( const int lv, const double dTime_dt )
{
   const double Jet_SrcGamma = sqrt(1.0 + SQR(Jet_SrcVel));
   const double Jet_Src3Vel = Jet_SrcVel / Jet_SrcGamma;

   const double dh  = amr->dh[MAX_LEVEL];

   const double Cs = 0.182574; // 1.0/sqrt(3); // TODO: not enough digit for double
   double dt_user  = DT__FLUID * dh / (Jet_Src3Vel+3.0*Cs);

   return dt_user;
} // FUNCTION : Mis_GetTimeStep_User_Template



void AddNewField_Jet()
{
#  if ( NCOMP_PASSIVE_USER > 0)
   if ( Passive_0000 == 5 ) Passive_0000 = AddField( "Passive_0000", FIXUP_FLUX_YES, FIXUP_REST_YES, NORMALIZE_NO, INTERP_FRAC_NO );
   if ( Passive_0001 == 6 ) Passive_0001 = AddField( "Passive_0001", FIXUP_FLUX_YES, FIXUP_REST_YES, NORMALIZE_NO, INTERP_FRAC_NO );
   if ( Passive_0002 == 7 ) Passive_0002 = AddField( "Passive_0002", FIXUP_FLUX_YES, FIXUP_REST_YES, NORMALIZE_NO, INTERP_FRAC_NO );
#  endif
} // FUNCTION : AddNewField_Jet



//-------------------------------------------------------------------------------------------------------
// Function    :  Init_TestProb_Hydro_Jet
// Description :  Test problem initializer
//
// Note        :  None
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Init_TestProb_Hydro_Jet()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );


// validate the compilation flags and runtime parameters
   Validate();


// set the problem-specific runtime parameters
   SetParameter();


// get enclosed mass
   Init_Function_User_Ptr   = SetGridIC;
   Flag_User_Ptr            = Flag_User;
   Flag_Region_Ptr          = Flag_Region;

   //if (Jet_Fire > 0)
   //{
   //  OPT__DT_USER             = 1;
   //  Mis_GetTimeStep_User_Ptr = Mis_GetTimeStep_User;
   //}
   //else
   //{
   //  OPT__DT_USER             = 0;
   //  Mis_GetTimeStep_User_Ptr = NULL;
   //}

   BC_User_Ptr              = NULL;
   Flu_ResetByUser_Func_Ptr = Flu_ResetByUser_Jets;
   Output_User_Ptr          = NULL;
   Aux_Record_User_Ptr      = NULL;
   End_User_Ptr             = NULL;
#  ifdef GRAVITY
   Init_ExtPot_Ptr          = Init_ExtPot_IsothermalSlab;
#  endif

#  if ( NCOMP_PASSIVE_USER > 0 )
   Init_Field_User_Ptr      = AddNewField_Jet;
#  endif


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : Init_TestProb_Hydro_Jets

#endif // #if ( MODEL == HYDRO )
