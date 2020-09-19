#include <random>
#include "GAMER.h"
#include "TestProb.h"

void Cartesian2Spherical( const double x[], double r[] );
void Spherical2Cartesian( const double r[], double x[] );
void CartesianRotate( double x[], double theta, double phi, bool inverse );

#if ( MODEL == SR_HYDRO )

// problem-specific global variables
// =======================================================================================

// options
static int      Jet_Ambient;                                    // [0/1/2/3] : uniform/iothermal/singular isothermal/beta model
static bool     Jet_DiffPrecession;                             //
static bool     Jet_Precession;                                 //
static bool     Jet_TimeDependentSrc;                           //
static int      Jet_Exist;                                      // [0/1/2/3]: no jet/jet1/jet2/full jet

// general parameters
static double   Jet_ParticleMassSrc;                            // particle mass in jet source [g]
static double   Jet_ParticleMassAmbient;                        // particle mass in ambient    [g]

// uniform background parameters
static double   Jet_Rd_UpperBound;                              // amplitude of random number
static double   Jet_UniformDens;                                // uniform ambient density
static double   Jet_UniformVel[3];                              // uniform ambient 4-velocity
static double   Jet_UniformTemp;                                // uniform ambient temperature

// isothermal sphere parameters
static double   Jet_HSE_Dx                         = NULL_REAL; // for Jet_Ambient=1: distance between box center and jet source (along x)
static double   Jet_HSE_Dy                         = NULL_REAL; // for Jet_Ambient=1: distance between box center and jet source (along y)
static double   Jet_HSE_Dz                         = NULL_REAL; // for Jet_Ambient=1: distance between box center and jet source (along z)
static char     Jet_HSE_BgTable_File[MAX_STRING];               // for Jet_Ambient=1: filename of the background gas table
static double  *Jet_HSE_BgTable_Data               = NULL;      // for Jet_Ambient=1: background gas table [radius/density/temperature]
static int      Jet_HSE_BgTable_NBin;                           // for Jet_Ambient=1: number of bins in Jet_HSE_BgTable_Data[]
static double   Jet_HSE_Radius                     = NULL_REAL; // for Jet_Ambient=1: radius of halo
static double   Jet_HSE_Mass                       = NULL_REAL; // for Jet_Ambient=1: enclosed mass of halo

// beta model
static double   Jet_Beta_Beta;
static double   Jet_Beta_Rcore;
static double   Jet_Beta_PeakDens;

// jet fluid parameters
static double   Jet_SrcVel;                                     // jet 4-velocity
static double   Jet_SrcDens;                                    // jet density
static double   Jet_SrcTemp;                                    // jet temperature
static bool     Jet_SmoothVel;


// sound speed
static double Cs;
static double CrossingTime;

// =======================================================================================
//        G       A       C
//          _____________
//          \     |     /
//           \   E|    /        z
//            \  /|\  /         ^
//             \/_|_\/          |
//             /\O| /\B         |
//            /  \|/  \
//           /   D|    \
//          /_____|_____\
//                F
// =======================================================================================
//
// jet geometry parameters
static double   Jet_Radius;                                     // OB
static double   Jet_HalfHeight;                                 // OA
static double   Jet_HalfOpeningAngle;                           // half-opening angle
static double   Jet_CenOffset[3];                               // jet central coordinates offset
static double   Jet_Cen[3];                                     // jet central coordinates
static double   Jet_MaxDis;                                     // maximum distance between the cylinder-shape jet source and the jet center

// precession parameters
static double   Jet_Angular_Velocity;                           // precession angular velocity (degree per code_time)
static double   Jet_PrecessionAngle;                            // precession angle in degree
static double   Jet_PrecessionAxis[3];                          // cone orientation vector (x,y,z) (NOT necessary to be a unit vector)

// time-depent source
static double   Jet_BurstStartTime;                             // start burst time in jet
static double   Jet_BurstEndTime;                               // end burst time in jet
static double   Jet_Burst4Vel;                                  // burst 4-velocity
static double   Jet_BurstDens;                                  // burst proper density
static double   Jet_BurstTemp;                                  // burst temperature
static bool     Flag_Burst4Vel;
static bool     Flag_BurstDens;
static bool     Flag_BurstTemp;

// pressure-balanced sphere
static double   Sphere_Radius;
static double   Sphere_CoreDens;
static double   Sphere_CoreRadius;
static double   Sphere_Center_x;
static double   Sphere_Center_y;
static double   Sphere_Center_z;

static double *Table_R, *Table_D, *Table_T;

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
#  if ( MODEL != SR_HYDRO )
   Aux_Error( ERROR_INFO, "MODEL must be SR_HYDRO !!\n" );
#  endif


#  ifdef PARTICLE
   Aux_Error( ERROR_INFO, "PARTICLE must be disabled !!\n" );
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
// ********************************************************************************************************************************
// ReadPara->Add( "KEY_IN_THE_FILE",   &VARIABLE,              DEFAULT,       MIN,              MAX               );
// ********************************************************************************************************************************

// load options
   ReadPara->Add( "Jet_Exist",               &Jet_Exist,                  3,                     0,              3    );
   ReadPara->Add( "Jet_Ambient",             &Jet_Ambient,                0,                     0,              3    );
   ReadPara->Add( "Jet_DiffPrecession",      &Jet_DiffPrecession,     false,          Useless_bool,   Useless_bool    );
   ReadPara->Add( "Jet_Precession",          &Jet_Precession,         false,          Useless_bool,   Useless_bool    );
   ReadPara->Add( "Jet_TimeDependentSrc",    &Jet_TimeDependentSrc,   false,          Useless_bool,   Useless_bool    );

// load jet fluid parameters
   ReadPara->Add( "Jet_SrcVel",              &Jet_SrcVel    ,          -1.0,          NoMin_double,   NoMax_double    );
   ReadPara->Add( "Jet_SmoothVel",           &Jet_SmoothVel ,          true,          Useless_bool,   Useless_bool    );
   ReadPara->Add( "Jet_SrcDens",             &Jet_SrcDens   ,          -1.0,          Eps_double,     NoMax_double    );
   ReadPara->Add( "Jet_SrcTemp",             &Jet_SrcTemp   ,          -1.0,          Eps_double,     NoMax_double    );
   ReadPara->Add( "Jet_ParticleMassSrc",     &Jet_ParticleMassSrc,     -1.0,          Eps_double,     NoMax_double    );
   ReadPara->Add( "Jet_ParticleMassAmbient", &Jet_ParticleMassAmbient, -1.0,          Eps_double,     NoMax_double    );

// load source geometry parameters
   ReadPara->Add( "Jet_Radius",              &Jet_Radius,              -1.0,          Eps_double,     NoMax_double    );
   ReadPara->Add( "Jet_HalfHeight",          &Jet_HalfHeight,          -1.0,          Eps_double,     NoMax_double    );
   ReadPara->Add( "Jet_HalfOpeningAngle",    &Jet_HalfOpeningAngle,    -1.0,                   0.0,           90.0    );
   ReadPara->Add( "Jet_CenOffset_x",         &Jet_CenOffset [0],        NoDef_double, NoMin_double,   NoMax_double    );
   ReadPara->Add( "Jet_CenOffset_y",         &Jet_CenOffset [1],        NoDef_double, NoMin_double,   NoMax_double    );
   ReadPara->Add( "Jet_CenOffset_z",         &Jet_CenOffset [2],        NoDef_double, NoMin_double,   NoMax_double    );

// load precission parameters
   ReadPara->Add( "Jet_Angular_Velocity",    &Jet_Angular_Velocity,     NoDef_double,          0.0,   NoMax_double    );
   ReadPara->Add( "Jet_PrecessionAngle",     &Jet_PrecessionAngle,      NoDef_double, NoMin_double,           90.0    );
   ReadPara->Add( "Jet_PrecessionAxis_x",    &Jet_PrecessionAxis[0],    NoDef_double, NoMin_double,   NoMax_double    );               
   ReadPara->Add( "Jet_PrecessionAxis_y",    &Jet_PrecessionAxis[1],    NoDef_double, NoMin_double,   NoMax_double    );  
   ReadPara->Add( "Jet_PrecessionAxis_z",    &Jet_PrecessionAxis[2],    NoDef_double, NoMin_double,   NoMax_double    );  

// load uniform background parameters
   ReadPara->Add( "Jet_UniformDens",         &Jet_UniformDens,         -1.0,          Eps_double,     NoMax_double    );
   ReadPara->Add( "Jet_Rd_UpperBound",       &Jet_Rd_UpperBound,        0.0,                   0.0,   NoMax_double    );
   ReadPara->Add( "Jet_UniformVel_x",        &Jet_UniformVel[0],        0.0,          NoMin_double,   NoMax_double    );
   ReadPara->Add( "Jet_UniformVel_y",        &Jet_UniformVel[1],        0.0,          NoMin_double,   NoMax_double    );
   ReadPara->Add( "Jet_UniformVel_z",        &Jet_UniformVel[2],        0.0,          NoMin_double,   NoMax_double    );
   ReadPara->Add( "Jet_UniformTemp",         &Jet_UniformTemp,         -1.0,          Eps_double,     NoMax_double    );

// load isothermal sphere parameters
   ReadPara->Add( "Jet_HSE_Radius",          &Jet_HSE_Radius,          -1.0,          NoMin_double,   NoMax_double    );
   ReadPara->Add( "Jet_HSE_Dx",              &Jet_HSE_Dx,               0.0,          NoMin_double,   NoMax_double    );
   ReadPara->Add( "Jet_HSE_Dy",              &Jet_HSE_Dy,               0.0,          NoMin_double,   NoMax_double    );
   ReadPara->Add( "Jet_HSE_Dz",              &Jet_HSE_Dz,               0.0,          NoMin_double,   NoMax_double    );
   ReadPara->Add( "Jet_HSE_BgTable_File",     Jet_HSE_BgTable_File,     Useless_str,  Useless_str,    Useless_str     );

// load beta model parameters
   ReadPara->Add( "Jet_Beta_Beta",           &Jet_Beta_Beta,            0.0,          NoMin_double,   NoMax_double    );
   ReadPara->Add( "Jet_Beta_Rcore",          &Jet_Beta_Rcore,           0.0,          NoMin_double,   NoMax_double    );
   ReadPara->Add( "Jet_Beta_PeakDens",       &Jet_Beta_PeakDens,        0.0,          NoMin_double,   NoMax_double    );

// load pressured-balanced sphere
   ReadPara->Add( "Sphere_Radius",           &Sphere_Radius,           -1.0,          NoMin_double,   NoMax_double    );
   ReadPara->Add( "Sphere_CoreDens",         &Sphere_CoreDens,          0.0,          Eps_double,     NoMax_double    );
   ReadPara->Add( "Sphere_CoreRadius",       &Sphere_CoreRadius,        0.0,          Eps_double,     NoMax_double    );
   ReadPara->Add( "Sphere_Center_x",         &Sphere_Center_x,          0.0,          NoMin_double,   NoMax_double    );
   ReadPara->Add( "Sphere_Center_y",         &Sphere_Center_y,          0.0,          NoMin_double,   NoMax_double    );
   ReadPara->Add( "Sphere_Center_z",         &Sphere_Center_z,          0.0,          NoMin_double,   NoMax_double    );

// load time-dependent source varibles
   ReadPara->Add( "Jet_BurstStartTime",      &Jet_BurstStartTime,      -1.0,          NoMin_double,   NoMax_double    );
   ReadPara->Add( "Jet_BurstEndTime",        &Jet_BurstEndTime,        -1.0,          NoMin_double,   NoMax_double    );
   ReadPara->Add( "Jet_Burst4Vel",           &Jet_Burst4Vel,           -1.0,          NoMin_double,   NoMax_double    );
   ReadPara->Add( "Jet_BurstDens",           &Jet_BurstDens,           -1.0,          NoMin_double,   NoMax_double    );
   ReadPara->Add( "Jet_BurstTemp",           &Jet_BurstTemp,           -1.0,          NoMin_double,   NoMax_double    );
   ReadPara->Add( "Flag_Burst4Vel",          &Flag_Burst4Vel,           false,        Useless_bool,   Useless_bool    );
   ReadPara->Add( "Flag_BurstDens",          &Flag_BurstDens,           false,        Useless_bool,   Useless_bool    );
   ReadPara->Add( "Flag_BurstTemp",          &Flag_BurstTemp,           false,        Useless_bool,   Useless_bool    );

   ReadPara->Read( FileName );

   delete ReadPara;

// replace useless parameters with NaN
   if ( Jet_Ambient != 0 )
   {
     Jet_UniformDens      = NAN;
     Jet_UniformVel[0]    = NAN;
     Jet_UniformVel[1]    = NAN;
     Jet_UniformVel[2]    = NAN;
   }


   // uniform precession
   if ( Jet_Precession )       
   {
   }

   if ( Sphere_Radius < 0.0 )
   {
      Sphere_Radius     = NAN;
      Sphere_CoreDens   = NAN;
      Sphere_CoreRadius = NAN;
      Sphere_Center_x   = NAN; 
      Sphere_Center_y   = NAN;
      Sphere_Center_z   = NAN;   
   }

   // differential precession
   if ( Jet_DiffPrecession )  
   {
   }

   // no precession
   if ( !Jet_Precession && !Jet_DiffPrecession ) 
   {
   }
   
   if ( !Jet_TimeDependentSrc )
   {
     Jet_BurstDens       = NAN;
     Jet_Burst4Vel       = NAN;
     Jet_BurstTemp       = NAN; 
     Jet_BurstStartTime  = NAN;
     Jet_BurstEndTime    = NAN;
   }

// (1-2) check runtime parameters

// check time-dependent source
   if ( Jet_TimeDependentSrc )
   {
     if ( !Flag_Burst4Vel && !Flag_BurstDens && !Flag_BurstTemp )
       Aux_Error( ERROR_INFO, "One of Flag_Burst4Vel, Flag_BurstDens or Flag_BurstTemp must be enabled !!\n" );

     if ( Jet_BurstEndTime <= Jet_BurstStartTime)
       Aux_Error( ERROR_INFO, "Jet_BurstEndTime <= Jet_BurstStartTime !!\n" );
  
     if ( Jet_BurstEndTime >= END_T )
       Aux_Error( ERROR_INFO, "Jet_BurstEndTime >= END_T !!\n" );
  
     if ( Flag_Burst4Vel && Jet_Burst4Vel <= Eps_double )
       Aux_Error( ERROR_INFO, "Jet_Burst4Vel <= Eps_double !!\n" );
  
     if ( Flag_BurstDens && Jet_BurstDens <= Eps_double )
       Aux_Error( ERROR_INFO, "Jet_BurstDens <= Eps_double !!\n" );
  
     if ( Flag_BurstTemp && Jet_BurstTemp <= Eps_double )
       Aux_Error( ERROR_INFO, "Jet_BurstTemp <= Eps_double !!\n" );
   }


// check isothermal sphere
   if ( Jet_HSE_Radius > 0.0 )
   {
      if ( !OPT__FLAG_REGION && OPT__FLAG_LOHNER_DENS && OPT__FLAG_LOHNER_TEMP && OPT__FLAG_LOHNER_ENGY ) 
        Aux_Error( ERROR_INFO, "OPT__FLAG_REGION must be enabled !!\n" );

      if ( Jet_HSE_Radius == NULL_REAL ) 
        Aux_Error( ERROR_INFO, "Jet_HSE_Radius = %e !!\n" );

      if ( Jet_HSE_Radius >= 0.5*amr->BoxSize[0] || 
           Jet_HSE_Radius >= 0.5*amr->BoxSize[1] || 
           Jet_HSE_Radius >= 0.5*amr->BoxSize[2]   )
        Aux_Error( ERROR_INFO, "halo size is greater or equal to box size !!\n" );

#      ifndef GRAVITY
       Aux_Error( ERROR_INFO, "GRAVITY must be enabled !!\n" );
#      endif
   }


// check uniform ambient
   else if ( Jet_Ambient == 0 )
   {
#     ifdef GRAVITY
      Aux_Error( ERROR_INFO, "GRAVITY must be disabled !!\n" );
#     endif
   }

// check beta model
   else if ( Jet_Ambient == 3 )
   {
	  if ( Jet_Beta_Beta == 0.0 )
      Aux_Error( ERROR_INFO, "Jet_Beta_Beta = 0.0 !!\n" );

	  if ( Jet_Beta_Rcore == 0.0 )
      Aux_Error( ERROR_INFO, "Jet_Beta_Rcore = 0.0 !!\n" );

	  if ( Jet_Beta_PeakDens == 0.0 )
      Aux_Error( ERROR_INFO, "Jet_Beta_PeakDens = 0.0 !!\n" );
   }

// check precession parameters
   if ( Jet_Precession && Jet_DiffPrecession )
      Aux_Error( ERROR_INFO, "Jet_Precession and Jet_DiffPrecession are not allowed to enable simultaneously !!\n" );

// (1-2) convert to code unit
   Jet_SrcVel               *= Const_c   / UNIT_V;
   Jet_SrcTemp              *= Const_GeV / ( Jet_ParticleMassSrc * SQR(Const_c) );
   Jet_SrcDens              *= 1.0       / UNIT_D;

   Jet_Radius               *= Const_pc / UNIT_L;
   Jet_HalfHeight           *= Const_pc / UNIT_L;
   Jet_HalfOpeningAngle     *= M_PI / 180.0;
   Jet_PrecessionAngle      *= M_PI / 180.0;

   Jet_CenOffset[0]         *= Const_pc / UNIT_L;
   Jet_CenOffset[1]         *= Const_pc / UNIT_L;
   Jet_CenOffset[2]         *= Const_pc / UNIT_L;

   Jet_HSE_Dx               *= Const_pc / UNIT_L;  
   Jet_HSE_Dy               *= Const_pc / UNIT_L;  
   Jet_HSE_Dz               *= Const_pc / UNIT_L;  

   if ( Jet_HSE_Radius > 0.0 )
   {
     Jet_HSE_Radius         *= Const_pc / UNIT_L; 
   }

   if ( Jet_Ambient == 0 )
   {
     Jet_UniformDens        *= 1.0       / UNIT_D;
     Jet_UniformVel[0]      *= Const_c   / UNIT_V;
     Jet_UniformVel[1]      *= Const_c   / UNIT_V;
     Jet_UniformVel[2]      *= Const_c   / UNIT_V;
   }
   
   if ( Jet_Ambient == 3 )
   {
     Jet_Beta_Rcore         *= Const_pc  / UNIT_L;
     Jet_Beta_PeakDens      *= 1.0       / UNIT_D;
   }

   Jet_UniformTemp          *= Const_GeV / ( Jet_ParticleMassAmbient * SQR(Const_c) );
   Jet_Angular_Velocity     *= 1.0;    // the unit of Jet_Angular_Velocity is UNIT_T

   
   if ( Sphere_Radius > 0.0 )
   {
      Sphere_Radius         *= Const_pc / UNIT_L;
      Sphere_CoreRadius     *= Const_pc / UNIT_L;
      Sphere_CoreDens       *=      1.0 / UNIT_D;
      Sphere_Center_x       *= Const_pc / UNIT_L; 
      Sphere_Center_y       *= Const_pc / UNIT_L;
      Sphere_Center_z       *= Const_pc / UNIT_L;   
   }

   if ( Jet_Precession ) // uniform precession
   {
   }

   if ( Jet_DiffPrecession ) // differential precession
   {
   }

   if ( !Jet_Precession && !Jet_DiffPrecession ) // no precession
   {
   }

   if ( Jet_TimeDependentSrc )
   {
     Jet_BurstStartTime     *= 1e3*Const_yr / UNIT_T;
     Jet_BurstEndTime       *= 1e3*Const_yr / UNIT_T;
     Jet_Burst4Vel          *= Const_c      / UNIT_V;
     Jet_BurstDens          *= 1.0          / UNIT_D;
     Jet_BurstTemp          *= Const_GeV    / ( Jet_ParticleMassSrc * SQR(Const_c) );
   }

// (2) set the problem-specific derived parameters
   double SecAngle = 1.0 / cos(0.5*Jet_HalfOpeningAngle);
   double TanAngle = sin(0.5*Jet_HalfOpeningAngle) * SecAngle;

   Jet_MaxDis  = sqrt( SQR( Jet_Radius ) + SQR( Jet_HalfHeight * SecAngle ) + 2.0 * Jet_Radius * Jet_HalfHeight * TanAngle );

   for (int d=0; d<3; d++)    Jet_Cen[d] = 0.5*amr->BoxSize[d] + Jet_CenOffset[d];



// (3) load the hydrostatic equilibrium table
   if ( Jet_Ambient == 1  &&  OPT__INIT != INIT_BY_RESTART )
   {
      const bool RowMajor_No  = false;       // load data into the column-major order
      const bool AllocMem_Yes = true;        // allocate memory for Merger_Prof1/2
      const int  NCol         = 3;           // total number of columns to load
      const int  Col[NCol]    = {1, 2, 3};   // target columns: (radius, density, temperature)

      Jet_HSE_BgTable_NBin = Aux_LoadTable( Jet_HSE_BgTable_Data, Jet_HSE_BgTable_File, NCol, Col, RowMajor_No, AllocMem_Yes );

      if ( Jet_HSE_BgTable_Data == NULL )
      {
         Aux_Error( ERROR_INFO, "Jet_HSE_BgTable_Data == NULL !!\n" );
      }

//    convert to code units
      Table_R = Jet_HSE_BgTable_Data + 0*Jet_HSE_BgTable_NBin;
      Table_D = Jet_HSE_BgTable_Data + 1*Jet_HSE_BgTable_NBin;
      Table_T = Jet_HSE_BgTable_Data + 2*Jet_HSE_BgTable_NBin;

      for (int b=0; b<Jet_HSE_BgTable_NBin; b++)
      {
         Table_R[b] *= Const_pc / UNIT_L;
         Table_D[b] *= 1.0       / UNIT_D;
         Table_T[b] *= 1.0;
         Table_T[b] *= Const_GeV / ( Jet_ParticleMassAmbient * SQR(Const_c) );
      }
   } // if ( Jet_Ambient  &&  OPT__INIT != INIT_BY_RESTART )


// computing sound speed and half-crossing time
   double AmbientTemp, Distance;

   if      ( Jet_Ambient == 1 )       AmbientTemp = Table_T[0];
   else                               AmbientTemp = Jet_UniformTemp;

   Cs = sqrt( AmbientTemp * SQR(Const_c) );

   ( Jet_HSE_Radius <= 0.0 ) ? Distance = BOX_SIZE : Distance = Jet_HSE_Radius;

   CrossingTime = ( Distance * Const_pc ) / Cs;

// (4) reset other general-purpose parameters
//     --> a helper macro PRINT_WARNING is defined in TestProb.h
   const long   End_Step_Default = __INT_MAX__;
   const double End_T_Default    = BOX_SIZE / Jet_SrcVel / UNIT_T;

   if ( END_STEP < 0 ) {
      END_STEP = End_Step_Default;
      PRINT_WARNING( "END_STEP", END_STEP, FORMAT_LONG );
   }

   if ( END_T < 0.0 ) {
      END_T = End_T_Default;
      PRINT_WARNING( "END_T", END_T, FORMAT_REAL );
   }


// (4) make a note
   if ( MPI_Rank == 0 )
   {
     Aux_Message( stdout, "=============================================================================\n" );
     Aux_Message( stdout, "  test problem ID         = %d\n",                TESTPROB_ID                                     );
     Aux_Message( stdout, "  Jet_Exist               = %d\n",                Jet_Exist                                       );
     Aux_Message( stdout, "  Jet_SmoothVel           = %d\n",                Jet_SmoothVel                                   );
     Aux_Message( stdout, "  Jet_Ambient             = %d\n",                Jet_Ambient                                     );
     Aux_Message( stdout, "  Jet_DiffPrecession      = %d\n",                Jet_DiffPrecession                              );
     Aux_Message( stdout, "  Jet_Precession          = %d\n",                Jet_Precession                                  );
     Aux_Message( stdout, "  Jet_TimeDependentSrc    = %d\n",                Jet_TimeDependentSrc                            );
     Aux_Message( stdout, "  Jet_ParticleMassSrc     = %14.7e g\n",          Jet_ParticleMassSrc                             );
     Aux_Message( stdout, "  Jet_ParticleMassAmbient = %14.7e g\n",          Jet_ParticleMassAmbient                         );
     Aux_Message( stdout, "  Jet_SrcVel              = %14.7e c\n",          Jet_SrcVel                                      );
     Aux_Message( stdout, "  Jet_SrcDens             = %14.7e g/cm^3\n",     Jet_SrcDens*UNIT_D                              );
     Aux_Message( stdout, "  Jet_SrcTemp             = %14.7e GeV\n",        Jet_SrcTemp*Jet_ParticleMassSrc*SQR(Const_c)/Const_GeV     );
     Aux_Message( stdout, "  Jet_NumDensSrc          = %14.7e per cc\n",     Jet_SrcDens*UNIT_D/Jet_ParticleMassSrc          );
     Aux_Message( stdout, "  Jet_CenOffset[x]        = %14.7e pc\n",         Jet_CenOffset [0]*UNIT_L/Const_pc               );
     Aux_Message( stdout, "  Jet_CenOffset[y]        = %14.7e pc\n",         Jet_CenOffset [1]*UNIT_L/Const_pc               );
     Aux_Message( stdout, "  Jet_CenOffset[z]        = %14.7e pc\n",         Jet_CenOffset [2]*UNIT_L/Const_pc               );
     Aux_Message( stdout, "  Jet_Angular_Velocity    = %14.7e degree/kyr\n", Jet_Angular_Velocity                            );
     Aux_Message( stdout, "  Jet_PrecessionAngle     = %14.7e degree\n",     Jet_PrecessionAngle*180.0/M_PI                  );
     Aux_Message( stdout, "  Jet_HalfOpeningAngle    = %14.7e degree\n",     Jet_HalfOpeningAngle*180.0/M_PI                 );
     Aux_Message( stdout, "  Jet_Radius              = %14.7e pc\n",         Jet_Radius*UNIT_L/Const_pc                      );
     Aux_Message( stdout, "  Jet_HalfHeight          = %14.7e pc\n",         Jet_HalfHeight*UNIT_L/Const_pc                  );
     Aux_Message( stdout, "  Jet_MaxDis              = %14.7e pc\n",         Jet_MaxDis*UNIT_L/Const_pc                      );
   }

   if ( Jet_Ambient == 0 && MPI_Rank == 0 )
   {
     Aux_Message( stdout, "  Jet_UniformDens         = %14.7e g/cm^3\n",     Jet_UniformDens*UNIT_D                          );
     Aux_Message( stdout, "  Jet_Rd_UpperBound       = %14.7e \n",           Jet_Rd_UpperBound                               );
     Aux_Message( stdout, "  Jet_UniformTemp         = %14.7e GeV\n",        Jet_UniformTemp*Jet_ParticleMassAmbient*SQR(Const_c)/Const_GeV );
     Aux_Message( stdout, "  Jet_UniformVel[x]       = %14.7e c\n",          Jet_UniformVel[0]                               );
     Aux_Message( stdout, "  Jet_UniformVel[y]       = %14.7e c\n",          Jet_UniformVel[1]                               );
     Aux_Message( stdout, "  Jet_UniformVel[z]       = %14.7e c\n",          Jet_UniformVel[2]                               );
     Aux_Message( stdout, "  Jet_UniformNumDens      = %14.7e per cc\n",     Jet_UniformDens*UNIT_D/Jet_ParticleMassAmbient  );
   }

   if ( MPI_Rank == 0 )
   {
     Aux_Message( stdout, "  Cs                      = %14.7e c\n",          Cs / Const_c                                    );
     Aux_Message( stdout, "  CrossingTime            = %14.7e pc/c\n",      CrossingTime / UNIT_T                            );
   }

   if ( Jet_Ambient == 1 && MPI_Rank == 0 )
   {
     Aux_Message( stdout, "  Jet_HSE_BgTable_File    = %s\n",                Jet_HSE_BgTable_File                            );
   }

   if ( Jet_HSE_Radius > 0.0 && MPI_Rank == 0 )
   {
     Aux_Message( stdout, "  Jet_HSE_Radius          = %14.7e pc\n",        Jet_HSE_Radius*UNIT_L/Const_pc                   );
     Aux_Message( stdout, "  Jet_HSE_Dx              = %14.7e pc\n",        Jet_HSE_Dx*UNIT_L/Const_pc                       );
     Aux_Message( stdout, "  Jet_HSE_Dy              = %14.7e pc\n",        Jet_HSE_Dy*UNIT_L/Const_pc                       );
     Aux_Message( stdout, "  Jet_HSE_Dz              = %14.7e pc\n",        Jet_HSE_Dz*UNIT_L/Const_pc                       );
   }

   if ( Jet_Ambient == 3 && MPI_Rank == 0 )
   {
     Aux_Message( stdout, "  Jet_Beta_Rcore          = %14.7e pc\n",        Jet_Beta_Rcore*UNIT_L/Const_pc                   );
     Aux_Message( stdout, "  Jet_Beta_PeakDens       = %14.7e g/cm^3\n",     Jet_Beta_PeakDens*UNIT_D                        );
     Aux_Message( stdout, "  Jet_Beta_Beta           = %14.7e\n",            Jet_Beta_Beta                                   );
   }

   if ( MPI_Rank == 0 )
   {
     Aux_Message( stdout, "  Jet_PrecessionAxis[x]   = %14.7e\n",            Jet_PrecessionAxis[0]                           );
     Aux_Message( stdout, "  Jet_PrecessionAxis[y]   = %14.7e\n",            Jet_PrecessionAxis[1]                           );
     Aux_Message( stdout, "  Jet_PrecessionAxis[z]   = %14.7e\n",            Jet_PrecessionAxis[2]                           );
   }

   if ( Sphere_Radius > 0.0 && MPI_Rank == 0 )
   {
     Aux_Message( stdout, "  Sphere_Radius           = %14.7e pc\n",        Sphere_Radius*UNIT_L/Const_pc                    );
     Aux_Message( stdout, "  Sphere_CoreRadius       = %14.7e pc\n",        Sphere_CoreRadius*UNIT_L/Const_pc                );
     Aux_Message( stdout, "  Sphere_CoreDens         = %14.7e g/cm^3\n",    Sphere_CoreDens*UNIT_D                           );
     Aux_Message( stdout, "  Sphere_DensSurface      = %14.7e g/cm^3\n",    Sphere_CoreDens / ( 1.0 + SQR( Sphere_Radius / Sphere_CoreRadius) )*UNIT_D );
     Aux_Message( stdout, "  Sphere_Center_x         = %14.7e pc\n",        Sphere_Center_x*UNIT_L/Const_pc                  );
     Aux_Message( stdout, "  Sphere_Center_y         = %14.7e pc\n",        Sphere_Center_y*UNIT_L/Const_pc                  );
     Aux_Message( stdout, "  Sphere_Center_z         = %14.7e pc\n",        Sphere_Center_z*UNIT_L/Const_pc                  );

   }

   if ( Jet_TimeDependentSrc && MPI_Rank == 0 )
   {
     Aux_Message( stdout, "  Jet_BurstStartTime      = %14.7e kyr \n",       Jet_BurstStartTime*UNIT_T/(1e3*Const_yr)        );
     Aux_Message( stdout, "  Jet_BurstEndTime        = %14.7e kyr \n",       Jet_BurstEndTime*UNIT_T/(1e3*Const_yr)          );
     Aux_Message( stdout, "  Jet_Burst4Vel           = %14.7e c \n",         Jet_Burst4Vel                                   );
     Aux_Message( stdout, "  Jet_BurstDens           = %14.7e g/cm^3\n",     Jet_BurstDens*UNIT_D                            );
     Aux_Message( stdout, "  Jet_BurstTemp           = %14.7e GeV\n",        Jet_BurstTemp*Jet_ParticleMassSrc*SQR(Const_c)/Const_GeV   );
     Aux_Message( stdout, "  Flag_Burst4Vel          = %d\n",                Flag_Burst4Vel                                  );
     Aux_Message( stdout, "  Flag_BurstDens          = %d\n",                Flag_BurstDens                                  );
     Aux_Message( stdout, "  Flag_BurstTemp          = %d\n",                Flag_BurstTemp                                  );
   }

   if ( MPI_Rank == 0 )
   {
     Aux_Message( stdout, "=============================================================================\n"                  );
   }


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Setting runtime parameters ... done\n"                                   );

} // FUNCTION : SetParameter


void GetEnclosedMass( double *Jet_HSE_Mass )
{
  double *Table_D, *Table_R;
  int N = 4096;
  double dr = Jet_HSE_Radius/(double)N;
  double r[N];
  const bool RowMajor_No  = false;       // load data into the column-major order
  const bool AllocMem_Yes = true;        // allocate memory for Merger_Prof1/2
  const int  NCol         = 3;           // total number of columns to load
  const int  Col[NCol]    = {1, 2, 3};   // target columns: (radius, density, temperature)


  if ( Jet_HSE_BgTable_Data == NULL )
  {
     Jet_HSE_BgTable_NBin = Aux_LoadTable( Jet_HSE_BgTable_Data, Jet_HSE_BgTable_File, NCol, Col, RowMajor_No, AllocMem_Yes );
     Aux_Error( ERROR_INFO, "Jet_HSE_BgTable_Data == NULL !!\n" );
  }
 
  Table_D = Jet_HSE_BgTable_Data + 1*Jet_HSE_BgTable_NBin;

  *Jet_HSE_Mass = 0.0;  

  Table_R = Jet_HSE_BgTable_Data + 0*Jet_HSE_BgTable_NBin;


// To reduce rounding error, we should accumulate mass from core
  for (int i=0; i<N; i++)
  {
    r[i] = dr * (double)(i+1);
    *Jet_HSE_Mass += dr*r[i]*r[i]*Mis_InterpolateFromTable( Jet_HSE_BgTable_NBin, Table_R, Table_D, r[i] );
  }

  *Jet_HSE_Mass *= 4.0*M_PI;
  
}


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

// variables for Jet_Ambient
//   const double *Table_R=NULL, *Table_D=NULL, *Table_T=NULL;
   double dx, dy, dz, r;
   const int  NCol         = 3;           // total number of columns to load
   const int  Col[NCol]    = {1, 2, 3};   // target columns: (radius, density, temperature)

   double DropRatio = 1e-2;

// variables for jet
   double Pri4Vel[NCOMP_FLUID];

   dx = x - amr->BoxCenter[0] - Jet_HSE_Dx;
   dy = y - amr->BoxCenter[1] - Jet_HSE_Dy;
   dz = z - amr->BoxCenter[2] - Jet_HSE_Dz;

   r = sqrt( dx*dx + dy*dy + dz*dz );

// random number generator
   std::random_device rd;

// Mersenne twister
   std::mt19937 generator( rd() );

// amplitude of random number

   std::uniform_real_distribution<real> unif(-Jet_Rd_UpperBound, Jet_Rd_UpperBound);
   double RandNumber = unif(generator);
   RandNumber = 0.0;

// uniform ambient
   if ( Jet_Ambient == 0 )
   {
	  Pri4Vel[0] = Jet_UniformDens * ( 1.0 + RandNumber / Jet_UniformDens );
      Pri4Vel[1] = Jet_UniformVel[0];
      Pri4Vel[2] = Jet_UniformVel[1];
      Pri4Vel[3] = Jet_UniformVel[2];
	  Pri4Vel[4] = Jet_UniformTemp * Jet_UniformDens;
   }

// isotherml sphere
   if ( Jet_Ambient == 1 )
   {
      double Temp, RhoSurface;

      Pri4Vel[1] = 0.0;
      Pri4Vel[2] = 0.0;
      Pri4Vel[3] = 0.0;


      if ( Jet_HSE_Radius > 0.0 )
      {
         if ( r <= Jet_HSE_Radius )
         {
           Pri4Vel[0] = Mis_InterpolateFromTable( Jet_HSE_BgTable_NBin, Table_R, Table_D, r );
           Temp = Mis_InterpolateFromTable( Jet_HSE_BgTable_NBin, Table_R, Table_T, r );
           Pri4Vel[4] = Temp * Pri4Vel[0];
         }
         else
         {
           Temp = Mis_InterpolateFromTable( Jet_HSE_BgTable_NBin, Table_R, Table_T, Jet_HSE_Radius );

           // uniform ambient outside halo
           Pri4Vel[0] =              DropRatio * Mis_InterpolateFromTable( Jet_HSE_BgTable_NBin, Table_R, Table_D, Jet_HSE_Radius );
           Pri4Vel[4] = Pri4Vel[0] / DropRatio * Temp;
         }
      }
      else
      {
         Pri4Vel[0] = Mis_InterpolateFromTable( Jet_HSE_BgTable_NBin, Table_R, Table_D, r );
         Temp = Mis_InterpolateFromTable( Jet_HSE_BgTable_NBin, Table_R, Table_T, r );
         Pri4Vel[4] = Temp * Pri4Vel[0];
      }
   }

   // singular isothermal sphere or beta model
   if ( Jet_Ambient == 2 ) 
   {
#     ifdef GRAVITY
      Pri4Vel[1] = 0.0;
      Pri4Vel[2] = 0.0;
      Pri4Vel[3] = 0.0;

      if ( Jet_HSE_Radius > 0.0 )
      {
        if ( r <= Jet_HSE_Radius )
        {
	      Pri4Vel[0]  = Jet_UniformTemp / (2.0*M_PI*NEWTON_G);
	      Pri4Vel[0] /= SQR( r );
	      Pri4Vel[4]  = Jet_UniformTemp * Pri4Vel[0];
        }
        else
        {
	      Pri4Vel[0]  = Jet_UniformTemp / (2.0*M_PI*NEWTON_G);
	      Pri4Vel[0] /= SQR( Jet_HSE_Radius );
          Pri4Vel[0] *= DropRatio;
          Pri4Vel[4]  = Pri4Vel[0] * Jet_UniformTemp / DropRatio;
        }
      }
      else
      {
	      Pri4Vel[0]  = Jet_UniformTemp / (2.0*M_PI*NEWTON_G);
	      Pri4Vel[0] /= SQR( r );
	      Pri4Vel[4]  = Jet_UniformTemp * Pri4Vel[0];
      }
#     endif
   }

   if ( Jet_Ambient == 3 ) 
   {
#     ifdef GRAVITY
      Pri4Vel[1] = 0.0;
      Pri4Vel[2] = 0.0;
      Pri4Vel[3] = 0.0;

      if ( Jet_HSE_Radius > 0.0 )
      {
        if ( r <= Jet_HSE_Radius )
        {
	      Pri4Vel[0]  = Jet_Beta_PeakDens * pow( 1.0 + (r/Jet_Beta_Rcore)*(r/Jet_Beta_Rcore), -1.5*Jet_Beta_Beta );
	      Pri4Vel[4]  = Jet_UniformTemp * Pri4Vel[0];
        }
        else
        {
	      Pri4Vel[0]  = Jet_Beta_PeakDens * pow( 1.0 + SQR( Jet_HSE_Radius/Jet_Beta_Rcore ), -1.5*Jet_Beta_Beta );
          Pri4Vel[0] *= DropRatio;
          Pri4Vel[4]  = Pri4Vel[0] * Jet_UniformTemp / DropRatio;
        }
      }
      else
      {
	      Pri4Vel[0]  = Jet_Beta_PeakDens * pow( 1.0 + (r/Jet_Beta_Rcore)*(r/Jet_Beta_Rcore), -1.5*Jet_Beta_Beta );
	      Pri4Vel[4]  = Jet_UniformTemp * Pri4Vel[0];
      }
#     endif
   }

   dx = x - amr->BoxCenter[0] - Sphere_Center_x;
   dy = y - amr->BoxCenter[1] - Sphere_Center_y;
   dz = z - amr->BoxCenter[2] - Sphere_Center_z;

   r = sqrt( dx*dx + dy*dy + dz*dz );

   if ( 0.0 < Sphere_Radius && Sphere_Radius > r )
   {
      Pri4Vel[0] = Sphere_CoreDens / ( 1.0 + SQR( r / Sphere_CoreRadius) );
   }



// cast double to real
#  ifndef FLOAT8
   double Out[NCOMP_FLUID];

   SRHydro_Pri2Con(Pri4Vel, Out, GAMMA);

   fluid [0] = (real) Out[0];
   fluid [1] = (real) Out[1];
   fluid [2] = (real) Out[2];
   fluid [3] = (real) Out[3];
   fluid [4] = (real) Out[4];
#  else
   SRHydro_Pri2Con(Pri4Vel, fluid, GAMMA);
#  endif

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
//        G       A       C
//          _____________
//          \     |     /
//           \   E|    /        z
//            \  /|\  /         ^
//             \/_|_\/          |
//             /\O| /\B         |
//            /  \|/  \
//           /   D|    \
//          /_____|_____\
//                F
//-------------------------------------------------------------------------------------------------------
bool Flu_ResetByUser_Jets( real fluid[], const double x, const double y, const double z, const double Time,
                                         const int lv, double AuxArray[] )
{
  if ( Jet_Exist == 0 ) return false;

  double xp[3], rp[3];
  double Prim[5], Cons[5], Vel[3];
  double PrecessionAxis_Spherical[3], Omega_t;
  bool InsideUpperCone, InsideLowerCone;
  double Jet_SrcVelSmooth;

  Omega_t = Jet_Angular_Velocity * Time * M_PI / 180.0;



  // shift the coordinate origin to the source center (the point O)
  xp[0] = x - Jet_Cen[0];
  xp[1] = y - Jet_Cen[1];
  xp[2] = z - Jet_Cen[2];

  if ( Jet_PrecessionAxis[0] != 0.0 || Jet_PrecessionAxis[1] != 0.0 ||  Jet_PrecessionAxis[2] == 0.0 )
  {
    // get theta, phi for the first rotation
    Cartesian2Spherical( Jet_PrecessionAxis, PrecessionAxis_Spherical );

    // rotate coordinate to align z-axis with fixed precession axis
	CartesianRotate(xp, PrecessionAxis_Spherical[1], PrecessionAxis_Spherical[2], false);
  }

  // rotate coordinate to align z-axis with rotating symmetric axis
  CartesianRotate(xp, Jet_PrecessionAngle, Omega_t, false);

  // determine whether or not the point is inside of source
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
    InsideUpperCone &= 0.0 <= xp[2] && xp[2] <= Jet_HalfHeight;
    InsideLowerCone &= -Jet_HalfHeight <= xp[2] && xp[2] <= 0.0;
  }



  // set fluid variable inside source
  if ( ( InsideUpperCone && ( Jet_Exist == 1 || Jet_Exist == 3 ) ) 
	|| ( InsideLowerCone && ( Jet_Exist == 2 || Jet_Exist == 3 ) ) )
  {
    if ( Jet_HalfOpeningAngle == 0.0 )
  	{
  	  Vel[0] = 0.0;
  	  Vel[1] = 0.0;
   
  	  if ( InsideUpperCone == true ) Vel[2] = +Jet_SrcVel;
      else                           Vel[2] = -Jet_SrcVel;
  
  	  CartesianRotate(Vel, Jet_PrecessionAngle, Omega_t, true);
  
      if ( Jet_PrecessionAxis[0] != 0.0 || Jet_PrecessionAxis[1] != 0.0 ||  Jet_PrecessionAxis[2] == 0.0 )
           CartesianRotate(Vel, PrecessionAxis_Spherical[1], PrecessionAxis_Spherical[2], true);
  
      Prim[0] = Jet_SrcDens;
  	  Prim[1] = Vel[0];
  	  Prim[2] = Vel[1];
  	  Prim[3] = Vel[2];
      Prim[4] = Jet_SrcTemp*Jet_SrcDens;
  	}
  	else
  	{
      // shift origin to the point D/E
	  if ( InsideUpperCone == true )   xp[2] += Jet_Radius/tan(Jet_HalfOpeningAngle);
	  else                             xp[2] -= Jet_Radius/tan(Jet_HalfOpeningAngle);

      CartesianRotate(xp, Jet_PrecessionAngle, Omega_t, true);

      // smooth velocity  
      Cartesian2Spherical(xp, rp);

      if ( InsideLowerCone == true ) rp[1] -= M_PI;

      if ( Jet_SmoothVel ) Jet_SrcVelSmooth = Jet_SrcVel*SQR(cos( 0.5 * M_PI * rp[1] / Jet_HalfOpeningAngle ));
      else                 Jet_SrcVelSmooth = Jet_SrcVel;

  
      if ( Jet_PrecessionAxis[0] != 0.0 || Jet_PrecessionAxis[1] != 0.0 ||  Jet_PrecessionAxis[2] == 0.0 )
         CartesianRotate(xp, PrecessionAxis_Spherical[1], PrecessionAxis_Spherical[2], true);
  
      Cartesian2Spherical(xp, rp);

      Prim[0] = Jet_SrcDens;
      Prim[1] = Jet_SrcVelSmooth*sin(rp[1])*cos(rp[2]);
      Prim[2] = Jet_SrcVelSmooth*sin(rp[1])*sin(rp[2]);
      Prim[3] = Jet_SrcVelSmooth*cos(rp[1]);
      Prim[4] = Jet_SrcTemp*Jet_SrcDens;
  	}

#   ifdef FLOAT8
	SRHydro_Pri2Con( Prim, fluid, GAMMA );
#   else
	SRHydro_Pri2Con( Prim, Cons, GAMMA );

	fluid[0] = (real)Cons[0];
	fluid[1] = (real)Cons[1];
	fluid[2] = (real)Cons[2];
	fluid[3] = (real)Cons[3];
	fluid[4] = (real)Cons[4];
#   endif

	return true;
  }



  return false;


} // FUNCTION : Flu_ResetByUser_Jets


// (true/false): if the target cell (is/is not) within the region to be refined
static bool Flag_Region( const int i, const int j, const int k, const int lv, const int PID )
{

   if ( Jet_HSE_Radius > 0.0 && OPT__FLAG_LOHNER_DENS == 1 )
   {
      const double dh     = amr->dh[lv];                                                         // grid size
      const double Pos[3] = { amr->patch[0][lv][PID]->EdgeL[0] + (i+0.5)*dh,  // x,y,z position
                              amr->patch[0][lv][PID]->EdgeL[1] + (j+0.5)*dh,
                              amr->patch[0][lv][PID]->EdgeL[2] + (k+0.5)*dh  };
   
      bool Flag = false;  
   
      const double Center[3]      = { 0.5*amr->BoxSize[0] + Jet_HSE_Dx, 
                                      0.5*amr->BoxSize[1] + Jet_HSE_Dy, 
                                      0.5*amr->BoxSize[2] + Jet_HSE_Dz };
   
      const double dR[3]          = { Pos[0]-Center[0]-Jet_CenOffset[0], 
                                      Pos[1]-Center[1]-Jet_CenOffset[1], 
                                      Pos[2]-Center[2]-Jet_CenOffset[2] };
   
      const double R              = sqrt( SQR(dR[0]) + SQR(dR[1]) + SQR(dR[2]) );
   
      const double ShellThickness = 16*amr->dh[3];
     
   
   
      if ( R < Jet_HSE_Radius - ShellThickness )   return true;
      else                                         return false;
   }

   return true;

} // FUNCTION : Flag_Region



bool Flag_User( const int i, const int j, const int k, const int lv, const int PID, const double Threshold )
{
   const double dh     = amr->dh[lv];                                                  // grid size
   const double Pos[3] = { amr->patch[0][lv][PID]->EdgeL[0] + (i+0.5)*dh,              // x,y,z position
                           amr->patch[0][lv][PID]->EdgeL[1] + (j+0.5)*dh,
                           amr->patch[0][lv][PID]->EdgeL[2] + (k+0.5)*dh  };

   const double Center[3]      = { Jet_Cen[0], 
                                   Jet_Cen[1], 
                                   Jet_Cen[2] };

   const double dR[3]          = { Pos[0]-Center[0], Pos[1]-Center[1], Pos[2]-Center[2] };
   const double R              = sqrt( SQR(dR[0]) + SQR(dR[1]) + SQR(dR[2]) );


  
   if ( Jet_HSE_Radius*0.9 < R && R < Jet_HSE_Radius && Jet_HSE_Radius > 0.0 && Step > 1 )
   {
      if ( lv == MAX_LEVEL-1 )
      {
        if ( MPI_Rank == 0 ) 
        #pragma omp master
        {
          system("touch STOP_GAMER_STOP");
        }
      }
     
   }





   bool Flag = false;

   Flag |= R < dh*1.8;
   Flag |= R < Jet_MaxDis;

   //Flag |= fabs(dR[0]) < 2.5 && fabs(dR[1]) < 0.3 && fabs(dR[2]) < 0.3;

   if ( Flag ) return true;
   else        return false;

} // FUNCTION : Flag_User


#  ifdef GRAVITY
void Init_ExternalAcc()
{
  ExtAcc_AuxArray[0] = 0.5*amr->BoxSize[0] + Jet_HSE_Dx;
  ExtAcc_AuxArray[1] = 0.5*amr->BoxSize[1] + Jet_HSE_Dy;
  ExtAcc_AuxArray[2] = 0.5*amr->BoxSize[2] + Jet_HSE_Dz;

  if ( Jet_Ambient == 2 )
  {
    ExtAcc_AuxArray[3] = 2.0 * Jet_UniformTemp;
    ExtAcc_AuxArray[4] = 0.0;
    ExtAcc_AuxArray[5] = 2.0;
  }
  else if ( Jet_Ambient == 3 )
  {
    ExtAcc_AuxArray[3] = 3.0 * Jet_UniformTemp * Jet_Beta_Beta;
    ExtAcc_AuxArray[4] = Jet_Beta_Rcore;
    ExtAcc_AuxArray[5] = 4.0;
  }
  
}
#endif

//-------------------------------------------------------------------------------------------------------
// Function    :  Init_TestProb_SRHydro_Jets
// Description :  Test problem initializer
//
// Note        :  None
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Init_TestProb_SRHydro_Jets()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );


// validate the compilation flags and runtime parameters
   Validate();


// set the problem-specific runtime parameters
   SetParameter();


// get enclosed mass
//   if ( Jet_Ambient == 1 )
//   GetEnclosedMass( &Jet_HSE_Mass );

   Init_Function_User_Ptr   = SetGridIC;
   Flag_User_Ptr            = Flag_User;
   Flag_Region_Ptr          = Flag_Region;
   Mis_GetTimeStep_User_Ptr = NULL;
   BC_User_Ptr              = NULL;
   Flu_ResetByUser_Func_Ptr = Flu_ResetByUser_Jets;
   Output_User_Ptr          = NULL;
   Aux_Record_User_Ptr      = NULL;
   End_User_Ptr             = NULL;
#  ifdef GRAVITY
   Init_ExternalAcc_Ptr     = Init_ExternalAcc;
#  endif


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : Init_TestProb_SRHydro_Jets

#endif
