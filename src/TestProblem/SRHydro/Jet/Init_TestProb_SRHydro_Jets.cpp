#include "GAMER.h"
#include "TestProb.h"


#if ( MODEL == SR_HYDRO )

// problem-specific global variables
// =======================================================================================

// options
static bool     Jet_Ambient;                                    // [0/1] : uniform/hydro-static equilbrium ambient
static bool     Jet_DiffPrecession;                             //
static bool     Jet_Precession;                                 //
static bool     Jet_TimeDependentSrc;                           //
static bool     Jet_Exist;                                      //

// general parameters
static double   Jet_ParticleMassSrc;                            // particle mass in jet source [g]
static double   Jet_ParticleMassAmbient;                        // particle mass in ambient    [g]

// uniform background parameters
static double   Jet_UniformDens;                                // uniform ambient density
static double   Jet_UniformVel[3];                              // uniform ambient 4-velocity
static double   Jet_UniformTemp;                                // uniform ambient temperature

// halo parameters
static double   Jet_HSE_Dx                         = NULL_REAL; // for Jet_Ambient=1: distance between box and cluster centre (along x)
static double   Jet_HSE_Dy                         = NULL_REAL; // for Jet_Ambient=1: distance between box and cluster centre (along y)
static double   Jet_HSE_Dz                         = NULL_REAL; // for Jet_Ambient=1: distance between box and cluster centre (along z)
static char     Jet_HSE_BgTable_File[MAX_STRING];               // for Jet_Ambient=1: filename of the background gas table
static double  *Jet_HSE_BgTable_Data               = NULL;      // for Jet_Ambient=1: background gas table [radius/density/temperature]
static int      Jet_HSE_BgTable_NBin;                           // for Jet_Ambient=1: number of bins in Jet_HSE_BgTable_Data[]
static double   Jet_HSE_Radius                     = NULL_REAL; // for Jet_Ambient=1: radius of halo
static double   Jet_HSE_Mass                       = NULL_REAL; // for Jet_Ambient=1: enclosed mass of halo

// jet parameters
static double   Jet_Radius;                                     // radius of the cylinder-shape jet source
static double   Jet_HalfHeight;                                 // half height of the cylinder-shape jet source
static double   Jet_SrcVel;                                     // jet 4-velocity
static double   Jet_SrcDens;                                    // jet density
static double   Jet_Vec[3];                                     // jet orientation vector (x,y,z) (NOT necessary to be a unit vector)
static double   Jet_CenOffset[3];                               // jet central coordinates offset
static double   Jet_SrcTemp;                                    // jet temperature
static double   Jet_Cen[3];                                     // jet central coordinates
static double   Jet_WaveK;                                      // jet wavenumber used in the sin() function to have smooth bidirectional jets
static double   Jet_MaxDis;                                     // maximum distance between the cylinder-shape jet source and the jet center
                                                                // --> ( Jet_Radius^2 + Jet_HalfHeight^2 )^0.5
// precession parameters
static double   Jet_Angular_Velocity;                           // precession angular velocity (degree per code_time)
static double   Jet_Angle;                                      // precession angle in degree
static    int   Jet_NumShell;                                   // number of shells
static double   BH_Radius;                                      // radius of black hole
static double   Jet_Cone_Vec[3];                                // cone orientation vector (x,y,z) (NOT necessary to be a unit vector)

// time-depent source
static double   Jet_BurstStartTime;                             // start burst time in jet
static double   Jet_BurstEndTime;                               // end burst time in jet
static double   Jet_Burst4Vel;                                  // burst 4-velocity
static double   Jet_BurstDens;                                  // burst proper density
static double   Jet_BurstTemp;                                  // burst temperature
static bool     Flag_Burst4Vel;
static bool     Flag_BurstDens;
static bool     Flag_BurstTemp;


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
   ReadPara->Add( "Jet_Ambient",             &Jet_Ambient,            false,          Useless_bool,   Useless_bool    );
   ReadPara->Add( "Jet_DiffPrecession",      &Jet_DiffPrecession,     false,          Useless_bool,   Useless_bool    );
   ReadPara->Add( "Jet_Precession",          &Jet_Precession,         false,          Useless_bool,   Useless_bool    );
   ReadPara->Add( "Jet_TimeDependentSrc",    &Jet_TimeDependentSrc,   false,          Useless_bool,   Useless_bool    );
   ReadPara->Add( "Jet_Exist",               &Jet_Exist,              false,          Useless_bool,   Useless_bool    );

// load general parameters
   ReadPara->Add( "Jet_ParticleMassSrc",     &Jet_ParticleMassSrc,     -1.0,          Eps_double,     NoMax_double    );
   ReadPara->Add( "Jet_ParticleMassAmbient", &Jet_ParticleMassAmbient, -1.0,          Eps_double,     NoMax_double    );
 

// load hydrostatic equilibrium parameters
   ReadPara->Add( "Jet_HSE_Radius",          &Jet_HSE_Radius,          -1.0,          Eps_double,     NoMax_double    );
   ReadPara->Add( "Jet_HSE_Dx",              &Jet_HSE_Dx,               0.0,          NoMin_double,   NoMax_double    );
   ReadPara->Add( "Jet_HSE_Dy",              &Jet_HSE_Dy,               0.0,          NoMin_double,   NoMax_double    );
   ReadPara->Add( "Jet_HSE_Dz",              &Jet_HSE_Dz,               0.0,          NoMin_double,   NoMax_double    );
   ReadPara->Add( "Jet_HSE_BgTable_File",     Jet_HSE_BgTable_File,     Useless_str,  Useless_str,    Useless_str     );

// load background parameters
   ReadPara->Add( "Jet_UniformDens",         &Jet_UniformDens,         -1.0,          Eps_double,     NoMax_double    );
   ReadPara->Add( "Jet_UniformVel_x",        &Jet_UniformVel[0],        0.0,          NoMin_double,   NoMax_double    );
   ReadPara->Add( "Jet_UniformVel_y",        &Jet_UniformVel[1],        0.0,          NoMin_double,   NoMax_double    );
   ReadPara->Add( "Jet_UniformVel_z",        &Jet_UniformVel[2],        0.0,          NoMin_double,   NoMax_double    );
   ReadPara->Add( "Jet_UniformTemp",         &Jet_UniformTemp,         -1.0,          Eps_double,     NoMax_double    );

// load jet parameters
   ReadPara->Add( "Jet_Radius",              &Jet_Radius    ,          -1.0,          Eps_double,     NoMax_double    );
   ReadPara->Add( "Jet_HalfHeight",          &Jet_HalfHeight,          -1.0,          Eps_double,     NoMax_double    );
   ReadPara->Add( "Jet_SrcVel",              &Jet_SrcVel    ,          -1.0,          NoMin_double,   NoMax_double    );
   ReadPara->Add( "Jet_SrcDens",             &Jet_SrcDens   ,          -1.0,          Eps_double,     NoMax_double    );
   ReadPara->Add( "Jet_SrcTemp",             &Jet_SrcTemp   ,          -1.0,          Eps_double,     NoMax_double    );
   ReadPara->Add( "Jet_Vec_x",               &Jet_Vec       [0],        NoDef_double, NoMin_double,   NoMax_double    );
   ReadPara->Add( "Jet_Vec_y",               &Jet_Vec       [1],        NoDef_double, NoMin_double,   NoMax_double    );
   ReadPara->Add( "Jet_Vec_z",               &Jet_Vec       [2],        NoDef_double, NoMin_double,   NoMax_double    );
   ReadPara->Add( "Jet_CenOffset_x",         &Jet_CenOffset [0],        NoDef_double, NoMin_double,   NoMax_double    );
   ReadPara->Add( "Jet_CenOffset_y",         &Jet_CenOffset [1],        NoDef_double, NoMin_double,   NoMax_double    );
   ReadPara->Add( "Jet_CenOffset_z",         &Jet_CenOffset [2],        NoDef_double, NoMin_double,   NoMax_double    );

// load precission parameters
   ReadPara->Add( "Jet_Angular_Velocity",    &Jet_Angular_Velocity,     NoDef_double, NoMin_double,   NoMax_double    );
   ReadPara->Add( "Jet_Angle",               &Jet_Angle,                NoDef_double, NoMin_double,   NoMax_double    );
   ReadPara->Add( "Jet_NumShell",            &Jet_NumShell,             1           , NoMin_int,      NoMax_int       );
   ReadPara->Add( "BH_Radius",               &BH_Radius,                NoDef_double, Eps_double,     NoMax_double    );
   ReadPara->Add( "Jet_Cone_Vec_x",          &Jet_Cone_Vec[0],          NoDef_double, NoMin_double,   NoMax_double    );               
   ReadPara->Add( "Jet_Cone_Vec_y",          &Jet_Cone_Vec[1],          NoDef_double, NoMin_double,   NoMax_double    );  
   ReadPara->Add( "Jet_Cone_Vec_z",          &Jet_Cone_Vec[2],          NoDef_double, NoMin_double,   NoMax_double    );  


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

// check hydro-static equilibrium
   if ( Jet_Ambient )
   {
      if ( !OPT__FLAG_REGION ) 
        Aux_Error( ERROR_INFO, "OPT__FLAG_REGION must be enabled !!\n" );

      if ( Jet_HSE_Radius == NULL_REAL ) 
        Aux_Error( ERROR_INFO, "Jet_HSE_Radius = %e !!\n" );

      if ( Jet_HSE_Radius >= 0.5*amr->BoxSize[0] || 
           Jet_HSE_Radius >= 0.5*amr->BoxSize[1] || 
           Jet_HSE_Radius >= 0.5*amr->BoxSize[2]   )
        Aux_Error( ERROR_INFO, "halo size is greater or equal to box size !!\n" );

#     ifndef GRAVITY
      Aux_Error( ERROR_INFO, "GRAVITY must be enabled !!\n" );
#     endif
   }
// check uniform ambient
   else
   {
#     ifdef GRAVITY
      Aux_Error( ERROR_INFO, "GRAVITY must be disabled !!\n" );
#     endif
   }

// check precession parameters
   if ( Jet_Precession && Jet_DiffPrecession )
      Aux_Error( ERROR_INFO, "Jet_Precession and Jet_DiffPrecession are not allowed to enable simultaneously !!\n" );


// replace useless parameters with NaN
   if ( Jet_Ambient )
   {
     Jet_UniformDens      = NAN;
     Jet_UniformTemp      = NAN;
     Jet_UniformVel[0]    = NAN;
     Jet_UniformVel[1]    = NAN;
     Jet_UniformVel[2]    = NAN;
   }
   else
   {
     Jet_HSE_Radius       = NAN;
     Jet_HSE_Dx           = NAN;
     Jet_HSE_Dy           = NAN;
     Jet_HSE_Dz           = NAN;
   }

   if ( Jet_Precession )       // uniform precession
   {
     Jet_NumShell = 1;
     BH_Radius    = NAN;
   }
   else if ( Jet_DiffPrecession )  // differential precession
   {
   }
   else if ( !Jet_Precession && !Jet_DiffPrecession ) // no precession
   {
     Jet_NumShell = 1;
     BH_Radius    = NAN;
     Jet_Angular_Velocity = 0.0;
     Jet_Angle   = 0.0;
     Jet_Cone_Vec[0] = Jet_Vec[0];
     Jet_Cone_Vec[1] = Jet_Vec[1];
     Jet_Cone_Vec[2] = Jet_Vec[2];
   }
   
   if ( !Jet_TimeDependentSrc )
   {
     Jet_BurstDens       = NAN;
     Jet_Burst4Vel       = NAN;
     Jet_BurstTemp       = NAN; 
     Jet_BurstStartTime  = NAN;
     Jet_BurstEndTime    = NAN;
   }

// (1-2) convert to code unit
   Jet_Radius               *= Const_kpc / UNIT_L;
   Jet_HalfHeight           *= Const_kpc / UNIT_L;
   Jet_SrcVel               *= Const_c   / UNIT_V;
   Jet_SrcTemp              *= Const_keV / ( Jet_ParticleMassSrc * SQR(Const_c) );
   Jet_SrcDens              *= 1.0       / UNIT_D;

   Jet_CenOffset[0]         *= Const_kpc / UNIT_L;
   Jet_CenOffset[1]         *= Const_kpc / UNIT_L;
   Jet_CenOffset[2]         *= Const_kpc / UNIT_L;

   if ( Jet_Ambient )
   {
     Jet_HSE_Radius         *= Const_kpc / UNIT_L; 
     Jet_HSE_Dx             *= Const_kpc / UNIT_L;  
     Jet_HSE_Dy             *= Const_kpc / UNIT_L;  
     Jet_HSE_Dz             *= Const_kpc / UNIT_L;  
   }
   else
   {
     Jet_UniformDens        *= 1.0       / UNIT_D;
     Jet_UniformTemp        *= Const_keV / ( Jet_ParticleMassAmbient * SQR(Const_c) );
     Jet_UniformVel[0]      *= Const_c   / UNIT_V;
     Jet_UniformVel[1]      *= Const_c   / UNIT_V;
     Jet_UniformVel[2]      *= Const_c   / UNIT_V;
   }

   Jet_Angular_Velocity   *= 1e3*Const_yr / UNIT_T;  // the unit of Jet_Angular_Velocity is 1/kyr

   if ( Jet_Precession ) // uniform precession
   {
   }

   if ( Jet_DiffPrecession ) // differential precession
   {
     BH_Radius              *= Const_kpc / UNIT_L;
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
     Jet_BurstTemp          *= Const_keV    / ( Const_mp * SQR(Const_c) );
   }

// (2) set the problem-specific derived parameters
   Jet_WaveK   = 0.5*M_PI/Jet_HalfHeight;
   Jet_MaxDis  = sqrt( SQR(Jet_HalfHeight) + SQR(Jet_Radius) );

   for (int d=0; d<3; d++)    Jet_Cen[d] = 0.5*amr->BoxSize[d] + Jet_CenOffset[d];


// (3) load the hydrostatic equilibrium table
   if ( Jet_Ambient  &&  OPT__INIT != INIT_BY_RESTART )
   {
      const bool RowMajor_No  = false;       // load data into the column-major order
      const bool AllocMem_Yes = true;        // allocate memory for Merger_Prof1/2
      const int  NCol         = 3;           // total number of columns to load
      const int  Col[NCol]    = {1, 2, 3};   // target columns: (radius, density, temperature)
      double *Table_R, *Table_D, *Table_T;

      Jet_HSE_BgTable_NBin = Aux_LoadTable( Jet_HSE_BgTable_Data, Jet_HSE_BgTable_File, NCol, Col, RowMajor_No, AllocMem_Yes );


//    convert to code units
      Table_R = Jet_HSE_BgTable_Data + 0*Jet_HSE_BgTable_NBin;
      Table_D = Jet_HSE_BgTable_Data + 1*Jet_HSE_BgTable_NBin;
      Table_T = Jet_HSE_BgTable_Data + 2*Jet_HSE_BgTable_NBin;

      for (int b=0; b<Jet_HSE_BgTable_NBin; b++)
      {
         Table_R[b] *= Const_kpc / UNIT_L;
         Table_D[b] *= 1.0       / UNIT_D;
         Table_T[b] *= Const_keV / ( Const_mp * SQR(Const_c) );
      }
   } // if ( Jet_Ambient  &&  OPT__INIT != INIT_BY_RESTART )


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
     Aux_Message( stdout, "  Jet_Ambient             = %d\n",                Jet_Ambient                                     );
     Aux_Message( stdout, "  Jet_DiffPrecession      = %d\n",                Jet_DiffPrecession                              );
     Aux_Message( stdout, "  Jet_Precession          = %d\n",                Jet_Precession                                  );
     Aux_Message( stdout, "  Jet_TimeDependentSrc    = %d\n",                Jet_TimeDependentSrc                            );
     Aux_Message( stdout, "  Jet_ParticleMassSrc     = %14.7e g\n",          Jet_ParticleMassSrc                             );
     Aux_Message( stdout, "  Jet_ParticleMassAmbient = %14.7e g\n",          Jet_ParticleMassAmbient                         );
     Aux_Message( stdout, "  Jet_Radius              = %14.7e kpc\n",        Jet_Radius*UNIT_L/Const_kpc                     );
     Aux_Message( stdout, "  Jet_HalfHeight          = %14.7e kpc\n",        Jet_HalfHeight*UNIT_L/Const_kpc                 );
     Aux_Message( stdout, "  Jet_SrcVel              = %14.7e c\n",          Jet_SrcVel                                      );
     Aux_Message( stdout, "  Jet_SrcDens             = %14.7e g/cm^3\n",     Jet_SrcDens*UNIT_D                              );
     Aux_Message( stdout, "  Jet_SrcTemp             = %14.7e keV\n",        Jet_SrcTemp*Const_mp*SQR(Const_c)/Const_keV     );
     Aux_Message( stdout, "  Jet_NumDensSrc          = %14.7e per cc\n",     Jet_SrcDens*UNIT_D/Jet_ParticleMassSrc          );
     Aux_Message( stdout, "  Jet_CenOffset[x]        = %14.7e kpc\n",        Jet_CenOffset [0]*UNIT_L/Const_kpc              );
     Aux_Message( stdout, "  Jet_CenOffset[y]        = %14.7e kpc\n",        Jet_CenOffset [1]*UNIT_L/Const_kpc              );
     Aux_Message( stdout, "  Jet_CenOffset[z]        = %14.7e kpc\n",        Jet_CenOffset [2]*UNIT_L/Const_kpc              );
     Aux_Message( stdout, "  Jet_Angular_Velocity    = %14.7e degree/kyr\n", Jet_Angular_Velocity* UNIT_T / (1e3*Const_yr)   );
     Aux_Message( stdout, "  Jet_Angle               = %14.7e\n",            Jet_Angle                                       );
     Aux_Message( stdout, "  Jet_NumShell            = %d\n",                Jet_NumShell                                    );
   }

   if ( Jet_Ambient && MPI_Rank == 0 ) 
   {
     Aux_Message( stdout, "  Jet_HSE_Radius          = %e kpc\n",            Jet_HSE_Radius*UNIT_L/Const_kpc                 );
     Aux_Message( stdout, "  Jet_HSE_Dx              = %e kpc\n",            Jet_HSE_Dx*UNIT_L/Const_kpc                     );
     Aux_Message( stdout, "  Jet_HSE_Dy              = %e kpc\n",            Jet_HSE_Dy*UNIT_L/Const_kpc                     );
     Aux_Message( stdout, "  Jet_HSE_Dz              = %e kpc\n",            Jet_HSE_Dz*UNIT_L/Const_kpc                     );
     Aux_Message( stdout, "  Jet_HSE_BgTable_File    = %s\n",                Jet_HSE_BgTable_File                            );
   }
   else if ( !Jet_Ambient && MPI_Rank == 0 )
   {
     Aux_Message( stdout, "  Jet_UniformDens         = %14.7e g/cm^3\n",     Jet_UniformDens*UNIT_D                          );
     Aux_Message( stdout, "  Jet_UniformTemp         = %14.7e keV\n",        Jet_UniformTemp*Const_mp*SQR(Const_c)/Const_keV );
     Aux_Message( stdout, "  Jet_UniformVel[x]       = %14.7e c\n",          Jet_UniformVel[0]                               );
     Aux_Message( stdout, "  Jet_UniformVel[y]       = %14.7e c\n",          Jet_UniformVel[1]                               );
     Aux_Message( stdout, "  Jet_UniformVel[z]       = %14.7e c\n",          Jet_UniformVel[2]                               );
     Aux_Message( stdout, "  Jet_UniformNumDens      = %14.7e per cc\n",     Jet_UniformDens*UNIT_D/Jet_ParticleMassAmbient  );
   }

   if ( ( Jet_Precession || Jet_DiffPrecession ) && MPI_Rank == 0 )
   {
     Aux_Message( stdout, "  Jet_Cone_Vec[x]         = %14.7e\n",            Jet_Cone_Vec[0]                                 );
     Aux_Message( stdout, "  Jet_Cone_Vec[y]         = %14.7e\n",            Jet_Cone_Vec[1]                                 );
     Aux_Message( stdout, "  Jet_Cone_Vec[z]         = %14.7e\n",            Jet_Cone_Vec[2]                                 );
   }

   if ( Jet_DiffPrecession && MPI_Rank == 0 )
   {
     Aux_Message( stdout, "  BH_Radius               = %14.7e kpc\n",        BH_Radius*UNIT_L/Const_kpc                      );
   }

   if ( !Jet_Precession && !Jet_DiffPrecession && MPI_Rank == 0 )
   {
     Aux_Message( stdout, "  Jet_Vec[x]              = %14.7e\n",            Jet_Vec[0]                                      );
     Aux_Message( stdout, "  Jet_Vec[y]              = %14.7e\n",            Jet_Vec[1]                                      );
     Aux_Message( stdout, "  Jet_Vec[z]              = %14.7e\n",            Jet_Vec[2]                                      );
   }

   if ( Jet_TimeDependentSrc && MPI_Rank == 0 )
   {
     Aux_Message( stdout, "  Jet_BurstStartTime      = %14.7e kyr \n",       Jet_BurstStartTime*UNIT_T/(1e3*Const_yr)        );
     Aux_Message( stdout, "  Jet_BurstEndTime        = %14.7e kyr \n",       Jet_BurstEndTime*UNIT_T/(1e3*Const_yr)          );
     Aux_Message( stdout, "  Jet_Burst4Vel           = %14.7e c \n",         Jet_Burst4Vel                                   );
     Aux_Message( stdout, "  Jet_BurstDens           = %14.7e g/cm^3\n",     Jet_BurstDens*UNIT_D                            );
     Aux_Message( stdout, "  Jet_BurstTemp           = %14.7e keV\n",        Jet_BurstTemp*Const_mp*SQR(Const_c)/Const_keV   );
     Aux_Message( stdout, "  Flag_Burst4Vel          = %d\n",                Flag_Burst4Vel                                  );
     Aux_Message( stdout, "  Flag_BurstDens          = %d\n",                Flag_BurstDens                                  );
     Aux_Message( stdout, "  Flag_BurstTemp          = %d\n",                Flag_BurstTemp                                  );
   }

   if ( MPI_Rank == 0 )
   {
     Aux_Message( stdout, "=============================================================================\n" );
   }


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Setting runtime parameters ... done\n" );

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
   const double *Table_R=NULL, *Table_D=NULL, *Table_T=NULL;
   double dx, dy, dz, r;
   const int  NCol         = 3;           // total number of columns to load
   const int  Col[NCol]    = {1, 2, 3};   // target columns: (radius, density, temperature)

   double DropRatio = 1e-2;

// variables for jet
   double Pri4Vel[NCOMP_FLUID];

   Pri4Vel[1] = Jet_UniformVel[0];
   Pri4Vel[2] = Jet_UniformVel[1];
   Pri4Vel[3] = Jet_UniformVel[2];

   if ( Jet_Ambient )
   {
      double Temp, RhoSurface;

      Table_R = Jet_HSE_BgTable_Data + (Col[0]-1)*Jet_HSE_BgTable_NBin;
      Table_D = Jet_HSE_BgTable_Data + (Col[1]-1)*Jet_HSE_BgTable_NBin;
      Table_T = Jet_HSE_BgTable_Data + (Col[2]-1)*Jet_HSE_BgTable_NBin;

      dx = x - amr->BoxCenter[0] + Jet_HSE_Dx;
      dy = y - amr->BoxCenter[1] + Jet_HSE_Dy;
      dz = z - amr->BoxCenter[2] + Jet_HSE_Dz;

      r = sqrt( dx*dx + dy*dy + dz*dz );

      Pri4Vel[0] = Mis_InterpolateFromTable( Jet_HSE_BgTable_NBin, Table_R, Table_D, r );
      Pri4Vel[4] = Mis_InterpolateFromTable( Jet_HSE_BgTable_NBin, Table_R, Table_T, r ) * Pri4Vel[0];

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

        // exponential ambient outside halo
        //RhoSurface = DropRatio * Mis_InterpolateFromTable( Jet_HSE_BgTable_NBin, Table_R, Table_D, Jet_HSE_Radius );

		//Temp /= DropRatio;

        //Pri4Vel[0] = RhoSurface * exp( ( NEWTON_G * Jet_HSE_Mass / Temp ) * ( 1.0/r - 1.0/Jet_HSE_Radius ) );
        //Pri4Vel[4] = Pri4Vel[0] * Temp;
      }
   }
   else
   {
	  Pri4Vel[0] = Jet_UniformDens;
	  Pri4Vel[4] = Jet_UniformTemp * Jet_UniformDens;
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
//-------------------------------------------------------------------------------------------------------
bool Flu_ResetByUser_Jets( real fluid[], const double x, const double y, const double z, const double Time,
                                    const int lv, double AuxArray[] )
{
   if ( !Jet_Exist ) return false;

   const double r[3] = { x, y, z };

   double Jet_dr, Jet_dh, S, Area;
   double Dis_c2m, Dis_c2v, Dis_v2m, Vec_c2m[3], Vec_v2m[3];
   double Jet_SrcVel_xyz[3],Jet_BurstVel_xyz[3], Dis_Cone_Vec, Dis_Cone_Vec_2;
   double Cos_Phi, Sin_Phi, Sin_Angle, Cos_Angle, Sin_Omega, Cos_Omega;
   double Omega_t, ShellOmega, ShellRadius, ShellThick, Pri4Vel[NCOMP_FLUID];
   real   MomSin;
   bool BurstDens_Yes, Burst4Vel_Yes, BurstTemp_Yes;
   int state = 0;

// distance: jet center to mesh
   for (int d=0; d<3; d++)    Vec_c2m[d] = r[d] - Jet_Cen[d];
   Dis_c2m = sqrt( SQR(Vec_c2m[0]) + SQR(Vec_c2m[1]) + SQR(Vec_c2m[2]) );

// exclude cells far away from the jet source
   if ( Dis_c2m > Jet_MaxDis )   return false;

   Dis_Cone_Vec_2 = sqrt(SQR(Jet_Cone_Vec[0]) + SQR(Jet_Cone_Vec[1])); 
   Dis_Cone_Vec   = sqrt(SQR(Jet_Cone_Vec[0]) + SQR(Jet_Cone_Vec[1]) + SQR(Jet_Cone_Vec[2]));

   ShellThick = Jet_Radius / (double) Jet_NumShell;

   Sin_Angle = sin(Jet_Angle*M_PI / 180.0);
   Cos_Angle = cos(Jet_Angle*M_PI / 180.0);

       
   for ( int i = Jet_NumShell ; i >= 1; i-- ) 
   {
       ShellRadius = i * ShellThick;
   
       if ( Jet_DiffPrecession )
         ShellOmega = Jet_Angular_Velocity * pow( 1.0 + SQR( ( ShellRadius - 0.5*ShellThick ) / BH_Radius ), -1.50 );
       else
         ShellOmega = Jet_Angular_Velocity;

   
       Omega_t = ShellOmega * Time * M_PI / 180.0;
   
       Sin_Omega = sin(Omega_t);
       Cos_Omega = cos(Omega_t);

//     coordinate transformation
       if ( Dis_Cone_Vec_2 >  __DBL_EPSILON__ )
       {
          Cos_Phi = Jet_Cone_Vec[0] / Dis_Cone_Vec_2;
          Sin_Phi = Jet_Cone_Vec[1] / Dis_Cone_Vec_2;

   
          Jet_Vec[0] = - Dis_Cone_Vec   * Sin_Phi   * Sin_Angle * Cos_Omega
                       - Jet_Cone_Vec[2]* Cos_Phi   * Sin_Angle * Sin_Omega
                       + Dis_Cone_Vec_2 * Cos_Phi   * Cos_Angle;
   
          Jet_Vec[1] =   Dis_Cone_Vec   * Cos_Phi   * Sin_Angle * Cos_Omega
                       - Jet_Cone_Vec[2]* Sin_Phi   * Sin_Angle * Sin_Omega
                       + Dis_Cone_Vec_2 * Sin_Phi   * Cos_Angle;
   
          Jet_Vec[2] =   Dis_Cone_Vec_2 * Sin_Angle * Sin_Omega 
                       + Jet_Cone_Vec[2]* Cos_Angle;
       }
       else
       {
          Jet_Vec[0] = Dis_Cone_Vec * Sin_Angle * Cos_Omega;
          Jet_Vec[1] = Dis_Cone_Vec * Sin_Angle * Sin_Omega;
          Jet_Vec[2] = Dis_Cone_Vec * Cos_Angle;
       }


//     distance: jet center to temporary vector
       Dis_c2v = sqrt( SQR(Jet_Vec[0]) + SQR(Jet_Vec[1]) + SQR(Jet_Vec[2]) );

//     velocity components along different directions
       for (int d=0; d<3; d++)
       {
         Jet_SrcVel_xyz[d]   = Jet_SrcVel    * Jet_Vec[d] / Dis_c2v;
         Jet_BurstVel_xyz[d] = Jet_Burst4Vel * Jet_Vec[d] / Dis_c2v;
       }

//     distance: temporary vector to mesh
       for (int d=0; d<3; d++)    Vec_v2m[d] = r[d] - Jet_Cen[d] - Jet_Vec[d];
       Dis_v2m = sqrt( SQR(Vec_v2m[0]) + SQR(Vec_v2m[1]) + SQR(Vec_v2m[2]) );


//     check whether or not the target cell is within the jet source
       S      = 0.5*( Dis_c2m + Dis_v2m + Dis_c2v );
       Area   = sqrt( S*(S-Dis_c2m)*(S-Dis_v2m)*(S-Dis_c2v) );

       Jet_dr = 2.0*Area/Dis_c2v;
       Jet_dh = sqrt( Dis_c2m*Dis_c2m - Jet_dr*Jet_dr );

//     reset the fluid variables within the jet source
       if ( Jet_dh <= Jet_HalfHeight  &&  Jet_dr <= ShellRadius )
       {
//        use a sine function to make the velocity smooth within the jet from +Jet_SrcVel to -Jet_SrcVel
          MomSin      = sin( Jet_WaveK*Jet_dh );
          MomSin     *= SIGN( Vec_c2m[0]*Jet_Vec[0] + Vec_c2m[1]*Jet_Vec[1] + Vec_c2m[2]*Jet_Vec[2] );
         
          BurstDens_Yes = ( Flag_BurstDens && Jet_BurstStartTime < Time && Time < Jet_BurstEndTime ) ? true : false;
          Burst4Vel_Yes = ( Flag_Burst4Vel && Jet_BurstStartTime < Time && Time < Jet_BurstEndTime ) ? true : false;
          BurstTemp_Yes = ( Flag_BurstTemp && Jet_BurstStartTime < Time && Time < Jet_BurstEndTime ) ? true : false;

          Pri4Vel[0] = BurstDens_Yes ? Jet_BurstDens               : Jet_SrcDens;
          Pri4Vel[1] = Burst4Vel_Yes ? Jet_BurstVel_xyz[0]*MomSin  : Jet_SrcVel_xyz[0]*MomSin;
          Pri4Vel[2] = Burst4Vel_Yes ? Jet_BurstVel_xyz[1]*MomSin  : Jet_SrcVel_xyz[1]*MomSin;
          Pri4Vel[3] = Burst4Vel_Yes ? Jet_BurstVel_xyz[2]*MomSin  : Jet_SrcVel_xyz[2]*MomSin;
          Pri4Vel[4] = BurstTemp_Yes ? Jet_BurstTemp*Jet_BurstDens : Jet_SrcTemp*Jet_SrcDens;

//        cast double to real
#         ifndef FLOAT8
          double Out[NCOMP_FLUID];

          SRHydro_Pri2Con(Pri4Vel, Out,  GAMMA);

          fluid [0] = (real) Out[0]*pow(10,i);
          fluid [1] = (real) Out[1];
          fluid [2] = (real) Out[2];
          fluid [3] = (real) Out[3];
          fluid [4] = (real) Out[4];
#         else
          SRHydro_Pri2Con(Pri4Vel, fluid, GAMMA);
#         endif
   
          state++;
       } // if (  Jet_dh <= Jet_HalfHeight  &&  Jet_dr <= Jet_Radius )

   }

   if ( state > 0 ) return true;
   else             return false;





} // FUNCTION : Flu_ResetByUser_Jets


// (true/false): if the target cell (is/is not) within the region to be refined
bool Flag_Region( const int i, const int j, const int k, const int lv, const int PID )
{
   if ( !Jet_Ambient ) return true;


   const double dh     = amr->dh[lv];                                                  // grid size
   const double Pos[3] = { amr->patch[0][lv][PID]->EdgeL[0] + (i+0.5)*dh,              // x,y,z position
                           amr->patch[0][lv][PID]->EdgeL[1] + (j+0.5)*dh,
                           amr->patch[0][lv][PID]->EdgeL[2] + (k+0.5)*dh  };

   bool Flag = false;  

   const double Center[3]      = { 0.5*amr->BoxSize[0], 0.5*amr->BoxSize[1], 0.5*amr->BoxSize[2] };
   const double dR[3]          = { Pos[0]-Center[0], Pos[1]-Center[1], Pos[2]-Center[2] };
   const double R              = sqrt( SQR(dR[0]) + SQR(dR[1]) + SQR(dR[2]) );

   const double ShellThickness = 16*amr->dh[0];
  


   if ( R < Jet_HSE_Radius - ShellThickness )   return true;
//   else if ( Jet_HSE_Radius - ShellThickness <= R && R <= Jet_HSE_Radius + ShellThickness && lv <= 0 )
//   {
//      Flag = true;
//   }
   else                                         return false;

} // FUNCTION : Flag_Region



bool Flag_User( const int i, const int j, const int k, const int lv, const int PID, const double Threshold )
{
   const double dh     = amr->dh[lv];                                                  // grid size
   const double Pos[3] = { amr->patch[0][lv][PID]->EdgeL[0] + (i+0.5)*dh,              // x,y,z position
                           amr->patch[0][lv][PID]->EdgeL[1] + (j+0.5)*dh,
                           amr->patch[0][lv][PID]->EdgeL[2] + (k+0.5)*dh  };

//   const real (*Dens)[PS1][PS1] = amr->patch[ amr->FluSg[lv] ][lv][PID]->fluid[DENS];  // density
//   const real (*MomX)[PS1][PS1] = amr->patch[ amr->FluSg[lv] ][lv][PID]->fluid[MOMX];  // momentum x
//   const real (*MomY)[PS1][PS1] = amr->patch[ amr->FluSg[lv] ][lv][PID]->fluid[MOMY];  // momentum y
//   const real (*MomZ)[PS1][PS1] = amr->patch[ amr->FluSg[lv] ][lv][PID]->fluid[MOMZ];  // momentum z
//   const real (*Engy)[PS1][PS1] = amr->patch[ amr->FluSg[lv] ][lv][PID]->fluid[ENGY];  // total energy
//#  ifdef GRAVITY
//   const real (*Pot )[PS1][PS1] = amr->patch[ amr->PotSg[lv] ][lv][PID]->pot;          // potential
//#  endif


   const double Center[3]      = { Jet_Cen[0], Jet_Cen[1], Jet_Cen[2] };
   const double dR[3]          = { Pos[0]-Center[0], Pos[1]-Center[1], Pos[2]-Center[2] };
   const double R              = sqrt( SQR(dR[0]) + SQR(dR[1]) + SQR(dR[2]) );


   if ( R <= Jet_MaxDis )   return true;
   else                     return false;

} // FUNCTION : Flag_User



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
//   if ( Jet_Ambient )
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
   Init_ExternalAcc_Ptr     = NULL;
#  endif


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : Init_TestProb_SRHydro_Jets

#endif
