#include "GAMER.h"
#include "TestProb.h"

#if  ( MODEL == SR_HYDRO )

// problem-specific global variables
// =======================================================================================
static double Pulsar_DensBg;           // background mass density
static double Pulsar_PresBg;           // background pressure

static double Pulsar_WindDensRatio;     // density ratio of center to background
static double Pulsar_WindVelocity;      //
static double Pulsar_WindPresRatio;     // density ratio of center to background

static double Pulsar_JetDensRatio;      // pressure ratio of center to background
static double Pulsar_JetVelocity;       //
static double Pulsar_JetPresRatio;      // pressure ratio of center to background

static double Pulsar_DiskRadius;        // explosion radius
static double Pulsar_DiskHeight;        // explosion radius
static double Pulsar_Center[3];     // explosion center
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


#  if ( MODEL != SR_HYDRO )
   Aux_Error( ERROR_INFO, "MODEL != SR_HYDRO !!\n" );
#  endif

#  ifdef GRAVITY
   Aux_Error( ERROR_INFO, "GRAVITY must be disabled !!\n" );
#  endif

#  ifdef COMOVING
   Aux_Error( ERROR_INFO, "COMOVING must be disabled !!\n" );
#  endif

#  ifdef PARTICLE
   Aux_Error( ERROR_INFO, "PARTICLE must be disabled !!\n" );
#  endif

   if ( !OPT__INIT_RESTRICT )
      Aux_Error( ERROR_INFO, "OPT__INIT_RESTRICT must be enabled !!\n" );


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Validating test problem %d ... done\n", TESTPROB_ID );

} // FUNCTION : Validate



#if ( MODEL == SR_HYDRO )
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

// add parameters in the following format:
// --> note that VARIABLE, DEFAULT, MIN, and MAX must have the same data type
// --> some handy constants (e.g., NoMin_int, Eps_float, ...) are defined in "include/ReadPara.h"
// ********************************************************************************************************************************
// ReadPara->Add( "KEY_IN_THE_FILE",       &VARIABLE_ADDRESS,     DEFAULT,         MIN,                MAX    );
// ********************************************************************************************************************************
   ReadPara->Add( "Pulsar_DensBg",         &Pulsar_DensBg,           -1.0,  Eps_double,       NoMax_double    );
   ReadPara->Add( "Pulsar_PresBg",         &Pulsar_PresBg,           -1.0,  Eps_double,       NoMax_double    );
   ReadPara->Add( "Pulsar_WindDensRatio",  &Pulsar_WindDensRatio,    -1.0,  Eps_double,       NoMax_double    );
   ReadPara->Add( "Pulsar_WindVelocity",   &Pulsar_WindVelocity,     -1.0,  NoMin_double,     NoMax_double    );
   ReadPara->Add( "Pulsar_WindPresRatio",  &Pulsar_WindPresRatio,    -1.0,  Eps_double,       NoMax_double    );
   ReadPara->Add( "Pulsar_JetDensRatio",   &Pulsar_JetDensRatio,     -1.0,  Eps_double,       NoMax_double    );
   ReadPara->Add( "Pulsar_JetVelocity",    &Pulsar_JetVelocity,      -1.0,  NoMin_double,     NoMax_double    );
   ReadPara->Add( "Pulsar_JetPresRatio",   &Pulsar_JetPresRatio,     -1.0,  Eps_double,       NoMax_double    );
   ReadPara->Add( "Pulsar_DiskRadius",     &Pulsar_DiskRadius,       -1.0,  Eps_double,       NoMax_double    );
   ReadPara->Add( "Pulsar_DiskHeight",     &Pulsar_DiskHeight,       -1.0,  Eps_double,       NoMax_double    );
   ReadPara->Add( "Pulsar_Center_X",       &Pulsar_Center[0],        -1.0,  NoMin_double,     amr->BoxSize[0] );
   ReadPara->Add( "Pulsar_Center_Y",       &Pulsar_Center[1],        -1.0,  NoMin_double,     amr->BoxSize[1] );
   ReadPara->Add( "Pulsar_Center_Z",       &Pulsar_Center[2],        -1.0,  NoMin_double,     amr->BoxSize[2] );

   ReadPara->Read( FileName );

   delete ReadPara;

// set the default explosion center
   for (int d=0; d<3; d++)
      if ( Pulsar_Center[d] < 0.0 )  Pulsar_Center[d] = 0.5*amr->BoxSize[d];


// (2) set the problem-specific derived parameters


// (3) reset other general-purpose parameters
//     --> a helper macro PRINT_WARNING is defined in TestProb.h
   const double End_T_Default    = 5.0e-3;
   const long   End_Step_Default = __INT_MAX__;

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
      Aux_Message( stdout, "  test problem ID                 = %d\n",     TESTPROB_ID );
      Aux_Message( stdout, "  ambient density                 = %13.7e\n", Pulsar_DensBg );
      Aux_Message( stdout, "  ambient pressure                = %13.7e\n", Pulsar_PresBg );
      Aux_Message( stdout, "  wind density/ambient density    = %13.7e\n", Pulsar_WindDensRatio );
      Aux_Message( stdout, "  wind velocity                   = %13.7e\n", Pulsar_WindVelocity );
      Aux_Message( stdout, "  wind pressure/ambient pressure  = %13.7e\n", Pulsar_WindPresRatio );
      Aux_Message( stdout, "  jet density/ambient density     = %13.7e\n", Pulsar_JetDensRatio );
      Aux_Message( stdout, "  jet velocity                    = %13.7e\n", Pulsar_JetVelocity );
      Aux_Message( stdout, "  jet pressure/ambient pressure   = %13.7e\n", Pulsar_JetPresRatio );
      Aux_Message( stdout, "  disk radius                     = %13.7e\n", Pulsar_DiskRadius );
      Aux_Message( stdout, "  disk height                     = %13.7e\n", Pulsar_DiskHeight );
      Aux_Message( stdout, "                                  = (%13.7e, %13.7e, %13.7e)\n",
                                                               Pulsar_Center[0], Pulsar_Center[1], Pulsar_Center[2] );
      Aux_Message( stdout, "=============================================================================\n" );
   }


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Setting runtime parameters ... done\n" );

} // FUNCTION : SetParameter



//-------------------------------------------------------------------------------------------------------
// Function    :  SetGridIC
// Description :  Set the problem-specific initial condition on grids
//
// Note        :  1. This function may also be used to estimate the numerical errors when OPT__OUTPUT_USER is enabled
//                   --> In this case, it should provide the analytical solution at the given "Time"
//                2. This function will be invoked by multiple OpenMP threads when OPENMP is enabled
//                   --> Please ensure that everything here is thread-safe
//                3. Even when DUAL_ENERGY is adopted for HYDRO, one does NOT need to set the dual-energy variable here
//                   --> It will be calculated automatically
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
   // primitive variables
   double Prim[NCOMP_FLUID];
   
   // shift origin to center of pulsar
   const double X  = x - Pulsar_Center[0];
   const double Y  = y - Pulsar_Center[1];
   const double Z  = z - Pulsar_Center[2];

   // transform to cylindrical coordinate 
   const double R   = SQRT( SQR( X ) + SQR( Y ) );
   const double Phi = ATAN2( X, Y ); 

   if ( R < Pulsar_DiskRadius && FABS( Z - 0.5*Pulsar_DiskHeight ) )
   {
      Prim[0] = Pulsar_WindDensRatio*Pulsar_DensBg;
      Prim[1] = Pulsar_WindVelocity*SIN(Phi);
      Prim[2] = Pulsar_WindVelocity*COS(Phi); 
      Prim[3] = 0.0;
      Prim[4] = Pulsar_WindPresRatio*Pulsar_PresBg;
   }
   else
   {
      Prim[0] = Pulsar_DensBg;
      Prim[1] = 0.0;
      Prim[2] = 0.0;
      Prim[3] = 0.0;
      Prim[4] = Pulsar_PresBg;
   }

   SRHydro_Pri2Con( Prim, fluid, GAMMA );

} // FUNCTION : SetGridIC
#endif // #if ( MODEL == SR_HYDRO )


// (true/false): if the target cell (is/is not) within the region to be refined
static bool Flag_Region( const int i, const int j, const int k, const int lv, const int PID )
{
   const double dh     = amr->dh[lv];                                                  // grid size
   const double Pos[3] = { amr->patch[0][lv][PID]->EdgeL[0] + (i+0.5)*dh,              // x,y,z position
                           amr->patch[0][lv][PID]->EdgeL[1] + (j+0.5)*dh,
                           amr->patch[0][lv][PID]->EdgeL[2] + (k+0.5)*dh  };

   bool Flag = false;  

   const double Center[3]      = { 0.5*amr->BoxSize[0], 0.5*amr->BoxSize[1], 0.5*amr->BoxSize[2] };
   const double dR[3]          = { Pos[0]-Center[0], Pos[1]-Center[1], Pos[2]-Center[2] };
   const double R              = sqrt( SQR(dR[0]) + SQR(dR[1]) + SQR(dR[2]) );



   if ( R < 1.1*Pulsar_DiskRadius )   return true;
   else                               return false;

} // FUNCTION : Flag_Region


//-------------------------------------------------------------------------------------------------------
// Function    :  Init_TestProb_Hydro_PulsarWave
// Description :  Test problem initializer
//
// Note        :  None
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Init_TestProb_SRHydro_PulsarWind()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );


// validate the compilation flags and runtime parameters
   Validate();


#  if ( MODEL == SR_HYDRO )
// set the problem-specific runtime parameters
   SetParameter();


// set the function pointers of various problem-specific routines
   Init_Function_User_Ptr   = SetGridIC;
   Output_User_Ptr          = NULL;
   Flag_Region_Ptr          = Flag_Region;
   Mis_GetTimeStep_User_Ptr = NULL;
   Aux_Record_User_Ptr      = NULL;
   BC_User_Ptr              = NULL;
   Flu_ResetByUser_Func_Ptr = NULL;
   End_User_Ptr             = NULL;
#  endif // #if ( MODEL == SR_HYDRO )


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : Init_TestProb_Hydro_PulsarWave

#endif
