#include "GAMER.h"
#include "TestProb.h"



// problem-specific global variables
// =======================================================================================
static double InclinedAngle;           // incline angle of shock
static double DensUpStream;            // proper mass density in up-stream
static double UxUpStream;              // x component of 4-velocity in up-stream
static double UyUpStream;              // y component of 4-velocity in up-stream
static double UzUpStream;              // z component of 4-velocity in up-stream
static double UxDownStream;            // x component of 4-velocity in down-stream
static double UyDownStream;            // y component of 4-velocity in down-stream
static double UzDownStream;            // z component of 4-velocity in down-stream
static double PresUpstream;            // pressure in up-stream
static double DensDownStream;          // proper mass density in down-stream
static double PresDownStream;          // pressure in down-stream
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
   Aux_Error( ERROR_INFO, "MODEL != SR_HYDRO !!\n" );
#  endif

#  ifdef GRAVITY
   Aux_Error( ERROR_INFO, "GRAVITY must be disabled !!\n" );
#  endif

#  ifdef PARTICLE
   Aux_Error( ERROR_INFO, "PARTICLE must be disabled !!\n" );
#  endif

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
// ReadPara->Add( "KEY_IN_THE_FILE",        &VARIABLE,             DEFAULT,                  MIN,                MAX      );
// ********************************************************************************************************************************
   ReadPara->Add( "InclinedAngle",          &InclinedAngle,           -1.0,         NoMin_double,                   +90.0 );
   ReadPara->Add( "DensUpStream",           &DensUpStream,            -1.0,           Eps_double,            NoMax_double );
   ReadPara->Add( "UxUpStream",             &UxUpStream,              -1.0,         NoMin_double,            NoMax_double );
   ReadPara->Add( "UyUpStream",             &UyUpStream,              -1.0,         NoMin_double,            NoMax_double );
   ReadPara->Add( "UzUpStream",             &UzUpStream,              -1.0,         NoMin_double,            NoMax_double );
   ReadPara->Add( "PresUpstream",           &PresUpstream,            -1.0,           Eps_double,            NoMax_double );
   ReadPara->Add( "DensDownStream",         &DensDownStream,          -1.0,           Eps_double,            NoMax_double );
   ReadPara->Add( "UxDownStream",           &UxDownStream,            -1.0,         NoMin_double,            NoMax_double );
   ReadPara->Add( "UyDownStream",           &UyDownStream,            -1.0,         NoMin_double,            NoMax_double );
   ReadPara->Add( "UzDownStream",           &UzDownStream,            -1.0,         NoMin_double,            NoMax_double );
   ReadPara->Add( "PresDownStream",         &PresDownStream,          -1.0,           Eps_double,            NoMax_double );

   ReadPara->Read( FileName );

   delete ReadPara;

// (1-2) set the default values
   if ( InclinedAngle == -1.0 )
        InclinedAngle = 45.0;

// (1-3) check the runtime parameters


// (2) set the problem-specific derived parameters


// (3) reset other general-purpose parameters
//     --> a helper macro PRINT_WARNING is defined in TestProb.h
   const long   End_Step_Default = __INT_MAX__;
   const double End_T_Default    = __FLT_MAX__;

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
      Aux_Message( stdout, "  test problem ID                          = %d\n",       TESTPROB_ID       );
      Aux_Message( stdout, "  inclined angle of shock                  = %13.7e\n",   InclinedAngle     );
      Aux_Message( stdout, "  proper mass density in up-stream         = %13.7e\n",   DensUpStream      );
      Aux_Message( stdout, "  x component of 4-velocity in up-stream   = %13.7e\n",   UxUpStream        );
      Aux_Message( stdout, "  y component of 4-velocity in up-stream   = %13.7e\n",   UyUpStream        );
      Aux_Message( stdout, "  z component of 4-velocity in up-stream   = %13.7e\n",   UzUpStream        );
      Aux_Message( stdout, "  pressure in up-stream                    = %13.7e\n",   PresUpstream      );
      Aux_Message( stdout, "  proper mass density in down-stream       = %13.7e\n",   DensDownStream    );
      Aux_Message( stdout, "  x component of 4-velocity in down-stream = %13.7e\n",   UxDownStream      );
      Aux_Message( stdout, "  y component of 4-velocity in down-stream = %13.7e\n",   UyDownStream      );
      Aux_Message( stdout, "  z component of 4-velocity in down-stream = %13.7e\n",   UzDownStream      );
      Aux_Message( stdout, "  pressure in down-stream                  = %13.7e\n",   PresDownStream    );
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
//                3. Even when DUAL_ENERGY is adopted for SR_HYDRO, one does NOT need to set the dual-energy variable here
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
   double Prim[NCOMP_FLUID];
   double BoxSizeX = amr->BoxSize[0];

   //if ( y >= tan( InclinedAngle*M_PI/180.0 )*x - 1e-2*BoxSizeX ) // down-stream
   if ( y >= tan( InclinedAngle*M_PI/180.0 )*x ) // down-stream
   {
      Prim[0] = DensDownStream;
      Prim[1] = UxDownStream;
      Prim[2] = UyDownStream;
      Prim[3] = UzDownStream;
      Prim[4] = PresDownStream;
   }
   else                                                            // up-stream
   {
      Prim[0] = DensUpStream;
      Prim[1] = UxUpStream;
      Prim[2] = UyUpStream;
      Prim[3] = UzUpStream;
      Prim[4] = PresUpstream;
   }


// cast double to real
#  ifndef FLOAT8
   double Out[NCOMP_FLUID];

   SRHydro_Pri2Con (Prim, Out, GAMMA);

   fluid [0] = (real) Out[0];
   fluid [1] = (real) Out[1];
   fluid [2] = (real) Out[2];
   fluid [3] = (real) Out[3];
   fluid [4] = (real) Out[4];
#  else
   SRHydro_Pri2Con (Prim, fluid, GAMMA);
#  endif

} // FUNCTION : SetGridIC

//-------------------------------------------------------------------------------------------------------
// Function    :  BC_User
// Description :  User-specified boundary condition
//
// Note        :  1. Invoked by "Flu_BoundaryCondition_User" using the function pointer "BC_User_Ptr"
//                   --> The function pointer may be reset by various test problem initializers, in which case
//                       this funtion will become useless
//                2. Always return NCOMP_TOTAL variables
//                3. Enabled by the runtime options "OPT__BC_FLU_* == 4"
//
// Parameter   :  fluid    : Fluid field to be set
//                x/y/z    : Physical coordinates
//                Time     : Physical time
//                lv       : Refinement level
//                AuxArray : Auxiliary array
//
// Return      :  fluid
//-------------------------------------------------------------------------------------------------------
void BC_User( real fluid[], const double x, const double y, const double z, const double Time,
              const int lv, double AuxArray[] )
{


}

//-------------------------------------------------------------------------------------------------------
// Function    :  Init_TestProb_SRHydro_DoubleMachReflection
// Description :  Test problem initializer
//
// Note        :  None
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Init_TestProb_SRHydro_DoubleMachReflection()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );


// validate the compilation flags and runtime parameters
   Validate();


// replace SR_HYDRO by the target model (e.g., MHD/ELBDM) and also check other compilation flags if necessary (e.g., GRAVITY/PARTICLE)
#  if ( MODEL == SR_HYDRO )
// set the problem-specific runtime parameters
   SetParameter();


// procedure to enable a problem-specific function:
// 1. define a user-specified function (example functions are given below)
// 2. declare its function prototype on the top of this file
// 3. set the corresponding function pointer below to the new problem-specific function
// 4. enable the corresponding runtime option in "Input__Parameter"
//    --> for instance, enable OPT__OUTPUT_USER for Output_User_Ptr
   Init_Function_User_Ptr      = SetGridIC;
   Init_Field_User_Ptr         = NULL;    // set NCOMP_PASSIVE_USER;        example: TestProblem/Hydro/Plummer/Init_TestProb_Hydro_Plummer.cpp --> AddNewField()
   Flag_User_Ptr               = NULL;    // option: OPT__FLAG_USER;        example: Refine/Flag_User.cpp
   Mis_GetTimeStep_User_Ptr    = NULL;    // option: OPT__DT_USER;          example: Miscellaneous/Mis_GetTimeStep_User.cpp
   BC_User_Ptr                 = NULL;    // option: OPT__BC_FLU_*=4;       example: TestProblem/ELBDM/ExtPot/Init_TestProb_ELBDM_ExtPot.cpp --> BC()
   Flu_ResetByUser_Func_Ptr    = NULL;    // option: OPT__RESET_FLUID;      example: Fluid/Flu_ResetByUser.cpp
   Output_User_Ptr             = NULL;    // option: OPT__OUTPUT_USER;      example: TestProblem/Hydro/AcousticWave/Init_TestProb_Hydro_AcousticWave.cpp --> OutputError()
   Aux_Record_User_Ptr         = NULL;    // option: OPT__RECORD_USER;      example: Auxiliary/Aux_Record_User.cpp
   End_User_Ptr                = NULL;    // option: none;                  example: TestProblem/Hydro/ClusterMerger_vs_Flash/Init_TestProb_ClusterMerger_vs_Flash.cpp --> End_ClusterMerger()
#  endif // #if ( MODEL == SR_HYDRO )


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : Init_TestProb_SRHydro_DoubleMachReflection
