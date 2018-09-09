#include "GAMER.h"
#include "TestProb.h"

# define PI 3.14159265

// problem-specific global variables
// =======================================================================================
static double Theta;     // incline angle of shock
static double Height;    // vertical shift
static double DensUp;    // primitive density in up-stream
static double VelyUp;    // magnitude of 3-velocity in up-stream
static double PresUp;    // pressure in up-stream
static double DensDown;  // primitive density in down-stream
static double VelyDown;  // magnitude of 3-velocity in down-stream
static double PresDown;  // pressure in down-stream
// =======================================================================================

void CPU_Pri2Con( const real In[], real Out[], const real Gamma);
void CPU_3Velto4Vel( const real In[], real Out[] );


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

// (1-1) add parameters in the following format:
// --> note that VARIABLE, DEFAULT, MIN, and MAX must have the same data type
// --> some handy constants (e.g., Useless_bool, Eps_double, NoMin_int, ...) are defined in "include/ReadPara.h"
// ********************************************************************************************************************************
// ReadPara->Add( "KEY_IN_THE_FILE",  &VARIABLE,          DEFAULT,        MIN,              MAX               );
// ********************************************************************************************************************************
   ReadPara->Add( "Theta",            &Theta,             -1.0,           Eps_double,       +90.0             );
//   ReadPara->Add( "Height",           &Height,            -1.0,           Eps_double,       NoMax_double      );
   ReadPara->Add( "DensUp",           &DensUp,            -1.0,           Eps_double,       NoMax_double      );
   ReadPara->Add( "VelyUp",           &VelyUp,            -1.0,           -1.0 ,            +1.0              );
   ReadPara->Add( "PresUp",           &PresUp,            -1.0,           Eps_double,       NoMax_double      );
   ReadPara->Add( "DensDown",         &DensDown,          -1.0,           Eps_double,       NoMax_double      );
   ReadPara->Add( "VelyDown",         &VelyDown,          -1.0,           -1.0      ,       +1.0              );
   ReadPara->Add( "PresDown",         &PresDown,          -1.0,           Eps_double,       NoMax_double      );

   ReadPara->Read( FileName );

   delete ReadPara;

// (1-2) set the default values
   Theta = 45;
   Height = 0.05*amr->BoxSize[1];

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
      Aux_Message( stdout, "  test problem ID                        = %d\n",     TESTPROB_ID );
      Aux_Message( stdout, "  incline angle of shock                 = %13.7e\n", Theta );
      Aux_Message( stdout, "  vertical shift                         = %13.7e\n", Height );
      Aux_Message( stdout, "  primitive density in up-stream         = %13.7e\n", DensUp );
      Aux_Message( stdout, "  magnitude of 3-velocity in up-stream   = %13.7e\n", VelyUp );
      Aux_Message( stdout, "  pressure in up-stream                  = %13.7e\n", PresUp );
      Aux_Message( stdout, "  primitive density in down-stream       = %13.7e\n", DensDown );
      Aux_Message( stdout, "  magnitude of 3-velocity in down-stream = %13.7e\n", VelyDown );
      Aux_Message( stdout, "  pressure in down-stream                = %13.7e\n", PresDown );
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
   real Prim1[NCOMP_FLUID]; // store 3-velocity
   real Prim2[NCOMP_FLUID]; // store 4-velocity

   double VelyUp_x = VelyUp*SIN(Theta*PI/180.0);
   double VelyUp_y = VelyUp*COS(Theta*PI/180.0);
   double VelyDown_x = VelyDown*SIN(Theta*PI/180.0);
   double VelyDown_y = VelyDown*COS(Theta*PI/180.0);
   double slope = tan(Theta*PI/180.0); // shock slope

   if ( y >= slope*x +  Height) // down-stream
   {
      Prim1[0] = DensDown;
      Prim1[1] = VelyDown_x;
      Prim1[2] = VelyDown_y;
      Prim1[3] = 0.0;
      Prim1[4] = PresDown;

      CPU_3Velto4Vel (Prim1, Prim2);
      CPU_Pri2Con (Prim2, fluid, GAMMA);
   }else{ // up-stream
      Prim1[0] = DensUp;
      Prim1[1] = VelyUp_x;
      Prim1[2] = VelyUp_y;
      Prim1[3] = 0.0;
      Prim1[4] = PresUp;

      CPU_3Velto4Vel (Prim1, Prim2);
      CPU_Pri2Con (Prim2, fluid, GAMMA);
   }

} // FUNCTION : SetGridIC
#endif // #if ( MODEL == SR_HYDRO )



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
