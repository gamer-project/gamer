#include "GAMER.h"
#include "TestProb.h"



// problem-specific global variables
// =======================================================================================
static double var_double;
static int    var_int;
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


// examples
/*
#  if ( MODEL != HYDRO )
   Aux_Error( ERROR_INFO, "MODEL != HYDRO !!\n" );
#  endif

   if ( amr->BoxSize[0] != amr->BoxSize[1]  ||  amr->BoxSize[0] != amr->BoxSize[2] )
      Aux_Error( ERROR_INFO, "simulation domain must be CUBIC !!\n" );
*/


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Validating test problem %d ... done\n", TESTPROB_ID );

} // FUNCTION : Validate



// replace XXX by the target model (e.g., HYDRO/ELBDM)
#if ( MODEL == XXX )
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

// add parameters in the following format (some handy constants are defined in TestProb.h):
// --> note that VARIABLE, DEFAULT, MIN, and MAX must have the same data type
// --> some handy constants (e.g., NoMin_int, Eps_float, ...) are defined in "ReadPara.h"
// ********************************************************************************************************************************
// ReadPara->Add( "KEY_IN_THE_FILE",   &VARIABLE,              DEFAULT,       MIN,              MAX               );
// ********************************************************************************************************************************
   ReadPara->Add( "var_double",        &var_double,            -1.0,          NoMin_double,     NoMax_double      );
   ReadPara->Add( "var_int",           &var_int,               -1,            NoMin_int,        NoMax_int         );

   ReadPara->Read( FileName );

   delete ReadPara;

// set the default values


// (2) set the problem-specific derived parameters


// (3) reset other general-purpose parameters
//     --> a helper macro PRINT_WARNING is defined in TestProb.h
   const long   End_Step_Default = 10;
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
      Aux_Message( stdout, "  test problem ID           = %d\n",     TESTPROB_ID );
      Aux_Message( stdout, "  var_double                = %13.7e\n", var_double );
      Aux_Message( stdout, "  var_int                   = %d\n",     var_int );
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
// Parameter   :  fluid : Fluid field to be initialized
//                x/y/z : Physical coordinates
//                Time  : Physical time
//
// Return      :  fluid
//-------------------------------------------------------------------------------------------------------
void SetGridIC( real *fluid, const double x, const double y, const double z, const double Time )
{

   /*
// HYDRO example
   fluid[DENS] = 1.0;
   fluid[MOMX] = 0.0;
   fluid[MOMY] = 0.0;
   fluid[MOMZ] = 0.0;
   fluid[ENGY] = 1.0 + 0.5*( SQR(fluid[MOMX]) + SQR(fluid[MOMY]) + SQR(fluid[MOMZ]) ) / fluid[DENS];
   */

   /*
// ELBDM example
   fluid[REAL] = 1.0;
   fluid[IMAG] = 2.0;
   fluid[DENS] = SQR(fluid[REAL]) + SQR(fluid[IMAG])
   */

} // FUNCTION : SetGridIC
#endif // #if ( MODEL == XXX )



//-------------------------------------------------------------------------------------------------------
// Function    :  Init_TestProb_Template
// Description :  Test problem initializer
//
// Note        :  None
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Init_TestProb_Template()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );


// validate the compilation flags and runtime parameters
   Validate();


// replace XXX by the target model (e.g., HYDRO/ELBDM)
#  if ( MODEL == XXX )
// set the problem-specific runtime parameters
   SetParameter();


// set the function pointers of various problem-specific routines
   Init_Function_User_Ptr   = SetGridIC;
   Output_User_Ptr          = NULL;       // example: Hydro/AcousticWave/Init_TestProb_AcousticWave.cpp --> OutputError()
   Flag_User_Ptr            = NULL;       // example: Hydro/Bondi/Flag_Bondi.cpp
   Mis_GetTimeStep_User_Ptr = NULL;
   Aux_Record_User_Ptr      = NULL;       // example: Hydro/Bondi/Record_Bondi.cpp
   BC_User_Ptr              = NULL;
   Flu_ResetByUser_Func_Ptr = NULL;
   End_User_Ptr             = NULL;
#  ifdef GRAVITY
   Init_ExternalAcc_Ptr     = NULL;       // example: Hydro/Bondi/Init_ExternalAcc_Bondi.cpp
   Init_ExternalPot_Ptr     = NULL;
#  endif
#  endif // #if ( MODEL == XXX )


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : Init_TestProb_Template
