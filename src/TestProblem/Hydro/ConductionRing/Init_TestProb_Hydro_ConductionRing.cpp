#include "GAMER.h"


// problem-specific global variables
// =======================================================================================
static double Ring_U1;               // internal energy in ring
static double Ring_U2;               // internal energy outside ring
static double Ring_R1;               // inner radius of ring
static double Ring_R2;               // outer radius of ring
static double Ring_Rho;              // background mass density
static double Ring_Angle;            // angle defining initial region of U2

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


#  if ( MODEL != HYDRO )
   Aux_Error( ERROR_INFO, "MODEL != HYDRO !!\n" );
#  endif

#  ifndef MHD
   Aux_Error( ERROR_INFO, "MHD must be enabled !!\n" );
#  endif

#  ifndef CONDUCTION
   Aux_Error( ERROR_INFO, "CONDUCTION must be enabled !!\n" );
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

   if ( !OPT__FREEZE_HYDRO )
      Aux_Error( ERROR_INFO, "OPT__FREEZE_HYDRO must be enabled !!\n" );

   if ( CONDUCTION_FLUX_TYPE != ANISOTROPIC_CONDUCTION )
      Aux_Error( ERROR_INFO, "please set \"CONDUCTION_FLUX_TYPE = 2\" (i.e., anisotropic conduction) !!\n" );
   
   for (int f=0; f<6; f++)
   if ( OPT__BC_FLU[f] != BC_FLU_OUTFLOW )
      Aux_Error( ERROR_INFO, "please set \"OPT__BC_FLU_* = 2\" (i.e., outflow BC) !!\n" );

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Validating test problem %d ... done\n", TESTPROB_ID );

} // FUNCTION : Validate



#if ( MODEL == HYDRO )
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
// ReadPara->Add( "KEY_IN_THE_FILE",     &VARIABLE_ADDRESS,      DEFAULT,       MIN,              MAX               );
// ********************************************************************************************************************************
   ReadPara->Add( "Ring_U1",             &Ring_U1,               10.0,          Eps_double,       NoMax_double      );
   ReadPara->Add( "Ring_U2",             &Ring_U2,               12.0,          Eps_double,       NoMax_double      );
   ReadPara->Add( "Ring_Rho",            &Ring_Rho,               1.0,          Eps_double,       NoMax_double      );
   ReadPara->Add( "Ring_R1",             &Ring_R1,                0.5,          0.0,              NoMax_double      );
   ReadPara->Add( "Ring_R2",             &Ring_R2,                0.7,          0.0,              NoMax_double      );
   ReadPara->Add( "Ring_Angle",          &Ring_Angle,            15.0,          0.0,              360.0             );

   ReadPara->Read( FileName );

   delete ReadPara;

// (2) set the problem-specific derived parameters


// (3) reset other general-purpose parameters
//     --> a helper macro PRINT_RESET_PARA is defined in Macro.h
   const double End_T_Default    = 200.0;
   const long   End_Step_Default = __INT_MAX__;

   if ( END_STEP < 0 ) {
      END_STEP = End_Step_Default;
      PRINT_RESET_PARA( END_STEP, FORMAT_LONG, "" );
   }

   if ( END_T < 0.0 ) {
      END_T = End_T_Default;
      PRINT_RESET_PARA( END_T, FORMAT_REAL, "" );
   }


// (4) make a note
   if ( MPI_Rank == 0 )
   {
      Aux_Message( stdout, "=============================================================================\n" );
      Aux_Message( stdout, "  Ring_U1              = % 14.7e\n", Ring_U1 );
      Aux_Message( stdout, "  Ring_U2              = % 14.7e\n", Ring_U2 );
      Aux_Message( stdout, "  Ring_R1              = % 14.7e\n", Ring_R1 );
      Aux_Message( stdout, "  Ring_R2              = % 14.7e\n", Ring_R2 );
      Aux_Message( stdout, "  Ring_Rho             = % 14.7e\n", Ring_Rho );
      Aux_Message( stdout, "  Ring_Angle           = % 14.7e\n", Ring_Angle );
      Aux_Message( stdout, "=============================================================================\n" );
   }

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Setting runtime parameters ... done\n" );

} // FUNCTION : SetParameter



//-------------------------------------------------------------------------------------------------------
// Function    :  SetGridIC
// Description :  Set the problem-specific initial condition on grids
//
// Note        :  This function may also be used to estimate the numerical errors when OPT__OUTPUT_USER is enabled
//                --> In this case, it should provide the analytical solution at the given "Time"
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

   const double x0         = 0.5*amr->BoxSize[0];
   const double y0         = 0.5*amr->BoxSize[0];
   const double angle_rads = (real)Ring_Angle*M_PI/180.0;
   const double r          = sqrt( SQR(x-x0) + SQR(y-y0) );
   const double phi        = atan2( y-y0, x-x0 );

   double Dens, MomX, MomY, MomZ, Eint, Etot;

// inside region
   if ( 0.5 < r && r < 0.7 && FABS(phi) < angle_rads )
   {
      Eint = Ring_U1*Ring_Rho;
   }
// outside region
   else
   {
      Eint = Ring_U2*Ring_Rho;
   }

   Dens = Ring_Rho;
   MomX = 0.0;
   MomY = 0.0;
   MomZ = 0.0;

// compute the total gas energy
   Etot = Hydro_ConEint2Etot( Dens, MomX, MomY, MomZ, Eint, 0.0 );     // do NOT include magnetic energy here

// set the output array
   fluid[DENS] = Dens;
   fluid[MOMX] = MomX;
   fluid[MOMY] = MomY;
   fluid[MOMZ] = MomZ;
   fluid[ENGY] = Etot;

} // FUNCTION : SetGridIC

//-------------------------------------------------------------------------------------------------------
// Function    :  SetBFieldIC
// Description :  Set the problem-specific initial condition of magnetic field
//
// Note        :  1. This function will be invoked by multiple OpenMP threads when OPENMP is enabled
//                   (unless OPT__INIT_GRID_WITH_OMP is disabled)
//                   --> Please ensure that everything here is thread-safe
//
// Parameter   :  magnetic : Array to store the output magnetic field
//                x/y/z    : Target physical coordinates
//                Time     : Target physical time
//                lv       : Target refinement level
//                AuxArray : Auxiliary array
//
// Return      :  magnetic
//-------------------------------------------------------------------------------------------------------
void SetBFieldIC( real magnetic[], const double x, const double y, const double z, const double Time,
                  const int lv, double AuxArray[] )
{

   const double x0         = 0.5*amr->BoxSize[0];
   const double y0         = 0.5*amr->BoxSize[0];
   const double r          = sqrt( SQR(x-x0) + SQR(y-y0) );
   const double phi        = atan2( y-y0, x-x0 );

   magnetic[MAGX] = - ( y - y0 ) * 0.5/sqrt(M_PI) / r;
   magnetic[MAGY] =   ( x - x0 ) * 0.5/sqrt(M_PI) / r;
   magnetic[MAGZ] = 0.0;

} // FUNCTION : SetBFieldIC
#endif // #if ( MODEL == HYDRO )


//-------------------------------------------------------------------------------------------------------
// Function    :  Init_TestProb_Hydro_ConductionRing
// Description :  Test problem initializer
//
// Note        :  None
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Init_TestProb_Hydro_ConductionRing()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );


// validate the compilation flags and runtime parameters
   Validate();


#  if ( MODEL == HYDRO )
// set the problem-specific runtime parameters
   SetParameter();


// set the function pointers of various problem-specific routines
   Init_Function_User_Ptr = SetGridIC;
#  ifdef MHD
   Init_Function_BField_User_Ptr = SetBFieldIC;
#  endif
#  endif // #if ( MODEL == HYDRO )


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : Init_TestProb_Hydro_ConductionRing
