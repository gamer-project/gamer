#include "GAMER.h"
#include "TestProb.h"


static void BC( real fluid[], const double x, const double y, const double z, const double Time,
                const int lv, double AuxArray[] );
static void Init_ExtPot();


// problem-specific global variables
// =======================================================================================
static double ExtPot_Amp;        // initial wave function amplitude
static double ExtPot_M;          // point source mass
static double ExtPot_Cen[3];     // point source position
// =======================================================================================

// external potential routines
void SetCPUExtPot_PointMass( ExtPot_t &CPUExtPot_Ptr );
# ifdef GPU
void SetGPUExtPot_PointMass( ExtPot_t &GPUExtPot_Ptr );
# endif




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


#  if ( MODEL != ELBDM )
   Aux_Error( ERROR_INFO, "MODEL != ELBDM !!\n" );
#  endif

#  ifndef GRAVITY
   Aux_Error( ERROR_INFO, "GRAVITY must be enabled !!\n" );
#  endif

#  ifdef COMOVING
   Aux_Error( ERROR_INFO, "COMOVING must be disabled !!\n" );
#  endif

#  ifdef PARTICLE
   Aux_Error( ERROR_INFO, "PARTICLE must be disabled !!\n" );
#  endif

#  ifdef GRAVITY
   if ( !OPT__EXT_POT )
   Aux_Error( ERROR_INFO, "OPT__EXT_POT must be enabled !!\n" );
#  endif


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Validating test problem %d ... done\n", TESTPROB_ID );

} // FUNCTION : Validate



#if ( MODEL == ELBDM  &&  defined GRAVITY )
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
// ReadPara->Add( "KEY_IN_THE_FILE",   &VARIABLE,              DEFAULT,       MIN,              MAX               );
// ********************************************************************************************************************************
   ReadPara->Add( "ExtPot_Amp",        &ExtPot_Amp,            -1.0,          Eps_double,       NoMax_double      );
   ReadPara->Add( "ExtPot_M",          &ExtPot_M,              -1.0,          Eps_double,       NoMax_double      );
   ReadPara->Add( "ExtPot_Cen_X",      &ExtPot_Cen[0],         -1.0,          NoMin_double,     NoMax_double      );
   ReadPara->Add( "ExtPot_Cen_Y",      &ExtPot_Cen[1],         -1.0,          NoMin_double,     NoMax_double      );
   ReadPara->Add( "ExtPot_Cen_Z",      &ExtPot_Cen[2],         -1.0,          NoMin_double,     NoMax_double      );

   ReadPara->Read( FileName );

   delete ReadPara;

// set the default center
   for (int d=0; d<3; d++)
      if ( ExtPot_Cen[d] < 0.0 )    ExtPot_Cen[d] = 0.5*amr->BoxSize[d];


// (2) set the problem-specific derived parameters


// (3) reset other general-purpose parameters
//     --> a helper macro PRINT_RESET_PARA is defined in Macro.h
   const long   End_Step_Default = __INT_MAX__;
   const double End_T_Default    = 1.0e-2;

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
      Aux_Message( stdout, "  test problem ID         = %d\n",                       TESTPROB_ID );
      Aux_Message( stdout, "  wave function amplitude = %13.7e\n",                   ExtPot_Amp );
      Aux_Message( stdout, "  point source mass       = %13.7e\n",                   ExtPot_M );
      Aux_Message( stdout, "  point source position   = (%13.7e, %13.7e, %13.7e)\n", ExtPot_Cen[0],
                                                                                     ExtPot_Cen[1],
                                                                                     ExtPot_Cen[2] );
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

   const double r     = sqrt( SQR(x-ExtPot_Cen[0]) + SQR(y-ExtPot_Cen[1]) + SQR(z-ExtPot_Cen[2]) );
   const double Coeff = 2.0*SQR(ELBDM_ETA)*NEWTON_G*ExtPot_M;
   const double R     = sqrt( Coeff*r );

   fluid[REAL] = ExtPot_Amp*j1( 2.0*R )/R;
   fluid[IMAG] = 0.0;                                       // imaginary part is always zero --> no initial velocity
   fluid[DENS] = SQR( fluid[REAL] ) + SQR( fluid[IMAG] );

} // FUNCTION : SetGridIC



//-------------------------------------------------------------------------------------------------------
// Function    :  BC
// Description :  Set the extenral boundary condition to the analytical solution
//
// Note        :  1. Linked to the function pointer "BC_User_Ptr"
//
// Parameter   :  fluid    : Fluid field to be set
//                x/y/z    : Physical coordinates
//                Time     : Physical time
//                lv       : Refinement level
//                AuxArray : Auxiliary array
//
// Return      :  fluid
//-------------------------------------------------------------------------------------------------------
void BC( real fluid[], const double x, const double y, const double z, const double Time,
         const int lv, double AuxArray[] )
{

   SetGridIC( fluid, x, y, z, Time, lv, AuxArray );

} // FUNCTION : BC



//-------------------------------------------------------------------------------------------------------
// Function    :  Init_ExtPotAuxArray_ExtPotTest
// Description :  Set the auxiliary array ExtPot_AuxArray[] used by the external potential routine
//
// Note        :  1. External potential can be enabled by the runtime option "OPT__EXT_POT"
//                2. Link to the function pointer "Init_ExtPotAuxArray_Ptr"
//                3. AuxArray[] has the size of EXT_POT_NAUX_MAX defined in Macro.h (default = 20)
//
// Parameter   :  AuxArray : Array to be filled up
//
// Return      :  AuxArray[]
//-------------------------------------------------------------------------------------------------------
void Init_ExtPotAuxArray_ExtPotTest( double AuxArray[] )
{

   AuxArray[0] = ExtPot_Cen[0];
   AuxArray[1] = ExtPot_Cen[1];
   AuxArray[2] = ExtPot_Cen[2];
   AuxArray[3] = ExtPot_M*NEWTON_G;

} // FUNCTION : Init_ExtPotAuxArray_ExtPotTest
#endif // #if ( MODEL == ELBDM  &&  defined GRAVITY )



//-------------------------------------------------------------------------------------------------------
// Function    :  Init_TestProb_ELBDM_ExtPot
// Description :  Test problem initializer
//
// Note        :  None
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Init_TestProb_ELBDM_ExtPot()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );


// validate the compilation flags and runtime parameters
   Validate();


#  if ( MODEL == ELBDM  &&  defined GRAVITY )
// set the problem-specific runtime parameters
   SetParameter();


// set the function pointers of various problem-specific routines
   Init_Function_User_Ptr  = SetGridIC;
   BC_User_Ptr             = BC;
   Init_ExtPotAuxArray_Ptr = Init_ExtPotAuxArray_ExtPotTest;
   SetCPUExtPot_Ptr        = SetCPUExtPot_PointMass;
#  ifdef GPU
   SetGPUExtPot_Ptr        = SetGPUExtPot_PointMass;
#  endif
#  endif // #if ( MODEL == ELBDM  &&  defined GRAVITY )


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : Init_TestProb_ELBDM_ExtPot
