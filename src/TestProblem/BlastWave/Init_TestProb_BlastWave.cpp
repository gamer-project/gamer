#include "GAMER.h"
#include "TestProb.h"



// problem-specific global variables
// =======================================================================================
static double Blast_Dens_Bg;       // background mass density
static double Blast_Engy_Bg;       // background energy density
static double Blast_Engy_Exp;      // total explosion energy
static double Blast_Radius;        // explosion radius
static double Blast_Center[3];     // explosion center
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

#  ifdef GRAVITY
   Aux_Error( ERROR_INFO, "GRAVITY must be disabled !!\n" );
#  endif

#  ifdef COMOVING
   Aux_Error( ERROR_INFO, "COMOVING must be disabled !!\n" );
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

// add parameters in the following format (some handy constants are defined in TestProb.h):
// ********************************************************************************************************************************
// ReadPara->Add( "KEY_IN_THE_FILE",   &VARIABLE_ADDRESS,      DEFAULT,       MIN,              MAX               );
// ********************************************************************************************************************************
   ReadPara->Add( "Blast_Dens_Bg",     &Blast_Dens_Bg,         -1.0,          NoZero_double,    NoMax_double      );
   ReadPara->Add( "Blast_Engy_Bg",     &Blast_Engy_Bg,         -1.0,          NoZero_double,    NoMax_double      );
   ReadPara->Add( "Blast_Engy_Exp",    &Blast_Engy_Exp,        -1.0,          NoZero_double,    NoMax_double      );
   ReadPara->Add( "Blast_Radius",      &Blast_Radius,          -1.0,          NoZero_double,    NoMax_double      );
   ReadPara->Add( "Blast_Center_X",    &Blast_Center[0],       -1.0,          NoMin_double,     amr->BoxSize[0]   );
   ReadPara->Add( "Blast_Center_Y",    &Blast_Center[1],       -1.0,          NoMin_double,     amr->BoxSize[1]   );
   ReadPara->Add( "Blast_Center_Z",    &Blast_Center[2],       -1.0,          NoMin_double,     amr->BoxSize[2]   );

   ReadPara->Read( FileName );

   delete ReadPara;

// set the default explosion center
   for (int d=0; d<3; d++)
      if ( Blast_Center[d] < 0.0 )  Blast_Center[d] = 0.5*amr->BoxSize[d];


// (2) set the problem-specific derived parameters


// (3) reset other general-purpose parameters
//     --> a helper macro PRINT_WARNING is defined in TestProb.h
   const double End_T_Default    = 3.0e-2;
   const long   End_Step_Default = __INT_MAX__;

   if ( END_STEP < 0 ) {
      END_STEP = End_Step_Default;
      PRINT_WARNING( "END_STEP", END_STEP, FORMAT_LONG );
   }

   if ( END_T < 0.0 ) {
      END_T = End_T_Default;
      PRINT_WARNING( "END_T", END_T, FORMAT_REAL );
   }

   if ( !OPT__INIT_RESTRICT ) {
      OPT__INIT_RESTRICT = true;
      PRINT_WARNING( "OPT__INIT_RESTRICT", OPT__INIT_RESTRICT, FORMAT_BOOL );
   }


// (4) make a note
   if ( MPI_Rank == 0 )
   {
      Aux_Message( stdout, "=============================================================================\n" );
      Aux_Message( stdout, "  test problem ID           = %d\n",     TESTPROB_ID );
      Aux_Message( stdout, "  background mass density   = %13.7e\n", Blast_Dens_Bg );
      Aux_Message( stdout, "  background energy density = %13.7e\n", Blast_Engy_Bg);
      Aux_Message( stdout, "  total explosion energy    = %13.7e\n", Blast_Engy_Exp );
      Aux_Message( stdout, "  explosion energy density  = %13.7e\n", Blast_Engy_Exp/(4.0*M_PI/3.0*CUBE(Blast_Radius) ) );
      Aux_Message( stdout, "  explosion radius          = %13.7e\n", Blast_Radius );
      Aux_Message( stdout, "  explosion center          = (%13.7e, %13.7e, %13.7e)\n", Blast_Center[0], Blast_Center[1],
                                                                                       Blast_Center[2] );
      Aux_Message( stdout, "=============================================================================\n" );
   }


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Setting runtime parameters ... done\n" );

} // FUNCTION : SetParameter



//-------------------------------------------------------------------------------------------------------
// Function    :  SetGridIC
// Description :  Set the problem-specific initial condition on grids
//
// Note        :  1. This function will also be used to estimate the numerical errors when OPT__OUTPUT_TEST_ERROR is enabled
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

   const double Blast_Engy_Exp_Density = Blast_Engy_Exp/(4.0*M_PI/3.0*Blast_Radius*Blast_Radius*Blast_Radius);
   const double r = SQRT( SQR(x-Blast_Center[0]) + SQR(y-Blast_Center[1]) + SQR(z-Blast_Center[2]) );

   fluid[DENS] = Blast_Dens_Bg;
   fluid[MOMX] = 0.0;
   fluid[MOMY] = 0.0;
   fluid[MOMZ] = 0.0;

   if ( r <= Blast_Radius )   fluid[ENGY] = Blast_Engy_Exp_Density;
   else                       fluid[ENGY] = Blast_Engy_Bg;

} // FUNCTION : SetGridIC



//-------------------------------------------------------------------------------------------------------
// Function    :  Init_TestProb_BlastWave
// Description :  Test problem initializer
//
// Note        :  None
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Init_TestProb_BlastWave()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );


// validate the compilation flags and runtime parameters
   Validate();


// set the problem-specific runtime parameters
   SetParameter();


// set various problem-specific routines
   Init_Function_User_Ptr = SetGridIC;
   Output_User_Ptr        = NULL;


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : Init_TestProb_BlastWave
