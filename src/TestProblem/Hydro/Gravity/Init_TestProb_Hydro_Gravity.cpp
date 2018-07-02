#include "GAMER.h"
#include "TestProb.h"



// problem-specific global variables
// =======================================================================================
typedef int DensProf_t;
const DensProf_t
   DENSPROF_NFW       = 1,
   DENSPROF_HERNQUIST = 2;

static DensProf_t Gra_DensProf;  // density profile
static double     Gra_Radius0;   // radius  parameter for the adopted density profile
static double     Gra_Dens0;     // density parameter for the adopted density profile
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
#  if ( MODEL != HYDRO )
   Aux_Error( ERROR_INFO, "MODEL != HYDRO !!\n" );
#  endif

#  ifndef GRAVITY
   Aux_Error( ERROR_INFO, "GRAVITY must be enabled !!\n" );
#  endif

#  ifdef PARTICLE
   Aux_Error( ERROR_INFO, "PARTICLE must be disabled !!\n" );
#  endif

#  ifdef COMOVING
   Aux_Error( ERROR_INFO, "COMOVING must be disabled !!\n" );
#  endif

#  ifdef DUAL_ENERGY
   Aux_Error( ERROR_INFO, "DUAL_ENERGY must be disabled\n" );
#  endif

   if ( END_STEP > 0 )  Aux_Error( ERROR_INFO, "do not set END_STEP for this test !!\n" );

   if ( END_T > 0.0 )   Aux_Error( ERROR_INFO, "do not set END_T for this test !!\n" );


// warnings
   if ( MPI_Rank == 0 )
   {
      if ( OPT__BC_POT == BC_POT_PERIODIC )
         Aux_Message( stderr, "WARNING : periodic BC for gravity (OPT__BC_POT=1) is not recommended for this test !!\n" );

      if ( !OPT__OUTPUT_POT )
         Aux_Message( stderr, "WARNING : turn on OPT__OUTPUT_POT to output gravitational potential !!\n" );
   }


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

// (1-1) add parameters in the following format:
// --> note that VARIABLE, DEFAULT, MIN, and MAX must have the same data type
// --> some handy constants (e.g., Useless_bool, Eps_double, NoMin_int, ...) are defined in "include/ReadPara.h"
// ********************************************************************************************************************************
// ReadPara->Add( "KEY_IN_THE_FILE",   &VARIABLE,              DEFAULT,       MIN,              MAX               );
// ********************************************************************************************************************************
   ReadPara->Add( "Gra_DensProf",      &Gra_DensProf,          2,             1,                2                 );
   ReadPara->Add( "Gra_Radius0",       &Gra_Radius0,           1.0,           Eps_double,       NoMax_double      );
   ReadPara->Add( "Gra_Dens0",         &Gra_Dens0,             1.0,           Eps_double,       NoMax_double      );

   ReadPara->Read( FileName );

   delete ReadPara;

// (1-2) set the default values

// (1-3) check the runtime parameters


// (2) set the problem-specific derived parameters


// (3) reset other general-purpose parameters
//     --> a helper macro PRINT_WARNING is defined in TestProb.h
   const long   End_Step_Default = 0;
   const double End_T_Default    = 0.0;

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
      Aux_Message( stdout, "  test problem ID = %d\n",     TESTPROB_ID  );
      Aux_Message( stdout, "  Gra_DensProf    = %d\n",     Gra_DensProf );
      Aux_Message( stdout, "  Gra_Radius0     = %13.7e\n", Gra_Radius0  );
      Aux_Message( stdout, "  Gra_Dens0       = %13.7e\n", Gra_Dens0    );
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

   const double r = sqrt( SQR(x-amr->BoxCenter[0]) + SQR(y-amr->BoxCenter[1]) + SQR(z-amr->BoxCenter[2]) );
   const double s = r / Gra_Radius0;

// only need to set density here
   switch ( Gra_DensProf )
   {
      case DENSPROF_NFW :
         fluid[DENS] = Gra_Dens0 / ( s*SQR(1.0+s) );
         break;

      case DENSPROF_HERNQUIST :
         fluid[DENS] = Gra_Dens0 / ( s*CUBE(1.0+s) );
         break;

      default :
         Aux_Error( ERROR_INFO, "unsupported density profile (Gra_DensProf = %d) !!\n", Gra_DensProf );
         exit( -1 );
   }

   fluid[MOMX] = NULL_REAL;
   fluid[MOMY] = NULL_REAL;
   fluid[MOMZ] = NULL_REAL;
   fluid[ENGY] = NULL_REAL;

} // FUNCTION : SetGridIC
#endif // #if ( MODEL == HYDRO )



//-------------------------------------------------------------------------------------------------------
// Function    :  Init_TestProb_Hydro_Gravity
// Description :  Test problem initializer
//
// Note        :  None
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Init_TestProb_Hydro_Gravity()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );


// validate the compilation flags and runtime parameters
   Validate();


#  if ( MODEL == HYDRO )
// set the problem-specific runtime parameters
   SetParameter();


   Init_Function_User_Ptr = SetGridIC;
   Output_User_Ptr        = NULL;
#  endif // #if ( MODEL == HYDRO )


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : Init_TestProb_Hydro_Gravity
