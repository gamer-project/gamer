#include "GAMER.h"
#include "Par_EquilibriumIC.h"
#include "string"

using namespace std;

// negligibly small uniform density and energy
double ParEqmIC_SmallGas;

// problem-specific function prototypes
#ifdef MASSIVE_PARTICLES
void Par_Init_ByFunction_ParEqmIC( const long NPar_ThisRank, const long NPar_AllRank,
                                   real_par *ParMass, real_par *ParPosX, real_par *ParPosY, real_par *ParPosZ,
                                   real_par *ParVelX, real_par *ParVelY, real_par *ParVelZ, real_par *ParTime,
                                   long_par *ParType, real_par *AllAttributeFlt[PAR_NATT_FLT_TOTAL],
                                   long_par *AllAttributeInt[PAR_NATT_INT_TOTAL] );
#endif

// external potential routines
void Init_ExtPot_ParEqmIC();




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

#  ifdef COMOVING
   Aux_Error( ERROR_INFO, "COMOVING must be disabled !!\n" );
#  endif

#  ifdef PARTICLE
   if ( OPT__INIT == INIT_BY_FUNCTION  &&  amr->Par->Init != PAR_INIT_BY_FUNCTION )
      Aux_Error( ERROR_INFO, "please set PAR_INIT = 1 (by FUNCTION) !!\n" );
#  else
   Aux_Error( ERROR_INFO, "PARTICLE must be enabled !!\n" );
#  endif

#  ifndef SUPPORT_GSL
   Aux_Error( ERROR_INFO, "SUPPORT_GSL must be enabled !!\n" );
#  endif


// warnings
   if ( MPI_Rank == 0 )
   {
      for (int f=0; f<6; f++)
      if ( OPT__BC_FLU[f] == BC_FLU_PERIODIC )
         Aux_Message( stderr, "WARNING : periodic BC for fluid is not recommended for this test !!\n" );

#     ifdef GRAVITY
      if ( OPT__BC_POT == BC_POT_PERIODIC )
         Aux_Message( stderr, "WARNING : periodic BC for gravity is not recommended for this test !!\n" );
#     endif
   } // if ( MPI_Rank == 0 )


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
//     --> a helper macro PRINT_RESET_PARA is defined in Macro.h
   const long   End_Step_Default = __INT_MAX__;
   const double End_T_Default    =  10;

   if ( END_STEP < 0 ) {
      END_STEP = End_Step_Default;
      PRINT_RESET_PARA( END_STEP, FORMAT_LONG, "" );
   }

   if ( END_T < 0.0 ) {
      END_T = End_T_Default;
      PRINT_RESET_PARA( END_T, FORMAT_REAL, "" );
   }

// load run-time parameters
   const char* FileName = "Input__TestProb";
   ReadPara_t *ReadPara  = new ReadPara_t;
   // ********************************************************************************************************************************
   // ReadPara->Add( "KEY_IN_THE_FILE",      &VARIABLE,              DEFAULT,       MIN,              MAX               );
   // ********************************************************************************************************************************
   ReadPara->Add( "ParEqmIC_SmallGas",       &ParEqmIC_SmallGas,     1e-3,          0.,               NoMax_double      );
   ReadPara->Read( FileName );
   delete ReadPara;

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

   fluid[DENS] = ParEqmIC_SmallGas;
   fluid[MOMX] = 0;
   fluid[MOMY] = 0;
   fluid[MOMZ] = 0;
#  ifdef GRAVITY
   fluid[ENGY] = ParEqmIC_SmallGas;
#  endif

// just set all passive scalars as zero
   for (int v=NCOMP_FLUID; v<NCOMP_TOTAL; v++)  fluid[v] = 0.0;

} // FUNCTION : SetGridIC
#endif // #if ( MODEL == HYDRO )



//-------------------------------------------------------------------------------------------------------
// Function    :  Init_TestProb_Hydro_ParEqmIC
// Description :  Test problem initializer
//
// Note        :  None
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Init_TestProb_Hydro_ParEqmIC()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );


// validate the compilation flags and runtime parameters
   Validate();

#  if ( MODEL == HYDRO )
// set the problem-specific runtime parameters
   SetParameter();


   Init_Function_User_Ptr  = SetGridIC;
#  ifdef MASSIVE_PARTICLES
   Par_Init_ByFunction_Ptr = Par_Init_ByFunction_ParEqmIC;
#  endif
#  ifdef GRAVITY
   if ( OPT__EXT_POT == EXT_POT_FUNC )
   Init_ExtPot_Ptr         = Init_ExtPot_ParEqmIC;
#  endif
#  endif // #if ( MODEL == HYDRO )


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : Init_TestProb_Hydro_ParEqmIC
