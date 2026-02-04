#include "GAMER.h"



// problem-specific global variables
// =======================================================================================
static double StarFormationThreshold_MassDensity_Min;    // Minimum mass density in the box
static double StarFormationThreshold_MassDensity_Max;    // Maximum mass density in the box
static double StarFormationThreshold_Temperature_Min;    // Minimum temperature in the box
static double StarFormationThreshold_Temperature_Max;    // Maximum temperature in the box

static double StarFormationThreshold_logDens_Min;        // Minimum log( mass density ) in the box
static double StarFormationThreshold_logDens_Max;        // Maximum log( mass density ) in the box
static double StarFormationThreshold_logDens_Range;      // Range of log ( mass density )
static double StarFormationThreshold_logTemp_Min;        // Minimum log( temperature ) in the box
static double StarFormationThreshold_logTemp_Max;        // Maximum log( temperature ) in the box
static double StarFormationThreshold_logTemp_Range;      // Range of log ( temperature )
static double StarFormationThreshold_FreeFallTime_Min;   // Minimum free-fall time in the box
static double StarFormationThreshold_FreeFallTime_Max;   // Maximum free-fall time in the box
// =======================================================================================


// problem-specific function prototypes
#ifdef PARTICLE
void Par_Init_ByFunction_StarFormationThreshold( const long NPar_ThisRank, const long NPar_AllRank,
                                                 real_par *ParMass, real_par *ParPosX, real_par *ParPosY, real_par *ParPosZ,
                                                 real_par *ParVelX, real_par *ParVelY, real_par *ParVelZ, real_par *ParTime,
                                                 long_par *ParType, real_par *AllAttributeFlt[PAR_NATT_FLT_TOTAL],
                                                 long_par *AllAttributeInt[PAR_NATT_INT_TOTAL] );
#endif




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

#  ifndef PARTICLE
   Aux_Error( ERROR_INFO, "PARTICLE must be enabled !!\n" );
#  endif

#  ifndef GRAVITY
   Aux_Error( ERROR_INFO, "GRAVITY must be enabled !!\n" );
#  endif

#  ifndef STAR_FORMATION
   Aux_Error( ERROR_INFO, "STAR_FORMATION must be enabled !!\n" );
#  endif

#  ifdef COMOVING
   Aux_Error( ERROR_INFO, "COMOVING must be disabled !!\n" );
#  endif

   if ( !OPT__UNIT )
      Aux_Error( ERROR_INFO, "OPT__UNIT must be enabled !!\n" );

   if ( !OPT__FREEZE_FLUID )
      Aux_Error( ERROR_INFO, "OPT__FREEZE_FLUID must be enabled !!\n" );

   for (int f=0; f<6; f++)
   if ( OPT__BC_FLU[f] == BC_FLU_PERIODIC )
      Aux_Error( ERROR_INFO, "do not use periodic BC (OPT__BC_FLU* = 1) for this test !!\n" );

#  ifdef GRAVITY
   if ( OPT__BC_POT == BC_POT_PERIODIC )
      Aux_Error( ERROR_INFO, "do not use periodic BC (OPT__BC_POT = 1) for this test !!\n" );
#  endif

#  ifdef PARTICLE
   if ( OPT__INIT == INIT_BY_FUNCTION  &&  amr->Par->Init != PAR_INIT_BY_FUNCTION )
      Aux_Error( ERROR_INFO, "please set PAR_INIT = 1 (by FUNCTION) !!\n" );

   if ( !OPT__FREEZE_PAR )
      Aux_Error( ERROR_INFO, "OPT__FREEZE_PAR must be enabled !!\n" );
#  endif


// warnings
   if ( MPI_Rank == 0 )
   {

      if ( !OPT__RESET_FLUID )
         Aux_Message( stderr, "WARNING : it's recommended to enable OPT__RESET_FLUID to have a constant environment for this test !!\n" );

   } // if ( MPI_Rank == 0 )


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Validating test problem %d ... done\n", TESTPROB_ID );

} // FUNCTION : Validate



#if ( MODEL == HYDRO  &&  defined MASSIVE_PARTICLES )
//-------------------------------------------------------------------------------------------------------
// Function    :  LoadInputTestProb
// Description :  Read problem-specific runtime parameters from Input__TestProb and store them in HDF5 snapshots (Data_*)
//
// Note        :  1. Invoked by SetParameter() to read parameters
//                2. Invoked by Output_DumpData_Total_HDF5() using the function pointer Output_HDF5_InputTest_Ptr to store parameters
//                3. If there is no problem-specific runtime parameter to load, add at least one parameter
//                   to prevent an empty structure in HDF5_Output_t
//                   --> Example:
//                       LOAD_PARA( load_mode, "TestProb_ID", &TESTPROB_ID, TESTPROB_ID, TESTPROB_ID, TESTPROB_ID );
//
// Parameter   :  load_mode      : Mode for loading parameters
//                                 --> LOAD_READPARA    : Read parameters from Input__TestProb
//                                     LOAD_HDF5_OUTPUT : Store parameters in HDF5 snapshots
//                ReadPara       : Data structure for reading parameters (used with LOAD_READPARA)
//                HDF5_InputTest : Data structure for storing parameters in HDF5 snapshots (used with LOAD_HDF5_OUTPUT)
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void LoadInputTestProb( const LoadParaMode_t load_mode, ReadPara_t *ReadPara, HDF5_Output_t *HDF5_InputTest )
{

#  ifndef SUPPORT_HDF5
   if ( load_mode == LOAD_HDF5_OUTPUT )   Aux_Error( ERROR_INFO, "please turn on SUPPORT_HDF5 in the Makefile for load_mode == LOAD_HDF5_OUTPUT !!\n" );
#  endif

   if ( load_mode == LOAD_READPARA     &&  ReadPara       == NULL )   Aux_Error( ERROR_INFO, "load_mode == LOAD_READPARA and ReadPara == NULL !!\n" );
   if ( load_mode == LOAD_HDF5_OUTPUT  &&  HDF5_InputTest == NULL )   Aux_Error( ERROR_INFO, "load_mode == LOAD_HDF5_OUTPUT and HDF5_InputTest == NULL !!\n" );

// add parameters in the following format:
// --> note that VARIABLE, DEFAULT, MIN, and MAX must have the same data type
// --> some handy constants (e.g., NoMin_int, Eps_float, ...) are defined in "include/ReadPara.h"
// --> LOAD_PARA() is defined in "include/TestProb.h"
// ********************************************************************************************************************************
// LOAD_PARA( load_mode, "KEY_IN_THE_FILE",                        &VARIABLE,                                     DEFAULT,       MIN,              MAX               );
// ********************************************************************************************************************************
   LOAD_PARA( load_mode, "StarFormationThreshold_MassDensity_Min", &StarFormationThreshold_MassDensity_Min,       1.0e-29,       Eps_double,       NoMax_double      );
   LOAD_PARA( load_mode, "StarFormationThreshold_MassDensity_Max", &StarFormationThreshold_MassDensity_Max,       1.0e-21,       Eps_double,       NoMax_double      );
   LOAD_PARA( load_mode, "StarFormationThreshold_Temperature_Min", &StarFormationThreshold_Temperature_Min,       1.0e+00,       Eps_double,       NoMax_double      );
   LOAD_PARA( load_mode, "StarFormationThreshold_Temperature_Max", &StarFormationThreshold_Temperature_Max,       1.0e+08,       Eps_double,       NoMax_double      );

} // FUNCITON : LoadInputTestProb



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
//                3. Must call EoS_Init() before calling any other EoS routine
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void SetParameter()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Setting runtime parameters ...\n" );


// (1) load the problem-specific runtime parameters
// (1-1) read parameters from Input__TestProb
   const char FileName[] = "Input__TestProb";
   ReadPara_t *ReadPara  = new ReadPara_t;

   LoadInputTestProb( LOAD_READPARA, ReadPara, NULL );

   ReadPara->Read( FileName );

   delete ReadPara;

// (1-2) set the default values

// (1-3) check the runtime parameters
   if ( StarFormationThreshold_MassDensity_Min > StarFormationThreshold_MassDensity_Max )
      Aux_Error( ERROR_INFO, "MassDensity_Min = %14.7e > MassDensity_Max = %14.7e !!\n",
                 StarFormationThreshold_MassDensity_Min, StarFormationThreshold_MassDensity_Max );

   if ( StarFormationThreshold_Temperature_Min > StarFormationThreshold_Temperature_Max )
      Aux_Error( ERROR_INFO, "Temperature_Min = %14.7e > Temperature_Max = %14.7e !!\n",
                 StarFormationThreshold_Temperature_Min, StarFormationThreshold_Temperature_Max );


// (2) set the problem-specific derived parameters
// convert to code units
   StarFormationThreshold_MassDensity_Min /= UNIT_D;
   StarFormationThreshold_MassDensity_Max /= UNIT_D;

   StarFormationThreshold_logDens_Min   = log10( StarFormationThreshold_MassDensity_Min );
   StarFormationThreshold_logDens_Max   = log10( StarFormationThreshold_MassDensity_Max );
   StarFormationThreshold_logDens_Range = StarFormationThreshold_logDens_Max - StarFormationThreshold_logDens_Min;
   StarFormationThreshold_logTemp_Min   = log10( StarFormationThreshold_Temperature_Min );
   StarFormationThreshold_logTemp_Max   = log10( StarFormationThreshold_Temperature_Max );
   StarFormationThreshold_logTemp_Range = StarFormationThreshold_logTemp_Max - StarFormationThreshold_logTemp_Min;

   StarFormationThreshold_FreeFallTime_Max = sqrt( (3.0*M_PI) / (32.0*NEWTON_G*StarFormationThreshold_MassDensity_Min) );
   StarFormationThreshold_FreeFallTime_Min = sqrt( (3.0*M_PI) / (32.0*NEWTON_G*StarFormationThreshold_MassDensity_Max) );


// (3) reset other general-purpose parameters
//     --> a helper macro PRINT_RESET_PARA is defined in Macro.h
   const long   End_Step_Default = __INT_MAX__;
   const double End_T_Default    = 1.0*sqrt( StarFormationThreshold_FreeFallTime_Max * StarFormationThreshold_FreeFallTime_Min );

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
      Aux_Message( stdout, "  test problem ID                             = %d\n",             TESTPROB_ID                                              );
      Aux_Message( stdout, "  StarFormationThreshold_MassDensity_Min      = %13.7e UNIT_D\n",  StarFormationThreshold_MassDensity_Min                   );
      Aux_Message( stdout, "                                              = %13.7e g/cm^3\n",  StarFormationThreshold_MassDensity_Min*UNIT_D            );
      Aux_Message( stdout, "                                              = %13.7e mH/cm^3\n", StarFormationThreshold_MassDensity_Min*UNIT_D/Const_mH   );
      Aux_Message( stdout, "                           -> Free-fall time  = %13.7e UNIT_T\n",  StarFormationThreshold_FreeFallTime_Max                  );
      Aux_Message( stdout, "                                              = %13.7e Myr\n",     StarFormationThreshold_FreeFallTime_Max*UNIT_T/Const_Myr );
      Aux_Message( stdout, "  StarFormationThreshold_MassDensity_Max      = %13.7e UNIT_D\n",  StarFormationThreshold_MassDensity_Max                   );
      Aux_Message( stdout, "                                              = %13.7e g/cm^3\n",  StarFormationThreshold_MassDensity_Max*UNIT_D            );
      Aux_Message( stdout, "                                              = %13.7e mH/cm^3\n", StarFormationThreshold_MassDensity_Max*UNIT_D/Const_mH   );
      Aux_Message( stdout, "                           -> Free-fall time  = %13.7e UNIT_T\n",  StarFormationThreshold_FreeFallTime_Min                  );
      Aux_Message( stdout, "                                              = %13.7e Myr\n",     StarFormationThreshold_FreeFallTime_Min*UNIT_T/Const_Myr );
      Aux_Message( stdout, "  StarFormationThreshold_Temperature_Min      = %13.7e K\n",       StarFormationThreshold_Temperature_Min                   );
      Aux_Message( stdout, "  StarFormationThreshold_Temperature_Max      = %13.7e K\n",       StarFormationThreshold_Temperature_Max                   );
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
//                   (unless OPT__INIT_GRID_WITH_OMP is disabled)
//                   --> Please ensure that everything here is thread-safe
//                3. Even when DUAL_ENERGY is adopted for HYDRO, one does NOT need to set the dual-energy variable here
//                   --> It will be calculated automatically
//                4. For MHD, do NOT add magnetic energy (i.e., 0.5*B^2) to fluid[ENGY] here
//                   --> It will be added automatically later
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


// check
#  ifdef GAMER_DEBUG
   if ( EoS_DensTemp2Pres_CPUPtr == NULL )
      Aux_Error( ERROR_INFO, "EoS_DensTemp2Pres_CPUPtr == NULL !!\n" );
#  endif

// compute the gas log density     by linear interpolation in x-direction
   const double logDens = StarFormationThreshold_logDens_Min +
                          StarFormationThreshold_logDens_Range*( x / amr->BoxSize[0] );
   const double Dens    = pow( 10, logDens );

// compute the gas log temperature by linear interpolation in y-direction
   const double logTemp = StarFormationThreshold_logTemp_Min +
                          StarFormationThreshold_logTemp_Range*( y / amr->BoxSize[1] );
   const double Temp    = pow( 10, logTemp );

// compute the gas pressure
   const double Pres = EoS_DensTemp2Pres_CPUPtr( Dens, Temp, NULL, EoS_AuxArray_Flt, EoS_AuxArray_Int,
                                                 h_EoS_Table ); // assuming EoS requires no passive scalars

// assume no momentum
   const double MomX = 0.0;
   const double MomY = 0.0;
   const double MomZ = 0.0;

// compute the total gas energy
   const double Eint = EoS_DensPres2Eint_CPUPtr( Dens, Pres, NULL, EoS_AuxArray_Flt,
                                                 EoS_AuxArray_Int, h_EoS_Table );   // assuming EoS requires no passive scalars
   const double Etot = Hydro_ConEint2Etot( Dens, MomX, MomY, MomZ, Eint, 0.0 );     // do NOT include magnetic energy here


// set the output array
   fluid[DENS] = Dens;
   fluid[MOMX] = MomX;
   fluid[MOMY] = MomY;
   fluid[MOMZ] = MomZ;
   fluid[ENGY] = Etot;

} // FUNCTION : SetGridIC



//-------------------------------------------------------------------------------------------------------
// Function    :  Flu_ResetByUser_StarFormationThreshold
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
//                Emag     : Magnetic energy (MHD only)
//                x/y/z    : Target physical coordinates
//                Time     : Target physical time
//                dt       : Time interval to advance solution
//                lv       : Target refinement level
//                AuxArray : Auxiliary array
//
// Return      :  true  : This cell has been reset
//                false : This cell has not been reset
//-------------------------------------------------------------------------------------------------------
int Flu_ResetByUser_StarFormationThreshold( real fluid[], const double Emag, const double x, const double y, const double z, const double Time,
                                            const double dt, const int lv, double AuxArray[] )
{

   SetGridIC( fluid, x, y, z, Time, lv, AuxArray );

   return true;

} // FUNCTION : Flu_ResetByUser_StarFormationThreshold



//-------------------------------------------------------------------------------------------------------
// Function    :  Mis_GetTimeStep_StarFormationThreshold
// Description :  0.1 times minimum free-fall time
//
// Note        :  1. This function should be applied to both physical and comoving coordinates and always
//                   return the evolution time-step (dt) actually used in various solvers
//                   --> Physical coordinates : dt = physical time interval
//                       Comoving coordinates : dt = delta(scale_factor) / ( Hubble_parameter*scale_factor^3 )
//                   --> We convert dt back to the physical time interval, which equals "delta(scale_factor)"
//                       in the comoving coordinates, in Mis_GetTimeStep()
//                2. Invoked by Mis_GetTimeStep() using the function pointer "Mis_GetTimeStep_User_Ptr",
//                   which must be set by a test problem initializer
//                3. Enabled by the runtime option "OPT__DT_USER"
//
// Parameter   :  lv       : Target refinement level
//                dTime_dt : dTime/dt (== 1.0 if COMOVING is off)
//
// Return      :  dt
//-------------------------------------------------------------------------------------------------------
double Mis_GetTimeStep_StarFormationThreshold( const int lv, const double dTime_dt )
{

   double dt_user = 0.1*StarFormationThreshold_FreeFallTime_Min;

   return dt_user;

} // FUNCTION : Mis_GetTimeStep_StarFormationThreshold
#endif // #if ( MODEL == HYDRO  &&  defined MASSIVE_PARTICLES )



//-------------------------------------------------------------------------------------------------------
// Function    :  Init_TestProb_Hydro_StarFormationThreshold
// Description :  Test problem initializer
//
// Note        :  None
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Init_TestProb_Hydro_StarFormationThreshold()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );


// validate the compilation flags and runtime parameters
   Validate();


#  if ( MODEL == HYDRO  &&  defined MASSIVE_PARTICLES )
// set the problem-specific runtime parameters
   SetParameter();


// set the function pointers of various problem-specific routines
   Init_Function_User_Ptr            = SetGridIC;
   Par_Init_ByFunction_Ptr           = Par_Init_ByFunction_StarFormationThreshold;
   Flu_ResetByUser_Func_Ptr          = Flu_ResetByUser_StarFormationThreshold;
   Mis_GetTimeStep_User_Ptr          = Mis_GetTimeStep_StarFormationThreshold;
#  ifdef SUPPORT_HDF5
   Output_HDF5_InputTest_Ptr         = LoadInputTestProb;
#  endif
#  endif // #if ( MODEL == HYDRO  &&  defined MASSIVE_PARTICLES )


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : Init_TestProb_Hydro_StarFormationThreshold
