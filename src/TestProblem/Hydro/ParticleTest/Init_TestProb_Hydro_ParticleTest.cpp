#include "GAMER.h"



// problem-specific global variables
// =======================================================================================
static double ParTest_Dens_Bg;        // background mass density
static double ParTest_Pres_Bg;        // background pressure
static double ParTest_Ang_Freq;       // gas angular frequency
       int    ParTest_NPar[3];        // particles on a side
       double ParTest_Point_Mass;     // the mass of the active particles
       double ParTest_Par_Sep;        // the separation between the active particles
       bool   ParTest_Use_Tracers;    // whether or not to include tracers
       bool   ParTest_Use_Massive;    // whether or not to include massive particles

// =======================================================================================

// problem-specific function prototypes
#ifdef PARTICLE
void Par_Init_ByFunction_ParticleTest( const long NPar_ThisRank, const long NPar_AllRank,
                                       real_par *ParMass, real_par *ParPosX, real_par *ParPosY, real_par *ParPosZ,
                                       real_par *ParVelX, real_par *ParVelY, real_par *ParVelZ, real_par *ParTime,
                                       long_par *ParType, real_par *AllAttributeFlt[PAR_NATT_FLT_TOTAL],
                                       long_par *AllAttributeInt[PAR_NATT_INT_TOTAL] );
#endif

bool Flag_ParticleTest( const int i, const int j, const int k, const int lv,
                        const int PID, const double *Threshold );




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

#  ifndef PARTICLE
   Aux_Error( ERROR_INFO, "PARTICLE must be enabled !!\n" );
#  endif

#  ifdef COMOVING
   Aux_Error( ERROR_INFO, "COMOVING must be disabled !!\n" );
#  endif

   if ( !OPT__FREEZE_FLUID )
      Aux_Error( ERROR_INFO, "OPT__FREEZE_FLUID must be enabled !!\n" );


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Validating test problem %d ... done\n", TESTPROB_ID );

} // FUNCTION : Validate



#if ( MODEL == HYDRO )
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
// **************************************************************************************************************************
// LOAD_PARA( load_mode, "KEY_IN_THE_FILE",     &VARIABLE,               DEFAULT,      MIN,              MAX               );
// **************************************************************************************************************************
   LOAD_PARA( load_mode, "ParTest_Dens_Bg",     &ParTest_Dens_Bg,        1.0e-2,       Eps_double,       NoMax_double      );
   LOAD_PARA( load_mode, "ParTest_Pres_Bg",     &ParTest_Pres_Bg,        1.0e-2,       Eps_double,       NoMax_double      );
   LOAD_PARA( load_mode, "ParTest_Ang_Freq",    &ParTest_Ang_Freq,       0.00051668,   Eps_double,       1.0e-3            );
   LOAD_PARA( load_mode, "ParTest_NParX",       &ParTest_NPar[0],        32,           2,                128               );
   LOAD_PARA( load_mode, "ParTest_NParY",       &ParTest_NPar[1],        32,           2,                128               );
   LOAD_PARA( load_mode, "ParTest_NParZ",       &ParTest_NPar[2],        32,           2,                128               );
   LOAD_PARA( load_mode, "ParTest_Par_Sep",     &ParTest_Par_Sep,        0.5,          Eps_double,       NoMax_double      );
   LOAD_PARA( load_mode, "ParTest_Point_Mass",  &ParTest_Point_Mass,     1.0,          Eps_double,       NoMax_double      );
   LOAD_PARA( load_mode, "ParTest_Use_Tracers", &ParTest_Use_Tracers,    true,         Useless_bool,     Useless_bool      );
   LOAD_PARA( load_mode, "ParTest_Use_Massive", &ParTest_Use_Massive,    true,         Useless_bool,     Useless_bool      );

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


// (2) check the runtime parameters
   if ( !ParTest_Use_Tracers  &&  !ParTest_Use_Massive )
      Aux_Error( ERROR_INFO,
                 "either ParTest_Use_Tracer, ParTest_Use_Massive, or both must be true !!\n" );

#  ifndef TRACER
   if ( ParTest_Use_Tracers )    Aux_Error( ERROR_INFO, "must enable TRACER for ParTest_Use_Tracers !!\n" );
#  endif

#  ifndef GRAVITY
   if ( ParTest_Use_Massive )    Aux_Error( ERROR_INFO, "must enable GRAVITY for ParTest_Use_Massive !!\n" );
#  endif


// (3) reset other general-purpose parameters
//     --> a helper macro PRINT_RESET_PARA is defined in Macro.h
   const double End_T_Default    = 2*M_PI/ParTest_Ang_Freq;
   const long   End_Step_Default = __INT_MAX__;

   if ( END_STEP < 0 ) {
      END_STEP = End_Step_Default;
      PRINT_RESET_PARA( END_STEP, FORMAT_LONG, "" );
   }

   if ( END_T < 0.0 ) {
      END_T = End_T_Default;
      PRINT_RESET_PARA( END_T, FORMAT_REAL, "" );
   }

// overwrite the total number of particles
#  ifdef PARTICLE
   amr->Par->NPar_Active_AllRank = 0;
   if ( ParTest_Use_Massive )    amr->Par->NPar_Active_AllRank += 2;
   if ( ParTest_Use_Tracers )    amr->Par->NPar_Active_AllRank += ParTest_NPar[0]*ParTest_NPar[1]*ParTest_NPar[2];
   PRINT_RESET_PARA( amr->Par->NPar_Active_AllRank, FORMAT_LONG, "(PAR_NPAR in Input__Parameter)" );
#  endif


// (4) make a note
   if ( MPI_Rank == 0 )
   {
      Aux_Message( stdout, "=============================================================================\n" );
      Aux_Message( stdout, "  test problem ID            = %d\n",     TESTPROB_ID         );
      Aux_Message( stdout, "  background mass density    = %13.7e\n", ParTest_Dens_Bg     );
      Aux_Message( stdout, "  background pressure        = %13.7e\n", ParTest_Pres_Bg     );
      Aux_Message( stdout, "  angular frequency          = %13.7e\n", ParTest_Ang_Freq    );
      Aux_Message( stdout, "  number of particles (x)    = %d\n",     ParTest_NPar[0]     );
      Aux_Message( stdout, "                      (y)    = %d\n",     ParTest_NPar[1]     );
      Aux_Message( stdout, "                      (z)    = %d\n",     ParTest_NPar[2]     );
      Aux_Message( stdout, "  active particle separation = %13.7e\n", ParTest_Par_Sep     );
      Aux_Message( stdout, "  active particle mass       = %13.7e\n", ParTest_Point_Mass  );
      Aux_Message( stdout, "  include tracer particles   = %d\n",     ParTest_Use_Tracers );
      Aux_Message( stdout, "  include massive particles  = %d\n",     ParTest_Use_Massive );
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


   const double dr_Bg [2]  = { x - 0.25*amr->BoxSize[0], y - 0.25*amr->BoxSize[1] };
   const double dr_Mom[2]  = { x - 0.50*amr->BoxSize[0], y - 0.50*amr->BoxSize[1] };
   const double Radius_Bg  = sqrt( dr_Bg [0]*dr_Bg [0] + dr_Bg [1]*dr_Bg [1] );
   const double Radius_Mom = sqrt( dr_Mom[0]*dr_Mom[0] + dr_Mom[1]*dr_Mom[1] );

   const double Velocity = ParTest_Ang_Freq*Radius_Mom;

   const double Cos_theta = dr_Mom[0]/Radius_Mom;
   const double Sin_theta = dr_Mom[1]/Radius_Mom;

   double Dens, MomX, MomY, MomZ, Pres, Eint, Etot;

   Dens = ParTest_Dens_Bg * ( 1.0 + 5.0 * Radius_Bg / amr->BoxSize[0] );
   Pres = ParTest_Pres_Bg * ( 1.0 + 5.0 * Radius_Bg / amr->BoxSize[0] );
   MomX = -Dens*Velocity*Sin_theta;
   MomY = Dens*Velocity*Cos_theta;
   MomZ = 0.0;

   Eint = EoS_DensPres2Eint_CPUPtr( Dens, Pres, NULL, EoS_AuxArray_Flt,
                                    EoS_AuxArray_Int, h_EoS_Table );    // assuming EoS requires no passive scalars

   Etot = Hydro_ConEint2Etot( Dens, MomX, MomY, MomZ, Eint, 0.0 );      // do NOT include magnetic energy here

   fluid[DENS] = Dens;
   fluid[MOMX] = MomX;
   fluid[MOMY] = MomY;
   fluid[MOMZ] = MomZ;
   fluid[ENGY] = Etot;

} // FUNCTION : SetGridIC
#endif // #if ( MODEL == HYDRO )



//-------------------------------------------------------------------------------------------------------
// Function    :  Init_TestProb_Hydro_ParticleTest
// Description :  Test problem initializer
//
// Note        :  None
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Init_TestProb_Hydro_ParticleTest()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );


// validate the compilation flags and runtime parameters
   Validate();


#  if ( MODEL == HYDRO )
// set the problem-specific runtime parameters
   SetParameter();


// set the function pointers of various problem-specific routines
   Init_Function_User_Ptr        = SetGridIC;
   Output_User_Ptr               = NULL;
   Flag_User_Ptr                 = Flag_ParticleTest;
   Mis_GetTimeStep_User_Ptr      = NULL;
   Aux_Record_User_Ptr           = NULL;
   BC_User_Ptr                   = NULL;
   Flu_ResetByUser_Func_Ptr      = NULL;
   End_User_Ptr                  = NULL;
#  ifdef PARTICLE
   Par_Init_ByFunction_Ptr       = Par_Init_ByFunction_ParticleTest;
#  endif
#  ifdef SUPPORT_HDF5
   Output_HDF5_InputTest_Ptr     = LoadInputTestProb;
#  endif
#  endif // #if ( MODEL == HYDRO )


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : Init_TestProb_Hydro_ParticleTest



bool Flag_ParticleTest( const int i, const int j, const int k, const int lv,
                        const int PID, const double *Threshold )
{

   // Refine a rectanglar solid region in the center just to test behavior in
   // non-uniform regions

   const double dh     = amr->dh[lv];                                                  // grid size
   const double Pos[3] = { amr->patch[0][lv][PID]->EdgeL[0] + (i+0.5)*dh,              // x,y,z position
                           amr->patch[0][lv][PID]->EdgeL[1] + (j+0.5)*dh,
                           amr->patch[0][lv][PID]->EdgeL[2] + (k+0.5)*dh  };

   const double Center[3] = { 0.5*amr->BoxSize[0], 0.5*amr->BoxSize[1], 0.5*amr->BoxSize[2] };
   const double dr[3]     = { Pos[0]-Center[0], Pos[1]-Center[1], Pos[2]-Center[2] };

   bool Flag = (FABS(dr[0]) < Threshold[0]) && (FABS(dr[1]) < 2.0*Threshold[0]) && (FABS(dr[2]) < Threshold[0]);

   return Flag;

}
