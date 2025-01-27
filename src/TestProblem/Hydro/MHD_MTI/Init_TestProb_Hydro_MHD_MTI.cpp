#include "GAMER.h"



// problem-specific global variables
// =======================================================================================
static double MHD_MTI_Rho0;          // scale density
static double MHD_MTI_P0;            // scale pressure
static double MHD_MTI_v0;            // velocity perturbation amplitude
static double MHD_MTI_z0;            // scale height
static double MHD_MTI_B0;            // background magnetic field
static int    MHD_MTI_Dir;           // magnetic field direction (0/1/2) --> (x/y/z)
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

#  ifndef MHD
   Aux_Error( ERROR_INFO, "MHD must be enabled !!\n" );
#  endif

#  ifndef CONDUCTION
   Aux_Error( ERROR_INFO, "CONDUCTION must be enabled !!\n" );
#  endif

#  ifndef GRAVITY
   Aux_Error( ERROR_INFO, "GRAVITY must be enabled !!\n" );
#  endif

#  ifdef PARTICLE
   Aux_Error( ERROR_INFO, "PARTICLE must be disabled !!\n" );
#  endif

#  if ( EOS != EOS_GAMMA )
   Aux_Error( ERROR_INFO, "EOS != EOS_GAMMA !!\n" );
#  endif

   for (int f=0; f<4; f++)
   if ( OPT__BC_FLU[f] != BC_FLU_PERIODIC )
      Aux_Error( ERROR_INFO, "must adopt periodic BC (i.e., \"OPT__BC_FLU_* = 1\") for the x-y directions!!\n" );

// warnings
   if ( MPI_Rank == 0 )
   {
#     ifndef FLOAT8
      Aux_Message( stderr, "WARNING : it's recommended to enable FLOAT8 for this test !!\n" );
#     endif
   }

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Validating test problem %d ... done\n", TESTPROB_ID );

} // FUNCTION : Validate



#if ( MODEL == HYDRO  &&  defined GRAVITY )
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
   ReadPara->Add( "MHD_MTI_Rho0",    &MHD_MTI_Rho0,        -1.0,          Eps_double,       NoMax_double      );
   ReadPara->Add( "MHD_MTI_P0",      &MHD_MTI_P0,          -1.0,          Eps_double,       NoMax_double      );
   ReadPara->Add( "MHD_MTI_v0",      &MHD_MTI_v0,           0.0,          NoMin_double,     NoMax_double      );
   ReadPara->Add( "MHD_MTI_z0",      &MHD_MTI_z0,          -1.0,          0.0,              NoMax_double      );
   ReadPara->Add( "MHD_MTI_Dir",     &MHD_MTI_Dir,          0,            0,                2                 );
   ReadPara->Add( "MHD_MTI_B0",      &MHD_MTI_B0,          -1.0,          0.0,              NoMax_double      );

   ReadPara->Read( FileName );

   delete ReadPara;


// (2) set the problem-specific derived parameters
   MHD_MTI_WaveLength = ( MHD_MTI_Dir == 3 ) ? amr->BoxSize[0]/sqrt(3.0) : amr->BoxSize[MHD_MTI_Dir];

// assuming EOS_GAMMA
   if ( MHD_MTI_Mode == 1 )
      MHD_MTI_WaveSpeed = sqrt( GAMMA*MHD_MTI_P0/MHD_MTI_Rho0 + SQR(MHD_MTI_B0)/MHD_MTI_Rho0 );
   else
      Aux_Error( ERROR_INFO, "unsupported MHD_MTI_Mode = %d !!\n", MHD_MTI_Mode );


// (3) reset other general-purpose parameters
//     --> a helper macro PRINT_RESET_PARA is defined in Macro.h
   const double End_T_Default    = MHD_MTI_WaveLength / MHD_MTI_WaveSpeed;
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
      Aux_Message( stdout, "  test problem ID      = %d\n",      TESTPROB_ID        );
      Aux_Message( stdout, "  mode                 = %d\n",      MHD_MTI_Mode       );
      Aux_Message( stdout, "  background density   = % 14.7e\n", MHD_MTI_Rho0       );
      Aux_Message( stdout, "  background pressure  = % 14.7e\n", MHD_MTI_P0         );
      Aux_Message( stdout, "  velocity amplitude   = % 14.7e\n", MHD_MTI_v0         );
      Aux_Message( stdout, "  background B field   = % 14.7e\n", MHD_MTI_B0         );
      Aux_Message( stdout, "  B-field direction    = %d\n",      MHD_MTI_Dir        );
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

   double Dens, MomX, MomY, MomZ, Pres, Eint, Etot;

   Dens = MHD_MTI_Rho0 *  SQR( 1.0 - z / MHD_MTI_z0 );
   Pres = MHD_MTI_P0   * CUBE( 1.0 - z / MHD_MTI_z0 );
   MomX = 0.0;
   MomY = 0.0;
   MomZ = Dens * MHD_MTI_v0 * sin ( 2.0*M_PI*z / amr->BoxSize[2] );

// compute the total gas energy
   Eint = EoS_DensPres2Eint_CPUPtr( Dens, Pres, NULL, EoS_AuxArray_Flt,
                                    EoS_AuxArray_Int, h_EoS_Table );   // assuming EoS requires no passive scalars
   Etot = Hydro_ConEint2Etot( Dens, MomX, MomY, MomZ, Eint, 0.0 );     // do NOT include magnetic energy here

// set the output array
   fluid[DENS] = Dens;
   fluid[MOMX] = MomX;
   fluid[MOMY] = MomY;
   fluid[MOMZ] = MomZ;
   fluid[ENGY] = Etot;

} // FUNCTION : SetGridIC



#ifdef MHD
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

   magnetic[MAGX]        = 0.0;
   magnetic[MAGY]        = 0.0;
   magnetic[MAGZ]        = 0.0;
   magnetic[MHD_MTI_Dir] = MHD_MTI_B0;

} // FUNCTION : SetBFieldIC


#endif // #ifdef MHD

//-------------------------------------------------------------------------------------------------------
// Function    :  MTI_BC
// Description :  Set the external boundary condition to the IC condition
//
// Note        :  1. Linked to the function pointer "BC_User_Ptr"
//
// Parameter   :  Array          : Array to store the prepared data including ghost zones
//                ArraySize      : Size of Array including the ghost zones on each side
//                fluid          : Fluid fields to be set
//                NVar_Flu       : Number of fluid variables to be prepared
//                GhostSize      : Number of ghost zones
//                idx            : Array indices
//                pos            : Physical coordinates
//                Time           : Physical time
//                lv             : Refinement level
//                TFluVarIdxList : List recording the target fluid variable indices ( = [0 ... NCOMP_TOTAL-1] )
//                AuxArray       : Auxiliary array
//
// Return      :  fluid
//-------------------------------------------------------------------------------------------------------
void MTI_BC( real Array[], const int ArraySize[], real fluid[], const int NVar_Flu,
             const int GhostSize, const int idx[], const double pos[], const double Time,
             const int lv, const int TFluVarIdxList[], double AuxArray[] )
{

   const real g0 = -3.0*MHD_MTI_P0/MHD_MTI_Rho0/MHD_MTI_z0;
   const real T0 = MHD_MTI_P0 / MHD_MTI_Rho0;
   const real zbnd = ;
   const real Temp = T0 * ( 1.0 - zbnd / MHD_MTI_z0 );
   const real Dens = ;
   const real Pres = Dens * Temp;
   const real MomX = 0.0;
   const real MomY = 0.0;
   const real MomZ = 0.0;

// compute the total gas energy
   Eint = EoS_DensPres2Eint_CPUPtr( Dens, Pres, NULL, EoS_AuxArray_Flt,
                                    EoS_AuxArray_Int, h_EoS_Table );   // assuming EoS requires no passive scalars
   Etot = Hydro_ConEint2Etot( Dens, MomX, MomY, MomZ, Eint, 0.0 );     // do NOT include magnetic energy here

// set the output array
   fluid[DENS] = Dens;
   fluid[MOMX] = MomX;
   fluid[MOMY] = MomY;
   fluid[MOMZ] = MomZ;
   fluid[ENGY] = Etot;

} // FUNCTION : MTI_BC

#endif // #if ( MODEL == HYDRO  &&  defined GRAVITY )

//-------------------------------------------------------------------------------------------------------
// Function    :  Init_TestProb_Hydro_MHD_MTI
// Description :  Test problem initializer
//
// Note        :  None
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Init_TestProb_Hydro_MHD_MTI()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );


// validate the compilation flags and runtime parameters
   Validate();


#  if ( MODEL == HYDRO )
// set the problem-specific runtime parameters
   SetParameter();


// set the function pointers of various problem-specific routines
   Init_Function_User_Ptr        = SetGridIC;
   BC_User_Ptr                   = MTI_BC;
#  ifdef MHD
   Init_Function_BField_User_Ptr = SetBFieldIC;
#  endif
#  ifdef GRAVITY
   Init_ExtAcc_Ptr               = Init_ExtAcc_MTI;
#  endif
#  endif // #if ( MODEL == HYDRO )


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : Init_TestProb_Hydro_MHD_MTI
