#include "GAMER.h"



// problem-specific global variables
// =======================================================================================
static double CR_Shocktube_Rho_R;       // density on the Right state
static double CR_Shocktube_Rho_L;       // density on the Left state
static double CR_Shocktube_Pres_R;      // gas pressure on the Right state
static double CR_Shocktube_Pres_L;      // gas pressure on the Left state
static double CR_Shocktube_PresCR_R;    // cosmic ray pressure on the Right state
static double CR_Shocktube_PresCR_L;    // cosmic ray pressure on the Left state
static int    CR_Shocktube_Dir;         // shock tube direction
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

#  ifndef COSMIC_RAY
   Aux_Error( ERROR_INFO, "COSMIC_RAY must be enabled !!\n" );
#  endif // #ifndef COSMIC_RAY

#  if ( EOS != EOS_COSMIC_RAY )
   Aux_Error( ERROR_INFO, "EOS != EOS_COSMIC_RAY when enable COSMIC_RAY !!\n" );
#  endif

#  ifdef CR_DIFFUSION
   Aux_Error( ERROR_INFO, "CR_DIFFUSION must be disabled !!\n" );
#  endif

// warnings


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Validating test problem %d ... done\n", TESTPROB_ID );

} // FUNCTION : Validate



#if ( MODEL == HYDRO  &&  defined COSMIC_RAY )
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
// LOAD_PARA( load_mode, "KEY_IN_THE_FILE",       &VARIABLE,                  DEFAULT,       MIN,              MAX               );
// ********************************************************************************************************************************
   LOAD_PARA( load_mode, "CR_Shocktube_Rho_R",    &CR_Shocktube_Rho_R,        0.0,           0.0,              NoMax_double      );
   LOAD_PARA( load_mode, "CR_Shocktube_Rho_L",    &CR_Shocktube_Rho_L,        0.0,           0.0,              NoMax_double      );
   LOAD_PARA( load_mode, "CR_Shocktube_Pres_R",   &CR_Shocktube_Pres_R,       0.0,           0.0,              NoMax_double      );
   LOAD_PARA( load_mode, "CR_Shocktube_Pres_L",   &CR_Shocktube_Pres_L,       0.0,           0.0,              NoMax_double      );
   LOAD_PARA( load_mode, "CR_Shocktube_PresCR_R", &CR_Shocktube_PresCR_R,     0.0,           0.0,              NoMax_double      );
   LOAD_PARA( load_mode, "CR_Shocktube_PresCR_L", &CR_Shocktube_PresCR_L,     0.0,           0.0,              NoMax_double      );
   LOAD_PARA( load_mode, "CR_Shocktube_Dir",      &CR_Shocktube_Dir,          0,             0,                2                 );

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

// (1-2) set the default values

// (1-3) check the runtime parameters


// (2) set the problem-specific derived parameters

// (3) reset other general-purpose parameters
//     --> a helper macro PRINT_RESET_PARA is defined in Macro.h
   const long   End_Step_Default = __INT_MAX__;
   const double End_T_Default    = __FLT_MAX__;

   if ( END_STEP < 0 ) {
      END_STEP = End_Step_Default;
      PRINT_RESET_PARA( END_STEP, FORMAT_LONG, "" );
   }

   if ( END_T < 0.0 ) {
      END_T = End_T_Default;
      PRINT_RESET_PARA( END_T, FORMAT_REAL, "" );
   }

   if (  ( CR_Shocktube_Dir == 0 && OPT__OUTPUT_PART != OUTPUT_X )  ||
         ( CR_Shocktube_Dir == 1 && OPT__OUTPUT_PART != OUTPUT_Y )  ||
         ( CR_Shocktube_Dir == 2 && OPT__OUTPUT_PART != OUTPUT_Z )    )
   {
      OPT__OUTPUT_PART = ( CR_Shocktube_Dir == 0 ) ? OUTPUT_X : ( CR_Shocktube_Dir == 1 ) ? OUTPUT_Y : OUTPUT_Z;
      PRINT_RESET_PARA( OPT__OUTPUT_PART, FORMAT_INT, "" );
   }


// (4) make a note
   if ( MPI_Rank == 0 )
   {
      Aux_Message( stdout, "=============================================================================\n" );
      Aux_Message( stdout, "  test problem ID       = %d\n",     TESTPROB_ID           );
      Aux_Message( stdout, "  CR_Shocktube_Rho_R    = %14.7e\n", CR_Shocktube_Rho_R    );
      Aux_Message( stdout, "  CR_Shocktube_RhoR_L   = %14.7e\n", CR_Shocktube_Rho_L    );
      Aux_Message( stdout, "  CR_Shocktube_Pres_R   = %14.7e\n", CR_Shocktube_Pres_R   );
      Aux_Message( stdout, "  CR_Shocktube_Pres_L   = %14.7e\n", CR_Shocktube_Pres_L   );
      Aux_Message( stdout, "  CR_Shocktube_PresCR_R = %14.7e\n", CR_Shocktube_PresCR_R );
      Aux_Message( stdout, "  CR_Shocktube_PresCR_L = %14.7e\n", CR_Shocktube_PresCR_L );
      Aux_Message( stdout, "  CR_Shocktube_Dir      = %d\n",     CR_Shocktube_Dir      );
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

   double Dens, MomX, MomY, MomZ, Pres, Eint, Etot, P_cr, CRay;
   double r;

   switch ( CR_Shocktube_Dir )
   {
      case 0: r=x; break;
      case 1: r=y; break;
      case 2: r=z; break;
   } // switch ( CR_Shocktube_Dir )

   if      ( r < amr->BoxCenter[CR_Shocktube_Dir] )
   {
      Dens = CR_Shocktube_Rho_L;
      Pres = CR_Shocktube_Pres_L;
      P_cr = CR_Shocktube_PresCR_L;
   }

   else if ( r > amr->BoxCenter[CR_Shocktube_Dir] )
   {
      Dens = CR_Shocktube_Rho_R;
      Pres = CR_Shocktube_Pres_R;
      P_cr = CR_Shocktube_PresCR_R;
   }

   else
   {
      Dens = 0.5*( CR_Shocktube_Rho_L + CR_Shocktube_Rho_R );
      Pres = 0.5*( CR_Shocktube_Pres_L + CR_Shocktube_Pres_R );
      P_cr = 0.5*( CR_Shocktube_PresCR_L + CR_Shocktube_PresCR_R );
   }

   MomX = 0.0;
   MomY = 0.0;
   MomZ = 0.0;

   const double GAMMA_CR_m1_inv = 1.0 / (GAMMA_CR - 1.0);
   Pres = Pres + P_cr;
   CRay = GAMMA_CR_m1_inv*P_cr;

// set the output array of passive scaler
   fluid[CRAY] = CRay;

   Eint = EoS_DensPres2Eint_CPUPtr( Dens, Pres, fluid+NCOMP_FLUID, EoS_AuxArray_Flt, EoS_AuxArray_Int, h_EoS_Table);
   Etot = Hydro_ConEint2Etot( Dens, MomX, MomY, MomZ, Eint, 0.0 );      // do NOT include magnetic energy here

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

   magnetic[MAGX] = 0.0;
   magnetic[MAGY] = 0.0;
   magnetic[MAGZ] = 0.0;

} // FUNCTION : SetBFieldIC
#endif // #ifdef MHD
#endif // #if ( MODEL == HYDRO  &&  defined COSMIC_RAY )



//-------------------------------------------------------------------------------------------------------
// Function    :  Init_TestProb_Hydro_CR_ShockTube
// Description :  Test problem initializer
//
// Note        :  None
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Init_TestProb_Hydro_CR_ShockTube()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );


// validate the compilation flags and runtime parameters
   Validate();


#  if ( MODEL == HYDRO  &&  defined COSMIC_RAY )
// set the problem-specific runtime parameters
   SetParameter();


// procedure to enable a problem-specific function:
   Init_Function_User_Ptr        = SetGridIC;
#  ifdef MHD
   Init_Function_BField_User_Ptr = SetBFieldIC;
#  endif
#  ifdef SUPPORT_HDF5
   Output_HDF5_InputTest_Ptr     = LoadInputTestProb;
#  endif
#  endif // #if ( MODEL == HYDRO  &&  defined COSMIC_RAY )


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : Init_TestProb_Hydro_CR_ShockTube
