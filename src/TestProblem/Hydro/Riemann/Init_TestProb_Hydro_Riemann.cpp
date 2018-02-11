#include "GAMER.h"
#include "TestProb.h"



// problem-specific global variables
// =======================================================================================
typedef int Riemann_t;
const Riemann_t
   SOD_SHOCK_TUBE = 0,
   STRONG_SHOCK   = 1,
   TWO_SHOCKS     = 2,
   EINFELDT_1203  = 3,
   EINFELDT_1125  = 4,
   SONIC_RARE     = 5;

static Riemann_t Riemann_Prob;         // target Riemann problem
static char      Riemann_Name[100];    // name of the target Riemann problem
static real      Riemann_RhoL;         // left-state density
static real      Riemann_VelL;         // left-state velocity
static real      Riemann_VelL_T;       // left-state transverse velocity
static real      Riemann_PreL;         // left-state pressure
static real      Riemann_RhoR;         // right-state density
static real      Riemann_VelR;         // right-state velocity
static real      Riemann_VelR_T;       // right-state transverse velocity
static real      Riemann_PreR;         // right-state pressure
static double    Riemann_EndT;         // end physical time
static int       Riemann_LR;           // wave propagation direction (>0/<0 --> positive/negative direction)
static int       Riemann_XYZ;          // wave propagation direction (0/1/2 --> x/y/z)
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

   if ( OPT__BC_FLU[0] != BC_FLU_OUTFLOW )
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

// (1-1) add parameters in the following format:
// --> note that VARIABLE, DEFAULT, MIN, and MAX must have the same data type
// --> some handy constants (e.g., NoMin_int, Eps_float, ...) are defined in "include/ReadPara.h"
// ********************************************************************************************************************************
// ReadPara->Add( "KEY_IN_THE_FILE",   &VARIABLE,              DEFAULT,       MIN,              MAX               );
// ********************************************************************************************************************************
   ReadPara->Add( "Riemann_Prob",      &Riemann_Prob,          -1,            0,                5                 );
   ReadPara->Add( "Riemann_LR",        &Riemann_LR,             1,            NoMin_int,        NoMax_int         );
   ReadPara->Add( "Riemann_XYZ",       &Riemann_XYZ,            0,            0,                2                 );

   ReadPara->Read( FileName );

   delete ReadPara;

// (1-2) set the default values
   switch ( Riemann_Prob )
   {
      case SOD_SHOCK_TUBE : Riemann_RhoL = 1.0;    Riemann_VelL = 0.0;  Riemann_PreL = 1.0;  Riemann_VelL_T = 0.0;
                            Riemann_RhoR = 0.125;  Riemann_VelR = 0.0;  Riemann_PreR = 0.1;  Riemann_VelR_T = 0.0;
                            Riemann_EndT = 0.1;
                            sprintf( Riemann_Name, "Sod's shock tube" );
                            break;

      case STRONG_SHOCK   : Riemann_RhoL = 1250.0;  Riemann_VelL = 0.0;  Riemann_PreL = 500.0;  Riemann_VelL_T = 0.0;
                            Riemann_RhoR =  125.0;  Riemann_VelR = 0.0;  Riemann_PreR =   5.0;  Riemann_VelR_T = 0.0;
                            Riemann_EndT = 0.4;
                            sprintf( Riemann_Name, "strong shock" );
                            break;

      case TWO_SHOCKS     : Riemann_RhoL = 1.0;  Riemann_VelL = 3.0;  Riemann_PreL = 1.0;  Riemann_VelL_T = 0.0;
                            Riemann_RhoR = 2.0;  Riemann_VelR = 1.0;  Riemann_PreR = 1.0;  Riemann_VelR_T = 0.0;
                            Riemann_EndT = 0.1;
                            sprintf( Riemann_Name, "two shocks" );
                            break;

      case EINFELDT_1203  : Riemann_RhoL = 1.0;  Riemann_VelL = -2.0;  Riemann_PreL = GAMMA-1.0;  Riemann_VelL_T = 0.0;
                            Riemann_RhoR = 1.0;  Riemann_VelR = +2.0;  Riemann_PreR = GAMMA-1.0;  Riemann_VelR_T = 0.0;
                            Riemann_EndT = 0.1;
                            sprintf( Riemann_Name, "Einfeldt's 1-2-0-3" );
                            break;

      case EINFELDT_1125  : Riemann_RhoL = 1.0;  Riemann_VelL = -1.0;  Riemann_PreL = 2.5*(GAMMA-1.0);  Riemann_VelL_T = -2.0;
                            Riemann_RhoR = 1.0;  Riemann_VelR = +1.0;  Riemann_PreR = 2.5*(GAMMA-1.0);  Riemann_VelR_T = +2.0;
                            Riemann_EndT = 0.1;
                            sprintf( Riemann_Name, "Einfeldt's 1-1-2-5" );
                            break;

      case SONIC_RARE     : Riemann_RhoL = 1.0;    Riemann_VelL = 0.75;  Riemann_PreL = 1.0;  Riemann_VelL_T = 0.0;
                            Riemann_RhoR = 0.125;  Riemann_VelR = 0.0;   Riemann_PreR = 0.1;  Riemann_VelR_T = 0.0;
                            Riemann_EndT = 0.1;
                            sprintf( Riemann_Name, "sonic rarefaction wave" );
                            break;

      default : Aux_Error( ERROR_INFO, "unsupported Riemann problem (%d) !!\n", Riemann_Prob );
   } // switch ( Riemann_Prob )

// (1-3) check the runtime parameters


// (2) set the problem-specific derived parameters


// (3) reset other general-purpose parameters
//     --> a helper macro PRINT_WARNING is defined in TestProb.h
   const long   End_Step_Default = __INT_MAX__;
   const double End_T_Default    = Riemann_EndT;

   if ( END_STEP < 0 ) {
      END_STEP = End_Step_Default;
      PRINT_WARNING( "END_STEP", END_STEP, FORMAT_LONG );
   }

   if ( END_T < 0.0 ) {
      END_T = End_T_Default;
      PRINT_WARNING( "END_T", END_T, FORMAT_REAL );
   }

   if (  ( Riemann_XYZ == 0 && OPT__OUTPUT_PART != OUTPUT_X )  ||
         ( Riemann_XYZ == 1 && OPT__OUTPUT_PART != OUTPUT_Y )  ||
         ( Riemann_XYZ == 2 && OPT__OUTPUT_PART != OUTPUT_Z )    )
   {
      OPT__OUTPUT_PART = ( Riemann_XYZ == 0 ) ? OUTPUT_X : ( Riemann_XYZ == 1 ) ? OUTPUT_Y : OUTPUT_Z;
      PRINT_WARNING( "OPT__OUTPUT_PART", OPT__OUTPUT_PART, FORMAT_INT );
   }


// (4) make a note
   if ( MPI_Rank == 0 )
   {
      Aux_Message( stdout, "=============================================================================\n" );
      Aux_Message( stdout, "  test problem ID                 = %d\n",     TESTPROB_ID );
      Aux_Message( stdout, "  target Riemann problem          = %s\n",     Riemann_Name );
      Aux_Message( stdout, "  left-state density              = %13.7e\n", Riemann_RhoL );
      Aux_Message( stdout, "  left-state velocity             = %13.7e\n", Riemann_VelL );
      Aux_Message( stdout, "  left-state transverse velocity  = %13.7e\n", Riemann_VelL_T );
      Aux_Message( stdout, "  left-state pressure             = %13.7e\n", Riemann_PreL );
      Aux_Message( stdout, "  right-state density             = %13.7e\n", Riemann_RhoR );
      Aux_Message( stdout, "  right-state velocity            = %13.7e\n", Riemann_VelR );
      Aux_Message( stdout, "  right-state transverse velocity = %13.7e\n", Riemann_VelR_T );
      Aux_Message( stdout, "  right-state pressure            = %13.7e\n", Riemann_PreR );
      Aux_Message( stdout, "  propagation direction           = %s%s\n",   ( Riemann_LR > 0 ) ? "+" : "-",
                                                                           ( Riemann_XYZ == 0 ) ? "x" :
                                                                           ( Riemann_XYZ == 1 ) ? "y" : "z" );
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

   const double _Gamma_m1 = 1.0/(GAMMA-1.0);

   double r, BoxCen;
   int    TVar[NCOMP_FLUID];

   switch ( Riemann_XYZ )
   {
      case 0 : r=x; TVar[0]=DENS; TVar[1]=MOMX; TVar[2]=MOMY; TVar[3]=MOMZ; TVar[4]=ENGY; BoxCen=0.5*amr->BoxSize[0]; break;
      case 1 : r=y; TVar[0]=DENS; TVar[1]=MOMY; TVar[2]=MOMZ; TVar[3]=MOMX; TVar[4]=ENGY; BoxCen=0.5*amr->BoxSize[1]; break;
      case 2 : r=z; TVar[0]=DENS; TVar[1]=MOMZ; TVar[2]=MOMX; TVar[3]=MOMY; TVar[4]=ENGY; BoxCen=0.5*amr->BoxSize[2]; break;
      default : Aux_Error( ERROR_INFO, "incorrect parameter %s = %d !!\n", "Riemann_XYZ", Riemann_XYZ );
   }

   if (  ( Riemann_LR > 0 && r < BoxCen )  ||  ( Riemann_LR < 0 && r > BoxCen )  )
   {
      fluid[ TVar[0] ] = Riemann_RhoL;
      fluid[ TVar[1] ] = Riemann_RhoL*Riemann_VelL;
      fluid[ TVar[2] ] = Riemann_RhoL*Riemann_VelL_T;
      fluid[ TVar[3] ] = 0.0;
      fluid[ TVar[4] ] = 0.5*( SQR(fluid[MOMX]) + SQR(fluid[MOMY]) + SQR(fluid[MOMZ]) )/fluid[DENS] + Riemann_PreL*_Gamma_m1;
   }

   else
   {
      fluid[ TVar[0] ] = Riemann_RhoR;
      fluid[ TVar[1] ] = Riemann_RhoR*Riemann_VelR;
      fluid[ TVar[2] ] = Riemann_RhoR*Riemann_VelR_T;
      fluid[ TVar[3] ] = 0.0;
      fluid[ TVar[4] ] = 0.5*( SQR(fluid[MOMX]) + SQR(fluid[MOMY]) + SQR(fluid[MOMZ]) )/fluid[DENS] + Riemann_PreR*_Gamma_m1;
   }

   if ( Riemann_LR < 0 )
   {
      fluid[ TVar[1] ] = -fluid[ TVar[1] ];
      fluid[ TVar[2] ] = -fluid[ TVar[2] ];
   }

} // FUNCTION : SetGridIC
#endif // #if ( MODEL == HYDRO )



//-------------------------------------------------------------------------------------------------------
// Function    :  Init_TestProb_Hydro_Riemann
// Description :  Test problem initializer
//
// Note        :  None
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Init_TestProb_Hydro_Riemann()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );


// validate the compilation flags and runtime parameters
   Validate();


#  if ( MODEL == HYDRO )
// set the problem-specific runtime parameters
   SetParameter();


// set the function pointers of various problem-specific routines
   Init_Function_User_Ptr   = SetGridIC;
   Output_User_Ptr          = NULL;       // example: Hydro/AcousticWave/Init_TestProb_Hydro_AcousticWave.cpp --> OutputError()
   Flag_User_Ptr            = NULL;       // example: AGORA_IsolatedGalaxy/Flag_AGORA.cpp
   Mis_GetTimeStep_User_Ptr = NULL;
   Aux_Record_User_Ptr      = NULL;
   BC_User_Ptr              = NULL;       // example: ELBDM/ExtPot/Init_TestProb_ELBDM_ExtPot.cpp --> BC()
   Flu_ResetByUser_Func_Ptr = NULL;
   End_User_Ptr             = NULL;       // example: Hydro/ClusterMerger_vs_Flash/Init_TestProb_Hydro_ClusterMerger_vs_Flash.cpp --> End_ClusterMerger()
#  endif // #if ( MODEL == HYDRO )


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : Init_TestProb_Hydro_Riemann
