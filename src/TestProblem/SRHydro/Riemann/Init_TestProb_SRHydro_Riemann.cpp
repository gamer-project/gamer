#include "GAMER.h"
#include "TestProb.h"


#if  ( MODEL == SR_HYDRO )

// problem-specific global variables
// =======================================================================================
typedef int Riemann_t;
const Riemann_t
   SOD_SHOCK_TUBE = 0,
   STRONG_SHOCK   = 1,
   TWO_SHOCKS     = 2,
   EINFELDT_1203  = 3,
   EINFELDT_1125  = 4,
   SONIC_RARE     = 5,
   USER_DEFINED   = 6;

static Riemann_t Riemann_Prob;         // target Riemann problem
static char      Riemann_Name[100];    // name of the target Riemann problem
static real      Riemann_RhoL;         // left-state density
static real      Riemann_VelL;         // left-state 4-velocity
static real      Riemann_VelL_T;       // left-state transverse 4-velocity
static real      Riemann_PreL;         // left-state pressure
static real      Riemann_RhoR;         // right-state density
static real      Riemann_VelR;         // right-state 4-velocity
static real      Riemann_VelR_T;       // right-state transverse 4-velocity
static real      Riemann_PreR;         // right-state pressure
static real      Riemann_EndT;         // end physical time
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

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Validating test problem %d ... done\n", TESTPROB_ID );

#  ifdef GRAVITY
   Aux_Error( ERROR_INFO, "GRAVITY must be disabled!!\n" );
#  endif

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

// (1-1) add parameters in the following format:
// --> note that VARIABLE, DEFAULT, MIN, and MAX must have the same data type
// --> some handy constants (e.g., NoMin_int, Eps_float, ...) are defined in "include/ReadPara.h"
// ********************************************************************************************************************************
// ReadPara->Add( "KEY_IN_THE_FILE",   &VARIABLE,              DEFAULT,       MIN,              MAX               );
// ********************************************************************************************************************************
   ReadPara->Add( "Riemann_Prob",      &Riemann_Prob,          -1,            0,                6                 );
   ReadPara->Add( "Riemann_LR",        &Riemann_LR,             1,            NoMin_int,        NoMax_int         );
   ReadPara->Add( "Riemann_XYZ",       &Riemann_XYZ,            0,            0,                2                 );
   ReadPara->Add( "Riemann_RhoL",      &Riemann_RhoL,           HUGE_NUMBER,  TINY_NUMBER,      HUGE_NUMBER       );
   ReadPara->Add( "Riemann_RhoR",      &Riemann_RhoR,           HUGE_NUMBER,  TINY_NUMBER,      HUGE_NUMBER       );
   ReadPara->Add( "Riemann_VelL",      &Riemann_VelL,           HUGE_NUMBER, -HUGE_NUMBER,      HUGE_NUMBER       );
   ReadPara->Add( "Riemann_VelR",      &Riemann_VelR,           HUGE_NUMBER, -HUGE_NUMBER,      HUGE_NUMBER       );
   ReadPara->Add( "Riemann_PreL",      &Riemann_PreL,           HUGE_NUMBER,  TINY_NUMBER,      HUGE_NUMBER       );
   ReadPara->Add( "Riemann_PreR",      &Riemann_PreR,           HUGE_NUMBER,  TINY_NUMBER,      HUGE_NUMBER       );
   ReadPara->Add( "Riemann_VelL_T",    &Riemann_VelL_T,         HUGE_NUMBER, -HUGE_NUMBER,      HUGE_NUMBER       );
   ReadPara->Add( "Riemann_VelR_T",    &Riemann_VelR_T,         HUGE_NUMBER, -HUGE_NUMBER,      HUGE_NUMBER       );
   ReadPara->Add( "Riemann_EndT",      &Riemann_EndT,           TINY_NUMBER,  TINY_NUMBER,      HUGE_NUMBER       );

   ReadPara->Read( FileName );

   delete ReadPara;

// (1-2) set the default values
   switch ( Riemann_Prob )
   {
      case SOD_SHOCK_TUBE : Riemann_RhoL =  5.0;  Riemann_VelL = 0.1;  Riemann_PreL = 10.0;  Riemann_VelL_T = 0.0;
                            Riemann_RhoR = 10.0;  Riemann_VelR = 0.1;  Riemann_PreR =  8.0;  Riemann_VelR_T = 0.0;
                            Riemann_EndT = 0.1;
                            sprintf( Riemann_Name, "Sod's shock tube" );
                            break;
      case STRONG_SHOCK   : Riemann_RhoL = 0.9   ;  Riemann_VelL = 0.1;  Riemann_PreL = 5.0e+10;  Riemann_VelL_T = 0.0;
                            Riemann_RhoR = 1.0;     Riemann_VelR = 0.1;  Riemann_PreR = 1.0e-5;  Riemann_VelR_T = 0.0;
                            Riemann_EndT = 0.1;
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

      case USER_DEFINED   : sprintf( Riemann_Name, "user defined" );
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
      Aux_Message( stdout, "  test problem ID                   = %d\n",     TESTPROB_ID );
      Aux_Message( stdout, "  target Riemann problem            = %s\n",     Riemann_Name );
      Aux_Message( stdout, "  left-state density                = %13.7e\n", Riemann_RhoL );
      Aux_Message( stdout, "  left-state 4-velocity             = %13.7e\n", Riemann_VelL );
      Aux_Message( stdout, "  left-state transverse 4-velocity  = %13.7e\n", Riemann_VelL_T );
      Aux_Message( stdout, "  left-state pressure               = %13.7e\n", Riemann_PreL );
      Aux_Message( stdout, "  right-state density               = %13.7e\n", Riemann_RhoR );
      Aux_Message( stdout, "  right-state 4-velocity            = %13.7e\n", Riemann_VelR );
      Aux_Message( stdout, "  right-state transverse 4-velocity = %13.7e\n", Riemann_VelR_T );
      Aux_Message( stdout, "  right-state pressure              = %13.7e\n", Riemann_PreR );
      Aux_Message( stdout, "  propagation direction             = %s%s\n",   ( Riemann_LR  >  0 ) ? "+" : "-",
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
//                3. Even when DUAL_ENERGY is adopted for SR_HYDRO, one does NOT need to set the dual-energy variable here
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
   double ConVarL[NCOMP_FLUID];
   double ConVarR[NCOMP_FLUID];
   double PriVarL[NCOMP_FLUID];
   double PriVarR[NCOMP_FLUID];

   PriVarL[0] = Riemann_RhoL;
   PriVarL[1] = Riemann_VelL;
   PriVarL[2] = Riemann_VelL_T;
   PriVarL[3] = 0.0;
   PriVarL[4] = Riemann_PreL;
   
   PriVarR[0] = Riemann_RhoR;
   PriVarR[1] = Riemann_VelR;
   PriVarR[2] = Riemann_VelR_T;
   PriVarR[3] = 0.0;
   PriVarR[4] = Riemann_PreR;

#  ifdef USE_3_VELOCITY
   SRHydro_4Velto3Vel( PriVarL, PriVarL );
   SRHydro_4Velto3Vel( PriVarR, PriVarR );
#  endif

// left-state
   SRHydro_Pri2Con(PriVarL, ConVarL, GAMMA);

// right-state
   SRHydro_Pri2Con(PriVarR, ConVarR, GAMMA);

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
      fluid[ TVar[0] ] = (real) ConVarL[0];
      fluid[ TVar[1] ] = (real) ConVarL[1];
      fluid[ TVar[2] ] = (real) ConVarL[2];
      fluid[ TVar[3] ] = (real) ConVarL[3];
      fluid[ TVar[4] ] = (real) ConVarL[4];
   }

   else
   {
      fluid[ TVar[0] ] = (real) ConVarR[0];
      fluid[ TVar[1] ] = (real) ConVarR[1];
      fluid[ TVar[2] ] = (real) ConVarR[2];
      fluid[ TVar[3] ] = (real) ConVarR[3];
      fluid[ TVar[4] ] = (real) ConVarR[4];
   }

   if ( Riemann_LR < 0 )
   {
      fluid[ TVar[1] ] = -fluid[ TVar[1] ];
      fluid[ TVar[2] ] = -fluid[ TVar[2] ];
   }


} // FUNCTION : SetGridIC


bool Flag_User_Riemann( const int i, const int j, const int k, const int lv, const int PID, const double Threshold )
{
   const double dh     = amr->dh[lv];                                                  // grid size
   const double Pos[3] = { amr->patch[0][lv][PID]->EdgeL[0] + (i+0.5)*dh,              // x,y,z position
                           amr->patch[0][lv][PID]->EdgeL[1] + (j+0.5)*dh,
                           amr->patch[0][lv][PID]->EdgeL[2] + (k+0.5)*dh  };


  
   if ( Step == 0 && 0.49*amr->BoxSize[0] < Pos[0] && Pos[0] < 0.51*amr->BoxSize[0] )
               return true;
   else        return false;

} // FUNCTION : Flag_User

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
void Init_TestProb_SRHydro_Riemann()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );


// validate the compilation flags and runtime parameters
   Validate();


// set the problem-specific runtime parameters
   SetParameter();


// set the function pointers of various problem-specific routines
   Init_Function_User_Ptr   = SetGridIC;
   Output_User_Ptr          = NULL;       // example: Hydro/AcousticWave/Init_TestProb_Hydro_AcousticWave.cpp --> OutputError()
   Flag_User_Ptr            = Flag_User_Riemann;  // example: AGORA_IsolatedGalaxy/Flag_AGORA.cpp
   Mis_GetTimeStep_User_Ptr = NULL;
   Aux_Record_User_Ptr      = NULL;
   BC_User_Ptr              = NULL;       // example: ELBDM/ExtPot/Init_TestProb_ELBDM_ExtPot.cpp --> BC()
   Flu_ResetByUser_Func_Ptr = NULL;
   End_User_Ptr             = NULL;       // example: Hydro/ClusterMerger_vs_Flash/Init_TestProb_Hydro_ClusterMerger_vs_Flash.cpp --> End_ClusterMerger()


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : Init_TestProb_Hydro_Riemann


#endif
