#include "GAMER.h"
#include "TestProb.h"



// problem-specific global variables
// =======================================================================================
       double Bondi_MassBH;         // black hole mass
static double Bondi_Rho0;           // background density
static double Bondi_T0;             // background temperature
       double Bondi_RefineRadius0;  // refinement radius at the base level
                                    // NOTE: refinement radius at Lv is set to Bondi_RefineRadius0*2^(-Lv)
                                    // --> all refinement shells have roughly the same number of cells at its level
                                    //     (except Lv=MAX_LEVEL, which will have twice the number of cells along the radius
                                    //     unless Bondi_HalfMaxLvRefR is on)
       bool   Bondi_HalfMaxLvRefR;  // halve the refinement radius at the maximum level
       double Bondi_InBC_Rho;       // density     inside the void region
static double Bondi_InBC_T;         // temperature inside the void region
static double Bondi_InBC_NCell;     // number of finest cells for the radius of the void region
static double Bondi_Soften_NCell;   // number of finest cells for the soften length (<=0.0 ==> disable)

       double Bondi_InBC_R;         // radius of the void region (=Bondi_InBC_NCell*dh[MAX_LEVEL])
       double Bondi_InBC_P;         // pressure inside the void region
       double Bondi_InBC_E;         // energy inside the void region
       double Bondi_Soften_R;       // soften length (=Bondi_Soften_NCell*dh[MAX_LEVEL])
static double Bondi_P0;             // background pressure
static double Bondi_Cs;             // background sound speed
static double Bondi_RS;             // Schwarzschild radius
static double Bondi_RB;             // Bondi radius
static double Bondi_TimeB;          // Bondi time

       double Bondi_SinkMass;       // total mass             in the void region removed in one global time-step
       double Bondi_SinkMomX;       // total x-momentum       ...
       double Bondi_SinkMomY;       // total y-momentum       ...
       double Bondi_SinkMomZ;       // total z-momentum       ...
       double Bondi_SinkMomXAbs;    // total |x-momentum|     ...
       double Bondi_SinkMomYAbs;    // total |y-momentum|     ...
       double Bondi_SinkMomZAbs;    // total |z-momentum|     ...
       double Bondi_SinkEk;         // total kinematic energy ...
       double Bondi_SinkEt;         // total thermal   energy ...
       int    Bondi_SinkNCell;      // total number of finest cells within the void region

// external units in cgs
const double UnitExt_L = Const_kpc;
const double UnitExt_D = 1.0;
const double UnitExt_M = Const_Msun;
const double UnitExt_E = Const_keV;


// hydrostatic equilibrium (HSE)
static bool   Bondi_HSE;            // enable HSE
static int    Bondi_HSE_Mode;       // initial configuration (1:T=Bondi_T0, 2:rho~1/r, 3:beta model)
static double Bondi_HSE_Dens_NormR; // normalize the density profile to density(r=NormR)=NormD
static double Bondi_HSE_Dens_NormD; // see Bondi_HSE_Dens_NormR

// parameters for Bondi_HSE_Mode=1
static int    Bondi_HSE_Dens_NBin;  // number of bins in the density profile table
static double Bondi_HSE_Dens_MinR;  // minimum radius in the density profile
static double Bondi_HSE_Dens_MaxR;  // maximum ...
static bool   Bondi_HSE_Truncate;   // truncate density within r<TrunR to density=TrunD
static double Bondi_HSE_TrunR;      // see Bondi_HSE_Truncate
static double Bondi_HSE_TrunD;      // see Bondi_HSE_Truncate
static double Bondi_HSE_TrunSmoothR;// smooth out density within TrunR-SmoothR<r<TrunR+SmoothR

static double *Bondi_HSE_DensProf[2] = { NULL, NULL };   // density profile table: [0/1] = [radius/density]

// parameters for Bondi_HSE_Mode=2
static bool   Bondi_HSE_Pres_NormT;    // true --> adjust P2 in the pressure profile such that T(r=NormR)=Bondi_T0
static double Bondi_HSE_Pres_NormP1;   // P=P1*r^-2+P2 (P2=0 when Bondi_HSE_Pres_NormT=false)
static double Bondi_HSE_Pres_NormP2;

// parameters for Bondi_HSE_Mode=3 (beta model)
const  double Bondi_HSE_Beta = 2.0/3.0;   // beta (must be 2/3 for now)
static double Bondi_HSE_Beta_Rcore;       // core radius (input parameter)
static double Bondi_HSE_Beta_Rho0;        // peak density (set by Bondi_HSE_Dens_NormR/D)
static double Bondi_HSE_Beta_P1;          // P(r) = P1*( 1/x + atan(x) ) + P2 assuming beta=2/3, where x=r/Rcore,
static double Bondi_HSE_Beta_P2;          // P1=G*MassBH*Rho0/Rcore, and P2 currently fixed to -0.5*pi*P1 so that P(inf)=0
// =======================================================================================


// problem-specific function prototypes
void Init_ExtAcc_Bondi();
void Record_Bondi();
bool Flag_Bondi( const int i, const int j, const int k, const int lv, const int PID, const double *Threshold );
int Flu_ResetByUser_Func_Bondi( real fluid[], const double Emag, const double x, const double y, const double z, const double Time,
                                const double dt, const int lv, double AuxArray[] );
void Flu_ResetByUser_API_Bondi( const int lv, const int FluSg, const int MagSg, const double TimeNew, const double dt );
static void HSE_SetDensProfileTable();

// this test problem needs to reset both Flu_ResetByUser_API_Ptr and Flu_ResetByUser_Func_Ptr, while
// the former is not defined in TestProb.h (because it's rarely required)
extern void (*Flu_ResetByUser_API_Ptr)( const int lv, const int FluSg, const int MagSg, const double TimeNew, const double dt );




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

#  ifndef GRAVITY
   Aux_Error( ERROR_INFO, "GRAVITY must be enabled !!\n" );
#  endif

#  ifdef COMOVING
   Aux_Error( ERROR_INFO, "COMOVING must be disabled !!\n" );
#  endif

#  ifdef PARTICLE
   Aux_Error( ERROR_INFO, "PARTICLE must be disabled !!\n" );
#  endif

   if ( !OPT__RESET_FLUID || !OPT__RESET_FLUID_INIT )
      Aux_Error( ERROR_INFO, "must enable both OPT__RESET_FLUID and OPT__RESET_FLUID_INIT !!\n" );

#  ifdef GRAVITY
   if ( OPT__BC_POT != BC_POT_ISOLATED )
      Aux_Error( ERROR_INFO, "please set \"OPT__BC_POT = 2\" (i.e., isolated gravity) !!\n" );

   if ( OPT__SELF_GRAVITY  &&  MPI_Rank == 0 )
      Aux_Message( stderr, "WARNING : are you sure about enabling \"OPT__SELF_GRAVITY\" !?\n" );

   if ( OPT__EXT_ACC != EXT_ACC_FUNC  &&  MPI_Rank == 0 )
      Aux_Message( stderr, "WARNING : OPT__EXT_ACC != EXT_ACC_FUNC (1) !?\n" );
#  endif

#  ifndef DUAL_ENERGY
   if ( MPI_Rank == 0 )
   {
      Aux_Message( stderr, "WARNING : it's recommended to enable DUAL_ENERGY, espeically for small GAMMA and large MAX_LEVEL !!\n" );

      if ( amr->BoxSize[0] != amr->BoxSize[1]  ||  amr->BoxSize[0] != amr->BoxSize[2] )
         Aux_Message( stderr, "WARNING : simulation domain is not cubic !?\n" );
   }
#  endif


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
   const char FileName[] = "Input__TestProb";
   ReadPara_t *ReadPara  = new ReadPara_t;

// (1-1) add parameters in the following format:
// add parameters in the following format (some handy constants are defined in TestProb.h):
// --> note that VARIABLE, DEFAULT, MIN, and MAX must have the same data type
// ********************************************************************************************************************************
// ReadPara->Add( "KEY_IN_THE_FILE",      &VARIABLE,                  DEFAULT,      MIN,              MAX               );
// ********************************************************************************************************************************
   ReadPara->Add( "Bondi_MassBH",         &Bondi_MassBH,             -1.0,          Eps_double,       NoMax_double      );
   ReadPara->Add( "Bondi_Rho0",           &Bondi_Rho0,               -1.0,          Eps_double,       NoMax_double      );
   ReadPara->Add( "Bondi_T0",             &Bondi_T0,                 -1.0,          Eps_double,       NoMax_double      );
   ReadPara->Add( "Bondi_RefineRadius0",  &Bondi_RefineRadius0,      -1.0,          Eps_double,       NoMax_double      );
   ReadPara->Add( "Bondi_HalfMaxLvRefR",  &Bondi_HalfMaxLvRefR,       true,         Useless_bool,     Useless_bool      );
   ReadPara->Add( "Bondi_InBC_Rho",       &Bondi_InBC_Rho,           -1.0,          Eps_double,       NoMax_double      );
   ReadPara->Add( "Bondi_InBC_T",         &Bondi_InBC_T,             -1.0,          Eps_double,       NoMax_double      );
   ReadPara->Add( "Bondi_InBC_NCell",     &Bondi_InBC_NCell,         -1.0,          Eps_double,       NoMax_double      );
   ReadPara->Add( "Bondi_Soften_NCell",   &Bondi_Soften_NCell,       -1.0,          NoMin_double,     NoMax_double      );

   ReadPara->Add( "Bondi_HSE",            &Bondi_HSE,                 false,        Useless_bool,     Useless_bool      );
   ReadPara->Add( "Bondi_HSE_Mode",       &Bondi_HSE_Mode,            1,            1,                3                 );
   ReadPara->Add( "Bondi_HSE_Dens_NBin",  &Bondi_HSE_Dens_NBin,       10000,        2,                NoMax_int         );
   ReadPara->Add( "Bondi_HSE_Dens_MinR",  &Bondi_HSE_Dens_MinR,      -1.0,          NoMin_double,     NoMax_double      );
   ReadPara->Add( "Bondi_HSE_Dens_MaxR",  &Bondi_HSE_Dens_MaxR,      -1.0,          NoMin_double,     NoMax_double      );
   ReadPara->Add( "Bondi_HSE_Dens_NormR", &Bondi_HSE_Dens_NormR,     -1.0,          NoMin_double,     NoMax_double      );
   ReadPara->Add( "Bondi_HSE_Dens_NormD", &Bondi_HSE_Dens_NormD,     -1.0,          Eps_double,       NoMax_double      );
   ReadPara->Add( "Bondi_HSE_Truncate",   &Bondi_HSE_Truncate,        true,         Useless_bool,     Useless_bool      );
   ReadPara->Add( "Bondi_HSE_TrunR",      &Bondi_HSE_TrunR,          -1.0,          NoMin_double,     NoMax_double      );
   ReadPara->Add( "Bondi_HSE_TrunD",      &Bondi_HSE_TrunD,          -1.0,          Eps_double,       NoMax_double      );
   ReadPara->Add( "Bondi_HSE_TrunSmoothR",&Bondi_HSE_TrunSmoothR,    -1.0,          NoMin_double,     NoMax_double      );
   ReadPara->Add( "Bondi_HSE_Pres_NormT", &Bondi_HSE_Pres_NormT,      false,        Useless_bool,     Useless_bool      );
   ReadPara->Add( "Bondi_HSE_Beta_Rcore", &Bondi_HSE_Beta_Rcore,     -1.0,          Eps_double,       NoMax_double      );

   ReadPara->Read( FileName );

   delete ReadPara;

// (1-2) set the default values

// (1-3) convert from external to internal units
   Bondi_MassBH        *= UnitExt_M/UNIT_M;
   Bondi_Rho0          *= UnitExt_D/UNIT_D;
   Bondi_T0            *= UnitExt_E/UNIT_E;
   Bondi_RefineRadius0 *= UnitExt_L/UNIT_L;
   Bondi_InBC_Rho      *= UnitExt_D/UNIT_D;
   Bondi_InBC_T        *= UnitExt_E/UNIT_E;

// (1-4) check the runtime parameters
   if ( Bondi_HSE )
   {
      for (int s=0; s<6; s++)
         if ( OPT__BC_FLU[s] != BC_FLU_USER  &&  MPI_Rank == 0 )
            Aux_Message( stderr, "WARNING : OPT__BC_FLU[%d] != BC_FLU_USER for HSE setup !?\n", s );
   }

   else
   {
      for (int s=0; s<6; s++)
         if ( OPT__BC_FLU[s] != BC_FLU_OUTFLOW  &&  MPI_Rank == 0 )
            Aux_Message( stderr, "WARNING : OPT__BC_FLU[%d] != BC_FLU_OUTFLOW for non-HSE setup !?\n", s );
   }


// (2) set the problem-specific derived parameters
// must initialize EoS first
   EoS_Init();

   if ( EoS_DensTemp2Pres_CPUPtr == NULL )
      Aux_Error( ERROR_INFO, "EoS_DensTemp2Pres_CPUPtr == NULL !!\n" );

// the following assumes that EoS requires no passive scalars
   Bondi_InBC_R   = Bondi_InBC_NCell*amr->dh[MAX_LEVEL];
   Bondi_InBC_P   = EoS_DensTemp2Pres_CPUPtr( Bondi_InBC_Rho, Bondi_InBC_T*UNIT_E/Const_kB, NULL,
                                              EoS_AuxArray_Flt, EoS_AuxArray_Int, h_EoS_Table );
   Bondi_InBC_E   = EoS_DensPres2Eint_CPUPtr( Bondi_InBC_Rho, Bondi_InBC_P, NULL, EoS_AuxArray_Flt, EoS_AuxArray_Int,
                                              h_EoS_Table );
   Bondi_Soften_R = Bondi_Soften_NCell*amr->dh[MAX_LEVEL];
   Bondi_P0       = EoS_DensTemp2Pres_CPUPtr( Bondi_Rho0, Bondi_T0*UNIT_E/Const_kB, NULL,
                                              EoS_AuxArray_Flt, EoS_AuxArray_Int, h_EoS_Table );
   Bondi_Cs       = sqrt(  EoS_DensPres2CSqr_CPUPtr( Bondi_Rho0, Bondi_P0, NULL, EoS_AuxArray_Flt, EoS_AuxArray_Int,
                                                     h_EoS_Table )  );
#  ifdef GRAVITY
   Bondi_RS       = 2.0*NEWTON_G*Bondi_MassBH/SQR(Const_c/UNIT_V);
   Bondi_RB       =     NEWTON_G*Bondi_MassBH/SQR(Bondi_Cs);
#  endif
   Bondi_TimeB    = Bondi_RB/Bondi_Cs;


// (3) initialize the HSE setup
   if ( Bondi_HSE )
   {
//    (3-1) set the default values
      if ( Bondi_HSE_Dens_MinR < 0.0 ) {
         Bondi_HSE_Dens_MinR = 0.1*amr->dh[MAX_LEVEL]*UNIT_L/UnitExt_L;
         if ( MPI_Rank == 0 )    Aux_Message( stdout, "NOTE : parameter [%-25s] is set to the default value [%- 21.14e]\n",
                                              "Bondi_HSE_Dens_MinR", Bondi_HSE_Dens_MinR );
      }

      if ( Bondi_HSE_Dens_MaxR < 0.0 ) {
         Bondi_HSE_Dens_MaxR = (0.5*sqrt(3.0)*amr->BoxSize[0]+PS1*amr->dh[0])*UNIT_L/UnitExt_L;
         if ( MPI_Rank == 0 )    Aux_Message( stdout, "NOTE : parameter [%-25s] is set to the default value [%- 21.14e]\n",
                                              "Bondi_HSE_Dens_MaxR", Bondi_HSE_Dens_MaxR );
      }

      if ( Bondi_HSE_Dens_NormR < 0.0 ) {
         Bondi_HSE_Dens_NormR = Bondi_RB*UNIT_L/UnitExt_L;
         if ( MPI_Rank == 0 )    Aux_Message( stdout, "NOTE : parameter [%-25s] is set to the default value [%- 21.14e]\n",
                                              "Bondi_HSE_Dens_NormR", Bondi_HSE_Dens_NormR );
      }

      if ( Bondi_HSE_TrunR < 0.0 ) {
         Bondi_HSE_TrunR = Bondi_RB*UNIT_L/UnitExt_L;
         if ( MPI_Rank == 0 )    Aux_Message( stdout, "NOTE : parameter [%-25s] is set to the default value [%- 21.14e]\n",
                                              "Bondi_HSE_TrunR", Bondi_HSE_TrunR );
      }

//    (3-2) convert from external to internal units
      Bondi_HSE_Dens_MinR   *= UnitExt_L/UNIT_L;
      Bondi_HSE_Dens_MaxR   *= UnitExt_L/UNIT_L;
      Bondi_HSE_Dens_NormR  *= UnitExt_L/UNIT_L;
      Bondi_HSE_Dens_NormD  *= UnitExt_D/UNIT_D;
      Bondi_HSE_TrunR       *= UnitExt_L/UNIT_L;
      Bondi_HSE_TrunD       *= UnitExt_D/UNIT_D;
      Bondi_HSE_TrunSmoothR *= UnitExt_L/UNIT_L;
      Bondi_HSE_Beta_Rcore  *= UnitExt_L/UNIT_L;

//    (3-3) set the density profile table for HSE
      if      ( Bondi_HSE_Mode == 1 )
         HSE_SetDensProfileTable();

      else if ( Bondi_HSE_Mode == 2 )
      {
         Bondi_HSE_Pres_NormP1 = 0.5*NEWTON_G*Bondi_MassBH*Bondi_HSE_Dens_NormD*Bondi_HSE_Dens_NormR;
         Bondi_HSE_Pres_NormP2 = ( Bondi_HSE_Pres_NormT ) ?
                                 Bondi_T0*Bondi_HSE_Dens_NormD/(MOLECULAR_WEIGHT*Const_mH/UNIT_M)
                                 - Bondi_HSE_Pres_NormP1/SQR(Bondi_HSE_Dens_NormR)
                                 : 0.0;
      }

      else if ( Bondi_HSE_Mode == 3 )
      {
//       Rho0 is for arbitrary beta, but P1/P2 is for beta=2/3 only
         Bondi_HSE_Beta_Rho0 = Bondi_HSE_Dens_NormD*
                               pow( 1.0+SQR(Bondi_HSE_Dens_NormR/Bondi_HSE_Beta_Rcore), 1.5*Bondi_HSE_Beta );
         Bondi_HSE_Beta_P1   = NEWTON_G*Bondi_MassBH*Bondi_HSE_Beta_Rho0/Bondi_HSE_Beta_Rcore;
         Bondi_HSE_Beta_P2   = -0.5*M_PI*Bondi_HSE_Beta_P1;    // such that P(inf)=0
      }
   } // if ( Bondi_HSE )


// (4) reset other general-purpose parameters
//     --> a helper macro PRINT_RESET_PARA is defined in Macro.h
   const long   End_Step_Default = __INT_MAX__;
   const double End_T_Default    = 1.0e1*Bondi_TimeB;    // 10 Bondi time

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
      Aux_Message( stdout, "  test problem ID       = %d\n",                     TESTPROB_ID                                                   );
      Aux_Message( stdout, "  Bondi_MassBH          = %13.7e (%13.7e Msun)\n",   Bondi_MassBH, Bondi_MassBH*UNIT_M/Const_Msun                  );
      Aux_Message( stdout, "  Bondi_Rho0            = %13.7e (%13.7e g/cm^3)\n", Bondi_Rho0, Bondi_Rho0*UNIT_D                                 );
      Aux_Message( stdout, "  Bondi_T0              = %13.7e (%13.7e keV)\n",    Bondi_T0, Bondi_T0*UNIT_E/Const_keV                           );
      Aux_Message( stdout, "  Bondi_RefineRadius0   = %13.7e (%13.7e kpc)\n",    Bondi_RefineRadius0, Bondi_RefineRadius0*UNIT_L/Const_kpc     );
      Aux_Message( stdout, "  Bondi_HalfMaxLvRefR   = %s\n",                     (Bondi_HalfMaxLvRefR)?"YES":"NO"                              );
      Aux_Message( stdout, "  Bondi_InBC_Rho        = %13.7e (%13.7e g/cm^3)\n", Bondi_InBC_Rho, Bondi_InBC_Rho*UNIT_D                         );
      Aux_Message( stdout, "  Bondi_InBC_T          = %13.7e (%13.7e keV)\n",    Bondi_InBC_T, Bondi_InBC_T*UNIT_E/Const_keV                   );
      Aux_Message( stdout, "  Bondi_InBC_NCell      = %13.7e\n",                 Bondi_InBC_NCell                                              );
      Aux_Message( stdout, "  Bondi_InBC_R          = %13.7e (%13.7e kpc)\n",    Bondi_InBC_R, Bondi_InBC_R*UNIT_L/Const_kpc                   );
      Aux_Message( stdout, "  Bondi_InBC_E          = %13.7e\n",                 Bondi_InBC_E                                                  );
      Aux_Message( stdout, "  Bondi_Soften_NCell    = %13.7e\n",                 Bondi_Soften_NCell                                            );
      Aux_Message( stdout, "  Bondi_Soften_R        = %13.7e (%13.7e kpc)\n",    Bondi_Soften_R, Bondi_Soften_R*UNIT_L/Const_kpc               );
      Aux_Message( stdout, "  Bondi_Cs              = %13.7e (%13.7e km/s)\n",   Bondi_Cs, Bondi_Cs*UNIT_V/Const_km                            );
      Aux_Message( stdout, "  Schwarzschild radius  = %13.7e (%13.7e kpc)\n",    Bondi_RS, Bondi_RS*UNIT_L/Const_kpc                           );
      Aux_Message( stdout, "  Bondi         radius  = %13.7e (%13.7e kpc)\n",    Bondi_RB, Bondi_RB*UNIT_L/Const_kpc                           );
      Aux_Message( stdout, "  Bondi         time    = %13.7e (%13.7e Myr)\n",    Bondi_TimeB, Bondi_TimeB*UNIT_T/Const_Myr                     );

      Aux_Message( stdout, "  Bondi_HSE             = %s\n",                     (Bondi_HSE)?"YES":"NO"                                        );
      if ( Bondi_HSE ) {
      Aux_Message( stdout, "  Bondi_HSE_Mode        = %d\n",                     Bondi_HSE_Mode                                                );
      Aux_Message( stdout, "  Bondi_HSE_Dens_NBin   = %d\n",                     Bondi_HSE_Dens_NBin                                           );
      Aux_Message( stdout, "  Bondi_HSE_Dens_MinR   = %13.7e (%13.7e kpc)\n",    Bondi_HSE_Dens_MinR, Bondi_HSE_Dens_MinR*UNIT_L/Const_kpc     );
      Aux_Message( stdout, "  Bondi_HSE_Dens_MaxR   = %13.7e (%13.7e kpc)\n",    Bondi_HSE_Dens_MaxR, Bondi_HSE_Dens_MaxR*UNIT_L/Const_kpc     );
      Aux_Message( stdout, "  Bondi_HSE_Dens_NormR  = %13.7e (%13.7e kpc)\n",    Bondi_HSE_Dens_NormR, Bondi_HSE_Dens_NormR*UNIT_L/Const_kpc   );
      Aux_Message( stdout, "  Bondi_HSE_Dens_NormD  = %13.7e (%13.7e g/cm^3)\n", Bondi_HSE_Dens_NormD, Bondi_HSE_Dens_NormD*UNIT_D             );
      Aux_Message( stdout, "  Bondi_HSE_Truncate    = %s\n",                     (Bondi_HSE_Truncate)?"YES":"NO"                               );
      Aux_Message( stdout, "  Bondi_HSE_TrunR       = %13.7e (%13.7e kpc)\n",    Bondi_HSE_TrunR, Bondi_HSE_TrunR*UNIT_L/Const_kpc             );
      Aux_Message( stdout, "  Bondi_HSE_TrunD       = %13.7e (%13.7e g/cm^3)\n", Bondi_HSE_TrunD, Bondi_HSE_TrunD*UNIT_D                       );
      Aux_Message( stdout, "  Bondi_HSE_TrunSmoothR = %13.7e (%13.7e kpc)\n",    Bondi_HSE_TrunSmoothR, Bondi_HSE_TrunSmoothR*UNIT_L/Const_kpc );
      Aux_Message( stdout, "  Bondi_HSE_Pres_NormT  = %s\n",                     (Bondi_HSE_Pres_NormT)?"YES":"NO"                             );
      Aux_Message( stdout, "  Bondi_HSE_Beta        = %13.7e\n",                 Bondi_HSE_Beta                                                );
      Aux_Message( stdout, "  Bondi_HSE_Beta_Rho0   = %13.7e (%13.7e g/cm^3)\n", Bondi_HSE_Beta_Rho0, Bondi_HSE_Beta_Rho0*UNIT_D               );
      Aux_Message( stdout, "  Bondi_HSE_Beta_Rcore  = %13.7e (%13.7e kpc)\n",    Bondi_HSE_Beta_Rcore, Bondi_HSE_Beta_Rcore*UNIT_L/Const_kpc   ); }
      Aux_Message( stdout, "=============================================================================\n" );
   } // if ( MPI_Rank == 0 )


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

   double Dens, MomX, MomY, MomZ, Pres, Eint, Etot;

   MomX = 0.0;
   MomY = 0.0;
   MomZ = 0.0;

// note that the central void region will be reset by calling Flu_ResetByUser_Func_Bondi()
// hydrostatic equilibrium
   if ( Bondi_HSE )
   {
      const double r = sqrt( SQR(x-amr->BoxCenter[0]) + SQR(y-amr->BoxCenter[1]) + SQR(z-amr->BoxCenter[2]) );

      if ( Bondi_HSE_Mode == 1 )
      {
         const double *Table_R = Bondi_HSE_DensProf[0];
         const double *Table_D = Bondi_HSE_DensProf[1];

         Dens = Mis_InterpolateFromTable( Bondi_HSE_Dens_NBin, Table_R, Table_D, r );

         if ( Dens == NULL_REAL )
            Aux_Error( ERROR_INFO, "interpolation failed at radius %13.7e --> reset Bondi_HSE_Dens_MinR/MaxR !!\n", r );

         Pres = EoS_DensTemp2Pres_CPUPtr( Dens, Bondi_T0*UNIT_E/Const_kB, NULL,
                                          EoS_AuxArray_Flt, EoS_AuxArray_Int, h_EoS_Table );
      }

      else if ( Bondi_HSE_Mode == 2 )
      {
         Dens = Bondi_HSE_Dens_NormD*Bondi_HSE_Dens_NormR/r;
         Pres = Bondi_HSE_Pres_NormP1/SQR(r) + Bondi_HSE_Pres_NormP2;
      }

      else if ( Bondi_HSE_Mode == 3 )
      {
         const double x = r / Bondi_HSE_Beta_Rcore;

         Dens = Bondi_HSE_Beta_Rho0*pow( 1.0+SQR(x), -1.5*Bondi_HSE_Beta );   // for arbitrary beta
         Pres = Bondi_HSE_Beta_P1*( 1.0/x + atan(x) ) + Bondi_HSE_Beta_P2;    // for beta=2/3 only
      }

      else
         Aux_Error( ERROR_INFO, "unsupported Bondi_HSE_Mode (%d) !!\n", Bondi_HSE_Mode );
   } // if ( Bondi_HSE )


// uniform background
   else
   {
      Dens = Bondi_Rho0;
      Pres = Bondi_P0;
   }


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



//-------------------------------------------------------------------------------------------------------
// Function    :  HSE_SetDensProfileTable
// Description :  Set up the density profile table for HSE
//
// Note        :  1. Assume isotheral (for Bondi_HSE_Mode=1)
//                   --> log(rho2/rho1) = G*MassBH/A*( 1/r2 - 1/r1 ), where A = kB*T/(mu*mH)
//                   --> r1 and rho1 are for normalization and are set by Bondi_HSE_Dens_NormR/D
//                2. Invoked by SetParameter()
//                3. Truncate density when enabling Bondi_HSE_Truncate
//                4. Smooth out density around the trucation radius
//
// Parameter   :  None
//-------------------------------------------------------------------------------------------------------
void HSE_SetDensProfileTable()
{

// allocate table --> deallocated by End_Bondi()
   for (int v=0; v<2; v++)    Bondi_HSE_DensProf[v] = new double [Bondi_HSE_Dens_NBin];


// set radius
// --> r_min should be smaller than 0.5*dh_min for proper interpolation
   const int    NBin  = Bondi_HSE_Dens_NBin;
   const double r_min = Bondi_HSE_Dens_MinR;
   const double r_max = Bondi_HSE_Dens_MaxR;
   const double dr    = pow( r_max/r_min, 1.0/(NBin-1.0) );

   for (int b=0; b<NBin; b++)    Bondi_HSE_DensProf[0][b] = r_min*pow( dr, (double)b );


// set density
   const double A = Bondi_T0/(MOLECULAR_WEIGHT*Const_mH/UNIT_M);
   const double B = NEWTON_G*Bondi_MassBH/A;

   for (int b=0; b<NBin; b++)
      Bondi_HSE_DensProf[1][b] = Bondi_HSE_Dens_NormD*exp( B*(1.0/Bondi_HSE_DensProf[0][b] - 1.0/Bondi_HSE_Dens_NormR) );


// truncate and smooth out density
   if ( Bondi_HSE_Truncate )
   {
//    truncation
      for (int b=0; b<NBin; b++)
      {
         if ( Bondi_HSE_DensProf[0][b] < Bondi_HSE_TrunR )
            Bondi_HSE_DensProf[1][b] = Bondi_HSE_TrunD;
         else
            break;
      }

//    smoothing
      if ( Bondi_HSE_TrunSmoothR > 0.0 )
      {
         double r1, r2;
         long   n;

//       backup the original profile
         double *DensProf_Ori = new double [Bondi_HSE_Dens_NBin];
         memcpy( DensProf_Ori, Bondi_HSE_DensProf[1], Bondi_HSE_Dens_NBin*sizeof(double) );


         for (int b1=0; b1<NBin; b1++)
         {
            r1 = Bondi_HSE_DensProf[0][b1];

//          smooth out density within TrunR-SmoothR<r<TrunR+SmoothR
            if ( r1 > Bondi_HSE_TrunR-Bondi_HSE_TrunSmoothR  &&  r1 < Bondi_HSE_TrunR+Bondi_HSE_TrunSmoothR )
            {
               Bondi_HSE_DensProf[1][b1] = 0.0;
               n = 0;

               for (int b2=0; b2<NBin; b2++)
               {
                  r2 = Bondi_HSE_DensProf[0][b2];

//                we also adopt Bondi_HSE_TrunSmoothR as the smoothing kernel width
                  if ( fabs(r1-r2) < Bondi_HSE_TrunSmoothR )
                  {
                     Bondi_HSE_DensProf[1][b1] += DensProf_Ori[b2];
                     n++;
                  }
               }

               Bondi_HSE_DensProf[1][b1] /= (double)n;

//             do not allow the smoothed out data to exceed the original values
               Bondi_HSE_DensProf[1][b1] = MIN( Bondi_HSE_DensProf[1][b1], DensProf_Ori[b1] );
            }
         } // for (int b1=0; b1<NBin; b1++)

         delete [] DensProf_Ori;
      } // if ( Bondi_HSE_TrunSmoothR > 0.0 )
   } // Bondi_HSE_Truncate

} // FUNCTION : HSE_SetDensProfileTable


//-------------------------------------------------------------------------------------------------------
// Function    :  End_Bondi
// Description :  Free memory before terminating the program
//
// Note        :  1. Linked to the function pointer "End_User_Ptr" to replace "End_User()"
//
// Parameter   :  None
//-------------------------------------------------------------------------------------------------------
void End_Bondi()
{

   for (int v=0; v<2; v++)
   {
      delete [] Bondi_HSE_DensProf[v];
      Bondi_HSE_DensProf[v] = NULL;
   }

} // FUNCTION : End_Bondi
#endif // #if ( MODEL == HYDRO  &&  defined GRAVITY )



//-------------------------------------------------------------------------------------------------------
// Function    :  Init_TestProb_Hydro_Bondi
// Description :  Test problem initializer
//
// Note        :  None
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Init_TestProb_Hydro_Bondi()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );


// validate the compilation flags and runtime parameters
   Validate();


#  if ( MODEL == HYDRO  &&  defined GRAVITY )
// set the problem-specific runtime parameters
   SetParameter();


// set the function pointers of various problem-specific routines
   Init_Function_User_Ptr   = SetGridIC;
   Flag_User_Ptr            = Flag_Bondi;
   Aux_Record_User_Ptr      = Record_Bondi;
   BC_User_Ptr              = SetGridIC;
   Flu_ResetByUser_Func_Ptr = Flu_ResetByUser_Func_Bondi;
   Flu_ResetByUser_API_Ptr  = Flu_ResetByUser_API_Bondi;
   End_User_Ptr             = End_Bondi;
#  ifdef GRAVITY
   Init_ExtAcc_Ptr         = Init_ExtAcc_Bondi;
#  endif
#  endif // #if ( MODEL == HYDRO  &&  defined GRAVITY )


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : Init_TestProb_Hydro_Bondi
