#include "GAMER.h"
#include "TestProb.h"



// problem-specific global variables
// =======================================================================================
static double Jeans_Rho0;           // background density
static double Jeans_Rho1;           // amplitude of the density perturbation
static double Jeans_P0;             // background pressure
static double Jeans_v0;             // background velocity
static double Jeans_Sign;           // (+1/-1) --> (stable: right/left-moving wave; unstable: growing/decaying mode)
static double Jeans_Phase0;         // initial phase shift
#ifdef MHD
static double Jeans_B0;             // background magnetic field
#endif

static double Jeans_WaveLength;     // wavelength
static double Jeans_WaveK;          // wavenumber
static double Jeans_WaveKj;         // critical wavenumber
static double Jeans_WaveW;          // wave angular frequency
static double Jeans_WaveSpeed;      // propagation speed (sound speed in hydro or fast wave speed in MHD)
static bool   Jeans_Stable;         // true/false --> Jeans stable/unstable
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

#  ifdef COMOVING
   Aux_Error( ERROR_INFO, "COMOVING must be disabled !!\n" );
#  endif

#  ifdef PARTICLE
   Aux_Error( ERROR_INFO, "PARTICLE must be disabled !!\n" );
#  endif

   if ( amr->BoxSize[0] != amr->BoxSize[1]  ||  amr->BoxSize[0] != amr->BoxSize[2] )
      Aux_Error( ERROR_INFO, "simulation domain must be cubic !!\n" );

   for (int f=0; f<6; f++)
   if ( OPT__BC_FLU[f] != BC_FLU_PERIODIC )
      Aux_Error( ERROR_INFO, "must adopt periodic BC for the fluid (i.e., \"OPT__BC_FLU_* = 1\") !!\n" );

#  ifdef GRAVITY
   if ( OPT__BC_POT != BC_POT_PERIODIC )
      Aux_Error( ERROR_INFO, "must adopt periodic BC for the gravity (i.e., \"OPT__BC_POT = 1\") !!\n" );
#  endif


// warnings
   if ( MPI_Rank == 0 )
   {
      if ( !OPT__OUTPUT_USER )   Aux_Message( stdout, "WARNING : OPT__OUTPUT_USER is off !!\n" );
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
   ReadPara->Add( "Jeans_Rho0",        &Jeans_Rho0,           -1.0,           Eps_double,       NoMax_double      );
   ReadPara->Add( "Jeans_Rho1",        &Jeans_Rho1,           -1.0,           0.0,              NoMax_double      );
   ReadPara->Add( "Jeans_P0",          &Jeans_P0,             -1.0,           Eps_double,       NoMax_double      );
   ReadPara->Add( "Jeans_v0",          &Jeans_v0,              0.0,           NoMin_double,     NoMax_double      );
   ReadPara->Add( "Jeans_Sign",        &Jeans_Sign,            1.0,           NoMin_double,     NoMax_double      );
   ReadPara->Add( "Jeans_Phase0",      &Jeans_Phase0,          0.0,           NoMin_double,     NoMax_double      );
#  ifdef MHD
   ReadPara->Add( "Jeans_B0",          &Jeans_B0,             -1.0,           0.0,              NoMax_double      );
#  endif

   ReadPara->Read( FileName );

   delete ReadPara;


// (2) set the problem-specific derived parameters
   Jeans_WaveLength = amr->BoxSize[0] / sqrt(3.0);   // 3 wavelengths along the diagonal
#  ifdef MHD
   Jeans_WaveSpeed  = sqrt( GAMMA*Jeans_P0/Jeans_Rho0 + SQR(Jeans_B0)/Jeans_Rho0 );
#  else
   Jeans_WaveSpeed  = sqrt( GAMMA*Jeans_P0/Jeans_Rho0 );
#  endif
   Jeans_WaveK      = 2.0*M_PI/Jeans_WaveLength;
   Jeans_WaveKj     = sqrt( 4.0*M_PI*NEWTON_G*Jeans_Rho0/SQR(Jeans_WaveSpeed) );
   Jeans_Stable     = ( Jeans_WaveK > Jeans_WaveKj );
   Jeans_WaveW      = Jeans_WaveSpeed*sqrt(  fabs( SQR(Jeans_WaveK) - SQR(Jeans_WaveKj) )  );

// force Jeans_Sign to be +1.0/-1.0
   if ( Jeans_Sign >= 0.0 )  Jeans_Sign = +1.0;
   else                      Jeans_Sign = -1.0;


// (3) reset other general-purpose parameters
//     --> a helper macro PRINT_WARNING is defined in TestProb.h
   const double End_T_Default    = ( Jeans_Stable ) ? 2.0*M_PI/Jeans_WaveW : log(50.0)/Jeans_WaveW;
   const long   End_Step_Default = __INT_MAX__;

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
      Aux_Message( stdout, "  test problem ID        = %d\n",     TESTPROB_ID      );
      Aux_Message( stdout, "  background density     = %14.7e\n", Jeans_Rho0       );
      Aux_Message( stdout, "  density perturbation   = %14.7e\n", Jeans_Rho1       );
      Aux_Message( stdout, "  background pressure    = %14.7e\n", Jeans_P0         );
      Aux_Message( stdout, "  background velocity    = %14.7e\n", Jeans_v0         );
#     ifdef MHD
      Aux_Message( stdout, "  background B field     = %14.7e\n", Jeans_B0         );
#     endif
      Aux_Message( stdout, "  sign (grow/decay;R/L)  = %14.7e\n", Jeans_Sign       );
      Aux_Message( stdout, "  initial phase shift    = %14.7e\n", Jeans_Phase0     );
      Aux_Message( stdout, "  wave speed             = %14.7e\n", Jeans_WaveSpeed  );
      Aux_Message( stdout, "  wavelength             = %14.7e\n", Jeans_WaveLength );
      Aux_Message( stdout, "  wavenumber             = %14.7e\n", Jeans_WaveK      );
      Aux_Message( stdout, "  critial wavenumber     = %14.7e\n", Jeans_WaveKj     );
      Aux_Message( stdout, "  stable                 = %s\n",     (Jeans_Stable)?"YES":"NO" );
      Aux_Message( stdout, "  wave angular frequency = %14.7e\n", Jeans_WaveW      );
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

// short names
   const double Rho0   = Jeans_Rho0;
   const double Rho1   = Jeans_Rho1;
   const double P0     = Jeans_P0;
   const double v0     = Jeans_v0;
   const double Sign   = Jeans_Sign;
   const double Phase0 = Jeans_Phase0;
   const double WaveK  = Jeans_WaveK;
   const double WaveW  = Jeans_WaveW;

   double r, v1, P1, CosPhase, SinPhase, ExpPhase;

   r  = 1.0/sqrt(3.0)*( x + y + z ) - v0*Time;
   v1 = Sign*Rho1/Rho0*WaveW/WaveK;
   P1 = GAMMA*P0/Rho0*Rho1;

   if ( Jeans_Stable )
   {
      CosPhase = cos( WaveK*r - Sign*WaveW*Time + Phase0 );

      fluid[DENS] = Rho0 + Rho1*CosPhase;
      fluid[MOMX] = fluid[DENS]*( v0 + v1*CosPhase )/sqrt(3.0);
      fluid[MOMY] = fluid[MOMX];
      fluid[MOMZ] = fluid[MOMX];
      fluid[ENGY] = 0.5*( SQR(fluid[MOMX]) + SQR(fluid[MOMY]) + SQR(fluid[MOMZ]) )/fluid[DENS]
                    + ( P0 + P1*CosPhase )/(GAMMA-1.0);
   }

   else
   {
      CosPhase = cos( WaveK*r + Phase0 );
      SinPhase = sin( WaveK*r + Phase0 );
      ExpPhase = exp( Jeans_Sign*Jeans_WaveW*Time );

      fluid[DENS] = Rho0 + Rho1*CosPhase*ExpPhase;
      fluid[MOMX] = fluid[DENS]*( v0 - v1*SinPhase*ExpPhase )/sqrt(3.0);
      fluid[MOMY] = fluid[MOMX];
      fluid[MOMZ] = fluid[MOMX];
      fluid[ENGY] = 0.5*( SQR(fluid[MOMX]) + SQR(fluid[MOMY]) + SQR(fluid[MOMZ]) )/fluid[DENS]
                    + ( P0 + P1*CosPhase*ExpPhase )/(GAMMA-1.0);
   }

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

// short names
   const double Rho0   = Jeans_Rho0;
   const double Rho1   = Jeans_Rho1;
   const double v0     = Jeans_v0;
   const double B0     = Jeans_B0;
   const double Sign   = Jeans_Sign;
   const double Phase0 = Jeans_Phase0;
   const double WaveK  = Jeans_WaveK;
   const double WaveW  = Jeans_WaveW;

   double r, B1, CosPhase, ExpPhase;

   r  = 1.0/sqrt(3.0)*( x + y + z ) - v0*Time;
   B1 = B0/Rho0*Rho1;

// assuming transverse component only
   if ( Jeans_Stable )
   {
      CosPhase = cos( WaveK*r - Sign*WaveW*Time + Phase0 );

      magnetic[MAGZ] = ( B0 + B1*CosPhase ) / sqrt(1.5);
   }

   else
   {
      CosPhase = cos( WaveK*r + Phase0 );
      ExpPhase = exp( Jeans_Sign*Jeans_WaveW*Time );

      magnetic[MAGZ] = ( B0 + B1*CosPhase*ExpPhase ) / sqrt(1.5);
   }

   magnetic[MAGX] = -0.5*magnetic[MAGZ];
   magnetic[MAGY] = -0.5*magnetic[MAGZ];

} // FUNCTION : SetBFieldIC
#endif // #ifdef MHD



//-------------------------------------------------------------------------------------------------------
// Function    :  OutputError
// Description :  Output the L1 error
//
// Note        :  1. Invoke Output_L1Error()
//                2. Use SetGridIC() to provide the analytical solution at any given time
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
static void OutputError()
{

   const char Prefix[100]     = "Jeans";
   const OptOutputPart_t Part = OUTPUT_DIAG;

#  ifdef MHD
   Output_L1Error( SetGridIC, SetBFieldIC, Prefix, Part, NULL_REAL, NULL_REAL, NULL_REAL );
#  else
   Output_L1Error( SetGridIC, NULL,        Prefix, Part, NULL_REAL, NULL_REAL, NULL_REAL );
#  endif

} // FUNCTION : OutputError
#endif // #if ( MODEL == HYDRO  &&  defined GRAVITY )



//-------------------------------------------------------------------------------------------------------
// Function    :  Init_TestProb_Hydro_JeansInstability
// Description :  Test problem initializer
//
// Note        :  None
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Init_TestProb_Hydro_JeansInstability()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );


// validate the compilation flags and runtime parameters
   Validate();


#  if ( MODEL == HYDRO  &&  defined GRAVITY )
// set the problem-specific runtime parameters
   SetParameter();


// set the function pointers of various problem-specific routines
   Init_Function_User_Ptr        = SetGridIC;
   Output_User_Ptr               = OutputError;
   Init_Field_User_Ptr           = NULL;
   Init_User_Ptr                 = NULL;
   Flag_User_Ptr                 = NULL;
   Mis_GetTimeStep_User_Ptr      = NULL;
   Aux_Record_User_Ptr           = NULL;
   BC_User_Ptr                   = NULL;
   Flu_ResetByUser_Func_Ptr      = NULL;
   End_User_Ptr                  = NULL;
#  ifdef MHD
   Init_Function_BField_User_Ptr = SetBFieldIC;
   BC_BField_User_Ptr            = NULL;
#  endif
#  endif // if ( MODEL == HYDRO  &&  defined GRAVITY )


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : Init_TestProb_Hydro_JeansInstability
