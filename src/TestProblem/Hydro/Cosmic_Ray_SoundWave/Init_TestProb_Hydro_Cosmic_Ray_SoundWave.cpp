#include "GAMER.h"
#include "TestProb.h"



static void OutputError();


// problem-specific global variables
// =======================================================================================
static double CR_Acoustic_Delta;        // amplitude of the velocity perturbation
static double CR_Acoustic_Rho0;         // background density
static double CR_Acoustic_Pres0;        // background pressure
static double CR_Acoustic_Pres_CR0;     // background pressure of cosmic rays
static double CR_Acoustic_V0;           // background velocity
static double CR_Acoustic_Sign;         // (+1/-1) --> (right/left-moving wave)
static double CR_Acoustic_Phase;        // initial phase shift
static int    CR_Acoustic_Dir;          // wave direction (0/1) --> (x/diagonal)
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

#  if ( defined COSMIC_ARY && EOS != EOS_COSMIC_RAY )
   Aux_Error( ERROR_INFO, "EOS != EOS_COSMIC_RAY when enable COSMIC_RAY!!\n" );
#  endif

#  endif // #ifndef COSMIC_RAY

// warnings


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Validating test problem %d ... done\n", TESTPROB_ID );

} // FUNCTION : Validate



// replace HYDRO by the target model (e.g., MHD/ELBDM) and also check other compilation flags if necessary (e.g., GRAVITY/PARTICLE)
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
// ReadPara->Add( "KEY_IN_THE_FILE",      &VARIABLE,              DEFAULT,       MIN,              MAX               );
// ********************************************************************************************************************************
//   ReadPara->Add( "GAMMA_CR",             &GAMMA_CR,                  0.0,           0.0,              NoMax_double);
   ReadPara->Add( "CR_Acoustic_Delta",    &CR_Acoustic_Delta,         0.0,           0.0,              NoMax_double);
   ReadPara->Add( "CR_Acoustic_Rho0",     &CR_Acoustic_Rho0,          0.0,           0.0,              NoMax_double);
   ReadPara->Add( "CR_Acoustic_Pres0",    &CR_Acoustic_Pres0,         0.0,           0.0,              NoMax_double);
   ReadPara->Add( "CR_Acoustic_Pres_CR0", &CR_Acoustic_Pres_CR0,      0.0,           0.0,              NoMax_double);
   ReadPara->Add( "CR_Acoustic_V0",       &CR_Acoustic_V0,            0.0,           0.0,              NoMax_double);
   ReadPara->Add( "CR_Acoustic_Sign",     &CR_Acoustic_Sign,          0.0,           0.0,              NoMax_double);
   ReadPara->Add( "CR_Acoustic_Phase",    &CR_Acoustic_Phase,         0.0,           0.0,              NoMax_double);
   ReadPara->Add( "CR_Acoustic_Dir",      &CR_Acoustic_Dir,             0,             0,                 NoMax_int);

   ReadPara->Read( FileName );

   delete ReadPara;

// (1-2) set the default values
// force Acoustic_Sign to be +1.0/-1.0 
   if ( CR_Acoustic_Sign >= 0.0 )   CR_Acoustic_Sign = +1.0;
   else                             CR_Acoustic_Sign = -1.0;

// (1-3) check the runtime parameters


// (2) set the problem-specific derived parameters

// (3) reset other general-purpose parameters
//     --> a helper macro PRINT_WARNING is defined in TestProb.h
   const long   End_Step_Default = __INT_MAX__;
   const double End_T_Default    = __FLT_MAX__;

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
      Aux_Message( stdout, "  test problem ID           = %d\n",         TESTPROB_ID );
      Aux_Message( stdout, "  CR_Acoustic_Delta         = %14.7e\n",     CR_Acoustic_Delta );
      Aux_Message( stdout, "  CR_Acoustic_Rho0          = %14.7e\n",     CR_Acoustic_Rho0 );
      Aux_Message( stdout, "  CR_Acoustic_Pres0         = %14.7e\n",     CR_Acoustic_Pres0 );
      Aux_Message( stdout, "  CR_Acoustic_Pres_CR0      = %14.7e\n",     CR_Acoustic_Pres_CR0 );
      Aux_Message( stdout, "  CR_Acoustic_V0            = %14.7e\n",     CR_Acoustic_V0 );
      Aux_Message( stdout, "  CR_Acoustic_Sign          = %14.7e\n",     CR_Acoustic_Sign );
      Aux_Message( stdout, "  CR_Acoustic_Phase         = %14.7e\n",     CR_Acoustic_Phase );
      Aux_Message( stdout, "  CR_Acoustic_Dir           = %d\n",         CR_Acoustic_Dir );
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
#  ifdef COSMIC_RAY
   double cs = SQRT( ( GAMMA * CR_Acoustic_Pres0 + GAMMA_CR * CR_Acoustic_Pres_CR0 ) / CR_Acoustic_Rho0 );
#  else
   double cs = SQRT( ( GAMMA * CR_Acoustic_Pres0 ) / CR_Acoustic_Rho0 );
#  endif
   double delta_cs = CR_Acoustic_Delta / cs;
   
   double wavelength, WaveK, WaveW, r, wave;
   if ( CR_Acoustic_Dir == 0 )
   {
      wavelength = 0.5;
      WaveK = 2.0*M_PI/wavelength;
      WaveW = WaveK * cs;
      r = x - CR_Acoustic_V0 * Time;
      wave = SIN( WaveK * r - CR_Acoustic_Sign * WaveW * Time + CR_Acoustic_Phase );

      Dens = ( 1.0 + delta_cs * wave ) * CR_Acoustic_Rho0;
      MomX = Dens * ( CR_Acoustic_Sign * CR_Acoustic_Delta * wave + CR_Acoustic_V0);
      MomY = 1.0 + Dens * 1.e-16 * sin( 8. * M_PI * y );
      MomZ = 0.0;
      Pres = ( 1.0 + delta_cs * wave * GAMMA ) * CR_Acoustic_Pres0;
   }
   else if ( CR_Acoustic_Dir == 1 )
   {
      wavelength = SQRT(3.0) / 3.;
      WaveK = 2.0*M_PI/wavelength;
      WaveW = WaveK * cs;
      r = (x + y + z) / SQRT(3.) - CR_Acoustic_V0 * Time;
      wave = SIN( WaveK * r - CR_Acoustic_Sign * WaveW * Time + CR_Acoustic_Phase );

      Dens = ( 1.0 + delta_cs * wave ) * CR_Acoustic_Rho0;
      MomX = Dens * ( CR_Acoustic_Sign * CR_Acoustic_Delta * wave + CR_Acoustic_V0) / SQRT(3.0);
      MomY = MomX;
      MomZ = MomX;
      Pres = ( 1.0 + delta_cs * wave * GAMMA ) * CR_Acoustic_Pres0;
   }
   else
   {
      Aux_Error( ERROR_INFO, "CR_Acoustic_Dir = %d is NOT supported [0/1] !!\n", CR_Acoustic_Dir );
   }
   
#  ifdef COSMIC_RAY
   double GAMMA_CR_m1_inv = 1.0 / (GAMMA_CR - 1.0);
   P_cr = ( 1.0 + delta_cs * wave * GAMMA_CR ) * CR_Acoustic_Pres_CR0;
   Pres = Pres + P_cr;
   CRay = GAMMA_CR_m1_inv * P_cr;

// set the output array of passive scaler
   fluid[CRAY] = CRay;
#  endif
   
   Eint = EoS_DensPres2Eint_CPUPtr( Dens, Pres, fluid+NCOMP_FLUID, EoS_AuxArray_Flt, EoS_AuxArray_Int, h_EoS_Table );
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
void OutputError()
{

   const char Prefix[100]     = "CosmicRay_Acousticwave";
   if ( CR_Acoustic_Dir == 0 ){
      const OptOutputPart_t Part = OUTPUT_X;
#     ifdef MHD
      Output_L1Error( SetGridIC, SetBFieldIC, Prefix, Part, OUTPUT_PART_X, OUTPUT_PART_Y, OUTPUT_PART_Z );
#     else
      Output_L1Error( SetGridIC, NULL,        Prefix, Part, OUTPUT_PART_X, OUTPUT_PART_Y, OUTPUT_PART_Z );
#     endif
   } else if ( CR_Acoustic_Dir == 1 ){
      const OptOutputPart_t Part = OUTPUT_DIAG;
#     ifdef MHD
      Output_L1Error( SetGridIC, SetBFieldIC, Prefix, Part, NULL_REAL, NULL_REAL, NULL_REAL );
#     else
      Output_L1Error( SetGridIC, NULL,        Prefix, Part, NULL_REAL, NULL_REAL, NULL_REAL );
#     endif
   }

} // FUNCTION : OutputError
#endif // #if ( MODEL == HYDRO )



//-------------------------------------------------------------------------------------------------------
// Function    :  Init_TestProb_Hydro_Cosmic_Ray_SoundWave
// Description :  Test problem initializer
//
// Note        :  None
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Init_TestProb_Hydro_Cosmic_Ray_SoundWave()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );


// validate the compilation flags and runtime parameters
   Validate();


#  if ( MODEL == HYDRO )
// set the problem-specific runtime parameters
   SetParameter();


// set the function pointers of various problem-specific routines
   Init_Function_User_Ptr         = SetGridIC;

#  ifdef MHD
   Init_Function_BField_User_Ptr  = SetBFieldIC;
#  endif

   Output_User_Ptr                = OutputError;

#  endif // #if ( MODEL == HYDRO )

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : Init_TestProb_Hydro_Cosmic_Ray_SoundWave
