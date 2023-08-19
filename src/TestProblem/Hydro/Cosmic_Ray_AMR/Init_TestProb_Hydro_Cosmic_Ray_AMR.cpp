#include "GAMER.h"
#include "TestProb.h"





// problem-specific global variables
// =======================================================================================
static double GAMMA_CR;                 // adiabatic index of cosmic rays
static double CR_Diffusion_CenterX;
static double CR_Diffusion_CenterY;
static double CR_Diffusion_CenterZ;
static double CR_Diffusion_Vx;
static double CR_Diffusion_Vy;
static double CR_Diffusion_Vz;
static double CR_Diffusion_Rho0;
static double CR_Diffusion_PGas0;
static double CR_Diffusion_E0_CR;
static double CR_Diffusion_BG_CR;
static double CR_Diffusion_R02_CR;
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

#  if ( EOS != EOS_COSMIC_RAY )
   Aux_Error( ERROR_INFO, "EOS != EOS_COSMIC_RAY !!\n" );
#  endif

// warnings
   if ( MPI_Rank == 0 )
   {
#     ifndef DUAL_ENERGY
         Aux_Message( stderr, "WARNING : it's recommended to enable DUAL_ENERGY for this test !!\n" );
#     endif
      // What is this for?@@@
      if ( FLAG_BUFFER_SIZE < 5 )
         Aux_Message( stderr, "WARNING : it's recommended to set FLAG_BUFFER_SIZE >= 5 for this test !!\n" );
   } // if ( MPI_Rank == 0 )



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
   ReadPara->Add( "GAMMA_CR",             &GAMMA_CR,                  0.0,           0.0,              NoMax_double);
   ReadPara->Add( "CR_Diffusion_CenterX",      &CR_Diffusion_CenterX,           0.5,           0.0,              NoMax_double);
   ReadPara->Add( "CR_Diffusion_CenterY",      &CR_Diffusion_CenterY,           0.5,           0.0,              NoMax_double);
   ReadPara->Add( "CR_Diffusion_CenterZ",      &CR_Diffusion_CenterZ,           0.5,           0.0,              NoMax_double);
   ReadPara->Add( "CR_Diffusion_Vx",           &CR_Diffusion_Vx,                0.0,           NoMin_double,     NoMax_double);
   ReadPara->Add( "CR_Diffusion_Vy",           &CR_Diffusion_Vy,                0.0,           NoMin_double,     NoMax_double);
   ReadPara->Add( "CR_Diffusion_Vz",           &CR_Diffusion_Vz,                0.0,           NoMin_double,     NoMax_double);
   ReadPara->Add( "CR_Diffusion_Rho0",         &CR_Diffusion_Rho0,              1.0,           0.0,              NoMax_double);
   ReadPara->Add( "CR_Diffusion_PGas0",        &CR_Diffusion_PGas0,             1.0,           0.0,              NoMax_double);
   ReadPara->Add( "CR_Diffusion_E0_CR",        &CR_Diffusion_E0_CR,             0.1,           0.0,              NoMax_double);
   ReadPara->Add( "CR_Diffusion_BG_CR",        &CR_Diffusion_BG_CR,             0.1,           0.0,              NoMax_double);
   ReadPara->Add( "CR_Diffusion_R02_CR",       &CR_Diffusion_R02_CR,            0.05,          0.0,              NoMax_double);

   ReadPara->Read( FileName );

   delete ReadPara;

// (1-2) set the default values
//
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
      Aux_Message( stdout, "  GAMMA_CR                  = %14.7e\n",     GAMMA_CR );
      Aux_Message( stdout, "  CR_Diffusion_CenterX      = %14.7e\n",     CR_Diffusion_CenterX );
      Aux_Message( stdout, "  CR_Diffusion_CenterY      = %14.7e\n",     CR_Diffusion_CenterY );
      Aux_Message( stdout, "  CR_Diffusion_CenterZ      = %14.7e\n",     CR_Diffusion_CenterZ );
      Aux_Message( stdout, "  CR_Diffusion_Vx           = %14.7e\n",     CR_Diffusion_Vx );
      Aux_Message( stdout, "  CR_Diffusion_Vy           = %14.7e\n",     CR_Diffusion_Vy );
      Aux_Message( stdout, "  CR_Diffusion_Vz           = %14.7e\n",     CR_Diffusion_Vz );
      Aux_Message( stdout, "  CR_Diffusion_Rho0         = %14.7e\n",     CR_Diffusion_Rho0 );
      Aux_Message( stdout, "  CR_Diffusion_PGas0        = %14.7e\n",     CR_Diffusion_PGas0 );
      Aux_Message( stdout, "  CR_Diffusion_E0_CR        = %14.7e\n",     CR_Diffusion_E0_CR );
      Aux_Message( stdout, "  CR_Diffusion_BG_CR        = %14.7e\n",     CR_Diffusion_BG_CR );
      Aux_Message( stdout, "  CR_Diffusion_R02_CR       = %14.7e\n",     CR_Diffusion_R02_CR );
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

   double Dens, MomX, MomY, MomZ, Pres, Eint, Etot;
   double P_cr, CRay;
   double GAMMA_CR_m1_inv = 1.0 / (GAMMA_CR - 1.0);
   
   
   Dens = CR_Diffusion_Rho0;
   MomX = Dens * CR_Diffusion_Vx;
   MomY = Dens * CR_Diffusion_Vy + Dens * 1.e-16 * sin( 8. * M_PI * y);
   MomZ = Dens * CR_Diffusion_Vz;
   Pres = CR_Diffusion_PGas0;
   
   CRay = CR_Diffusion_E0_CR;
   double D1, D2, D3, r;
   D1 = x - CR_Diffusion_CenterX - CR_Diffusion_Vx*Time;
   D2 = y - CR_Diffusion_CenterY - CR_Diffusion_Vy*Time;
   D3 = z - CR_Diffusion_CenterZ - CR_Diffusion_Vz*Time;
   CRay *= EXP( -CR_Diffusion_R02_CR * D1 * D1 );
   //CRay *= EXP( -CR_Diffusion_R02_CR * D2 * D2 );
   //CRay *= EXP( -CR_Diffusion_R02_CR * D3 * D3 );
   r     = (D1 + D2) / sqrt(2.0);
   // CRay *= EXP( -CR_Diffusion_R02_CR * r * r );
   CRay += CR_Diffusion_BG_CR;

// set the output array of passive scaler
#ifdef COSMIC_RAY
   P_cr = CRay * ( GAMMA_CR - 1.0 );
   Pres = Pres + P_cr;
   
   fluid[CRAY] = CRay;
#endif
   
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


#endif // #if ( MODEL == HYDRO )



//-------------------------------------------------------------------------------------------------------
// Function    :  Init_TestProb_Hydro_Cosmic_Ray_AMR
// Description :  Test problem initializer
//
// Note        :  None
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Init_TestProb_Hydro_Cosmic_Ray_AMR()
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

#  endif // #if ( MODEL == HYDRO )

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : Init_TestProb_Hydro_Cosmic_Ray_SoundWave
