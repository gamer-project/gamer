#include "GAMER.h"



// =======================================================================================
// CR_Streaming_Test : test problems (Jiang & Oh 2018, Sec 4)
// =======================================================================================
#define CR_TEST_TRIANGULAR_1D     0     // Sec 4.1.1 : 1D triangular Ec profile             (streaming)
#define CR_TEST_GAUSSIAN_1D       1     // Sec 4.1.1 / 4.1.4 : 1D Gaussian Ec               (streaming or diffusion)
#define CR_TEST_GAUSSIAN_2D       2     // Sec 4.1.2 / 4.1.5 : 2D Gaussian Ec, diagonal B   (streaming or anisotropic diffusion)
#define CR_TEST_CIRCLE_2D         3     // Sec 4.1.5 : 2D ring Ec, circular B               (anisotropic diffusion)
#define CR_TEST_BOTTLENECK_1D     4     // Sec 4.1.3 : 1D density cloud, CR injected from the -x boundary
#define CR_TEST_WAVE_1D           5     // Sec 4.2.1 : 1D CR-driven wave                    (full MHD + CR coupling)
#define CR_TEST_BLAST_2D          6     // Sec 4.2.3 : 2D CR-driven blast wave              (full MHD + CR coupling)
// -- ports of the classic (COSMIC_RAY) tests for comparing CR_E with CRay --
#define CR_TEST_CLASSIC_DIFFUSION 7     // classic CR_Diffusion : Gaussian Ec, uniform B||x
#define CR_TEST_CLASSIC_SHOCKTUBE 8     // classic CR_ShockTube : CR-hydro shock, CR locked to gas (high opacity + CR_SOURCE)
#define CR_TEST_CLASSIC_SOUNDWAVE 9     // classic CR_SoundWave : CR-modified acoustic wave, CR locked to gas
// -- Jiang & Oh 2018 tests (cont.) --
#define CR_TEST_PAPER_SHOCK      10     // Sec 4.2.2 : colliding-flow CR-modified shock (Fig. 11)
#define CR_TEST_BLAST_3D         11     // 3D version of Sec 4.2.3; B direction set by CR_Streaming_B_theta/B_phi


// problem-specific global variables
// =======================================================================================
static int    CR_Streaming_Test;           // test selector (one of the CR_TEST_* macros above)
static int    CR_Streaming_Dir;            // streaming/wave axis for the 1D tests
                                           // (0/1/2 = x/y/z; 3 = (1,1,1) for CR_TEST_CLASSIC_SOUNDWAVE only)
static double CR_Streaming_FlowV;          // uniform flow velocity along CR_Streaming_Dir (Sec 4.1.4)
                                           // --> inflow speed for CR_TEST_PAPER_SHOCK (paper: 10)
static double CR_Streaming_Ec0;            // uniform Ec for CR_TEST_PAPER_SHOCK (paper: 1, 50, or 200)
static int    CR_Streaming_FcInit;         // initial Fc for CR_TEST_PAPER_SHOCK: 1 = advective (4/3)*v*Ec (paper),
                                           // 0 = zero (no flux jump at the collision at t = 0)
static double CR_Streaming_B_theta;        // CR_TEST_BLAST_3D: polar angle of B from +z (deg)
static double CR_Streaming_B_phi;          // CR_TEST_BLAST_3D: azimuthal angle of B from +x (deg)
static int    CR_Streaming_GradOutflowBC;  // gradient-preserving CR outflow BC on the x faces (see GradOutflowBC())
// =======================================================================================


// function prototypes shared within this file
#if ( MODEL == HYDRO  &&  defined CR_TWOMOMENT )
void SetGridIC( real fluid[], const double x, const double y, const double z, const double Time,
                const int lv, double AuxArray[] );
void ShockBC( real Array[], const int ArraySize[], real fluid[], const int NVar_Flu,
              const int GhostSize, const int idx[], const double pos[], const double Time,
              const int lv, const int TFluVarIdxList[], double AuxArray[] );
void GradOutflowBC( real Array[], const int ArraySize[], real fluid[], const int NVar_Flu,
                    const int GhostSize, const int idx[], const double pos[], const double Time,
                    const int lv, const int TFluVarIdxList[], double AuxArray[] );
bool Flag_CR_Streaming( const int i, const int j, const int k, const int lv, const int PID, const double *Threshold );
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

#  ifndef CR_TWOMOMENT
   Aux_Error( ERROR_INFO, "CR_TWOMOMENT must be enabled !!\n" );
#  endif

#  ifndef MHD
   Aux_Error( ERROR_INFO, "MHD must be enabled (CR streaming/diffusion is field-aligned) !!\n" );
#  endif

#  if ( EOS != EOS_GAMMA )
   Aux_Error( ERROR_INFO, "EOS must be EOS_GAMMA !!\n" );
#  endif

#  ifdef CR_DIFFUSION
   Aux_Error( ERROR_INFO, "CR_DIFFUSION must be disabled !!\n" );
#  endif


// warnings


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Validating test problem %d ... done\n", TESTPROB_ID );

} // FUNCTION : Validate



#if ( MODEL == HYDRO  &&  defined CR_TWOMOMENT )
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
// --> note that VARIABLE, DEFAULT, MIN, and MAX must have the same data type
// --> some handy constants (e.g., Useless_bool, Eps_double, NoMin_int, ...) are defined in "include/ReadPara.h"
// *****************************************************************************************************************************
// ReadPara->Add( "KEY_IN_THE_FILE",            &VARIABLE,                   DEFAULT,     MIN,              MAX               );
// *****************************************************************************************************************************
   ReadPara->Add( "CR_Streaming_Test",          &CR_Streaming_Test,          0,           0,                11                );
   ReadPara->Add( "CR_Streaming_Dir",           &CR_Streaming_Dir,           0,           0,                3                 );
   ReadPara->Add( "CR_Streaming_FlowV",         &CR_Streaming_FlowV,         0.0,         NoMin_double,     NoMax_double      );
   ReadPara->Add( "CR_Streaming_Ec0",           &CR_Streaming_Ec0,           1.0,         Eps_double,       NoMax_double      );
   ReadPara->Add( "CR_Streaming_FcInit",        &CR_Streaming_FcInit,        1,           0,                1                 );
   ReadPara->Add( "CR_Streaming_GradOutflowBC", &CR_Streaming_GradOutflowBC, 0,           0,                1                 );
   ReadPara->Add( "CR_Streaming_B_theta",       &CR_Streaming_B_theta,       90.0,        NoMin_double,     NoMax_double      );
   ReadPara->Add( "CR_Streaming_B_phi",         &CR_Streaming_B_phi,         0.0,         NoMin_double,     NoMax_double      );

   ReadPara->Read( FileName );

   delete ReadPara;

// (1-2) set the default values

// (1-3) check the runtime parameters
   if ( CR_Streaming_Test != CR_TEST_BLAST_3D  &&
        ( CR_Streaming_B_theta != 90.0  ||  CR_Streaming_B_phi != 0.0 ) )
      Aux_Error( ERROR_INFO, "CR_Streaming_B_theta/B_phi are only supported for CR_Streaming_Test = %d (3D blast) !!\n",
                 CR_TEST_BLAST_3D );

// TODO: Reject CR_Streaming_Dir = 3 except for CR_TEST_CLASSIC_SOUNDWAVE (SetGridIC() reads amr->BoxCenter[3])
// TODO: Reject CR_Streaming_GradOutflowBC with CR_TEST_BOTTLENECK_1D/PAPER_SHOCK (overrides their BC_User_Ptr)


// (2) set the problem-specific derived parameters
   const bool is_1D = ( CR_Streaming_Test == CR_TEST_TRIANGULAR_1D     ||
                        CR_Streaming_Test == CR_TEST_GAUSSIAN_1D       ||
                        CR_Streaming_Test == CR_TEST_BOTTLENECK_1D     ||
                        CR_Streaming_Test == CR_TEST_WAVE_1D           ||
                        CR_Streaming_Test == CR_TEST_CLASSIC_DIFFUSION ||
                        CR_Streaming_Test == CR_TEST_CLASSIC_SHOCKTUBE ||
                        CR_Streaming_Test == CR_TEST_PAPER_SHOCK         );


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

// for the 1D tests, output a line along the streaming direction
   if ( is_1D )
   {
      const int target = ( CR_Streaming_Dir == 0 ) ? OUTPUT_X : ( CR_Streaming_Dir == 1 ) ? OUTPUT_Y : OUTPUT_Z;
      if ( OPT__OUTPUT_PART != target )
      {
         OPT__OUTPUT_PART = target;
         PRINT_RESET_PARA( OPT__OUTPUT_PART, FORMAT_INT, "" );
      }
   }

// CR_Streaming_GradOutflowBC requires the user BC on the x faces
   if ( CR_Streaming_GradOutflowBC )
   {
      if ( OPT__BC_FLU[0] != BC_FLU_USER )
      {
         OPT__BC_FLU[0] = BC_FLU_USER;
         PRINT_RESET_PARA( OPT__BC_FLU[0], FORMAT_INT, "(gradient-preserving CR outflow BC)" );
      }
      if ( OPT__BC_FLU[1] != BC_FLU_USER )
      {
         OPT__BC_FLU[1] = BC_FLU_USER;
         PRINT_RESET_PARA( OPT__BC_FLU[1], FORMAT_INT, "(gradient-preserving CR outflow BC)" );
      }
   }


// (4) make a note
   if ( MPI_Rank == 0 )
   {
      Aux_Message( stdout, "=============================================================================\n" );
      Aux_Message( stdout, "  test problem ID       = %d\n",     TESTPROB_ID );
      Aux_Message( stdout, "  CR_Streaming_Test     = %d\n",     CR_Streaming_Test  );
      Aux_Message( stdout, "  CR_Streaming_Dir      = %d\n",     CR_Streaming_Dir   );
      Aux_Message( stdout, "  CR_Streaming_FlowV    = %14.7e\n", CR_Streaming_FlowV );
      Aux_Message( stdout, "  CR_Streaming_Ec0      = %14.7e\n", CR_Streaming_Ec0   );
      Aux_Message( stdout, "  CR_Streaming_FcInit   = %d\n",     CR_Streaming_FcInit );
      Aux_Message( stdout, "  CR_Streaming_GradOutflowBC = %d\n", CR_Streaming_GradOutflowBC );
      Aux_Message( stdout, "  CR_Streaming_B_theta  = %14.7e\n", CR_Streaming_B_theta );
      Aux_Message( stdout, "  CR_Streaming_B_phi    = %14.7e\n", CR_Streaming_B_phi   );
      Aux_Message( stdout, "=============================================================================\n" );
   }


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Setting runtime parameters ... done\n" );

} // FUNCTION : SetParameter



//-------------------------------------------------------------------------------------------------------
// Function    :  CR_Bottleneck_Dens
// Description :  Density profile of the cold cloud used in the bottleneck test (Sec 4.1.3, Eq. 24)
//
// Note        :  1. Assume the paper domain x in [0,1000] with the cloud at x0 = 200
//
// Parameter   :  x : x coordinate
//
// Return      :  Gas density
//-------------------------------------------------------------------------------------------------------
static double CR_Bottleneck_Dens( const double x )
{
   const double rho_c = 1.0, rho_h = 0.1, x0 = 200.0, dx0 = 25.0;
   const double t1 = 1.0 + std::tanh( (x  - x0)/dx0 );
   const double t2 = 1.0 + std::tanh( (x0 - x )/dx0 );
   return rho_h + (rho_c - rho_h)*t1*t2;
} // FUNCTION : CR_Bottleneck_Dens



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
//                5. ADV_* are overwritten by CR_UpdateOpacity() before use
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

   const double xc = amr->BoxCenter[0];
   const double yc = amr->BoxCenter[1];

   double Dens = 1.0, vx = 0.0, vy = 0.0, vz = 0.0, Pgas = 1.0, cr_E = 0.0;
   double cr_F1 = 0.0;   // reduced CR flux Fc/CR_VMAX along x

   switch ( CR_Streaming_Test )
   {
      case CR_TEST_TRIANGULAR_1D :
      {
//       Sec 4.1.1: Ec = 2 - |r| for |r| < 1, else 1   (r relative to the box center)
         const double r = ( (CR_Streaming_Dir==0)?x : (CR_Streaming_Dir==1)?y : z ) - amr->BoxCenter[CR_Streaming_Dir];
         const double d = std::fabs( r );
         cr_E = ( d < 1.0 ) ? (2.0 - d) : 1.0;
         vx   = ( CR_Streaming_Dir==0 ) ? CR_Streaming_FlowV : 0.0;
         vy   = ( CR_Streaming_Dir==1 ) ? CR_Streaming_FlowV : 0.0;
         vz   = ( CR_Streaming_Dir==2 ) ? CR_Streaming_FlowV : 0.0;
         break;
      }

      case CR_TEST_GAUSSIAN_1D :
      {
//       Sec 4.1.1 / 4.1.4: Ec = exp(-40*r^2)   (r relative to the box center)
         const double r = ( (CR_Streaming_Dir==0)?x : (CR_Streaming_Dir==1)?y : z ) - amr->BoxCenter[CR_Streaming_Dir];
         cr_E = std::exp( -40.0 * r*r );
         vx   = ( CR_Streaming_Dir==0 ) ? CR_Streaming_FlowV : 0.0;
         vy   = ( CR_Streaming_Dir==1 ) ? CR_Streaming_FlowV : 0.0;
         vz   = ( CR_Streaming_Dir==2 ) ? CR_Streaming_FlowV : 0.0;
         break;
      }

      case CR_TEST_GAUSSIAN_2D :
      {
//       Sec 4.1.2 / 4.1.5: Ec = exp(-40*(dx^2 + dy^2)), B along the diagonal (set in SetBFieldIC())
         const double dx = x - xc, dy = y - yc;
         cr_E = std::exp( -40.0 * (dx*dx + dy*dy) );
         break;
      }

      case CR_TEST_CIRCLE_2D :
      {
//       Sec 4.1.5: a ring of Ec; B is circular (set in SetBFieldIC())
//       Ec = 12 in 0.5 < r < 0.7 and |phi| < pi/12, else 10 (analytic solution: Eq. 28)
//       --> rho = 1e8 emulates the paper's fixed MHD variables under the unbalanced hoop stress of the circular B
//           (rho is unused by the CR module without streaming)
         const double dx = x - xc, dy = y - yc;
         const double r   = std::sqrt( dx*dx + dy*dy );
         const double phi = std::atan2( dy, dx );
         cr_E = ( r > 0.5  &&  r < 0.7  &&  std::fabs(phi) < M_PI/12.0 ) ? 12.0 : 10.0;
         Dens = 1.0e8;
         break;
      }

      case CR_TEST_BOTTLENECK_1D :
      {
//       Sec 4.1.3: a cold dense cloud; CRs are injected from the -x boundary by BottleneckBC()
         Dens = CR_Bottleneck_Dens( x );
         cr_E = 1.0e-6;
         break;
      }

      case CR_TEST_WAVE_1D :
      {
//       Sec 4.2.1: Ec = 20 + 10*sin(pi*(x-xc)); uniform gas; full MHD + CR coupling (CR_SOURCE = 1)
         cr_E = 20.0 + 10.0*std::sin( M_PI*(x - xc) );
         break;
      }

      case CR_TEST_BLAST_2D :
      {
//       Sec 4.2.3: Ec = 100 inside r < 0.02, else 0.1; uniform background (rho = 1, Eint = 2.5)
         const double dx = x - xc, dy = y - yc;
         const double r  = std::sqrt( dx*dx + dy*dy );
         cr_E = ( r < 0.02 ) ? 100.0 : 0.1;
         Pgas = (GAMMA - 1.0)*2.5;
         break;
      }

      case CR_TEST_BLAST_3D :
      {
//       3D version of Sec 4.2.3: Ec = 100 inside the sphere r < 0.02, else 0.1; uniform background (rho = 1, Eint = 2.5)
//       --> B direction set by CR_Streaming_B_theta/B_phi (see SetBFieldIC())
         const double dx = x - xc, dy = y - yc, dz = z - amr->BoxCenter[2];
         const double r  = std::sqrt( dx*dx + dy*dy + dz*dz );
         cr_E = ( r < 0.02 ) ? 100.0 : 0.1;
         Pgas = (GAMMA - 1.0)*2.5;
         break;
      }

      case CR_TEST_CLASSIC_DIFFUSION :
      {
//       classic CR_Diffusion (Type 0, Mag_Type 0): Ec = E0*exp(-R02*r^2) + BG, uniform B||x, rho = 1, Pgas = 5/3
//       --> set CR_SIGMA = 1/(3*kappa_classic) and a large CR_SIGMA_PERP for kappa_perp = 0
//       --> CR_SOURCE = 1: live gas as in the classic module; 0: static gas (pure diffusion, Eq. 25)
         const double E0 = 1.0, R02 = 40.0, BG = 0.1;
         const double dxc = x - xc, dyc = y - yc, dzc = z - amr->BoxCenter[2];
         cr_E = E0*std::exp( -R02*( dxc*dxc + dyc*dyc + dzc*dzc ) ) + BG;
         Pgas = 1.666666666667;                  // classic CR_Diffusion_PGas0
         break;
      }

      case CR_TEST_CLASSIC_SHOCKTUBE :
      {
//       classic CR_ShockTube (gamma_cr = 4/3, no diffusion): Ec = Pcr/(gamma_cr-1) = 3*Pcr
//       --> CR locked to the gas by a large CR_SIGMA with CR_STREAM = 0 and CR_SOURCE = 1
         const double Pcr_L = 1.3e5, Pcr_R = 2.4e2;
         const double _gcrm1 = 3.0;   // = 1/(gamma_cr-1) with gamma_cr = 4/3
         if      ( x < xc ) { Dens = 1.0;            Pgas = 6.7e4;             cr_E = Pcr_L*_gcrm1;             }
         else if ( x > xc ) { Dens = 0.2;            Pgas = 2.4e2;             cr_E = Pcr_R*_gcrm1;             }
         else               { Dens = 0.5*(1.0+0.2); Pgas = 0.5*(6.7e4+2.4e2); cr_E = 0.5*(Pcr_L+Pcr_R)*_gcrm1; }
         break;
      }

      case CR_TEST_CLASSIC_SOUNDWAVE :
      {
//       classic CR_SoundWave (gamma_cr = 4/3): cs^2 = (gamma*Pg0 + gamma_cr*Pcr0)/rho0
//       --> same CR settings as CR_TEST_CLASSIC_SHOCKTUBE
//       --> CR_Streaming_Dir = 0/1/2: along x/y/z; 3: along (1,1,1)
         const double Rho0 = 1.0, Pg0 = 1.0, Pcr0 = 1.0, Delta = 1.0e-6, Sign = 1.0;
         const double gamma_cr  = 4.0/3.0;
         const double cs       = std::sqrt( GAMMA*Pg0 + gamma_cr*Pcr0 );   // rho0 = 1
         const double delta_cs = Delta/cs;
         const bool    diag    = ( CR_Streaming_Dir == 3 );
         const double WaveL    = diag ? amr->BoxSize[0]/std::sqrt(3.0) : amr->BoxSize[CR_Streaming_Dir];
         const double WaveK    = 2.0*M_PI/WaveL;
         const double r        = diag ? ( x + y + z )/std::sqrt(3.0)
                                      : ( (CR_Streaming_Dir==0)?x : (CR_Streaming_Dir==1)?y : z );
         const double wave     = std::sin( WaveK*r );                      // t = 0
         const double vel      = Sign*Delta*wave;                          // velocity perturbation
         Dens = ( 1.0 + delta_cs*wave )*Rho0;
         if ( diag ) {
            vx = vy = vz = vel/std::sqrt(3.0);
         } else {
            vx = ( CR_Streaming_Dir==0 ) ? vel : 0.0;
            vy = ( CR_Streaming_Dir==1 ) ? vel : 0.0;
            vz = ( CR_Streaming_Dir==2 ) ? vel : 0.0;
         }
         Pgas = ( 1.0 + delta_cs*wave*GAMMA    )*Pg0;
         cr_E = ( 1.0 + delta_cs*wave*gamma_cr )*Pcr0/(gamma_cr - 1.0);
         break;
      }

      case CR_TEST_PAPER_SHOCK :
      {
//       Sec 4.2.2 (Fig. 11): gas streams at +/-FlowV collide at the box center and drive a CR-modified shock to each side
//       --> rho = Pgas = 1, Ec = CR_Streaming_Ec0, B||x (v_A = 1 upstream); x ghost zones fixed by ShockBC()
//       --> paper: FlowV = 10, Ec0 = 1/50/200, CR_SIGMA = 10, CR_VMAX = 100
//       --> CR_F1 stores the reduced flux Fc/CR_VMAX, so the advective Fc = (4/3)*v*Ec is divided by CR_VMAX
         const double v = ( x < xc ) ? CR_Streaming_FlowV : -CR_Streaming_FlowV;
         Dens  = 1.0;
         Pgas  = 1.0;
         vx    = v;
         cr_E  = CR_Streaming_Ec0;
         cr_F1 = ( CR_Streaming_FcInit == 0 ) ? 0.0 : (4.0/3.0)*v*cr_E/CR_VMAX;
         break;
      }

      default :
         Aux_Error( ERROR_INFO, "unsupported CR_Streaming_Test (%d) !!\n", CR_Streaming_Test );
   } // switch ( CR_Streaming_Test )


// momentum
   const double MomX = Dens*vx;
   const double MomY = Dens*vy;
   const double MomZ = Dens*vz;
   const double Pres = Pgas;  // gas pressure only

   const double Eint = EoS_DensPres2Eint_CPUPtr( Dens, Pres, fluid+NCOMP_FLUID, EoS_AuxArray_Flt,
                                                 EoS_AuxArray_Int, h_EoS_Table );
   const double Etot = Hydro_ConEint2Etot( Dens, MomX, MomY, MomZ, Eint, 0.0 );   // do NOT include magnetic energy here

// CR fields
   fluid[CR_E ] = cr_E;
   fluid[CR_F1] = cr_F1;
   fluid[CR_F2] = 0.0;
   fluid[CR_F3] = 0.0;

// streaming opacity and velocity (see Note 5)
   fluid[ADV_SIGMA] = CR_MAX_OPACITY;
   fluid[ADV_VX   ] = 0.0;
   fluid[ADV_VY   ] = 0.0;
   fluid[ADV_VZ   ] = 0.0;

// gas fields
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

   switch ( CR_Streaming_Test )
   {
      case CR_TEST_TRIANGULAR_1D :
      case CR_TEST_GAUSSIAN_1D :
//       uniform field along the streaming direction, |B| = 1  --> v_A = 1
         magnetic[MAGX] = ( CR_Streaming_Dir==0 ) ? 1.0 : 0.0;
         magnetic[MAGY] = ( CR_Streaming_Dir==1 ) ? 1.0 : 0.0;
         magnetic[MAGZ] = ( CR_Streaming_Dir==2 ) ? 1.0 : 0.0;
         break;

      case CR_TEST_GAUSSIAN_2D :
//       uniform field along the x-y diagonal, |B| = 1  --> v_A = 1   (Sec 4.1.2 / 4.1.5)
         magnetic[MAGX] = 1.0/std::sqrt(2.0);
         magnetic[MAGY] = 1.0/std::sqrt(2.0);
         magnetic[MAGZ] = 0.0;
         break;

      case CR_TEST_CIRCLE_2D :
      {
//       circular field B = ( -(y-yc), (x-xc) ) / r, |B| = 1   (Sec 4.1.5)
//       --> set as the discrete curl of Az = -r so that div(B) = 0 to round-off, as in Athena++ src/pgen/cr_diffusion.cpp
         const double xc = amr->BoxCenter[0];
         const double yc = amr->BoxCenter[1];
         const double dh = amr->dh[lv];                 // each B component is set at its own face center
#        define CR_CIRCLE_AZ( X, Y )  ( -std::sqrt( SQR((X)-xc) + SQR((Y)-yc) ) )
         magnetic[MAGX] =  ( CR_CIRCLE_AZ( x, y+0.5*dh ) - CR_CIRCLE_AZ( x, y-0.5*dh ) ) / dh;
         magnetic[MAGY] = -( CR_CIRCLE_AZ( x+0.5*dh, y ) - CR_CIRCLE_AZ( x-0.5*dh, y ) ) / dh;
         magnetic[MAGZ] = 0.0;
#        undef CR_CIRCLE_AZ
         break;
      }

      case CR_TEST_CLASSIC_SOUNDWAVE :
//       B parallel to the wave vector (no perpendicular CR terms)
         if ( CR_Streaming_Dir == 3 ) {
            magnetic[MAGX] = 1.0/std::sqrt(3.0);
            magnetic[MAGY] = 1.0/std::sqrt(3.0);
            magnetic[MAGZ] = 1.0/std::sqrt(3.0);
         } else {
            magnetic[MAGX] = ( CR_Streaming_Dir==0 ) ? 1.0 : 0.0;
            magnetic[MAGY] = ( CR_Streaming_Dir==1 ) ? 1.0 : 0.0;
            magnetic[MAGZ] = ( CR_Streaming_Dir==2 ) ? 1.0 : 0.0;
         }
         break;

      case CR_TEST_BLAST_3D :
      {
//       uniform field along (theta, phi) in degrees, |B| = 1  --> v_A = 1  (rho = 1)
//       --> the default (90, 0) recovers B||x of CR_TEST_BLAST_2D
         const double theta = CR_Streaming_B_theta * M_PI/180.0;
         const double phi   = CR_Streaming_B_phi   * M_PI/180.0;
         magnetic[MAGX] = std::sin(theta)*std::cos(phi);
         magnetic[MAGY] = std::sin(theta)*std::sin(phi);
         magnetic[MAGZ] = std::cos(theta);
         break;
      }

      case CR_TEST_BOTTLENECK_1D :
      case CR_TEST_WAVE_1D :
      case CR_TEST_BLAST_2D :
      case CR_TEST_CLASSIC_DIFFUSION :
      case CR_TEST_CLASSIC_SHOCKTUBE :
      case CR_TEST_PAPER_SHOCK :
      default :
//       uniform field along x, |B| = 1
         magnetic[MAGX] = 1.0;
         magnetic[MAGY] = 0.0;
         magnetic[MAGZ] = 0.0;
         break;
   } // switch ( CR_Streaming_Test )

} // FUNCTION : SetBFieldIC
#endif // #ifdef MHD



//-------------------------------------------------------------------------------------------------------
// Function    :  BottleneckBC
// Description :  User boundary condition for the bottleneck test (Sec 4.1.3): inject CRs from the
//                -x boundary by fixing Ec = 3 while keeping the background gas profile
//
// Note        :  1. Linked to the function pointer "BC_User_Ptr"
//                2. Only for the -x face (OPT__BC_FLU_XM = 4)
//                3. Ec_ghost = 3 and Fc_ghost = -Fc_active (reflecting)
//                4. Fc_active is read from "Array", whose interior is already filled (as in JetICMWall)
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
//                TFluVarIdxList : List recording the target fluid variable indices
//                AuxArray       : Auxiliary array
//
// Return      :  fluid
//
// Reference   :  Athena++ src/pgen/cr_diffusion.cpp BottleneckCRInnerX1()
//-------------------------------------------------------------------------------------------------------
void BottleneckBC( real Array[], const int ArraySize[], real fluid[], const int NVar_Flu,
                   const int GhostSize, const int idx[], const double pos[], const double Time,
                   const int lv, const int TFluVarIdxList[], double AuxArray[] )
{

// start from the IC
   SetGridIC( fluid, pos[0], pos[1], pos[2], Time, lv, AuxArray );

// fixed Ec reservoir and zero transverse CR fluxes
   fluid[CR_E ] = 3.0;
   fluid[CR_F2] = 0.0;
   fluid[CR_F3] = 0.0;

// reflecting CR flux: Fc_ghost = -Fc_active
// --> slot v of "Array" stores the field TFluVarIdxList[v]
   typedef real (*vla)[ ArraySize[2] ][ ArraySize[1] ][ ArraySize[0] ];
   vla Array3D = ( vla )Array;

   const int i_ref = GhostSize;   // active cell adjacent to the -x boundary
   const int jg    = idx[1];
   const int kg    = idx[2];

   real CRF1_active = 0.0;
   for (int v=0; v<NVar_Flu; v++)
      if ( TFluVarIdxList[v] == CR_F1 )   CRF1_active = Array3D[v][kg][jg][i_ref];

   fluid[CR_F1] = -CRF1_active;

} // FUNCTION : BottleneckBC



//-------------------------------------------------------------------------------------------------------
// Function    :  ShockBC
// Description :  User boundary condition for CR_TEST_PAPER_SHOCK: fix the x ghost zones to the initial
//                condition (Jiang & Oh 2018, Sec 4.2.2)
//
// Note        :  1. Linked to the function pointer "BC_User_Ptr" for both the -x and +x faces (OPT__BC_FLU_XM/XP = 4)
//                2. Valid until the shocks reach the boundaries
//
// Parameter   :  See BottleneckBC()
//
// Return      :  fluid
//-------------------------------------------------------------------------------------------------------
void ShockBC( real Array[], const int ArraySize[], real fluid[], const int NVar_Flu,
              const int GhostSize, const int idx[], const double pos[], const double Time,
              const int lv, const int TFluVarIdxList[], double AuxArray[] )
{

// fixed ghost zones = initial condition
   SetGridIC( fluid, pos[0], pos[1], pos[2], Time, lv, AuxArray );

} // FUNCTION : ShockBC



//-------------------------------------------------------------------------------------------------------
// Function    :  GradOutflowBC
// Description :  Gradient-preserving CR outflow boundary condition on the -x and +x faces (streaming tests)
//
// Note        :  1. Linked to the function pointer "BC_User_Ptr" for both the -x and +x faces
//                   --> OPT__BC_FLU_XM/XP are reset in SetParameter()
//                2. Gas fields are set to the initial condition
//                   --> valid only for a uniform gas decoupled from CRs (CR_SOURCE = 0)
//                3. Ec is linearly extrapolated from the last two active cells and Fc is copied from the last active cell
//                   --> keeps grad(Pc) at the boundary so that CRs stream out at (4/3)*v_A*Ec (Jiang & Oh 2018, Fig. 3),
//                       whereas copying Ec (Athena++ outflow BC) flattens Ec near the boundary
//
// Parameter   :  See BottleneckBC()
//
// Return      :  fluid
//-------------------------------------------------------------------------------------------------------
void GradOutflowBC( real Array[], const int ArraySize[], real fluid[], const int NVar_Flu,
                    const int GhostSize, const int idx[], const double pos[], const double Time,
                    const int lv, const int TFluVarIdxList[], double AuxArray[] )
{

// start from the IC
   SetGridIC( fluid, pos[0], pos[1], pos[2], Time, lv, AuxArray );

// map the prepared 1D "Array" to 4D [NVar_Flu][k][j][i]
   typedef real (*vla)[ ArraySize[2] ][ ArraySize[1] ][ ArraySize[0] ];
   vla Array3D = ( vla )Array;

   const int ig = idx[0];
   const int jg = idx[1];
   const int kg = idx[2];
   const int Nx = ArraySize[0];

// identify the x-face and the two nearest active cells along the outward normal
   int i_last, i_lastm1, ndepth;
   if ( ig < GhostSize )                     // -x ghost zone
   {
      i_last   = GhostSize;                  // last active cell (adjacent to the -x boundary)
      i_lastm1 = GhostSize + 1;              // one cell further into the interior
      ndepth   = GhostSize - ig;             // 1, 2, ... outward from the boundary
   }
   else                                      // +x ghost zone ( ig >= Nx - GhostSize )
   {
      i_last   = Nx - GhostSize - 1;         // last active cell (adjacent to the +x boundary)
      i_lastm1 = Nx - GhostSize - 2;
      ndepth   = ig - i_last;                // 1, 2, ...
   }

// read Ec (two active cells) and Fc1/2/3 (last active cell) from the prepared interior
   real Ec_last = 0.0, Ec_lastm1 = 0.0, Fc1 = 0.0, Fc2 = 0.0, Fc3 = 0.0;
   for (int v=0; v<NVar_Flu; v++)
   {
      const int f = TFluVarIdxList[v];
      if      ( f == CR_E  ) { Ec_last = Array3D[v][kg][jg][i_last]; Ec_lastm1 = Array3D[v][kg][jg][i_lastm1]; }
      else if ( f == CR_F1 )   Fc1     = Array3D[v][kg][jg][i_last];
      else if ( f == CR_F2 )   Fc2     = Array3D[v][kg][jg][i_last];
      else if ( f == CR_F3 )   Fc3     = Array3D[v][kg][jg][i_last];
   }

// gradient-preserving Ec (linear extrapolation, positive-clamped) + outflow-copy Fc
   const real CR_EC_TINY = (real)1.0e-30;
   const real ec_extrap  = Ec_last + (real)ndepth*( Ec_last - Ec_lastm1 );
   fluid[CR_E ] = FMAX( ec_extrap, CR_EC_TINY );
   fluid[CR_F1] = Fc1;
   fluid[CR_F2] = Fc2;
   fluid[CR_F3] = Fc3;

} // FUNCTION : GradOutflowBC
#endif // #if ( MODEL == HYDRO  &&  defined CR_TWOMOMENT )



//-------------------------------------------------------------------------------------------------------
// Function    :  Init_TestProb_Hydro_CR_Streaming
// Description :  Test problem initializer
//
// Note        :  None
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Init_TestProb_Hydro_CR_Streaming()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );


// validate the compilation flags and runtime parameters
   Validate();


#  if ( MODEL == HYDRO  &&  defined CR_TWOMOMENT )
// set the problem-specific runtime parameters
   SetParameter();


// procedure to enable a problem-specific function:
// 1. define a user-specified function (example functions are given below)
// 2. declare its function prototype on the top of this file
// 3. set the corresponding function pointer below to the new problem-specific function
// 4. enable the corresponding runtime option in "Input__Parameter"
//    --> for instance, enable OPT__OUTPUT_USER for Output_User_Ptr
   Init_Function_User_Ptr            = SetGridIC;
#  ifdef MHD
   Init_Function_BField_User_Ptr     = SetBFieldIC;
#  endif

// spherical-shell refinement (see Flag_CR_Streaming())
   Flag_User_Ptr                     = Flag_CR_Streaming;

// the bottleneck test injects CRs through a user boundary condition on the -x face
   if ( CR_Streaming_Test == CR_TEST_BOTTLENECK_1D )
   {
      BC_User_Ptr                    = BottleneckBC;
//    uniform B --> reuse SetBFieldIC() for the B ghost zones
#     ifdef MHD
      BC_BField_User_Ptr             = SetBFieldIC;
#     endif
   }

// the Sec 4.2.2 shock test fixes the x ghost zones to the initial (inflow) state
   if ( CR_Streaming_Test == CR_TEST_PAPER_SHOCK )
   {
      BC_User_Ptr                    = ShockBC;
//    uniform B --> reuse SetBFieldIC() for the B ghost zones
#     ifdef MHD
      BC_BField_User_Ptr             = SetBFieldIC;
#     endif
   }

// gradient-preserving CR outflow BC on the x faces
   if ( CR_Streaming_GradOutflowBC )
   {
      BC_User_Ptr                    = GradOutflowBC;
//    assume uniform B --> reuse SetBFieldIC() for the B ghost zones
#     ifdef MHD
      BC_BField_User_Ptr             = SetBFieldIC;
#     endif
   }
#  endif // #if ( MODEL == HYDRO  &&  defined CR_TWOMOMENT )


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : Init_TestProb_Hydro_CR_Streaming
