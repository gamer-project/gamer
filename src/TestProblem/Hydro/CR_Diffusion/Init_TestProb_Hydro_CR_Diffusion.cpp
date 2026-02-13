#include "GAMER.h"



// problem-specific global variables
// =======================================================================================
static double CR_Diffusion_CenterX;     // Center location on x-axis
static double CR_Diffusion_CenterY;     // Center location on y-axis
static double CR_Diffusion_CenterZ;     // Center location on z-axis
static double CR_Diffusion_Vx;          // Background x-velocity
static double CR_Diffusion_Vy;          // Background y-velocity
static double CR_Diffusion_Vz;          // Background z-velocity
static double CR_Diffusion_Rho0;        // Background density
static double CR_Diffusion_PGas0;       // Background gas pressure
static double CR_Diffusion_E0_CR;       // Amplitude of cosmic ray energy density
static double CR_Diffusion_BG_CR;       // Background cosmic ray energy density
static int    CR_Diffusion_Type;        // The initial condition type:
                                        // 0: Gaussian distribution ball CR
                                        // 1: Step function ring distribution CR
                                        // 2: Gaussian distribution ring CR
                                        // 3: Gaussian distribution 1D
                                        // 4: CR blast wave
static int    CR_Diffusion_Mag_Type;    // The magnetic field type:
                                        // 0: Uniform
                                        // 1: Circular
                                        // 2: Random within (-1.0,+1.0) (Not divergence free)
                                        // 3: Radiative (Not divergence free)
                                        // 4: Suwa+ 2007
static double CR_Diffusion_MagX;        // Magnitude of x component
static double CR_Diffusion_MagY;        // Magnitude of y component
static double CR_Diffusion_MagZ;        // Magnitude of z component
// simulation dimension
static int    CR_Diffusion_GX;          // Include x direction or not [0/1]
static int    CR_Diffusion_GY;          // Include y direction or not [0/1]
static int    CR_Diffusion_GZ;          // Include z direction or not [0/1]
// random magnetic field
static int    CR_Diffusion_Seed;        // Random seed of magnetic field
static RandomNumber_t *RNG = NULL;
// gaussian ball
static double CR_Diffusion_R02_CR;      // Inverse of variance of the Gaussian distribution
// step function ring
static double CR_Diffusion_R_In;        // Inner radius
static double CR_Diffusion_R_Out;       // Outer radius
// gaussian ring
static double CR_Diffusion_CenterR;     // Center radius
static double CR_Diffusion_delR;        // Standard deviation of the Gaussian distribution on radial direction
static double CR_Diffusion_delPhi;      // Standard deviation of the Gaussian distribution on azimuthal direction
// blast wave
static double CR_Diffusion_R0_CR;       // CR source radius
static double CR_Diffusion_R0_B;        // Magnetic filed characteristic length. See Suwa+ 2007.


static int    CR_Diffusion_Space;       // Simulation space
static int    CR_Diffusion_Dim;         // Simulation dimensions
// =======================================================================================


// problem-specific function prototypes
static void OutputError();
static double RandomNumber( RandomNumber_t *RNG, const double Min, const double Max );




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

#  ifndef COSMIC_RAY
   Aux_Error( ERROR_INFO, "COSMIC_RAY must be enabled !!\n" );
#  endif

#  ifndef CR_DIFFUSION
   Aux_Error( ERROR_INFO, "CR_DIFFUSION must be enabled !!\n" );
#  endif

#  if ( EOS != EOS_COSMIC_RAY )
   Aux_Error( ERROR_INFO, "EOS != EOS_COSMIC_RAY !!\n" );
#  endif


// warnings


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Validating test problem %d ... done\n", TESTPROB_ID );

} // FUNCTION : Validate



#if ( MODEL == HYDRO  &&  defined CR_DIFFUSION )
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
// *****************************************************************************************************************************************
// LOAD_PARA( load_mode, "KEY_IN_THE_FILE",           &VARIABLE,                       DEFAULT,       MIN,              MAX               );
// *****************************************************************************************************************************************
   LOAD_PARA( load_mode, "CR_Diffusion_CenterX",      &CR_Diffusion_CenterX,          -1.0,           NoMin_double,     amr->BoxSize[0]   );
   LOAD_PARA( load_mode, "CR_Diffusion_CenterY",      &CR_Diffusion_CenterY,          -1.0,           NoMin_double,     amr->BoxSize[1]   );
   LOAD_PARA( load_mode, "CR_Diffusion_CenterZ",      &CR_Diffusion_CenterZ,          -1.0,           NoMin_double,     amr->BoxSize[2]   );
   LOAD_PARA( load_mode, "CR_Diffusion_Vx",           &CR_Diffusion_Vx,                0.0,           NoMin_double,     NoMax_double      );
   LOAD_PARA( load_mode, "CR_Diffusion_Vy",           &CR_Diffusion_Vy,                0.0,           NoMin_double,     NoMax_double      );
   LOAD_PARA( load_mode, "CR_Diffusion_Vz",           &CR_Diffusion_Vz,                0.0,           NoMin_double,     NoMax_double      );
   LOAD_PARA( load_mode, "CR_Diffusion_Rho0",         &CR_Diffusion_Rho0,              1.0,           Eps_double,       NoMax_double      );
   LOAD_PARA( load_mode, "CR_Diffusion_PGas0",        &CR_Diffusion_PGas0,             1.0,           Eps_double,       NoMax_double      );
   LOAD_PARA( load_mode, "CR_Diffusion_E0_CR",        &CR_Diffusion_E0_CR,             0.1,           0.0,              NoMax_double      );
   LOAD_PARA( load_mode, "CR_Diffusion_BG_CR",        &CR_Diffusion_BG_CR,             0.1,           0.0,              NoMax_double      );
   LOAD_PARA( load_mode, "CR_Diffusion_R02_CR",       &CR_Diffusion_R02_CR,            0.05,          0.0,              NoMax_double      );
   LOAD_PARA( load_mode, "CR_Diffusion_Type",         &CR_Diffusion_Type,              0,             0,                4                 );
   LOAD_PARA( load_mode, "CR_Diffusion_Mag_Type",     &CR_Diffusion_Mag_Type,          0,             0,                4                 );
   LOAD_PARA( load_mode, "CR_Diffusion_MagX",         &CR_Diffusion_MagX,              0.0,           NoMin_double,     NoMax_double      );
   LOAD_PARA( load_mode, "CR_Diffusion_MagY",         &CR_Diffusion_MagY,              0.0,           NoMin_double,     NoMax_double      );
   LOAD_PARA( load_mode, "CR_Diffusion_MagZ",         &CR_Diffusion_MagZ,              0.0,           NoMin_double,     NoMax_double      );
   LOAD_PARA( load_mode, "CR_Diffusion_Seed",         &CR_Diffusion_Seed,              0,             NoMin_int,        NoMax_int         );
   LOAD_PARA( load_mode, "CR_Diffusion_R_In",         &CR_Diffusion_R_In,              0.5,           0.0,              NoMax_double      );
   LOAD_PARA( load_mode, "CR_Diffusion_R_Out",        &CR_Diffusion_R_Out,             0.7,           0.0,              NoMax_double      );
   LOAD_PARA( load_mode, "CR_Diffusion_CenterR",      &CR_Diffusion_CenterR,           0.6,           0.0,              NoMax_double      );
   LOAD_PARA( load_mode, "CR_Diffusion_delR",         &CR_Diffusion_delR,              0.05,          0.0,              NoMax_double      );
   LOAD_PARA( load_mode, "CR_Diffusion_delPhi",       &CR_Diffusion_delPhi,            0.5,           0.0,              NoMax_double      );
   LOAD_PARA( load_mode, "CR_Diffusion_R0_CR",        &CR_Diffusion_R0_CR,             0.02,          0.0,              NoMax_double      );
   LOAD_PARA( load_mode, "CR_Diffusion_R0_B",         &CR_Diffusion_R0_B,              1.0,           0.0,              NoMax_double      );
   LOAD_PARA( load_mode, "CR_Diffusion_GX",           &CR_Diffusion_GX,                0,             0,                1                 );
   LOAD_PARA( load_mode, "CR_Diffusion_GY",           &CR_Diffusion_GY,                0,             0,                1                 );
   LOAD_PARA( load_mode, "CR_Diffusion_GZ",           &CR_Diffusion_GZ,                0,             0,                1                 );

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

// get the number of OpenMP threads
   int NT;
#  ifdef OPENMP
#  pragma omp parallel
#  pragma omp master
   {  NT = omp_get_num_threads();  }
#  else
   {  NT = 1;                      }
#  endif

// allocate RNG
   RNG = new RandomNumber_t( NT );

// set the random seed of each MPI rank
   for (int t=0; t<NT; t++) {
      RNG->SetSeed(t, CR_Diffusion_Seed+MPI_Rank*1000+t);
   }

// (1-2) set the default values
   if ( CR_Diffusion_CenterX < 0.0 )  CR_Diffusion_CenterX = 0.5*amr->BoxSize[0];
   if ( CR_Diffusion_CenterY < 0.0 )  CR_Diffusion_CenterY = 0.5*amr->BoxSize[1];
   if ( CR_Diffusion_CenterZ < 0.0 )  CR_Diffusion_CenterZ = 0.5*amr->BoxSize[2];

// (1-3) check the runtime parameters
   if ( CR_Diffusion_Type == 0  ||  CR_Diffusion_Type == 1  ||  CR_Diffusion_Type == 2  ||  CR_Diffusion_Type == 3 )
      Aux_Error( ERROR_INFO, "CR_Diffusion_Type (%d) is not supported yet !! (Unless you have fixed the all the fluid except cosmic ray.)\n", CR_Diffusion_Type );

   if ( (CR_Diffusion_Type == 0  ||  CR_Diffusion_Type == 3)  &&  CR_Diffusion_Mag_Type != 0 )
      Aux_Error( ERROR_INFO, "We only support the uniform magnetic field ( CR_Diffusion_Mag_Type = 0 ) for this test problem type (%d).\n", CR_Diffusion_Type );

   if ( CR_Diffusion_Type == 3  &&  CR_DIFF_PERP != 0.0 )
      Aux_Error( ERROR_INFO, "CR_DIFF_PERP must be 0.0 for this test problem type (%d).\n", CR_Diffusion_Type );


   if ( CR_Diffusion_Type == 0  ||  CR_Diffusion_Type == 1  ||  CR_Diffusion_Type == 2  ||  CR_Diffusion_Type == 3 )
      Aux_Message( stderr, "WARNING : All the fluid except cosmic ray should be fixed to follow the analytical solution. The adiabatic work done term should be disable also.\n" );

   if ( CR_Diffusion_Type == 1 )
      Aux_Message( stderr, "WARNING : The analytical solution is only worked for the short time evolution." );

   if ( CR_Diffusion_Mag_Type == 1 )
      Aux_Message( stderr, "WARNING : CR_Diffusion_MagY and CR_Diffusion_MagZ is useless in this type of magnetic field.\n" );

   if ( CR_Diffusion_Mag_Type == 2 )
      Aux_Message( stderr, "WARNING : This type of magnetic field is not divergence free.\n" );

   if ( CR_Diffusion_Mag_Type == 3 )
   {
      Aux_Message( stderr, "WARNING : CR_Diffusion_MagY and CR_Diffusion_MagZ is useless in this type of magnetic field.\n" );
      Aux_Message( stderr, "WARNING : This type of magnetic field is not divergence free.\n" );
   }

   if ( CR_Diffusion_Mag_Type == 4 )
      Aux_Message( stderr, "WARNING : CR_Diffusion_MagY and CR_Diffusion_MagZ is useless in this type of magnetic field.\n" );


// (2) set the problem-specific derived parameters
   CR_Diffusion_Dim   = CR_Diffusion_GX + CR_Diffusion_GY + CR_Diffusion_GZ;
   CR_Diffusion_Space = CR_Diffusion_GX + 2 * CR_Diffusion_GY + 4 * CR_Diffusion_GZ;


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


// (4) make a note
   if ( MPI_Rank == 0 )
   {
      Aux_Message( stdout, "=============================================================================\n" );
      Aux_Message( stdout, "  test problem ID       = %d\n",     TESTPROB_ID           );
      Aux_Message( stdout, "  CR_Diffusion_CenterX  = %14.7e\n", CR_Diffusion_CenterX  );
      Aux_Message( stdout, "  CR_Diffusion_CenterY  = %14.7e\n", CR_Diffusion_CenterY  );
      Aux_Message( stdout, "  CR_Diffusion_CenterZ  = %14.7e\n", CR_Diffusion_CenterZ  );
      Aux_Message( stdout, "  CR_Diffusion_Vx       = %14.7e\n", CR_Diffusion_Vx       );
      Aux_Message( stdout, "  CR_Diffusion_Vy       = %14.7e\n", CR_Diffusion_Vy       );
      Aux_Message( stdout, "  CR_Diffusion_Vz       = %14.7e\n", CR_Diffusion_Vz       );
      Aux_Message( stdout, "  CR_Diffusion_Rho0     = %14.7e\n", CR_Diffusion_Rho0     );
      Aux_Message( stdout, "  CR_Diffusion_PGas0    = %14.7e\n", CR_Diffusion_PGas0    );
      Aux_Message( stdout, "  CR_Diffusion_E0_CR    = %14.7e\n", CR_Diffusion_E0_CR    );
      Aux_Message( stdout, "  CR_Diffusion_BG_CR    = %14.7e\n", CR_Diffusion_BG_CR    );
      Aux_Message( stdout, "  CR_Diffusion_R02_CR   = %14.7e\n", CR_Diffusion_R02_CR   );
      Aux_Message( stdout, "  CR_Diffusion_Type     = %d\n",     CR_Diffusion_Type     );
      Aux_Message( stdout, "  CR_Diffusion_Mag_Type = %d\n",     CR_Diffusion_Mag_Type );
      Aux_Message( stdout, "  CR_Diffusion_MagX     = %14.7e\n", CR_Diffusion_MagX     );
      Aux_Message( stdout, "  CR_Diffusion_MagY     = %14.7e\n", CR_Diffusion_MagY     );
      Aux_Message( stdout, "  CR_Diffusion_MagZ     = %14.7e\n", CR_Diffusion_MagZ     );
      Aux_Message( stdout, "  CR_Diffusion_Seed     = %d\n",     CR_Diffusion_Seed     );
      Aux_Message( stdout, "  CR_Diffusion_R_In     = %14.7e\n", CR_Diffusion_R_In     );
      Aux_Message( stdout, "  CR_Diffusion_R_Out    = %14.7e\n", CR_Diffusion_R_Out    );
      Aux_Message( stdout, "  CR_Diffusion_CenterR  = %14.7e\n", CR_Diffusion_CenterR  );
      Aux_Message( stdout, "  CR_Diffusion_delR     = %14.7e\n", CR_Diffusion_delR     );
      Aux_Message( stdout, "  CR_Diffusion_delPhi   = %14.7e\n", CR_Diffusion_delPhi   );
      Aux_Message( stdout, "  CR_Diffusion_R0_CR    = %14.7e\n", CR_Diffusion_R0_CR    );
      Aux_Message( stdout, "  CR_Diffusion_R0_B     = %14.7e\n", CR_Diffusion_R0_B     );
      Aux_Message( stdout, "  CR_Diffusion_GX       = %d\n",     CR_Diffusion_GX       );
      Aux_Message( stdout, "  CR_Diffusion_GY       = %d\n",     CR_Diffusion_GY       );
      Aux_Message( stdout, "  CR_Diffusion_GZ       = %d\n",     CR_Diffusion_GZ       );
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
   double f;

   Dens = CR_Diffusion_Rho0;
   MomX = Dens * CR_Diffusion_Vx;
   MomY = Dens * CR_Diffusion_Vy;
   MomZ = Dens * CR_Diffusion_Vz;
   Pres = CR_Diffusion_PGas0;

   CRay = CR_Diffusion_E0_CR;
   double magD1, magD2, magD3, D1, D2, D3, fD1, fD2, fD3;

   switch ( CR_Diffusion_Space )
   {
//    1D simulation
      case 1: magD1 = CR_Diffusion_MagX; magD2 = CR_Diffusion_MagY; magD3 = CR_Diffusion_MagZ;
              D1 = x - CR_Diffusion_CenterX - CR_Diffusion_Vx*Time;
              D2 = y - CR_Diffusion_CenterY - CR_Diffusion_Vy*Time;
              D3 = z - CR_Diffusion_CenterZ - CR_Diffusion_Vz*Time;
              break;
      case 2: magD1 = CR_Diffusion_MagY; magD2 = CR_Diffusion_MagZ; magD3 = CR_Diffusion_MagX;
              D1 = y - CR_Diffusion_CenterY - CR_Diffusion_Vy*Time;
              D2 = z - CR_Diffusion_CenterZ - CR_Diffusion_Vz*Time;
              D3 = x - CR_Diffusion_CenterX - CR_Diffusion_Vx*Time;
              break;
      case 4: magD1 = CR_Diffusion_MagZ; magD2 = CR_Diffusion_MagX; magD3 = CR_Diffusion_MagY;
              D1 = z - CR_Diffusion_CenterZ - CR_Diffusion_Vz*Time;
              D2 = x - CR_Diffusion_CenterX - CR_Diffusion_Vx*Time;
              D3 = y - CR_Diffusion_CenterY - CR_Diffusion_Vy*Time;
              break;

//    2D simulation
      case 3: magD1 = CR_Diffusion_MagX; magD2 = CR_Diffusion_MagY; magD3 = CR_Diffusion_MagZ;
              D1 = x - CR_Diffusion_CenterX - CR_Diffusion_Vx*Time;
              D2 = y - CR_Diffusion_CenterY - CR_Diffusion_Vy*Time;
              D3 = z - CR_Diffusion_CenterZ - CR_Diffusion_Vz*Time;
              break;
      case 5: magD1 = CR_Diffusion_MagY; magD2 = CR_Diffusion_MagZ; magD3 = CR_Diffusion_MagX;
              D1 = y - CR_Diffusion_CenterY - CR_Diffusion_Vy*Time;
              D2 = z - CR_Diffusion_CenterZ - CR_Diffusion_Vz*Time;
              D3 = x - CR_Diffusion_CenterX - CR_Diffusion_Vx*Time;
              break;
      case 6: magD1 = CR_Diffusion_MagZ; magD2 = CR_Diffusion_MagX; magD3 = CR_Diffusion_MagY;
              D1 = z - CR_Diffusion_CenterZ - CR_Diffusion_Vz*Time;
              D2 = x - CR_Diffusion_CenterX - CR_Diffusion_Vx*Time;
              D3 = y - CR_Diffusion_CenterY - CR_Diffusion_Vy*Time;
              break;

//    3D simulation
      case 7: magD1 = CR_Diffusion_MagX; magD2 = CR_Diffusion_MagY; magD3 = CR_Diffusion_MagZ;
              D1 = x - CR_Diffusion_CenterX - CR_Diffusion_Vx*Time;
              D2 = y - CR_Diffusion_CenterY - CR_Diffusion_Vy*Time;
              D3 = z - CR_Diffusion_CenterZ - CR_Diffusion_Vz*Time;
              break;
   } // switch ( CR_Diffusion_Space )

   if ( CR_Diffusion_Type == 0 )
   {
      if ( CR_Diffusion_Dim == 1 )
      {
//       Analytical solution:
//       f(k, t) = 1 + 4 * R_0**2 * k * t
//       E(x, y) = E_0 * EXP[ -R_0**2 * x**2 / f(k_x, t) ] / SQRT[ f(k_x, t) ]
//       Check the Magnetic field to define the time evolution factor
         if      ( magD1 != 0  &&   magD2 == 0  &&  magD3 == 0  )   fD1 = 1.0 + 4.0 * CR_Diffusion_R02_CR * CR_DIFF_PARA * Time;
         else if ( magD1 == 0  &&  (magD2 != 0  ||  magD3 != 0) )   fD1 = 1.0 + 4.0 * CR_Diffusion_R02_CR * CR_DIFF_PERP * Time;
         else     Aux_Error( ERROR_INFO, "This kind of magnetic field is not supported for 1D diffusion case.\n" );

         CRay *= EXP( -CR_Diffusion_R02_CR * D1 * D1 / fD1 ) / SQRT(fD1);
      }

      else if ( CR_Diffusion_Dim == 2 )
      {
//       Check the Magnetic field to define the diffusion coefficient
         if ( magD3 != 0  &&  (magD1 != 0 || magD2 != 0) )
            Aux_Error( ERROR_INFO, "This kind of magnetic field is not supported for 2D diffusion case.\n" );

         else if ( magD3 != 0 )
         {
            fD1 = 1.0 + 4.0 * CR_Diffusion_R02_CR * CR_DIFF_PERP * Time;
            fD2 = 1.0 + 4.0 * CR_Diffusion_R02_CR * CR_DIFF_PERP * Time;
         }

         else
         {
//          Analytical solution:
//          f(k, t) = 1 + 4 * R_0**2 * k * t
//          k_x = k_parallel, k_y = k_perpendicular
//          E(x, y) = E_0 * EXP[ -R_0**2 * x**2 / f(k_x, t) ] / SQRT[ f(k_x, t) ]
//                        * EXP[ -R_0**2 * y**2 / f(k_y, t) ] / SQRT[ f(k_y, t) ]
//
//          We bascially apply a rotation matrix to rotate the magnetic field to (1, 0, 0) direction, so the simulation will be matched the solution above.
//          C = A \cross B, where A and B are the unit vector.
//          sin = C, cos = A \dot B
//          R = [  cos | sin ]
//              [ -sin | cos ]
//
//          B = RA, B = (1, 0), A = (B_x, B_y)
            double Mag_len2 = magD1 * magD1 + magD2 * magD2;
            double A[2] = {          magD1, magD2 };
            double B[2] = { SQRT(Mag_len2),   0.0 };
            double sin  = -A[1] * B[0] / Mag_len2;
            double cos  =  A[0] * B[0] / Mag_len2;
            double R[2][2] = { {  cos, sin },
                               { -sin, cos } };

            double pos_new[2] = {0.0, 0.0};
            double pos_old[2] = { D1, D2 };
            for (int i=0; i<2; i++)   for (int j = 0; j < 2; j++ )   pos_new[i] += R[i][j] * pos_old[j]; // new coordinate move to origin
            D1 = SQR( pos_new[0] );
            fD1 = 1.0 + 4.0 * CR_Diffusion_R02_CR * CR_DIFF_PARA * Time;
            D2 = SQR( pos_new[1] );
            fD2 = 1.0 + 4.0 * CR_Diffusion_R02_CR * CR_DIFF_PERP * Time;
         }
         CRay *= EXP( -CR_Diffusion_R02_CR * D1 * D1 / fD1 ) / SQRT(fD1);
         CRay *= EXP( -CR_Diffusion_R02_CR * D2 * D2 / fD2 ) / SQRT(fD2);
      }

      else if ( CR_Diffusion_Dim == 3 )
      {
//       Analytical solution:
//       f(k, t) = 1 + 4 * R_0**2 * k * t
//       k_x = k_parallel, k_y = k_z = k_perpendicular
//       E(x, y, z) = E_0 * EXP[ -R_0**2 * x**2 / f(k_x, t) ] / SQRT[ f(k_x, t) ]
//                        * EXP[ -R_0**2 * y**2 / f(k_y, t) ] / SQRT[ f(k_y, t) ]
//                        * EXP[ -R_0**2 * z**2 / f(k_z, t) ] / SQRT[ f(k_z, t) ]
//
//       We bascially apply a rotation matrix to rotate the magnetic field to (1, 0, 0) direction, so the simulation will be matched the solution above.
//       C = A \cross B, where A and B are the unit vector.
//       sin = abs(C), cos = A \dot B
//       R = I + v + v**2 / (1+cos), **NOTICE**: this would be invaliad when cos = -1.
//           [    1 - (C_2**2 + C_3**2) / (1+cos) | -C_3 +         (C_1*C_2) / (1+cos) |  C_2 +         (C_1*C_3) / (1+cos) ]
//       R = [  C_3 +         (C_1*C_2) / (1+cos) | 1    - (C_3**2 + C_1**2) / (1+cos) | -C_1 +         (C_2*C_3) / (1+cos) ]
//           [ -C_2 +         (C_1*C_3) / (1+cos) |  C_1 +         (C_2*C_3) / (1+cos) |    1 - (C_1**2 + C_2**2) / (1+cos) ]
//
//       B = RA, B = (1, 0, 0), A = (B_x, B_y, B_z)

         double Mag_len2 = magD1*magD1 + magD2*magD2 + magD3*magD3;
         double A[3]     = {          magD1,                 magD2,                 magD3 };
         double B[3]     = { SQRT(Mag_len2),                   0.0,                   0.0 };
         double C[3]     = {            0.0, A[2]*B[0] / Mag_len2 , -A[1]*B[0] / Mag_len2 };
         double cosP1    = A[0] * B[0] / Mag_len2 + 1.0;
         double R[3][3]  = { {   1.0 - (C[1]*C[1] + C[2]*C[2]) / cosP1, -C[2] +             (C[0]*C[1]) / cosP1,  C[1] +             (C[0]*C[2]) / cosP1 },
                             {  C[2] +             (C[0]*C[1]) / cosP1,   1.0 - (C[2]*C[2] + C[0]*C[0]) / cosP1, -C[0] +             (C[1]*C[2]) / cosP1 },
                             { -C[1] +             (C[0]*C[2]) / cosP1,  C[0] +             (C[1]*C[2]) / cosP1,   1.0 - (C[0]*C[0] + C[1]*C[1]) / cosP1 } };

         double pos_new[3] = {0.0, 0.0, 0.0};
         double pos_old[3] = { D1,  D2,  D3};
         for (int i=0; i<3; i++)   for (int j=0; j<3; j++)   pos_new[i] += R[i][j] * pos_old[j]; // new coordinate move to origin

         double x2 = SQR( pos_new[0] );
         double fx = 1.0 + 4.0 * CR_Diffusion_R02_CR * CR_DIFF_PARA * Time;
         CRay *= EXP( -CR_Diffusion_R02_CR * x2 / fx ) / SQRT(fx);
         double y2 = SQR( pos_new[1] );
         double fy = 1.0 + 4.0 * CR_Diffusion_R02_CR * CR_DIFF_PERP * Time;
         CRay *= EXP( -CR_Diffusion_R02_CR * y2 / fy ) / SQRT(fy);
         double z2 = SQR( pos_new[2] );
         double fz = 1.0 + 4.0 * CR_Diffusion_R02_CR * CR_DIFF_PERP * Time;
         CRay *= EXP( -CR_Diffusion_R02_CR * z2 / fz ) / SQRT(fz);
      }

      else
         Aux_Error( ERROR_INFO, "At least one of the {CR_Diffusion_GX, CR_Diffusion_GY, CR_Diffusion_GZ} should not be zero.\n" );
   }

   else if ( CR_Diffusion_Type == 1 )
   {
      if ( CR_Diffusion_Dim != 2 )   Aux_Error( ERROR_INFO, "This is 2-D test type. Only one of the CR_Diffusion_G{X,Y,Z} can be zero.\n" );

      const double r   = SQRT( SQR(D1) + SQR(D2) );
      const double phi = ATAN2( D2, D1 );
      const double D   = SQRT( 4. * Time * CR_DIFF_PARA );

      if ( Time == 0.0 )
      {
         CRay = (r > CR_Diffusion_R_In  &&  r < CR_Diffusion_R_Out && FABS(phi) < M_PI/12.) ? 2.0 : 0.0;
      }

      else
      {
         if ( r > CR_Diffusion_R_In  &&  r < CR_Diffusion_R_Out && FABS(phi) < M_PI/12. ) {
            CRay = erfc( (phi - M_PI/12.) * r / D ) + erfc( (phi + M_PI/12.) * r / D );
         } else {
            CRay = 0.0;
         }
      }
   }

   else if ( CR_Diffusion_Type == 2 )
   {
      if ( CR_Diffusion_Dim != 2 )   Aux_Error( ERROR_INFO, "This is 2-D test type. Only one of the CR_Diffusion_G{X,Y,Z} can be zero.\n" );
      const double r        = SQRT( SQR(D1) + SQR(D2) );
      const double phi      = ATAN2( D2, D1 );
      const double del_r2   = SQR(CR_Diffusion_delR);
      const double del_phi2 = SQR(CR_Diffusion_delPhi) + 2. * CR_DIFF_PARA * Time / SQR(r);

      CRay *= SQRT( SQR(CR_Diffusion_delPhi) / del_phi2 );
      CRay *= EXP( -0.5 * SQR(r-CR_Diffusion_CenterR) / del_r2 );
      CRay *= EXP( -0.5 * SQR(phi) / del_phi2 );
   }
   else if ( CR_Diffusion_Type == 3 )
   {
      if ( CR_Diffusion_Space == 0 )   Aux_Error( ERROR_INFO, "At least one of the CR_Diffusion_G{X,Y,Z} should not be zero.\n" );
      double r;
      switch ( CR_Diffusion_Dim )
      {
         case 1: r = D1;
                 if ( magD1 == 0.0  ||  magD2 != 0.0  ||  magD3 != 0.0 )
                    Aux_Error( ERROR_INFO, "The magnetic field must align to the simulation axis (%d,%d,%d).\n",
                               CR_Diffusion_GX, CR_Diffusion_GY, CR_Diffusion_GZ );
                 break;
         case 2: r = D1 + D2;
                 if ( magD1 == 0.0  ||  magD1 != magD2  ||  magD3 != 0.0 )
                    Aux_Error( ERROR_INFO, "The magnetic field must align to the diagonal direction (%d,%d,%d).\n",
                               CR_Diffusion_GX, CR_Diffusion_GY, CR_Diffusion_GZ );
                 break;
         case 3: r = D1 + D2 + D3;
                 if ( magD1 == 0.0  ||  magD1 != magD2  ||  magD1 != magD3 )
                    Aux_Error( ERROR_INFO, "The magnetic field must align to the diagonal direction (%d,%d,%d).\n",
                               CR_Diffusion_GX, CR_Diffusion_GY, CR_Diffusion_GZ );
                 break;
      } // switch ( CR_Diffusion_Dim )

      const double fr = 1.0 + 4.0 * CR_Diffusion_R02_CR * CR_DIFF_PARA * Time;
      CRay *= EXP( -CR_Diffusion_R02_CR * r * r / fr ) / SQRT(fr);
   }
   else if ( CR_Diffusion_Type == 4 )
   {
      const double r = SQRT( SQR(D1) + SQR(D2) + SQR(D3) );
      if ( r > CR_Diffusion_R0_CR )   CRay = 0.0;
   }
   else
   {
      Aux_Error( ERROR_INFO, "CR_Diffusion_Type = %d is NOT supported [0/1/2/3/4] !!\n", CR_Diffusion_Type );
   } // if ( CR_Diffusion_Type )

   CRay = CRay + CR_Diffusion_BG_CR;

   P_cr = EoS_CREint2CRPres_CPUPtr( CRay, EoS_AuxArray_Flt, EoS_AuxArray_Int, h_EoS_Table );
   Pres = Pres + P_cr;

// set the output array of passive scaler
   fluid[CRAY] = CRay;

   Eint = EoS_DensPres2Eint_CPUPtr( Dens, Pres, fluid+NCOMP_FLUID, EoS_AuxArray_Flt, EoS_AuxArray_Int, h_EoS_Table );
   Etot = Hydro_ConEint2Etot( Dens, MomX, MomY, MomZ, Eint, 0.0 );      // do NOT include magnetic energy here

// set the output array
   fluid[DENS] = Dens;
   fluid[MOMX] = MomX;
   fluid[MOMY] = MomY;
   fluid[MOMZ] = MomZ;
   fluid[ENGY] = Etot;

} // FUNCTION : SetGridIC



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

   double D1, D2, D3, _one = 1.0;
   if ( CR_Diffusion_Mag_Type == 1 ) _one = -1.0;

   const double pos[3] = {x - CR_Diffusion_CenterX, y - CR_Diffusion_CenterY, z - CR_Diffusion_CenterZ};

   switch ( CR_Diffusion_Space )
   {
      case 3: D1 = pos[0];        D2 = pos[1] * _one; D3 = 0.0;           break;
      case 5: D1 = pos[0] * _one; D2 = 0.0;           D3 = pos[2];        break;
      case 6: D1 = 0.0;           D2 = pos[1];        D3 = pos[2] * _one; break;
   } // switch ( CR_Diffusion_Space )

   if      ( CR_Diffusion_Mag_Type == 0 )
   {
      magnetic[MAGX] = CR_Diffusion_MagX;
      magnetic[MAGY] = CR_Diffusion_MagY;
      magnetic[MAGZ] = CR_Diffusion_MagZ;
   }

   else if ( CR_Diffusion_Mag_Type == 1 )
   {
      const double r = SQRT( SQR(D1) + SQR(D2) + SQR(D3) );
      magnetic[MAGX] = (D1 != 0.0) ? CR_Diffusion_MagX * (D2+D3)/r : 0.0;
      magnetic[MAGY] = (D2 != 0.0) ? CR_Diffusion_MagX * (D1+D3)/r : 0.0;
      magnetic[MAGZ] = (D3 != 0.0) ? CR_Diffusion_MagX * (D1+D2)/r : 0.0;
   }

   else if ( CR_Diffusion_Mag_Type == 2 )
   {
      magnetic[MAGX] = RandomNumber( RNG, -1.0, +1.0 );
      magnetic[MAGY] = RandomNumber( RNG, -1.0, +1.0 );
      magnetic[MAGZ] = RandomNumber( RNG, -1.0, +1.0 );
   }

   else if ( CR_Diffusion_Mag_Type == 3 )
   {
      const double r = SQRT( SQR(D1) + SQR(D2) + SQR(D3) );
      magnetic[MAGX] = CR_Diffusion_MagX * D1 / r;
      magnetic[MAGY] = CR_Diffusion_MagX * D2 / r;
      magnetic[MAGZ] = CR_Diffusion_MagX * D3 / r;
   }

   else if ( CR_Diffusion_Mag_Type == 4 )
   {
//    Reference: Suwa+ 2007, PASJ, 59, 771
//    A_phi = 0.5 * B0 * ( R0^3 / (r^3 + R0^3) ) * r * sin(theta), A_r = A_theta = 0
      const double r = SQRT( SQR(pos[0]) + SQR(pos[1]) + SQR(pos[2]) );
      const double _r0_cub = 1. / CUBE(CR_Diffusion_R0_B);
      const double _factor = 1. / SQR( 1. + CUBE(r) * _r0_cub );
      magnetic[MAGX] =  1.5 * CR_Diffusion_MagX * pos[0] * pos[2] * r * _r0_cub * _factor;
      magnetic[MAGY] =  1.5 * CR_Diffusion_MagX * pos[1] * pos[2] * r * _r0_cub * _factor;
      magnetic[MAGZ] =  SQRT(_factor)
                       -1.5 * CR_Diffusion_MagX * (pos[0]*pos[0] + pos[1]*pos[1]) * r * _r0_cub * _factor;
   }

   else
      Aux_Error( ERROR_INFO, "CR_Diffusion_Mag_Type = %d is NOT supported [0/1/2/3/4] !!\n", CR_Diffusion_Mag_Type );

} // FUNCTION : SetBFieldIC



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

   const char Prefix[100] = "CR_Diffusion";
   OptOutputPart_t Part;

   switch ( CR_Diffusion_Space )
   {
      case 1: Part = OUTPUT_X;    break;
      case 2: Part = OUTPUT_Y;    break;
      case 3: Part = OUTPUT_XY;   break;
      case 4: Part = OUTPUT_Z;    break;
      case 5: Part = OUTPUT_XZ;   break;
      case 6: Part = OUTPUT_YZ;   break;
      case 7: Part = OUTPUT_DIAG; break;
   } // switch ( CR_Diffusion_Space )

   Output_L1Error( SetGridIC, SetBFieldIC, Prefix, Part, OUTPUT_PART_X, OUTPUT_PART_Y, OUTPUT_PART_Z );

} // FUNCTION : OutputError
#endif // #if ( MODEL == HYDRO  &&  defined CR_DIFFUSION )



//-------------------------------------------------------------------------------------------------------
// Function    :  RandomNumber
// Description :  Generate a single random number
//
// Parameter   :  RNG : Random number generator
//                Min : Lower limit of the random number
//                Max : Upper limit of the random number
//
// Return      :  Random number
//-------------------------------------------------------------------------------------------------------
double RandomNumber( RandomNumber_t *RNG, const double Min, const double Max )
{

// thread-private variables
#  ifdef OPENMP
   const int TID = omp_get_thread_num();
#  else
   const int TID = 0;
#  endif

   return RNG->GetValue( TID, Min, Max );

} // FUNCTION : RandomNumber



//-------------------------------------------------------------------------------------------------------
// Function    :  Init_TestProb_Hydro_CR_Diffusion
// Description :  Test problem initializer
//
// Note        :  None
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Init_TestProb_Hydro_CR_Diffusion()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );


// validate the compilation flags and runtime parameters
   Validate();


#  if ( MODEL == HYDRO  &&  defined CR_DIFFUSION )
// set the problem-specific runtime parameters
   SetParameter();


   Init_Function_User_Ptr        = SetGridIC;
#  ifdef MHD
   Init_Function_BField_User_Ptr = SetBFieldIC;
#  endif
   Output_User_Ptr               = OutputError;
#  ifdef SUPPORT_HDF5
   Output_HDF5_InputTest_Ptr     = LoadInputTestProb;
#  endif
#  endif // #if ( MODEL == HYDRO  &&  defined CR_DIFFUSION )


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : Init_TestProb_Hydro_CR_Diffusion
