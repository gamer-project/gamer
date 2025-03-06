#include "GAMER.h"



// problem-specific global variables
// =======================================================================================
static int    MHDLinear_Mode;          // (1/2/3) --> (fast/slow/Alfven wave)
static double MHDLinear_Rho0;          // background density
static double MHDLinear_Rho1;          // amplitude of the density perturbation
static double MHDLinear_P0;            // background pressure
static double MHDLinear_v0;            // background velocity
static double MHDLinear_B0;            // background magnetic field
static int    MHDLinear_Dir;           // wave direction: (0/1/2/3) --> (x/y/z/diagonal)
static double MHDLinear_Sign;          // (+1/-1) --> (right/left-moving wave)
static double MHDLinear_Phase0;        // initial phase shift

static double MHDLinear_WaveSpeed;     // propagation speed
static double MHDLinear_WaveLength;    // wavelength
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

#  ifdef GRAVITY
   Aux_Error( ERROR_INFO, "GRAVITY must be disabled !!\n" );
#  endif

#  ifdef PARTICLE
   Aux_Error( ERROR_INFO, "PARTICLE must be disabled !!\n" );
#  endif

#  if ( EOS != EOS_GAMMA )
   Aux_Error( ERROR_INFO, "EOS != EOS_GAMMA !!\n" );
#  endif

   if ( MHDLinear_Dir == 3  &&  ( amr->BoxSize[0] != amr->BoxSize[1] || amr->BoxSize[0] != amr->BoxSize[2] )  )
      Aux_Error( ERROR_INFO, "simulation domain must be cubic for MHDLinear_Dir = %d !!\n", MHDLinear_Dir );

   for (int f=0; f<6; f++)
   if ( OPT__BC_FLU[f] != BC_FLU_PERIODIC )
      Aux_Error( ERROR_INFO, "must adopt periodic BC (i.e., \"OPT__BC_FLU_* = 1\") !!\n" );


// warnings
   if ( MPI_Rank == 0 )
   {
#     ifndef FLOAT8
      Aux_Message( stderr, "WARNING : it's recommended to enable FLOAT8 for this test !!\n" );
#     endif

      if ( !OPT__OUTPUT_USER )   Aux_Message( stdout, "WARNING : OPT__OUTPUT_USER is off !!\n" );
   }


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Validating test problem %d ... done\n", TESTPROB_ID );

} // FUNCTION : Validate



#if ( MODEL == HYDRO )
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
// ************************************************************************************************************************
// LOAD_PARA( load_mode, "KEY_IN_THE_FILE",   &VARIABLE,               DEFAULT,      MIN,              MAX               );
// ************************************************************************************************************************
   LOAD_PARA( load_mode, "MHDLinear_Mode",    &MHDLinear_Mode,        -1,            1,                1                 );
   LOAD_PARA( load_mode, "MHDLinear_Rho0",    &MHDLinear_Rho0,        -1.0,          Eps_double,       NoMax_double      );
   LOAD_PARA( load_mode, "MHDLinear_Rho1",    &MHDLinear_Rho1,        -1.0,          0.0,              NoMax_double      );
   LOAD_PARA( load_mode, "MHDLinear_P0",      &MHDLinear_P0,          -1.0,          Eps_double,       NoMax_double      );
   LOAD_PARA( load_mode, "MHDLinear_v0",      &MHDLinear_v0,           0.0,          NoMin_double,     NoMax_double      );
   LOAD_PARA( load_mode, "MHDLinear_Dir",     &MHDLinear_Dir,          3,            0,                3                 );
   LOAD_PARA( load_mode, "MHDLinear_B0",      &MHDLinear_B0,          -1.0,          0.0,              NoMax_double      );
   LOAD_PARA( load_mode, "MHDLinear_Sign",    &MHDLinear_Sign,         1.0,          NoMin_double,     NoMax_double      );
   LOAD_PARA( load_mode, "MHDLinear_Phase0",  &MHDLinear_Phase0,       0.0,          NoMin_double,     NoMax_double      );

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


// (2) set the problem-specific derived parameters
   MHDLinear_WaveLength = ( MHDLinear_Dir == 3 ) ? amr->BoxSize[0]/sqrt(3.0) : amr->BoxSize[MHDLinear_Dir];

// assuming EOS_GAMMA
   if ( MHDLinear_Mode == 1 )
      MHDLinear_WaveSpeed = sqrt( GAMMA*MHDLinear_P0/MHDLinear_Rho0 + SQR(MHDLinear_B0)/MHDLinear_Rho0 );
   else
      Aux_Error( ERROR_INFO, "unsupported MHDLinear_Mode = %d !!\n", MHDLinear_Mode );

// force MHDLinear_Sign to be +1.0/-1.0
   if ( MHDLinear_Sign >= 0.0 )  MHDLinear_Sign = +1.0;
   else                          MHDLinear_Sign = -1.0;


// (3) reset other general-purpose parameters
//     --> a helper macro PRINT_RESET_PARA is defined in Macro.h
   const double End_T_Default    = MHDLinear_WaveLength / MHDLinear_WaveSpeed;
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
      Aux_Message( stdout, "  test problem ID      = %d\n",      TESTPROB_ID          );
      Aux_Message( stdout, "  mode                 = %d\n",      MHDLinear_Mode       );
      Aux_Message( stdout, "  background density   = % 14.7e\n", MHDLinear_Rho0       );
      Aux_Message( stdout, "  density perturbation = % 14.7e\n", MHDLinear_Rho1       );
      Aux_Message( stdout, "  background pressure  = % 14.7e\n", MHDLinear_P0         );
      Aux_Message( stdout, "  background velocity  = % 14.7e\n", MHDLinear_v0         );
      Aux_Message( stdout, "  background B field   = % 14.7e\n", MHDLinear_B0         );
      Aux_Message( stdout, "  direction            = %d\n",      MHDLinear_Dir        );
      Aux_Message( stdout, "  sign (R/L)           = % 14.7e\n", MHDLinear_Sign       );
      Aux_Message( stdout, "  initial phase shift  = % 14.7e\n", MHDLinear_Phase0     );
      Aux_Message( stdout, "  wave speed           = % 14.7e\n", MHDLinear_WaveSpeed  );
      Aux_Message( stdout, "  wavelength           = % 14.7e\n", MHDLinear_WaveLength );
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
   const double Rho0       = MHDLinear_Rho0;
   const double Rho1       = MHDLinear_Rho1;
   const double P0         = MHDLinear_P0;
   const double v0         = MHDLinear_v0;
   const double B0         = MHDLinear_B0;
   const double Sign       = MHDLinear_Sign;
   const double Phase0     = MHDLinear_Phase0;
   const double WaveSpeed  = MHDLinear_WaveSpeed;
   const double WaveLength = MHDLinear_WaveLength;

   double r, v1, B1, P1, WaveK, WaveW, SinPhase;
   double Dens, Mom, MomX, MomY, MomZ, Pres, Eint, Etot;

   switch ( MHDLinear_Dir ) {
      case 0:  r = x;                        break;
      case 1:  r = y;                        break;
      case 2:  r = z;                        break;
      case 3:  r = ( x + y + z )/sqrt(3.0);  break;
   }
   r -= v0*Time;

   if ( MHDLinear_Mode == 1 )
   {
      v1       = Sign*Rho1*WaveSpeed/Rho0;
      B1       = Sign*v1*B0/WaveSpeed;
      P1       = Sign*v1*P0*GAMMA/WaveSpeed;    // assuming EOS_GAMMA
      WaveK    = 2.0*M_PI/WaveLength;
      WaveW    = WaveK*WaveSpeed;
      SinPhase = sin( WaveK*r - Sign*WaveW*Time + Phase0 );

      Dens = Rho0 + Rho1*SinPhase;
      Mom  = Dens*( v1*SinPhase + v0 );

      switch ( MHDLinear_Dir ) {
         case 0:  MomX = Mom;            MomY = 0.0;   MomZ = 0.0;   break;
         case 1:  MomX = 0.0;            MomY = Mom;   MomZ = 0.0;   break;
         case 2:  MomX = 0.0;            MomY = 0.0;   MomZ = Mom;   break;
         case 3:  MomX = Mom/sqrt(3.0);  MomY = MomX;  MomZ = MomX;  break;
      }

      Pres = P0 + P1*SinPhase;
   } // if ( MHDLinear_Mode == 1 )

   else
      Aux_Error( ERROR_INFO, "unsupported MHDLinear_Mode = %d !!\n", MHDLinear_Mode );

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

// short names
   const double Rho0       = MHDLinear_Rho0;
   const double Rho1       = MHDLinear_Rho1;
   const double v0         = MHDLinear_v0;
   const double B0         = MHDLinear_B0;
   const double Sign       = MHDLinear_Sign;
   const double Phase0     = MHDLinear_Phase0;
   const double WaveSpeed  = MHDLinear_WaveSpeed;
   const double WaveLength = MHDLinear_WaveLength;

   double r, v1, B1, B, WaveK, WaveW, SinPhase;

   switch ( MHDLinear_Dir ) {
      case 0:  r = x;                        break;
      case 1:  r = y;                        break;
      case 2:  r = z;                        break;
      case 3:  r = ( x + y + z )/sqrt(3.0);  break;
   }
   r -= v0*Time;

   if ( MHDLinear_Mode == 1 )
   {
      v1       = Sign*Rho1*WaveSpeed/Rho0;
      B1       = Sign*v1*B0/WaveSpeed;
      WaveK    = 2.0*M_PI/WaveLength;
      WaveW    = WaveK*WaveSpeed;
      SinPhase = sin( WaveK*r - Sign*WaveW*Time + Phase0 );
      B        = B0 + B1*SinPhase;

      switch ( MHDLinear_Dir ) {
         case 0: magnetic[MAGX] = 0.0;
                 magnetic[MAGY] = +B / sqrt(2.0);
                 magnetic[MAGZ] = +B / sqrt(2.0);
                 break;

         case 1: magnetic[MAGX] = +B / sqrt(2.0);
                 magnetic[MAGY] = 0.0;
                 magnetic[MAGZ] = +B / sqrt(2.0);
                 break;

         case 2: magnetic[MAGX] = +B / sqrt(2.0);
                 magnetic[MAGY] = +B / sqrt(2.0);
                 magnetic[MAGZ] = 0.0;
                 break;

         case 3: magnetic[MAGX] = -B / ( 2.0*sqrt(1.5) );
                 magnetic[MAGY] = -B / ( 2.0*sqrt(1.5) );
                 magnetic[MAGZ] = +B / sqrt(1.5);
                 break;
      }
   } // if ( MHDLinear_Mode == 1 )

   else
      Aux_Error( ERROR_INFO, "unsupported MHDLinear_Mode = %d !!\n", MHDLinear_Mode );

} // FUNCTION : SetBFieldIC
#endif // #ifdef MHD



#ifdef MHD
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

   const char Prefix[100]     = "MHDLinearWave";
   const OptOutputPart_t Part = OUTPUT_X + MHDLinear_Dir;

   Output_L1Error( SetGridIC, SetBFieldIC, Prefix, Part, OUTPUT_PART_X, OUTPUT_PART_Y, OUTPUT_PART_Z );

} // FUNCTION : OutputError
#endif // #ifdef MHD
#endif // #if ( MODEL == HYDRO )



//-------------------------------------------------------------------------------------------------------
// Function    :  Init_TestProb_Hydro_MHD_LinearWave
// Description :  Test problem initializer
//
// Note        :  None
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Init_TestProb_Hydro_MHD_LinearWave()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );


// validate the compilation flags and runtime parameters
   Validate();


#  if ( MODEL == HYDRO )
// set the problem-specific runtime parameters
   SetParameter();


// set the function pointers of various problem-specific routines
   Init_Function_User_Ptr        = SetGridIC;
#  ifdef MHD
   Init_Function_BField_User_Ptr = SetBFieldIC;
   Output_User_Ptr               = OutputError;
#  endif
#  ifdef SUPPORT_HDF5
   Output_HDF5_InputTest_Ptr     = LoadInputTestProb;
#  endif
#  endif // #if ( MODEL == HYDRO )


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : Init_TestProb_Hydro_MHD_LinearWave
