#include "GAMER.h"



// problem-specific global variables
// =======================================================================================
static double CR_Acoustic_Delta;        // amplitude of the velocity perturbation
static double CR_Acoustic_Rho0;         // background density
static double CR_Acoustic_Pres0;        // background pressure
static double CR_Acoustic_Pres_CR0;     // background pressure of cosmic rays
static double CR_Acoustic_V0;           // background velocity
static double CR_Acoustic_Sign;         // (+1/-1) --> (right/left-moving wave)
static double CR_Acoustic_Phase;        // initial phase shift
static int    CR_Acoustic_Dir;          // wave direction (0/1/2/3) --> (x/y/z/diagonal)

static double CR_Acoustic_WaveSpeed;    // wave speed
static double CR_Acoustic_WaveL;        // wavelength
// =======================================================================================


// problem-specific function prototypes
static void OutputError();




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

#  if ( defined COSMIC_ARY  &&  EOS != EOS_COSMIC_RAY )
   Aux_Error( ERROR_INFO, "EOS != EOS_COSMIC_RAY when enable COSMIC_RAY!!\n" );
#  endif

   if ( CR_Acoustic_Dir == 3  &&  ( amr->BoxSize[0] != amr->BoxSize[1] || amr->BoxSize[0] != amr->BoxSize[2] )  )
      Aux_Error( ERROR_INFO, "simulation domain must be cubic for CR_Acoustic_Dir = %d !!\n", CR_Acoustic_Dir );

   for (int f=0; f<6; f++)
   if ( OPT__BC_FLU[f] != BC_FLU_PERIODIC )
      Aux_Error( ERROR_INFO, "please set \"OPT__BC_FLU_* = 1\" (i.e., periodic BC) !!\n" );


// warnings
   if ( MPI_Rank == 0 )
   {
      if ( !OPT__OUTPUT_USER )   Aux_Message( stdout, "WARNING : OPT__OUTPUT_USER is off !!\n" );
   }


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
// *******************************************************************************************************************************
// LOAD_PARA( load_mode, "KEY_IN_THE_FILE",      &VARIABLE,                  DEFAULT,       MIN,              MAX               );
// *******************************************************************************************************************************
   LOAD_PARA( load_mode, "CR_Acoustic_Delta",    &CR_Acoustic_Delta,         0.0,           0.0,              NoMax_double      );
   LOAD_PARA( load_mode, "CR_Acoustic_Rho0",     &CR_Acoustic_Rho0,          0.0,           0.0,              NoMax_double      );
   LOAD_PARA( load_mode, "CR_Acoustic_Pres0",    &CR_Acoustic_Pres0,         0.0,           0.0,              NoMax_double      );
   LOAD_PARA( load_mode, "CR_Acoustic_Pres_CR0", &CR_Acoustic_Pres_CR0,      0.0,           0.0,              NoMax_double      );
   LOAD_PARA( load_mode, "CR_Acoustic_V0",       &CR_Acoustic_V0,            0.0,           NoMin_double,     NoMax_double      );
   LOAD_PARA( load_mode, "CR_Acoustic_Sign",     &CR_Acoustic_Sign,          0.0,           NoMin_double,     NoMax_double      );
   LOAD_PARA( load_mode, "CR_Acoustic_Phase",    &CR_Acoustic_Phase,         0.0,           NoMin_double,     NoMax_double      );
   LOAD_PARA( load_mode, "CR_Acoustic_Dir",      &CR_Acoustic_Dir,           0,             0,                3                 );

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
// force Acoustic_Sign to be +1.0/-1.0
   if ( CR_Acoustic_Sign >= 0.0 )   CR_Acoustic_Sign = +1.0;
   else                             CR_Acoustic_Sign = -1.0;

// (1-3) check the runtime parameters


// (2) set the problem-specific derived parameters
   CR_Acoustic_WaveSpeed = SQRT(  ( GAMMA*CR_Acoustic_Pres0 + GAMMA_CR*CR_Acoustic_Pres_CR0 ) / CR_Acoustic_Rho0  );
   CR_Acoustic_WaveL     = ( CR_Acoustic_Dir == 3 ) ? amr->BoxSize[0]/sqrt(3.0) : amr->BoxSize[CR_Acoustic_Dir];

// (3) reset other general-purpose parameters
//     --> a helper macro PRINT_RESET_PARA is defined in Macro.h
   const long   End_Step_Default = __INT_MAX__;
   const double End_T_Default    = CR_Acoustic_WaveL / CR_Acoustic_WaveSpeed;

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
      Aux_Message( stdout, "  CR_Acoustic_Delta     = %14.7e\n", CR_Acoustic_Delta     );
      Aux_Message( stdout, "  CR_Acoustic_Rho0      = %14.7e\n", CR_Acoustic_Rho0      );
      Aux_Message( stdout, "  CR_Acoustic_Pres0     = %14.7e\n", CR_Acoustic_Pres0     );
      Aux_Message( stdout, "  CR_Acoustic_Pres_CR0  = %14.7e\n", CR_Acoustic_Pres_CR0  );
      Aux_Message( stdout, "  CR_Acoustic_V0        = %14.7e\n", CR_Acoustic_V0        );
      Aux_Message( stdout, "  CR_Acoustic_Sign      = %14.7e\n", CR_Acoustic_Sign      );
      Aux_Message( stdout, "  CR_Acoustic_Phase     = %14.7e\n", CR_Acoustic_Phase     );
      Aux_Message( stdout, "  CR_Acoustic_Dir       = %d\n",     CR_Acoustic_Dir       );
      Aux_Message( stdout, "  CR_Acoustic_WaveSpeed = %14.7e\n", CR_Acoustic_WaveSpeed );
      Aux_Message( stdout, "  CR_Acoustic_WaveL     = %14.7e\n", CR_Acoustic_WaveL     );
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

   const double cs         = CR_Acoustic_WaveSpeed;
   const double delta_cs   = CR_Acoustic_Delta / cs;
   const double wavelength = CR_Acoustic_WaveL;

   double Dens, MomX, MomY, MomZ, Pres, Eint, Etot, P_cr, CRay;
   double WaveK, WaveW, r, wave, Mom;


   WaveK = 2.0*M_PI/wavelength;
   WaveW = WaveK * cs;

   switch ( CR_Acoustic_Dir ) {
      case 0:  r = x;                        break;
      case 1:  r = y;                        break;
      case 2:  r = z;                        break;
      case 3:  r = ( x + y + z )/sqrt(3.0);  break;
   }

   r    -= CR_Acoustic_V0 * Time;
   wave  = sin( WaveK*r - CR_Acoustic_Sign*WaveW*Time + CR_Acoustic_Phase );


   Dens = ( 1.0 + delta_cs*wave )*CR_Acoustic_Rho0;
   Mom  = Dens*( CR_Acoustic_Sign*CR_Acoustic_Delta*wave + CR_Acoustic_V0 );

   switch ( CR_Acoustic_Dir ) {
      case 0:  MomX = Mom;            MomY = 0.0;   MomZ = 0.0;   break;
      case 1:  MomX = 0.0;            MomY = Mom;   MomZ = 0.0;   break;
      case 2:  MomX = 0.0;            MomY = 0.0;   MomZ = Mom;   break;
      case 3:  MomX = Mom/sqrt(3.0);  MomY = MomX;  MomZ = MomX;  break;
   }

   Pres = ( 1.0 + delta_cs*wave*GAMMA )*CR_Acoustic_Pres0;

   double GAMMA_CR_m1_inv = 1.0 / (GAMMA_CR - 1.0);
   P_cr = ( 1.0 + delta_cs*wave*GAMMA_CR )*CR_Acoustic_Pres_CR0;
   Pres = Pres + P_cr;
   CRay = GAMMA_CR_m1_inv*P_cr;

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

   const char Prefix[100]     = "CR_SoundWave";
   const OptOutputPart_t Part = OUTPUT_X + CR_Acoustic_Dir;

#  ifdef MHD
   Output_L1Error( SetGridIC, SetBFieldIC, Prefix, Part, OUTPUT_PART_X, OUTPUT_PART_Y, OUTPUT_PART_Z );
#  else
   Output_L1Error( SetGridIC, NULL,        Prefix, Part, OUTPUT_PART_X, OUTPUT_PART_Y, OUTPUT_PART_Z );
#  endif

} // FUNCTION : OutputError
#endif // #if ( MODEL == HYDRO  &&  defined COSMIC_RAY )



//-------------------------------------------------------------------------------------------------------
// Function    :  Init_TestProb_Hydro_CR_SoundWave
// Description :  Test problem initializer
//
// Note        :  None
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Init_TestProb_Hydro_CR_SoundWave()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );


// validate the compilation flags and runtime parameters
   Validate();


#  if ( MODEL == HYDRO  &&  defined COSMIC_RAY )
// set the problem-specific runtime parameters
   SetParameter();


// set the function pointers of various problem-specific routines
   Init_Function_User_Ptr        = SetGridIC;
#  ifdef MHD
   Init_Function_BField_User_Ptr = SetBFieldIC;
#  endif
   Output_User_Ptr               = OutputError;
#  ifdef SUPPORT_HDF5
   Output_HDF5_InputTest_Ptr     = LoadInputTestProb;
#  endif
#  endif // #if ( MODEL == HYDRO  &&  defined COSMIC_RAY )


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : Init_TestProb_Hydro_CR_SoundWave
