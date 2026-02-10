#include "GAMER.h"


// problem-specific global variables
// =======================================================================================
static double Ring_U;                // background internal energy
static double Ring_R1;               // inner radius of ring
static double Ring_R2;               // outer radius of ring
static double Ring_Rho;              // background mass density
static double Ring_Angle;            // angle defining initial region of U2
static double chi;                   // factor for conduction coefficient
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

#  ifndef MHD
   Aux_Error( ERROR_INFO, "MHD must be enabled !!\n" );
#  endif

#  ifndef CONDUCTION
   Aux_Error( ERROR_INFO, "CONDUCTION must be enabled !!\n" );
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

   if ( !OPT__FREEZE_HYDRO )
      Aux_Error( ERROR_INFO, "OPT__FREEZE_HYDRO must be enabled !!\n" );

#  ifdef CONDUCTION
   if ( CONDUCTION_FLUX_TYPE != ANISOTROPIC_CONDUCTION )
      Aux_Error( ERROR_INFO, "please set \"CONDUCTION_FLUX_TYPE = 2\" (i.e., anisotropic conduction) !!\n" );
#  endif

   for (int f=0; f<6; f++)
   if ( OPT__BC_FLU[f] != BC_FLU_OUTFLOW )
      Aux_Error( ERROR_INFO, "please set \"OPT__BC_FLU_* = 2\" (i.e., outflow BC) !!\n" );

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
// *****************************************************************************************************************
// LOAD_PARA( load_mode, "KEY_IN_THE_FILE",  &VARIABLE,        DEFAULT,       MIN,              MAX               );
// *****************************************************************************************************************
   LOAD_PARA( load_mode, "Ring_U",           &Ring_U,          10.0,          Eps_double,       NoMax_double      );
   LOAD_PARA( load_mode, "Ring_Rho",         &Ring_Rho,        1.0,           Eps_double,       NoMax_double      );
   LOAD_PARA( load_mode, "Ring_R1",          &Ring_R1,         0.5,           0.0,              NoMax_double      );
   LOAD_PARA( load_mode, "Ring_R2",          &Ring_R2,         0.7,           0.0,              NoMax_double      );
   LOAD_PARA( load_mode, "Ring_Angle",       &Ring_Angle,      15.0,          0.0,              360.0             );

} // FUNCTION : LoadInputTestProb

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

   LoadInputTestProb( LOAD_READPARA, ReadPara, NULL );

   ReadPara->Read( FileName );

   delete ReadPara;

// (2) set the problem-specific derived parameters

#  ifdef CONDUCTION
   chi = CONDUCTION_CONSTANT_COEFF*(GAMMA-1.0)*MOLECULAR_WEIGHT/Ring_Rho;
#  endif

// (3) reset other general-purpose parameters
//     --> a helper macro PRINT_RESET_PARA is defined in Macro.h
   const double End_T_Default    = 200.0;
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
      Aux_Message( stdout, "  Ring_U1              = % 14.7e\n", Ring_U  );
      Aux_Message( stdout, "  Ring_R1              = % 14.7e\n", Ring_R1 );
      Aux_Message( stdout, "  Ring_R2              = % 14.7e\n", Ring_R2 );
      Aux_Message( stdout, "  Ring_Rho             = % 14.7e\n", Ring_Rho );
      Aux_Message( stdout, "  Ring_Angle           = % 14.7e\n", Ring_Angle );
      Aux_Message( stdout, "=============================================================================\n" );
   }

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Setting runtime parameters ... done\n" );

} // FUNCTION : SetParameter



//-------------------------------------------------------------------------------------------------------
// Function    :  SetGridIC
// Description :  Set the problem-specific initial condition on grids
//
// Note        :  This function may also be used to estimate the numerical errors when OPT__OUTPUT_USER is enabled
//                --> In this case, it should provide the analytical solution at the given "Time"
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

   const double x0         = 0.5*amr->BoxSize[0];
   const double y0         = 0.5*amr->BoxSize[0];
   const double angle_rads = Ring_Angle*M_PI/180.0;
   const double r          = SQRT( SQR(x-x0) + SQR(y-y0) );
   const double phi        = SATAN2( y-y0, x-x0 );
   const double D          = 2.0*SQRT(chi*Time);

   double Dens, MomX, MomY, MomZ, Eint, Etot;

// inside region
   if ( 0.5 < r && r < 0.7 )
      Eint = Ring_U*Ring_Rho + erfc((phi-angle_rads)*r/D)-erfc((phi+angle_rads)*r/D);
// outside region
   else
      Eint = Ring_U*Ring_Rho;

   Dens = Ring_Rho;
   MomX = 0.0;
   MomY = 0.0;
   MomZ = 0.0;

// compute the total gas energy
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

   const double x0         = 0.5*amr->BoxSize[0];
   const double y0         = 0.5*amr->BoxSize[0];
   const double r          = SQRT( SQR(x-x0) + SQR(y-y0) );
   const double phi        = SATAN2( y-y0, x-x0 );

   magnetic[MAGX] = - ( y - y0 ) * 0.5/SQRT(M_PI) / r;
   magnetic[MAGY] =   ( x - x0 ) * 0.5/SQRT(M_PI) / r;
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
static void OutputError()
{

   const char Prefix[100]     = "ConductionRing";
   const OptOutputPart_t Part = OUTPUT_XY;

#  ifdef MHD
   Output_L1Error( SetGridIC, SetBFieldIC, Prefix, Part, OUTPUT_PART_X, OUTPUT_PART_Y, 0.5*amr->BoxSize[2] );
#  else
   Output_L1Error( SetGridIC, NULL,        Prefix, Part, OUTPUT_PART_X, OUTPUT_PART_Y, 0.5*amr->BoxSize[2] );
#  endif

} // FUNCTION : OutputError

#endif // #if ( MODEL == HYDRO )


//-------------------------------------------------------------------------------------------------------
// Function    :  Init_TestProb_Hydro_ConductionRing
// Description :  Test problem initializer
//
// Note        :  None
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Init_TestProb_Hydro_ConductionRing()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );


// validate the compilation flags and runtime parameters
   Validate();


#  if ( MODEL == HYDRO )
// set the problem-specific runtime parameters
   SetParameter();


// set the function pointers of various problem-specific routines
   Init_Function_User_Ptr = SetGridIC;
   Output_User_Ptr               = OutputError;
#  ifdef MHD
   Init_Function_BField_User_Ptr = SetBFieldIC;
#  endif
#  ifdef SUPPORT_HDF5
   Output_HDF5_InputTest_Ptr     = LoadInputTestProb;
#  endif
#  endif // #if ( MODEL == HYDRO )


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : Init_TestProb_Hydro_ConductionRing
