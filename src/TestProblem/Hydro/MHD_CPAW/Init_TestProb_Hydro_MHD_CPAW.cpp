#include "GAMER.h"



// problem-specific global variables
// =======================================================================================
static double MHD_CPAW_Rho0;          // background density
static double MHD_CPAW_A;             // perturbation amplitude
static double MHD_CPAW_P0;            // background pressure
static double MHD_CPAW_B0;            // background magnetic field
static double MHD_CPAW_Angle2;        // 
static double MHD_CPAW_Angle3;        // 

static double MHD_CPAW_WaveSpeed;     // propagation speed
static double MHD_CPAW_WaveLength;    // wavelength
static double MHD_CPAW_WaveNumber;    // wavenumber
static double MHD_CPAW_WaveFrequency; // wave frequency
static double MHD_CPAW_CosAngle2;     // cosine of angle 2
static double MHD_CPAW_SinAngle2;     // sine of angle 2
static double MHD_CPAW_CosAngle3;     // cosine of angle 3
static double MHD_CPAW_SinAngle3;     // sine of angle 3
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

#  ifdef MHD
   if ( OPT__INIT_BFIELD_BYVECPOT != INIT_MAG_BYVECPOT_FUNC )
      Aux_Error( ERROR_INFO, "Must set OPT__INIT_BFIELD_BYVECPOT != %d !!\n", INIT_MAG_BYVECPOT_FUNC );
#  endif
   
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
// **********************************************************************************************************************
// LOAD_PARA( load_mode, "KEY_IN_THE_FILE",  &VARIABLE,              DEFAULT,      MIN,              MAX               );
// **********************************************************************************************************************
   LOAD_PARA( load_mode, "MHD_CPAW_Rho0",    &MHD_CPAW_Rho0,         1.0,          Eps_double,       NoMax_double      );
   LOAD_PARA( load_mode, "MHD_CPAW_A",       &MHD_CPAW_A,            1.0e-3,       Eps_double,       NoMax_double      );
   LOAD_PARA( load_mode, "MHD_CPAW_P0",      &MHD_CPAW_P0,           0.1,          Eps_double,       NoMax_double      );
   LOAD_PARA( load_mode, "MHD_CPAW_B0",      &MHD_CPAW_B0,           1.0,          Eps_double,       NoMax_double      );
   LOAD_PARA( load_mode, "MHD_CPAW_Angle2",  &MHD_CPAW_Angle2,      -1.0,          NoMin_double,     180.0             );
   LOAD_PARA( load_mode, "MHD_CPAW_Angle3",  &MHD_CPAW_Angle3,      -1.0,          NoMin_double,     180.0             );

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
   // (1-1) read parameters from Input__TestProb
   const char FileName[] = "Input__TestProb";
   ReadPara_t *ReadPara  = new ReadPara_t;

   LoadInputTestProb( LOAD_READPARA, ReadPara, NULL );

   ReadPara->Read( FileName );

   delete ReadPara;

// (2) set the problem-specific derived parameters
   if ( MHD_CPAW_Angle3 == -1.0) 
   {
      MHD_CPAW_Angle3  = ATAN( amr->BoxSize[0]/amr->BoxSize[1] );
      MHD_CPAW_Angle3 *= 180.0/M_PI;
   }
   MHD_CPAW_CosAngle3     = cos( MHD_CPAW_Angle3*M_PI/180.0 );
   MHD_CPAW_SinAngle3     = sin( MHD_CPAW_Angle3*M_PI/180.0 );

   if ( MHD_CPAW_Angle2 == -1.0) 
   {
      MHD_CPAW_Angle2 = ATAN( 0.5*( amr->BoxSize[0] * MHD_CPAW_CosAngle3 + 
                                    amr->BoxSize[1] * MHD_CPAW_SinAngle3 ) /
                                    amr->BoxSize[2] );
      MHD_CPAW_Angle2 *= 180.0/M_PI;
   }
   MHD_CPAW_CosAngle2     = cos( MHD_CPAW_Angle2*M_PI/180.0 );
   MHD_CPAW_SinAngle2     = sin( MHD_CPAW_Angle2*M_PI/180.0 );

   const double x1 = amr->BoxSize[0]*MHD_CPAW_CosAngle2*MHD_CPAW_CosAngle3;
   const double x2 = amr->BoxSize[1]*MHD_CPAW_CosAngle2*MHD_CPAW_SinAngle3;
   const double x3 = amr->BoxSize[2]*MHD_CPAW_SinAngle2;

   MHD_CPAW_WaveLength    = x1;
   if ( MHD_CPAW_Angle3 != 0.0 ) MHD_CPAW_WaveLength = FMIN( MHD_CPAW_WaveLength, x2 );
   if ( MHD_CPAW_Angle2 != 0.0 ) MHD_CPAW_WaveLength = FMIN( MHD_CPAW_WaveLength, x3 );
   MHD_CPAW_WaveNumber    = 2.0*M_PI/MHD_CPAW_WaveLength;
   MHD_CPAW_WaveSpeed     = MHD_CPAW_B0/SQRT(MHD_CPAW_Rho0);
   MHD_CPAW_WaveFrequency = MHD_CPAW_WaveNumber*MHD_CPAW_WaveSpeed;

// (3) reset other general-purpose parameters
//     --> a helper macro PRINT_RESET_PARA is defined in Macro.h
   const double End_T_Default    = 2.0*M_PI/MHD_CPAW_WaveFrequency;
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
      Aux_Message( stdout, "  test problem ID        = %d\n",      TESTPROB_ID            );
      Aux_Message( stdout, "  background density     = % 14.7e\n", MHD_CPAW_Rho0          );
      Aux_Message( stdout, "  perturbation amplitude = % 14.7e\n", MHD_CPAW_A             );
      Aux_Message( stdout, "  background pressure    = % 14.7e\n", MHD_CPAW_P0            );
      Aux_Message( stdout, "  background B field     = % 14.7e\n", MHD_CPAW_B0            );
      Aux_Message( stdout, "  wave speed             = % 14.7e\n", MHD_CPAW_WaveSpeed     );
      Aux_Message( stdout, "  wavelength             = % 14.7e\n", MHD_CPAW_WaveLength    );
      Aux_Message( stdout, "  wavenumber             = % 14.7e\n", MHD_CPAW_WaveNumber    );
      Aux_Message( stdout, "  wave frequency         = % 14.7e\n", MHD_CPAW_WaveFrequency );
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
   const double Rho0       = MHD_CPAW_Rho0;
   const double A          = MHD_CPAW_A;
   const double P0         = MHD_CPAW_P0;
   const double WaveSpeed  = MHD_CPAW_WaveSpeed;
   const double WaveLength = MHD_CPAW_WaveLength;
   const double WaveK      = MHD_CPAW_WaveNumber;
   const double WaveW      = MHD_CPAW_WaveFrequency;
   const double sin_a3     = MHD_CPAW_SinAngle3;
   const double cos_a3     = MHD_CPAW_CosAngle3;
   const double sin_a2     = MHD_CPAW_SinAngle2;
   const double cos_a2     = MHD_CPAW_CosAngle2;

   double kr, r, WaveForm1, WaveForm2, Gamma, DampForm;
   double Dens, Mom0, Mom1, Mom2, MomX, MomY, MomZ, Pres, Eint, Etot;

   r  = cos_a2 * ( x*cos_a3 + y*sin_a3 ) + z*sin_a2;
   kr = WaveK*r;

   Dens = Rho0;
   Pres = P0;

   WaveForm1 =  sin( kr - WaveW*Time );
   WaveForm2 =  cos( kr - WaveW*Time );
   Mom0      =  Dens * WaveSpeed;
   Mom1      = -Mom0 * A * WaveForm1;
   Mom2      = -Mom0 * A * WaveForm2;

   MomX = Mom0*cos_a2*cos_a3 - Mom1*sin_a3 - Mom2*sin_a2*cos_a3;  
   MomY = Mom0*cos_a2*sin_a3 + Mom1*cos_a3 - Mom2*sin_a2*sin_a3;     
   MomZ = Mom0*sin_a2                      + Mom2*cos_a2;       
   
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
// Function    :  SetAFieldIC
// Description :  Function template to initialize the magnetic vector potential
//
// Note        :  1. Invoked by MHD_Init_BField_ByVecPot_Function() using the function pointer
//                   "Init_BField_ByVecPot_User_Ptr", which must be set by a test problem initializer
//                2. This function will be invoked by multiple OpenMP threads when OPENMP is enabled
//                   (unless OPT__INIT_GRID_WITH_OMP is disabled)
//                   --> Please ensure that everything here is thread-safe
//
// Parameter   :  x/y/z     : Target physical coordinates
//                Time      : Target physical time
//                lv        : Target refinement level
//                Component : Component of the output magnetic vector potential
//                            --> Supported components: 'x', 'y', 'z'
//                AuxArray  : Auxiliary array
//                            --> Useless since it is currently fixed to NULL
//
// Return      :  "Component"-component of the magnetic vector potential at (x, y, z, Time)
//-------------------------------------------------------------------------------------------------------
double SetAFieldIC( const double x, const double y, const double z, const double Time,
                    const int lv, const char Component, double AuxArray[] )
{

// short names
   const double A          = MHD_CPAW_A;
   const double B0         = MHD_CPAW_B0;
   const double WaveK      = MHD_CPAW_WaveNumber;
   const double WaveW      = MHD_CPAW_WaveFrequency;
   const double sin_a3     = MHD_CPAW_SinAngle3;
   const double cos_a3     = MHD_CPAW_CosAngle3;
   const double sin_a2     = MHD_CPAW_SinAngle2;
   const double cos_a2     = MHD_CPAW_CosAngle2;

   double kr, r, xx, WaveForm1, WaveForm2;
   double e1x, e1y, e1z, e2x, e2y, e2z;
   double A1, A2, mag_vecpot;
   double xcomp, ycomp, zcomp;

   r  = cos_a2 * ( x*cos_a3 + y*sin_a3 ) + sin_a2*z;
   kr = WaveK*r;
   xx = -x*sin_a3 + y*cos_a3;

   WaveForm1 =  sin( kr - WaveW*Time );
   WaveForm2 =  cos( kr - WaveW*Time );
   A1        =  B0 * A * WaveForm1 / WaveK;
   A2        =  B0 * A * WaveForm2 / WaveK + B0*xx;

   switch ( Component )
   {
      case 'x':  mag_vecpot = -A1*sin_a3-A2*sin_a2*cos_a3;  break;
      case 'y':  mag_vecpot =  A1*cos_a3-A2*sin_a2*sin_a3;  break;
      case 'z':  mag_vecpot =  A2*cos_a2;                   break;
      default :  Aux_Error( ERROR_INFO, "unsupported component (%c) !!\n", Component );
   }

   return mag_vecpot;

} // FUNCTION : SetAFieldIC

//-------------------------------------------------------------------------------------------------------
// Function    :  CheckBField
// Description :  Check the problem-specific magnetic field
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
void CheckBField( real magnetic[], const double x, const double y, const double z, const double Time,
                    const int lv, double AuxArray[] )
{


// short names
   const double A          = MHD_CPAW_A;
   const double B0         = MHD_CPAW_B0;
   const double WaveK      = MHD_CPAW_WaveNumber;
   const double WaveW      = MHD_CPAW_WaveFrequency;
   const double sin_a3     = MHD_CPAW_SinAngle3;
   const double cos_a3     = MHD_CPAW_CosAngle3;
   const double sin_a2     = MHD_CPAW_SinAngle2;
   const double cos_a2     = MHD_CPAW_CosAngle2;

   double kr, r, B1, B2, WaveForm1, WaveForm2;

   r  = cos_a2 * ( x*cos_a3 + y*sin_a3 ) + z*sin_a2;
   kr = WaveK*r;

   WaveForm1 =  sin( kr - WaveW*Time );
   WaveForm2 =  cos( kr - WaveW*Time );
   B1        =  B0 * A * WaveForm1;
   B2        =  B0 * A * WaveForm2;

   magnetic[MAGX] = B0*cos_a2*cos_a3 - B1*sin_a3 - B2*sin_a2*cos_a3;  
   magnetic[MAGY] = B0*cos_a2*sin_a3 + B1*cos_a3 - B2*sin_a2*sin_a3;     
   magnetic[MAGZ] = B0*sin_a2                    + B2*cos_a2;       

} // FUNCTION : CheckBField


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

   const char Prefix[100]     = "MHD_CPAW";
   const OptOutputPart_t Part = OUTPUT_XY;

   Output_L1Error( SetGridIC, CheckBField, Prefix, Part, OUTPUT_PART_X, OUTPUT_PART_Y, OUTPUT_PART_Z );

} // FUNCTION : OutputError
#endif // #ifdef MHD
#endif // #if ( MODEL == HYDRO )



//-------------------------------------------------------------------------------------------------------
// Function    :  Init_TestProb_Hydro_MHD_CPAW
// Description :  Test problem initializer
//
// Note        :  None
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Init_TestProb_Hydro_MHD_CPAW()
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
   Init_BField_ByVecPot_User_Ptr = SetAFieldIC;
   Output_User_Ptr               = OutputError;
#  endif
#  ifdef SUPPORT_HDF5
   Output_HDF5_InputTest_Ptr     = LoadInputTestProb;
#  endif
#  endif // #if ( MODEL == HYDRO )


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : Init_TestProb_Hydro_MHD_CPAW
