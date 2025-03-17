#include "GAMER.h"



// problem-specific global variables
// =======================================================================================
static double Blast_Dens_Bg;        // background mass density
static double Blast_Pres_Bg;        // background pressure
static double Blast_Pres_Exp;       // explosion pressure
static double Blast_Radius;         // explosion radius
static double Blast_Center[3];      // explosion center
#ifdef MHD
static double Blast_BField;         // magnetic field strength along the diagonal direction
       double Blast_ResetB_amp;     // amplitude (for resetting magnetic field)
       double Blast_ResetB_r0;      // scale radius
       double Blast_ResetB_tmin;    // starting time
       double Blast_ResetB_tmax;    // ending time
static bool   Blast_ResetB_VecPot;  // use vector potential to reset magnetic field (recommended)
#endif
// =======================================================================================

// problem-specific function prototypes
bool Flag_BlastWave( const int i, const int j, const int k, const int lv, const int PID, const double *Threshold );
#ifdef MHD
double MHD_ResetByUser_VecPot_BlastWave( const double x, const double y, const double z, const double Time,
                                         const double dt, const int lv, const char Component, double AuxArray[] );
double MHD_ResetByUser_BField_BlastWave( const double x, const double y, const double z, const double Time,
                                         const double dt, const int lv, const char Component, double AuxArray[], const double B_in,
                                         const bool UseVecPot, const real *Ax, const real *Ay, const real *Az,
                                         const int i, const int j, const int k );
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


#  if ( MODEL != HYDRO )
   Aux_Error( ERROR_INFO, "MODEL != HYDRO !!\n" );
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

   if ( !OPT__INIT_RESTRICT )
      Aux_Error( ERROR_INFO, "OPT__INIT_RESTRICT must be enabled !!\n" );


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
// ***************************************************************************************************************************
// LOAD_PARA( load_mode, "KEY_IN_THE_FILE",      &VARIABLE,               DEFAULT,      MIN,              MAX               );
// ***************************************************************************************************************************
   LOAD_PARA( load_mode, "Blast_Dens_Bg",        &Blast_Dens_Bg,         -1.0,          Eps_double,       NoMax_double      );
   LOAD_PARA( load_mode, "Blast_Pres_Bg",        &Blast_Pres_Bg,         -1.0,          Eps_double,       NoMax_double      );
   LOAD_PARA( load_mode, "Blast_Pres_Exp",       &Blast_Pres_Exp,        -1.0,          Eps_double,       NoMax_double      );
   LOAD_PARA( load_mode, "Blast_Radius",         &Blast_Radius,          -1.0,          Eps_double,       NoMax_double      );
   LOAD_PARA( load_mode, "Blast_Center_X",       &Blast_Center[0],       -1.0,          NoMin_double,     amr->BoxSize[0]   );
   LOAD_PARA( load_mode, "Blast_Center_Y",       &Blast_Center[1],       -1.0,          NoMin_double,     amr->BoxSize[1]   );
   LOAD_PARA( load_mode, "Blast_Center_Z",       &Blast_Center[2],       -1.0,          NoMin_double,     amr->BoxSize[2]   );
#  ifdef MHD
   LOAD_PARA( load_mode, "Blast_BField",         &Blast_BField,           5.0e-2,       NoMin_double,     NoMax_double      );
   LOAD_PARA( load_mode, "Blast_ResetB_amp",     &Blast_ResetB_amp,       1.0e2,        0.0,              NoMax_double      );
   LOAD_PARA( load_mode, "Blast_ResetB_r0",      &Blast_ResetB_r0,        1.0e-2,       Eps_double,       NoMax_double      );
   LOAD_PARA( load_mode, "Blast_ResetB_tmin",    &Blast_ResetB_tmin,      1.0e-3,       0.0,              NoMax_double      );
   LOAD_PARA( load_mode, "Blast_ResetB_tmax",    &Blast_ResetB_tmax,      2.0e-3,       0.0,              NoMax_double      );
   LOAD_PARA( load_mode, "Blast_ResetB_VecPot",  &Blast_ResetB_VecPot,    true,         Useless_bool,     Useless_bool      );
#  endif

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

// set the default explosion center
   for (int d=0; d<3; d++)
      if ( Blast_Center[d] < 0.0 )  Blast_Center[d] = 0.5*amr->BoxSize[d];


// (2) set the problem-specific derived parameters


// (3) reset other general-purpose parameters
//     --> a helper macro PRINT_RESET_PARA is defined in Macro.h
   const double End_T_Default    = 5.0e-3;
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
//    assuming EOS_GAMMA (must not invoke any EoS routine here since it has not been initialized)
      const double ExpVol  = 4.0*M_PI/3.0*CUBE(Blast_Radius);
      const double ExpEngy = Blast_Pres_Exp/(GAMMA-1.0)*ExpVol;
#     if ( EOS != EOS_GAMMA )
      Aux_Message( stderr, "WARNING : the total explosion energy below assumes EOS_GAMMA !!\n" );
#     endif

      Aux_Message( stdout, "=============================================================================\n" );
      Aux_Message( stdout, "  test problem ID           = %d\n",     TESTPROB_ID );
      Aux_Message( stdout, "  background mass density   = %13.7e\n", Blast_Dens_Bg );
      Aux_Message( stdout, "  background pressure       = %13.7e\n", Blast_Pres_Bg );
      Aux_Message( stdout, "  explosion pressure        = %13.7e\n", Blast_Pres_Exp );
      Aux_Message( stdout, "  total explosion energy    = %13.7e (assuming constant-gamma EoS)\n", ExpEngy );
      Aux_Message( stdout, "  explosion radius          = %13.7e\n", Blast_Radius );
      Aux_Message( stdout, "  explosion center          = (%13.7e, %13.7e, %13.7e)\n", Blast_Center[0], Blast_Center[1],
                                                                                       Blast_Center[2] );
#     ifdef MHD
      Aux_Message( stdout, "  magnetic field strength   = %13.7e\n", Blast_BField );
      Aux_Message( stdout, "  resetting magnetic field  = %d\n",     OPT__RESET_FLUID );
      if ( OPT__RESET_FLUID ) {
      Aux_Message( stdout, "     amplitude              = %13.7e\n", Blast_ResetB_amp );
      Aux_Message( stdout, "     scale radius           = %13.7e\n", Blast_ResetB_r0 );
      Aux_Message( stdout, "     starting time          = %13.7e\n", Blast_ResetB_tmin );
      Aux_Message( stdout, "     ending time            = %13.7e\n", Blast_ResetB_tmax );
      Aux_Message( stdout, "     use vector potential   = %s\n",     (Blast_ResetB_VecPot)?"YES":"NO" );
      }
#     endif
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

   const double r = SQRT( SQR(x-Blast_Center[0]) + SQR(y-Blast_Center[1]) + SQR(z-Blast_Center[2]) );
   double Dens, MomX, MomY, MomZ, Pres, Eint, Etot;

   Dens = Blast_Dens_Bg;
   Pres = ( r <= Blast_Radius ) ? Blast_Pres_Exp : Blast_Pres_Bg;
#  ifdef SRHD
   real Prim[NCOMP_TOTAL];
   Prim[0] = Dens;
   Prim[1] = 0.0;
   Prim[2] = 0.0;
   Prim[3] = 0.0;
   Prim[4] = Pres;
   Hydro_Pri2Con( Prim, fluid, NULL_BOOL, NULL_INT, NULL,
                  EoS_DensPres2Eint_CPUPtr, EoS_Temp2HTilde_CPUPtr, EoS_HTilde2Temp_CPUPtr,
                  EoS_AuxArray_Flt, EoS_AuxArray_Int, h_EoS_Table, NULL );

#  else // #ifdef SRHD

   MomX = 0.0;
   MomY = 0.0;
   MomZ = 0.0;
   Eint = EoS_DensPres2Eint_CPUPtr( Dens, Pres, NULL, EoS_AuxArray_Flt,
                                    EoS_AuxArray_Int, h_EoS_Table );   // assuming EoS requires no passive scalars
   Etot = Hydro_ConEint2Etot( Dens, MomX, MomY, MomZ, Eint, 0.0 );     // do NOT include magnetic energy here

// set the output array
   fluid[DENS] = Dens;
   fluid[MOMX] = MomX;
   fluid[MOMY] = MomY;
   fluid[MOMZ] = MomZ;
   fluid[ENGY] = Etot;
#  endif // #ifdef SRHD ... else ...

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

   magnetic[MAGX] = Blast_BField / sqrt(3.0);
   magnetic[MAGY] = Blast_BField / sqrt(3.0);
   magnetic[MAGZ] = Blast_BField / sqrt(3.0);

} // FUNCTION : SetBFieldIC
#endif // #ifdef MHD
#endif // #if ( MODEL == HYDRO )



//-------------------------------------------------------------------------------------------------------
// Function    :  Init_TestProb_Hydro_BlastWave
// Description :  Test problem initializer
//
// Note        :  None
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Init_TestProb_Hydro_BlastWave()
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
   MHD_ResetByUser_BField_Ptr    = MHD_ResetByUser_BField_BlastWave;
   MHD_ResetByUser_VecPot_Ptr    = (Blast_ResetB_VecPot) ? MHD_ResetByUser_VecPot_BlastWave : NULL;
#  endif
   Flag_User_Ptr                 = Flag_BlastWave;
#  ifdef SUPPORT_HDF5
   Output_HDF5_InputTest_Ptr     = LoadInputTestProb;
#  endif
#  endif // #if ( MODEL == HYDRO )


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : Init_TestProb_Hydro_BlastWave
