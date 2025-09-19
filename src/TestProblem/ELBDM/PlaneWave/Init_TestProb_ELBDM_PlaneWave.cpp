#include "GAMER.h"



// problem-specific global variables
// =======================================================================================
static double PWave_NWavelength;   // number of plane wave wavelength (will be reset to 3 times input value if PWave_XYZ == 3)
static double PWave_Amp;           // plane wave amplitude
static double PWave_Phase0;        // plane wave phase constant
static int    PWave_XYZ;           // plane wave direction (0/1/2/3 --> x/y/z/diagonal)
static int    PWave_LSR;           // plane wave direction (<0/0/>0 --> Left-moving/Standing/Right-moving)

static double PWave_Lambda;        // plane wave wavelength
static double PWave_Period;        // plane wave period
static double PWave_WaveK;         // plane wave wavenumber
static double PWave_WaveW;         // plane wave angular frequency
static double PWave_PhaseV;        // plane wave phase velocity
static double PWave_GroupV;        // plane wave group velocity

static FieldIdx_t PWave_Idx_WrappedPhase = Idx_Undefined;    // field index for the wrapped phase
// =======================================================================================

static void OutputError();
static void BC( real Array[], const int ArraySize[], real fluid[], const int NVar_Flu,
                const int GhostSize, const int idx[], const double pos[], const double Time,
                const int lv, const int TFluVarIdxList[], double AuxArray[] );




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
#  if ( MODEL != ELBDM )
   Aux_Error( ERROR_INFO, "MODEL != ELBDM !!\n" );
#  endif

#  ifdef GRAVITY
   Aux_Error( ERROR_INFO, "GRAVITY must be disabled !!\n" );
#  endif

#  ifdef PARTICLE
   Aux_Error( ERROR_INFO, "PARTICLE must be disabled !!\n" );
#  endif

   for (int f=0; f<6; f++)
   if ( OPT__BC_FLU[f] != BC_FLU_PERIODIC  &&  OPT__BC_FLU[f] != BC_FLU_USER )
      Aux_Error( ERROR_INFO, "must adopt periodic or user BC for fluid --> reset OPT__BC_FLU* !!\n" );

   if ( NCOMP_PASSIVE_USER != 1 )
      Aux_Error( ERROR_INFO, "please set NCOMP_PASSIVE_USER to 1 !!\n" );


// warnings
   if ( MPI_Rank == 0 )
   {
#     ifndef FLOAT8
      Aux_Message( stderr, "WARNING : it's recommended to enable FLOAT8 for this test !!\n" );
#     endif

      if ( !OPT__OUTPUT_USER )
         Aux_Message( stderr, "WARNING : it's recommended to enable OPT__OUTPUT_USER !!\n" );
   }


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Validating test problem %d ... done\n", TESTPROB_ID );

} // FUNCTION : Validate



#if ( MODEL == ELBDM )
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
// LOAD_PARA( load_mode, "KEY_IN_THE_FILE",   &VARIABLE,              DEFAULT,       MIN,              MAX               );
// ************************************************************************************************************************
   LOAD_PARA( load_mode, "PWave_NWavelength", &PWave_NWavelength,     2.0,           Eps_double,       NoMax_double      );
   LOAD_PARA( load_mode, "PWave_Amp",         &PWave_Amp,             1.0,           Eps_double,       NoMax_double      );
   LOAD_PARA( load_mode, "PWave_Phase0",      &PWave_Phase0,          0.0,           NoMin_double,     NoMax_double      );
   LOAD_PARA( load_mode, "PWave_XYZ",         &PWave_XYZ,             0,             0,                3                 );
   LOAD_PARA( load_mode, "PWave_LSR",         &PWave_LSR,             1,             NoMin_int,        NoMax_int         );

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

// (1-3) check and reset the runtime parameters
   if ( PWave_XYZ == 3  &&  ( amr->BoxSize[0] != amr->BoxSize[1]  ||  amr->BoxSize[0] != amr->BoxSize[2] ) )
      Aux_Error( ERROR_INFO, "simulation domain must be CUBIC in %s test if PWave_XYZ == 3 !!\n", "ELBDM PlaneWave" );

   if ( PWave_XYZ == 3 ) {
      PWave_NWavelength *= 3;
      PRINT_RESET_PARA( PWave_NWavelength, FORMAT_REAL, "");
   }

   if ( ! Mis_CompareRealValue( PWave_NWavelength, round(PWave_NWavelength), NULL, false ) )
      for (int f=((PWave_XYZ < 3)?(2*PWave_XYZ):0); f<((PWave_XYZ < 3)?(2*PWave_XYZ+2):6); f++)
         if ( OPT__BC_FLU[f] != BC_FLU_USER )
            Aux_Error( ERROR_INFO, "OPT__BC_FLU[%d] must be %d for non-integer PWave_NWavelength (=%13.7e) !!\n",
                                   f, BC_FLU_USER, PWave_NWavelength );

#  if ( ELBDM_SCHEME == ELBDM_HYBRID )
   if ( ELBDM_FIRST_WAVE_LEVEL > MAX_LEVEL )
   {
      if ( PWave_LSR == 0 )
      {
//       there are nodes in a standing wave
         Aux_Error( ERROR_INFO, "Standing wave (PWave_LSR = %d) cannot work in the fluid scheme !!\n", PWave_LSR );
      }
      else // if ( PWave_LSR == 0 )
      {
//       the unwrapped phase of a travelling wave is not periodic
         for (int f=((PWave_XYZ < 3)?(2*PWave_XYZ):0); f<((PWave_XYZ < 3)?(2*PWave_XYZ+2):6); f++)
            if ( OPT__BC_FLU[f] != BC_FLU_USER )
               Aux_Error( ERROR_INFO, "OPT__BC_FLU[%d] must be %d for travelling wave in fluid scheme !!\n",
                                      f, BC_FLU_USER );
      } // if ( PWave_LSR == 0 ) ... else
   }
   else // if ( ELBDM_FIRST_WAVE_LEVEL > MAX_LEVEL )
   {
//    the phase on the wave level is wrapped
      if ( !ELBDM_MATCH_PHASE )
         Aux_Error( ERROR_INFO, "ELBDM_MATCH_PHASE should be enabled to make the phase on fluid levels continuous !!\n" );
   } // if ( ELBDM_FIRST_WAVE_LEVEL > MAX_LEVEL ) ... else
#  endif // # if ( ELBDM_SCHEME == ELBDM_HYBRID )


// (2) set the problem-specific derived parameters
   PWave_Lambda = ( PWave_XYZ == 3 ) ? amr->BoxSize[0]*sqrt(3.0)/PWave_NWavelength : amr->BoxSize[PWave_XYZ]/PWave_NWavelength;
   PWave_WaveK  = 2.0*M_PI/PWave_Lambda;
   PWave_WaveW  = 0.5*SQR( PWave_WaveK )/ELBDM_ETA;
   PWave_Period = 2.0*M_PI/PWave_WaveW;
   PWave_PhaseV = 0.5*PWave_WaveK/ELBDM_ETA;
   PWave_GroupV = PWave_WaveK/ELBDM_ETA;



// (3) reset other general-purpose parameters
//     --> a helper macro PRINT_RESET_PARA is defined in Macro.h
   const long   End_Step_Default = __INT_MAX__;
   const double End_T_Default    = 6.0*PWave_Period; // 6 periods

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
      Aux_Message( stdout, "  test problem ID                  = %d\n",     TESTPROB_ID                );
      Aux_Message( stdout, "  number of plane wave wavelength  = %13.7e\n", PWave_NWavelength          );
      Aux_Message( stdout, "  plane wave wavelength            = %13.7e\n", PWave_Lambda               );
      Aux_Message( stdout, "  plane wave amplitude             = %13.7e\n", PWave_Amp                  );
      Aux_Message( stdout, "  plane wave phase constant        = %13.7e\n", PWave_Phase0               );
      Aux_Message( stdout, "  plane wave wavenumber            = %13.7e\n", PWave_WaveK                );
      Aux_Message( stdout, "  plane wave angular frequency     = %13.7e\n", PWave_WaveW                );
      Aux_Message( stdout, "  plane wave period                = %13.7e\n", PWave_Period               );
      Aux_Message( stdout, "  standing wave                    = %s\n",    ( PWave_LSR == 0 ) ? "YES" : "NO" );
      Aux_Message( stdout, "  plane wave direction             = %s%s\n",  ( PWave_LSR == 0 ) ? ""  :
                                                                           ( PWave_LSR >  0 ) ? "+" : "-",
                                                                           ( PWave_XYZ == 0 ) ? "x" :
                                                                           ( PWave_XYZ == 1 ) ? "y" :
                                                                           ( PWave_XYZ == 2 ) ? "z" : "diagonal" );
      if ( PWave_LSR != 0 ) {
      Aux_Message( stdout, "  plane wave phase velocity        = %13.7e\n", PWave_PhaseV               );
      Aux_Message( stdout, "  plane wave group velocity        = %13.7e\n", PWave_GroupV               );
      }
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
//                3. fluid[Idx_WrappedPhase] is used to store the wrapped phase
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

   double r, PhaseR, PhaseL, Real, Imag, Phase;
   switch ( PWave_XYZ )
   {
      case 0 : r = x;                              break;
      case 1 : r = y;                              break;
      case 2 : r = z;                              break;
      case 3 : r = 1.0/sqrt(3.0)*( x + y + z );    break;
      default : Aux_Error( ERROR_INFO, "incorrect parameter \"%s = %d [0/1/2/3]\" !!\n", "PWave_XYZ", PWave_XYZ );
      break;
   }

   PhaseR  =  PWave_WaveK*r - PWave_WaveW*Time + PWave_Phase0;
   PhaseL  = -PWave_WaveK*r - PWave_WaveW*Time + PWave_Phase0;

// set the real and imaginary parts
   if      ( PWave_LSR > 0 ) { // right-moving wave
      Real  = PWave_Amp*cos( PhaseR );
      Imag  = PWave_Amp*sin( PhaseR );
      Phase = PhaseR;
   }
   else if ( PWave_LSR < 0 ) { // left-moving wave
      Real  = PWave_Amp*cos( PhaseL );
      Imag  = PWave_Amp*sin( PhaseL );
      Phase = PhaseL;
   }
   else {                      // standing wave (PWave_LSR == 0 )
      Real  = 0.5*(  PWave_Amp*cos( PhaseR ) + PWave_Amp*cos( PhaseL )  );
      Imag  = 0.5*(  PWave_Amp*sin( PhaseR ) + PWave_Amp*sin( PhaseL )  );
      Phase = ( cos( PWave_WaveK*r ) > 0 ) ? 0.5*( PhaseR + PhaseL ) : M_PI+0.5*( PhaseR + PhaseL );
   }

// set the density
   fluid[DENS] = SQR( Real ) + SQR( Imag );

#  if ( ELBDM_SCHEME == ELBDM_HYBRID )
   if ( amr->use_wave_flag[lv] ) {
#  endif
   fluid[REAL] = Real;
   fluid[IMAG] = Imag;
// set the wrapped phase
   fluid[PWave_Idx_WrappedPhase] = SATAN2( Imag, Real );
#  if ( ELBDM_SCHEME == ELBDM_HYBRID )
   } else {
   fluid[PHAS] = Phase;
   fluid[STUB] = 0.0;
// set the wrapped phase
   fluid[PWave_Idx_WrappedPhase] = SATAN2( SIN(Phase), COS(Phase) );
   }
#  endif

} // FUNCTION : SetGridIC



//-------------------------------------------------------------------------------------------------------
// Function    :  AddNewField_PlaneWave
// Description :  Add the wrapped phase as a problem-specific field
//
// Note        :  1. Ref: https://github.com/gamer-project/gamer/wiki/Adding-New-Simulations#v-add-problem-specific-grid-fields-and-particle-attributes
//                2. Invoke AddField() for each of the problem-specific field:
//                   --> Field label sent to AddField() will be used as the output name of the field
//                   --> Field index returned by AddField() can be used to access the field data
//                3. Pre-declared field indices are put in Field.h
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void AddNewField_PlaneWave()
{

   if ( PWave_Idx_WrappedPhase == Idx_Undefined )
      PWave_Idx_WrappedPhase = AddField( "WrappedPhase", FIXUP_FLUX_NO, FIXUP_REST_NO, NORMALIZE_NO, INTERP_FRAC_NO );

} // FUNCTION : AddNewField_PlaneWave



//-------------------------------------------------------------------------------------------------------
// Function    :  Output_UserWorkBeforeOutput_PlaneWave
// Description :  Calculate and update the wrapped phase field before dumping data
//
// Note        :  1. Invoked by Output_DumpData() using the function pointer "Output_UserWorkBeforeOutput_Ptr"
//
// Parameter   :
//
// Return      :
//-------------------------------------------------------------------------------------------------------
void Output_UserWorkBeforeOutput_PlaneWave()
{

   for (int lv=0; lv<NLEVEL; lv++)
   {
      int FluSg = amr->FluSg[lv];

      for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
      {
         for (int k=0; k<PATCH_SIZE; k++)
         for (int j=0; j<PATCH_SIZE; j++)
         for (int i=0; i<PATCH_SIZE; i++)
         {
//          record the wrapped phase
            double Phase;

#           if ( ELBDM_SCHEME == ELBDM_HYBRID )
            if ( amr->use_wave_flag[lv] )
            {
#           endif
               Phase = SATAN2( amr->patch[FluSg][lv][PID]->fluid[IMAG][k][j][i],
                               amr->patch[FluSg][lv][PID]->fluid[REAL][k][j][i] );
#           if ( ELBDM_SCHEME == ELBDM_HYBRID )
            } // if ( amr->use_wave_flag[lv] )
            else
            {
               Phase = SATAN2( SIN(amr->patch[FluSg][lv][PID]->fluid[PHAS][k][j][i]),
                               COS(amr->patch[FluSg][lv][PID]->fluid[PHAS][k][j][i]) );
            } // if ( amr->use_wave_flag[lv] ) ... else
#           endif

            amr->patch[FluSg][lv][PID]->fluid[PWave_Idx_WrappedPhase][k][j][i] = Phase;

         } // i,j,k
      } // for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
   } // for (int lv=0; lv<NLEVEL; lv++)

} // FUNCTION : Output_UserWorkBeforeOutput_PlaneWave



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

   const char Prefix[MAX_STRING] = "PlaneWave";
   const OptOutputPart_t Part    = ( PWave_XYZ == 3 ) ? OUTPUT_DIAG : (OUTPUT_X + PWave_XYZ);

   Output_L1Error( SetGridIC, NULL, Prefix, Part, OUTPUT_PART_X, OUTPUT_PART_Y, OUTPUT_PART_Z );

} // FUNCTION : OutputError



//-------------------------------------------------------------------------------------------------------
// Function    :  BC
// Description :  Set the extenral boundary condition to the analytical solution
//
// Note        :  1. Linked to the function pointer "BC_User_Ptr"
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
//                TFluVarIdxList : List recording the target fluid variable indices ( = [0 ... NCOMP_TOTAL-1] )
//                AuxArray       : Auxiliary array
//
// Return      :  fluid
//-------------------------------------------------------------------------------------------------------
void BC( real Array[], const int ArraySize[], real fluid[], const int NVar_Flu,
         const int GhostSize, const int idx[], const double pos[], const double Time,
         const int lv, const int TFluVarIdxList[], double AuxArray[] )
{

// simply call the IC function
   SetGridIC( fluid, pos[0], pos[1], pos[2], Time, lv, AuxArray );

} // FUNCTION : BC
#endif // #if ( MODEL == ELBDM )



//-------------------------------------------------------------------------------------------------------
// Function    :  Init_TestProb_ELBDM_PlaneWave
// Description :  Test problem initializer
//
// Note        :  None
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Init_TestProb_ELBDM_PlaneWave()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );


// validate the compilation flags and runtime parameters
   Validate();


#  if ( MODEL == ELBDM )
// set the problem-specific runtime parameters
   SetParameter();


   Init_Function_User_Ptr          = SetGridIC;
   Init_Field_User_Ptr             = AddNewField_PlaneWave;
   Output_User_Ptr                 = OutputError;
   Output_UserWorkBeforeOutput_Ptr = Output_UserWorkBeforeOutput_PlaneWave;
#  ifdef SUPPORT_HDF5
   Output_HDF5_InputTest_Ptr       = LoadInputTestProb;
#  endif
   BC_User_Ptr                     = BC;
#  endif // #if ( MODEL == ELBDM )


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : Init_TestProb_ELBDM_PlaneWave
