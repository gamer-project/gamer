#include "GAMER.h"



// problem-specific global variables
// =======================================================================================
static int    PWave_NWavelength;   // number of plane wave wavelength (will be reset to 3 times input value if PWave_XYZ == 3)
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

static FieldIdx_t PWave_Idx_Phase = Idx_Undefined;    // field index for unwrapped phase
// =======================================================================================

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
   if ( OPT__BC_FLU[f] != BC_FLU_PERIODIC )
      Aux_Error( ERROR_INFO, "must adopt periodic BC for fluid --> reset OPT__BC_FLU* !!\n" );

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
// ReadPara->Add( "KEY_IN_THE_FILE",   &VARIABLE,              DEFAULT,       MIN,              MAX               );
// ********************************************************************************************************************************
   ReadPara->Add( "PWave_NWavelength", &PWave_NWavelength,     2,             1,                NoMax_int         );
   ReadPara->Add( "PWave_Amp",         &PWave_Amp,             1.0,           Eps_double,       NoMax_double      );
   ReadPara->Add( "PWave_Phase0",      &PWave_Phase0,          0.0,           NoMin_double,     NoMax_double      );
   ReadPara->Add( "PWave_XYZ",         &PWave_XYZ,             0,             0,                3                 );
   ReadPara->Add( "PWave_LSR",         &PWave_LSR,             1,             NoMin_int,        NoMax_int         );

   ReadPara->Read( FileName );

   delete ReadPara;

// (1-2) set the default values

// (1-3) check and reset the runtime parameters
   if ( PWave_XYZ == 3  &&  ( amr->BoxSize[0] != amr->BoxSize[1]  ||  amr->BoxSize[0] != amr->BoxSize[2] ) )
      Aux_Error( ERROR_INFO, "simulation domain must be CUBIC in %s test if PWave_XYZ == 3 !!\n", "ELBDM PlaneWave" );

   if ( PWave_XYZ == 3 ) {
      PWave_NWavelength *= 3;
      PRINT_RESET_PARA( PWave_NWavelength, FORMAT_INT, "");
   }


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
      Aux_Message( stdout, "  number of plane wave wavelength  = %d\n",     PWave_NWavelength          );
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
//                3. fluid[Idx_Phase] is used to store the unwrapped phase
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
      Phase = SATAN2 ( Imag, Real );
   }

// set the density
   fluid[DENS] = SQR( Real ) + SQR( Imag );

#  if ( ELBDM_SCHEME == ELBDM_HYBRID )
   if ( amr->use_wave_flag[lv] ) {
#  endif
   fluid[REAL] = Real;
   fluid[IMAG] = Imag;
// set the unwrapped phase
   fluid[PWave_Idx_Phase] = SATAN2( Imag, Real );
#  if ( ELBDM_SCHEME == ELBDM_HYBRID )
   } else {
   fluid[PHAS] = Phase;
   fluid[STUB] = 0.0;
   }
#  endif

} // FUNCTION : SetGridIC



//-------------------------------------------------------------------------------------------------------
// Function    :  AddNewField_PlaneWave
// Description :  Add the unwrapped phase as a problem-specific field
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

   if ( PWave_Idx_Phase == Idx_Undefined )
      PWave_Idx_Phase = AddField( "Phase", FIXUP_FLUX_NO, FIXUP_REST_NO, NORMALIZE_NO, INTERP_FRAC_NO );

} // FUNCTION : AddNewField_PlaneWave



//-------------------------------------------------------------------------------------------------------
// Function    :  Output_UserWorkBeforeOutput_PlaneWave
// Description :  Calculate and update the unwrapped phase field before dumping data
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
//          record the unwrapped phase
            double Phase;

            Phase = ATAN2( amr->patch[FluSg][lv][PID]->fluid[IMAG][k][j][i],
                           amr->patch[FluSg][lv][PID]->fluid[REAL][k][j][i] );

            amr->patch[FluSg][lv][PID]->fluid[PWave_Idx_Phase][k][j][i] = Phase;

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
#  endif // #if ( MODEL == ELBDM )


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : Init_TestProb_ELBDM_PlaneWave
