#include "GAMER.h"



// problem-specific global variables
// =======================================================================================
static int LSS_InitMode;   // initialization mode: 1=density-only, 2=real and imaginary parts or phase and density
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
#  if ( MODEL != ELBDM)
   Aux_Error( ERROR_INFO, "MODEL != ELBDM !!\n" );
#  endif

#  ifndef GRAVITY
   Aux_Error( ERROR_INFO, "GRAVITY must be enabled !!\n" );
#  endif

#  ifndef COMOVING
   Aux_Error( ERROR_INFO, "COMOVING must be enabled !!\n" );
#  endif

#  ifdef GRAVITY
   if ( OPT__BC_FLU[0] != BC_FLU_PERIODIC  ||  OPT__BC_POT != BC_POT_PERIODIC )
      Aux_Error( ERROR_INFO, "must adopt periodic BC for this test !!\n" );
#  endif

   if ( OPT__INIT == INIT_BY_FUNCTION )
      Aux_Error( ERROR_INFO, "OPT__INIT=INIT_BY_FUNCTION (1) is not supported for this test !!\n" );


// warnings


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
   LOAD_PARA( load_mode, "LSS_InitMode",      &LSS_InitMode,          1,             1,                2                 );

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

// (1-3) check the runtime parameters


// (2) set the problem-specific derived parameters


// (3) reset other general-purpose parameters
//     --> a helper macro PRINT_RESET_PARA is defined in Macro.h
   const long   End_Step_Default = __INT_MAX__;
   const double End_T_Default    = 1.0;

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
      Aux_Message( stdout, "  test problem ID     = %d\n", TESTPROB_ID  );
      Aux_Message( stdout, "  initialization mode = %d\n", LSS_InitMode );
      Aux_Message( stdout, "=============================================================================\n" );
   }


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Setting runtime parameters ... done\n" );

} // FUNCTION : SetParameter



//-------------------------------------------------------------------------------------------------------
// Function    :  Init_ByFile_ELBDM_LSS
// Description :  Function to actually set the fluid field from the input uniform-mesh array
//
// Note        :  1. Invoked by Init_ByFile_AssignData() using the function pointer Init_ByFile_User_Ptr()
//                   --> The function pointer may be reset by various test problem initializers, in which case
//                       this funtion will become useless
//                2. One can use LSS_InitMode to support different data formats
//                3. For ELBDM_SCHEME == ELBDM_WAVE this function expects:
//                       LSS_InitMode == 1: Density
//                       LSS_InitMode == 2: Real and imaginary part
//                   For ELBDM_SCHEME == ELBDM_HYBRID this function expects:
//                       LSS_InitMode == 1: Density
//                       LSS_InitMode == 2: Density and phase
// Parameter   :  fluid_out : Fluid field to be set
//                fluid_in  : Fluid field loaded from the uniform-mesh array (UM_IC)
//                nvar_in   : Number of variables in fluid_in
//                x/y/z     : Target physical coordinates
//                Time      : Target physical time
//                lv        : Target AMR level
//                AuxArray  : Auxiliary array
//
// Return      :  fluid_out
//-------------------------------------------------------------------------------------------------------
void Init_ByFile_ELBDM_LSS( real fluid_out[], const real fluid_in[], const int nvar_in,
                            const double x, const double y, const double z, const double Time,
                            const int lv, double AuxArray[] )
{

   double Re, Im, De;

#  if ( ELBDM_SCHEME == ELBDM_HYBRID )
   double Ph;
#  endif

   switch ( LSS_InitMode )
   {
      case 1:
      {
         if ( nvar_in != 1 )  Aux_Error( ERROR_INFO, "nvar_in (%d) != 1 for LSS_InitMode 1 !!\n", nvar_in );

         const double AveDens     = 1.0;        // assuming background density = 1.0
         const double GrowingFrac = 3.0/5.0;    // growing-mode amplitude = total amplitude * 3/5

         Re = sqrt( (fluid_in[0]-AveDens )/GrowingFrac + AveDens );
         Im = 0.0;   // constant phase
         De = SQR( Re ) + SQR( Im );

#        if ( ELBDM_SCHEME == ELBDM_HYBRID )
         Ph = 0.0;
#        endif

         break;
      }

      case 2:
      {
         if ( nvar_in != 2 )  Aux_Error( ERROR_INFO, "nvar_in (%d) != 2 for LSS_InitMode 2 !!\n", nvar_in );

//       ELBDM_WAVE   expects real and imaginary parts
//       ELBDM_HYBRID expects density and phase
#        if   ( ELBDM_SCHEME == ELBDM_WAVE )
         Re = fluid_in[0];
         Im = fluid_in[1];
         De = SQR( Re ) + SQR( Im );
#        elif ( ELBDM_SCHEME == ELBDM_HYBRID )
         De = fluid_in[0];
         Ph = fluid_in[1];
         Re = sqrt( De )*cos( Ph );
         Im = sqrt( De )*sin( Ph );
#        else
#        error : ERROR : unsupported ELBDM_SCHEME !!
#        endif
         break;
      }

      default:
         Aux_Error( ERROR_INFO, "unsupported initialization mode (%s = %d) !!\n",
                    "LSS_InitMode", LSS_InitMode );
   } // switch ( LSS_InitMode )

#  if ( ELBDM_SCHEME == ELBDM_HYBRID )
   if ( amr->use_wave_flag[lv] ) {
#  endif
   fluid_out[DENS] = De;
   fluid_out[REAL] = Re;
   fluid_out[IMAG] = Im;
#  if ( ELBDM_SCHEME == ELBDM_HYBRID )
   } else {
   fluid_out[DENS] = De;
   fluid_out[PHAS] = Ph;
   fluid_out[STUB] = 0.0;
   }
#  endif

} // FUNCTION : Init_ByFile_ELBDM_LSS
#endif // #if ( MODEL == ELBDM )



//-------------------------------------------------------------------------------------------------------
// Function    :  Init_TestProb_ELBDM_LSS
// Description :  Test problem initializer
//
// Note        :  None
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Init_TestProb_ELBDM_LSS()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );


// validate the compilation flags and runtime parameters
   Validate();


#  if ( MODEL == ELBDM )
// set the problem-specific runtime parameters
   SetParameter();


   Init_ByFile_User_Ptr      = Init_ByFile_ELBDM_LSS;
#  ifdef SUPPORT_HDF5
   Output_HDF5_InputTest_Ptr = LoadInputTestProb;
#  endif
#  endif // #if ( MODEL == ELBDM )


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : Init_TestProb_ELBDM_LSS
