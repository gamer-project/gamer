#include "GAMER.h"



// problem-specific global variables
// =======================================================================================
static double Gau_v0;       // mean velocity
static double Gau_Width;    // Gaussian width
static double Gau_Center;   // Gaussian center
static int    Gau_XYZ;      // wave propagation direction (0/1/2 --> x/y/z)
static int    Gau_PeriodicN;// periodic boundary condition
                            // (0 = non-periodic, >0 = number of periodic images each side)
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


// warnings
   if ( MPI_Rank == 0 )
   {
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
   ReadPara->Add( "Gau_v0",            &Gau_v0,                1.0,           NoMin_double,     NoMax_double      );
   ReadPara->Add( "Gau_Width",         &Gau_Width,             0.1,           Eps_double,       NoMax_double      );
   ReadPara->Add( "Gau_Center",        &Gau_Center,            NoDef_double,  NoMin_double,     NoMax_double      );
   ReadPara->Add( "Gau_XYZ",           &Gau_XYZ,               0,             0,                2                 );
   ReadPara->Add( "Gau_PeriodicN",     &Gau_PeriodicN,         0,             0,                NoMax_int         );

   ReadPara->Read( FileName );

   delete ReadPara;

// (1-2) set the default values
   if ( Gau_Center == NoDef_double )   Gau_Center = amr->BoxCenter[Gau_XYZ];

// (1-3) check the runtime parameters
   if ( MPI_Rank == 0 )
   {
      if ( OPT__BC_FLU[2*Gau_XYZ] == BC_FLU_PERIODIC  &&  Gau_PeriodicN == 0 )
         Aux_Message( stderr, "WARNING : Gau_PeriodicN should be >0 when adopting periodic BC for fluid !!\n" );
   }


// (2) set the problem-specific derived parameters


// (3) reset other general-purpose parameters
//     --> a helper macro PRINT_RESET_PARA is defined in Macro.h
   const long   End_Step_Default = __INT_MAX__;
   const double End_T_Default    = 0.50*amr->BoxSize[Gau_XYZ]/fabs( Gau_v0 );

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
      Aux_Message( stdout, "  test problem ID                    = %d\n",     TESTPROB_ID   );
      Aux_Message( stdout, "  mean velocity                      = %14.7e\n", Gau_v0        );
      Aux_Message( stdout, "  Gaussian width                     = %14.7e\n", Gau_Width     );
      Aux_Message( stdout, "  Gaussian center                    = %14.7e\n", Gau_Center    );
      Aux_Message( stdout, "  propagation direction              = %d\n",     Gau_XYZ       );
      Aux_Message( stdout, "  number of periodic image each side = %d\n",     Gau_PeriodicN );
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

   double r;
   switch ( Gau_XYZ )
   {
      case 0: r = x;    break;
      case 1: r = y;    break;
      case 2: r = z;    break;

      default : Aux_Error( ERROR_INFO, "incorrect Gau_XYZ (%d) !!\n", Gau_XYZ );
      break;
   }

   const double Gau_Const1 = 1.0 + pow(  Time / ( ELBDM_ETA*SQR(Gau_Width) ), 2.0  );
   const double Gau_Theta1 = -0.5*acos(  pow( Gau_Const1, -0.5 )  );
   double Re=0.0, Im=0.0;

// n=0, m=0: original wave packet
// n>0, m=0/1: images for periodic BC on the plus(+)/minus(-) direction
   for (int n=0; n<Gau_PeriodicN+1; n++) {
   for (int m=0; m<((n==0)?1:2); m++) {
      const double Center     = Gau_Center + n*(1-2*m)*amr->BoxSize[Gau_XYZ];
      const double dr1        = r -     Gau_v0*Time - Center;
      const double dr2        = r - 0.5*Gau_v0*Time - Center;
      const double Gau_Const2 = pow( SQR(Gau_Width)*M_PI*Gau_Const1, -0.25 )
                                *exp(  -0.5*pow( dr1/Gau_Width, 2.0 )/Gau_Const1  );
      const double Gau_Theta2 = 0.5*pow( dr1, 2.0 )*ELBDM_ETA*Time/(  pow( ELBDM_ETA*SQR(Gau_Width), 2.0) + SQR(Time)  )
                                + Gau_v0*ELBDM_ETA*dr2;

      Re += Gau_Const2*cos( Gau_Theta1 + Gau_Theta2 );
      Im += Gau_Const2*sin( Gau_Theta1 + Gau_Theta2 );
   }}

   fluid[DENS] = SQR( Re ) + SQR( Im );

#  if ( ELBDM_SCHEME == ELBDM_HYBRID )
   if ( amr->use_wave_flag[lv] ) {
#  endif
   fluid[REAL] = Re;
   fluid[IMAG] = Im;
#  if  ( ELBDM_SCHEME == ELBDM_HYBRID )
   } else {
   fluid[PHAS] = SATAN2( Im, Re );
   fluid[STUB] = 0.0;
   }
#  endif

} // FUNCTION : SetGridIC



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

   const char Prefix[100]     = "Gaussian";
   const OptOutputPart_t Part = OUTPUT_X + Gau_XYZ;

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
// Function    :  Init_TestProb_ELBDM_GaussianWavePacket
// Description :  Test problem initializer
//
// Note        :  None
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Init_TestProb_ELBDM_GaussianWavePacket()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );


// validate the compilation flags and runtime parameters
   Validate();


#  if ( MODEL == ELBDM )
// set the problem-specific runtime parameters
   SetParameter();


   Init_Function_User_Ptr = SetGridIC;
   BC_User_Ptr            = BC;
   Output_User_Ptr        = OutputError;
#  endif // #if ( MODEL == ELBDM )


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : Init_TestProb_ELBDM_GaussianWavePacket
