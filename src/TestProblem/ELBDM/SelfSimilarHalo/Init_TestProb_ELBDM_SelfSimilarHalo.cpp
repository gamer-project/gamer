#include "GAMER.h"



// problem-specific global variables
// =======================================================================================
static char    SelSimHalo_Filename[MAX_STRING]; // filename of the profile table

static double *SelSimHalo_Prof = NULL;          // radial profiles [radius/mass_density/phase]
static int     SelSimHalo_NBin;                 // number of radial bins in the profile table
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
#  if ( MODEL != ELBDM )
   Aux_Error( ERROR_INFO, "MODEL != ELBDM !!\n" );
#  endif

#  ifndef GRAVITY
   Aux_Error( ERROR_INFO, "GRAVITY must be enabled !!\n" );
#  endif

#  ifndef COMOVING
   Aux_Error( ERROR_INFO, "COMOVING must be enabled !!\n" );
#  endif

#  ifdef PARTICLE
   Aux_Error( ERROR_INFO, "PARTICLE must be disabled !!\n" );
#  endif


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Validating test problem %d ... done\n", TESTPROB_ID );

} // FUNCTION : Validate



#if ( MODEL == ELBDM  &&  defined GRAVITY  &&  defined COMOVING )
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
// ReadPara->Add( "KEY_IN_THE_FILE",      &VARIABLE,              DEFAULT,       MIN,              MAX               );
// ********************************************************************************************************************************
   ReadPara->Add( "SelSimHalo_Filename",   SelSimHalo_Filename,   NoDef_str,     Useless_str,      Useless_str       );

   ReadPara->Read( FileName );

   delete ReadPara;

// (1-2) set the default values

// (1-3) check the runtime parameters


// (2) load the radial profiles
   if ( OPT__INIT != INIT_BY_RESTART )
   {
//    load table (note that the 3rd column in the table is velocity and we need to convert it to phase later)
      const bool RowMajor_No  = false;       // load data into the column-major order
      const bool AllocMem_Yes = true;        // allocate memory for SelSimHalo_Prof
      const int  NCol         = 3;           // total number of columns to be loaded
      const int  Col[NCol]    = {0, 1, 2};   // target columns
      double *Prof_Radius, *Prof_Dens, *Prof_Phase;

      SelSimHalo_NBin = Aux_LoadTable( SelSimHalo_Prof, SelSimHalo_Filename, NCol, Col, RowMajor_No, AllocMem_Yes );

      Prof_Radius = SelSimHalo_Prof + 0*SelSimHalo_NBin;
      Prof_Dens   = SelSimHalo_Prof + 1*SelSimHalo_NBin;
      Prof_Phase  = SelSimHalo_Prof + 2*SelSimHalo_NBin;


//    backup the velocity data
      double *Prof_Velocity = new double [SelSimHalo_NBin];

      memcpy( Prof_Velocity, Prof_Phase, SelSimHalo_NBin*sizeof(double) );


//    integrate velocity to get phase
      Prof_Phase[0] = 0.0;

      for (int b=1; b<SelSimHalo_NBin; b++)
      {
         const int bm1 = b - 1;

         Prof_Phase[b] = Prof_Phase[bm1] + 0.5*( Prof_Velocity[b] + Prof_Velocity[bm1] )*( Prof_Radius[b] - Prof_Radius[bm1] );
      }

      delete [] Prof_Velocity;


//    convert to code units (assuming the input units are dimensionless)
      const double H0    = 100.0*Const_km/Const_s/Const_Mpc*HUBBLE0*UNIT_T;      // Hubble parameter at z=0 in the code units
      const double Eta2x = pow( 1.5*H0*ELBDM_ETA, -0.5 )*pow( Time[0], -0.25 );  // conversion factor between eta and the comoving distance

      for (int b=0; b<SelSimHalo_NBin; b++)
      {
         Prof_Radius[b] *= Eta2x;   // length_in_gamer  = comoving_distance
//       Prof_Dens  [b] *= 1.0;     // density_in_gamer = density_in_self_similar_solution
//       Prof_Phase [b] *= 1.0;     // phase_in_gamer   = phase_in_self_similar_solution
      }
   } // if ( OPT__INIT != INIT_BY_RESTART )


// (3) reset other general-purpose parameters
//     --> a helper macro PRINT_RESET_PARA is defined in Macro.h
   const long   End_Step_Default = __INT_MAX__;
   const double End_T_Default    = 1.0/3.0;     // z = 2.0

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
      Aux_Message( stdout, "  test problem ID = %d\n", TESTPROB_ID         );
      Aux_Message( stdout, "  profile table   = %s\n", SelSimHalo_Filename );
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

   const double *Prof_Radius = SelSimHalo_Prof + 0*SelSimHalo_NBin;
   const double *Prof_Dens   = SelSimHalo_Prof + 1*SelSimHalo_NBin;
   const double *Prof_Phase  = SelSimHalo_Prof + 2*SelSimHalo_NBin;

   double r, Dens, Phase;

   r     = sqrt( SQR(x-amr->BoxCenter[0]) + SQR(y-amr->BoxCenter[1]) + SQR(z-amr->BoxCenter[2]) );
   Dens  = Mis_InterpolateFromTable( SelSimHalo_NBin, Prof_Radius, Prof_Dens,  r );
   Phase = Mis_InterpolateFromTable( SelSimHalo_NBin, Prof_Radius, Prof_Phase, r );

   if ( Dens == NULL_REAL  ||  Phase == NULL_REAL )
      Aux_Error( ERROR_INFO, "interpolation failed at radius %13.7e (probably outside the input table) !!\n", r );

   fluid[DENS] = Dens;

#  if ( ELBDM_SCHEME == ELBDM_HYBRID )
   if ( amr->use_wave_flag[lv] ) {
#  endif
   fluid[REAL] = sqrt( Dens )*cos( Phase );
   fluid[IMAG] = sqrt( Dens )*sin( Phase );
#  if ( ELBDM_SCHEME == ELBDM_HYBRID )
   } else {
   fluid[PHAS] = Phase;
   fluid[STUB] = 0.0;
   }
#  endif

} // FUNCTION : SetGridIC



//-------------------------------------------------------------------------------------------------------
// Function    :  End_SelfSimilarHalo
// Description :  Free memory before terminating the program
//
// Note        :  1. Linked to the function pointer "End_User_Ptr" to replace "End_User()"
//
// Parameter   :  None
//-------------------------------------------------------------------------------------------------------
void End_SelfSimilarHalo()
{

   delete [] SelSimHalo_Prof;

} // FUNCTION : End_SelfSimilarHalo
#endif // #if ( MODEL == ELBDM  &&  defined GRAVITY  &&  defined COMOVING )



//-------------------------------------------------------------------------------------------------------
// Function    :  Init_TestProb_ELBDM_SelfSimilarHalo
// Description :  Test problem initializer
//
// Note        :  None
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Init_TestProb_ELBDM_SelfSimilarHalo()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );


// validate the compilation flags and runtime parameters
   Validate();


#  if ( MODEL == ELBDM  &&  defined GRAVITY  &&  defined COMOVING )
// set the problem-specific runtime parameters
   SetParameter();


   Init_Function_User_Ptr = SetGridIC;
   End_User_Ptr           = End_SelfSimilarHalo;
#  endif // if ( MODEL == ELBDM  &&  defined GRAVITY  &&  defined COMOVING )


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : Init_TestProb_ELBDM_SelfSimilarHalo
