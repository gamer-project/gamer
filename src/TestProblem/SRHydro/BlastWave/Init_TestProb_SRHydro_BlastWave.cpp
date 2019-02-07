#include "GAMER.h"
#include "TestProb.h"


// problem-specific global variables
// =======================================================================================
static double Blast_Dens_Bg;       // background mass density
static double Blast_Pres_Bg;       // background pressure
static double Blast_Dens_Ratio;    // density ratio of center to background
static double Blast_Pres_Ratio;    // pressure ratio of center to background
static double Blast_Radius;        // explosion radius
static double Blast_Center[3];     // explosion center
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


#  if ( MODEL != SR_HYDRO )
   Aux_Error( ERROR_INFO, "MODEL != SR_HYDRO !!\n" );
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



#if ( MODEL == SR_HYDRO )
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

// add parameters in the following format:
// --> note that VARIABLE, DEFAULT, MIN, and MAX must have the same data type
// --> some handy constants (e.g., NoMin_int, Eps_float, ...) are defined in "include/ReadPara.h"
// ********************************************************************************************************************************
// ReadPara->Add( "KEY_IN_THE_FILE",   &VARIABLE_ADDRESS,      DEFAULT,       MIN,              MAX               );
// ********************************************************************************************************************************
   ReadPara->Add( "Blast_Dens_Bg",     &Blast_Dens_Bg,         -1.0,          Eps_double,       NoMax_double      );
   ReadPara->Add( "Blast_Pres_Bg",     &Blast_Pres_Bg,         -1.0,          Eps_double,       NoMax_double      );
   ReadPara->Add( "Blast_Dens_Ratio",  &Blast_Dens_Ratio,      -1.0,          Eps_double,       NoMax_double      );
   ReadPara->Add( "Blast_Pres_Ratio",  &Blast_Pres_Ratio,      -1.0,          Eps_double,       NoMax_double      );
   ReadPara->Add( "Blast_Radius",      &Blast_Radius,          -1.0,          Eps_double,       NoMax_double      );
   ReadPara->Add( "Blast_Center_X",    &Blast_Center[0],       -1.0,          NoMin_double,     amr->BoxSize[0]   );
   ReadPara->Add( "Blast_Center_Y",    &Blast_Center[1],       -1.0,          NoMin_double,     amr->BoxSize[1]   );
   ReadPara->Add( "Blast_Center_Z",    &Blast_Center[2],       -1.0,          NoMin_double,     amr->BoxSize[2]   );

   ReadPara->Read( FileName );

   delete ReadPara;

// set the default explosion center
   for (int d=0; d<3; d++)
      if ( Blast_Center[d] < 0.0 )  Blast_Center[d] = 0.5*amr->BoxSize[d];


// (2) set the problem-specific derived parameters


// (3) reset other general-purpose parameters
//     --> a helper macro PRINT_WARNING is defined in TestProb.h
   const double End_T_Default    = 5.0e-3;
   const long   End_Step_Default = __INT_MAX__;

   if ( END_STEP < 0 ) {
      END_STEP = End_Step_Default;
      PRINT_WARNING( "END_STEP", END_STEP, FORMAT_LONG );
   }

   if ( END_T < 0.0 ) {
      END_T = End_T_Default;
      PRINT_WARNING( "END_T", END_T, FORMAT_REAL );
   }


// (4) make a note
   if ( MPI_Rank == 0 )
   {
      Aux_Message( stdout, "=============================================================================\n" );
      Aux_Message( stdout, "  test problem ID                         = %d\n",     TESTPROB_ID );
      Aux_Message( stdout, "  background mass density                 = %13.7e\n", Blast_Dens_Bg );
      Aux_Message( stdout, "  background pressure                     = %13.7e\n", Blast_Pres_Bg);
      Aux_Message( stdout, "  density ratio of center to background   = %13.7e\n", Blast_Dens_Ratio );
      Aux_Message( stdout, "  pressure ratio of center to background  = %13.7e\n", Blast_Pres_Ratio );
      Aux_Message( stdout, "  explosion radius                        = %13.7e\n", Blast_Radius );
      Aux_Message( stdout, "  explosion center                        = (%13.7e, %13.7e, %13.7e)\n",
                                                          Blast_Center[0], Blast_Center[1], Blast_Center[2] );
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

   double Prim_BG1[5] = { Blast_Dens_Bg, 0, 0, 0, Blast_Pres_Bg }; // store 3-velocity
   double Prim_BG2[5] = {0};                                       // store 4-velocity
   double Cons_BG[5] = {0};

   double Prim_EXP1[5] = { Blast_Dens_Bg*Blast_Dens_Ratio, 0, 0, 0, Blast_Pres_Bg*Blast_Pres_Ratio }; // store 3-velocity
   double Prim_EXP2[5] = {0};                                                                         // store 4-velocity
   double Cons_EXP[5] = {0};

   SRHydro_3Velto4Vel ( Prim_BG1, Prim_BG2);
   SRHydro_3Velto4Vel ( Prim_EXP1, Prim_EXP2);

   SRHydro_Pri2Con(Prim_BG2, Cons_BG, GAMMA);
   SRHydro_Pri2Con(Prim_EXP2, Cons_EXP, GAMMA);

   if ( r <= Blast_Radius )
   {
     fluid[DENS] = Cons_EXP[0];
     fluid[MOMX] = Cons_EXP[1];
     fluid[MOMY] = Cons_EXP[2];
     fluid[MOMZ] = Cons_EXP[3];
     fluid[ENGY] = Cons_EXP[4];
   }
   else
   { 
     fluid[DENS] = Cons_BG[0];
     fluid[MOMX] = Cons_BG[1];
     fluid[MOMY] = Cons_BG[2];
     fluid[MOMZ] = Cons_BG[3];
     fluid[ENGY] = Cons_BG[4];
   }

} // FUNCTION : SetGridIC
#endif // #if ( MODEL == SR_HYDRO )



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
void Init_TestProb_SRHydro_BlastWave()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );


// validate the compilation flags and runtime parameters
   Validate();


#  if ( MODEL == SR_HYDRO )
// set the problem-specific runtime parameters
   SetParameter();


// set the function pointers of various problem-specific routines
   Init_Function_User_Ptr   = SetGridIC;
   Output_User_Ptr          = NULL;
   Flag_User_Ptr            = NULL;
   Mis_GetTimeStep_User_Ptr = NULL;
   Aux_Record_User_Ptr      = NULL;
   BC_User_Ptr              = NULL;
   Flu_ResetByUser_Func_Ptr = NULL;
   End_User_Ptr             = NULL;
#  endif // #if ( MODEL == SR_HYDRO )


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : Init_TestProb_Hydro_BlastWave

