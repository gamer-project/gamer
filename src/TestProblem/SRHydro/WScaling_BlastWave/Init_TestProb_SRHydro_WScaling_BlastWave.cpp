#include "GAMER.h"
#include "TestProb.h"

#if  ( MODEL == SR_HYDRO )

// problem-specific global variables
// =======================================================================================
static double Blast_Dens_Bg;         // background mass density
static double Blast_Pres_Bg;         // background pressure
static double Blast_Dens_Ratio;      // density ratio of center to background
static double Blast_Pres_Ratio;      // pressure ratio of center to background
static double Blast_Radius;          // initial explosion radius
static double (*Blast_Center)[3];    // explosion center
static int (*Blast_Center_Temp)[3];  // explosion center
static int    Number_BlastWave_X;    // number of blast wave in x-direction
static int    Number_BlastWave_Y;    // number of blast wave in y-direction
static int    Number_BlastWave_Z;    // number of blast wave in z-direction
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
// ReadPara->Add( "KEY_IN_THE_FILE",     &VARIABLE_ADDRESS,      DEFAULT,       MIN,              MAX               );
// ********************************************************************************************************************************
   ReadPara->Add( "Blast_Dens_Bg",       &Blast_Dens_Bg,         -1.0,          Eps_double,       NoMax_double      );
   ReadPara->Add( "Blast_Pres_Bg",       &Blast_Pres_Bg,         -1.0,          Eps_double,       NoMax_double      );
   ReadPara->Add( "Blast_Dens_Ratio",    &Blast_Dens_Ratio,      -1.0,          Eps_double,       NoMax_double      );
   ReadPara->Add( "Blast_Pres_Ratio",    &Blast_Pres_Ratio,      -1.0,          Eps_double,       NoMax_double      );
   ReadPara->Add( "Blast_Radius",        &Blast_Radius,          -1.0,          Eps_double,       NoMax_double      );
   ReadPara->Add( "Number_BlastWave_X",  &Number_BlastWave_X,       1,                   1,       NoMax_int         );
   ReadPara->Add( "Number_BlastWave_Y",  &Number_BlastWave_Y,       1,                   1,       NoMax_int         );
   ReadPara->Add( "Number_BlastWave_Z",  &Number_BlastWave_Z,       1,                   1,       NoMax_int         );

   ReadPara->Read( FileName );

   delete ReadPara;

   int Total_BlastWave = Number_BlastWave_X * Number_BlastWave_Y * Number_BlastWave_Z;

   Blast_Center      = new double [Total_BlastWave][3];
   Blast_Center_Temp = new    int [Total_BlastWave][3];

   double dX[3];

   dX[0] = 0.5 * amr->BoxSize[0] / (double) Number_BlastWave_X;
   dX[1] = 0.5 * amr->BoxSize[1] / (double) Number_BlastWave_Y;
   dX[2] = 0.5 * amr->BoxSize[2] / (double) Number_BlastWave_Z;

// set the explosion centers
   for (int i=0;i<Total_BlastWave;i++)
    {
       Blast_Center_Temp[i][2] = i % Number_BlastWave_Z;
       Blast_Center_Temp[i][1] = ( ( i-Blast_Center_Temp[i][2] ) / Number_BlastWave_Z ) % Number_BlastWave_Y;
       Blast_Center_Temp[i][0] = ( i-Blast_Center_Temp[i][2]-Blast_Center_Temp[i][1] ) / (Number_BlastWave_Y*Number_BlastWave_Z);

       Blast_Center[i][2] = ( 2 * Blast_Center_Temp[i][2] + 1 ) * dX[2];
       Blast_Center[i][1] = ( 2 * Blast_Center_Temp[i][1] + 1 ) * dX[1];
       Blast_Center[i][0] = ( 2 * Blast_Center_Temp[i][0] + 1 ) * dX[0];
    }
// (2) reset other general-purpose parameters
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


// (3) make a note
   if ( MPI_Rank == 0 )
   {
      Aux_Message( stdout, "=============================================================================\n" );
      Aux_Message( stdout, "  test problem ID                         = %d\n",     TESTPROB_ID );
      Aux_Message( stdout, "  background mass density                 = %13.7e\n", Blast_Dens_Bg );
      Aux_Message( stdout, "  background pressure                     = %13.7e\n", Blast_Pres_Bg);
      Aux_Message( stdout, "  density ratio of center to background   = %13.7e\n", Blast_Dens_Ratio );
      Aux_Message( stdout, "  pressure ratio of center to background  = %13.7e\n", Blast_Pres_Ratio );
      Aux_Message( stdout, "  number of blast wave in x-direction     = %13d\n"  , Number_BlastWave_X );
      Aux_Message( stdout, "  number of blast wave in y-direction     = %13d\n"  , Number_BlastWave_Y );
      Aux_Message( stdout, "  number of blast wave in z-direction     = %13d\n"  , Number_BlastWave_Z );
      Aux_Message( stdout, "  explosion radius                        = %13.7e\n", Blast_Radius );
      for (int i=0; i<Total_BlastWave; i++)
        for (int d=0; d<3; d++)
            Aux_Message( stdout, "  blast wave center[%d][%d]                 = %13.7e\n", i, d,Blast_Center[i][d] );
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
   double Prim_BG[5] = { Blast_Dens_Bg, 0, 0, 0, Blast_Pres_Bg };
   double Cons_BG[5] = { 0 };

   double Prim_EXP[5] = { Blast_Dens_Bg*Blast_Dens_Ratio, 0, 0, 0, Blast_Pres_Bg*Blast_Pres_Ratio };
   double Cons_EXP [5] = { 0 };

   int Total = Number_BlastWave_X * Number_BlastWave_Y * Number_BlastWave_Z;

   double r[Total];

   SRHydro_Pri2Con ( Prim_BG , Cons_BG,  GAMMA);
   SRHydro_Pri2Con ( Prim_EXP, Cons_EXP, GAMMA);

   for (int i=0; i<Total; i++)
    {
      r[i] = SQRT( SQR(x-Blast_Center[i][0]) + SQR(y-Blast_Center[i][1]) + SQR(z-Blast_Center[i][2]) );
    
       if ( r[i] <= Blast_Radius )
       {
         fluid[DENS] = (real) Cons_EXP[0];
         fluid[MOMX] = (real) Cons_EXP[1];
         fluid[MOMY] = (real) Cons_EXP[2];
         fluid[MOMZ] = (real) Cons_EXP[3];
         fluid[ENGY] = (real) Cons_EXP[4];
         i = Total; // break loop immediately
       }
       else
       { 
         fluid[DENS] = (real) Cons_BG[0];
         fluid[MOMX] = (real) Cons_BG[1];
         fluid[MOMY] = (real) Cons_BG[2];
         fluid[MOMZ] = (real) Cons_BG[3];
         fluid[ENGY] = (real) Cons_BG[4];
       }
    }
} // FUNCTION : SetGridIC
#endif // #if ( MODEL == SR_HYDRO )

void End_WSBlastWave()
{
   delete [] Blast_Center;
   delete [] Blast_Center_Temp;
}

//-------------------------------------------------------------------------------------------------------
// Function    :  Init_TestProb_SRHydro_BlastWave
// Description :  Test problem initializer
//
// Note        :  None
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Init_TestProb_SRHydro_WScaling_BlastWave()
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
   End_User_Ptr             = End_WSBlastWave;
#  endif // #if ( MODEL == SR_HYDRO )


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : Init_TestProb_SRHydro_BlastWave

#endif
