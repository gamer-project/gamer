#include <random>
#include "GAMER.h"
#include "TestProb.h"

#if ( MODEL == HYDRO )

// problem-specific global variables
// =======================================================================================

// background parameters
static double   ICM_Density;             // ICM density
static double   ICM_Position;            // position of interface
static double   ICM_Angle;               // inclination of interface
static double   ICM_Tangent;             // tangent of interface angle
static double   Amb_Pressure;            // ambient pressure
static double   Lobe_ICM_Ratio;          // ratio of lobe/ICM densities
static double   Lobe_Density;            // lobe density

// jet parameters
static bool     Jet_Fire;                // [true/false]: jet on/off
static double   Jet_Velocity;            // jet y-velocity (units of c)
static double   Jet_Radius;              // radius of jet
static double   Jet_Position;            // position of jet
static double   Jet_Lobe_Ratio;          // ratio of jet/lobe densities
static double   Jet_Center[2];           // jet central coordinates
static double   Jet_Density;             // jet density

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
#  ifdef PARTICLE
   Aux_Error( ERROR_INFO, "PARTICLE must be disabled !!\n" );
#  endif

#  ifdef GRAVITY
   Aux_Error( ERROR_INFO, "GRAVITY must be disabled !!\n" );
#  endif

   if ( !OPT__UNIT )
       Aux_Error( ERROR_INFO, "OPT__UNIT must be enabled !!\n" );

// warnings
   if ( MPI_Rank == 0 )
   {
      for (int s=0; s<6; s++)
         if ( OPT__BC_FLU[s] != BC_FLU_OUTFLOW )
            Aux_Message( stderr, "WARNING : it's recommended to use the outflow BC (currently OPT__BC_FLU[%d] = %d != 2)\n",
                         s, OPT__BC_FLU[s] );
   }


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Validating test problem %d ... done\n", TESTPROB_ID );

} // FUNCTION : Validate



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
// ************************************************************************************************************************
// ReadPara->Add( "KEY_IN_THE_FILE",  &VARIABLE,        DEFAULT,       MIN,            MAX    );
// ************************************************************************************************************************

// background parameters
   ReadPara->Add( "ICM_Density",      &ICM_Density,     NoDef_double,  NoMin_double,   NoMax_double );
   ReadPara->Add( "ICM_Position",     &ICM_Position,    NoDef_double,  NoMin_double,   NoMax_double );
   ReadPara->Add( "ICM_Angle",        &ICM_Angle,       NoDef_double,  NoMin_double,   NoMax_double );
   ReadPara->Add( "Amb_Pressure",     &Amb_Pressure,    NoDef_double,  NoMin_double,   NoMax_double );
   ReadPara->Add( "Lobe_ICM_Ratio",   &Lobe_ICM_Ratio,  NoDef_double,  NoMin_double,   NoMax_double );

// jet parameters
   ReadPara->Add( "Jet_Fire",         &Jet_Fire,        false,         Useless_bool,   Useless_bool );
   ReadPara->Add( "Jet_Lobe_Ratio",   &Jet_Lobe_Ratio,  NoDef_double,  NoMin_double,   NoMax_double );
   ReadPara->Add( "Jet_Radius",       &Jet_Radius,      NoDef_double,  NoMin_double,   NoMax_double );
   ReadPara->Add( "Jet_Position",     &Jet_Position,    NoDef_double,  NoMin_double,   NoMax_double );
   ReadPara->Add( "Jet_Velocity",     &Jet_Velocity,    NoDef_double,  NoMin_double,   NoMax_double );

   ReadPara->Read( FileName );

   delete ReadPara;

// (1-2) convert to code unit
   Jet_Velocity *= Const_c   / UNIT_V;
   ICM_Density  *= 1.0       / UNIT_D;
   Jet_Radius   *= Const_kpc / UNIT_L;
   Jet_Position *= Const_kpc / UNIT_L;
   ICM_Position *= Const_kpc / UNIT_L;
   Amb_Pressure *= 1.0       / UNIT_P;

// (2) set the problem-specific derived parameters

   ICM_Tangent   = tan(ICM_Angle*M_PI/180.0);
   Lobe_Density  = Lobe_ICM_Ratio*ICM_Density;
   Jet_Density   = Jet_Lobe_Ratio*Lobe_Density;
   Jet_Center[0] = Jet_Position;
   Jet_Center[1] = amr->BoxCenter[2];

// (3) reset other general-purpose parameters
//     --> a helper macro PRINT_WARNING is defined in TestProb.h
   const long   End_Step_Default = __INT_MAX__;
   const double End_T_Default    = 100.0*Const_kyr / UNIT_T;

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
     Aux_Message( stdout, "  test problem ID          = %d\n",                TESTPROB_ID                   );
     Aux_Message( stdout, "  ICM_Density              = %14.7e g/cm^3\n",     ICM_Density*UNIT_D            );
     Aux_Message( stdout, "  ICM_Position             = %14.7e kpc\n",        ICM_Position*UNIT_L/Const_kpc );
     Aux_Message( stdout, "  ICM_Angle                = %14.7e degree\n",     ICM_Angle                     );
     Aux_Message( stdout, "  Amb_Pressure             = %14.7e erg/cm^3\n",   Amb_Pressure*UNIT_P           );
     Aux_Message( stdout, "  Lobe_ICM_Ratio           = %14.7e\n",            Lobe_ICM_Ratio                );
     Aux_Message( stdout, "  Lobe_Density             = %14.7e g/cm^3\n",     Lobe_Density*UNIT_D           );
   }

   if ( Jet_Fire && MPI_Rank == 0 )
   {
     Aux_Message( stdout, "  Jet_Radius               = %14.7e kpc\n",        Jet_Radius*UNIT_L/Const_kpc   );
     Aux_Message( stdout, "  Jet_Position             = %14.7e kpc\n",        Jet_Position*UNIT_L/Const_kpc );
     Aux_Message( stdout, "  Jet_Velocity             = %14.7e c\n",          Jet_Velocity*UNIT_V/Const_c   );
     Aux_Message( stdout, "  Jet_Density              = %14.7e g/cm^3\n",     Jet_Density*UNIT_D            );
     Aux_Message( stdout, "  Jet_Lobe_Ratio           = %14.7e\n",            Jet_Lobe_Ratio                );
   }

   if ( MPI_Rank == 0 )
   {
     Aux_Message( stdout, "=============================================================================\n"                   );
   }

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Setting runtime parameters ... done\n"                                    );

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
// variables for jet
   real PriReal[NCOMP_FLUID];

   if ( x < ICM_Position + ICM_Tangent*(y-0.5*amr->BoxSize[1]) ) {
       PriReal[0] = (real)ICM_Density;
   } else {
       PriReal[0] = (real)Lobe_Density;
   }
   PriReal[1] = 0.0;
   PriReal[2] = 0.0;
   PriReal[3] = 0.0;
   PriReal[4] = (real)Amb_Pressure;

   Hydro_Pri2Con( PriReal, fluid, NULL_BOOL, NULL_INT, NULL, EoS_DensPres2Eint_CPUPtr,
                  EoS_Temp2HTilde_CPUPtr, EoS_HTilde2Temp_CPUPtr,
                  EoS_AuxArray_Flt, EoS_AuxArray_Int, h_EoS_Table, NULL );


} // FUNCTION : SetGridIC


//-------------------------------------------------------------------------------------------------------
// Function    :  BC
// Description :  Set the extenral boundary condition to the analytical solution
//
// Note        :  1. Linked to the function pointer "BC_User_Ptr"
//
// Parameter   :  fluid    : Fluid field to be set
//                x/y/z    : Physical coordinates
//                Time     : Physical time
//                lv       : Refinement level
//                AuxArray : Auxiliary array
//
// Return      :  fluid
//-------------------------------------------------------------------------------------------------------
void BC( real Array[], const int ArraySize[], real fluid[], const int NVar_Flu,
	 const int GhostSize, const int idx[], const double pos[], const double Time,
	 const int lv, const int TFluVarIdxList[], double AuxArray[] )
{

    if ( !Jet_Fire ) return;

    real PriReal[NCOMP_FLUID];

    const int j_ref = GhostSize;  // reference j index

    double rad = sqrt( SQR(pos[0]-Jet_Center[0]) + SQR(pos[2]-Jet_Center[1]) );

    // 1D array -> 3D array
    real (*Array3D)[ArraySize[0]][ArraySize[1]][ArraySize[2]] = ( real (*)[ArraySize[0]][ArraySize[1]][ArraySize[2]] )Array;

    if ( rad <= Jet_Radius )
    {

        // set fluid variable inside source

        PriReal[0] = (real)Jet_Density;
        PriReal[1] = 0.0;
        PriReal[2] = (real)Jet_Velocity;
        PriReal[3] = 0.0;
        PriReal[4] = (real)Amb_Pressure;

    } else {

        for (int v=0; v<NVar_Flu; v++)
            PriReal[ TFluVarIdxList[v] ] = Array3D[v][idx[2]][j_ref][idx[0]];

    }

    Hydro_Pri2Con( PriReal, fluid, NULL_BOOL, NULL_INT, NULL, EoS_DensPres2Eint_CPUPtr,
                   EoS_Temp2HTilde_CPUPtr, EoS_HTilde2Temp_CPUPtr,
                   EoS_AuxArray_Flt, EoS_AuxArray_Int, h_EoS_Table, NULL );

} // FUNCTION : BC


//-------------------------------------------------------------------------------------------------------
// Function    :  Init_TestProb_SRHydro_JetICMWall
// Description :  Test problem initializer
//
// Note        :  None
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Init_TestProb_SRHydro_JetICMWall()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );


// validate the compilation flags and runtime parameters
   Validate();


// set the problem-specific runtime parameters
   SetParameter();


// get enclosed mass
   Init_Function_User_Ptr   = SetGridIC;
   Flag_User_Ptr            = NULL;
   Flag_Region_Ptr          = NULL;
   Mis_GetTimeStep_User_Ptr = NULL;
   BC_User_Ptr              = BC;
   Flu_ResetByUser_Func_Ptr = NULL;
   Output_User_Ptr          = NULL;
   Aux_Record_User_Ptr      = NULL;
   End_User_Ptr             = NULL;


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : Init_TestProb_SRHydro_JetICMWall

#endif
