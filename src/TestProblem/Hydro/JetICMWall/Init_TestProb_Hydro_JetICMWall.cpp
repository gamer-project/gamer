#include "GAMER.h"



static FieldIdx_t JetFieldIdx  = Idx_Undefined;
static FieldIdx_t ICMFieldIdx  = Idx_Undefined;
static FieldIdx_t LobeFieldIdx = Idx_Undefined;
static FieldIdx_t IntFieldIdx  = Idx_Undefined;

// problem-specific global variables
// =======================================================================================
// background parameters
static double   ICM_Density;             // ICM density
static double   Jump_Position_x;         // position of interface
static double   Jump_Position_y;         // position of interface
static double   Jump_Angle;              // inclination of interface
static double   Jump_Tangent;            // tangent of interface angle
static double   Jump_Width;              // width of jump
static double   Amb_Pressure;            // ambient pressure
static double   Lobe_ICM_Ratio;          // ratio of lobe/ICM densities
static double   Lobe_Density;            // lobe density

// jet parameters
static bool     Jet_Fire;                // [true/false]: jet on/off
static double   Jet_Velocity;            // jet y-velocity (units of c)
static double   Jet_VelSlope;            // Slope of velocity gradient across jet
static double   Jet_VelCenter;           // jet central velocity
static double   Jet_Radius;              // radius of jet
static double   Jet_Position;            // position of jet
static double   Jet_Lobe_Ratio;          // ratio of jet/lobe densities
static double   Jet_Center[2];           // jet central coordinates
static double   Jet_Density;             // jet density
static double   Jet_Gamma;               // jet relativistic gamma
static double   Jet_PrecessAngle;        // jet relativistic gamma
static double   Jet_PrecessPeriod;       // jet relativistic gamma
static double   Jet_PrecessOmega;        // jet relativistic gamma
static double   Jet_Cosine;              // jet relativistic gamma
static double   Jet_Sine;                // jet relativistic gamma
static double   Jump_Sine;
static double   Jump_Cosine;
// =======================================================================================

static void JetBC( real Array[], const int ArraySize[], real fluid[], const int NVar_Flu,
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
#  ifndef SRHD
   Aux_Error( ERROR_INFO, "SRHD must be enabled !!\n" );
#  endif

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



#if ( MODEL == HYDRO )
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
// ReadPara->Add( "KEY_IN_THE_FILE",    &VARIABLE,          DEFAULT,       MIN,            MAX          );
// ************************************************************************************************************************

// background parameters
   ReadPara->Add( "ICM_Density",        &ICM_Density,       NoDef_double,  NoMin_double,   NoMax_double );
   ReadPara->Add( "Jump_Position_x",    &Jump_Position_x,   NoDef_double,  NoMin_double,   NoMax_double );
   ReadPara->Add( "Jump_Position_y",    &Jump_Position_y,   NoDef_double,  NoMin_double,   NoMax_double );
   ReadPara->Add( "Jump_Angle",         &Jump_Angle,        NoDef_double,  NoMin_double,   NoMax_double );
   ReadPara->Add( "Amb_Pressure",       &Amb_Pressure,      NoDef_double,  NoMin_double,   NoMax_double );
   ReadPara->Add( "Lobe_ICM_Ratio",     &Lobe_ICM_Ratio,    NoDef_double,  NoMin_double,   NoMax_double );
   ReadPara->Add( "Jump_Width",         &Jump_Width,        NoDef_double,  NoMin_double,   NoMax_double );

// jet parameters
   ReadPara->Add( "Jet_Fire",           &Jet_Fire,          false,         Useless_bool,   Useless_bool );
   ReadPara->Add( "Jet_Lobe_Ratio",     &Jet_Lobe_Ratio,    NoDef_double,  NoMin_double,   NoMax_double );
   ReadPara->Add( "Jet_Radius",         &Jet_Radius,        NoDef_double,  NoMin_double,   NoMax_double );
   ReadPara->Add( "Jet_Position",       &Jet_Position,      NoDef_double,  NoMin_double,   NoMax_double );
   ReadPara->Add( "Jet_Velocity",       &Jet_Velocity,      NoDef_double,  NoMin_double,   NoMax_double );
   ReadPara->Add( "Jet_VelSlope",       &Jet_VelSlope,      0.0,           NoMin_double,   0.0          );
   ReadPara->Add( "Jet_VelCenter",      &Jet_VelCenter,     NoDef_double,  0.0,            NoMax_double );
   ReadPara->Add( "Jet_PrecessAngle",   &Jet_PrecessAngle,  0.0,           0.0,            NoMax_double );
   ReadPara->Add( "Jet_PrecessPeriod",  &Jet_PrecessPeriod, 0.0,           0.0,            NoMax_double );

   ReadPara->Read( FileName );

   delete ReadPara;

   Jet_Gamma = sqrt(1.0+SQR(Jet_Velocity));

// (1-2) convert to code unit
   Jet_Velocity      *= Const_c   / UNIT_V;
   ICM_Density       *= 1.0       / UNIT_D;
   Jet_Radius        *= Const_kpc / UNIT_L;
   Jet_Position      *= Const_kpc / UNIT_L;
   Jump_Position_x   *= Const_kpc / UNIT_L;
   Jump_Position_y   *= Const_kpc / UNIT_L;
   Amb_Pressure      *= 1.0       / UNIT_P;
   Jump_Width        *= Const_kpc / UNIT_L;
   Jet_PrecessPeriod *= 1000.0*Const_yr / UNIT_T;


// (2) set the problem-specific derived parameters
   Jump_Tangent     = tan( Jump_Angle*M_PI/180.0 );
   Jump_Sine        = sin( Jump_Angle*M_PI/180.0 );
   Jump_Cosine      = cos( Jump_Angle*M_PI/180.0 );
   Lobe_Density     = Lobe_ICM_Ratio*ICM_Density;
   Jet_Density      = Jet_Lobe_Ratio*Lobe_Density;
   Jet_Center[0]    = Jet_Position;
   Jet_Center[1]    = amr->BoxCenter[2];
   Jet_PrecessOmega = 2.0*M_PI/Jet_PrecessPeriod;
   Jet_Cosine       = cos( Jet_PrecessAngle*M_PI/180.0 );
   Jet_Sine         = sin( Jet_PrecessAngle*M_PI/180.0 );


// (3) reset other general-purpose parameters
   const long   End_Step_Default = __INT_MAX__;
   const double End_T_Default    = 100.0*Const_kyr / UNIT_T;

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
     Aux_Message( stdout, "=============================================================================\n " );
     Aux_Message( stdout, "  test problem ID = %d\n",              TESTPROB_ID                      );
     Aux_Message( stdout, "  ICM_Density     = %14.7e g/cm^3\n",   ICM_Density*UNIT_D               );
     Aux_Message( stdout, "  Jump_Position_x = %14.7e kpc\n",      Jump_Position_x*UNIT_L/Const_kpc );
     Aux_Message( stdout, "  Jump_Position_y = %14.7e kpc\n",      Jump_Position_y*UNIT_L/Const_kpc );
     Aux_Message( stdout, "  Jump_Angle      = %14.7e degree\n",   Jump_Angle                       );
     Aux_Message( stdout, "  Jump_Width      = %14.7e kpc\n",      Jump_Width*UNIT_L/Const_kpc      );
     Aux_Message( stdout, "  Amb_Pressure    = %14.7e erg/cm^3\n", Amb_Pressure*UNIT_P              );
     Aux_Message( stdout, "  Lobe_ICM_Ratio  = %14.7e\n",          Lobe_ICM_Ratio                   );
     Aux_Message( stdout, "  Lobe_Density    = %14.7e g/cm^3\n",   Lobe_Density*UNIT_D              );
   }

   if ( Jet_Fire && MPI_Rank == 0 )
   {
     Aux_Message( stdout, "  Jet_Radius      = %14.7e kpc\n",      Jet_Radius*UNIT_L/Const_kpc      );
     Aux_Message( stdout, "  Jet_Position    = %14.7e kpc\n",      Jet_Position*UNIT_L/Const_kpc    );
     Aux_Message( stdout, "  Jet_Velocity    = %14.7e c\n",        Jet_Velocity*UNIT_V/Const_c      );
     Aux_Message( stdout, "  Jet_VelSlope    = %14.7e\n",          Jet_VelSlope                     );
     Aux_Message( stdout, "  Jet_Density     = %14.7e g/cm^3\n",   Jet_Density*UNIT_D               );
     Aux_Message( stdout, "  Jet_Lobe_Ratio  = %14.7e\n",          Jet_Lobe_Ratio                   );
   }

   if ( MPI_Rank == 0 )
   {
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

// variables for jet
   real PriReal[NCOMP_TOTAL];

   double xx = (x - Jump_Position_x)*Jump_Cosine - (y - Jump_Position_y)*Jump_Sine;
   double d;

   double xw = xx/Jump_Width;
   if (xw > 200.0)
     d = Lobe_Density;
   else
     d = (ICM_Density + Lobe_Density*exp(xw)) / (1.0 + exp(xw));

   PriReal[0] = (real)d;
   PriReal[1] = 0.0;
   PriReal[2] = 0.0;
   PriReal[3] = 0.0;
   PriReal[4] = (real)Amb_Pressure;

   Hydro_Pri2Con( PriReal, fluid, false, PassiveNorm_NVar, PassiveNorm_VarIdx,
                  EoS_DensPres2Eint_CPUPtr, EoS_Temp2HTilde_CPUPtr, EoS_HTilde2Temp_CPUPtr,
                  EoS_AuxArray_Flt, EoS_AuxArray_Int, h_EoS_Table, NULL );

   double ICM_x = (d-0.4*ICM_Density)/(0.4*ICM_Density);
   double lobe_x = (2.5*Lobe_Density-d)/(1.25*Lobe_Density);
   if ( ICM_x  < 0.0 ) ICM_x  = 0.0;
   if ( lobe_x < 0.0 ) lobe_x = 0.0;
   if ( ICM_x  > 1.0 ) ICM_x  = 1.0;
   if ( lobe_x > 1.0 ) lobe_x = 1.0;

   fluid[JetFieldIdx ] = 0.0;
   fluid[ICMFieldIdx ] = (real)(ICM_x*d);
   fluid[LobeFieldIdx] = (real)(lobe_x*d);
   fluid[IntFieldIdx ] = (real)((1-ICM_x-lobe_x)*d);

} // FUNCTION : SetGridIC



//-------------------------------------------------------------------------------------------------------
// Function    :  JetBC
// Description :  Set the external boundary condition for the JetICMWall problem. On the -y
//                boundary, inside the jet region, the jet inflows with the jet density and a
//                velocity profile that peaks in the jet center and linearly falls off to the
//                edge of the jet. The jet passive scalar density is set to the jet density
//                within this region. Outside the jet region, the boundary condition is "diode",
//                meaning outflow-only.
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
// Return                        : fluid
//-------------------------------------------------------------------------------------------------------
void JetBC( real Array[], const int ArraySize[], real BVal[], const int NVar_Flu,
            const int GhostSize, const int idx[], const double pos[], const double Time,
            const int lv, const int TFluVarIdxList[], double AuxArray[] )
{

   int i, j, k;

   i = idx[0];
   j = idx[1];
   k = idx[2];

   real PriReal[NCOMP_TOTAL];

   const int j_ref = GhostSize;  // reference j index

   int TFluVarIdx;

   double rad = sqrt( SQR(pos[0]-Jet_Center[0]) + SQR(pos[2]-Jet_Center[1]) );

// 1D array -> 3D array
   typedef real (*vla)[ ArraySize[2] ][ ArraySize[1] ][ ArraySize[0] ];
   vla Array3D = ( vla )Array;

   double x = rad/Jet_Radius;

   if ( Jet_Fire  &&  x <= 1.0 )
   {
      double u_jet    = Jet_Velocity*( Jet_VelSlope*x+Jet_VelCenter );
      double cos_phi  = cos( Jet_PrecessOmega*Time );
      double sin_phi  = sin( Jet_PrecessOmega*Time );
      double LntzFact = sqrt( 1.0 + u_jet*u_jet );
      double u_jet_x  = 0.0;
      double u_jet_y  = u_jet;
      double u_jet_z  = 0.0;

//    set fluid variable inside source
      PriReal[0           ] = (real)Jet_Density;
      PriReal[1           ] = (real)0.0;
      PriReal[2           ] = (real)u_jet_y;
      PriReal[3           ] = (real)0.0;
      PriReal[4           ] = (real)Amb_Pressure;
      PriReal[JetFieldIdx ] = (real)Jet_Density;
      PriReal[ICMFieldIdx ] = (real)0.0;
      PriReal[LobeFieldIdx] = (real)0.0;
      PriReal[IntFieldIdx ] = (real)0.0;

      Hydro_Pri2Con( PriReal, BVal, false, PassiveNorm_NVar, PassiveNorm_VarIdx,
                     EoS_DensPres2Eint_CPUPtr, EoS_Temp2HTilde_CPUPtr, EoS_HTilde2Temp_CPUPtr,
                     EoS_AuxArray_Flt, EoS_AuxArray_Int, h_EoS_Table, NULL );
   } // if ( Jet_Fire  &&  x <= 1.0 )

   else
   {
      for (int v=0; v<NVar_Flu; v++) {
         TFluVarIdx = TFluVarIdxList[v];

         if ( TFluVarIdx == MOMY )
            BVal[TFluVarIdx] = MIN( Array3D[v][k][j_ref][i], 0.0 );
         else
            BVal[TFluVarIdx] = Array3D[v][k][j_ref][i];
     }
   } // if ( Jet_Fire  &&  x <= 1.0 ) ... else ...

} // FUNCTION : JetBC



//-------------------------------------------------------------------------------------------------------
// Function    :  AddNewField_JetICMWall
// Description :  Add the problem-specific fields
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
void AddNewField_JetICMWall()
{

  if ( JetFieldIdx == Idx_Undefined )
    JetFieldIdx = AddField( "JetField", FIXUP_FLUX_YES, FIXUP_REST_YES,
                            NORMALIZE_YES, INTERP_FRAC_NO );
  if ( ICMFieldIdx == Idx_Undefined )
    ICMFieldIdx = AddField( "ICMField", FIXUP_FLUX_YES, FIXUP_REST_YES,
                            NORMALIZE_YES, INTERP_FRAC_NO );
  if ( LobeFieldIdx == Idx_Undefined )
    LobeFieldIdx = AddField( "LobeField", FIXUP_FLUX_YES, FIXUP_REST_YES,
                             NORMALIZE_YES, INTERP_FRAC_NO );
  if ( IntFieldIdx == Idx_Undefined )
    IntFieldIdx = AddField( "IntField", FIXUP_FLUX_YES, FIXUP_REST_YES,
                            NORMALIZE_YES, INTERP_FRAC_NO );

} // FUNCTION : AddNewField_JetICMWall
#endif // #if ( MODEL == HYDRO )



//-------------------------------------------------------------------------------------------------------
// Function    :  Init_TestProb_Hydro_JetICMWall
// Description :  Test problem initializer
//
// Note        :  None
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Init_TestProb_Hydro_JetICMWall()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );


// validate the compilation flags and runtime parameters
   Validate();


#  if ( MODEL == HYDRO )
// set the problem-specific runtime parameters
   SetParameter();


// set the function pointers of various problem-specific routines
   Init_Function_User_Ptr   = SetGridIC;
   Init_Field_User_Ptr      = AddNewField_JetICMWall;
   Flag_User_Ptr            = NULL;
   Flag_Region_Ptr          = NULL;
   Mis_GetTimeStep_User_Ptr = NULL;
   BC_User_Ptr              = JetBC;
   Flu_ResetByUser_Func_Ptr = NULL;
   Output_User_Ptr          = NULL;
   Aux_Record_User_Ptr      = NULL;
   End_User_Ptr             = NULL;
#  endif // #if ( MODEL == HYDRO )


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : Init_TestProb_SRHydro_JetICMWall
