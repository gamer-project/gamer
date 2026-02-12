#include "GAMER.h"



// problem-specific global variables
// =======================================================================================
static double SNBlast_Dens_Bg;          // background mass density (in g/cm^3)
static double SNBlast_Temp_Bg;          // background temperature  (in K)
static double SNBlast_Pres_Bg;          // background pressure
static double SNBlast_Eint_Bg;          // backgroudn internal energy
static double SNBlast_MetalMassFrac_Bg; // background metal mass fraction (metal_mass / gas_mass)
       double SNBlast_ParMass;          // supernova progenitor particle mass (in Msun)
       double SNBlast_ParCenter[3];     // supernova progenitor initial center
       double SNBlast_ParVelocity[3];   // supernova progenitor initial velocity
       double SNBlast_ParMetalMassFrac; // supernova progenitor initial metal mass fraction (metal_mass / particle_mass)
// =======================================================================================

// problem-specific function prototypes
#ifdef MASSIVE_PARTICLES
void Par_Init_ByFunction_SNFeedbackBlastWave( const long NPar_ThisRank, const long NPar_AllRank,
                                              real_par *ParMass, real_par *ParPosX, real_par *ParPosY, real_par *ParPosZ,
                                              real_par *ParVelX, real_par *ParVelY, real_par *ParVelZ, real_par *ParTime,
                                              long_par *ParType, real_par *AllAttributeFlt[PAR_NATT_FLT_TOTAL],
                                              long_par *AllAttributeInt[PAR_NATT_INT_TOTAL] );
bool Flag_SNFeedbackBlastWave( const int i, const int j, const int k, const int lv, const int PID, const double *Threshold );
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

#  ifndef GRAVITY
   Aux_Error( ERROR_INFO, "GRAVITY must be enabled !!\n" );
#  endif

#  ifdef COMOVING
   Aux_Error( ERROR_INFO, "COMOVING must be disabled !!\n" );
#  endif

#  ifndef PARTICLE
   Aux_Error( ERROR_INFO, "PARTICLE must be enabled !!\n" );
#  endif

#  ifdef STAR_FORMATION
   Aux_Error( ERROR_INFO, "STAR_FORMATION must be disabled !!\n" );
#  endif

#  ifndef FEEDBACK
   Aux_Error( ERROR_INFO, "FEEDBACK must be enabled !!\n" );
#  endif

   if ( !OPT__UNIT )
      Aux_Error( ERROR_INFO, "OPT__UNIT must be enabled !!\n" );

#  ifdef GRAVITY
   if ( OPT__SELF_GRAVITY )
      Aux_Error( ERROR_INFO, "OPT__SELF_GRAVITY must be disabled !!\n" );
#  endif

#  ifdef PARTICLE
   if ( OPT__INIT == INIT_BY_FUNCTION  &&  amr->Par->Init != PAR_INIT_BY_FUNCTION )
      Aux_Error( ERROR_INFO, "please set PAR_INIT = 1 (by FUNCTION) !!\n" );

   if ( amr->Par->NPar_Active_AllRank != 1 )
      Aux_Error( ERROR_INFO, "PAR_NPAR must be 1 !!\n" );
#  endif

#  ifdef FEEDBACK
   if ( !FB_RESOLVED_SNEII )
      Aux_Error( ERROR_INFO, "FB_RESOLVED_SNEII must be enabled !!\n" );

   if ( !OPT__FLAG_NPAR_CELL )
      Aux_Error( ERROR_INFO, "OPT__FLAG_NPAR_CELL must be enabled to make sure the particle is in the MAX_LEVEL for feedback !!\n" );
#  endif


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Validating test problem %d ... done\n", TESTPROB_ID );

} // FUNCTION : Validate



#if ( MODEL == HYDRO  &&  defined FEEDBACK )
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
// ********************************************************************************************************************************
// LOAD_PARA( load_mode, "KEY_IN_THE_FILE",          &VARIABLE,                  DEFAULT,      MIN,              MAX               );
// ********************************************************************************************************************************
   LOAD_PARA( load_mode, "SNBlast_Dens_Bg",          &SNBlast_Dens_Bg,           1.0e-22,      Eps_double,       NoMax_double      );
   LOAD_PARA( load_mode, "SNBlast_Temp_Bg",          &SNBlast_Temp_Bg,           1.0e+01,      Eps_double,       NoMax_double      );
   LOAD_PARA( load_mode, "SNBlast_MetalMassFrac_Bg", &SNBlast_MetalMassFrac_Bg,  0.0,          0.0,              1.0               );
   LOAD_PARA( load_mode, "SNBlast_ParMass",          &SNBlast_ParMass,           1.0e+02,      Eps_double,       NoMax_double      );
   LOAD_PARA( load_mode, "SNBlast_ParCenter_X",      &SNBlast_ParCenter[0],     -1.0,          NoMin_double,     amr->BoxSize[0]   );
   LOAD_PARA( load_mode, "SNBlast_ParCenter_Y",      &SNBlast_ParCenter[1],     -1.0,          NoMin_double,     amr->BoxSize[1]   );
   LOAD_PARA( load_mode, "SNBlast_ParCenter_Z",      &SNBlast_ParCenter[2],     -1.0,          NoMin_double,     amr->BoxSize[2]   );
   LOAD_PARA( load_mode, "SNBlast_ParVelocity_X",    &SNBlast_ParVelocity[0],    0.0,          NoMin_double,     NoMax_double      );
   LOAD_PARA( load_mode, "SNBlast_ParVelocity_Y",    &SNBlast_ParVelocity[1],    0.0,          NoMin_double,     NoMax_double      );
   LOAD_PARA( load_mode, "SNBlast_ParVelocity_Z",    &SNBlast_ParVelocity[2],    0.0,          NoMin_double,     NoMax_double      );
   LOAD_PARA( load_mode, "SNBlast_ParMetalMassFrac", &SNBlast_ParMetalMassFrac,  0.0,          0.0,              1.0               );

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
      if ( SNBlast_ParCenter[d] < 0.0 )  SNBlast_ParCenter[d] = amr->BoxCenter[d];

// convert to code units
   SNBlast_Dens_Bg /= UNIT_D;
   SNBlast_ParMass *= Const_Msun / UNIT_M;


// (2) set the problem-specific derived parameters
// must initialize EoS first
   EoS_Init();

   if ( EoS_DensTemp2Pres_CPUPtr == NULL )
      Aux_Error( ERROR_INFO, "EoS_DensTemp2Pres_CPUPtr == NULL !!\n" );

   if ( EoS_DensPres2Eint_CPUPtr == NULL )
      Aux_Error( ERROR_INFO, "EoS_DensPres2Eint_CPUPtr == NULL !!\n" );

   if ( EoS_DensEint2Pres_CPUPtr == NULL )
      Aux_Error( ERROR_INFO, "EoS_DensEint2Pres_CPUPtr == NULL !!\n" );

   SNBlast_Pres_Bg      = EoS_DensTemp2Pres_CPUPtr( SNBlast_Dens_Bg, SNBlast_Temp_Bg, NULL,
                                                    EoS_AuxArray_Flt, EoS_AuxArray_Int, h_EoS_Table ); // assuming EoS requires no passive scalars
   SNBlast_Eint_Bg      = EoS_DensPres2Eint_CPUPtr( SNBlast_Dens_Bg, SNBlast_Pres_Bg, NULL,
                                                    EoS_AuxArray_Flt, EoS_AuxArray_Int, h_EoS_Table ); // assuming EoS requires no passive scalars

   const double ExpEngy = FB_RESOLVED_SNEII_EJECT_ENGY;
   const double ExpDens = SNBlast_Dens_Bg + FB_RESOLVED_SNEII_EJECT_MASS/CUBE(amr->dh[MAX_LEVEL]);
   const double ExpEint = SNBlast_Eint_Bg + ExpEngy/CUBE(amr->dh[MAX_LEVEL]);
   const double ExpPres = EoS_DensEint2Pres_CPUPtr( ExpDens, ExpEint, NULL,
                                                    EoS_AuxArray_Flt, EoS_AuxArray_Int, h_EoS_Table ); // assuming EoS requires no passive scalars


// (3) reset other general-purpose parameters
//     --> a helper macro PRINT_RESET_PARA is defined in Macro.h
   const double End_T_Default    = MAX( 10.0*FB_RESOLVED_SNEII_DELAY_TIME, 10.0*Const_Myr/UNIT_T );
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

      Aux_Message( stdout, "=============================================================================\n" );
      Aux_Message( stdout, "  test problem ID           = %d\n",                        TESTPROB_ID              );
      Aux_Message( stdout, "  background mass density   = %13.7e\n",                    SNBlast_Dens_Bg          );
      Aux_Message( stdout, "                            = %13.7e g/cm^3\n",             SNBlast_Dens_Bg*UNIT_D   );
      Aux_Message( stdout, "  background temperature    = %13.7e K\n",                  SNBlast_Temp_Bg          );
      Aux_Message( stdout, "  background pressure       = %13.7e\n",                    SNBlast_Pres_Bg          );
      Aux_Message( stdout, "  background energy density = %13.7e\n",                    SNBlast_Eint_Bg          );
      Aux_Message( stdout, "  background metal fraction = %13.7e\n",                    SNBlast_MetalMassFrac_Bg );
      Aux_Message( stdout, "  SN explosion energy       = %13.7e\n",                    ExpEngy                  );
      Aux_Message( stdout, "                            = %13.7e erg\n",                ExpEngy*UNIT_E           );
      Aux_Message( stdout, "  SN explosion pressure     = %13.7e\n",                    ExpPres                  );
      Aux_Message( stdout, "  SN initial mass           = %13.7e\n",                    SNBlast_ParMass          );
      Aux_Message( stdout, "  SN initial center         = (%13.7e, %13.7e, %13.7e)\n",  SNBlast_ParCenter[0],
                                                                                        SNBlast_ParCenter[1],
                                                                                        SNBlast_ParCenter[2]     );
      Aux_Message( stdout, "  SN initial velocity       = (%13.7e, %13.7e, %13.7e)\n",  SNBlast_ParVelocity[0],
                                                                                        SNBlast_ParVelocity[1],
                                                                                        SNBlast_ParVelocity[2]   );
      Aux_Message( stdout, "  SN initial metal fracion  = %13.7e\n",                    SNBlast_ParMetalMassFrac );
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

// check
#  ifdef GAMER_DEBUG
   if ( Idx_Metal == Idx_Undefined )
      Aux_Error( ERROR_INFO, "Idx_Metal is undefined !!\n" );
#  endif

   double Dens, MomX, MomY, MomZ, Eint, Etot, Metal=NULL_REAL;

   Dens  = SNBlast_Dens_Bg;
   MomX  = 0.0;
   MomY  = 0.0;
   MomZ  = 0.0;
   Eint  = SNBlast_Eint_Bg;
   Etot  = Hydro_ConEint2Etot( Dens, MomX, MomY, MomZ, Eint, 0.0 );     // do NOT include magnetic energy here
   Metal = Dens*SNBlast_MetalMassFrac_Bg;

// set the output array
   fluid[DENS]      = Dens;
   fluid[MOMX]      = MomX;
   fluid[MOMY]      = MomY;
   fluid[MOMZ]      = MomZ;
   fluid[ENGY]      = Etot;
   fluid[Idx_Metal] = Metal;

} // FUNCTION : SetGridIC



//-------------------------------------------------------------------------------------------------------
// Function    :  AddNewField_SNFeedbackBlastWave
// Description :  Add the problem-specific grid fields
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
void AddNewField_SNFeedbackBlastWave()
{

// add the metallicity field only if it has not been done
// --> since Grackle may already add this field automatically when GRACKLE_METAL is enabled
// --> also note that "Idx_Metal" has been predefined in Field.h
   if ( Idx_Metal == Idx_Undefined )
      Idx_Metal = AddField( "Metal", FIXUP_FLUX_YES, FIXUP_REST_YES, FLOOR_YES, NORMALIZE_NO, INTERP_FRAC_YES );

} // FUNCTION : AddNewField_SNFeedbackBlastWave



#ifdef MASSIVE_PARTICLES
//-------------------------------------------------------------------------------------------------------
// Function    :  AddNewParticleAttribute_SNFeedbackBlastWave
// Description :  Add the problem-specific particle attributes
//
// Note        :  1. Ref: https://github.com/gamer-project/gamer/wiki/Adding-New-Simulations#v-add-problem-specific-grid-fields-and-particle-attributes
//                2. Invoke AddParticleField() for each of the problem-specific particle attribute:
//                   --> Attribute label sent to AddParticleField() will be used as the output name of the attribute
//                   --> Attribute index returned by AddParticleField() can be used to access the particle attribute data
//                3. Pre-declared attribute indices are put in Field.h
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void AddNewParticleAttribute_SNFeedbackBlastWave()
{

// "Idx_ParMetalFrac" has been predefined in Field.h
   if ( Idx_ParMetalFrac == Idx_Undefined )
      Idx_ParMetalFrac = AddParticleAttributeFlt( "ParMetalFrac" );

} // FUNCTION : AddNewParticleAttribute_SNFeedbackBlastWave
#endif // #ifdef MASSIVE_PARTICLES
#endif // #if ( MODEL == HYDRO  &&  defined FEEDBACK )



//-------------------------------------------------------------------------------------------------------
// Function    :  Init_TestProb_Hydro_SNFeedbackBlastWave
// Description :  Test problem initializer
//
// Note        :  None
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Init_TestProb_Hydro_SNFeedbackBlastWave()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );


// validate the compilation flags and runtime parameters
   Validate();


#  if ( MODEL == HYDRO  &&  defined FEEDBACK )
// set the problem-specific runtime parameters
   SetParameter();


// set the function pointers of various problem-specific routines
   Init_Function_User_Ptr        = SetGridIC;
   Flag_User_Ptr                 = Flag_SNFeedbackBlastWave;
   Init_Field_User_Ptr           = AddNewField_SNFeedbackBlastWave;
#  ifdef MASSIVE_PARTICLES
   Par_Init_ByFunction_Ptr       = Par_Init_ByFunction_SNFeedbackBlastWave;
   Par_Init_Attribute_User_Ptr   = AddNewParticleAttribute_SNFeedbackBlastWave;
#  ifdef SUPPORT_HDF5
   Output_HDF5_InputTest_Ptr     = LoadInputTestProb;
#  endif
#  endif
#  endif // #if ( MODEL == HYDRO  &&  defined FEEDBACK )


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : Init_TestProb_Hydro_SNFeedbackBlastWave
