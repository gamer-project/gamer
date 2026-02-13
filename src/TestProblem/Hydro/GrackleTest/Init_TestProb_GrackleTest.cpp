#include "GAMER.h"



// problem-specific global variables
// =======================================================================================
static int    GrackleTest_DefaultTestMode;    // Default test mode
static double GrackleTest_MassDensity_Min;    // Minimum total mass density in the box (in g cm^-3) [1.0e-29]
static double GrackleTest_MassDensity_Max;    // Maximum total mass density in the box (in g cm^-3) [1.0e-21]
static double GrackleTest_TempOverMMW_Min;    // Minimum temperature over mean molecular weight (T/mu) in the box (in K) [1.0e+00]
static double GrackleTest_TempOverMMW_Max;    // Maximum temperature over mean molecular weight (T/mu) in the box (in K) [1.0e+08]
static double GrackleTest_MFrac_Metal;        // Metal mass fraction    (GRACKLE_METAL only)           [0.01295]
static double GrackleTest_MFrac_e;            // Electron mass fraction (GRACKLE_PRIMORDIAL >= 1 only) [0.0]
static double GrackleTest_MFrac_HI;           // HI mass fraction       (GRACKLE_PRIMORDIAL >= 1 only) [0.750158]
static double GrackleTest_MFrac_HII;          // HII mass fraction      (GRACKLE_PRIMORDIAL >= 1 only) [0.0]
static double GrackleTest_MFrac_HeI;          // HeI mass fraction      (GRACKLE_PRIMORDIAL >= 1 only) [0.236892]
static double GrackleTest_MFrac_HeII;         // HeII mass fraction     (GRACKLE_PRIMORDIAL >= 1 only) [0.0]
static double GrackleTest_MFrac_HeIII;        // HeIII mass fraction    (GRACKLE_PRIMORDIAL >= 1 only) [0.0]
static double GrackleTest_MFrac_HM;           // HM mass fraction       (GRACKLE_PRIMORDIAL >= 2 only) [0.0]
static double GrackleTest_MFrac_H2I;          // H2I mass fraction      (GRACKLE_PRIMORDIAL >= 2 only) [0.0]
static double GrackleTest_MFrac_H2II;         // H2II mass fraction     (GRACKLE_PRIMORDIAL >= 2 only) [0.0]
static double GrackleTest_MFrac_DI;           // DI mass fraction       (GRACKLE_PRIMORDIAL >= 3 only) [0.0]
static double GrackleTest_MFrac_DII;          // DII mass fraction      (GRACKLE_PRIMORDIAL >= 3 only) [0.0]
static double GrackleTest_MFrac_HDI;          // HDI mass fraction      (GRACKLE_PRIMORDIAL >= 3 only) [0.0]
static double GrackleTest_HeatingRate;        // User-provided heating rate (in erg cm^-3 s^-1 n_H^-1) [0.0]
static double GrackleTest_CoolingRate;        // User-provided cooling rate (in erg cm^-3 s^-1 n_H^-2) [0.0]

static double GrackleTest_logDens_Min;        // Minimum log( mass density ) in the box
static double GrackleTest_logDens_Max;        // Maximum log( mass density ) in the box
static double GrackleTest_logDens_Range;      // Range of log ( mass density )
static double GrackleTest_logTemp_Min;        // Minimum log( temperature ) in the box
static double GrackleTest_logTemp_Max;        // Maximum log( temperature ) in the box
static double GrackleTest_logTemp_Range;      // Range of log ( temperature )
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
#  if ( MODEL != HYDRO )
   Aux_Error( ERROR_INFO, "MODEL != HYDRO !!\n" );
#  endif

#  ifndef SUPPORT_GRACKLE
   Aux_Error( ERROR_INFO, "SUPPORT_GRACKLE must be enabled !!\n" );
#  endif

#  ifdef PARTICLE
   Aux_Error( ERROR_INFO, "PARTICLE must be disabled !!\n" );
#  endif

#  ifdef GRAVITY
   Aux_Error( ERROR_INFO, "GRAVITY must be disabled !!\n" );
#  endif

#  ifdef COMOVING
   Aux_Error( ERROR_INFO, "COMOVING must be disabled !!\n" );
#  endif

   if ( !OPT__UNIT )
      Aux_Error( ERROR_INFO, "OPT__UNIT must be enabled for this test !!\n" );

   if ( !OPT__FREEZE_FLUID )
      Aux_Error( ERROR_INFO, "OPT__FREEZE_FLUID must be enabled for this test !!\n" );

#  ifdef SUPPORT_GRACKLE
   if ( !OPT__UNFREEZE_GRACKLE )
      Aux_Error( ERROR_INFO, "OPT__UNFREEZE_GRACKLE must be enabled for this test !!\n" );

   if ( !GRACKLE_ACTIVATE )
      Aux_Error( ERROR_INFO, "GRACKLE_ACTIVATE must be enabled for this test !!\n" );

   if ( !OPT__OUTPUT_GRACKLE_TEMP )
      Aux_Error( ERROR_INFO, "OPT__OUTPUT_GRACKLE_TEMP must be enabled for this test !!\n" );

   if ( !OPT__OUTPUT_GRACKLE_MU )
      Aux_Error( ERROR_INFO, "OPT__OUTPUT_GRACKLE_MU must be enabled for this test !!\n" );

   if ( !OPT__OUTPUT_GRACKLE_TCOOL )
      Aux_Error( ERROR_INFO, "OPT__OUTPUT_GRACKLE_TCOOL must be enabled for this test !!\n" );

   if ( END_STEP != 0  &&  END_T != 0.0  &&  !GRACKLE_COOLING )
      Aux_Error( ERROR_INFO, "GRACKLE_COOLING must be enabled for time evolution in this test !!\n" );


   if ( MPI_Rank == 0 )
   {
      if ( DT__GRACKLE_COOLING < 0.0 )
         Aux_Message( stderr, "WARNING : it's recommended to set DT__GRACKLE_COOLING for this test !!\n" );

      if ( !OPT__FLAG_COOLING_LEN )
         Aux_Message( stderr, "WARNING : it's recommended to set OPT__FLAG_COOLING_LEN for this test !!\n" );
   } // if ( MPI_Rank == 0 )
#  endif // #ifdef SUPPORT_GRACKLE


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Validating test problem %d ... done\n", TESTPROB_ID );

} // FUNCTION : Validate



#if ( MODEL == HYDRO  &&  defined SUPPORT_GRACKLE )
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
// LOAD_PARA( load_mode, "KEY_IN_THE_FILE",             &VARIABLE,                          DEFAULT,       MIN,              MAX               );
// ********************************************************************************************************************************
   LOAD_PARA( load_mode, "GrackleTest_DefaultTestMode", &GrackleTest_DefaultTestMode,       0,             0,                3                 );
   LOAD_PARA( load_mode, "GrackleTest_MassDensity_Min", &GrackleTest_MassDensity_Min,       1.0e-29,       Eps_double,       NoMax_double      );
   LOAD_PARA( load_mode, "GrackleTest_MassDensity_Max", &GrackleTest_MassDensity_Max,       1.0e-21,       Eps_double,       NoMax_double      );
   LOAD_PARA( load_mode, "GrackleTest_TempOverMMW_Min", &GrackleTest_TempOverMMW_Min,       1.0e+00,       Eps_double,       NoMax_double      );
   LOAD_PARA( load_mode, "GrackleTest_TempOverMMW_Max", &GrackleTest_TempOverMMW_Max,       1.0e+08,       Eps_double,       NoMax_double      );
   LOAD_PARA( load_mode, "GrackleTest_MFrac_Metal",     &GrackleTest_MFrac_Metal,           0.01295,       0.0,              1.0               );
   LOAD_PARA( load_mode, "GrackleTest_MFrac_e",         &GrackleTest_MFrac_e,               0.0,           0.0,              1.0               );
   LOAD_PARA( load_mode, "GrackleTest_MFrac_HI",        &GrackleTest_MFrac_HI,              0.750158,      0.0,              1.0               );
   LOAD_PARA( load_mode, "GrackleTest_MFrac_HII",       &GrackleTest_MFrac_HII,             0.0,           0.0,              1.0               );
   LOAD_PARA( load_mode, "GrackleTest_MFrac_HeI",       &GrackleTest_MFrac_HeI,             0.236892,      0.0,              1.0               );
   LOAD_PARA( load_mode, "GrackleTest_MFrac_HeII",      &GrackleTest_MFrac_HeII,            0.0,           0.0,              1.0               );
   LOAD_PARA( load_mode, "GrackleTest_MFrac_HeIII",     &GrackleTest_MFrac_HeIII,           0.0,           0.0,              1.0               );
   LOAD_PARA( load_mode, "GrackleTest_MFrac_HM",        &GrackleTest_MFrac_HM,              0.0,           0.0,              1.0               );
   LOAD_PARA( load_mode, "GrackleTest_MFrac_H2I",       &GrackleTest_MFrac_H2I,             0.0,           0.0,              1.0               );
   LOAD_PARA( load_mode, "GrackleTest_MFrac_H2II",      &GrackleTest_MFrac_H2II,            0.0,           0.0,              1.0               );
   LOAD_PARA( load_mode, "GrackleTest_MFrac_DI",        &GrackleTest_MFrac_DI,              0.0,           0.0,              1.0               );
   LOAD_PARA( load_mode, "GrackleTest_MFrac_DII",       &GrackleTest_MFrac_DII,             0.0,           0.0,              1.0               );
   LOAD_PARA( load_mode, "GrackleTest_MFrac_HDI",       &GrackleTest_MFrac_HDI,             0.0,           0.0,              1.0               );
   LOAD_PARA( load_mode, "GrackleTest_HeatingRate",     &GrackleTest_HeatingRate,           0.0,           0.0,              NoMax_double      );
   LOAD_PARA( load_mode, "GrackleTest_CoolingRate",     &GrackleTest_CoolingRate,           0.0,           0.0,              NoMax_double      );

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
//                3. Must call EoS_Init() before calling any other EoS routine
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
   if ( ! GRACKLE_METAL )
   {
      GrackleTest_MFrac_Metal = 0.0;   PRINT_RESET_PARA( GrackleTest_MFrac_Metal, FORMAT_REAL, "for GRACKLE_METAL disabled" );
   }

   if ( GRACKLE_PRIMORDIAL < GRACKLE_PRI_CHE_NSPE12 )
   {
      GrackleTest_MFrac_DI    = 0.0;   PRINT_RESET_PARA( GrackleTest_MFrac_DI,  FORMAT_REAL, "for GRACKLE_PRIMORDIAL < GRACKLE_PRI_CHE_NSPE12" );
      GrackleTest_MFrac_DII   = 0.0;   PRINT_RESET_PARA( GrackleTest_MFrac_DII, FORMAT_REAL, "for GRACKLE_PRIMORDIAL < GRACKLE_PRI_CHE_NSPE12" );
      GrackleTest_MFrac_HDI   = 0.0;   PRINT_RESET_PARA( GrackleTest_MFrac_HDI, FORMAT_REAL, "for GRACKLE_PRIMORDIAL < GRACKLE_PRI_CHE_NSPE12" );
   }

   if ( GRACKLE_PRIMORDIAL < GRACKLE_PRI_CHE_NSPE9 )
   {
      GrackleTest_MFrac_HM    = 0.0;   PRINT_RESET_PARA( GrackleTest_MFrac_HM,   FORMAT_REAL, "for GRACKLE_PRIMORDIAL < GRACKLE_PRI_CHE_NSPE9" );
      GrackleTest_MFrac_H2I   = 0.0;   PRINT_RESET_PARA( GrackleTest_MFrac_H2I,  FORMAT_REAL, "for GRACKLE_PRIMORDIAL < GRACKLE_PRI_CHE_NSPE9" );
      GrackleTest_MFrac_H2II  = 0.0;   PRINT_RESET_PARA( GrackleTest_MFrac_H2II, FORMAT_REAL, "for GRACKLE_PRIMORDIAL < GRACKLE_PRI_CHE_NSPE9" );
   }

   if ( GRACKLE_PRIMORDIAL < GRACKLE_PRI_CHE_NSPE6 )
   {
      GrackleTest_MFrac_e     = 0.0;   PRINT_RESET_PARA( GrackleTest_MFrac_e,     FORMAT_REAL, "for GRACKLE_PRIMORDIAL < GRACKLE_PRI_CHE_NSPE6" );
      GrackleTest_MFrac_HI    = 0.0;   PRINT_RESET_PARA( GrackleTest_MFrac_HI,    FORMAT_REAL, "for GRACKLE_PRIMORDIAL < GRACKLE_PRI_CHE_NSPE6" );
      GrackleTest_MFrac_HII   = 0.0;   PRINT_RESET_PARA( GrackleTest_MFrac_HII,   FORMAT_REAL, "for GRACKLE_PRIMORDIAL < GRACKLE_PRI_CHE_NSPE6" );
      GrackleTest_MFrac_HeI   = 0.0;   PRINT_RESET_PARA( GrackleTest_MFrac_HeI,   FORMAT_REAL, "for GRACKLE_PRIMORDIAL < GRACKLE_PRI_CHE_NSPE6" );
      GrackleTest_MFrac_HeII  = 0.0;   PRINT_RESET_PARA( GrackleTest_MFrac_HeII,  FORMAT_REAL, "for GRACKLE_PRIMORDIAL < GRACKLE_PRI_CHE_NSPE6" );
      GrackleTest_MFrac_HeIII = 0.0;   PRINT_RESET_PARA( GrackleTest_MFrac_HeIII, FORMAT_REAL, "for GRACKLE_PRIMORDIAL < GRACKLE_PRI_CHE_NSPE6" );
   }

   if ( !GRACKLE_USE_V_HEATING_RATE )
   {
      GrackleTest_HeatingRate = 0.0;   PRINT_RESET_PARA( GrackleTest_HeatingRate, FORMAT_REAL, "for GRACKLE_USE_V_HEATING_RATE disabled" );
      GrackleTest_CoolingRate = 0.0;   PRINT_RESET_PARA( GrackleTest_CoolingRate, FORMAT_REAL, "for GRACKLE_USE_V_HEATING_RATE disabled" );
   }

   if ( GrackleTest_DefaultTestMode == 0 )
   {
//    keep the user's input
   }
   else if ( GrackleTest_DefaultTestMode == 1 )
   {
      GrackleTest_MassDensity_Min = 1.0e-29;   PRINT_RESET_PARA( GrackleTest_MassDensity_Min, FORMAT_REAL, "for GrackleTest_DefaultTestMode == 1" );
      GrackleTest_MassDensity_Max = 1.0e-21;   PRINT_RESET_PARA( GrackleTest_MassDensity_Max, FORMAT_REAL, "for GrackleTest_DefaultTestMode == 1" );
      GrackleTest_TempOverMMW_Min = 1.0e+00;   PRINT_RESET_PARA( GrackleTest_TempOverMMW_Min, FORMAT_REAL, "for GrackleTest_DefaultTestMode == 1" );
      GrackleTest_TempOverMMW_Max = 1.0e+08;   PRINT_RESET_PARA( GrackleTest_TempOverMMW_Max, FORMAT_REAL, "for GrackleTest_DefaultTestMode == 1" );
      GrackleTest_HeatingRate     = 0.0;       PRINT_RESET_PARA( GrackleTest_HeatingRate,     FORMAT_REAL, "for GrackleTest_DefaultTestMode == 1" );
      GrackleTest_CoolingRate     = 0.0;       PRINT_RESET_PARA( GrackleTest_CoolingRate,     FORMAT_REAL, "for GrackleTest_DefaultTestMode == 1" );
   }
   else if ( GrackleTest_DefaultTestMode == 2 )
   {
      GrackleTest_MassDensity_Min = 1.0e-24;   PRINT_RESET_PARA( GrackleTest_MassDensity_Min, FORMAT_REAL, "for GrackleTest_DefaultTestMode == 2" );
      GrackleTest_MassDensity_Max = 1.0e-24;   PRINT_RESET_PARA( GrackleTest_MassDensity_Max, FORMAT_REAL, "for GrackleTest_DefaultTestMode == 2" );
      GrackleTest_TempOverMMW_Min = 1.0e+04;   PRINT_RESET_PARA( GrackleTest_TempOverMMW_Min, FORMAT_REAL, "for GrackleTest_DefaultTestMode == 2" );
      GrackleTest_TempOverMMW_Max = 1.0e+04;   PRINT_RESET_PARA( GrackleTest_TempOverMMW_Max, FORMAT_REAL, "for GrackleTest_DefaultTestMode == 2" );
      GrackleTest_HeatingRate     = 0.0;       PRINT_RESET_PARA( GrackleTest_HeatingRate,     FORMAT_REAL, "for GrackleTest_DefaultTestMode == 2" );
      GrackleTest_CoolingRate     = 0.0;       PRINT_RESET_PARA( GrackleTest_CoolingRate,     FORMAT_REAL, "for GrackleTest_DefaultTestMode == 2" );
   }
   else if ( GrackleTest_DefaultTestMode == 3 )
   {
      if ( GRACKLE_PRIMORDIAL != GRACKLE_PRI_CHE_CLOUDY )
         Aux_Error( ERROR_INFO, "GRACKLE_PRIMORDIAL must be 0 for GrackleTest_DefaultTestMode = %d !!\n",
                    GrackleTest_DefaultTestMode );

      if ( GRACKLE_COOLING != 1 )
         Aux_Error( ERROR_INFO, "GRACKLE_COOLING must be 1 for GrackleTest_DefaultTestMode = %d !!\n",
                    GrackleTest_DefaultTestMode );

      if ( GRACKLE_USE_V_HEATING_RATE != 1 )
         Aux_Error( ERROR_INFO, "GRACKLE_USE_V_HEATING_RATE must be 1 for GrackleTest_DefaultTestMode = %d !!\n",
                    GrackleTest_DefaultTestMode );

      if ( GRACKLE_UV != 0 )
         Aux_Error( ERROR_INFO, "GRACKLE_UV must be 0 for GrackleTest_DefaultTestMode = %d !!\n",
                    GrackleTest_DefaultTestMode );

      if ( GRACKLE_CMB_FLOOR != 0 )
         Aux_Error( ERROR_INFO, "GRACKLE_CMB_FLOOR must be 0 for GrackleTest_DefaultTestMode = %d !!\n",
                    GrackleTest_DefaultTestMode );

      if ( GRACKLE_PE_HEATING != 0 )
         Aux_Error( ERROR_INFO, "GRACKLE_PE_HEATING must be 0 for GrackleTest_DefaultTestMode = %d !!\n",
                    GrackleTest_DefaultTestMode );

      GrackleTest_MassDensity_Min = 1.0e-28;   PRINT_RESET_PARA( GrackleTest_MassDensity_Min, FORMAT_REAL, "for GrackleTest_DefaultTestMode == 3" );
      GrackleTest_MassDensity_Max = 1.0e-28;   PRINT_RESET_PARA( GrackleTest_MassDensity_Max, FORMAT_REAL, "for GrackleTest_DefaultTestMode == 3" );
      GrackleTest_TempOverMMW_Min = 1.0e+06;   PRINT_RESET_PARA( GrackleTest_TempOverMMW_Min, FORMAT_REAL, "for GrackleTest_DefaultTestMode == 3" );
      GrackleTest_TempOverMMW_Max = 1.0e+06;   PRINT_RESET_PARA( GrackleTest_TempOverMMW_Max, FORMAT_REAL, "for GrackleTest_DefaultTestMode == 3" );
      GrackleTest_MFrac_Metal     = 0.0;       PRINT_RESET_PARA( GrackleTest_MFrac_Metal,     FORMAT_REAL, "for GrackleTest_DefaultTestMode == 3" );
      GrackleTest_HeatingRate     = 1.0e-24;   PRINT_RESET_PARA( GrackleTest_HeatingRate,     FORMAT_REAL, "for GrackleTest_DefaultTestMode == 3" );
      GrackleTest_CoolingRate     = 1.6e-20;   PRINT_RESET_PARA( GrackleTest_CoolingRate,     FORMAT_REAL, "for GrackleTest_DefaultTestMode == 3" );
   }
   else
   {
      Aux_Error( ERROR_INFO, "Unknown GrackleTest_DefaultTestMode = %d !!\n", GrackleTest_DefaultTestMode );
   }

// (1-3) check the runtime parameters
   if ( GrackleTest_MassDensity_Min > GrackleTest_MassDensity_Max )
      Aux_Error( ERROR_INFO, "MassDensity_Min = %14.7e > MassDensity_Max = %14.7e !!\n",
                 GrackleTest_MassDensity_Min, GrackleTest_MassDensity_Max );

   if ( GrackleTest_TempOverMMW_Min > GrackleTest_TempOverMMW_Max )
      Aux_Error( ERROR_INFO, "TempOverMMW_Min = %14.7e > TempOverMMW_Max = %14.7e !!\n",
                 GrackleTest_TempOverMMW_Min, GrackleTest_TempOverMMW_Max );

   if ( GRACKLE_PRIMORDIAL >= GRACKLE_PRI_CHE_NSPE6 )
   {
      const double MFrac_Total =  GrackleTest_MFrac_e   + GrackleTest_MFrac_HI   + GrackleTest_MFrac_HII
                                + GrackleTest_MFrac_HeI + GrackleTest_MFrac_HeII + GrackleTest_MFrac_HeIII
                                + GrackleTest_MFrac_HM  + GrackleTest_MFrac_H2I  + GrackleTest_MFrac_H2II
                                + GrackleTest_MFrac_DI  + GrackleTest_MFrac_DII  + GrackleTest_MFrac_HDI
                                + GrackleTest_MFrac_Metal;

      if ( ! Mis_CompareRealValue( MFrac_Total, 1.0, NULL, false ) )
         Aux_Error( ERROR_INFO, "Sum of mass fraction = %14.7e != 1.0 !!\n", MFrac_Total );
   }



// (2) set the problem-specific derived parameters
// convert to code units
   GrackleTest_MassDensity_Min /= UNIT_D;
   GrackleTest_MassDensity_Max /= UNIT_D;

   GrackleTest_logDens_Min   = log10( GrackleTest_MassDensity_Min );
   GrackleTest_logDens_Max   = log10( GrackleTest_MassDensity_Max );
   GrackleTest_logDens_Range = GrackleTest_logDens_Max - GrackleTest_logDens_Min;
   GrackleTest_logTemp_Min   = log10( GrackleTest_TempOverMMW_Min * MOLECULAR_WEIGHT );   // from Temp/mu to built-in Temp in GAMER
   GrackleTest_logTemp_Max   = log10( GrackleTest_TempOverMMW_Max * MOLECULAR_WEIGHT );
   GrackleTest_logTemp_Range = GrackleTest_logTemp_Max - GrackleTest_logTemp_Min;


// (3) reset other general-purpose parameters
//     --> a helper macro PRINT_RESET_PARA is defined in Macro.h
   const long   End_Step_Default = 10;                     // 10 * DT__GRACKLE_COOLING * cooling time
   const double End_T_Default    = 10.0*Const_Myr/UNIT_T;  // 10 Myr

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
      Aux_Message( stdout, "  test problem ID                             = %d\n",                           TESTPROB_ID                                 );
      Aux_Message( stdout, "  GrackleTest_DefaultTestMode                 = %d\n",                           GrackleTest_DefaultTestMode                 );
      Aux_Message( stdout, "  GrackleTest_MassDensity_Min                 = %13.7e UNIT_D\n",                GrackleTest_MassDensity_Min                 );
      Aux_Message( stdout, "                                              = %13.7e g/cm^3\n",                GrackleTest_MassDensity_Min*UNIT_D          );
      Aux_Message( stdout, "                                              = %13.7e mH/cm^3\n",               GrackleTest_MassDensity_Min*UNIT_D/Const_mH );
      Aux_Message( stdout, "  GrackleTest_MassDensity_Max                 = %13.7e UNIT_D\n",                GrackleTest_MassDensity_Max                 );
      Aux_Message( stdout, "                                              = %13.7e g/cm^3\n",                GrackleTest_MassDensity_Max*UNIT_D          );
      Aux_Message( stdout, "                                              = %13.7e mH/cm^3\n",               GrackleTest_MassDensity_Max*UNIT_D/Const_mH );
      Aux_Message( stdout, "  GrackleTest_TempOverMMW_Min                 = %13.7e K\n",                     GrackleTest_TempOverMMW_Min                 );
      Aux_Message( stdout, "  GrackleTest_TempOverMMW_Max                 = %13.7e K\n",                     GrackleTest_TempOverMMW_Max                 );
      Aux_Message( stdout, "  GrackleTest_MFrac_Metal                     = %13.7e\n",                       GrackleTest_MFrac_Metal                     );
      Aux_Message( stdout, "  GrackleTest_MFrac_e                         = %13.7e\n",                       GrackleTest_MFrac_e                         );
      Aux_Message( stdout, "  GrackleTest_MFrac_HI                        = %13.7e\n",                       GrackleTest_MFrac_HI                        );
      Aux_Message( stdout, "  GrackleTest_MFrac_HII                       = %13.7e\n",                       GrackleTest_MFrac_HII                       );
      Aux_Message( stdout, "  GrackleTest_MFrac_HeI                       = %13.7e\n",                       GrackleTest_MFrac_HeI                       );
      Aux_Message( stdout, "  GrackleTest_MFrac_HeII                      = %13.7e\n",                       GrackleTest_MFrac_HeII                      );
      Aux_Message( stdout, "  GrackleTest_MFrac_HeIII                     = %13.7e\n",                       GrackleTest_MFrac_HeIII                     );
      Aux_Message( stdout, "  GrackleTest_MFrac_HM                        = %13.7e\n",                       GrackleTest_MFrac_HM                        );
      Aux_Message( stdout, "  GrackleTest_MFrac_H2I                       = %13.7e\n",                       GrackleTest_MFrac_H2I                       );
      Aux_Message( stdout, "  GrackleTest_MFrac_H2II                      = %13.7e\n",                       GrackleTest_MFrac_H2II                      );
      Aux_Message( stdout, "  GrackleTest_MFrac_DI                        = %13.7e\n",                       GrackleTest_MFrac_DI                        );
      Aux_Message( stdout, "  GrackleTest_MFrac_DII                       = %13.7e\n",                       GrackleTest_MFrac_DII                       );
      Aux_Message( stdout, "  GrackleTest_MFrac_HDI                       = %13.7e\n",                       GrackleTest_MFrac_HDI                       );
      Aux_Message( stdout, "  GrackleTest_HeatingRate                     = %13.7e erg cm^-3 s^-1 n_H^-1\n", GrackleTest_HeatingRate                     );
      Aux_Message( stdout, "  GrackleTest_CoolingRate                     = %13.7e erg cm^-3 s^-1 n_H^-2\n", GrackleTest_CoolingRate                     );
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
//                   (unless OPT__INIT_GRID_WITH_OMP is disabled)
//                   --> Please ensure that everything here is thread-safe
//                3. Even when DUAL_ENERGY is adopted for HYDRO, one does NOT need to set the dual-energy variable here
//                   --> It will be calculated automatically
//                4. For MHD, do NOT add magnetic energy (i.e., 0.5*B^2) to fluid[ENGY] here
//                   --> It will be added automatically later
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
   if ( GRACKLE_PRIMORDIAL >= GRACKLE_PRI_CHE_NSPE6   &&  Idx_e     == Idx_Undefined )   Aux_Error( ERROR_INFO, "Idx_e is undefined !!\n" );
   if ( GRACKLE_PRIMORDIAL >= GRACKLE_PRI_CHE_NSPE6   &&  Idx_HI    == Idx_Undefined )   Aux_Error( ERROR_INFO, "Idx_HI is undefined !!\n" );
   if ( GRACKLE_PRIMORDIAL >= GRACKLE_PRI_CHE_NSPE6   &&  Idx_HII   == Idx_Undefined )   Aux_Error( ERROR_INFO, "Idx_HII is undefined !!\n" );
   if ( GRACKLE_PRIMORDIAL >= GRACKLE_PRI_CHE_NSPE6   &&  Idx_HeI   == Idx_Undefined )   Aux_Error( ERROR_INFO, "Idx_HeI is undefined !!\n" );
   if ( GRACKLE_PRIMORDIAL >= GRACKLE_PRI_CHE_NSPE6   &&  Idx_HeII  == Idx_Undefined )   Aux_Error( ERROR_INFO, "Idx_HeII is undefined !!\n" );
   if ( GRACKLE_PRIMORDIAL >= GRACKLE_PRI_CHE_NSPE6   &&  Idx_HeIII == Idx_Undefined )   Aux_Error( ERROR_INFO, "Idx_HeIII is undefined !!\n" );
   if ( GRACKLE_PRIMORDIAL >= GRACKLE_PRI_CHE_NSPE9   &&  Idx_HM    == Idx_Undefined )   Aux_Error( ERROR_INFO, "Idx_HM is undefined !!\n" );
   if ( GRACKLE_PRIMORDIAL >= GRACKLE_PRI_CHE_NSPE9   &&  Idx_H2I   == Idx_Undefined )   Aux_Error( ERROR_INFO, "Idx_H2I is undefined !!\n" );
   if ( GRACKLE_PRIMORDIAL >= GRACKLE_PRI_CHE_NSPE9   &&  Idx_H2II  == Idx_Undefined )   Aux_Error( ERROR_INFO, "Idx_H2II is undefined !!\n" );
   if ( GRACKLE_PRIMORDIAL >= GRACKLE_PRI_CHE_NSPE12  &&  Idx_DI    == Idx_Undefined )   Aux_Error( ERROR_INFO, "Idx_DI is undefined !!\n" );
   if ( GRACKLE_PRIMORDIAL >= GRACKLE_PRI_CHE_NSPE12  &&  Idx_DII   == Idx_Undefined )   Aux_Error( ERROR_INFO, "Idx_DII is undefined !!\n" );
   if ( GRACKLE_PRIMORDIAL >= GRACKLE_PRI_CHE_NSPE12  &&  Idx_HDI   == Idx_Undefined )   Aux_Error( ERROR_INFO, "Idx_HDI is undefined !!\n" );
   if ( GRACKLE_METAL                                 &&  Idx_Metal == Idx_Undefined )   Aux_Error( ERROR_INFO, "Idx_Metal is undefined !!\n" );
   if ( EoS_DensTemp2Pres_CPUPtr == NULL )                                               Aux_Error( ERROR_INFO, "EoS_DensTemp2Pres_CPUPtr == NULL !!\n" );
#  endif

// compute the gas log density     by linear interpolation in x-direction
   const double logDens = GrackleTest_logDens_Min +
                          GrackleTest_logDens_Range*( x / amr->BoxSize[0] );
   const double Dens    = pow( 10, logDens );

// compute the gas log temperature by linear interpolation in y-direction
   const double logTemp = GrackleTest_logTemp_Min +
                          GrackleTest_logTemp_Range*( y / amr->BoxSize[1] );
   const double Temp    = pow( 10, logTemp );

// compute the gas pressure
   const double Pres = EoS_DensTemp2Pres_CPUPtr( Dens, Temp, NULL, EoS_AuxArray_Flt, EoS_AuxArray_Int,
                                                 h_EoS_Table ); // assuming EoS requires no passive scalars

// assume no momentum
   const double MomX = 0.0;
   const double MomY = 0.0;
   const double MomZ = 0.0;

// compute the total gas energy
   const double Eint = EoS_DensPres2Eint_CPUPtr( Dens, Pres, NULL, EoS_AuxArray_Flt,
                                                 EoS_AuxArray_Int, h_EoS_Table );   // assuming EoS requires no passive scalars
   const double Etot = Hydro_ConEint2Etot( Dens, MomX, MomY, MomZ, Eint, 0.0 );     // do NOT include magnetic energy here


// set the output array
   fluid[DENS] = Dens;
   fluid[MOMX] = MomX;
   fluid[MOMY] = MomY;
   fluid[MOMZ] = MomZ;
   fluid[ENGY] = Etot;


// set passive scalars
// 6-species network
   if ( GRACKLE_PRIMORDIAL >= GRACKLE_PRI_CHE_NSPE6 )
   {
      fluid[Idx_e    ] = Dens * GrackleTest_MFrac_e;
      fluid[Idx_HI   ] = Dens * GrackleTest_MFrac_HI;
      fluid[Idx_HII  ] = Dens * GrackleTest_MFrac_HII;
      fluid[Idx_HeI  ] = Dens * GrackleTest_MFrac_HeI;
      fluid[Idx_HeII ] = Dens * GrackleTest_MFrac_HeII;
      fluid[Idx_HeIII] = Dens * GrackleTest_MFrac_HeIII;
   }

// 9-species network
   if ( GRACKLE_PRIMORDIAL >= GRACKLE_PRI_CHE_NSPE9 )
   {
      fluid[Idx_HM   ] = Dens * GrackleTest_MFrac_HM;
      fluid[Idx_H2I  ] = Dens * GrackleTest_MFrac_H2I;
      fluid[Idx_H2II ] = Dens * GrackleTest_MFrac_H2II;
   }

// 12-species network
   if ( GRACKLE_PRIMORDIAL >= GRACKLE_PRI_CHE_NSPE12 )
   {
      fluid[Idx_DI   ] = Dens * GrackleTest_MFrac_DI;
      fluid[Idx_DII  ] = Dens * GrackleTest_MFrac_DII;
      fluid[Idx_HDI  ] = Dens * GrackleTest_MFrac_HDI;
   }

// metallicity for metal cooling
   if ( GRACKLE_METAL )
      fluid[Idx_Metal] = Dens * GrackleTest_MFrac_Metal;

} // FUNCTION : SetGridIC



//-------------------------------------------------------------------------------------------------------
// Function    :  Grackle_vHeatingRate_GrackleTest
// Description :  Function to set Grackle's volumetric heating rate for GrackleTest
//
// Note        :  1. Invoked by Grackle_Prepare() using the function pointer
//                   "Grackle_vHeatingRate_User_Ptr", which must be set by a test problem initializer
//                2. This function will be invoked by multiple OpenMP threads when OPENMP is enabled
//                   --> Please ensure that everything here is thread-safe
//                3. Returned rate should be in the unit of erg s^-1 cm^-3
//
// Parameter   :  x/y/z    : Target physical coordinates
//                Time     : Target physical time
//                n_H      : Hydrogen number density, should be in the unit of cm^-3
//
// Return      :  volumetric_heating_rate
//-------------------------------------------------------------------------------------------------------
real_che Grackle_vHeatingRate_GrackleTest( const double x, const double y, const double z, const double Time, const double n_H )
{
   const double   volumetric_heating_rate_0 =       n_H * GrackleTest_HeatingRate; // GrackleTest_HeatingRate is in the unit of erg cm^-3 s^-1 n_H^-1
   const double   volumetric_cooling_rate_0 = n_H * n_H * GrackleTest_CoolingRate; // GrackleTest_CoolingRate is in the unit of erg cm^-3 s^-1 n_H^-2

// assume uniform distribution if it is not specified
   if ( GrackleTest_DefaultTestMode != 3 )
      return volumetric_heating_rate_0 - volumetric_cooling_rate_0;

// an example of spatial distribution only for GrackleTest_DefaultTestMode == 3
// two 2D Gaussian distributions: one heating and one cooling
   const double   Center_H[3]               = { amr->BoxCenter[0]+0.25*amr->BoxSize[0], amr->BoxCenter[1], amr->BoxCenter[2] }; // center of heating
   const double   Center_C[3]               = { amr->BoxCenter[0]-0.25*amr->BoxSize[0], amr->BoxCenter[1], amr->BoxCenter[2] }; // center of cooling
   const double   r_to_center_h             = sqrt( SQR(x-Center_H[0]) + SQR(y-Center_H[1]) );                                  // radius to heating
   const double   r_to_center_c             = sqrt( SQR(x-Center_C[0]) + SQR(y-Center_C[1]) );                                  // radius to cooling
   const double   distr_width               = 0.0625*amr->BoxSize[0];                                                           // width of Gaussian
   const double   distr_r_max               = 0.1250*amr->BoxSize[0];                                                           // cutoff radius
   const real_che volumetric_heating_rate_r = ( r_to_center_h <= distr_r_max ) ? volumetric_heating_rate_0 * exp( -0.5*SQR(r_to_center_h/distr_width) )
                                                                               : (real_che)0.0;
   const real_che volumetric_cooling_rate_r = ( r_to_center_c <= distr_r_max ) ? volumetric_cooling_rate_0 * exp( -0.5*SQR(r_to_center_c/distr_width) )
                                                                               : (real_che)0.0;

   return volumetric_heating_rate_r - volumetric_cooling_rate_r;

} // FUNCTION : Grackle_vHeatingRate_GrackleTest



//-------------------------------------------------------------------------------------------------------
// Function    :  Grackle_tempFloor_GrackleTest
// Description :  Function to set Grackle's temperature floor for GrackleTest
//
// Note        :  1. Invoked by Grackle_Prepare() using the function pointer
//                   "Grackle_tempFloor_User_Ptr", which must be set by a test problem initializer
//                2. This function will be invoked by multiple OpenMP threads when OPENMP is enabled
//                   --> Please ensure that everything here is thread-safe
//                3. Returned temperature should be in units of K
//
// Parameter   :  x/y/z     : Target physical coordinates
//                Time      : Target physical time
//                Dens_Gas  : Gas density, in code units
//                sEint_Gas : Gas specific internal energy, in code units
//
// Return      :  temperature_floor
//-------------------------------------------------------------------------------------------------------
real_che Grackle_tempFloor_GrackleTest( const double x, const double y, const double z, const double Time, const real_che Dens_Gas, const real_che sEint_Gas )
{
   const double  Dens_Gas_cgs =  Dens_Gas * UNIT_D;      // convert the unit to g cm^-3
   const double sEint_Gas_cgs = sEint_Gas * SQR(UNIT_V); // convert the unit to cm^2 s^-2

// arbitrary example:
// set a 1e10 K temperature floor to disable evolution for the high-density and high-temperature gases
   const real_che temperature_floor = ( Dens_Gas_cgs > 1.0e-24  &&  sEint_Gas_cgs > 2.0e+12 ) ? 1.0e10 : 0.0;

   return temperature_floor;

} // FUNCTION : Grackle_tempFloor_GrackleTest
#endif // #if ( MODEL == HYDRO  &&  defined SUPPORT_GRACKLE )



//-------------------------------------------------------------------------------------------------------
// Function    :  Init_TestProb_Hydro_GrackleTest
// Description :  Test problem initializer
//
// Note        :  None
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Init_TestProb_Hydro_GrackleTest()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );


// validate the compilation flags and runtime parameters
   Validate();


#  if ( MODEL == HYDRO  &&  defined SUPPORT_GRACKLE )
// set the problem-specific runtime parameters
   SetParameter();


// set the function pointers of various problem-specific routines
   Init_Function_User_Ptr            = SetGridIC;
   Grackle_vHeatingRate_User_Ptr     = Grackle_vHeatingRate_GrackleTest;
   Grackle_tempFloor_User_Ptr        = Grackle_tempFloor_GrackleTest;
#  ifdef SUPPORT_HDF5
   Output_HDF5_InputTest_Ptr         = LoadInputTestProb;
#  endif
#  endif // #if ( MODEL == HYDRO  &&  defined SUPPORT_GRACKLE )


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : Init_TestProb_Hydro_GrackleTest
