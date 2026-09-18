#include "GAMER.h"

#ifdef SUPPORT_GSL
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv2.h>
#endif // #ifdef SUPPORT_GSL



// problem-specific global variables
// =======================================================================================
static int    GrackleTest_DefaultTestMode;    // Default test mode
static double GrackleTest_MassDensity_Min;    // Minimum total mass density in the box (in g cm^-3) [1.0e-29]
static double GrackleTest_MassDensity_Max;    // Maximum total mass density in the box (in g cm^-3) [1.0e-21]
static double GrackleTest_TempOverMMW_Min;    // Minimum temperature over mean molecular weight (T/mu) in the box (in K) [1.0e+00]
static double GrackleTest_TempOverMMW_Max;    // Maximum temperature over mean molecular weight (T/mu) in the box (in K) [1.0e+08]
static double GrackleTest_MFrac_Metal;        // Metal mass fraction    (GRACKLE_METAL only)           [1.295e-2]
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
static double GrackleTest_DustToGasRatio;     // Dust-to-gas mass ratio (GRACKLE_DUST only)            [0.01]
static double GrackleTest_HeatingRate;        // User-provided heating rate (in erg cm^-3 s^-1 n_H^-1) [0.0]
static double GrackleTest_CoolingRate;        // User-provided cooling rate (in erg cm^-3 s^-1 n_H^-2) [0.0]
static double GrackleTest_ExpCoolCoeff;       // Coefficient k_cool for edot = -k_cool*e (Myr^-1; DefaultTestMode=5 only) [1.0]

static double GrackleTest_logDens_Min;        // Minimum log( mass density ) in the box
static double GrackleTest_logDens_Max;        // Maximum log( mass density ) in the box
static double GrackleTest_logDens_Range;      // Range of log ( mass density )
static double GrackleTest_logTemp_Min;        // Minimum log( temperature ) in the box
static double GrackleTest_logTemp_Max;        // Maximum log( temperature ) in the box
static double GrackleTest_logTemp_Range;      // Range of log ( temperature )
// =======================================================================================


#ifdef SUPPORT_GSL
static double DustSat_ComputeSaturationTime( const double T0_K, const double gas_rho_cgs, const double k_per_sec );
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

#  if ( NCOMP_PASSIVE != 14 )
   Aux_Error( ERROR_INFO, "NCOMP_PASSIVE must be 14 !!\n" );
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
   LOAD_PARA( load_mode, "GrackleTest_DefaultTestMode", &GrackleTest_DefaultTestMode,       0,             0,                5                 );
   LOAD_PARA( load_mode, "GrackleTest_MassDensity_Min", &GrackleTest_MassDensity_Min,       1.0e-29,       Eps_double,       NoMax_double      );
   LOAD_PARA( load_mode, "GrackleTest_MassDensity_Max", &GrackleTest_MassDensity_Max,       1.0e-21,       Eps_double,       NoMax_double      );
   LOAD_PARA( load_mode, "GrackleTest_TempOverMMW_Min", &GrackleTest_TempOverMMW_Min,       1.0e+00,       Eps_double,       NoMax_double      );
   LOAD_PARA( load_mode, "GrackleTest_TempOverMMW_Max", &GrackleTest_TempOverMMW_Max,       1.0e+08,       Eps_double,       NoMax_double      );
   LOAD_PARA( load_mode, "GrackleTest_MFrac_Metal",     &GrackleTest_MFrac_Metal,           1.295e-2,      0.0,              1.0               );
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
   LOAD_PARA( load_mode, "GrackleTest_DustToGasRatio",  &GrackleTest_DustToGasRatio,        0.01,          0.0,              1.0               );
   LOAD_PARA( load_mode, "GrackleTest_HeatingRate",     &GrackleTest_HeatingRate,           0.0,           0.0,              NoMax_double      );
   LOAD_PARA( load_mode, "GrackleTest_CoolingRate",     &GrackleTest_CoolingRate,           0.0,           0.0,              NoMax_double      );
   LOAD_PARA( load_mode, "GrackleTest_ExpCoolCoeff",    &GrackleTest_ExpCoolCoeff,          1.0,           0.0,              NoMax_double      );

} // FUNCTION : LoadInputTestProb



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
   if ( !GRACKLE_METAL )
   {
      GrackleTest_MFrac_Metal = 0.0;   PRINT_RESET_PARA( GrackleTest_MFrac_Metal, FORMAT_REAL, "for GRACKLE_METAL disabled" );
   }

   if ( !GRACKLE_DUST )
   {
      GrackleTest_DustToGasRatio = 0.0;PRINT_RESET_PARA( GrackleTest_DustToGasRatio, FORMAT_REAL, "for GRACKLE_DUST disabled" );
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

// note that we cannot reset the global GRACKLE_* parameters here since they have been used in Grackle_Init()
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

      if ( !GRACKLE_COOLING )
         Aux_Error( ERROR_INFO, "GRACKLE_COOLING must be 1 for GrackleTest_DefaultTestMode = %d !!\n",
                    GrackleTest_DefaultTestMode );

      if ( !GRACKLE_USE_V_HEATING_RATE )
         Aux_Error( ERROR_INFO, "GRACKLE_USE_V_HEATING_RATE must be 1 for GrackleTest_DefaultTestMode = %d !!\n",
                    GrackleTest_DefaultTestMode );

      if ( GRACKLE_UV )
         Aux_Error( ERROR_INFO, "GRACKLE_UV must be 0 for GrackleTest_DefaultTestMode = %d !!\n",
                    GrackleTest_DefaultTestMode );

      if ( GRACKLE_CMB_FLOOR )
         Aux_Error( ERROR_INFO, "GRACKLE_CMB_FLOOR must be 0 for GrackleTest_DefaultTestMode = %d !!\n",
                    GrackleTest_DefaultTestMode );

      if ( GRACKLE_PE_HEATING )
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
   else if ( GrackleTest_DefaultTestMode == 4 )
   {

#     ifdef COMOVING
      Aux_Error( ERROR_INFO, "COMOVING must be disabled for GrackleTest_DefaultTestMode = %d !!\n",
                 GrackleTest_DefaultTestMode );
#     endif

      if ( END_STEP != 0 )
      {
         END_STEP = 0;
         PRINT_RESET_PARA( END_STEP, FORMAT_LONG, "" );
      }

      if ( GRACKLE_PRIMORDIAL != GRACKLE_PRI_CHE_CLOUDY )
         Aux_Error( ERROR_INFO, "GRACKLE_PRIMORDIAL must be 0 for GrackleTest_DefaultTestMode = %d !!\n",
                    GrackleTest_DefaultTestMode );

      if ( GRACKLE_METAL )
         Aux_Error( ERROR_INFO, "GRACKLE_METAL must be 0 for GrackleTest_DefaultTestMode = %d !!\n",
                    GrackleTest_DefaultTestMode );

      if ( GRACKLE_COOLING )
         Aux_Error( ERROR_INFO, "GRACKLE_COOLING must be 0 for GrackleTest_DefaultTestMode = %d !!\n",
                    GrackleTest_DefaultTestMode );

      if ( !GRACKLE_UV )
         Aux_Error( ERROR_INFO, "GRACKLE_UV must be 1 for GrackleTest_DefaultTestMode = %d !!\n",
                    GrackleTest_DefaultTestMode );

      if ( GRACKLE_CMB_FLOOR )
         Aux_Error( ERROR_INFO, "GRACKLE_CMB_FLOOR must be 0 for GrackleTest_DefaultTestMode = %d !!\n",
                    GrackleTest_DefaultTestMode );

      if ( GRACKLE_PE_HEATING )
         Aux_Error( ERROR_INFO, "GRACKLE_PE_HEATING must be 0 for GrackleTest_DefaultTestMode = %d !!\n",
                    GrackleTest_DefaultTestMode );

      GrackleTest_MassDensity_Min = 1.0e-29;   PRINT_RESET_PARA( GrackleTest_MassDensity_Min, FORMAT_REAL, "for GrackleTest_DefaultTestMode == 4" );
      GrackleTest_MassDensity_Max = 1.0e-21;   PRINT_RESET_PARA( GrackleTest_MassDensity_Max, FORMAT_REAL, "for GrackleTest_DefaultTestMode == 4" );
      GrackleTest_TempOverMMW_Min = 1.0e+00;   PRINT_RESET_PARA( GrackleTest_TempOverMMW_Min, FORMAT_REAL, "for GrackleTest_DefaultTestMode == 4" );
      GrackleTest_TempOverMMW_Max = 1.0e+08;   PRINT_RESET_PARA( GrackleTest_TempOverMMW_Max, FORMAT_REAL, "for GrackleTest_DefaultTestMode == 4" );
      GrackleTest_HeatingRate     = 0.0;       PRINT_RESET_PARA( GrackleTest_HeatingRate,     FORMAT_REAL, "for GrackleTest_DefaultTestMode == 4" );
      GrackleTest_CoolingRate     = 0.0;       PRINT_RESET_PARA( GrackleTest_CoolingRate,     FORMAT_REAL, "for GrackleTest_DefaultTestMode == 4" );
   }
   else if ( GrackleTest_DefaultTestMode == 5 )
   {
      if ( !GRACKLE_METAL )
         Aux_Error( ERROR_INFO, "GRACKLE_METAL must be 1 for GrackleTest_DefaultTestMode = %d !!\n",
                  GrackleTest_DefaultTestMode );

      if ( !GRACKLE_DUST )
         Aux_Error( ERROR_INFO, "GRACKLE_DUST must be 1 for GrackleTest_DefaultTestMode = %d !!\n",
                  GrackleTest_DefaultTestMode );

      if ( GRACKLE_PRIMORDIAL != GRACKLE_PRI_CHE_CLOUDY )
         Aux_Error( ERROR_INFO, "GRACKLE_PRIMORDIAL must be 0 for GrackleTest_DefaultTestMode = %d !!\n",
                    GrackleTest_DefaultTestMode );

      if ( !GRACKLE_USE_V_HEATING_RATE )
         Aux_Error( ERROR_INFO, "GRACKLE_USE_V_HEATING_RATE must be 1 for GrackleTest_DefaultTestMode = %d !!\n",
                    GrackleTest_DefaultTestMode );

      if ( GRACKLE_UV )
         Aux_Error( ERROR_INFO, "GRACKLE_UV must be 0 for GrackleTest_DefaultTestMode = %d !!\n",
                    GrackleTest_DefaultTestMode );

      if ( GRACKLE_CMB_FLOOR )
         Aux_Error( ERROR_INFO, "GRACKLE_CMB_FLOOR must be 0 for GrackleTest_DefaultTestMode = %d !!\n",
                    GrackleTest_DefaultTestMode );

      if ( GRACKLE_PE_HEATING )
         Aux_Error( ERROR_INFO, "GRACKLE_PE_HEATING must be 0 for GrackleTest_DefaultTestMode = %d !!\n",
                    GrackleTest_DefaultTestMode );

      GrackleTest_MassDensity_Min = 1*( Const_amu / CUBE(Const_cm) );  PRINT_RESET_PARA( GrackleTest_MassDensity_Min, FORMAT_REAL, "for GrackleTest_DefaultTestMode == 5" );
      GrackleTest_MassDensity_Max = 1*( Const_amu / CUBE(Const_cm) );  PRINT_RESET_PARA( GrackleTest_MassDensity_Max, FORMAT_REAL, "for GrackleTest_DefaultTestMode == 5" );
      GrackleTest_TempOverMMW_Min = 1.0e+06/MOLECULAR_WEIGHT;          PRINT_RESET_PARA( GrackleTest_TempOverMMW_Min, FORMAT_REAL, "for GrackleTest_DefaultTestMode == 5" );
      GrackleTest_TempOverMMW_Max = 1.0e+06/MOLECULAR_WEIGHT;          PRINT_RESET_PARA( GrackleTest_TempOverMMW_Max, FORMAT_REAL, "for GrackleTest_DefaultTestMode == 5" );
      GrackleTest_MFrac_Metal     = 1.295e-2;                          PRINT_RESET_PARA( GrackleTest_MFrac_Metal,     FORMAT_REAL, "for GrackleTest_DefaultTestMode == 5" );
      GrackleTest_DustToGasRatio  = 0.01;                              PRINT_RESET_PARA( GrackleTest_DustToGasRatio,  FORMAT_REAL, "for GrackleTest_DefaultTestMode == 5" );
      GrackleTest_HeatingRate     = 0;                                 PRINT_RESET_PARA( GrackleTest_HeatingRate,     FORMAT_REAL, "for GrackleTest_DefaultTestMode == 5" );
      GrackleTest_CoolingRate     = 0;                                 PRINT_RESET_PARA( GrackleTest_CoolingRate,     FORMAT_REAL, "for GrackleTest_DefaultTestMode == 5" );

      if ( END_STEP < 0 )
      {
         END_STEP = __INT_MAX__;
         PRINT_RESET_PARA( END_STEP, FORMAT_LONG, "so that END_T controls the simulation end time for GrackleTest_DefaultTestMode == 5" );
      }

      if ( END_T < 0.0 )
      {
         if ( GrackleTest_ExpCoolCoeff <= 0.0 )
            Aux_Error( ERROR_INFO, "GrackleTest_ExpCoolCoeff must be > 0 to auto-compute END_T for GrackleTest_DefaultTestMode == 5 !!\n" );

         const double T0_K        = GrackleTest_TempOverMMW_Min*MOLECULAR_WEIGHT;
         const double gas_rho_cgs = GrackleTest_MassDensity_Min;
         const double k_per_sec   = GrackleTest_ExpCoolCoeff / Const_Myr;

//       The saturation time is the point at which the dust density stops changing appreciably. Using it as the default END_T
//       ensures the simulation runs long enough to capture the full dust-density decay without running unnecessarily longer.
         const double t_sat_sec = DustSat_ComputeSaturationTime( T0_K, gas_rho_cgs, k_per_sec );

         END_T = t_sat_sec / UNIT_T;
         PRINT_RESET_PARA( END_T, FORMAT_REAL, "to the auto-computed dust-sputtering saturation time for GrackleTest_DefaultTestMode == 5" );
      }

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

   if ( END_STEP != 0  &&  END_T != 0.0  &&  !GRACKLE_COOLING )
      Aux_Error( ERROR_INFO, "GRACKLE_COOLING must be enabled for time evolution (END_STEP = %ld, END_T = %14.7e) in this test !!\n", END_STEP, END_T );

// (4) make a note
   if ( MPI_Rank == 0 )
   {
      Aux_Message( stdout, "=============================================================================\n" );
      Aux_Message( stdout, "  test problem ID                             = %d\n",                           TESTPROB_ID                                 );
      Aux_Message( stdout, "  GrackleTest_DefaultTestMode                 = %d\n",                           GrackleTest_DefaultTestMode                 );
      Aux_Message( stdout, "  GrackleTest_MassDensity_Min                 = %13.7e UNIT_D\n",                GrackleTest_MassDensity_Min                 );
      Aux_Message( stdout, "                                              = %13.7e g/cm^3\n",                GrackleTest_MassDensity_Min*UNIT_D          );
      Aux_Message( stdout, "                                              = %13.7e mH/cm^3\n",               GrackleTest_MassDensity_Min*UNIT_D/Const_mH );
      Aux_Message( stdout, "                                              = %13.7e amu/cm^3\n",              GrackleTest_MassDensity_Min*UNIT_D/Const_amu);
      Aux_Message( stdout, "  GrackleTest_MassDensity_Max                 = %13.7e UNIT_D\n",                GrackleTest_MassDensity_Max                 );
      Aux_Message( stdout, "                                              = %13.7e g/cm^3\n",                GrackleTest_MassDensity_Max*UNIT_D          );
      Aux_Message( stdout, "                                              = %13.7e mH/cm^3\n",               GrackleTest_MassDensity_Max*UNIT_D/Const_mH );
      Aux_Message( stdout, "                                              = %13.7e amu/cm^3\n",              GrackleTest_MassDensity_Max*UNIT_D/Const_amu);
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
      Aux_Message( stdout, "  GrackleTest_DustToGasRatio                  = %13.7e\n",                       GrackleTest_DustToGasRatio                  );
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
   if ( Idx_e     == Idx_Undefined )         Aux_Error( ERROR_INFO, "Idx_e is undefined !!\n" );
   if ( Idx_HI    == Idx_Undefined )         Aux_Error( ERROR_INFO, "Idx_HI is undefined !!\n" );
   if ( Idx_HII   == Idx_Undefined )         Aux_Error( ERROR_INFO, "Idx_HII is undefined !!\n" );
   if ( Idx_HeI   == Idx_Undefined )         Aux_Error( ERROR_INFO, "Idx_HeI is undefined !!\n" );
   if ( Idx_HeII  == Idx_Undefined )         Aux_Error( ERROR_INFO, "Idx_HeII is undefined !!\n" );
   if ( Idx_HeIII == Idx_Undefined )         Aux_Error( ERROR_INFO, "Idx_HeIII is undefined !!\n" );
   if ( Idx_HM    == Idx_Undefined )         Aux_Error( ERROR_INFO, "Idx_HM is undefined !!\n" );
   if ( Idx_H2I   == Idx_Undefined )         Aux_Error( ERROR_INFO, "Idx_H2I is undefined !!\n" );
   if ( Idx_H2II  == Idx_Undefined )         Aux_Error( ERROR_INFO, "Idx_H2II is undefined !!\n" );
   if ( Idx_DI    == Idx_Undefined )         Aux_Error( ERROR_INFO, "Idx_DI is undefined !!\n" );
   if ( Idx_DII   == Idx_Undefined )         Aux_Error( ERROR_INFO, "Idx_DII is undefined !!\n" );
   if ( Idx_HDI   == Idx_Undefined )         Aux_Error( ERROR_INFO, "Idx_HDI is undefined !!\n" );
   if ( Idx_Metal == Idx_Undefined )         Aux_Error( ERROR_INFO, "Idx_Metal is undefined !!\n" );
   if ( Idx_Dust  == Idx_Undefined )         Aux_Error( ERROR_INFO, "Idx_Dust is undefined !!\n" );
   if ( EoS_DensTemp2Pres_CPUPtr == NULL )   Aux_Error( ERROR_INFO, "EoS_DensTemp2Pres_CPUPtr == NULL !!\n" );
#  endif

// compute the gas log10 density     by linear interpolation in x-direction
   const double logDens = GrackleTest_logDens_Min +
                          GrackleTest_logDens_Range*( x / amr->BoxSize[0] );
   const real   Dens    = pow( 10, logDens );

// compute the gas log10 temperature by linear interpolation in y-direction
   const double logTemp = GrackleTest_logTemp_Min +
                          GrackleTest_logTemp_Range*( y / amr->BoxSize[1] );
   const real   Temp    = pow( 10, logTemp );

// compute the gas pressure
   const real   Pres = EoS_DensTemp2Pres_CPUPtr( Dens, Temp, NULL, EoS_AuxArray_Flt, EoS_AuxArray_Int,
                                                 h_EoS_Table ); // assuming EoS requires no passive scalars

// assume no momentum
   const real   MomX = 0.0;
   const real   MomY = 0.0;
   const real   MomZ = 0.0;

// compute the total gas energy
   const real   Eint = EoS_DensPres2Eint_CPUPtr( Dens, Pres, NULL, EoS_AuxArray_Flt,
                                                 EoS_AuxArray_Int, h_EoS_Table );         // assuming EoS requires no passive scalars
   const real   Etot = Hydro_ConEint2Etot( Dens, MomX, MomY, MomZ, Eint, (real)0.0 );     // do NOT include magnetic energy here


// set the output array
   fluid[DENS] = Dens;
   fluid[MOMX] = MomX;
   fluid[MOMY] = MomY;
   fluid[MOMZ] = MomZ;
   fluid[ENGY] = Etot;


// set passive scalars
// 6-species network
   fluid[Idx_e    ] = Dens * (real)GrackleTest_MFrac_e;
   fluid[Idx_HI   ] = Dens * (real)GrackleTest_MFrac_HI;
   fluid[Idx_HII  ] = Dens * (real)GrackleTest_MFrac_HII;
   fluid[Idx_HeI  ] = Dens * (real)GrackleTest_MFrac_HeI;
   fluid[Idx_HeII ] = Dens * (real)GrackleTest_MFrac_HeII;
   fluid[Idx_HeIII] = Dens * (real)GrackleTest_MFrac_HeIII;

// 9-species network
   fluid[Idx_HM   ] = Dens * (real)GrackleTest_MFrac_HM;
   fluid[Idx_H2I  ] = Dens * (real)GrackleTest_MFrac_H2I;
   fluid[Idx_H2II ] = Dens * (real)GrackleTest_MFrac_H2II;

// 12-species network
   fluid[Idx_DI   ] = Dens * (real)GrackleTest_MFrac_DI;
   fluid[Idx_DII  ] = Dens * (real)GrackleTest_MFrac_DII;
   fluid[Idx_HDI  ] = Dens * (real)GrackleTest_MFrac_HDI;

// metallicity for metal cooling
   fluid[Idx_Metal] = Dens * (real)GrackleTest_MFrac_Metal;

// dust
   fluid[Idx_Dust ] = Dens * (real)GrackleTest_DustToGasRatio;

} // FUNCTION : SetGridIC



//-------------------------------------------------------------------------------------------------------
// Function    :  AddNewField_GrackleTest
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
void AddNewField_GrackleTest()
{

// add these fields only if they have not been initialized yet
// --> since Grackle may already add these fields automatically when GRACKLE_PRIMORDIAL or GRACKLE_METAL is enabled
//     in Init/Init_Field.cpp
// --> also note that "Idx_*" has been predefined in Field.h
// --> the purpose of adding these fields here when they were not added in Init/Init_Field.cpp is to
//     fix the total number of passive fields (NCOMP_PASSIVE) to 14, regardless of GRACKLE_PRIMORDIAL, GRACKLE_METAL, and GRACKLE_DUST
   if ( Idx_e     == Idx_Undefined )   Idx_e     = AddField( "Electron", FIXUP_FLUX_YES, FIXUP_REST_YES, FLOOR_YES, NORMALIZE_YES, INTERP_FRAC_YES );
   if ( Idx_HI    == Idx_Undefined )   Idx_HI    = AddField( "HI",       FIXUP_FLUX_YES, FIXUP_REST_YES, FLOOR_YES, NORMALIZE_YES, INTERP_FRAC_YES );
   if ( Idx_HII   == Idx_Undefined )   Idx_HII   = AddField( "HII",      FIXUP_FLUX_YES, FIXUP_REST_YES, FLOOR_YES, NORMALIZE_YES, INTERP_FRAC_YES );
   if ( Idx_HeI   == Idx_Undefined )   Idx_HeI   = AddField( "HeI",      FIXUP_FLUX_YES, FIXUP_REST_YES, FLOOR_YES, NORMALIZE_YES, INTERP_FRAC_YES );
   if ( Idx_HeII  == Idx_Undefined )   Idx_HeII  = AddField( "HeII",     FIXUP_FLUX_YES, FIXUP_REST_YES, FLOOR_YES, NORMALIZE_YES, INTERP_FRAC_YES );
   if ( Idx_HeIII == Idx_Undefined )   Idx_HeIII = AddField( "HeIII",    FIXUP_FLUX_YES, FIXUP_REST_YES, FLOOR_YES, NORMALIZE_YES, INTERP_FRAC_YES );
   if ( Idx_HM    == Idx_Undefined )   Idx_HM    = AddField( "HM",       FIXUP_FLUX_YES, FIXUP_REST_YES, FLOOR_YES, NORMALIZE_YES, INTERP_FRAC_YES );
   if ( Idx_H2I   == Idx_Undefined )   Idx_H2I   = AddField( "H2I",      FIXUP_FLUX_YES, FIXUP_REST_YES, FLOOR_YES, NORMALIZE_YES, INTERP_FRAC_YES );
   if ( Idx_H2II  == Idx_Undefined )   Idx_H2II  = AddField( "H2II",     FIXUP_FLUX_YES, FIXUP_REST_YES, FLOOR_YES, NORMALIZE_YES, INTERP_FRAC_YES );
   if ( Idx_DI    == Idx_Undefined )   Idx_DI    = AddField( "DI",       FIXUP_FLUX_YES, FIXUP_REST_YES, FLOOR_YES, NORMALIZE_YES, INTERP_FRAC_YES );
   if ( Idx_DII   == Idx_Undefined )   Idx_DII   = AddField( "DII",      FIXUP_FLUX_YES, FIXUP_REST_YES, FLOOR_YES, NORMALIZE_YES, INTERP_FRAC_YES );
   if ( Idx_HDI   == Idx_Undefined )   Idx_HDI   = AddField( "HDI",      FIXUP_FLUX_YES, FIXUP_REST_YES, FLOOR_YES, NORMALIZE_YES, INTERP_FRAC_YES );
   if ( Idx_Metal == Idx_Undefined )   Idx_Metal = AddField( "Metal",    FIXUP_FLUX_YES, FIXUP_REST_YES, FLOOR_YES, (GRACKLE_PRIMORDIAL==GRACKLE_PRI_CHE_CLOUDY)?NORMALIZE_NO:NORMALIZE_YES, INTERP_FRAC_YES );
   if ( Idx_Dust  == Idx_Undefined )   Idx_Dust  = AddField( "Dust",     FIXUP_FLUX_YES, FIXUP_REST_YES, FLOOR_YES, NORMALIZE_NO, INTERP_FRAC_YES  );

} // FUNCTION : AddNewField_GrackleTest



//-------------------------------------------------------------------------------------------------------
// Function    :  Mis_GetTimeStep_Dust
// Description :  Estimate the user-defined timestep from the cooling time of a
//                reference cell in the dust test problem.
//
// Note        :  1. This function computes the cooling time using Grackle .
//                2. The reference cell is currently fixed at patch[0] on level 0,
//                   with cell index [0][0][0].
//                3. The returned timestep is a user-defined multiple of the absolute cooling time.
//                4. The input arguments "lv" and "dTime_dt" are currently unused.
//                5. The cooling time fraction ("cool_frac") is chosen to mimic the internal
//                   time-step limit applied by Grackle for dust sputtering: 2% of the absolute
//                   cooling time above 3.0e5 K, and 10% at lower temperatures (see README Note 5).
//
// Parameter   :  lv          : Refinement level (unused here)
//                dTime_dt    : dTime/dt
//
// Return      :  User-defined timestep based on the cooling time
//-------------------------------------------------------------------------------------------------------
static double Mis_GetTimeStep_Dust( const int lv, const double dTime_dt )
{
// This additional timestep constraint is only relevant for dust sputtering
// and is only applied when GrackleTest_DefaultTestMode == 5
   if ( GrackleTest_DefaultTestMode != 5 )
      return HUGE_NUMBER;

// If no energy-decay coefficient is set, there is no physical driver for the
// temperature to change, so the cooling-time-based timestep below is undefined.
// Fall back to a fixed, reasonably small timestep in this case.
   if ( GrackleTest_ExpCoolCoeff == 0 )
      return 0.1*Const_Myr/UNIT_T;


   int FluSg = amr->FluSg[0];
   double Dens = amr->patch[FluSg][0][0]->fluid[DENS][0][0][0];
   double Eint = amr->patch[FluSg][0][0]->fluid[ENGY][0][0][0];

// Estimate gas temperature
   const double sEint_cgs = Eint/Dens * SQR( UNIT_V );
   const double Tgas = ( GAMMA - 1.0 ) * MOLECULAR_WEIGHT * Const_mH / Const_kB * sEint_cgs;

   const double cool_frac = ( Tgas > 3.0e5 ) ? 0.02 : 0.1;
   const double cooling_time = Grackle_GetTimeStep_CoolingTime( 0 ) / DT__GRACKLE_COOLING;

   return cool_frac * cooling_time;

} // FUNCTION : Mis_GetTimeStep_Dust



//-------------------------------------------------------------------------------------------------------
// Function    :  Grackle_vHeatingRate_GrackleTest
// Description :  Function to set Grackle's volumetric heating rate for GrackleTest
//
// Note        :  1. Invoked by Grackle_Prepare() using the function pointer
//                   "Grackle_vHeatingRate_User_Ptr", which must be set by a test problem initializer
//                2. This function will be invoked by multiple OpenMP threads when OPENMP is enabled
//                   --> Please ensure that everything here is thread-safe
//                3. Returned rate should be in unit of erg s^-1 cm^-3
//
// Parameter   :  x/y/z     : Target physical coordinates
//                Time      : Target physical time
//                n_H       : Hydrogen number density in units of cm^-3
//                sEint_Gas : Gas specific internal energy
//
// Return      :  volumetric_heating_rate
//-------------------------------------------------------------------------------------------------------
real_che Grackle_vHeatingRate_GrackleTest( const double x, const double y, const double z, const double Time, const double n_H, const real_che sEint_Gas )
{
   if ( GrackleTest_DefaultTestMode == 5 )
   {
      const real_che rho_cgs   = n_H * Const_mH / GRACKLE_HYDROGEN_MFRAC; // n_H is in cm^-3 ,rho_cgs is in g cm^-3
      const real_che sEint_cgs = sEint_Gas * SQR( UNIT_V );               // sEint_Gas is in code unit, sEint_cgs is in erg g^-1

//    apply the user-defined exponential cooling term:
//    d(rho*e)/dt = -GrackleTest_ExpCoolCoeff * rho * e, where rho = rho_cgs and e = sEint_cgs
      const real_che volumetric_heating_rate = - GrackleTest_ExpCoolCoeff / UNIT_T * sEint_cgs * rho_cgs;
      return volumetric_heating_rate;
   }

   const double   volumetric_heating_rate_0 =       n_H * GrackleTest_HeatingRate; // GrackleTest_HeatingRate has units of erg cm^-3 s^-1 n_H^-1
   const double   volumetric_cooling_rate_0 = n_H * n_H * GrackleTest_CoolingRate; // GrackleTest_CoolingRate has units of erg cm^-3 s^-1 n_H^-2

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
//                Dens_Gas  : Gas density in code units
//                sEint_Gas : Gas specific internal energy in code units
//
// Return      :  temperature_floor
//-------------------------------------------------------------------------------------------------------
real_che Grackle_tempFloor_GrackleTest( const double x, const double y, const double z, const double Time, const real_che Dens_Gas, const real_che sEint_Gas )
{
   const double  Dens_Gas_cgs =  Dens_Gas * UNIT_D;      // convert the unit to g cm^-3
   const double sEint_Gas_cgs = sEint_Gas * SQR(UNIT_V); // convert the unit to cm^2 s^-2

// arbitrary example:
// set a 1e10 K temperature floor to disable evolution for the high-density and high-temperature gas
   const real_che temperature_floor = ( Dens_Gas_cgs > 1.0e-24  &&  sEint_Gas_cgs > 2.0e+12 ) ? 1.0e10 : 0.0;

   return temperature_floor;

} // FUNCTION : Grackle_tempFloor_GrackleTest



// ============================================================
// dust-sputtering saturation-time estimator (GrackleTest_DefaultTestMode == 5 only)
// ============================================================
#ifdef SUPPORT_GSL

//-------------------------------------------------------------------------------------------------------
// Function    :  DustSat_InternalEnergy
// Description :  Evaluate the gas specific internal energy at time t under exponential decay
//
// Note        :  1. Follows e(t) = e0*exp(-k*t), consistent with the user-defined exponential
//                   cooling term applied in Grackle_vHeatingRate_GrackleTest()
//                2. Used by DustSat_ODE_RHS() to drive the dust-sputtering-time ODE
//
// Parameter   :  e0 : initial specific internal energy (erg/g)
//                k  : energy decay rate (s^-1)
//                t  : elapsed time (s)
//
// Return      :  specific internal energy at time t (erg/g)
//-------------------------------------------------------------------------------------------------------
static double DustSat_InternalEnergy( const double e0, const double k, const double t )
{
   return e0*exp( -k*t );
} // FUNCTION : DustSat_InternalEnergy



//-------------------------------------------------------------------------------------------------------
// Function    :  DustSat_SputteringTime
// Description :  Dust sputtering timescale formula, as given in Eq. (6) of
//                Richie et al. 2024, ApJ, 974, 81 ("Dust Survival in Galactic Winds")
//
// Note        :  1. t_sp = 0.17 Gyr * (a/0.1 um) * (1e-27 g/cm^3 / rho) * [(1e6.3 K / T)^omega + 1]
//                2. Grain radius and omega are fixed to the values adopted by Richie et al. 2024
//
// Parameter   :  energy_cgs  : gas specific internal energy (erg/g)
//                gas_rho_cgs : gas mass density (g/cm^3)
//
// Return      :  dust sputtering timescale (s)
//-------------------------------------------------------------------------------------------------------
static double DustSat_SputteringTime( const double energy_cgs, const double gas_rho_cgs )
{
   const double DustSat_GrainRadius_um = 0.1;     // grain radius (um)
   const double DustSat_Omega          = 2.5;     // exponent in the sputtering-time formula

   const double Coeff1 = ( 0.17*1.0e3*Const_Myr )*( DustSat_GrainRadius_um/0.1 )*( 1.0e-27/gas_rho_cgs );
   const double Coeff2 = pow( ( pow(10.0, 6.3)*Const_kB )/( (GAMMA-1.0)*MOLECULAR_WEIGHT*Const_mH * energy_cgs ), DustSat_Omega );
   return Coeff1*( Coeff2 + 1.0 );
} // FUNCTION : DustSat_SputteringTime



struct DustSat_ODEParams
{
   double e0;             // initial specific internal energy (erg/g)
   double k;              // energy decay rate (s^-1)
   double gas_rho_cgs;    // gas mass density (g/cm^3)
};



//-------------------------------------------------------------------------------------------------------
// Function    :  DustSat_ODE_RHS
// Description :  Right-hand side of the dust-density ODE d(rho_d)/dt = -3/tsp * rho_d, used by
//                DustSat_ComputeSaturationTime() to integrate the dust density over time, as given in
//                Eq. (3) of Richie et al. 2024, ApJ, 974, 81 ("Dust Survival in Galactic Winds")
//
// Note        :  1. GSL ODE right-hand-side callback; conforms to the gsl_odeiv2_system interface
//
// Parameter   :  t      : current time (s)
//                y      : current state, y[0] = normalized dust density
//                dydt   : output array for dy/dt
//                params : pointer to a DustSat_ODEParams struct
//
// Return      :  GSL_SUCCESS
//-------------------------------------------------------------------------------------------------------
static int DustSat_ODE_RHS( double t, const double y[], double dydt[], void *params )
{
   const DustSat_ODEParams *p = (const DustSat_ODEParams*) params;

   const double energy = DustSat_InternalEnergy( p->e0, p->k, t );
   const double tsp    = DustSat_SputteringTime( energy, p->gas_rho_cgs );

   dydt[0] = -3.0/tsp * y[0];

   return GSL_SUCCESS;
} // FUNCTION : DustSat_ODE_RHS



//-------------------------------------------------------------------------------------------------------
// Function    :  DustSat_ComputeSaturationTime
// Description :  Integrate the normalized dust-density ODE with GSL's rkf45 stepper and locate the
//                saturation time defined by |rho(t) - rho(t-dt_step)| / rho(t-dt_step) < DustSat_Tol
//
// Note        :  1. The step size dt_step mimics the timestep criterion used by Mis_GetTimeStep_Dust()
//                   in the actual simulation: 2% of the cooling time above 3.0e5 K and 10% at lower
//                   temperatures.
//                2. The saturation check compares the current step against the previous step
//                   (backward-looking), not an extrapolated future value
//
// Parameter   :  T0_K        : initial gas temperature (K)
//                gas_rho_cgs : gas mass density (g/cm^3)
//                k_per_sec   : GrackleTest_ExpCoolCoeff converted to s^-1
//
// Return      :  saturation time in seconds (falls back to DustSat_NCoolingTime*t_cool if not found)
//-------------------------------------------------------------------------------------------------------
static double DustSat_ComputeSaturationTime( const double T0_K, const double gas_rho_cgs, const double k_per_sec )
{
   const double DustSat_Tol            = 1.0e-3;  // saturation tolerance
   const double DustSat_NCoolingTime   = 5.0;     // safety cap (in units of t_cool) if no saturation is found

   const double e0          = Const_kB*T0_K / ( (GAMMA-1.0)*MOLECULAR_WEIGHT*Const_mH );  // specific internal energy (erg/g)
   const double t_cool_sec  = 1.0/k_per_sec;
   const double t_max_sec   = DustSat_NCoolingTime*t_cool_sec;

   DustSat_ODEParams params = { e0, k_per_sec, gas_rho_cgs };

   gsl_odeiv2_system sys = { DustSat_ODE_RHS, NULL, 1, &params };

   gsl_odeiv2_driver *driver = gsl_odeiv2_driver_alloc_y_new( &sys, gsl_odeiv2_step_rkf45, 1.0e-3*t_cool_sec, 1.0e-10, 1.0e-14 );

   double t          = 0.0;
   double y[1]       = { 1.0 };   // normalized dust density; y(0) = 1
   double y_prev     = 1.0;       // dust density at the previous step
   double t_sat_sec  = -1.0;

   while ( t < t_max_sec )
   {
      const double energy_now = DustSat_InternalEnergy( e0, k_per_sec, t );
      const double Tgas_now   = energy_now * (GAMMA-1.0) * MOLECULAR_WEIGHT / Const_kB * Const_mH;
      const double cool_frac  = ( Tgas_now > 3.0e5 ) ? 0.02 : 0.1;
      const double dt_step    = cool_frac * t_cool_sec;

      const double t_next = t + dt_step;
      const int status = gsl_odeiv2_driver_apply( driver, &t, t_next, y );
      if ( status != GSL_SUCCESS )
      {
         Aux_Message( stderr, "WARNING : GSL ODE integration failed (status = %d) in DustSat_ComputeSaturationTime !!\n", status );
         break;
      }

//    (rho_now - rho_prev)/rho_prev < DustSat_Tol, checked only after t > 0.5*t_cool
      if ( t > 0.5*t_cool_sec )
      {
         const double delta = fabs( y[0] - y_prev ) / y_prev;
         if ( delta < DustSat_Tol )
         {
            t_sat_sec = t;
            break;
         }
      }

      y_prev = y[0];
   }

   gsl_odeiv2_driver_free( driver );

   if ( t_sat_sec < 0.0 )
   {
      if ( MPI_Rank == 0 )
         Aux_Message( stderr, "WARNING : dust density does not saturate within %.1f cooling times "
                               "--> fall back to END_T = %.1f cooling times !!\n", DustSat_NCoolingTime, DustSat_NCoolingTime );
      t_sat_sec = t_max_sec;
   }

   return t_sat_sec;

} // FUNCTION : DustSat_ComputeSaturationTime
#endif // #ifdef SUPPORT_GSL
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
   Init_Function_User_Ptr        = SetGridIC;
   Init_Field_User_Ptr           = AddNewField_GrackleTest;
   Grackle_vHeatingRate_User_Ptr = Grackle_vHeatingRate_GrackleTest;
   Grackle_tempFloor_User_Ptr    = Grackle_tempFloor_GrackleTest;
   Mis_GetTimeStep_User_Ptr      = Mis_GetTimeStep_Dust;

#  ifdef SUPPORT_HDF5
   Output_HDF5_InputTest_Ptr     = LoadInputTestProb;
#  endif
#  endif // #if ( MODEL == HYDRO  &&  defined SUPPORT_GRACKLE )


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : Init_TestProb_Hydro_GrackleTest
