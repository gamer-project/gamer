#include "GAMER.h"
#include "Par_EquilibriumIC.h"



// problem-specific global variables
// =======================================================================================
static double   ParEqmIC_SmallGas;                                // negligibly small uniform density and energy
       int      ParEqmIC_NumCloud;                                // number of clouds
       char   (*ParEqmIC_Cloud_ParaFilenames)[MAX_STRING] = NULL; // filenames of the parameters for each cloud

       double (*ParEqmIC_Cloud_Center)[3]                 = NULL; // center coordinates of each cloud
       double (*ParEqmIC_Cloud_BulkVel)[3]                = NULL; // bulk velocity of each cloud
       char   (*ParEqmIC_Cloud_Type)[MAX_STRING]          = NULL; // type of each cloud
       double  *ParEqmIC_Cloud_Rho0                       = NULL; // scale density of each cloud
       double  *ParEqmIC_Cloud_R0                         = NULL; // scale radius of each cloud
       double  *ParEqmIC_Cloud_EinastoPowerFactor         = NULL; // Einasto power factor of each cloud
       char   (*ParEqmIC_Cloud_DensityTable)[MAX_STRING]  = NULL; // input density profile table of each cloud
       long    *ParEqmIC_Cloud_ParNum                     = NULL; // number of particles of each cloud
       double  *ParEqmIC_Cloud_MaxR                       = NULL; // maximum radius of each cloud
       int     *ParEqmIC_Cloud_NBin                       = NULL; // number of bins inside Cloud_MaxR in the density profile of each cloud
       int     *ParEqmIC_Cloud_RSeed                      = NULL; // random seed for particles of each cloud
       int     *ParEqmIC_Cloud_AddExtPotAnaly             = NULL; // whether adding analytical external potential for each cloud
       int     *ParEqmIC_Cloud_AddExtPotTable             = NULL; // whether adding external potential table for each cloud
       char   (*ParEqmIC_Cloud_ExtPotTable)[MAX_STRING]   = NULL; // input external potential table of each cloud
// =======================================================================================

// problem-specific function prototypes
#ifdef MASSIVE_PARTICLES
void Par_Init_ByFunction_ParEqmIC( const long NPar_ThisRank, const long NPar_AllRank,
                                   real_par *ParMass, real_par *ParPosX, real_par *ParPosY, real_par *ParPosZ,
                                   real_par *ParVelX, real_par *ParVelY, real_par *ParVelZ, real_par *ParTime,
                                   long_par *ParType, real_par *AllAttributeFlt[PAR_NATT_FLT_TOTAL],
                                   long_par *AllAttributeInt[PAR_NATT_INT_TOTAL] );
#endif

// external potential routines
void Init_ExtPot_ParEqmIC();




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

#  ifndef GRAVITY
   Aux_Error( ERROR_INFO, "GRAVITY must be enabled !!\n" );
#  endif

#  ifdef COMOVING
   Aux_Error( ERROR_INFO, "COMOVING must be disabled !!\n" );
#  endif

#  ifdef PARTICLE
   if ( OPT__INIT == INIT_BY_FUNCTION  &&  amr->Par->Init != PAR_INIT_BY_FUNCTION )
      Aux_Error( ERROR_INFO, "please set PAR_INIT = 1 (by FUNCTION) !!\n" );
#  else
   Aux_Error( ERROR_INFO, "PARTICLE must be enabled !!\n" );
#  endif

#  ifndef SUPPORT_GSL
   Aux_Error( ERROR_INFO, "SUPPORT_GSL must be enabled !!\n" );
#  endif


// warnings
   if ( MPI_Rank == 0 )
   {
      for (int f=0; f<6; f++)
      if ( OPT__BC_FLU[f] == BC_FLU_PERIODIC )
         Aux_Message( stderr, "WARNING : periodic BC for fluid is not recommended for this test !!\n" );

#     ifdef GRAVITY
      if ( OPT__BC_POT == BC_POT_PERIODIC )
         Aux_Message( stderr, "WARNING : periodic BC for gravity is not recommended for this test !!\n" );
#     endif
   } // if ( MPI_Rank == 0 )


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Validating test problem %d ... done\n", TESTPROB_ID );

} // FUNCTION : Validate



#if ( MODEL == HYDRO )
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
// ******************************************************************************************************************************
// LOAD_PARA( load_mode, "KEY_IN_THE_FILE",              &VARIABLE,                          DEFAULT,       MIN,           MAX               );
// ******************************************************************************************************************************
   LOAD_PARA( load_mode, "ParEqmIC_SmallGas",            &ParEqmIC_SmallGas,                 1e-3,          0.,            NoMax_double      );
   LOAD_PARA( load_mode, "ParEqmIC_NumCloud",            &ParEqmIC_NumCloud,                 1,             1,             NoMax_int         );
   for (int i=0; i<ParEqmIC_NumCloud; i++) {
   char ParEqmIC_Cloud_ParaFilename_i[MAX_STRING];
   sprintf( ParEqmIC_Cloud_ParaFilename_i, "ParEqmIC_Cloud_ParaFilename_%d", i+1 );
   LOAD_PARA( load_mode, ParEqmIC_Cloud_ParaFilename_i,   ParEqmIC_Cloud_ParaFilenames[i],   NoDef_str,     Useless_str,   Useless_str       );
   }

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
   const char* FileName = "Input__TestProb";

// load the number of clouds first
   ReadPara_t *ReadPara = new ReadPara_t;
   ReadPara->Add( "ParEqmIC_NumCloud", &ParEqmIC_NumCloud, 1, 1, NoMax_int );
   ReadPara->Read( FileName );
   delete ReadPara;

// load the remaining parameters
   ParEqmIC_Cloud_ParaFilenames = new char [ParEqmIC_NumCloud][MAX_STRING];
   ReadPara = new ReadPara_t;
   LoadInputTestProb( LOAD_READPARA, ReadPara, NULL );
   ReadPara->Read( FileName );

   delete ReadPara;

// (1-2) set the default values

// (1-3) check and reset the runtime parameters
   for (int i=0; i<ParEqmIC_NumCloud; i++)
      if ( !Aux_CheckFileExist( ParEqmIC_Cloud_ParaFilenames[i] ) )
         Aux_Error( ERROR_INFO, "ParEqmIC_Cloud_ParaFilename_%d file %s cannot be found !!\n", i+1, ParEqmIC_Cloud_ParaFilenames[i] );


// (2) set the parameters for each cloud
   ParEqmIC_Cloud_Center             = new double [ParEqmIC_NumCloud][3];
   ParEqmIC_Cloud_BulkVel            = new double [ParEqmIC_NumCloud][3];
   ParEqmIC_Cloud_Type               = new char   [ParEqmIC_NumCloud][MAX_STRING];
   ParEqmIC_Cloud_Rho0               = new double [ParEqmIC_NumCloud];
   ParEqmIC_Cloud_R0                 = new double [ParEqmIC_NumCloud];
   ParEqmIC_Cloud_EinastoPowerFactor = new double [ParEqmIC_NumCloud];
   ParEqmIC_Cloud_DensityTable       = new char   [ParEqmIC_NumCloud][MAX_STRING];
   ParEqmIC_Cloud_ParNum             = new long   [ParEqmIC_NumCloud];
   ParEqmIC_Cloud_MaxR               = new double [ParEqmIC_NumCloud];
   ParEqmIC_Cloud_NBin               = new int    [ParEqmIC_NumCloud];
   ParEqmIC_Cloud_RSeed              = new int    [ParEqmIC_NumCloud];
   ParEqmIC_Cloud_AddExtPotAnaly     = new int    [ParEqmIC_NumCloud];
   ParEqmIC_Cloud_AddExtPotTable     = new int    [ParEqmIC_NumCloud];
   ParEqmIC_Cloud_ExtPotTable        = new char   [ParEqmIC_NumCloud][MAX_STRING];

   for (int i=0; i<ParEqmIC_NumCloud; i++)
   {
      if ( MPI_Rank == 0 )   Aux_Message( stdout, "Reading the input ParEqmIC cloud parameters file for cloud_%d: %s\n", i+1, ParEqmIC_Cloud_ParaFilenames[i] );

//    (2-1) load the problem-specific runtime parameters
      ReadPara_t *ReadPara  = new ReadPara_t;

//    (2-1-1) add parameters in the following format:
//    --> note that VARIABLE, DEFAULT, MIN, and MAX must have the same data type
//    --> some handy constants (e.g., Useless_bool, Eps_double, NoMin_int, ...) are defined in "include/ReadPara.h"
//    ********************************************************************************************************************************
//    ReadPara->Add( "KEY_IN_THE_FILE",            &VARIABLE,                               DEFAULT,       MIN,              MAX               );
//    ********************************************************************************************************************************
      ReadPara->Add( "Cloud_CenterX",              &ParEqmIC_Cloud_Center[i][0],            NoDef_double,  NoMin_double,     NoMax_double      );
      ReadPara->Add( "Cloud_CenterY",              &ParEqmIC_Cloud_Center[i][1],            NoDef_double,  NoMin_double,     NoMax_double      );
      ReadPara->Add( "Cloud_CenterZ",              &ParEqmIC_Cloud_Center[i][2],            NoDef_double,  NoMin_double,     NoMax_double      );
      ReadPara->Add( "Cloud_BulkVelX",             &ParEqmIC_Cloud_BulkVel[i][0],           0.0,           NoMin_double,     NoMax_double      );
      ReadPara->Add( "Cloud_BulkVelY",             &ParEqmIC_Cloud_BulkVel[i][1],           0.0,           NoMin_double,     NoMax_double      );
      ReadPara->Add( "Cloud_BulkVelZ",             &ParEqmIC_Cloud_BulkVel[i][2],           0.0,           NoMin_double,     NoMax_double      );
      ReadPara->Add( "Cloud_Type",                  ParEqmIC_Cloud_Type[i],                 NoDef_str,     Useless_str,      Useless_str       );
      ReadPara->Add( "Cloud_Rho0",                 &ParEqmIC_Cloud_Rho0[i],                 1.0,           Eps_double,       NoMax_double      );
      ReadPara->Add( "Cloud_R0",                   &ParEqmIC_Cloud_R0[i],                   0.1,           Eps_double,       NoMax_double      );
      ReadPara->Add( "Cloud_EinastoPowerFactor",   &ParEqmIC_Cloud_EinastoPowerFactor[i],   1.0,           0.1,              10.0              );
      ReadPara->Add( "Cloud_DensityTable",          ParEqmIC_Cloud_DensityTable[i],         NoDef_str,     Useless_str,      Useless_str       );
      ReadPara->Add( "Cloud_ParNum",               &ParEqmIC_Cloud_ParNum[i],               1L,            1L,               NoMax_long        );
      ReadPara->Add( "Cloud_MaxR",                 &ParEqmIC_Cloud_MaxR[i],                 0.375,         Eps_double,       NoMax_double      );
      ReadPara->Add( "Cloud_NBin",                 &ParEqmIC_Cloud_NBin[i],                 1000,          1,                NoMax_int         );
      ReadPara->Add( "Cloud_RSeed",                &ParEqmIC_Cloud_RSeed[i],                123,           0,                NoMax_int         );
      ReadPara->Add( "Cloud_AddExtPotAnaly",       &ParEqmIC_Cloud_AddExtPotAnaly[i],       0,             0,                1                 );
      ReadPara->Add( "Cloud_AddExtPotTable",       &ParEqmIC_Cloud_AddExtPotTable[i],       0,             0,                1                 );
      ReadPara->Add( "Cloud_ExtPotTable",           ParEqmIC_Cloud_ExtPotTable[i],          NoDef_str,     Useless_str,      Useless_str       );

      ReadPara->Read( ParEqmIC_Cloud_ParaFilenames[i] );

      delete ReadPara;

//    (2-1-2) set the default values
      for (int d=0; d<3; d++)
         if ( ParEqmIC_Cloud_Center[i][d] == NoDef_double )  ParEqmIC_Cloud_Center[i][d] = amr->BoxCenter[d];

//    (2-2) warn against small Cloud_R0
      if ( MPI_Rank == 0  &&  ParEqmIC_Cloud_R0[i] < amr->dh[MAX_LEVEL] )
         Aux_Message( stdout, "WARNING : scale length R0 = %f of cloud_%d is smaller than the highest spatial resolution %f!\n", ParEqmIC_Cloud_R0[i], i+1, amr->dh[MAX_LEVEL] );

//    (2-3) check Cloud_Type
      if (  strcmp( ParEqmIC_Cloud_Type[i], "Plummer"   ) != 0  &&
            strcmp( ParEqmIC_Cloud_Type[i], "NFW"       ) != 0  &&
            strcmp( ParEqmIC_Cloud_Type[i], "Burkert"   ) != 0  &&
            strcmp( ParEqmIC_Cloud_Type[i], "Jaffe"     ) != 0  &&
            strcmp( ParEqmIC_Cloud_Type[i], "Hernquist" ) != 0  &&
            strcmp( ParEqmIC_Cloud_Type[i], "Einasto"   ) != 0  &&
            strcmp( ParEqmIC_Cloud_Type[i], "Table"     ) != 0  )
         Aux_Error( ERROR_INFO, "Incorrect ParEqmIC_Cloud_Type = %s for cloud_%d !!\n", ParEqmIC_Cloud_Type[i], i+1 );

//    (2-4) check Cloud_DensityTable
      if ( strcmp( ParEqmIC_Cloud_Type[i], "Table" ) == 0  &&  !Aux_CheckFileExist( ParEqmIC_Cloud_DensityTable[i] ) )
         Aux_Error( ERROR_INFO, "ParEqmIC_Cloud_DensityTable %s for cloud_%d cannot be found !!\n", ParEqmIC_Cloud_DensityTable[i], i+1 );

//    (2-5) check Cloud_ExtPotTable
      if ( ParEqmIC_Cloud_AddExtPotTable[i]  &&  !Aux_CheckFileExist( ParEqmIC_Cloud_ExtPotTable[i] ) )
         Aux_Error( ERROR_INFO, "ParEqmIC_Cloud_ExtPotTable %s for cloud_%d cannot be found !!\n", ParEqmIC_Cloud_ExtPotTable[i], i+1 );

   } // for (int i=0; i<ParEqmIC_NumCloud; i++)


// (3) reset other general-purpose parameters
//     --> a helper macro PRINT_RESET_PARA is defined in Macro.h
   const long   End_Step_Default = __INT_MAX__;
   const double End_T_Default    = 10.0;

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
      Aux_Message( stdout, "=============================================================================\n"                  );
      Aux_Message( stdout, "  test problem ID                              = %d\n",     TESTPROB_ID                          );
      Aux_Message( stdout, "  small gas                                    = %13.7e\n", ParEqmIC_SmallGas                    );
      Aux_Message( stdout, "  number of clouds                             = %d\n",     ParEqmIC_NumCloud                    );

      for (int i=0; i<ParEqmIC_NumCloud; i++)
      {

      Aux_Message( stdout, "\n" );
      Aux_Message( stdout, "  Cloud ID                                     = %d\n",     i+1                                  );

      Aux_Message( stdout, "     Cloud_Type                                = %s\n",     ParEqmIC_Cloud_Type[i]               );
      if ( strcmp( ParEqmIC_Cloud_Type[i], "Table" ) == 0 )
      Aux_Message( stdout, "     density profile file name                 = %s\n",     ParEqmIC_Cloud_DensityTable[i]       );

      for (int d=0; d<3; d++)
      Aux_Message( stdout, "     center coordinates [%3d]                  = %13.7e\n", d, ParEqmIC_Cloud_Center[i][d]       );
      for (int d=0; d<3; d++)
      Aux_Message( stdout, "     bulk velocity [%3d]                       = %13.7e\n", d, ParEqmIC_Cloud_BulkVel[i][d]      );

      Aux_Message( stdout, "     scale density                             = %13.7e\n", ParEqmIC_Cloud_Rho0[i]               );
      Aux_Message( stdout, "     scale radius                              = %13.7e\n", ParEqmIC_Cloud_R0[i]                 );
      if ( strcmp( ParEqmIC_Cloud_Type[i], "Einasto" ) == 0 )
      Aux_Message( stdout, "     Einasto power factor                      = %13.7e\n", ParEqmIC_Cloud_EinastoPowerFactor[i] );

      Aux_Message( stdout, "     number of particles                       = %ld\n",    ParEqmIC_Cloud_ParNum[i]             );
      Aux_Message( stdout, "     maximum radius of particles               = %13.7e\n", ParEqmIC_Cloud_MaxR[i]               );
      Aux_Message( stdout, "     number of radial bins inside Cloud_MaxR   = %d\n",     ParEqmIC_Cloud_NBin[i]               );
      Aux_Message( stdout, "     random seed for setting particle position = %d\n",     ParEqmIC_Cloud_RSeed[i]              );

      Aux_Message( stdout, "     adding external potential analytical      = %d\n",     ParEqmIC_Cloud_AddExtPotAnaly[i]     );
      Aux_Message( stdout, "     adding external potential table           = %d\n",     ParEqmIC_Cloud_AddExtPotTable[i]     );
      if ( ParEqmIC_Cloud_AddExtPotTable[i] )
      Aux_Message( stdout, "     external potential table file name        = %s\n",     ParEqmIC_Cloud_ExtPotTable[i]        );

      }

      Aux_Message( stdout, "=============================================================================\n"                  );
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

   fluid[DENS] = ParEqmIC_SmallGas;
   fluid[MOMX] = 0;
   fluid[MOMY] = 0;
   fluid[MOMZ] = 0;
   fluid[ENGY] = ParEqmIC_SmallGas;

// just set all passive scalars as zero
   for (int v=NCOMP_FLUID; v<NCOMP_TOTAL; v++)  fluid[v] = 0.0;

} // FUNCTION : SetGridIC



//-------------------------------------------------------------------------------------------------------
// Function    :  End_ParEqmIC
// Description :  Free memory before terminating the program
//
// Note        :  1. Linked to the function pointer "End_User_Ptr" to replace "End_User()"
//
// Parameter   :  None
//-------------------------------------------------------------------------------------------------------
void End_ParEqmIC()
{
   delete [] ParEqmIC_Cloud_ParaFilenames;
   delete [] ParEqmIC_Cloud_Center;
   delete [] ParEqmIC_Cloud_BulkVel;
   delete [] ParEqmIC_Cloud_Type;
   delete [] ParEqmIC_Cloud_Rho0;
   delete [] ParEqmIC_Cloud_R0;
   delete [] ParEqmIC_Cloud_EinastoPowerFactor;
   delete [] ParEqmIC_Cloud_DensityTable;
   delete [] ParEqmIC_Cloud_ParNum;
   delete [] ParEqmIC_Cloud_MaxR;
   delete [] ParEqmIC_Cloud_NBin;
   delete [] ParEqmIC_Cloud_RSeed;
   delete [] ParEqmIC_Cloud_AddExtPotAnaly;
   delete [] ParEqmIC_Cloud_AddExtPotTable;
   delete [] ParEqmIC_Cloud_ExtPotTable;

} // FUNCTION : End_ParEqmIC
#endif // #if ( MODEL == HYDRO )



//-------------------------------------------------------------------------------------------------------
// Function    :  Init_TestProb_Hydro_ParEqmIC
// Description :  Test problem initializer
//
// Note        :  None
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Init_TestProb_Hydro_ParEqmIC()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );


// validate the compilation flags and runtime parameters
   Validate();

#  if ( MODEL == HYDRO )
// set the problem-specific runtime parameters
   SetParameter();


   Init_Function_User_Ptr    = SetGridIC;
   End_User_Ptr              = End_ParEqmIC;
   Init_Function_User_Ptr    = SetGridIC;
#  ifdef MASSIVE_PARTICLES
   Par_Init_ByFunction_Ptr   = Par_Init_ByFunction_ParEqmIC;
#  endif
#  ifdef GRAVITY
   if ( OPT__EXT_POT == EXT_POT_FUNC )
   Init_ExtPot_Ptr           = Init_ExtPot_ParEqmIC;
#  endif
#  ifdef SUPPORT_HDF5
   Output_HDF5_InputTest_Ptr = LoadInputTestProb;
#  endif
#  endif // #if ( MODEL == HYDRO )


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : Init_TestProb_Hydro_ParEqmIC
