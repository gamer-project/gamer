#include "GAMER.h"
#include "Par_EquilibriumIC.h"
#include "string"

using namespace std;

// negligibly small uniform density and energy
double ParEqmIC_SmallGas;
int    ParEqmIC_CloudNum;
char   (*ParEqmIC_Params_Filenames)[MAX_STRING]        = NULL;

double (*ParEqmIC_Cloud_Center)[3]                     = NULL;
double (*ParEqmIC_Cloud_BulkVel)[3]                    = NULL;
char   (*ParEqmIC_Cloud_Type)[MAX_STRING]              = NULL;
double  *ParEqmIC_Cloud_Rho0                           = NULL;
double  *ParEqmIC_Cloud_R0                             = NULL;
double  *ParEqmIC_Cloud_Einasto_Power_Factor           = NULL;
char   (*ParEqmIC_Density_Table_Name)[MAX_STRING]      = NULL;
double  *ParEqmIC_Cloud_Par_Num_Ratio                  = NULL;
long    *ParEqmIC_Cloud_Par_Num                        = NULL;
double  *ParEqmIC_Cloud_MaxR                           = NULL;
int     *ParEqmIC_Cloud_MassProfNBin                   = NULL;
int     *ParEqmIC_Cloud_RSeed                          = NULL;
int     *ParEqmIC_AddExtPot                            = NULL;
char   (*ParEqmIC_ExtPot_Table_Name)[MAX_STRING]       = NULL;

// problem-specific function prototypes
#ifdef MASSIVE_PARTICLES
void Par_Init_ByFunction_ParEqmIC( const long NPar_ThisRank, const long NPar_AllRank,
                                   real *ParMass, real *ParPosX, real *ParPosY, real *ParPosZ,
                                   real *ParVelX, real *ParVelY, real *ParVelZ, real *ParTime,
                                   real *ParType, real *AllAttribute[PAR_NATT_TOTAL] );
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
//     --> a helper macro PRINT_RESET_PARA is defined in Macro.h
   const long   End_Step_Default = __INT_MAX__;
   const double End_T_Default    =  10;

   if ( END_STEP < 0 ) {
      END_STEP = End_Step_Default;
      PRINT_RESET_PARA( END_STEP, FORMAT_LONG, "" );
   }

   if ( END_T < 0.0 ) {
      END_T = End_T_Default;
      PRINT_RESET_PARA( END_T, FORMAT_REAL, "" );
   }

// load run-time parameters
   const char* FileName = "Input__TestProb";
   ReadPara_t *ReadPara  = new ReadPara_t;
   // ********************************************************************************************************************************
   // ReadPara->Add( "KEY_IN_THE_FILE",      &VARIABLE,              DEFAULT,       MIN,              MAX               );
   // ********************************************************************************************************************************
   ReadPara->Add( "ParEqmIC_SmallGas",       &ParEqmIC_SmallGas,     1e-3,          0.,               NoMax_double      );
   ReadPara->Add( "ParEqmIC_CloudNum",       &ParEqmIC_CloudNum,     1,             1,                NoMax_int         );
   ReadPara->Read( FileName );

   ParEqmIC_Params_Filenames = new char [ParEqmIC_CloudNum][MAX_STRING];

   char ParEqmIC_Params_Filename_i[MAX_STRING];

   for (int i=0; i<ParEqmIC_CloudNum; i++)
   {
      sprintf( ParEqmIC_Params_Filename_i, "ParEqmIC_Params_Filename_%d", i+1 );
      ReadPara->Add( ParEqmIC_Params_Filename_i,  ParEqmIC_Params_Filenames[i],  NoDef_str,            Useless_str,   Useless_str      );
   }

   ReadPara->Read( FileName );

   delete ReadPara;


   if ( MPI_Rank == 0 )
   {
      Aux_Message( stdout, "Checking Params_Filenames...\n" );

      fstream file;

      for (int i=0; i<ParEqmIC_CloudNum; i++)
      {
         file.open( ParEqmIC_Params_Filenames[i], ios::in );

         if ( !file )   Aux_Message( stdout, "ParEqmIC_Params_Filenames file %s cannot be found !!\n", ParEqmIC_Params_Filenames[i] );

         file.close();
      }

      Aux_Message( stdout, "Checking Params_Filenames... done\n" );
   }

   ParEqmIC_Cloud_Center               = new double [ParEqmIC_CloudNum][3];
   ParEqmIC_Cloud_BulkVel              = new double [ParEqmIC_CloudNum][3];
   ParEqmIC_Cloud_Type                 = new char   [ParEqmIC_CloudNum][MAX_STRING];
   ParEqmIC_Cloud_Rho0                 = new double [ParEqmIC_CloudNum];
   ParEqmIC_Cloud_R0                   = new double [ParEqmIC_CloudNum];
   ParEqmIC_Cloud_Einasto_Power_Factor = new double [ParEqmIC_CloudNum];
   ParEqmIC_Density_Table_Name         = new char   [ParEqmIC_CloudNum][MAX_STRING];
   ParEqmIC_Cloud_Par_Num_Ratio        = new double [ParEqmIC_CloudNum];
   ParEqmIC_Cloud_Par_Num              = new long   [ParEqmIC_CloudNum];
   ParEqmIC_Cloud_MaxR                 = new double [ParEqmIC_CloudNum];
   ParEqmIC_Cloud_MassProfNBin         = new int    [ParEqmIC_CloudNum];
   ParEqmIC_Cloud_RSeed                = new int    [ParEqmIC_CloudNum];
   ParEqmIC_AddExtPot                  = new int    [ParEqmIC_CloudNum];
   ParEqmIC_ExtPot_Table_Name          = new char   [ParEqmIC_CloudNum][MAX_STRING];

   for (int i=0; i<ParEqmIC_CloudNum; i++)
   {
      if ( MPI_Rank == 0 )   Aux_Message( stdout, "Reading physical parameters input file (%d): %s\n", i, ParEqmIC_Params_Filenames[i] );

      // (1) load the problem-specific runtime parameters
      ReadPara_t *ReadPara  = new ReadPara_t;

      // (1-1) add parameters in the following format:
      // --> note that VARIABLE, DEFAULT, MIN, and MAX must have the same data type
      // --> some handy constants (e.g., Useless_bool, Eps_double, NoMin_int, ...) are defined in "include/ReadPara.h"
      // ********************************************************************************************************************************
      // ReadPara->Add( "KEY_IN_THE_FILE",         &VARIABLE,                               DEFAULT,       MIN,              MAX               );
      // ********************************************************************************************************************************
      ReadPara->Add( "Cloud_CenterX",              &ParEqmIC_Cloud_Center[i][0],            NoDef_double,  NoMin_double,     NoMax_double      );
      ReadPara->Add( "Cloud_CenterY",              &ParEqmIC_Cloud_Center[i][1],            NoDef_double,  NoMin_double,     NoMax_double      );
      ReadPara->Add( "Cloud_CenterZ",              &ParEqmIC_Cloud_Center[i][2],            NoDef_double,  NoMin_double,     NoMax_double      );
      ReadPara->Add( "Cloud_BulkVelX",             &ParEqmIC_Cloud_BulkVel[i][0],           0.0,           NoMin_double,     NoMax_double      );
      ReadPara->Add( "Cloud_BulkVelY",             &ParEqmIC_Cloud_BulkVel[i][1],           0.0,           NoMin_double,     NoMax_double      );
      ReadPara->Add( "Cloud_BulkVelZ",             &ParEqmIC_Cloud_BulkVel[i][2],           0.0,           NoMin_double,     NoMax_double      );
      ReadPara->Add( "Cloud_Type",                  ParEqmIC_Cloud_Type[i],                 NoDef_str,     Useless_str,      Useless_str       );
      ReadPara->Add( "Cloud_Rho0",                 &ParEqmIC_Cloud_Rho0[i],                 1.0,           Eps_double,       NoMax_double      );
      ReadPara->Add( "Cloud_R0",                   &ParEqmIC_Cloud_R0[i],                   0.1,           Eps_double,       NoMax_double      );
      ReadPara->Add( "Cloud_Einasto_Power_Factor", &ParEqmIC_Cloud_Einasto_Power_Factor[i], 1.0,           0.1,              10.0              );
      ReadPara->Add( "Density_Table_Name",          ParEqmIC_Density_Table_Name[i],         NoDef_str,     Useless_str,      Useless_str       );
      ReadPara->Add( "Cloud_Par_Num_Ratio",        &ParEqmIC_Cloud_Par_Num_Ratio[i],        0.0,           0.0,              1.0               );
      ReadPara->Add( "Cloud_MaxR",                 &ParEqmIC_Cloud_MaxR[i],                 0.375,         Eps_double,       NoMax_double      );
      ReadPara->Add( "Cloud_MassProfNBin",         &ParEqmIC_Cloud_MassProfNBin[i],         1000,          2,                NoMax_int         );
      ReadPara->Add( "Cloud_RSeed",                &ParEqmIC_Cloud_RSeed[i],                123,           0,                NoMax_int         );
      ReadPara->Add( "AddExtPot",                  &ParEqmIC_AddExtPot[i],                  0,             0,                1                 );
      ReadPara->Add( "ExtPot_Table_Name",           ParEqmIC_ExtPot_Table_Name[i],          NoDef_str,     Useless_str,      Useless_str       );

      ReadPara->Read( ParEqmIC_Params_Filenames[i] );
      delete ReadPara;

      // (1-2) set the default values
      for (int d=0; d<3; d++)
         if ( ParEqmIC_Cloud_Center[i][d] == NoDef_double )  ParEqmIC_Cloud_Center[i][d] = amr->BoxCenter[d];

      // (2) make a note
      if ( MPI_Rank == 0 )
      {
         if ( i == 0 )
         {
         Aux_Message( stdout, "=============================================================================\n"          );
         Aux_Message( stdout, "  test problem ID                           = %d\n",     TESTPROB_ID                      );
         Aux_Message( stdout, "=============================================================================\n"          );
         }

         Aux_Message( stdout, "=============================================================================\n"          );
         Aux_Message( stdout, "  cloud ID                                  = %d\n",     i                                );
         Aux_Message( stdout, "  random seed for setting particle position = %d\n",     ParEqmIC_Cloud_RSeed[i]          );
         Aux_Message( stdout, "  peak density                              = %13.7e\n", ParEqmIC_Cloud_Rho0[i]           );
         Aux_Message( stdout, "  scale radius                              = %13.7e\n", ParEqmIC_Cloud_R0[i]             );
         Aux_Message( stdout, "  maximum radius of particles               = %13.7e\n", ParEqmIC_Cloud_MaxR[i]           );

         for (int d=0; d<3; d++)
         {
         Aux_Message( stdout, "  central coordinate [%d]                    = %13.7e\n", d, ParEqmIC_Cloud_Center[i][d]  );
         Aux_Message( stdout, "  bulk velocity [%d]                         = %13.7e\n", d, ParEqmIC_Cloud_BulkVel[i][d] );
         }

         if ( strcmp( ParEqmIC_Cloud_Type[i], "Table" ) != 0 )
         Aux_Message( stdout, "  number of radial bins in the mass profile = %d\n",     ParEqmIC_Cloud_MassProfNBin[i]   );

         Aux_Message( stdout, "  Cloud_Type                                = %s\n",     ParEqmIC_Cloud_Type[i]           );
         Aux_Message( stdout, "=============================================================================\n"          );

      }//if ( MPI_Rank == 0 )

      // (3) Warn against small R0
      if ( ParEqmIC_Cloud_R0[i] < amr->dh[MAX_LEVEL] )
         Aux_Message( stdout, "WARNING : Characteristic length R0:%f is smaller than spatial resolution %f!\n", ParEqmIC_Cloud_R0[i], amr->dh[MAX_LEVEL] );

      if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Setting runtime parameters ... done\n" );

      // (4) Check Cloud_Type and table filenames
      // Checking Cloud_Type
      if ( MPI_Rank == 0 )   Aux_Message( stdout, "Checking Cloud_Type\n" );

      int flag = 0;

      if      ( strcmp( ParEqmIC_Cloud_Type[i], "Plummer"   ) == 0 )   flag = 1;
      else if ( strcmp( ParEqmIC_Cloud_Type[i], "NFW"       ) == 0 )   flag = 1;
      else if ( strcmp( ParEqmIC_Cloud_Type[i], "Burkert"   ) == 0 )   flag = 1;
      else if ( strcmp( ParEqmIC_Cloud_Type[i], "Jaffe"     ) == 0 )   flag = 1;
      else if ( strcmp( ParEqmIC_Cloud_Type[i], "Hernquist" ) == 0 )   flag = 1;
      else if ( strcmp( ParEqmIC_Cloud_Type[i], "Einasto"   ) == 0 )   flag = 1;
      else if ( strcmp( ParEqmIC_Cloud_Type[i], "Table"     ) == 0 )   flag = 1;

      if ( flag == 0 )   Aux_Error( ERROR_INFO, "Error in the input of Cloud_Type: %s is not a Model Type !!\n", ParEqmIC_Cloud_Type[i] );

      // Checking Density_Table_Name
      if ( MPI_Rank == 0 )   Aux_Message( stdout, "Checking Density_Table_Name\n" );

      if ( strcmp( ParEqmIC_Cloud_Type[i], "Table" ) == 0 )
      {
         fstream file;
         file.open( ParEqmIC_Density_Table_Name[i], ios::in );

         if ( !file )   Aux_Error( ERROR_INFO, "Error in the input of Density_Table_Name: Density profile %s cannot be found !!\n", ParEqmIC_Density_Table_Name[i] );

         file.close();
      }

      // Checking ExtPot_Table_Name
      if ( MPI_Rank == 0 )   Aux_Message( stdout, "Checking ExtPot_Table_Name\n" );

      if ( ParEqmIC_AddExtPot[i] )
      {
         fstream file;
         file.open( ParEqmIC_ExtPot_Table_Name[i], ios::in );

         if ( !file )   Aux_Error( ERROR_INFO, "Error in the input of ExtPot_Table_Name: External potential profile %s cannot be found !!\n", ParEqmIC_ExtPot_Table_Name[i] );

         file.close();
      }

   } // for (int i=0; i<ParEqmIC_CloudNum; i++)

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
#  ifdef GRAVITY
   fluid[ENGY] = ParEqmIC_SmallGas;
#  endif

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
   delete [] ParEqmIC_Params_Filenames;
   delete [] ParEqmIC_Cloud_Center;
   delete [] ParEqmIC_Cloud_BulkVel;
   delete [] ParEqmIC_Cloud_Type;
   delete [] ParEqmIC_Cloud_Rho0;
   delete [] ParEqmIC_Cloud_R0;
   delete [] ParEqmIC_Cloud_Einasto_Power_Factor;
   delete [] ParEqmIC_Density_Table_Name;
   delete [] ParEqmIC_Cloud_Par_Num_Ratio;
   delete [] ParEqmIC_Cloud_Par_Num;
   delete [] ParEqmIC_Cloud_MaxR;
   delete [] ParEqmIC_Cloud_MassProfNBin;
   delete [] ParEqmIC_Cloud_RSeed;
   delete [] ParEqmIC_AddExtPot;
   delete [] ParEqmIC_ExtPot_Table_Name;

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


   Init_Function_User_Ptr  = SetGridIC;
   End_User_Ptr            = End_ParEqmIC;
#  ifdef MASSIVE_PARTICLES
   Par_Init_ByFunction_Ptr = Par_Init_ByFunction_ParEqmIC;
#  endif
#  ifdef GRAVITY
   if ( OPT__EXT_POT == EXT_POT_FUNC )
   Init_ExtPot_Ptr         = Init_ExtPot_ParEqmIC;
#  endif
#  endif // #if ( MODEL == HYDRO )


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : Init_TestProb_Hydro_ParEqmIC
