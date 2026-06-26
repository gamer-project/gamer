#include "GAMER.h"


static void Output_ExactCooling();

// problem-specific global variables
// =======================================================================================
static double EC_Temp;
static double EC_Dens;
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

#  if ( MODEL != HYDRO )
   Aux_Error( ERROR_INFO, "MODEL != HYDRO !!\n" );
#  endif

#  ifdef GRAVITY
   Aux_Error( ERROR_INFO, "GRAVITY must be disabled !!\n" );
#  endif

#  ifndef EXACT_COOLING
   Aux_Error( ERROR_INFO, "EXACT_COOLING must be enabled !!\n" );
#  endif

#  ifdef COMOVING
   Aux_Error( ERROR_INFO, "COMOVING must be disabled !!\n" );
#  endif

#  ifdef PARTICLE
   Aux_Error( ERROR_INFO, "PARTICLE must be disabled !!\n" );
#  endif
#  ifdef MHD
   Aux_Error( ERROR_INFO, "MHD must be disabled !!\n" );
#  endif
   
   if ( ! SrcTerms.ExactCooling )   Aux_Error( ERROR_INFO, "SRC_EXACTCOOLING must be enabled !!\n" );

// warnings
   if ( MPI_Rank == 0 )
   {
      if ( !OPT__OUTPUT_USER )   Aux_Message( stderr, "WARNING : OPT__OUTPUT_USER is off !!\n" );
   }

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
// ************************************************************************************************************************
// LOAD_PARA( load_mode, "KEY_IN_THE_FILE",   &VARIABLE,              DEFAULT,       MIN,              MAX               );
// ************************************************************************************************************************
   LOAD_PARA( load_mode, "EC_Temp",           &EC_Temp,               1000000.0,     Eps_double,       NoMax_double      );
   LOAD_PARA( load_mode, "EC_Dens",           &EC_Dens,               1.0,           Eps_double,       NoMax_double      );

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

// (1-3) check the runtime parameters


// (2) set the problem-specific derived parameters


// (3) reset other general-purpose parameters
//     --> a helper macro PRINT_WARNING is defined in TestProb.h
   const long   End_Step_Default = __INT_MAX__;
   const double End_T_Default    = 100.0*Const_Myr/UNIT_T;

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
      Aux_Message( stdout, "  test problem ID           = %d\n",     TESTPROB_ID );
      Aux_Message( stdout, "  EC_Temp                   = %13.7e\n", EC_Temp     );
      Aux_Message( stdout, "  EC_Dens                   = %13.7e\n", EC_Dens     );
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

   const double cl_X   = 0.7;                                      // mass-fraction of hydrogen
   const double cl_Z   = 0.018;                                    // metallicity (in Zsun)
   const double cl_mol = 1.0/(2*cl_X+0.75*(1-cl_X-cl_Z)+cl_Z*0.5); // mean (total) molecular weights

   double Dens, MomX, MomY, MomZ, Pres, Eint, Etot;
// Convert the input number density into mass density rho
   double cl_dens = (EC_Dens*MU_NORM*cl_mol) / UNIT_D;
   double cl_pres = EoS_DensTemp2Pres_CPUPtr( cl_dens, EC_Temp, NULL, EoS_AuxArray_Flt, EoS_AuxArray_Int, h_EoS_Table );

   Dens = cl_dens;
   MomX = 0.0;
   MomY = 0.0;
   MomZ = 0.0;
   Pres = cl_pres;
   Eint = EoS_DensPres2Eint_CPUPtr( Dens, Pres, NULL, EoS_AuxArray_Flt, EoS_AuxArray_Int, h_EoS_Table ); // assuming EoS requires no passive scalars
   Etot = Hydro_ConEint2Etot( Dens, MomX, MomY, MomZ, Eint, 0.0 ); // do NOT include magnetic energy here

// set the output array
   fluid[DENS] = Dens;
   fluid[MOMX] = MomX;
   fluid[MOMY] = MomY;
   fluid[MOMZ] = MomZ;
   fluid[ENGY] = Etot;

} // FUNCTION : SetGridIC



//-------------------------------------------------------------------------------------------------------
// Function    :  OutputExactCooling
// Description :  Output the temperature relative error in the exact cooling problem
//
// Note        :  1. Enabled by the runtime option "OPT__OUTPUT_USER"
//                2. Construct the analytical solution corresponding to the cooling function
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Output_ExactCooling()
{
   const char FileName[] = "Record__CoolingErr";
   static bool FirstTime = true;

// header
   if ( FirstTime )
   {
      if ( MPI_Rank == 0 )
      {
         if ( Aux_CheckFileExist( FileName ) )
            Aux_Message( stderr, "WARNING : file \"%s\" already exists !!\n", FileName );

         FILE *File_User = fopen( FileName, "a" );
         fprintf( File_User, "#%13s%10s ",  "Time", "DumpID" );
         fprintf( File_User, "%14s %14s %14s %14s %14s", "Temp_nume", "Temp_anal", "Err", "Tcool_nume", "Tcool_anal" );
         fprintf( File_User, "\n" );
         fclose( File_User );
      }
      FirstTime = false;
   } // if ( FirstTime )

   const double cl_X         = 0.7;                                      // mass-fraction of hydrogen
   const double cl_Z         = 0.018;                                    // metallicity (in Zsun)
   const double cl_mol       = 1.0/(2*cl_X+0.75*(1-cl_X-cl_Z)+cl_Z*0.5); // mean (total) molecular weights
   const double cl_mole      = 2.0/(1+cl_X);                             // mean electron molecular weights
   const double cl_moli      = 1.0/cl_X;                                 // mean proton molecular weights
   const double cl_moli_mole = cl_moli*cl_mole;                          // Assume the molecular weights are constant, mu_e*mu_i = 1.464
   const int    lv           = 0;

// get the numerical result
   real fluid[NCOMP_TOTAL];
   double Temp_nume     = 0.0;
   double Temp_nume_tmp = 0.0;
   double Tcool_nume    = 0.0;
   int    count         = 0;

   for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
   for (int k=0; k<PS1; k++) {
   for (int j=0; j<PS1; j++) {
   for (int i=0; i<PS1; i++) {
      for (int v=0; v<NCOMP_TOTAL; v++)   fluid[v] = amr->patch[ amr->FluSg[lv] ][lv][PID]->fluid[v][k][j][i];
      Temp_nume_tmp = Hydro_Con2Temp( fluid[0], fluid[1], fluid[2], fluid[3], fluid[4], fluid+NCOMP_FLUID,
                                      true, MIN_TEMP, PassiveFloorMask, 0.0, EoS_DensEint2Temp_CPUPtr, EoS_GuessHTilde_CPUPtr,
                                      EoS_HTilde2Temp_CPUPtr, EoS_AuxArray_Flt, EoS_AuxArray_Int, h_EoS_Table );
      Tcool_nume += 1.0/(GAMMA-1)*(Const_kB*cl_moli_mole*Temp_nume_tmp)/(fluid[0]*UNIT_D/MU_NORM*cl_mol*3.2217e-27*sqrt(Temp_nume_tmp))/Const_Myr;
      Temp_nume += Temp_nume_tmp;
      count += 1;
   }}} // i, j, k, PID
   Temp_nume  /= count;
   Tcool_nume /= count;

// compute the analytical solution for single branch cooling function
   double Temp_anal, Tcool_anal;
   if ( sqrt(EC_Temp) >= 3.2217e-27/2.0*(GAMMA-1)*EC_Dens*cl_mol*cl_mol/cl_mole/cl_moli/Const_kB*Time[0]*UNIT_T )
   {
      Temp_anal = pow(sqrt(EC_Temp) - 3.2217e-27/2.0*(GAMMA-1)*EC_Dens*cl_mol*cl_mol/cl_mole/cl_moli/Const_kB*Time[0]*UNIT_T, 2.0);
      if ( Temp_anal < MIN_TEMP )   Temp_anal = MIN_TEMP;
   }
   else   Temp_anal = MIN_TEMP;

   Tcool_anal = 1.0/(GAMMA-1)*EC_Dens*Const_kB*Temp_anal/((EC_Dens*cl_mol/cl_mole)*(EC_Dens*cl_mol/cl_moli)*3.2217e-27*sqrt(Temp_anal))/Const_Myr;

// record
   if ( MPI_Rank == 0 )
   {
      FILE *File_User = fopen( FileName, "a" );
      fprintf( File_User, "%14.7e%10d ", Time[0]*UNIT_T/Const_Myr, DumpID );
      fprintf( File_User, "%14.7e %14.7e %14.7e %14.7e %14.7e", Temp_nume, Temp_anal, (Temp_nume-Temp_anal)/Temp_anal, Tcool_nume, Tcool_anal );
      fprintf( File_User, "\n" );
      fclose( File_User );
   }

} // FUNCTION : Output_ExactCooling
#endif // #if ( MODEL == HYDRO )



//-------------------------------------------------------------------------------------------------------
// Function    :  Init_TestProb_Hydro_ExactCooling
// Description :  Test problem initializer
//
// Note        :  None
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Init_TestProb_Hydro_ExactCooling()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );


// validate the compilation flags and runtime parameters
   Validate();


#  if ( MODEL == HYDRO )
// set the problem-specific runtime parameters
   SetParameter();


   Init_Function_User_Ptr    = SetGridIC;
   Output_User_Ptr           = Output_ExactCooling;
#  ifdef SUPPORT_HDF5
   Output_HDF5_InputTest_Ptr = LoadInputTestProb;
#  endif
#  endif

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : Init_TestProb_Hydro_ExactCooling
