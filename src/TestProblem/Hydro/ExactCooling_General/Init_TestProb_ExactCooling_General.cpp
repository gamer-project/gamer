#include "GAMER.h"


static void Output_ExactCooling_General();

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

#  ifndef EXACT_COOLING_GENERAL
   Aux_Error( ERROR_INFO, "EXACT_COOLING_GENERAL must be enabled !!\n" );
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
   LOAD_PARA( load_mode, "EC_Temp",           &EC_Temp,               10000000.0,    Eps_double,       NoMax_double      );
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
//                4. Please run the simulation first then change END_T in Input__Parameter to about "t1" in Record__CoolingErr and modify data dump dt
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
   const double End_T_Default    = 8.0*Const_Myr/UNIT_T;

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
      Aux_Message( stdout, "  test problem ID = %d\n",     TESTPROB_ID );
      Aux_Message( stdout, "  EC_Temp         = %13.7e\n", EC_Temp     );
      Aux_Message( stdout, "  EC_Dens         = %13.7e\n", EC_Dens     );
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

   // gas settings
   const double cl_X   = 0.7;                                      // mass-fraction of hydrogen
   const double cl_Z   = 0.018;                                    // mass-fraction of metal
   const double cl_mol = 1.0/(2*cl_X+0.75*(1-cl_X-cl_Z)+cl_Z*0.5); // mean (total) molecular weights

   double Dens, MomX, MomY, MomZ, Pres, Eint, Etot;
   double cl_dens = (EC_Dens*MU_NORM*cl_mol) / UNIT_D;   // convert the input number density into mass density rho
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
// Function    :  OutputExactCooling_General
// Description :  Output the temperature relative error in the general exact cooling problem
//
// Note        :  1. Enabled by the runtime option "OPT__OUTPUT_USER"
//                2. Construct the analytical solution corresponding to the cooling function
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Output_ExactCooling_General()
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
         fprintf( File_User, "# t3        : time to cool down from initial temperature to 10^7.4K [Myr]\n");
         fprintf( File_User, "# t2        : time to cool down from initial temperature to 10^5K [Myr]\n");
         fprintf( File_User, "# t3        : time to cool down from initial temperature to 10^4K [Myr]\n");
         fprintf( File_User, "# Temp_nume : temperature updated by the exact-cooling solver [K]\n" );
         fprintf( File_User, "# Temp_anal : temperature obtained from the analytical integration [K]\n" );
         fprintf( File_User, "# Err       : relative error defined as (Temp_nume-Temp_anal)/Temp_anal\n" );
         fprintf( File_User, "# Tcool_nume: cooling time computed from the output density and output temperature Temp_nume [Myr]\n" );
         fprintf( File_User, "# Tcool_anal: cooling time computed from the input density and analytical temperature Temp_anal [Myr]\n" );
         fprintf( File_User, "# Lambda    : cooling efficiency at the current temperature [erg cm^3/s]\n");
         fprintf( File_User, "# ================================================================================================\n" );
         fprintf( File_User, "#%14s %14s %14s ", "t3", "t2", "t1" );
         fprintf( File_User, "#%13s%10s ",  "Time", "DumpID" );
         fprintf( File_User, "%14s %14s %14s %14s %14s %14s", "Temp_nume", "Temp_anal", "Err", "Tcool_nume", "Tcool_anal", "Lambda" );
         fprintf( File_User, "\n" );
         fclose( File_User );
      }

      FirstTime = false;
   } // if ( FirstTime )

   // gas settings
   const double cl_X         = 0.7;                                      // mass-fraction of hydrogen
   const double cl_Z         = 0.018;                                    // metallicity (in Zsun)
   const double cl_mol       = 1.0/(2*cl_X+0.75*(1-cl_X-cl_Z)+cl_Z*0.5); // mean (total) molecular weights
   const double cl_mole      = 2.0/(1+cl_X);                             // mean electron molecular weights
   const double cl_moli      = 1.0/cl_X;                                 // mean proton molecular weights
   const double cl_moli_mole = cl_moli*cl_mole;                          // Assume the molecular weights are constant, mu_e*mu_i = 1.464
   const int    lv           = 0;

   // get the numerical result
   real   fluid[NCOMP_TOTAL];
   double Temp_nume     = 0.0;
   double Temp_nume_tmp = 0.0;
   double Tcool_nume    = 0.0;
   double Lambda_nume   = 0.0;
   double Lambda_cell   = 0.0;
   double Etot          = 0.0;
   int    count         = 0;

   // cooling function constant (3 power law cooling function)
   const double T1   = 1.0e4;             // K
   const double T2   = 1.0e5;             // K
   const double T3   = pow(10.0, 7.4);    // K
   const double TN   = 1.0e9;             // K

   const double a1   =  0.74;
   const double a2   = -0.70;
   const double a3   =  0.50;

   const double L1   = 9.93e-23;               // erg cm^3 s^-1
   const double L2   = 5.51e-22;               // erg cm^3 s^-1
   const double L3   = 1.15e-23;               // erg cm^3 s^-1
   const double LN   = L3 * pow( TN / T3, a3); // erg cm^3 s^-1

   for (int k=1; k<PS1; k++) {
   for (int j=1; j<PS1; j++) {
   for (int i=1; i<PS1; i++) {
      for (int v=0; v<NCOMP_TOTAL; v++)   fluid[v] = amr->patch[ amr->FluSg[lv] ][lv][0]->fluid[v][k][j][i];
      Temp_nume_tmp = (real) Hydro_Con2Temp( fluid[0], fluid[1], fluid[2], fluid[3], fluid[4], fluid+NCOMP_FLUID,
                                             true, MIN_TEMP, 0, 0.0,
                                             EoS_DensEint2Temp_CPUPtr, EoS_GuessHTilde_CPUPtr, EoS_HTilde2Temp_CPUPtr,
                                             EoS_AuxArray_Flt, EoS_AuxArray_Int, h_EoS_Table );
                                     
      
      if (Temp_nume_tmp >= T1 && Temp_nume_tmp < T2)
      {
         Lambda_cell = L1 * pow( Temp_nume_tmp / T1, a1); // erg cm^3/s
         Tcool_nume += Const_kB * (fluid[0] * UNIT_D / Const_mp / cl_mol) * Temp_nume_tmp/ ((fluid[0] * UNIT_D / Const_mp / cl_mole) * (fluid[0] * UNIT_D / Const_mp / cl_moli) * (GAMMA - 1.0) * Lambda_cell)/ Const_Myr;
      }
      else if (Temp_nume_tmp >= T2 && Temp_nume_tmp < T3)
      {
         Lambda_cell = L2 * pow( Temp_nume_tmp / T2, a2); // erg cm^3/s
         Tcool_nume += Const_kB * (fluid[0] * UNIT_D / Const_mp / cl_mol) * Temp_nume_tmp/ ((fluid[0] * UNIT_D / Const_mp / cl_mole) * (fluid[0] * UNIT_D / Const_mp / cl_moli) * (GAMMA - 1.0) * Lambda_cell)/ Const_Myr;
      }
      else if (Temp_nume_tmp >= T3)
      {
         Lambda_cell = L3 * pow( Temp_nume_tmp / T3, a3); // erg cm^3/s
         Tcool_nume += Const_kB * (fluid[0] * UNIT_D / Const_mp / cl_mol) * Temp_nume_tmp/ ((fluid[0] * UNIT_D / Const_mp / cl_mole) * (fluid[0] * UNIT_D / Const_mp / cl_moli) * (GAMMA - 1.0) * Lambda_cell)/ Const_Myr;
      }

      Temp_nume   += Temp_nume_tmp;
      Lambda_nume += Lambda_cell;
      Etot        += fluid[ENGY];
      count       += 1;
   }}} // i, j, k
   Temp_nume  /= count;
   Tcool_nume /= count;
   Lambda_nume/= count;
   Etot       /= count;

// compute the analytical solution for cooling function
// please change to your analytical solution if you change the number of power law cooling function
   const double rho = EC_Dens * Const_mp * cl_mol;   // g cm^-3
   const double n   = EC_Dens;                       // cm^-3
   const double n_e = rho / (cl_mole * Const_mp);    // cm^-3
   const double n_H = rho / (cl_moli * Const_mp);    // cm^-3

   double Temp_anal, Tcool_anal;

   const double A1 = (GAMMA - 1.0) * n_e * n_H * L1 * pow( T1, -a1 ) / (n * Const_kB);
   const double A2 = (GAMMA - 1.0) * n_e * n_H * L2 * pow( T2, -a2 ) / (n * Const_kB);
   const double A3 = (GAMMA - 1.0) * n_e * n_H * L3 * pow( T3, -a3 ) / (n * Const_kB);   

   // physical time
   const double t_phys = Time[0] * UNIT_T;
   
   // time to cool down from initial temperature to T3, T2, T1
   double t3 = 0.0;
   double t2 = 0.0;
   double t1 = 0.0;

   // (1) T3 <= T0 < TN
   if ( EC_Temp >= T3 && EC_Temp < TN )
   {
      t3 = 1.0 / ((1.0 - a3) * A3) * ( pow( EC_Temp, 1.0 - a3 ) - pow( T3, 1.0 - a3 ) );
      t2 = t3 + 1.0 / ((1.0 - a2) * A2) * ( pow( T3, 1.0 - a2 ) - pow( T2, 1.0 - a2 ) );
      t1 = t2 + 1.0 / ((1.0 - a1) * A1) * ( pow( T2, 1.0 - a1 ) - pow( T1, 1.0 - a1 ) );

      if ( t_phys <= t3 )
      {
         Temp_anal  = pow( pow( EC_Temp, 1.0 - a3 ) - (1.0 - a3) * A3 * t_phys, 1.0 / (1.0 - a3) );
         Tcool_anal = Const_kB * n * Temp_anal / (n_e * n_H * (GAMMA - 1.0) * L3 * pow( Temp_anal / T3, a3 ) ) / Const_Myr;
      }
      else if ( t_phys <= t2 && t_phys > t3 )
      {
         Temp_anal  = pow( pow( T3, 1.0 - a2 ) - (1.0 - a2) * A2 * (t_phys - t3), 1.0 / (1.0 - a2) );
         Tcool_anal = Const_kB * n * Temp_anal / (n_e * n_H * (GAMMA - 1.0) * L2 * pow( Temp_anal / T2, a2 ) ) / Const_Myr;
      }
      else if ( t_phys <= t1 && t_phys > t2 )
      {
         Temp_anal  = pow( pow( T2, 1.0 - a1 ) - (1.0 - a1) * A1 * (t_phys - t2), 1.0 / (1.0 - a1) );
         Tcool_anal = Const_kB * n * Temp_anal / (n_e * n_H * (GAMMA - 1.0) * L1 * pow( Temp_anal / T1, a1 ) ) / Const_Myr;
      }
      else
      {
         Temp_anal  = T1;
         Tcool_anal = Const_kB * n * Temp_anal / (n_e * n_H * (GAMMA - 1.0) * L1 * pow( Temp_anal / T1, a1 ) ) / Const_Myr;
      }

   }

   // (2) T2 <= T0 < T3
   else if ( EC_Temp >= T2 && EC_Temp < T3 )
   {
      t2 = 1.0 / ((1.0 - a2) * A2) * ( pow( EC_Temp, 1.0 - a2 ) - pow( T2, 1.0 - a2 ) );
      t1 = t2 + 1.0 / ((1.0 - a1) * A1) * ( pow( T2, 1.0 - a1 ) - pow( T1, 1.0 - a1 ) );

      if ( t_phys <= t2 )
      {
         Temp_anal  = pow( pow( EC_Temp, 1.0 - a2 ) - (1.0 - a2) * A2 * t_phys, 1.0 / (1.0 - a2) );
         Tcool_anal = Const_kB * n * Temp_anal / (n_e * n_H * (GAMMA - 1.0) * L2 * pow( Temp_anal / T2, a2 ) ) / Const_Myr;
      }
      else if ( t_phys <= t1 && t_phys > t2 )
      {
         Temp_anal  = pow( pow( T2, 1.0 - a1 ) - (1.0 - a1) * A1 * (t_phys - t2), 1.0 / (1.0 - a1) );
         Tcool_anal = Const_kB * n * Temp_anal / (n_e * n_H * (GAMMA - 1.0) * L1 * pow( Temp_anal / T1, a1 ) ) / Const_Myr;
      }
      else
      {
         Temp_anal  = T1;
         Tcool_anal = Const_kB * n * Temp_anal / (n_e * n_H * (GAMMA - 1.0) * L1 * pow( Temp_anal / T1, a1 ) ) / Const_Myr;
      }

   }

   // (3) T1 <= T0 < T2
   else if ( EC_Temp >= T1 && EC_Temp < T2 )
   {
      t1 = 1.0 / ((1.0 - a1) * A1) * ( pow( EC_Temp, 1.0 - a1 ) - pow( T1, 1.0 - a1 ) );

      if ( t_phys <= t1 )
      {
         Temp_anal  = pow( pow( EC_Temp, 1.0 - a1 ) - (1.0 - a1) * A1 * t_phys, 1.0 / (1.0 - a1) );
         Tcool_anal = Const_kB * n * Temp_anal / (n_e * n_H * (GAMMA - 1.0) * L1 * pow( Temp_anal / T1, a1 ) ) / Const_Myr;
      }
      else
      {
         Temp_anal  = T1;
         Tcool_anal = Const_kB * n * Temp_anal / (n_e * n_H * (GAMMA - 1.0) * L1 * pow( Temp_anal / T1, a1 ) ) / Const_Myr;
      }

   }
   

// record
   if ( MPI_Rank == 0 )
   {
      FILE *File_User = fopen( FileName, "a" );
      fprintf( File_User, "%14.7e%14.7e%14.7e ", t3 / Const_Myr, t2 / Const_Myr, t1 / Const_Myr );
      fprintf( File_User, "%14.7e%10d ", Time[0]*UNIT_T/Const_Myr, DumpID );
      fprintf( File_User, "%14.7e %14.7e %14.7e %14.7e %14.7e %14.7e", Temp_nume, Temp_anal, (Temp_nume-Temp_anal)/Temp_anal, Tcool_nume, Tcool_anal, Lambda_nume );
      fprintf( File_User, "\n" );
      fclose( File_User );
   }

} // FUNCTION : Output_ExactCooling_General
#endif // #if ( MODEL == HYDRO )


//-------------------------------------------------------------------------------------------------------
// Function    :  Init_TestProb_Hydro_ExactCooling_General
// Description :  Test problem initializer
//
// Note        :  None
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Init_TestProb_Hydro_ExactCooling_General()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );


// validate the compilation flags and runtime parameters
   Validate();


#  if ( MODEL == HYDRO )
// set the problem-specific runtime parameters
   SetParameter();


   Init_Function_User_Ptr    = SetGridIC;
   Output_User_Ptr           = Output_ExactCooling_General;
#  ifdef SUPPORT_HDF5
   Output_HDF5_InputTest_Ptr = LoadInputTestProb;
#  endif
#  endif


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : Init_TestProb_Hydro_ExactCooling_General
