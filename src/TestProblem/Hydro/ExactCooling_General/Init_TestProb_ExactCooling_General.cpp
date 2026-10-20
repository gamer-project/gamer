#include "GAMER.h"


static void Output_ExactCooling_General();
void Cooling_efficiency_table( const double Z, double TEF_Temp[], double TEF_lambda[] );

// problem-specific global variables
static double EC_Temp;
static double EC_Dens;

static double CoolingTimeToTmin_ExactCooling_General()
{
   const int TEF_N = SrcTerms.ECG_TEF_N;

   if ( TEF_N != 201 )
      Aux_Error( ERROR_INFO, "SRC_ECG_TEF_N must be 201 for this cooling table !!\n" );

   double *TEF_Temp   = new double [TEF_N];
   double *TEF_lambda = new double [TEF_N];

   Cooling_efficiency_table( SrcTerms.ECG_Z, TEF_Temp, TEF_lambda );

   const double Tmin = TEF_Temp[0];
   const double Tmax = TEF_Temp[TEF_N-1];

   if ( EC_Temp < Tmin || EC_Temp > Tmax )
      Aux_Error( ERROR_INFO, "EC_Temp = %13.7e K is outside the cooling-table range [%13.7e, %13.7e] K !!\n",
                 EC_Temp, Tmin, Tmax );

   if ( EC_Temp == Tmin )
   {
      delete [] TEF_Temp;
      delete [] TEF_lambda;
      return 0.0;
   }

   const double X_H  = 0.7;
   const double X_Z  = 0.018 * SrcTerms.ECG_Z;
   const double mu_e = 2.0 / (1.0 + X_H);
   const double mu_H = 1.0 / X_H;
   const double mu   = 1.0 / ( 2.0*X_H + 0.75*(1.0-X_H-X_Z) + 0.5*X_Z );
   const double rho  = EC_Dens * Const_mp * mu;
   const double n    = EC_Dens;
   const double n_e  = rho / (mu_e * Const_mp);
   const double n_H  = rho / (mu_H * Const_mp);
   const double CoolingCoeff = (GAMMA - 1.0) * n_e * n_H / (n * Const_kB);

   double t_Tmin = 0.0;
   for (int i=0; i<TEF_N-1 && EC_Temp>TEF_Temp[i]; i++)
   {
      const double Tk = TEF_Temp[i];
      const double Tupper = fmin( EC_Temp, TEF_Temp[i+1] );
      const double Lk = TEF_lambda[i];
      const double alpha = log10( TEF_lambda[i+1] / Lk ) / log10( TEF_Temp[i+1] / Tk );

      if ( alpha == 1.0 )
         t_Tmin += Tk / Lk * log( Tupper / Tk ) / CoolingCoeff;
      else
         t_Tmin += Tk / Lk * ( pow( Tupper / Tk, 1.0-alpha ) - 1.0 ) / (1.0-alpha) / CoolingCoeff;

      if ( EC_Temp <= TEF_Temp[i+1] )
         break;
   }

   delete [] TEF_Temp;
   delete [] TEF_lambda;

   return t_Tmin;
}

static double CoolingEfficiency_ExactCooling_General( const double Temp )
{
   const int TEF_N = SrcTerms.ECG_TEF_N;

   if ( Temp <= h_SrcECG_TEF_Temp[0] )
      return h_SrcECG_TEF_lambda[0];

   if ( Temp >= h_SrcECG_TEF_Temp[TEF_N-1] )
      return h_SrcECG_TEF_lambda[TEF_N-1];

   for (int i=TEF_N-2; i>=0; i--)
   {
      if ( Temp >= h_SrcECG_TEF_Temp[i] )
         return h_SrcECG_TEF_lambda[i] * pow( Temp / h_SrcECG_TEF_Temp[i], h_SrcECG_TEF_alpha[i] );
   }

   Aux_Error( ERROR_INFO, "failed to find the cooling-efficiency interval for Temp = %13.7e K !!\n", Temp );
   return NULL_REAL;
}

static double CoolingTEF_ExactCooling_General( const double Temp )
{
   const int    TEF_N  = SrcTerms.ECG_TEF_N;
   const double TEF_TN = h_SrcECG_TEF_Temp[TEF_N-1];
   const double TEF_LN = h_SrcECG_TEF_lambda[TEF_N-1];

   if ( Temp <= h_SrcECG_TEF_Temp[0] )
      return h_SrcECG_TEF_Y[0];

   if ( Temp >= TEF_TN )
      return h_SrcECG_TEF_Y[TEF_N-1];

   for (int i=TEF_N-2; i>=0; i--)
   {
      if ( Temp >= h_SrcECG_TEF_Temp[i] )
      {
         const double Tk = h_SrcECG_TEF_Temp[i];
         const double Lk = h_SrcECG_TEF_lambda[i];
         const double ak = h_SrcECG_TEF_alpha[i];

         if ( ak == 1.0 )
            return h_SrcECG_TEF_Y[i] + TEF_LN / Lk * Tk / TEF_TN * log( Tk / Temp );
         else
            return h_SrcECG_TEF_Y[i] + 1.0 / (1.0-ak) * TEF_LN / Lk * Tk / TEF_TN *
                   ( 1.0 - pow( Tk / Temp, ak-1.0 ) );
      }
   }

   Aux_Error( ERROR_INFO, "failed to find the TEF interval for Temp = %13.7e K !!\n", Temp );
   return NULL_REAL;
}

static double CoolingTEFInverse_ExactCooling_General( const double Y )
{
   const int    TEF_N  = SrcTerms.ECG_TEF_N;
   const double TEF_TN = h_SrcECG_TEF_Temp[TEF_N-1];
   const double TEF_LN = h_SrcECG_TEF_lambda[TEF_N-1];

   if ( Y >= h_SrcECG_TEF_Y[0] )
      return h_SrcECG_TEF_Temp[0];

   if ( Y <= h_SrcECG_TEF_Y[TEF_N-1] )
      return TEF_TN;

   for (int i=1; i<TEF_N; i++)
   {
      if ( h_SrcECG_TEF_Y[i] <= Y && Y < h_SrcECG_TEF_Y[i-1] )
      {
         const int    k  = i - 1;
         const double Tk = h_SrcECG_TEF_Temp[k];
         const double Lk = h_SrcECG_TEF_lambda[k];
         const double ak = h_SrcECG_TEF_alpha[k];
         const double Yk = h_SrcECG_TEF_Y[k];

         if ( ak == 1.0 )
            return Tk * exp( -Lk / TEF_LN * TEF_TN / Tk * (Y-Yk) );
         else
            return Tk * pow( 1.0 - (1.0-ak) * Lk / TEF_LN * TEF_TN / Tk * (Y-Yk), 1.0/(1.0-ak) );
      }
   }

   Aux_Error( ERROR_INFO, "failed to find the inverse TEF interval for Y = %13.7e !!\n", Y );
   return NULL_REAL;
}

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
//                4. For a new run, END_T is set automatically to the exact cooling time from EC_Temp to Tmin
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
   const double End_T_Default   = 8.0*Const_Myr/UNIT_T;

   if ( END_STEP < 0 ) {
      END_STEP = End_Step_Default;
      PRINT_RESET_PARA( END_STEP, FORMAT_LONG, "" );
   }

   if ( OPT__INIT != INIT_BY_RESTART )
   {
      END_T = CoolingTimeToTmin_ExactCooling_General() / UNIT_T;
      PRINT_RESET_PARA( END_T, FORMAT_REAL, " (automatically set to t_Tmin)" );
   }
   else if ( END_T < 0.0 ) {
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
   const double cl_Z   = 0.018 * SrcTerms.ECG_Z;                   // mass-fraction of metal
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
         fprintf( File_User, "# t_Tmin     : time to cool from the initial temperature to Tmin [Myr]\n" );
         fprintf( File_User, "# Temp_nume  : temperature updated by the exact-cooling solver [K]\n" );
         fprintf( File_User, "# Temp_anal  : temperature obtained from the analytical integration [K]\n" );
         fprintf( File_User, "# Err        : relative error defined as (Temp_nume-Temp_anal)/Temp_anal\n" );
         fprintf( File_User, "# Tcool_nume : cooling time computed from the output density and output temperature Temp_nume [Myr]\n" );
         fprintf( File_User, "# Tcool_anal : cooling time computed from the input density and analytical temperature Temp_anal [Myr]\n" );
         fprintf( File_User, "# Lambda_nume: cooling efficiency at Temp_nume [erg cm^3/s]\n");
         fprintf( File_User, "# Lambda_anal: cooling efficiency at Temp_anal [erg cm^3/s]\n");
         fprintf( File_User, "# ================================================================================================\n" );
         fprintf( File_User, "#%14s ", "t_Tmin" );
         fprintf( File_User, "#%13s%10s ",  "Time", "DumpID" );
         fprintf( File_User, "%14s %14s %14s %14s %14s %14s %14s", "Temp_nume", "Temp_anal", "Err", "Tcool_nume", "Tcool_anal", "Lambda_nume", "Lambda_anal" );
         fprintf( File_User, "\n" );
         fclose( File_User );
      }

      FirstTime = false;
   } // if ( FirstTime )

   // gas settings
   const double cl_X         = 0.7;                                      // mass-fraction of hydrogen
   const double cl_Z         = 0.018 * SrcTerms.ECG_Z;                   // mass-fraction of metal
   const double cl_mol       = 1.0/(2*cl_X+0.75*(1-cl_X-cl_Z)+cl_Z*0.5); // mean (total) molecular weights
   const double cl_mole      = 2.0/(1+cl_X);                             // mean electron molecular weights
   const double cl_moli      = 1.0/cl_X;                                 // mean proton molecular weights
   const int    lv           = 0;

   // get the numerical result
   real   fluid[NCOMP_TOTAL];
   double Temp_nume     = 0.0;
   double Temp_nume_tmp = 0.0;
   double Tcool_nume    = 0.0;
   double Lambda_nume   = 0.0;
   double Lambda_cell   = 0.0;
   int    count         = 0;

   for (int k=1; k<PS1; k++) {
   for (int j=1; j<PS1; j++) {
   for (int i=1; i<PS1; i++) {
      for (int v=0; v<NCOMP_TOTAL; v++)   fluid[v] = amr->patch[ amr->FluSg[lv] ][lv][0]->fluid[v][k][j][i];
      Temp_nume_tmp = (real) Hydro_Con2Temp( fluid[0], fluid[1], fluid[2], fluid[3], fluid[4], fluid+NCOMP_FLUID,
                                             true, MIN_TEMP, 0, 0.0,
                                             EoS_DensEint2Temp_CPUPtr, EoS_GuessHTilde_CPUPtr, EoS_HTilde2Temp_CPUPtr,
                                             EoS_AuxArray_Flt, EoS_AuxArray_Int, h_EoS_Table );
                                     
      
      Lambda_cell = CoolingEfficiency_ExactCooling_General( Temp_nume_tmp );
      Tcool_nume += Const_kB * (fluid[0] * UNIT_D / Const_mp / cl_mol) * Temp_nume_tmp /
                    ( (fluid[0] * UNIT_D / Const_mp / cl_mole) * (fluid[0] * UNIT_D / Const_mp / cl_moli) *
                      (GAMMA - 1.0) * Lambda_cell ) / Const_Myr;

      Temp_nume   += Temp_nume_tmp;
      Lambda_nume += Lambda_cell;
      count       += 1;
   }}} // i, j, k
   Temp_nume  /= count;
   Tcool_nume /= count;
   Lambda_nume/= count;

// compute the analytical solution using all 200 power-law intervals
   const double rho = EC_Dens * Const_mp * cl_mol;   // g cm^-3
   const double n   = EC_Dens;                       // cm^-3
   const double n_e = rho / (cl_mole * Const_mp);    // cm^-3
   const double n_H = rho / (cl_moli * Const_mp);    // cm^-3

   const int    TEF_N        = SrcTerms.ECG_TEF_N;
   const double TEF_TN       = h_SrcECG_TEF_Temp[TEF_N-1];
   const double TEF_LN       = h_SrcECG_TEF_lambda[TEF_N-1];
   const double CoolingRate  = (GAMMA - 1.0) * n_e * n_H * TEF_LN / ( n * Const_kB * TEF_TN );
   const double Y_initial    = CoolingTEF_ExactCooling_General( EC_Temp );
   const double Y_Tmin       = h_SrcECG_TEF_Y[0];
   const double t_Tmin       = (Y_Tmin - Y_initial) / CoolingRate;
   const double t_phys       = Time[0] * UNIT_T;
   const double Y_anal       = Y_initial + CoolingRate * t_phys;
   const double Temp_anal    = CoolingTEFInverse_ExactCooling_General( Y_anal );
   const double Lambda_anal  = CoolingEfficiency_ExactCooling_General( Temp_anal );
   const double Tcool_anal   = Const_kB * n * Temp_anal / ( n_e * n_H * (GAMMA - 1.0) * Lambda_anal ) / Const_Myr;

// record
   if ( MPI_Rank == 0 )
   {
      FILE *File_User = fopen( FileName, "a" );
      fprintf( File_User, "%14.7e ", t_Tmin / Const_Myr );
      fprintf( File_User, "%14.7e%10d ", Time[0]*UNIT_T/Const_Myr, DumpID );
      fprintf( File_User, "%14.7e %14.7e %14.7e %14.7e %14.7e %14.7e %14.7e", Temp_nume, Temp_anal, (Temp_nume-Temp_anal)/Temp_anal, Tcool_nume, Tcool_anal, Lambda_nume, Lambda_anal );
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
