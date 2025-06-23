#include "GAMER.h"



// problem-specific global variables
// =======================================================================================
static double              GrackleComoving_InitialTemperature;
       double              GrackleComoving_InitialMetallicity;
#ifdef SUPPORT_GRACKLE
       grackle_field_data  my_fields;
       gr_float           *my_temperature;
       gr_float           *my_gamma;
       gr_float           *my_cooling_time;
#endif // #ifdef SUPPORT_GRACKLE
// =======================================================================================

// problem-specific function prototypes
#ifdef SUPPORT_GRACKLE
void Aux_CoolingCurve( const double z_value, const double T_min, const double T_max, const double dT );
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

#  ifndef COMOVING
   Aux_Error( ERROR_INFO, "COMOVING must be enabled !!\n" );
#  endif

#  ifndef SUPPORT_GRACKLE
   Aux_Error( ERROR_INFO, "SUPPORT_GRACKLE must be enabled !!\n" );
#  endif

#  ifdef GRAVITY
   if ( OPT__BC_FLU[0] != BC_FLU_PERIODIC  ||  OPT__BC_POT != BC_POT_PERIODIC )
      Aux_Error( ERROR_INFO, "must adopt periodic BC for this test !!\n" );
#  endif

   if ( OPT__INIT != INIT_BY_FUNCTION  &&  OPT__INIT != INIT_BY_RESTART )
      Aux_Error( ERROR_INFO, "OPT__INIT != FUNCTION (1) or RESTART (2) for this test !!\n" );

#  ifdef PARTICLE
   Aux_Error( ERROR_INFO, "PARTICLE must be disabled !!\n" );
#  endif

   if ( !OPT__RECORD_USER )
      Aux_Error( ERROR_INFO, "must turn on OPT__RECORD_USER for this test !!\n" );

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Validating test problem %d ... done\n", TESTPROB_ID );

} // FUNCTION : Validate



#if ( MODEL == HYDRO  &&  defined SUPPORT_GRACKLE  &&  defined COMOVING )
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
// **************************************************************************************************************************************
// ReadPara->Add( "KEY_IN_THE_FILE",                    &VARIABLE,                             DEFAULT,     MIN,          MAX          );
// **************************************************************************************************************************************
   ReadPara->Add( "GrackleComoving_InitialTemperature", &GrackleComoving_InitialTemperature,  -1.0,         Eps_double,   NoMax_double );
   ReadPara->Add( "GrackleComoving_InitialMetallicity", &GrackleComoving_InitialMetallicity,   0.0,         0.0,          1.0          );

   ReadPara->Read( FileName );

   delete ReadPara;

// (1-2) set the default values

// (1-3) check the runtime parameters


// (2) set the problem-specific derived parameters


// (3) reset other general-purpose parameters
//     --> a helper macro PRINT_RESET_PARA is defined in Macro.h
   const long   End_Step_Default = __INT_MAX__;
   const double End_T_Default    = 1.0;

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
      Aux_Message( stdout, "  test problem ID    = %d\n",       TESTPROB_ID );
      Aux_Message( stdout, "  InitialTemperature = %13.7e K\n", GrackleComoving_InitialTemperature );
      Aux_Message( stdout, "  InitialMetallicity = %13.7e\n",   GrackleComoving_InitialMetallicity );
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
   double Dens, MomX, MomY, MomZ, Eint, Etot;

   const double mu                = 4 / (8 - 5 * (1 - grackle_data->HydrogenFractionByMass));      // fully ionized gas
   const double temperature_units = get_temperature_units(&Che_Units);                             // set temperature units

// convert GracleComoving_InitialTemperature to comoving internal energy assuming mu above
// --> the mean molecular weight can be different from MOLECULAR_WEIGHT
//     the temperature in this test problem can be different the temperature field output by OUT__OUTPUT_TEMP
   const double u = GrackleComoving_InitialTemperature / (mu * (GAMMA - 1.) * temperature_units) * SQR(Time);


   Dens = 1.0;
   MomX = 0.0;
   MomY = 0.0;
   MomZ = 0.0;
   Eint = Dens * u;
   Etot = Hydro_ConEint2Etot( Dens, MomX, MomY, MomZ, Eint, 0.0 );     // do NOT include magnetic energy here

   fluid[DENS] = Dens;
   fluid[MOMX] = MomX;
   fluid[MOMY] = MomY;
   fluid[MOMZ] = MomZ;
   fluid[ENGY] = Etot;

// initialize all chemical species
// --> assuming fully ionized gas
   if ( GRACKLE_PRIMORDIAL >= GRACKLE_PRI_CHE_NSPE6 ) {
   fluid[Idx_HI   ] = 0.0;
   fluid[Idx_HII  ] = grackle_data->HydrogenFractionByMass * (1.0 - GrackleComoving_InitialMetallicity) * Dens;
   fluid[Idx_HeI  ] = 0.0;
   fluid[Idx_HeII ] = 0.0;
   fluid[Idx_HeIII] = (1.0 - grackle_data->HydrogenFractionByMass) * (1.0 - GrackleComoving_InitialMetallicity) * Dens;
   fluid[Idx_e    ] = (fluid[Idx_HII] + fluid[Idx_HeII] / 4.0 + 2.0 * fluid[Idx_HeIII] / 4.0) * Const_me / Const_mp;
   }
// 9-species network
   if ( GRACKLE_PRIMORDIAL >= GRACKLE_PRI_CHE_NSPE9 ) {
   fluid[Idx_HM   ] = 0.0;
   fluid[Idx_H2I  ] = 0.0;
   fluid[Idx_H2II ] = 0.0;
   }
// 12-species network
   if ( GRACKLE_PRIMORDIAL >= GRACKLE_PRI_CHE_NSPE12 ) {
   fluid[Idx_DI   ] = 0.0;
   fluid[Idx_DII  ] = grackle_data->DeuteriumToHydrogenRatio * fluid[Idx_HII];
   fluid[Idx_HDI  ] = 0.0;
   }
   if ( GRACKLE_METAL ) {
   fluid[Idx_Metal] = GrackleComoving_InitialMetallicity * Dens;
   }
} // FUNCTION : SetGridIC



//-------------------------------------------------------------------------------------------------------
// Function    :  Aux_Record_GrackleComoving
// Description :  Record the thermo-chemical properties of the first cell
//
// Note        :  1. Invoked by main() using the function pointer "Aux_Record_User_Ptr",
//                   which must be set by a test problem initializer
//                2. Enabled by the runtime option "OPT__RECORD_USER"
//                3. This function will be called both during the program initialization and after each full update
//
// Parameter   :  None
//-------------------------------------------------------------------------------------------------------
void Aux_Record_GrackleComoving()
{

   const char FileName[] = "Record__GrackleComoving";
   static bool FirstTime = true;

   if ( FirstTime )
   {
//    header
      if ( MPI_Rank == 0 )
      {
         if ( Aux_CheckFileExist(FileName) )    Aux_Message( stderr, "WARNING : file \"%s\" already exists !!\n", FileName );

         FILE *File_User = fopen( FileName, "a" );
         fprintf( File_User, "# a_scale              : scale factor\n");
         fprintf( File_User, "# Step                 : current step\n");
         fprintf( File_User, "# dt [s]               : physical time-step size in sec\n");
         fprintf( File_User, "# n [1/cc]             : Hydrogen number density in physical cm^-3\n");
         fprintf( File_User, "# mu                   : mean molecular weight\n");
         fprintf( File_User, "# Temp [K]             : temperature in Kelvin\n");
         fprintf( File_User, "# Edens [erg cm^-3]    : internal energy density in physical erg cm^-3\n");
         fprintf( File_User, "# Lcool[erg cm^3 s^-1] : cooling rate in physical erg cm^-3 s^-1\n");
         fprintf( File_User, "#%13s%14s%3s%22s%22s%22s%22s%22s%22s",  "a_scale", "Step", "", "dt [s]", "n [1/cc]", "mu", "Temp [K]", "Edens [erg cm^-3]", "Lcool[erg cm^3 s^-1]" );
         fprintf( File_User, "\n" );

         fclose( File_User );
      }

      FirstTime = false;
   }

// user-specified info
   if ( MPI_Rank == 0 )
   {
      FILE  *File_User = fopen( FileName, "a" );
      int    FluSg     = amr->FluSg[0];
      double Dens      = amr->patch[FluSg][0][0]->fluid[DENS][0][0][0];
      double Eint      = amr->patch[FluSg][0][0]->fluid[ENGY][0][0][0]; // assume no magnetic and kinetic energy

//    use the dual-energy variable to calculate the internal energy if applicable
#     ifdef DUAL_ENERGY
      double Dual = amr->patch[FluSg][0][0]->fluid[DUAL][0][0][0];

#     if   ( DUAL_ENERGY == DE_ENPY )
      const bool CheckMinPres_No  = false;
      double     Pres             = Hydro_DensDual2Pres( Dens, Dual, EoS_AuxArray_Flt[1], CheckMinPres_No, NULL_REAL );
//    EOS_GAMMA does not involve passive scalars
      Eint  = EoS_DensPres2Eint_CPUPtr( Dens, Pres, NULL, EoS_AuxArray_Flt, EoS_AuxArray_Int, h_EoS_Table );
#     elif ( DUAL_ENERGY == DE_EINT )
#     error : DE_EINT is NOT supported yet !!
#     endif

#     endif // #ifdef DUAL_ENERGY

      const double MassRatio_pe = Const_mp / Const_me;

      my_fields.density        [0] = Dens;
      my_fields.internal_energy[0] = Eint / Dens / SQR(Time[0]);

      if ( GRACKLE_PRIMORDIAL >= GRACKLE_PRI_CHE_NSPE6 ) {
         my_fields.e_density    [0] = amr->patch[FluSg][0][0]->fluid[Idx_e    ][0][0][0] * MassRatio_pe;
         my_fields.HI_density   [0] = amr->patch[FluSg][0][0]->fluid[Idx_HI   ][0][0][0];
         my_fields.HII_density  [0] = amr->patch[FluSg][0][0]->fluid[Idx_HII  ][0][0][0];
         my_fields.HeI_density  [0] = amr->patch[FluSg][0][0]->fluid[Idx_HeI  ][0][0][0];
         my_fields.HeII_density [0] = amr->patch[FluSg][0][0]->fluid[Idx_HeII ][0][0][0];
         my_fields.HeIII_density[0] = amr->patch[FluSg][0][0]->fluid[Idx_HeIII][0][0][0];
      }

      if ( GRACKLE_PRIMORDIAL >= GRACKLE_PRI_CHE_NSPE9 ) {
         my_fields.HM_density   [0] = amr->patch[FluSg][0][0]->fluid[Idx_HM   ][0][0][0];
         my_fields.H2I_density  [0] = amr->patch[FluSg][0][0]->fluid[Idx_H2I  ][0][0][0];
         my_fields.H2II_density [0] = amr->patch[FluSg][0][0]->fluid[Idx_H2II ][0][0][0];
      }

      if ( GRACKLE_PRIMORDIAL >= GRACKLE_PRI_CHE_NSPE12 ) {
         my_fields.DI_density   [0] = amr->patch[FluSg][0][0]->fluid[Idx_DI   ][0][0][0];
         my_fields.DII_density  [0] = amr->patch[FluSg][0][0]->fluid[Idx_DII  ][0][0][0];
         my_fields.HDI_density  [0] = amr->patch[FluSg][0][0]->fluid[Idx_HDI  ][0][0][0];
      }

      if ( GRACKLE_METAL ) {
         my_fields.metal_density[0] = amr->patch[FluSg][0][0]->fluid[Idx_Metal][0][0][0];
      }

      Che_Units.comoving_coordinates = 1;
      Che_Units.density_units        = UNIT_D / CUBE(Time[0]);
      Che_Units.length_units         = UNIT_L * Time[0];
      Che_Units.time_units           = UNIT_T;
      Che_Units.velocity_units       = UNIT_V;
      Che_Units.a_units              = 1.0;
      Che_Units.a_value              = Time[0];

//    calculate cooling time
      if ( calculate_cooling_time( &Che_Units, &my_fields, my_cooling_time ) == 0 )
         Aux_Error( ERROR_INFO, "Error in calculate_cooling_time.\n" );

//    calculate temperature
      if ( calculate_temperature( &Che_Units, &my_fields, my_temperature ) == 0 )
         Aux_Error( ERROR_INFO, "Error in calculate_temperature.\n" );

//    calculate gamma
      if ( calculate_gamma( &Che_Units, &my_fields, my_gamma ) == 0 )
         Aux_Error( ERROR_INFO, "Error in calculate_gamma.\n" );

      const double dt_SubStep        = Mis_dTime2dt( Time[0], dTime_Base ) * SQR(Time[0]) * UNIT_T;                                 // physical time-step size
      const double temperature_units = get_temperature_units(&Che_Units);                                                           // grackle temperature unit
      const double mu                = my_temperature[0] / (my_fields.internal_energy[0] * (my_gamma[0] - 1.) * temperature_units); // mean molecular weight
      const double n                 = Dens / CUBE(Time[0]) * UNIT_D / mu / Const_mp;                                               // Hydrogen number density
      const double Edens             = Eint * UNIT_P / SQR(Time[0]) / CUBE(Time[0]);                                                // internal energy density
      const double Temp              = my_temperature[0];                                                                           // temperature
      const double Lcool             = Edens / fabs(my_cooling_time[0] * UNIT_T) / n / n;                                           // cooling rate obtained from grackle

      fprintf( File_User, "%14.7e%14ld%3s%22.7e%22.7e%22.7e%22.7e%22.7e%22.7e", Time[0], Step, "", dt_SubStep, n, mu, Temp, Edens, Lcool );
      fprintf( File_User, "\n" );

      fclose( File_User );
   } // if ( MPI_Rank == 0 )

} // FUNCTION : Aux_Record_GrackleComoving



//-------------------------------------------------------------------------------------------------------
// Function    :  Init_GrackleComoving
// Description :  Initialization of GrackleComoving
//
// Note        :  1. Invoked by Init_GAMER() using the function pointer "Init_User_Ptr",
//                   which must be set by a test problem initializer
//                2. Gravitational potential on grids and particle acceleration have not been computed
//                   at this stage
//                   --> Use Init_User_AfterPoisson_Ptr() instead if this information is required
//                   --> It's OK to modify mass density on grids and particle mass/position here
//                       --> But grid distribution won't change unless you manually call the corresponding
//                           grid refinement routines here
//                       --> Modifying particle position requires special attention in order to ensure that
//                           all particles still reside in leaf patches. Do this only if you know what
//                           you are doing.
//                3. To add new particles, remember to call Par_AddParticleAfterInit()
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Init_GrackleComoving()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );

// set grid dimension and size
// grid_start and grid_end are used to ignore ghost zones
   const int field_size = 1;
   my_fields.grid_rank = 3;
   my_fields.grid_dimension = new int [3];
   my_fields.grid_start     = new int [3];
   my_fields.grid_end       = new int [3];
   for (int i=0; i<3; i++) {
      my_fields.grid_dimension[i] = 1; // the active dimension not including ghost zones.
      my_fields.grid_start    [i] = 0;
      my_fields.grid_end      [i] = 0;
   }
   my_fields.grid_dimension[0] = field_size;
   my_fields.grid_end[0]       = field_size - 1;
   my_fields.grid_dx           = 0.0; // used only for H2 self-shielding approximation

   my_fields.density         = new gr_float [field_size];
   my_fields.internal_energy = new gr_float [field_size];
   if ( GRACKLE_PRIMORDIAL >= GRACKLE_PRI_CHE_NSPE6 ) {
   my_fields.HI_density      = new gr_float [field_size];
   my_fields.HII_density     = new gr_float [field_size];
   my_fields.HeI_density     = new gr_float [field_size];
   my_fields.HeII_density    = new gr_float [field_size];
   my_fields.HeIII_density   = new gr_float [field_size];
   my_fields.e_density       = new gr_float [field_size];
   }
   if ( GRACKLE_PRIMORDIAL >= GRACKLE_PRI_CHE_NSPE9 ) {
   my_fields.HM_density      = new gr_float [field_size];
   my_fields.H2I_density     = new gr_float [field_size];
   my_fields.H2II_density    = new gr_float [field_size];
   }
   if ( GRACKLE_PRIMORDIAL >= GRACKLE_PRI_CHE_NSPE12 ) {
   my_fields.DI_density      = new gr_float [field_size];
   my_fields.DII_density     = new gr_float [field_size];
   my_fields.HDI_density     = new gr_float [field_size];
   }
   if ( GRACKLE_METAL ) {
   my_fields.metal_density   = new gr_float [field_size];
   }

   my_temperature            = new gr_float [field_size];
   my_gamma                  = new gr_float [field_size];
   my_cooling_time           = new gr_float [field_size];


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

}



//-------------------------------------------------------------------------------------------------------
// Function    :  End_GrackleComoving
// Description :  Free memory before terminating the program
//
// Note        :  1. Linked to the function pointer "End_User_Ptr" to replace "End_User()"
//
// Parameter   :  None
//-------------------------------------------------------------------------------------------------------
void End_GrackleComoving()
{
   // generate cooling rate table for comparison
   if ( MPI_Rank == 0 )
   {
      Aux_Message( stdout, "==========================================\n");
      Aux_Message( stdout, "Generating cooling rate tables ...\n");
      for ( double z=floor(1.0/END_T-1.0); z<=(1.0/A_INIT-1.0); z+=1.0 )
         Aux_CoolingCurve( z, 4, 8, 0.02);
      Aux_Message( stdout, "Generating cooling rate tables ... done\n");
      Aux_Message( stdout, "==========================================\n");
   }

   delete [] my_fields.grid_dimension;     my_fields.grid_dimension  = NULL;
   delete [] my_fields.grid_start;         my_fields.grid_start      = NULL;
   delete [] my_fields.grid_end;           my_fields.grid_end        = NULL;

   delete [] my_fields.density;            my_fields.density         = NULL;
   delete [] my_fields.internal_energy;    my_fields.internal_energy = NULL;
   if ( GRACKLE_PRIMORDIAL >= GRACKLE_PRI_CHE_NSPE6 ) {
   delete [] my_fields.HI_density;         my_fields.HI_density      = NULL;
   delete [] my_fields.HII_density;        my_fields.HII_density     = NULL;
   delete [] my_fields.HeI_density;        my_fields.HeI_density     = NULL;
   delete [] my_fields.HeII_density;       my_fields.HeII_density    = NULL;
   delete [] my_fields.HeIII_density;      my_fields.HeIII_density   = NULL;
   delete [] my_fields.e_density;          my_fields.e_density       = NULL;
   }
   if ( GRACKLE_PRIMORDIAL >= GRACKLE_PRI_CHE_NSPE9 ) {
   delete [] my_fields.HM_density;         my_fields.HM_density      = NULL;
   delete [] my_fields.H2I_density;        my_fields.H2I_density     = NULL;
   delete [] my_fields.H2II_density;       my_fields.H2II_density    = NULL;
   }
   if ( GRACKLE_PRIMORDIAL >= GRACKLE_PRI_CHE_NSPE12 ) {
   delete [] my_fields.DI_density;         my_fields.DI_density      = NULL;
   delete [] my_fields.DII_density;        my_fields.DII_density     = NULL;
   delete [] my_fields.HDI_density;        my_fields.HDI_density     = NULL;
   }
   if ( GRACKLE_METAL ) {
   delete [] my_fields.metal_density;      my_fields.metal_density   = NULL;
   }

   delete [] my_temperature;               my_temperature            = NULL;
   delete [] my_gamma;                     my_gamma                  = NULL;
   delete [] my_cooling_time;              my_cooling_time           = NULL;
} // FUNCTION : End_GrackleComoving



//-------------------------------------------------------------------------------------------------------
// Function    :  Mis_GetTimeStep_GrackleComoving
// Description :  returns 0.01 * cooling time
//
// Note        :  1. This function should be applied to both physical and comoving coordinates and always
//                   return the evolution time-step (dt) actually used in various solvers
//                   --> Physical coordinates : dt = physical time interval
//                       Comoving coordinates : dt = delta(scale_factor) / ( Hubble_parameter*scale_factor^3 )
//                   --> We convert dt back to the physical time interval, which equals "delta(scale_factor)"
//                       in the comoving coordinates, in Mis_GetTimeStep()
//                2. Invoked by Mis_GetTimeStep() using the function pointer "Mis_GetTimeStep_User_Ptr",
//                   which must be set by a test problem initializer
//                3. Enabled by the runtime option "OPT__DT_USER"
//
// Parameter   :  lv       : Target refinement level
//                dTime_dt : dTime/dt (== 1.0 if COMOVING is off)
//
// Return      :  dt
//-------------------------------------------------------------------------------------------------------
double Mis_GetTimeStep_GrackleComoving( const int lv, const double dTime_dt )
{

   double dt_user_phy = HUGE_NUMBER;

   int    FluSg = amr->FluSg[0];
   double Dens  = amr->patch[FluSg][0][0]->fluid[DENS][0][0][0];
   double Eint  = amr->patch[FluSg][0][0]->fluid[ENGY][0][0][0]; // assume no magnetic and kinetic energy

// use the dual-energy variable to calculate the internal energy if applicable
#  ifdef DUAL_ENERGY
   double Dual = amr->patch[FluSg][0][0]->fluid[DUAL][0][0][0];

#  if   ( DUAL_ENERGY == DE_ENPY )
   const bool CheckMinPres_No  = false;
   double     Pres             = Hydro_DensDual2Pres( Dens, Dual, EoS_AuxArray_Flt[1], CheckMinPres_No, NULL_REAL );
// EOS_GAMMA does not involve passive scalars
   Eint  = EoS_DensPres2Eint_CPUPtr( Dens, Pres, NULL, EoS_AuxArray_Flt, EoS_AuxArray_Int, h_EoS_Table );
#  elif ( DUAL_ENERGY == DE_EINT )
#  error : DE_EINT is NOT supported yet !!
#  endif

#  endif // #ifdef DUAL_ENERGY

   const double MassRatio_pe = Const_mp / Const_me;

   my_fields.density        [0] = Dens;
   my_fields.internal_energy[0] = Eint / Dens / SQR(Time[lv]);

   if ( GRACKLE_PRIMORDIAL >= GRACKLE_PRI_CHE_NSPE6 ) {
      my_fields.e_density    [0] = amr->patch[FluSg][0][0]->fluid[Idx_e    ][0][0][0] * MassRatio_pe;
      my_fields.HI_density   [0] = amr->patch[FluSg][0][0]->fluid[Idx_HI   ][0][0][0];
      my_fields.HII_density  [0] = amr->patch[FluSg][0][0]->fluid[Idx_HII  ][0][0][0];
      my_fields.HeI_density  [0] = amr->patch[FluSg][0][0]->fluid[Idx_HeI  ][0][0][0];
      my_fields.HeII_density [0] = amr->patch[FluSg][0][0]->fluid[Idx_HeII ][0][0][0];
      my_fields.HeIII_density[0] = amr->patch[FluSg][0][0]->fluid[Idx_HeIII][0][0][0];
   }

   if ( GRACKLE_PRIMORDIAL >= GRACKLE_PRI_CHE_NSPE9 ) {
      my_fields.HM_density   [0] = amr->patch[FluSg][0][0]->fluid[Idx_HM   ][0][0][0];
      my_fields.H2I_density  [0] = amr->patch[FluSg][0][0]->fluid[Idx_H2I  ][0][0][0];
      my_fields.H2II_density [0] = amr->patch[FluSg][0][0]->fluid[Idx_H2II ][0][0][0];
   }

   if ( GRACKLE_PRIMORDIAL >= GRACKLE_PRI_CHE_NSPE12 ) {
      my_fields.DI_density   [0] = amr->patch[FluSg][0][0]->fluid[Idx_DI   ][0][0][0];
      my_fields.DII_density  [0] = amr->patch[FluSg][0][0]->fluid[Idx_DII  ][0][0][0];
      my_fields.HDI_density  [0] = amr->patch[FluSg][0][0]->fluid[Idx_HDI  ][0][0][0];
   }

   if ( GRACKLE_METAL ) {
      my_fields.metal_density[0] = amr->patch[FluSg][0][0]->fluid[Idx_Metal][0][0][0];
   }

   Che_Units.comoving_coordinates = 1;
   Che_Units.density_units        = UNIT_D / CUBE(Time[lv]);
   Che_Units.length_units         = UNIT_L * Time[lv];
   Che_Units.time_units           = UNIT_T;
   Che_Units.velocity_units       = UNIT_V;
   Che_Units.a_units              = 1.0;
   Che_Units.a_value              = Time[lv];


// calculate cooling time
   if ( calculate_cooling_time( &Che_Units, &my_fields, my_cooling_time ) == 0 )
      Aux_Error( ERROR_INFO, "Error in calculate_cooling_time.\n" );

   dt_user_phy = FMIN(dt_user_phy, 0.01 * fabs(my_cooling_time[0]));


// recalculate cooling time with 10% lower internal energy
// --> to avoid overestimating the time-step size when the cooling time gets shorter as the temperature decreases
   my_fields.internal_energy[0] *= 0.9;

   if ( calculate_cooling_time( &Che_Units, &my_fields, my_cooling_time ) == 0 )
      Aux_Error( ERROR_INFO, "Error in calculate_cooling_time.\n" );

   dt_user_phy = FMIN(dt_user_phy, 0.01 * fabs(my_cooling_time[0]));


// convert the time-step size to comoving coordinates
   const double dt_user = dt_user_phy / SQR(Time[lv]);


   return dt_user;

} // FUNCTION : Mis_GetTimeStep_GrackleComoving
#endif // #if ( MODEL == HYDRO  &&  defined SUPPORT_GRACKLE  &&  defined COMOVING )



//-------------------------------------------------------------------------------------------------------
// Function    :  Init_TestProb_Hydro_Grackle_Comoving
// Description :  Test problem initializer
//
// Note        :  None
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Init_TestProb_Hydro_Grackle_Comoving()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );


// validate the compilation flags and runtime parameters
   Validate();

#  if ( MODEL == HYDRO  &&  defined SUPPORT_GRACKLE  &&  defined COMOVING )
// set the problem-specific runtime parameters
   SetParameter();

// set the function pointers of various problem-specific routines
   Init_Function_User_Ptr   = SetGridIC;
   Aux_Record_User_Ptr      = Aux_Record_GrackleComoving;
   Init_User_Ptr            = Init_GrackleComoving;
   End_User_Ptr             = End_GrackleComoving;
   Mis_GetTimeStep_User_Ptr = Mis_GetTimeStep_GrackleComoving;
#  endif // #if ( MODEL == HYDRO  &&  defined SUPPORT_GRACKLE  &&  defined COMOVING )


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : Init_TestProb_Hydro_Grackle_Comoving
