// #include "GAMER.h"


// // problem-specific global variables
// // =======================================================================================
// // static double  AGORA_DiskScaleLength;           // disk scale length
// // static double  AGORA_DiskScaleHeight;           // disk scale height
// // static double  AGORA_DiskTotalMass;             // disk total mass (gas + stars)
// static double  AGORA_DiskGasTemp;               // disk gas temperature
// static double  AGORA_Dens;           // halo atomic hydrogen number density (halo_gas_mass_density / atomic_hydrogen_mass)
//                                                 // --> necessary if one wants to enable metal_cooling in Grackle
// static double  AGORA_DiskMetalMassFrac;         // disk metal mass fraction (disk_metal_mass / disk_gas_mass)
// // =======================================================================================

// // problem-specific function prototypes
// double Mis_GetTimeStep_Dust( const int lv, const double dTime_dt )
// {

// // put your favorite time-step criteria here
// // ##########################################################################################################
//    double dt_user = HUGE_NUMBER;

// // Example 1 : set upper limit for the time interval to advance solution (per sub-step)
//    /*
//    dt_user = 1.0e-2;
//    */
// // Example 2 : set upper limit for the time interval to update the scale factor in COMOVING
//    /*
//    double dTime_user = 1.0e-5;
//    dt_user = dTime_user / dTime_dt;
//    */
// // ##########################################################################################################

//    double dTime_user = 2.0e0;

//    dt_user = dTime_user;
//    Aux_Message( stdout, "  dTime_dt     = %e \n",    dTime_dt  );
//    return dt_user;

// } // FUNCTION : Mis_GetTimeStep_Dust

// //-------------------------------------------------------------------------------------------------------
// // Function    :  Validate
// // Description :  Validate the compilation flags and runtime parameters for this test problem
// //
// // Note        :  None
// //
// // Parameter   :  None
// //
// // Return      :  None
// //-------------------------------------------------------------------------------------------------------
// void Validate()
// {

//    if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Validating test problem %d ...\n", TESTPROB_ID );

// #  if ( MODEL != HYDRO )
//    Aux_Error( ERROR_INFO, "MODEL != HYDRO !!\n" );
// #  endif

//    if ( !OPT__UNIT )
//       Aux_Error( ERROR_INFO, "OPT__UNIT must be enabled !!\n" );


// } // FUNCTION : Validate



// #if ( MODEL == HYDRO )
// //-------------------------------------------------------------------------------------------------------
// // Function    :  LoadInputTestProb
// // Description :  Read problem-specific runtime parameters from Input__TestProb and store them in HDF5 snapshots (Data_*)
// //
// // Note        :  1. Invoked by SetParameter() to read parameters
// //                2. Invoked by Output_DumpData_Total_HDF5() using the function pointer Output_HDF5_InputTest_Ptr to store parameters
// //                3. If there is no problem-specific runtime parameter to load, add at least one parameter
// //                   to prevent an empty structure in HDF5_Output_t
// //                   --> Example:
// //                       LOAD_PARA( load_mode, "TestProb_ID", &TESTPROB_ID, TESTPROB_ID, TESTPROB_ID, TESTPROB_ID );
// //
// // Parameter   :  load_mode      : Mode for loading parameters
// //                                 --> LOAD_READPARA    : Read parameters from Input__TestProb
// //                                     LOAD_HDF5_OUTPUT : Store parameters in HDF5 snapshots
// //                ReadPara       : Data structure for reading parameters (used with LOAD_READPARA)
// //                HDF5_InputTest : Data structure for storing parameters in HDF5 snapshots (used with LOAD_HDF5_OUTPUT)
// //
// // Return      :  None
// //-------------------------------------------------------------------------------------------------------
// void LoadInputTestProb( const LoadParaMode_t load_mode, ReadPara_t *ReadPara, HDF5_Output_t *HDF5_InputTest )
// {

// #  ifndef SUPPORT_HDF5
//    if ( load_mode == LOAD_HDF5_OUTPUT )   Aux_Error( ERROR_INFO, "please turn on SUPPORT_HDF5 in the Makefile for load_mode == LOAD_HDF5_OUTPUT !!\n" );
// #  endif

//    if ( load_mode == LOAD_READPARA     &&  ReadPara       == NULL )   Aux_Error( ERROR_INFO, "load_mode == LOAD_READPARA and ReadPara == NULL !!\n" );
//    if ( load_mode == LOAD_HDF5_OUTPUT  &&  HDF5_InputTest == NULL )   Aux_Error( ERROR_INFO, "load_mode == LOAD_HDF5_OUTPUT and HDF5_InputTest == NULL !!\n" );

// // add parameters in the following format:
// // --> note that VARIABLE, DEFAULT, MIN, and MAX must have the same data type
// // --> some handy constants (e.g., NoMin_int, Eps_float, ...) are defined in "include/ReadPara.h"
// // --> LOAD_PARA() is defined in "include/TestProb.h"
// // ********************************************************************************************************************************
// // LOAD_PARA( load_mode, "KEY_IN_THE_FILE",         &VARIABLE,                 DEFAULT,       MIN,              MAX               );
// // ********************************************************************************************************************************

//    // LOAD_PARA( load_mode, "AGORA_DiskScaleLength",   &AGORA_DiskScaleLength,   -1.0,           Eps_double,       NoMax_double      );
//    // LOAD_PARA( load_mode, "AGORA_DiskScaleHeight",   &AGORA_DiskScaleHeight,   -1.0,           Eps_double,       NoMax_double      );
//    // LOAD_PARA( load_mode, "AGORA_DiskTotalMass",     &AGORA_DiskTotalMass,     -1.0,           Eps_double,       NoMax_double      );
//    // LOAD_PARA( load_mode, "AGORA_DiskGasMassFrac",   &AGORA_DiskGasMassFrac,   -1.0,           Eps_double,       NoMax_double      );
//    LOAD_PARA( load_mode, "AGORA_DiskGasTemp",          &AGORA_DiskGasTemp,          -1.0,        Eps_double,    NoMax_double         );
//    LOAD_PARA( load_mode, "AGORA_Dens",                 &AGORA_Dens,                 -1.0,        Eps_double,    NoMax_double         );
//    LOAD_PARA( load_mode, "AGORA_DiskMetalMassFrac",    &AGORA_DiskMetalMassFrac,  -1.0,          Eps_double,      NoMax_double       );

// } // FUNCITON : LoadInputTestProb
// #  endif


// //-------------------------------------------------------------------------------------------------------
// // Function    :  SetParameter
// // Description :  Load and set the problem-specific runtime parameters
// //
// // Note        :  1. Filename is set to "Input__TestProb" by default
// //                2. Major tasks in this function:
// //                   (1) load the problem-specific runtime parameters
// //                   (2) set the problem-specific derived parameters
// //                   (3) reset other general-purpose parameters if necessary
// //                   (4) make a note of the problem-specific parameters
// //                3. Must call EoS_Init() before calling any other EoS routine
// //
// // Parameter   :  None
// //
// // Return      :  None
// //-------------------------------------------------------------------------------------------------------
// void SetParameter()
// {

//    if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Setting runtime parameters ...\n" );


// // (1) load the problem-specific runtime parameters
// // (1-1) read parameters from Input__TestProb
//    const char FileName[] = "Input__TestProb";
//    ReadPara_t *ReadPara  = new ReadPara_t;

//    LoadInputTestProb( LOAD_READPARA, ReadPara, NULL );

//    ReadPara->Read( FileName );

//    delete ReadPara;

// // convert to code units
//    // AGORA_DiskScaleLength *= Const_kpc  / UNIT_L;
//    // AGORA_DiskScaleHeight *= Const_kpc  / UNIT_L;
//    // AGORA_DiskTotalMass   *= Const_Msun / UNIT_M;
//    AGORA_Dens *= ( Const_amu / CUBE(Const_cm) ) / UNIT_D;


// // (3) reset other general-purpose parameters
// //     --> a helper macro PRINT_RESET_PARA is defined in Macro.h
//    const long   End_Step_Default = __INT_MAX__;
//    const double End_T_Default    =  500.0*Const_Myr/UNIT_T;

//    if ( END_STEP < 0 ) {
//       END_STEP = End_Step_Default;
//       PRINT_RESET_PARA( END_STEP, FORMAT_LONG, "" );
//    }

//    if ( END_T < 0.0 ) {
//       END_T = End_T_Default;
//       PRINT_RESET_PARA( END_T, FORMAT_REAL, "" );
//    }

// // (4) make a note
//    if ( MPI_Rank == 0 )
//    {
//       Aux_Message( stdout, "=============================================================================\n" );
//       Aux_Message( stdout, "  test problem ID      = %d\n",             TESTPROB_ID                               );
//       // Aux_Message( stdout, "  DiskScaleLength   = %13.7e kpc\n",     AGORA_DiskScaleLength * UNIT_L/Const_kpc  );
//       // Aux_Message( stdout, "  DiskScaleHeight   = %13.7e kpc\n",     AGORA_DiskScaleHeight * UNIT_L/Const_kpc  );
//       // Aux_Message( stdout, "  DiskTotalMass     = %13.7e Msun\n",    AGORA_DiskTotalMass   * UNIT_M/Const_Msun );
//       Aux_Message( stdout, "  AGORA_DiskMetalMassFrac     = %13.7e Msun\n",    AGORA_DiskMetalMassFrac  );
//       Aux_Message( stdout, "  DiskGasTemp          = %13.7e K\n",       AGORA_DiskGasTemp                         );
//       Aux_Message( stdout, "  AGORA_Dens           = %13.7e amu*cm^{-3}\n", AGORA_Dens * UNIT_D / ( Const_amu / CUBE(Const_cm) ) );
//       Aux_Message( stdout, "=============================================================================\n" );
//    }


//    if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Setting runtime parameters ... done\n" );

// } // FUNCTION : SetParameter



// //-------------------------------------------------------------------------------------------------------
// // Function    :  SetGridIC
// // Description :  Initialize the gas disk and halo for the AGORA isolated galaxy test
// //
// // Note        :  1. We do NOT truncate gas disk. Instead, we ensure pressure balance between gas disk and halo.
// //                2. This function will be invoked by multiple OpenMP threads when OPENMP is enabled
// //                   --> Please ensure that everything here is thread-safe
// //                3. Even when DUAL_ENERGY is adopted for HYDRO, one does NOT need to set the dual-energy variable here
// //                   --> It will be calculated automatically
// //
// // Parameter   :  fluid    : Fluid field to be initialized
// //                x/y/z    : Physical coordinates
// //                Time     : Physical time
// //                lv       : Target refinement level
// //                AuxArray : Auxiliary array
// //
// // Return      :  fluid
// //-------------------------------------------------------------------------------------------------------
// void SetGridIC( real fluid[], const double x, const double y, const double z, const double Time, const int lv, double AuxArray[] ){
// // check
// #  ifdef GAMER_DEBUG
//    if ( Idx_Metal == Idx_Undefined )
//       Aux_Error( ERROR_INFO, "Idx_Metal is undefined !!\n" );

//    if ( EoS_DensTemp2Pres_CPUPtr == NULL )
//       Aux_Error( ERROR_INFO, "EoS_DensTemp2Pres_CPUPtr == NULL !!\n" );
// #  endif
//    real MomX, MomY, MomZ, Pres, Eint, Etot;
//    real Metal=NULL_REAL;
//    real Passive[NCOMP_PASSIVE_USER];

//    // Aux_Message( stdout, "=============================================================================\n" );
//    // Aux_Message( stdout, "  NCOMP_PASSIVE_USER     = total number%d \n",    NCOMP_PASSIVE_USER  );
//    // Aux_Message( stdout, "  Idx_Metal              = No. %d \n",    Idx_Metal  );
//    // Aux_Message( stdout, "=============================================================================\n" );


//    MomX  = MomY = MomZ =  0.0;
//    Passive[NCOMP_PASSIVE_USER] = AGORA_Dens * AGORA_DiskMetalMassFrac;
//    Metal = Passive[NCOMP_PASSIVE_USER];

//    Pres = EoS_DensTemp2Pres_CPUPtr( AGORA_Dens, AGORA_DiskGasTemp, Passive, EoS_AuxArray_Flt, EoS_AuxArray_Int,
//                                            h_EoS_Table ); // assuming EoS requires no passive scalars
// // compute the total gas energy
//    Eint = EoS_DensPres2Eint_CPUPtr( AGORA_Dens, Pres, Passive, EoS_AuxArray_Flt,
//                                     EoS_AuxArray_Int, h_EoS_Table );   // assuming EoS requires no passive scalars
//    Etot = Hydro_ConEint2Etot( AGORA_Dens, MomX, MomY, MomZ, Eint, 0.0 );     // do NOT include magnetic energy here

// // set the output array
//    fluid[DENS] = AGORA_Dens;
//    fluid[MOMX] = MomX;
//    fluid[MOMY] = MomY;
//    fluid[MOMZ] = MomZ;
//    fluid[ENGY] = Etot;
//    fluid[Idx_Metal] = Metal;

// } // FUNCTION : SetGridIC



// //-------------------------------------------------------------------------------------------------------
// // Function    :  AddNewField_AGORA
// // Description :  Add the problem-specific grid fields
// //
// // Note        :  1. Ref: https://github.com/gamer-project/gamer/wiki/Adding-New-Simulations#v-add-problem-specific-grid-fields-and-particle-attributes
// //                2. Invoke AddField() for each of the problem-specific field:
// //                   --> Field label sent to AddField() will be used as the output name of the field
// //                   --> Field index returned by AddField() can be used to access the field data
// //                3. Pre-declared field indices are put in Field.h
// //
// // Parameter   :  None
// //
// // Return      :  None
// //-------------------------------------------------------------------------------------------------------
// void AddNewField_AGORA()
// {
// // add the metallicity field only if it has not been done
// // --> since Grackle may already add this field automatically when GRACKLE_METAL is enabled
// // --> also note that "Idx_Metal" has been predefined in Field.h
//    if ( Idx_Metal == Idx_Undefined )
//       Idx_Metal = AddField( "Metal", FIXUP_FLUX_YES, FIXUP_REST_YES, NORMALIZE_NO, INTERP_FRAC_YES );
// } // FUNCTION : AddNewField_AGORA


// //-------------------------------------------------------------------------------------------------------
// // Function    :  Init_TestProb_Hydro_Dust_Rho_Temp
// // Description :  Test problem initializer
// //
// // Note        :  None
// //
// // Parameter   :  None
// //
// // Return      :  None
// //-------------------------------------------------------------------------------------------------------
// void Init_TestProb_Hydro_Dust_Rho_Temp()
// {

//    if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );

// // validate the compilation flags and runtime parameters
//    Validate();

// #  if ( MODEL == HYDRO  )
//    // AddNewField_AGORA();
// // set the problem-specific runtime parameters
//    SetParameter();

// // set the function pointers of various problem-specific routines
//    Init_Function_User_Ptr      = SetGridIC;
//    Init_Field_User_Ptr         = AddNewField_AGORA;
//    Mis_GetTimeStep_User_Ptr    = Mis_GetTimeStep_Dust; // option: OPT__DT_USER;
// #  ifdef SUPPORT_HDF5
//    Output_HDF5_InputTest_Ptr   = LoadInputTestProb;
// #  endif
// #  endif // if ( MODEL == HYDRO )

//    if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );
// } // FUNCTION : Init_TestProb_Hydro_Dust_Rho_Temp
