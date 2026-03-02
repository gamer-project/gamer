#include "GAMER.h"

// problem-specific global variables
grackle_field_data my_fields;
static gr_float *my_cooling_time;
static gr_float *my_temperature;
static double  Gas_Temp;               // gas temperature
static double  Gas_Dens;               // gas density 
static double  Gas_MetalMassFrac;      // disk metal mass fraction (disk_metal_mass / disk_gas_mass)


// problem-specific function prototypes
double Mis_GetTimeStep_Dust( const int lv, const double dTime_dt )
{
   int FluSg = amr->FluSg[0];
   double Dens = amr->patch[FluSg][0][0]->fluid[DENS][0][0][0];
   double Eint = amr->patch[FluSg][0][0]->fluid[ENGY][0][0][0];

   my_fields.density[0] = Dens;
   my_fields.internal_energy[0] = Eint / Dens;

   // metal
   my_fields.metal_density[0] = amr->patch[FluSg][0][0]->fluid[Idx_Metal][0][0][0] ;
   Che_Units.density_units        = UNIT_D;
   Che_Units.length_units         = UNIT_L;
   Che_Units.time_units           = UNIT_T;
   Che_Units.velocity_units       = UNIT_V;

   // Calculate cooling time.
   if (calculate_cooling_time(&Che_Units, &my_fields, my_cooling_time) == 0) {
     Aux_Error( ERROR_INFO, "Error in calculate_cooling_time.\n");
   }

   double dTime_user = (0.1 * fabs(my_cooling_time[0])) * 1.0;
   Aux_Message( stdout, "  cooling_time   = %.15E,",    fabs(my_cooling_time[0])  );
   Aux_Message( stdout, "  dTime_user   = %.15E \n",    dTime_user  );
   return dTime_user;

} // FUNCTION : Mis_GetTimeStep_Dust


//-------------------------------------------------------------------------------------------------------
// Function    :  End_GrackleDust
// Description :  Free memory before terminating the program
//
// Note        :  1. Linked to the function pointer "End_User_Ptr" to replace "End_User()"
//
// Parameter   :  None
//-------------------------------------------------------------------------------------------------------
void End_GrackleDust()
{
   delete [] my_fields.grid_dimension;
   delete [] my_fields.grid_start;
   delete [] my_fields.grid_end;

   my_fields.grid_dimension = NULL;
   my_fields.grid_start = NULL;
   my_fields.grid_end = NULL;

   delete [] my_fields.density         ;
   delete [] my_fields.internal_energy ;
   delete [] my_fields.metal_density   ;

   my_fields.density         = NULL;
   my_fields.internal_energy = NULL;
   my_fields.metal_density   = NULL;
} // FUNCTION : End_GrackleDust


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

   if ( !OPT__UNIT )
      Aux_Error( ERROR_INFO, "OPT__UNIT must be enabled !!\n" );
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
// ********************************************************************************************************************************
// LOAD_PARA( load_mode, "KEY_IN_THE_FILE",         &VARIABLE,                 DEFAULT,       MIN,              MAX               );
// ********************************************************************************************************************************
   LOAD_PARA( load_mode, "Gas_Temp",             &Gas_Temp,           -1.0,        Eps_double,    NoMax_double       );
   LOAD_PARA( load_mode, "Gas_Dens",             &Gas_Dens,           -1.0,        Eps_double,    NoMax_double       );
   LOAD_PARA( load_mode, "Gas_MetalMassFrac",    &Gas_MetalMassFrac,  -1.0,        Eps_double,    NoMax_double       );

} // FUNCITON : LoadInputTestProb
#  endif


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

// convert to code units
   Gas_Dens *= ( Const_amu / CUBE(Const_cm) ) / UNIT_D;


// (3) reset other general-purpose parameters
//     --> a helper macro PRINT_RESET_PARA is defined in Macro.h
   const long   End_Step_Default = __INT_MAX__;
   const double End_T_Default    =  500.0*Const_Myr/UNIT_T;

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
      Aux_Message( stdout, "  test problem ID   = %d\n",                  TESTPROB_ID                        );
      Aux_Message( stdout, "  Gas_MetalMassFrac = %13.7e Msun\n",         Gas_MetalMassFrac                  );
      Aux_Message( stdout, "  Gas_Temp          = %13.7e K\n",            Gas_Temp                           );
      Aux_Message( stdout, "  Gas_Dens_kpc      = %13.7e amu*cm^{-3}\n",  Gas_Dens * UNIT_D / ( Const_amu / CUBE(Const_cm) ) );
      Aux_Message( stdout, "=============================================================================\n" );
   }


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Setting runtime parameters ... done\n" );

} // FUNCTION : SetParameter



//-------------------------------------------------------------------------------------------------------
// Function    :  SetGridIC
// Description :  Initialize the gas disk and halo for the AGORA isolated galaxy test
//
// Note        :  1. We do NOT truncate gas disk. Instead, we ensure pressure balance between gas disk and halo.
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
void SetGridIC( real fluid[], const double x, const double y, const double z, const double Time, const int lv, double AuxArray[] ){
// check
#  ifdef GAMER_DEBUG
   if ( Idx_Metal == Idx_Undefined )
      Aux_Error( ERROR_INFO, "Idx_Metal is undefined !!\n" );
   if ( Idx_Dust == Idx_Undefined )
      Aux_Error( ERROR_INFO, "Idx_Dust is undefined !!\n" );
   if ( EoS_DensTemp2Pres_CPUPtr == NULL )
      Aux_Error( ERROR_INFO, "EoS_DensTemp2Pres_CPUPtr == NULL !!\n" );
#  endif

   real MomX, MomY, MomZ, Pres, Eint, Etot;
   real Metal=NULL_REAL, Dust=NULL_REAL;
   real Passive[NCOMP_PASSIVE_USER];
   MomX  = MomY = MomZ =  0.0;
   Metal = Gas_Dens * Gas_MetalMassFrac;
   Dust = 0.01 * ( Const_amu / CUBE(Const_cm) ) / UNIT_D; // 0.1 amu
   Passive[0] = Metal;
   Passive[1] = Dust;

   Pres = EoS_DensTemp2Pres_CPUPtr( Gas_Dens, Gas_Temp, Passive, EoS_AuxArray_Flt, EoS_AuxArray_Int, h_EoS_Table ); // assuming EoS requires no passive scalars
   Eint = EoS_DensPres2Eint_CPUPtr( Gas_Dens, Pres, Passive, EoS_AuxArray_Flt, EoS_AuxArray_Int, h_EoS_Table );   // assuming EoS requires no passive scalars
   Etot = Hydro_ConEint2Etot( Gas_Dens, MomX, MomY, MomZ, Eint, 0.0 );    // compute the total gas energy  // do NOT include magnetic energy here

// set the output array
   fluid[DENS] = Gas_Dens;
   fluid[MOMX] = MomX;
   fluid[MOMY] = MomY;
   fluid[MOMZ] = MomZ;
   fluid[ENGY] = Etot;
   fluid[Idx_Metal] = Passive[0];
   fluid[Idx_Dust] = Passive[1];

} // FUNCTION : SetGridIC



//-------------------------------------------------------------------------------------------------------
// Function    :  AddNewField_DUST
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
void AddNewField_DUST()
{
// add the metallicity field only if it has not been done
// --> since Grackle may already add this field automatically when GRACKLE_METAL is enabled
// --> also note that "Idx_Metal" has been predefined in Field.h
   if ( Idx_Metal == Idx_Undefined )
      Idx_Metal = AddField( "Metal", FIXUP_FLUX_YES, FIXUP_REST_YES, NORMALIZE_NO, INTERP_FRAC_YES );
   if ( Idx_Dust == Idx_Undefined )
      Idx_Dust = AddField( "Dust", FIXUP_FLUX_YES, FIXUP_REST_YES, NORMALIZE_NO, INTERP_FRAC_YES );
} // FUNCTION : AddNewField_DUST


//-------------------------------------------------------------------------------------------------------
// Function    :  Init_TestProb_Hydro_Dust_Rho_Temp
// Description :  Test problem initializer
//
// Note        :  None
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Init_TestProb_Hydro_Dust_Rho_Temp()
{
   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );

   // Set grid dimension and size.
   // grid_start and grid_end are used to ignore ghost zones.
   int field_size = 1;
   my_fields.grid_rank = 3;
   my_fields.grid_dimension = new int[3];
   my_fields.grid_start = new int[3];
   my_fields.grid_end = new int[3];
   for (int i = 0;i < 3;i++) {
     my_fields.grid_dimension[i] = 1; // the active dimension not including ghost zones.
     my_fields.grid_start[i] = 0;
     my_fields.grid_end[i] = 0;
   }
   my_fields.grid_dimension[0] = field_size;
   my_fields.grid_end[0] = field_size - 1;
   my_fields.grid_dx = 0.0; // used only for H2 self-shielding approximation

   my_fields.density         = new gr_float[field_size];
   my_fields.internal_energy = new gr_float[field_size];
   my_fields.metal_density   = new gr_float[field_size];

   my_temperature  = new gr_float[field_size];
   my_cooling_time = new gr_float[field_size];


// validate the compilation flags and runtime parameters
   Validate();

#  if ( MODEL == HYDRO  )
   // AddNewField_DUST();
// set the problem-specific runtime parameters
   SetParameter();

// set the function pointers of various problem-specific routines
   Init_Function_User_Ptr      = SetGridIC;
   Init_Field_User_Ptr         = AddNewField_DUST;
   Mis_GetTimeStep_User_Ptr    = Mis_GetTimeStep_Dust; // option: OPT__DT_USER;
   End_User_Ptr                = End_GrackleDust;
#  ifdef SUPPORT_HDF5
   Output_HDF5_InputTest_Ptr   = LoadInputTestProb;
#  endif
#  endif // if ( MODEL == HYDRO )

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );
} // FUNCTION : Init_TestProb_Hydro_Dust_Rho_Temp
