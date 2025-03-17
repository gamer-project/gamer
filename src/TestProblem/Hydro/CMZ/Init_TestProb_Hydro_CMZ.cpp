#include "GAMER.h"



// problem-specific global variables
// =======================================================================================
       double BarredPot_V0;               // amplitude ( [km/s] ) of the lograthmic potential
       double BarredPot_q;                // axis ratio of the lograthmic potential
       double BarredPot_Rc;               // core radius [kpc] of the lograthmic potential
       double BarredPot_Omegabar;         // rotation speed [km/s/kpc] of the potential
       double BarredPot_fullBS;           // time [Gyr] for the barred potential to be fully turned-on
       double BarredPot_initT;            // initial gas temperature [K]
       double BarredPot_MeshCenter[3];    // central coordinates

// =======================================================================================

// problem-specific function prototypes
bool Flag_CMZ( const int i, const int j, const int k, const int lv, const int PID, const double *Threshold );
#ifdef PARTICLE
void Par_Init_ByFunction_BarredPot( const long NPar_ThisRank, const long NPar_AllRank,
                                    real_par *ParMass, real_par *ParPosX, real_par *ParPosY, real_par *ParPosZ,
                                    real_par *ParVelX, real_par *ParVelY, real_par *ParVelZ, real_par *ParTime,
                                    long_par *ParType, real_par *AllAttributeFlt[PAR_NATT_FLT_TOTAL],
                                    long_par *AllAttributeInt[PAR_NATT_INT_TOTAL] );
#endif
static void IsolatedBC( real Array[], const int ArraySize[], real fluid[], const int NVar_Flu,
                        const int GhostSize, const int idx[], const double pos[], const double Time,
                        const int lv, const int TFluVarIdxList[], double AuxArray[] );
//void Init_ExtAcc_BarredPot();
//void Init_ExtPot_BarredPot();
void Init_ExtPot_TabularP17();




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

   if ( !OPT__UNIT )
      Aux_Error( ERROR_INFO, "OPT__UNIT must be enabled !!\n" );


#  ifdef PARTICLE
   if ( OPT__INIT == INIT_BY_FUNCTION  &&  amr->Par->Init != PAR_INIT_BY_FUNCTION )
      Aux_Error( ERROR_INFO, "please set PAR_INIT = 1 (by FUNCTION) !!\n" );
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
// ***************************************************************************************************************************
// LOAD_PARA( load_mode, "KEY_IN_THE_FILE",      &VARIABLE,              DEFAULT,       MIN,              MAX               );
// ***************************************************************************************************************************
   LOAD_PARA( load_mode, "BarredPot_V0",         &BarredPot_V0,          200.,          Eps_double,       NoMax_double      );
   LOAD_PARA( load_mode, "BarredPot_q" ,         &BarredPot_q ,          0.9 ,          Eps_double,       NoMax_double      );
   LOAD_PARA( load_mode, "BarredPot_Rc",         &BarredPot_Rc,          0.05,          Eps_double,       NoMax_double      );
   LOAD_PARA( load_mode, "BarredPot_Omegabar",   &BarredPot_Omegabar,    40. ,          Eps_double,       NoMax_double      );
   LOAD_PARA( load_mode, "BarredPot_fullBS"  ,   &BarredPot_fullBS,      0.1 ,          Eps_double,       NoMax_double      );
   LOAD_PARA( load_mode, "BarredPot_initT"   ,   &BarredPot_initT,       1.0e4,         Eps_double,       NoMax_double      );

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

// convert to code units

   BarredPot_V0       *= Const_km             / UNIT_V;
   BarredPot_Rc       *= Const_kpc            / UNIT_L;
   BarredPot_Omegabar *= (Const_km/Const_kpc) / (1/UNIT_T);
   BarredPot_fullBS   *= Const_Gyr            / UNIT_T;


// (2) set the problem-specific derived parameters

   for (int d=0; d<3; d++)
       BarredPot_MeshCenter[d] = 0.5*amr->BoxSize[d];


// (3) reset other general-purpose parameters
//     --> a helper macro PRINT_RESET_PARA is defined in Macro.h
   const long   End_Step_Default = __INT_MAX__;
   const double End_T_Default    = 1.0*Const_Gyr/UNIT_T;

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
      Aux_Message( stdout, "==================================================================================\n"             );
      Aux_Message( stdout, "  test problem ID              = %d\n",              TESTPROB_ID                                  );
      Aux_Message( stdout, "  Log bar velocity             = %13.7e km/s\n",     BarredPot_V0*UNIT_V/Const_km                 );
      Aux_Message( stdout, "  Log bar axis ratio           = %13.7e\n",          BarredPot_q                                  );
      Aux_Message( stdout, "  Log bar pattern speed        = %13.7e km/s/kpc\n", BarredPot_Omegabar/UNIT_T*Const_kpc/Const_km );
      Aux_Message( stdout, "  Initial Gas disk temperature = %13.7e K\n",        BarredPot_initT                              );
      Aux_Message( stdout, "==================================================================================\n"             );
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

// check
#  ifdef GAMER_DEBUG
   if ( EoS_DensTemp2Pres_CPUPtr == NULL )
      Aux_Error( ERROR_INFO, "EoS_DensTemp2Pres_CPUPtr == NULL !!\n" );
#  endif


// Uniform gas disk, initially circular rotation

  const bool   CheckMinPres_Yes = true;
  const double dx               = x - BarredPot_MeshCenter[0];
  const double dy               = y - BarredPot_MeshCenter[1];
  const double dz               = z - BarredPot_MeshCenter[2];
  const double Rad              = sqrt( SQR(dx) + SQR(dy) );
  const double r                = sqrt( SQR(dx) + SQR(dy) + SQR(dz));

  double GasDens, GasPres, GasVel;
  double Eint,Etot,MomX,MomY,MomZ;
  double Metal=NULL_REAL;

  if (Rad <= 11.5 && fabs(dz) <= 0.2)
  {

     GasDens = 0.2735*exp(-Rad/4.8)*1/SQR(cosh(dz/0.13))*(Const_Msun/CUBE(Const_pc))/(UNIT_M/CUBE(UNIT_L));
     GasPres = EoS_DensTemp2Pres_CPUPtr( GasDens, BarredPot_initT, NULL, EoS_AuxArray_Flt, EoS_AuxArray_Int,
                                         h_EoS_Table ); // assuming EoS requires no passive scalars
//  GasVel  = BarredPot_V0*Rad/sqrt(SQR(Rad) + SQR(BarredPot_Rc));
     GasVel = -5.58683750e-05*pow(Rad,6) + 2.17357740e-03*pow(Rad,5) - 3.25132718e-02*pow(Rad,4)
               +2.32860976e-01*pow(Rad,3) - 8.14564481e-01*pow(Rad,2) + 1.35601708e+00*Rad + 1.06059808e+00;

     MomX  = -dy/r*GasVel*GasDens;
     MomY  = +dx/r*GasVel*GasDens;
     MomZ  = 0.0;

     Metal = GasDens*1.0e-2;

  }
  else
  {
     GasDens = 1.0e-6*(Const_Msun/CUBE(Const_pc))/(UNIT_M/CUBE(UNIT_L));
//   GasDens = 0.025 ;  // 1.0e-6
     GasPres = EoS_DensTemp2Pres_CPUPtr( GasDens, 100.*BarredPot_initT, NULL, EoS_AuxArray_Flt, EoS_AuxArray_Int,
                                         h_EoS_Table ); // assuming EoS requires no passive scalars
     MomX  = 0.0;
     MomY  = 0.0;
     MomZ  = 0.0;

     Metal = GasDens*1.0e-8;

  }


#  if ( EOS == EOS_ISOTHERMAL )

   GasPres = EoS_DensEint2Pres_CPUPtr( GasDens, 1.0,     NULL, EoS_AuxArray_Flt, EoS_AuxArray_Int, NULL );
   Eint    = EoS_DensPres2Eint_CPUPtr( GasDens, GasPres, NULL, EoS_AuxArray_Flt, EoS_AuxArray_Int, NULL );
   Etot    = Hydro_ConEint2Etot( GasDens, MomX, MomY, MomZ, Eint, 0.0 );

#  else

   Eint    = EoS_DensPres2Eint_CPUPtr( GasDens, GasPres, NULL, EoS_AuxArray_Flt,
                                       EoS_AuxArray_Int, h_EoS_Table );    // assuming EoS requires no passive scalars
   Etot    = Hydro_ConEint2Etot( GasDens, MomX, MomY, MomZ, Eint, 0.0 );   // do NOT include magnetic energy here
#  endif


  fluid[DENS] = GasDens ;
  fluid[MOMX] = MomX;
  fluid[MOMY] = MomY;
  fluid[MOMZ] = MomZ;
  fluid[ENGY] = Etot;
  fluid[Idx_Metal] = Metal;


} // FUNCTION : SetGridIC



//-------------------------------------------------------------------------------------------------------
// Function    :  IsolatedBC
// Description :  Isolated boundary condition for galaxies, only allow outflow velocities
//
// Note        :  1. Linked to the function pointer "BC_User_Ptr"
//
// Parameter   :  Array          : Array to store the prepared data including ghost zones
//                ArraySize      : Size of Array including the ghost zones on each side
//                fluid          : Fluid fields to be set
//                NVar_Flu       : Number of fluid variables to be prepared
//                GhostSize      : Number of ghost zones
//                idx            : Array indices
//                pos            : Physical coordinates
//                Time           : Physical time
//                lv             : Refinement level
//                TFluVarIdxList : List recording the target fluid variable indices ( = [0 ... NCOMP_TOTAL-1] )
//                AuxArray       : Auxiliary array
//
// Return      :  fluid
//-------------------------------------------------------------------------------------------------------
void IsolatedBC( real Array[], const int ArraySize[], real fluid[], const int NVar_Flu,
                 const int GhostSize, const int idx[], const double pos[], const double Time,
                 const int lv, const int TFluVarIdxList[], double AuxArray[] )
{

// simply call the IC function
   SetGridIC( fluid, pos[0], pos[1], pos[2], Time, lv, AuxArray );

} // FUNCTION : IsolatedBC



//-------------------------------------------------------------------------------------------------------
// Function    :  AddNewField_BarredPot
// Description :  Add the problem-specific fields
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
void AddNewField_BarredPot()
{
// add the metallicity field only if it has not been done
// --> since Grackle may already add this field automatically when GRACKLE_METAL is enabled
// --> also note that "Idx_Metal" has been predefined in Field.h
   if ( Idx_Metal == Idx_Undefined )
      Idx_Metal = AddField( "Metal", FIXUP_FLUX_YES, FIXUP_REST_YES, NORMALIZE_NO, INTERP_FRAC_YES );


} // FUNCTION : AddNewField_BarredPot

//-------------------------------------------------------------------------------------------------------
// Function    :  AddNewParticleAttribute_BarredPot
// Description :  Add the problem-specific particle attributes
//
// Note        :  1. Ref: https://github.com/gamer-project/gamer/wiki/Adding-New-Simulations#v-add-problem-specific-grid-fields-and-particle-attributes
//                2. Invoke AddParticleField() for each of the problem-specific particle attribute:
//                   --> Attribute label sent to AddParticleField() will be used as the output name of the attribute
//                   --> Attribute index returned by AddParticleField() can be used to access the particle attribute data
//                3. Pre-declared attribute indices are put in Field.h
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
//
#ifdef PARTICLE
void AddNewParticleAttribute_BarredPot()
{

// "Idx_ParMetalFrac" has been predefined in Field.h
   if ( Idx_ParMetalFrac == Idx_Undefined )
      Idx_ParMetalFrac = AddParticleAttributeFlt( "ParMetalFrac" );

} // FUNCTION : AddNewParticleAttribute_BarredPot
#endif
#endif // #if ( MODEL == HYDRO )


//-------------------------------------------------------------------------------------------------------
// Function    :  Init_TestProb_Hydro_BarredPot
// Description :  Test problem initializer
//
// Note        :  None
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Init_TestProb_Hydro_BarredPot()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );


// validate the compilation flags and runtime parameters
   Validate();


#  if ( MODEL == HYDRO )
// set the problem-specific runtime parameters
   SetParameter();

   Init_Function_User_Ptr      = SetGridIC;
   Init_Field_User_Ptr         = AddNewField_BarredPot;
   BC_User_Ptr                 = IsolatedBC;
   Flag_User_Ptr               = Flag_CMZ;
#  ifdef PARTICLE
   Par_Init_ByFunction_Ptr     = Par_Init_ByFunction_BarredPot;
   Par_Init_Attribute_User_Ptr = AddNewParticleAttribute_BarredPot;
#  endif
#  ifdef GRAVITY
//   Init_ExtAcc_Ptr              = Init_ExtAcc_BarredPot;
//   if ( OPT__EXT_POT == EXT_POT_FUNC )
   Init_ExtPot_Ptr             = Init_ExtPot_TabularP17;
#  endif // #ifdef GRAVITY
#  ifdef SUPPORT_HDF5
   Output_HDF5_InputTest_Ptr   = LoadInputTestProb;
#  endif
#  endif // #if ( MODEL == HYDRO )


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : Init_TestProb_Hydro_BarredPot
