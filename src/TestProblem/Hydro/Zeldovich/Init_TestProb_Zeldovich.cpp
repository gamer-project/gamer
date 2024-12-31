#include "GAMER.h"
#include <math.h>
#ifdef SUPPORT_GSL
#  include <gsl/gsl_errno.h>
#  include <gsl/gsl_roots.h>
#endif



// problem-specific global variables
// =======================================================================================
// =======================================================================================
static int    Gas_Par_Setup;     // 1=gas-only, 2=particle-only
static int    n_Pert_Wave_Len;   // "Pert_Wave_Len = BOX_SIZE_x/n_Pert_Wave_Len"
static int    NPar_X;            // 1D number of particles along the x-axis (axis of perturbation)
static int    NPar_YZ;           // 1D number of particles along the y/z-axis
static double Pert_Wave_Len;     // perturbation wavelength [Mpc/h]
static double k_Pert;            // perturbation wavenumber
static double z_Initial;         // redshift of simulation IC
static double z_Collapse;        // redshift at the onset of Zeldovich pancake collapse
static double GasTemp_Init;      // initial gas temperature in Kelvin
static double H0;                // Hubble parameter at redshift zero
static double rho_crit_i;        // comoving critical density = 3H0^2/8*pi*G
static double dhx;               // 1D inter-particle separation [Mpc/h]
// =======================================================================================


// problem-specific function prototypes
struct Zeldovich_coord_params
{
   double x_shift_input, z;
};
static double Zeldovich_coord_transformation( const double x_Lagrangian_Trial, void *params );
static void OutputError();
#ifdef SUPPORT_GSL
static double Zeldovich_coord_transformation_solver( const double x_shift_input, const double z );
static double Zeldovich_density_profile( const double x_Lagrangian, const double z );
static double Zeldovich_x_velocity_profile( const double x_Lagrangian, const double z );
#ifdef PARTICLE
void Par_Init_ByFunction_Zeldovich( const long NPar_ThisRank, const long NPar_AllRank,
                                    real_par *ParMass, real_par *ParPosX, real_par *ParPosY, real_par *ParPosZ,
                                    real_par *ParVelX, real_par *ParVelY, real_par *ParVelZ, real_par *ParTime,
                                    long_par *ParType, real_par *AllAttributeFlt[PAR_NATT_FLT_TOTAL],
                                    long_par *AllAttributeInt[PAR_NATT_INT_TOTAL] );
#endif // #ifdef PARTICLE
#endif // #ifdef SUPPORT_GSL




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

#  ifndef COMOVING
   Aux_Error( ERROR_INFO, "COMOVING must be enabled !!\n" );
#  endif

#  ifdef COMOVING
   if ( OMEGA_M0 != 1.0 )
      Aux_Error( ERROR_INFO, "OMEGA_M0 (%13.7e) != 1.0 (matter-only) for this test !!\n", OMEGA_M0 );
#  endif

#  ifdef GRAVITY
   if ( OPT__BC_FLU[0] != BC_FLU_PERIODIC  ||  OPT__BC_POT != BC_POT_PERIODIC )
      Aux_Error( ERROR_INFO, "must adopt periodic BC for this test !!\n" );
#  endif

   if ( OPT__INIT != INIT_BY_FUNCTION  &&  OPT__INIT != INIT_BY_RESTART )
      Aux_Error( ERROR_INFO, "OPT__INIT != FUNCTION (1) or RESTART (2) for this test !!\n" );

#  ifndef SUPPORT_GSL
   Aux_Error( ERROR_INFO, "SUPPORT_GSL must be enabled !!\n" );
#  endif

#  ifdef PARTICLE
   if ( amr->Par->Init != PAR_INIT_BY_FUNCTION  &&  amr->Par->Init != PAR_INIT_BY_RESTART )
      Aux_Error( ERROR_INFO, "PAR_INIT != FUNCTION (1) OR RESTART (2) for this test !!\n" );

   if ( NX0_TOT[1] != NX0_TOT[2] )
      Aux_Error( ERROR_INFO, "NX0_TOT[1] (%d) != NX0_TOT[2] (%d) !!\n", NX0_TOT[1], NX0_TOT[2] );
#  endif


// warnings
   if ( MPI_Rank == 0 )
   {
#     ifndef FLOAT8
      Aux_Message( stderr, "WARNING : it's recommended to enable FLOAT8 for gas-only setup to" );
      Aux_Message( stderr, " properly resolve the pressure and temperature of the cold flow\n" );
#     endif

#     ifndef DUAL_ENERGY
      Aux_Message( stderr, "WARNING : it's recommended to enable DUAL_ENERGY for this test\n" );
#     endif
   } // // if ( MPI_Rank == 0 )


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

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Setting runtime parameters ...\n" );


// (1) load the problem-specific runtime parameters
   const char FileName[] = "Input__TestProb";
   ReadPara_t *ReadPara  = new ReadPara_t;

// (1-1) add parameters in the following format:
// --> note that VARIABLE, DEFAULT, MIN, and MAX must have the same data type
// --> some handy constants (e.g., Useless_bool, Eps_double, NoMin_int, ...) are defined in "include/ReadPara.h"
// ********************************************************************************************************************************
// ReadPara->Add( "KEY_IN_THE_FILE",   &VARIABLE,               DEFAULT,      MIN,              MAX               );
// ********************************************************************************************************************************
   ReadPara->Add( "Gas_Par_Setup",     &Gas_Par_Setup,          1,            1,                2                 );
   ReadPara->Add( "n_Pert_Wave_Len",   &n_Pert_Wave_Len,        1,            1,                NoMax_int         );
   ReadPara->Add( "z_Collapse",        &z_Collapse,             1.0,          0.0,              NoMax_double      );
   ReadPara->Add( "GasTemp_Init",      &GasTemp_Init,           100.0,        0.0,              NoMax_double      );
   ReadPara->Add( "NPar_X",            &NPar_X,                 NX0_TOT[0],   0,                NoMax_int         );

   ReadPara->Read( FileName );

   delete ReadPara;


// (1-2) set the default values
   z_Initial = (1.0/Time[0]) - 1.0;


// (1-3) check the runtime parameters
   if ( z_Initial < z_Collapse )
      Aux_Error( ERROR_INFO, "z_Collapse (%13.7e) must be smaller than z_Initial (%13.7e) !!\n", z_Collapse, z_Initial );

   if      ( Gas_Par_Setup == 1 )   // for gas-only setup
   {
      if ( OPT__FREEZE_FLUID != 0 )
         Aux_Error( ERROR_INFO, "OPT__FREEZE_FLUID != 0 for Gas_Par_Setup == 1 !!\n");

      if ( OPT__OUTPUT_USER  != 1 )
         Aux_Error( ERROR_INFO, "OPT__OUTPUT_USER != 1 for Gas_Par_Setup == 1 !!\n");
   }

   else if ( Gas_Par_Setup == 2 )   // for particle-only setup
   {
      if ( OPT__FREEZE_FLUID != 1 )
         Aux_Error( ERROR_INFO, "OPT__FREEZE_FLUID != 1 for Gas_Par_Setup == 2 !!\n");

      if ( OPT__OUTPUT_USER  != 0  &&  MPI_Rank == 0 )
         Aux_Message( stderr, "WARNING : it's recommended to disable OPT__OUTPUT_USER for Gas_Par_Setup == 2\n" );
   }


// (2) set the problem-specific derived parameters
#  ifdef COMOVING
   H0            = HUBBLE0*(100.0*Const_km/(Const_s*Const_Mpc))*UNIT_T;
   Pert_Wave_Len = amr->BoxSize[0]/n_Pert_Wave_Len;
   rho_crit_i    = 1.0;            // [UNIT_D]; comoving critical density = 3H0^2/8*pi*G
   k_Pert        = 2.0*M_PI/Pert_Wave_Len;
   dhx           = amr->BoxSize[0] / NPar_X;
   NPar_YZ       = amr->BoxSize[1] / dhx;
#  endif


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

#  ifdef PARTICLE
   if      ( Gas_Par_Setup == 1 )   // overwrite the total number of particles in gas-only setup
      amr->Par->NPar_Active_AllRank = 0;

   else if ( Gas_Par_Setup == 2 )   // overwrite the total number of particles in particle-only setup
      amr->Par->NPar_Active_AllRank = (long)NPar_X*SQR((long)NPar_YZ);

   PRINT_RESET_PARA( amr->Par->NPar_Active_AllRank, FORMAT_LONG, "(PAR_NPAR in Input__Parameter)" );
#  endif


// (4) make a note
   if ( MPI_Rank == 0 )
   {
      Aux_Message( stdout, "=============================================================================\n" );
      Aux_Message( stdout, "  test problem ID                   = %d\n",     TESTPROB_ID     );
      Aux_Message( stdout, "  gas[1]/particle[2] setup          = %d\n",     Gas_Par_Setup   );
      Aux_Message( stdout, "  perturbation wavelength multiple  = %d\n",     n_Pert_Wave_Len );
      Aux_Message( stdout, "  Zeldovich collapse redshift       = %13.7e\n", z_Collapse      );
      Aux_Message( stdout, "  initial gas temperature in Kelvin = %13.7e\n", GasTemp_Init    );
      if ( Gas_Par_Setup == 2 )
      Aux_Message( stdout, "  particle number along the x-axis  = %d\n",     NPar_X          );
      Aux_Message( stdout, "=============================================================================\n" );
   }


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Setting runtime parameters ... done\n" );

} // FUNCTION : SetParameter



// compute conversion between Lagrangian (x_Lagrangian) and Eulerian/code coordinates (x_shift)
double Zeldovich_coord_transformation( const double x_Lagrangian_Trial, void *params )
{
   struct Zeldovich_coord_params *p = (struct Zeldovich_coord_params *) params;
   double x_shift_input = p->x_shift_input;
   double z             = p->z;

   return (x_Lagrangian_Trial - ((1.0+z_Collapse)/(1.0+z))*sin(k_Pert*x_Lagrangian_Trial)/k_Pert) - x_shift_input;
}



#ifdef SUPPORT_GSL
double Zeldovich_coord_transformation_solver( const double x_shift_input, const double z )
{
   int status;
   int iter = 0, max_iter = 100;
   const gsl_root_fsolver_type *T;
   gsl_root_fsolver *s;
   double x_Lagrangian_Trial = x_shift_input;
   double x_lower_bound = x_shift_input-Pert_Wave_Len/4.0, x_upper_bound = x_shift_input+Pert_Wave_Len/4.0;
   gsl_function F;
   struct Zeldovich_coord_params params = { x_shift_input, z };

   F.function = &Zeldovich_coord_transformation;
   F.params = &params;

   T = gsl_root_fsolver_brent;
   s = gsl_root_fsolver_alloc( T );
   gsl_root_fsolver_set( s, &F, x_lower_bound, x_upper_bound );

   do
   {
      iter++;
      status = gsl_root_fsolver_iterate( s );
      x_Lagrangian_Trial = gsl_root_fsolver_root( s );
      x_lower_bound = gsl_root_fsolver_x_lower( s );
      x_upper_bound = gsl_root_fsolver_x_upper( s );
      status = gsl_root_test_interval( x_lower_bound, x_upper_bound, Pert_Wave_Len/1.0e8, 1.0e-8 );
   }
   while ( status == GSL_CONTINUE  &&  iter < max_iter );

   gsl_root_fsolver_free( s );

   return x_Lagrangian_Trial;
}



double Zeldovich_density_profile( const double x_Lagrangian, const double z )
{
   return rho_crit_i/(1.0-((1.0+z_Collapse)/(1.0+z))*cos(k_Pert*x_Lagrangian));
}



double Zeldovich_x_velocity_profile( const double x_Lagrangian, const double z )
{
   return (-H0*(1.0+z_Collapse)/sqrt(1.0+z)/(1.0+z))*(sin(k_Pert*x_Lagrangian)/k_Pert);
}
#endif // #ifdef SUPPORT_GSL



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
#  if ( defined COMOVING  &&  defined SUPPORT_GSL )
// compute relevant physical quantities
   const double x_shift      = x - 0.5*amr->BoxSize[0]; // unperturbed simulation grids
   const double z_current    = 1.0/Time - 1.0;          // redshift at the current simulation time slice
   const double x_Lagrangian = Zeldovich_coord_transformation_solver( x_shift, z_current );

   if ( Gas_Par_Setup == 1 ) // gas-only setup
   {
      fluid[DENS] = Zeldovich_density_profile( x_Lagrangian, z_current );                // (comoving density)
      fluid[MOMX] = fluid[DENS]*Zeldovich_x_velocity_profile( x_Lagrangian, z_current ); // (comoving density)*(comoving velocity)
      fluid[MOMY] = 0.0;
      fluid[MOMZ] = 0.0;

      const real Pres = EoS_DensTemp2Pres_CPUPtr( fluid[DENS], SQR(A_INIT)*GasTemp_Init*SQR(cbrt(fluid[DENS]/rho_crit_i)),
                                                  NULL, EoS_AuxArray_Flt, EoS_AuxArray_Int, h_EoS_Table );
      const real Eint = EoS_DensPres2Eint_CPUPtr( fluid[DENS], Pres, NULL, EoS_AuxArray_Flt, EoS_AuxArray_Int, h_EoS_Table );
      fluid[ENGY] = Hydro_ConEint2Etot( fluid[DENS], fluid[MOMX], fluid[MOMY], fluid[MOMZ], Eint, 0.0 );
   }

   else if ( Gas_Par_Setup == 2 ) // particle-only setup
   {
      double Dens, MomX, MomY, MomZ, Eint, Etot;
//    gas density and energy cannot be zero so set them to an extremely small value
      Dens = 10.0*TINY_NUMBER;
      MomX = 0.0;
      MomY = 0.0;
      MomZ = 0.0;
      Eint = 10.0*TINY_NUMBER;
      Etot = Hydro_ConEint2Etot( Dens, MomX, MomY, MomZ, Eint, 0.0 );

      fluid[DENS] = Dens;
      fluid[MOMX] = MomX;
      fluid[MOMY] = MomY;
      fluid[MOMZ] = MomZ;
      fluid[ENGY] = Etot;
   }
#  endif // #if ( defined COMOVING  &&  defined SUPPORT_GSL )
} // FUNCTION : SetGridIC



#ifdef PARTICLE
//-------------------------------------------------------------------------------------------------------
// Function    :  Par_Init_ByFunction_Zeldovich
// Description :  User-specified function to initialize particle attributes
//
// Note        :  1. Invoked by Init_GAMER() using the function pointer "Par_Init_ByFunction_Ptr"
//                   --> This function pointer may be reset by various test problem initializers, in which case
//                       this funtion will become useless
//                2. Periodicity should be taken care of in this function
//                   --> No particles should lie outside the simulation box when the periodic BC is adopted
//                   --> However, if the non-periodic BC is adopted, particles are allowed to lie outside the box
//                       (more specifically, outside the "active" region defined by amr->Par->RemoveCell)
//                       in this function. They will later be removed automatically when calling Par_Aux_InitCheck()
//                       in Init_GAMER().
//                3. Particles set by this function are only temporarily stored in this MPI rank
//                   --> They will later be redistributed when calling Par_FindHomePatch_UniformGrid()
//                       and LB_Init_LoadBalance()
//                   --> Therefore, there is no constraint on which particles should be set by this function
//
// Parameter   :  NPar_ThisRank   : Number of particles to be set by this MPI rank
//                NPar_AllRank    : Total Number of particles in all MPI ranks
//                ParMass         : Particle mass     array with the size of NPar_ThisRank
//                ParPosX/Y/Z     : Particle position array with the size of NPar_ThisRank
//                ParVelX/Y/Z     : Particle velocity array with the size of NPar_ThisRank
//                ParTime         : Particle time     array with the size of NPar_ThisRank
//                ParType         : Particle type     array with the size of NPar_ThisRank
//                AllAttributeFlt : Pointer array for all particle floating-point attributes
//                                  --> Dimension = [PAR_NATT_FLT_TOTAL][NPar_ThisRank]
//                                  --> Use the attribute indices defined in Field.h (e.g., Idx_ParCreTime)
//                                      to access the data
//                AllAttributeInt : Pointer array for all particle integer attributes
//                                  --> Dimension = [PAR_NATT_INT_TOTAL][NPar_ThisRank]
//                                  --> Use the attribute indices defined in Field.h to access the data
//
// Return      :  ParMass, ParPosX/Y/Z, ParVelX/Y/Z, ParTime, ParType, AllAttributeFlt, AllAttributeInt
//-------------------------------------------------------------------------------------------------------
void Par_Init_ByFunction_Zeldovich( const long NPar_ThisRank, const long NPar_AllRank,
                                    real_par *ParMass, real_par *ParPosX, real_par *ParPosY, real_par *ParPosZ,
                                    real_par *ParVelX, real_par *ParVelY, real_par *ParVelZ, real_par *ParTime,
                                    long_par *ParType, real_par *AllAttributeFlt[PAR_NATT_FLT_TOTAL],
                                    long_par *AllAttributeInt[PAR_NATT_INT_TOTAL] )
{

#  ifdef SUPPORT_GSL
   if ( Gas_Par_Setup != 2 )   return;


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );


   real_par *ParFltData_AllRank[PAR_NATT_FLT_TOTAL];
   for (int v=0; v<PAR_NATT_FLT_TOTAL; v++)   ParFltData_AllRank[v] = NULL;
   long_par *ParIntData_AllRank[PAR_NATT_INT_TOTAL];
   for (int v=0; v<PAR_NATT_INT_TOTAL; v++)   ParIntData_AllRank[v] = NULL;

// only the master rank will construct the initial condition
   if ( MPI_Rank == 0 )
   {
      double PosVec[3], VelVec[3];
      int NPar_AllRank_Counter = 0;
//    determine the total mass enclosed within the simulation box boundaries
      const double TotM_Boundary = amr->BoxSize[0]*amr->BoxSize[1]*amr->BoxSize[2]*rho_crit_i;
//    determine the individual mass of identical particles
      const double ParM = TotM_Boundary / NPar_AllRank;

      ParFltData_AllRank[PAR_MASS] = new real_par [NPar_AllRank];
      ParFltData_AllRank[PAR_POSX] = new real_par [NPar_AllRank];
      ParFltData_AllRank[PAR_POSY] = new real_par [NPar_AllRank];
      ParFltData_AllRank[PAR_POSZ] = new real_par [NPar_AllRank];
      ParFltData_AllRank[PAR_VELX] = new real_par [NPar_AllRank];
      ParFltData_AllRank[PAR_VELY] = new real_par [NPar_AllRank];
      ParFltData_AllRank[PAR_VELZ] = new real_par [NPar_AllRank];

//    set particle attributes
      for (long px=0; px<NPar_X; px++)
      {
//       x-component position
         PosVec[0] = Zeldovich_coord_transformation_solver( px*dhx, z_Initial );
//       velocity
         VelVec[0] = -Zeldovich_x_velocity_profile( PosVec[0], z_Initial );
         VelVec[1] = 0.0;
         VelVec[2] = 0.0;

         for (long py=0; py<NPar_YZ; py++)
         {
//          y-component position
            PosVec[1] = py*dhx;

            for (long pz=0; pz<NPar_YZ; pz++)
            {
//             mass
               ParFltData_AllRank[PAR_MASS][NPar_AllRank_Counter] = (real_par)ParM;

//             z-component position
               PosVec[2] = pz*dhx;

               for (int d=0; d<3; d++)
               {
                  ParFltData_AllRank[PAR_POSX+d][NPar_AllRank_Counter] = (real_par)PosVec[d];
                  ParFltData_AllRank[PAR_VELX+d][NPar_AllRank_Counter] = (real_par)VelVec[d];
//                check periodicity
                  if ( OPT__BC_FLU[d*2] == BC_FLU_PERIODIC )
                     ParFltData_AllRank[PAR_POSX+d][NPar_AllRank_Counter]
                        = FMOD( ParFltData_AllRank[PAR_POSX+d][NPar_AllRank_Counter]+(real_par)amr->BoxSize[d], (real_par)amr->BoxSize[d] );
               }

               NPar_AllRank_Counter ++;
            } // for (long pz=0; pz<NPar_YZ; pz++)
         } // for (long py=0; py<NPar_YZ; py++)
      } // for (long px=0; px<NPar_X; px++)

      Aux_Message( stdout, "   Total mass within box boundaries = %13.7e\n",  TotM_Boundary   );
      Aux_Message( stdout, "   Total particle number            = %ld\n",     NPar_AllRank    );
      Aux_Message( stdout, "   Particle mass                    = %13.7e\n",  ParM            );
      Aux_Message( stdout, "   x_boundary_size                  = %13.7e\n",  amr->BoxSize[0] );
      Aux_Message( stdout, "   y/z_boundary_size                = %13.7e\n",  amr->BoxSize[1] );
      Aux_Message( stdout, "   z_Initial                        = %13.7e\n",  z_Initial       );
      Aux_Message( stdout, "   z_Collapse                       = %13.7e\n",  z_Collapse      );
      Aux_Message( stdout, "   Pert_Wave_Len                    = %13.7e\n",  Pert_Wave_Len   );
   } // if ( MPI_Rank == 0 )

// send particle attributes from the master rank to all ranks
   Par_ScatterParticleData( NPar_ThisRank, NPar_AllRank, _PAR_MASS|_PAR_POS|_PAR_VEL, _NONE,
                            ParFltData_AllRank, ParIntData_AllRank, AllAttributeFlt, AllAttributeInt );

// synchronize all particles to the physical time on the base level, and set generic particle type
   for (long p=0; p<NPar_ThisRank; p++)
   {
      ParTime[p] = (real_par)Time[0];
      ParType[p] = PTYPE_GENERIC_MASSIVE;
   }

// free memory
   if ( MPI_Rank == 0 )
   {
      for (int v=0; v<PAR_NATT_FLT_TOTAL; v++)   delete [] ParFltData_AllRank[v];
      for (int v=0; v<PAR_NATT_INT_TOTAL; v++)   delete [] ParIntData_AllRank[v];
   }


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );
#  endif // #ifdef SUPPORT_GSL

} // FUNCTION : Par_Init_ByFunction_Zeldovich
#endif // #ifdef PARTICLE



//-------------------------------------------------------------------------------------------------------
// Function    :  OutputError
// Description :  Output numerical error in density, x-velocity, and temperature profiles while "z > z_Collapse"
//
// Note        :  1. Invoke Output_L1Error()
//                2. Use SetGridIC() to provide the analytical solution at any given time
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void OutputError()
{

   const char Prefix[100]     = "Zeldovich";
   const OptOutputPart_t Part = OUTPUT_X + 0;  // along the x-direction

   Output_L1Error( SetGridIC, NULL, Prefix, Part, OUTPUT_PART_X, OUTPUT_PART_Y, OUTPUT_PART_Z );

} // FUNCTION : OutputError
#endif // #if ( MODEL == HYDRO )



//-------------------------------------------------------------------------------------------------------
// Function    :  Init_TestProb_Hydro_Zeldovich
// Description :  Test problem initializer
//
// Note        :  None
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Init_TestProb_Hydro_Zeldovich()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );


// validate the compilation flags and runtime parameters
   Validate();


#  if ( MODEL == HYDRO )
// set the problem-specific runtime parameters
   SetParameter();


// set the function pointers of various problem-specific routines
   Init_Function_User_Ptr  = SetGridIC;
#  ifdef PARTICLE
   Par_Init_ByFunction_Ptr = Par_Init_ByFunction_Zeldovich;
#  endif // #ifdef PARTICLE
   Output_User_Ptr         = ( Gas_Par_Setup == 1 ) ? OutputError : NULL; // output gas numerical and analytical profiles
#  endif // #if ( MODEL == HYDRO )


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : Init_TestProb_Hydro_Zeldovich
