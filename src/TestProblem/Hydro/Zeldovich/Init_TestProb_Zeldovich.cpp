#include "GAMER.h"
#include "TestProb.h"
#include <math.h>
#ifdef SUPPORT_GSL
   #include <gsl/gsl_errno.h>
   #include <gsl/gsl_roots.h>
#endif



// problem-specific global variables
// =======================================================================================
// =======================================================================================
static int Gas_Par_Setup;        // gas-only or particle-only setup
static int n_Pert_Wave_Len;      // "Pert_Wave_Len = BOXSIZE_x/n_Pert_Wave_Len"
static int NPar_X;               // 1D number of particles along the x-axis (axis of perturbation)
static int NPar_YZ;              // 1D number of particles along the y/z-axis
static double Pert_Wave_Len;     // perturbation wavelength [Mpc/h]
static double k_Pert;            // perturbation wavenumber
static double zi_Initial;        // redshfit of simulation IC
static double zc_Collapse;       // redshfit at the onset of Zeldovich pancake collapse
static double GasTemp_Init;      // initial gas temperature in Kelvin
static double H0;                // Hubble parameter at redshift zero
static double rho_crit_i;        // comoving critical density = 3H0^2/8*pi*G
static double dhx;               // 1D inter-particle separation [Mpc/h]
// =======================================================================================


// problem-specific function prototypes
struct Zeldovich_coord_params
{
   double x_shift_input, zi_Initial;
};
static double Zeldovich_coord_transformation(double x_Lagrangian_Trial, void *params);
static double Zeldovich_coord_transformation_solver(double x_shift_input, double zi_Initial);
static double Zeldovich_density_profile(double x_Lagrangian, double zi_Initial);
static double Zeldovich_x_velocity_profile(double x_Lagrangian, double zi_Initial);
static void OutputError();




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

#  ifdef GRAVITY
   if ( OPT__BC_FLU[0] != BC_FLU_PERIODIC  ||  OPT__BC_POT != BC_POT_PERIODIC )
      Aux_Error( ERROR_INFO, "must adopt periodic BC for this test !!\n" );
#  endif

   if ( OPT__INIT != INIT_BY_FUNCTION  &&  OPT__INIT != INIT_BY_RESTART )
      Aux_Error( ERROR_INFO, "OPT__INIT != FUNCTION (1) or RESTART (2) for this test !!\n" );

   if ( OMEGA_M0 != 1 )
      Aux_Error( ERROR_INFO, "OMEGA_M0 != 1 (matter-only) for this test !!\n");

#  ifndef FLOAT8
   Aux_Error( ERROR_INFO, "FLOAT8 must be enabled !!\n" );
#  endif

#  ifndef SUPPORT_GSL
   Aux_Error( ERROR_INFO, "SUPPORT_GSL must be enabled !!\n" );
#  endif

#  ifdef PARTICLE
   if ( amr->Par->Init != PAR_INIT_BY_FUNCTION )
      Aux_Error( ERROR_INFO, "PAR_INIT != FUNCTION (1) for this test !!\n" );
#  endif


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
   ReadPara->Add( "n_Pert_Wave_Len",   &n_Pert_Wave_Len,        2,            1,                NoMax_int         );
   ReadPara->Add( "zc_Collapse",       &zc_Collapse,            1.0,          0.0,              NoMax_double      );
   ReadPara->Add( "GasTemp_Init",      &GasTemp_Init,           100.0,        0.0,              NoMax_double      );
   ReadPara->Add( "NPar_X",            &NPar_X,                 NX0_TOT[0],    0,                NoMax_int         );

   ReadPara->Read( FileName );

   delete ReadPara;


// (1-2) set the default values

// (1-3) check the runtime parameters


// (2) set the problem-specific derived parameters
   H0            = HUBBLE0*(100*Const_km/(Const_s*Const_Mpc))*UNIT_T;
   Pert_Wave_Len = amr->BoxSize[0]/n_Pert_Wave_Len;
   rho_crit_i    = 1;             // [UNIT_D]; comoving critical density = 3H0^2/8*pi*G
   k_Pert        = 2.0*M_PI/Pert_Wave_Len;
   dhx           = amr->BoxSize[0] / NPar_X;
   NPar_YZ       = amr->BoxSize[1] / dhx;


// (3) reset other general-purpose parameters
//     --> a helper macro PRINT_WARNING is defined in TestProb.h
   const long   End_Step_Default = __INT_MAX__;
   const double End_T_Default    = 1.0;

   if ( END_STEP < 0 ) {
      END_STEP = End_Step_Default;
      PRINT_WARNING( "END_STEP", END_STEP, FORMAT_LONG );
   }

   if ( END_T < 0.0 ) {
      END_T = End_T_Default;
      PRINT_WARNING( "END_T", END_T, FORMAT_REAL );
   }

#  ifdef PARTICLE
   if ( Gas_Par_Setup == 1 ) // overwrite the total number of particles in gas-only setup
   {
      amr->Par->NPar_Active_AllRank = 0;
   }
   if ( Gas_Par_Setup == 2 ) // overwrite the total number of particles in particle-only setup
   {
      amr->Par->NPar_Active_AllRank = NPar_X*SQR(NPar_YZ);
   }
#  endif


// (4) make a note
   if ( MPI_Rank == 0 )
   {
      Aux_Message( stdout, "=============================================================================\n" );
      Aux_Message( stdout, "  test problem ID                   = %d\n",     TESTPROB_ID     );
      Aux_Message( stdout, "  gas[1]/particle[2] setup          = %d\n",     Gas_Par_Setup   );
      Aux_Message( stdout, "  perturbation wavelength multiple  = %14.7e\n", n_Pert_Wave_Len );
      Aux_Message( stdout, "  Zeldovich collapse redshift       = %14.7e\n", zc_Collapse     );
      Aux_Message( stdout, "  initial gas temperature in Kelvin = %14.7e\n", GasTemp_Init    );
      Aux_Message( stdout, "  particle number along the x-axis  = %d\n",     NPar_X          );
      Aux_Message( stdout, "=============================================================================\n" );
   }


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Setting runtime parameters ... done\n" );

} // FUNCTION : SetParameter



// compute conversion between Lagrangian (x_Lagrangian) and Eulerian/code coordinates (x_shift)
double Zeldovich_coord_transformation(double x_Lagrangian_Trial, void *params)
{
   struct Zeldovich_coord_params *p = (struct Zeldovich_coord_params *) params;
   double x_shift_input = p->x_shift_input;
   double zi_Initial = p->zi_Initial;

   return (x_Lagrangian_Trial - ((1+zc_Collapse)/(1+zi_Initial))*sin(k_Pert*x_Lagrangian_Trial)/k_Pert) - x_shift_input;
}
#ifdef SUPPORT_GSL
   double Zeldovich_coord_transformation_solver(double x_shift_input, double zi_Initial)
   {
      int status;
      int iter = 0, max_iter = 100;
      const gsl_root_fsolver_type *T;
      gsl_root_fsolver *s;
      double x_Lagrangian_Trial = x_shift_input;
      double x_lower_bound = x_shift_input-Pert_Wave_Len/4, x_upper_bound = x_shift_input+Pert_Wave_Len/4;
      gsl_function F;
      struct Zeldovich_coord_params params = {x_shift_input, zi_Initial};

      F.function = &Zeldovich_coord_transformation;
      F.params = &params;

      T = gsl_root_fsolver_brent;
      s = gsl_root_fsolver_alloc (T);
      gsl_root_fsolver_set (s, &F, x_lower_bound, x_upper_bound);

      do
      {
         iter++;
         status = gsl_root_fsolver_iterate (s);
         x_Lagrangian_Trial = gsl_root_fsolver_root (s);
         x_lower_bound = gsl_root_fsolver_x_lower (s);
         x_upper_bound = gsl_root_fsolver_x_upper (s);
         status = gsl_root_test_interval (x_lower_bound, x_upper_bound, Pert_Wave_Len/1e5, 1e-5);
      }
      while (status == GSL_CONTINUE && iter < max_iter);

      gsl_root_fsolver_free (s);

      return x_Lagrangian_Trial;
   }
double Zeldovich_density_profile(double x_Lagrangian, double zi_Initial)
{
   return rho_crit_i/(1-((1+zc_Collapse)/(1+zi_Initial))*cos(k_Pert*x_Lagrangian));
}
double Zeldovich_x_velocity_profile(double x_Lagrangian, double zi_Initial)
{
   return (-H0*(1+zc_Collapse)/sqrt(1+zi_Initial)/(1+zi_Initial))*(sin(k_Pert*x_Lagrangian)/k_Pert);
}
#endif



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
   #ifdef SUPPORT_GSL
//    compute relevant physical quantities
      zi_Initial = (1/Time) - 1;
      const double x_shift      = x - 0.5*amr->BoxSize[0]; // unperturbed simulation grids
      const double x_Lagrangian = Zeldovich_coord_transformation_solver(x_shift,zi_Initial);

   #  if   ( MODEL == HYDRO )
      if ( Gas_Par_Setup == 1 ) // gas-only setup
      {
      fluid[DENS] = Zeldovich_density_profile(x_Lagrangian, zi_Initial); // (comoving density)
      fluid[MOMX] = fluid[DENS]*Zeldovich_x_velocity_profile(x_Lagrangian, zi_Initial); // (comoving density)*(comoving velocity)
      fluid[MOMY] = 0.0;
      fluid[MOMZ] = 0.0;

      const double Pres = EoS_DensTemp2Pres_CPUPtr( fluid[DENS], SQR(A_INIT)*GasTemp_Init*SQR(cbrt(fluid[DENS]/rho_crit_i)),
                                                    NULL, EoS_AuxArray_Flt, EoS_AuxArray_Int, h_EoS_Table );
      const double Eint = EoS_DensPres2Eint_CPUPtr( fluid[DENS], Pres, NULL, EoS_AuxArray_Flt, EoS_AuxArray_Int, h_EoS_Table );
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
   #endif // #ifdef SUPPORT_GSL
      }
   #  endif //  #if ( MODEL == HYDRO )

   } // FUNCTION : SetGridIC



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
// Parameter   :  NPar_ThisRank : Number of particles to be set by this MPI rank
//                NPar_AllRank  : Total Number of particles in all MPI ranks
//                ParMass       : Particle mass     array with the size of NPar_ThisRank
//                ParPosX/Y/Z   : Particle position array with the size of NPar_ThisRank
//                ParVelX/Y/Z   : Particle velocity array with the size of NPar_ThisRank
//                ParTime       : Particle time     array with the size of NPar_ThisRank
//                ParType       : Particle type     array with the size of NPar_ThisRank
//                AllAttribute  : Pointer array for all particle attributes
//                                --> Dimension = [PAR_NATT_TOTAL][NPar_ThisRank]
//                                --> Use the attribute indices defined in Field.h (e.g., Idx_ParCreTime)
//                                    to access the data
//
// Return      :  ParMass, ParPosX/Y/Z, ParVelX/Y/Z, ParTime, ParType, AllAttribute
//-------------------------------------------------------------------------------------------------------
void Par_Init_ByFunction_Zeldovich( const long NPar_ThisRank, const long NPar_AllRank,
                                  real *ParMass, real *ParPosX, real *ParPosY, real *ParPosZ,
                                  real *ParVelX, real *ParVelY, real *ParVelZ, real *ParTime,
                                  real *ParType, real *AllAttribute[PAR_NATT_TOTAL] )
{
   if ( Gas_Par_Setup == 2 )
   {
      zi_Initial = (1/Time[0]) - 1;

      if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );


      real *ParData_AllRank[PAR_NATT_TOTAL];
      for (int v=0; v<PAR_NATT_TOTAL; v++)   ParData_AllRank[v] = NULL;

//    only the master rank will construct the initial condition
      if ( MPI_Rank == 0 )
      {
         double PosVec[3], VelVec[3];
         int NPar_AllRank_Counter = 0;
//       determine the total mass enclosed within the simulation box boundaries
         const double TotM_Boundary = amr->BoxSize[0]*amr->BoxSize[1]*amr->BoxSize[2]*rho_crit_i;
//       determine the individual mass of identical particles
         const double ParM = TotM_Boundary / NPar_AllRank;

         ParData_AllRank[PAR_MASS] = new real [NPar_AllRank];
         ParData_AllRank[PAR_POSX] = new real [NPar_AllRank];
         ParData_AllRank[PAR_POSY] = new real [NPar_AllRank];
         ParData_AllRank[PAR_POSZ] = new real [NPar_AllRank];
         ParData_AllRank[PAR_VELX] = new real [NPar_AllRank];
         ParData_AllRank[PAR_VELY] = new real [NPar_AllRank];
         ParData_AllRank[PAR_VELZ] = new real [NPar_AllRank];

//       set particle attributes
         for (long px=0; px<NPar_X; px++)
         {
//       x-component position
         PosVec[0] = Zeldovich_coord_transformation_solver(px*dhx,zi_Initial);
//       velocity
         VelVec[0] = -Zeldovich_x_velocity_profile(PosVec[0], zi_Initial);
         VelVec[1] = 0;
         VelVec[2] = 0;

            for (long py=0; py<NPar_YZ; py++)
            {
//          y-component position
            PosVec[1] = py*dhx;

               for (long pz=0; pz<NPar_YZ; pz++)
               {
//             mass
               ParData_AllRank[PAR_MASS][NPar_AllRank_Counter] = ParM;

//             z-component position
               PosVec[2] = pz*dhx;

                  for (int d=0; d<3; d++)
                  {
                     ParData_AllRank[PAR_POSX+d][NPar_AllRank_Counter] = PosVec[d];
                     ParData_AllRank[PAR_VELX+d][NPar_AllRank_Counter] = VelVec[d];
//                 check periodicity
                     if ( OPT__BC_FLU[d*2] == BC_FLU_PERIODIC )
                        ParData_AllRank[PAR_POSX+d][NPar_AllRank_Counter] = FMOD( ParData_AllRank[PAR_POSX+d][NPar_AllRank_Counter]+(real)amr->BoxSize[d], (real)amr->BoxSize[d] );
                  }

               NPar_AllRank_Counter += 1;
               }
         }
      }

      Aux_Message( stdout, "   Total mass within box boundaries = %13.7e\n",  TotM_Boundary );
      Aux_Message( stdout, "   Total particle number            = %d\n",      NPar_AllRank );
      Aux_Message( stdout, "   Particle mass                    = %13.7e\n",  ParM );
      Aux_Message( stdout, "   x_boundary_size                  = %13.7e\n",  amr->BoxSize[0] );
      Aux_Message( stdout, "   y/z_boundary_size                = %13.7e\n",  amr->BoxSize[1] );
      Aux_Message( stdout, "   zi_Initial                       = %13.7e\n",  zi_Initial );
      Aux_Message( stdout, "   zc_Collapse                      = %13.7e\n",  zc_Collapse );
      Aux_Message( stdout, "   Pert_Wave_Len                    = %13.7e\n",  Pert_Wave_Len );

   } // if ( MPI_Rank == 0 )

// send particle attributes from the master rank to all ranks
   Par_ScatterParticleData( NPar_ThisRank, NPar_AllRank, _PAR_MASS|_PAR_POS|_PAR_VEL, ParData_AllRank, AllAttribute );

// synchronize all particles to the physical time on the base level, and set generic particle type
   for (long p=0; p<NPar_ThisRank; p++) {
      ParTime[p] = Time[0];
      ParType[p] = PTYPE_GENERIC_MASSIVE;
   }

// free resource
   if ( MPI_Rank == 0 )
   {
      for (int v=0; v<PAR_NATT_TOTAL; v++)   delete [] ParData_AllRank[v];
   }


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );
   }
} // FUNCTION : Par_Init_ByFunction_Zeldovich



//-------------------------------------------------------------------------------------------------------
// Function    :  OutputError
// Description :  Output numerical error in density, x-velocity, and temperature profiles while "z > zc_Collapse"
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
   #ifdef SUPPORT_GSL
   if (zi_Initial > zc_Collapse)
   {
   const char Prefix[100]     = "Zeldovich";
   const OptOutputPart_t Part = OUTPUT_X + 0;  // along the x-direction

   Output_L1Error( SetGridIC, NULL, Prefix, Part, OUTPUT_PART_X, OUTPUT_PART_Y, OUTPUT_PART_Z );
   }
   #endif
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
   Init_Function_User_Ptr = SetGridIC;
   Par_Init_ByFunction_Ptr = Par_Init_ByFunction_Zeldovich; // Point mass function pointer -> initialization by function

   if ( Gas_Par_Setup == 1 )  // gas-only setup
   {
   Output_User_Ptr = OutputError; // output gas numerical and analytical profiles
   }


#  endif // #if ( MODEL == HYDRO )


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : Init_TestProb_Hydro_Zeldovich