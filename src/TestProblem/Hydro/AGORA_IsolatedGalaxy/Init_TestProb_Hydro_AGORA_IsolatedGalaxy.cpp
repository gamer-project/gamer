#include "GAMER.h"
#include "TestProb.h"



// problem-specific global variables
// =======================================================================================
static char    AGORA_VcProf_Filename[1000];     // filename of the circular velocity radial profile
       char    AGORA_HaloPar_Filename[1000];    // filename of the halo  particles
       char    AGORA_DiskPar_Filename[1000];    // filename of the disk  particles
       char    AGORA_BulgePar_Filename[1000];   // filename of the bulge particles
static double  AGORA_DiskScaleLength;           // disk scale length
static double  AGORA_DiskScaleHeight;           // disk scale height
static double  AGORA_DiskTotalMass;             // disk total mass (gas + stars)
static double  AGORA_DiskGasMassFrac;           // disk gas mass fraction (disk_gas_mass / disk_total_mass)
static double  AGORA_DiskGasTemp;               // disk gas temperature
static double  AGORA_HaloGasNumDensH;           // halo atomic hydrogen number density (halo_gas_mass_density / atomic_hydrogen_mass)
static double  AGORA_HaloGasTemp;               // halo gas temperature

       bool    AGORA_UseMetal = false;          // add and advect a metal density field
                                                // --> to enable this option, one must
                                                //     (1) set AGORA_(Disk/Halo)MetalMassFrac properly
                                                //     (2) set NCOMP_PASSIVE_MAKEFILE>=1 and PAR_NPASSIVE_MAKEFILE>=1 in the Makefile
                                                //     (3) define METAL and PAR_METAL_FRAC in Macro.h (hard coding, ugh!)
                                                // --> necessary if one wants to enable metal_cooling in Grackle
static double  AGORA_DiskMetalMassFrac;         // disk metal mass fraction (disk_metal_mass / disk_gas_mass)
static double  AGORA_HaloMetalMassFrac;         // halo metal mass fraction (halo_metal_mass / halo_gas_mass)


static double  AGORA_DiskGasDens0;              // disk gas mass density coefficient
static double  AGORA_HaloGasDens;               // halo gas mass density
static double  AGORA_HaloGasPres;               // halo gas pressure
static double *AGORA_VcProf = NULL;             // circular velocity radial profile [radius, velocity]
static int     AGORA_VcProf_NBin;               // number of radial bin in AGORA_VcProf
// =======================================================================================


// problem-specific function prototypes
double GaussianQuadratureIntegrate( const double dx, const double dy, const double dz, const double ds );
bool Flag_AGORA( const int i, const int j, const int k, const int lv, const int PID, const double Threshold );
void Par_Init_ByFunction_AGORA( const long NPar_ThisRank, const long NPar_AllRank,
                                real *ParMass, real *ParPosX, real *ParPosY, real *ParPosZ,
                                real *ParVelX, real *ParVelY, real *ParVelZ, real *ParTime,
                                real *ParPassive[PAR_NPASSIVE] );




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

#  ifndef PARTICLE
   Aux_Error( ERROR_INFO, "PARTICLE must be enabled !!\n" );
#  endif

#  ifndef GRAVITY
   Aux_Error( ERROR_INFO, "GRAVITY must be enabled !!\n" );
#  endif

#  ifdef COMOVING
   Aux_Error( ERROR_INFO, "COMOVING must be disabled !!\n" );
#  endif

   if ( !OPT__UNIT )
      Aux_Error( ERROR_INFO, "OPT__UNIT must be enabled !!\n" );

#  ifdef GRAVITY
   if ( OPT__BC_FLU[0] == BC_FLU_PERIODIC  ||  OPT__BC_POT == BC_POT_PERIODIC )
      Aux_Error( ERROR_INFO, "do not use periodic BC for this test !!\n" );
#  endif

#  ifdef PARTICLE
   if ( OPT__INIT == INIT_BY_FUNCTION  &&  amr->Par->Init != PAR_INIT_BY_FUNCTION )
      Aux_Error( ERROR_INFO, "please set PAR_INIT = 1 (by FUNCTION) !!\n" );
#  endif


   if ( MPI_Rank == 0 )
   {
#     if ( FLU_SCHEME != MHM )
         Aux_Message( stderr, "WARNING : it's recommended to adopt the MHM scheme for this test !!\n" );
#     endif

#     ifndef DUAL_ENERGY
         Aux_Message( stderr, "WARNING : it's recommended to enable DUAL_ENERGY for this test !!\n" );
#     endif

#     if ( MODEL == HYDRO  ||  MODEL == MHD )
      if ( MINMOD_COEFF > 1.5 )
         Aux_Message( stderr, "WARNING : it's recommended to set MINMOD_COEFF <= 1.5 for this test !!\n" );

      if ( !JEANS_MIN_PRES )
         Aux_Message( stderr, "WARNING : it's recommended to enable JEANS_MIN_PRES for this test !!\n" );
#     endif

      if ( FLAG_BUFFER_SIZE < 5 )
         Aux_Message( stderr, "WARNING : it's recommended to set FLAG_BUFFER_SIZE >= 5 for this test !!\n" );

      if ( FLAG_BUFFER_SIZE_MAXM1_LV < 2 )
         Aux_Message( stderr, "WARNING : it's recommended to set FLAG_BUFFER_SIZE_MAXM1_LV >= 2 for this test !!\n" );

      if ( amr->BoxSize[0] != amr->BoxSize[1]  ||  amr->BoxSize[0] != amr->BoxSize[2] )
         Aux_Message( stderr, "WARNING : non-cubic box (currently the flag routine \"Flag_AGORA()\" assumes a cubic box) !!\n" );

      if ( INIT_SUBSAMPLING_NCELL != 0 )
         Aux_Message( stderr, "WARNING : INIT_SUBSAMPLING_NCELL (%d) != 0 will lead to non-uniform initial disk temperature !!\n",
                      INIT_SUBSAMPLING_NCELL );
   } // if ( MPI_Rank == 0 )


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Validating test problem %d ... done\n", TESTPROB_ID );

} // FUNCTION : Validate



#if ( MODEL == HYDRO  &&  defined PARTICLE )
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

// add parameters in the following format (some handy constants are defined in TestProb.h):
// --> note that VARIABLE, DEFAULT, MIN, and MAX must have the same data type
// --> some handy constants (e.g., NoMin_int, Eps_float, ...) are defined in "ReadPara.h"
// ********************************************************************************************************************************
// ReadPara->Add( "KEY_IN_THE_FILE",         &VARIABLE,                 DEFAULT,       MIN,              MAX               );
// ********************************************************************************************************************************
   ReadPara->Add( "AGORA_VcProf_Filename",    AGORA_VcProf_Filename,    Useless_str,   Useless_str,      Useless_str       );
   ReadPara->Add( "AGORA_HaloPar_Filename",   AGORA_HaloPar_Filename,   Useless_str,   Useless_str,      Useless_str       );
   ReadPara->Add( "AGORA_DiskPar_Filename",   AGORA_DiskPar_Filename,   Useless_str,   Useless_str,      Useless_str       );
   ReadPara->Add( "AGORA_BulgePar_Filename",  AGORA_BulgePar_Filename,  Useless_str,   Useless_str,      Useless_str       );
   ReadPara->Add( "AGORA_DiskScaleLength",   &AGORA_DiskScaleLength,    -1.0,          Eps_double,       NoMax_double      );
   ReadPara->Add( "AGORA_DiskScaleHeight",   &AGORA_DiskScaleHeight,    -1.0,          Eps_double,       NoMax_double      );
   ReadPara->Add( "AGORA_DiskTotalMass",     &AGORA_DiskTotalMass,      -1.0,          Eps_double,       NoMax_double      );
   ReadPara->Add( "AGORA_DiskGasMassFrac",   &AGORA_DiskGasMassFrac,    -1.0,          Eps_double,       NoMax_double      );
   ReadPara->Add( "AGORA_DiskGasTemp",       &AGORA_DiskGasTemp,        -1.0,          Eps_double,       NoMax_double      );
   ReadPara->Add( "AGORA_HaloGasNumDensH",   &AGORA_HaloGasNumDensH,    -1.0,          Eps_double,       NoMax_double      );
   ReadPara->Add( "AGORA_HaloGasTemp",       &AGORA_HaloGasTemp,        -1.0,          Eps_double,       NoMax_double      );
   ReadPara->Add( "AGORA_UseMetal",          &AGORA_UseMetal,            false,        Useless_bool,     Useless_bool      );
   ReadPara->Add( "AGORA_DiskMetalMassFrac", &AGORA_DiskMetalMassFrac,   0.0,          0.0,              1.0               );
   ReadPara->Add( "AGORA_HaloMetalMassFrac", &AGORA_HaloMetalMassFrac,   0.0,          0.0,              1.0               );

   ReadPara->Read( FileName );

   delete ReadPara;

// check the runtime parameters
   if ( AGORA_UseMetal )
   {
#     if (  ( defined DUAL_ENERGY && NCOMP_PASSIVE < 2 )  ||  ( !defined DUAL_ENERGY && NCOMP_PASSIVE < 1 )  )
         Aux_Error( ERROR_INFO, "please set NCOMP_PASSIVE_MAKEFILE >= 1 in the Makefile for \"AGORA_UseMetal\" !!\n" );
#     endif

#     ifndef METAL
         Aux_Error( ERROR_INFO, "please define the symbolic constant \"METAL\" properly in Macro.h for \"AGORA_UseMetal\" !!\n" );
#     endif

#     if (  ( defined STAR_FORMATION && PAR_NPASSIVE < 2 )  ||  ( !defined STAR_FORMATION && PAR_NPASSIVE < 1 )  )
         Aux_Error( ERROR_INFO, "please set PAR_NPASSIVE_MAKEFILE >= 1 in the Makefile for \"AGORA_UseMetal\" !!\n" );
#     endif

#     ifndef PAR_METAL_FRAC
         Aux_Error( ERROR_INFO, "please define the symbolic constant \"PAR_METAL_FRAC\" properly in Macro.h for \"AGORA_UseMetal\" !!\n" );
#     endif
   }

   else
   {
#     ifdef SUPPORT_GRACKLE
      if ( GRACKLE_METAL )
         Aux_Error( ERROR_INFO, "please enable \"AGORA_UseMetal\" for \"GRACKLE_METAL\" !!\n" );
#     endif
   }

// convert to code units
   AGORA_DiskScaleLength *= Const_kpc  / UNIT_L;
   AGORA_DiskScaleHeight *= Const_kpc  / UNIT_L;
   AGORA_DiskTotalMass   *= Const_Msun / UNIT_M;
   AGORA_DiskGasTemp     *= Const_kB   / UNIT_E;
   AGORA_HaloGasNumDensH *= 1.0        * CUBE(UNIT_L);
   AGORA_HaloGasTemp     *= Const_kB   / UNIT_E;

// load the circular velocity radial profile
   if ( OPT__INIT != INIT_BY_RESTART )
   {
      const bool RowMajor_No = false;     // load data into the column-major order
      const bool AllocMem_Yes = true;     // allocate memory for AGORA_VcProf
      const int  NCol        = 2;         // total number of columns to load
      const int  Col[NCol]   = {0, 1};    // target columns: (radius, density)

      AGORA_VcProf_NBin = Aux_LoadTable( AGORA_VcProf, AGORA_VcProf_Filename, NCol, Col, RowMajor_No, AllocMem_Yes );

//    convert to code units (assuming the loaded radius and velocity are in kpc and km/s, respectively)
      double *Prof_R = AGORA_VcProf + 0*AGORA_VcProf_NBin;
      double *Prof_V = AGORA_VcProf + 1*AGORA_VcProf_NBin;

      for (int b=0; b<AGORA_VcProf_NBin; b++)
      {
         Prof_R[b] *= Const_kpc / UNIT_L;
         Prof_V[b] *= Const_km  / UNIT_V;
      }
   }


// (2) set the problem-specific derived parameters
   const bool CheckMinPres_Yes = true;
   AGORA_DiskGasDens0 = AGORA_DiskTotalMass*AGORA_DiskGasMassFrac / ( 4.0*M_PI*SQR(AGORA_DiskScaleLength)*AGORA_DiskScaleHeight );
   AGORA_HaloGasDens  = AGORA_HaloGasNumDensH*Const_mH/UNIT_M;
   AGORA_HaloGasPres  = CPU_Temperature2Pressure( AGORA_HaloGasDens, AGORA_HaloGasTemp, MOLECULAR_WEIGHT, Const_mH/UNIT_M,
                                                  CheckMinPres_Yes, MIN_PRES );


// (3) reset other general-purpose parameters
//     --> a helper macro PRINT_WARNING is defined in TestProb.h
   const long   End_Step_Default = __INT_MAX__;
   const double End_T_Default    =  500.0*Const_Myr/UNIT_T;

   if ( END_STEP < 0 ) {
      END_STEP = End_Step_Default;
      PRINT_WARNING( "END_STEP", END_STEP, FORMAT_LONG );
   }

   if ( END_T < 0.0 ) {
      END_T = End_T_Default;
      PRINT_WARNING( "END_T", END_T, FORMAT_REAL );
   }


// (4) make a note
   if ( MPI_Rank == 0 )
   {
      Aux_Message( stdout, "=============================================================================\n" );
      Aux_Message( stdout, "  test problem ID   = %d\n",             TESTPROB_ID                               );
      Aux_Message( stdout, "  VcProf_Filename   = %s\n",             AGORA_VcProf_Filename                     );
      Aux_Message( stdout, "  HaloPar_Filename  = %s\n",             AGORA_HaloPar_Filename                    );
      Aux_Message( stdout, "  DiskPar_Filename  = %s\n",             AGORA_DiskPar_Filename                    );
      Aux_Message( stdout, "  BulgePar_Filename = %s\n",             AGORA_BulgePar_Filename                   );
      Aux_Message( stdout, "  DiskScaleLength   = %13.7e kpc\n",     AGORA_DiskScaleLength * UNIT_L/Const_kpc  );
      Aux_Message( stdout, "  DiskScaleHeight   = %13.7e kpc\n",     AGORA_DiskScaleHeight * UNIT_L/Const_kpc  );
      Aux_Message( stdout, "  DiskTotalMass     = %13.7e Msun\n",    AGORA_DiskTotalMass   * UNIT_M/Const_Msun );
      Aux_Message( stdout, "  DiskGasMassFrac   = %13.7e\n",         AGORA_DiskGasMassFrac                     );
      Aux_Message( stdout, "  DiskGasTemp       = %13.7e K\n",       AGORA_DiskGasTemp     * UNIT_E/Const_kB   );
      Aux_Message( stdout, "  HaloGasNumDensH   = %13.7e cm^{-3}\n", AGORA_HaloGasNumDensH / CUBE(UNIT_L)      );
      Aux_Message( stdout, "  HaloGasTemp       = %13.7e K\n",       AGORA_HaloGasTemp     * UNIT_E/Const_kB   );
      Aux_Message( stdout, "  UseMetal          = %d\n",             AGORA_UseMetal                            );
      if ( AGORA_UseMetal ) {
      Aux_Message( stdout, "  DiskMetalMassFrac = %13.7e\n",         AGORA_DiskMetalMassFrac                   );
      Aux_Message( stdout, "  HaloMetalMassFrac = %13.7e\n",         AGORA_HaloMetalMassFrac                   ); }
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
void SetGridIC( real fluid[], const double x, const double y, const double z, const double Time,
                const int lv, double AuxArray[] )
{

   const bool   CheckMinPres_Yes = true;
   const double dx               = x - 0.5*amr->BoxSize[0];
   const double dy               = y - 0.5*amr->BoxSize[1];
   const double dz               = z - 0.5*amr->BoxSize[2];
   const double r                = sqrt( SQR(dx) + SQR(dy) );
   const double h                = fabs( dz );

   double DiskGasDens, DiskGasPres, DiskGasVel;

// use the 5-point Gaussian quadrature integration to improve the accuracy of the average cell density
// --> we don't want to use the sub-sampling for that purpose since it will break the initial uniform temperature
// DiskGasDens = AGORA_DiskGasDens0 * exp( -r/AGORA_DiskScaleLength ) * exp( -h/AGORA_DiskScaleHeight );
   DiskGasDens = GaussianQuadratureIntegrate( dx, dy, dz, amr->dh[lv] );
   DiskGasPres = CPU_Temperature2Pressure( DiskGasDens, AGORA_DiskGasTemp, MOLECULAR_WEIGHT, Const_mH/UNIT_M,
                                           CheckMinPres_Yes, MIN_PRES );

// disk component
   if ( DiskGasPres > AGORA_HaloGasPres )
   {
//    perform linear interpolation to get the gas circular velocity at r
      const double *Prof_R = AGORA_VcProf + 0*AGORA_VcProf_NBin;
      const double *Prof_V = AGORA_VcProf + 1*AGORA_VcProf_NBin;

      if (  ( DiskGasVel = Mis_InterpolateFromTable(AGORA_VcProf_NBin, Prof_R, Prof_V, r) ) == NULL_REAL  )
         Aux_Error( ERROR_INFO, "interpolation failed at radius %13.7e !!\n", r );

      fluid[DENS] = DiskGasDens;
      fluid[MOMX] = -dy/r*DiskGasVel*fluid[DENS];
      fluid[MOMY] = +dx/r*DiskGasVel*fluid[DENS];
      fluid[MOMZ] = 0.0;
      fluid[ENGY] = DiskGasPres / ( GAMMA - 1.0 )
                    + 0.5*( SQR(fluid[MOMX]) + SQR(fluid[MOMY]) + SQR(fluid[MOMZ]) ) / fluid[DENS];

//###: HARD-CODED FIELDS
#     if ( NCOMP_PASSIVE > 0 )
      if ( AGORA_UseMetal )
      fluid[METAL] = fluid[DENS]*AGORA_DiskMetalMassFrac;
#     endif
   }

// halo component
   else
   {
      fluid[DENS] = AGORA_HaloGasDens;
      fluid[MOMX] = 0.0;
      fluid[MOMY] = 0.0;
      fluid[MOMZ] = 0.0;
      fluid[ENGY] = AGORA_HaloGasPres / ( GAMMA - 1.0 )
                    + 0.5*( SQR(fluid[MOMX]) + SQR(fluid[MOMY]) + SQR(fluid[MOMZ]) ) / fluid[DENS];

//###: HARD-CODED FIELDS
#     if ( NCOMP_PASSIVE > 0 )
      if ( AGORA_UseMetal )
      fluid[METAL] = fluid[DENS]*AGORA_HaloMetalMassFrac;
#     endif
   } // if ( DiskPres > AGORA_HaloGasPres ) ... else ...

} // FUNCTION : SetGridIC


//-------------------------------------------------------------------------------------------------------
// Function    :  GaussianQuadratureIntegrate
// Description :  Use the 5-point Gaussian quadrature integration to calculate the average density of a given cell
//
// Note        :  1. This function mimics the AGORA isolated galaxy setup implemented in Enzo
//                2. Ref: http://mathworld.wolfram.com/Legendre-GaussQuadrature.html
//                        https://en.wikipedia.org/wiki/Gaussian_quadrature
//                3. 1D: Int[ f(x), {a,b} ] ~ 0.5*(b-a)*Sum_i[ weight_i*f( 0.5*(b-a)*abscissas_i + 0.5*(a+b) ) ]
//                   3D: Int[ f(x,y,z), { {xc-0.5*dx,xc+0.5*dx}, {yc-0.5*dy,yc+0.5*dy}, {zc-0.5*dz,zc+0.5*dz} } ]
//                       ~ dx*dy*dz/8 * Sum_i,j,k[ weight_i*weight_j*weight_k
//                                                 *f( xc+0.5*dx*abscissas_i, yc+0.5*dy*abscissas_j, zc+0.5*dz*abscissas_k ) ]
//
// Parameter   :  dx/dy/dz : Central coordinates of the target cell relative to the box center
//                ds       : Cell size
//
// Return      :  Average density
//-------------------------------------------------------------------------------------------------------
double GaussianQuadratureIntegrate( const double dx, const double dy, const double dz, const double ds )
{

   const int    NPoint            = 5;
   const double abscissas[NPoint] = { -0.90617984593866, -0.53846931010568, +0.0,              +0.53846931010568, +0.90617984593866 };
   const double weight   [NPoint] = { +0.23692688505619, +0.47862867049937, +0.56888888888889, +0.47862867049937, +0.23692688505619 };
   const double ds_2              = 0.5*ds;

   double dxx, dyy, dzz, r, h;
   double Dens = 0;

   for (int k=0; k<NPoint; k++)  {  dzz = dz + ds_2*abscissas[k];
                                    h   = fabs( dzz );
   for (int j=0; j<NPoint; j++)  {  dyy = dy + ds_2*abscissas[j];
   for (int i=0; i<NPoint; i++)  {  dxx = dx + ds_2*abscissas[i];
                                    r   = sqrt( SQR(dxx) + SQR(dyy) );

     Dens += weight[i]*weight[j]*weight[k] * exp( -r/AGORA_DiskScaleLength ) * exp( -h/AGORA_DiskScaleHeight );

   }}}

   Dens *= AGORA_DiskGasDens0/8.0;

   return Dens;

} // FUNCTION : GaussianQuadratureIntegrate



//-------------------------------------------------------------------------------------------------------
// Function    :  End_AGORA
// Description :  Free memory before terminating the program
//
// Note        :  1. Linked to the function pointer "End_User_Ptr" to replace "End_User()"
//
// Parameter   :  None
//-------------------------------------------------------------------------------------------------------
void End_AGORA()
{

   delete [] AGORA_VcProf;
   AGORA_VcProf = NULL;

} // FUNCTION : End_AGORA
#endif // #if ( MODEL == HYDRO  &&  defined PARTICLE )



//-------------------------------------------------------------------------------------------------------
// Function    :  Init_TestProb_Hydro_AGORA_IsolatedGalaxy
// Description :  Test problem initializer
//
// Note        :  None
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Init_TestProb_Hydro_AGORA_IsolatedGalaxy()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );


// validate the compilation flags and runtime parameters
   Validate();


#  if ( MODEL == HYDRO  &&  defined PARTICLE )
// set the problem-specific runtime parameters
   SetParameter();


// set the function pointers of various problem-specific routines
   Init_Function_User_Ptr   = SetGridIC;
   Output_User_Ptr          = NULL;
   Flag_User_Ptr            = Flag_AGORA;
   Mis_GetTimeStep_User_Ptr = NULL;
   Aux_Record_User_Ptr      = NULL;
   BC_User_Ptr              = NULL;
   Flu_ResetByUser_Func_Ptr = NULL;
   End_User_Ptr             = End_AGORA;
   Init_ExternalAcc_Ptr     = NULL;
   Init_ExternalPot_Ptr     = NULL;
   Par_Init_ByFunction_Ptr  = Par_Init_ByFunction_AGORA;
#  endif // if ( MODEL == HYDRO  &&  defined PARTICLE )


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : Init_TestProb_Hydro_AGORA_IsolatedGalaxy
