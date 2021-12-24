#include "GAMER.h"
#include "TestProb.h"



// problem-specific global variables
// =======================================================================================
static int      Soliton_N;                               // total number of solitons
static int      Soliton_RSeed;                           // random seed for setting Soliton_Center[]
                                                         //    (<0 --> hard coding each soliton)
static double   Soliton_CoreRadiusAll;                   // core radius for all solitons
                                                         //    (<=0.0 --> hard coding each soliton)
static double   Soliton_EmptyRegion;                     // soliton-free region from the boundary
                                                         //    (useful only when Soliton_RSeed>=0)
static char     Soliton_DensProf_Filename[MAX_STRING];   // filename of the reference soliton density profile
static char     Soliton_Distortion_Filename[MAX_STRING]; // filename for declaring the soliton distortion ratio along x, y, z

static int      Soliton_DensProf_NBin;                   // number of radial bins of the soliton density profile
static int      Soliton_Distortion;                      // distort the soliton or not (0: no distort; 1: distort)
static double  *Soliton_Distortion_Ratio = NULL;         // distortion ratio along x, y, z axis
static double  *Soliton_DensProf    = NULL;              // soliton density profile [radius/density]
static double  *Soliton_CoreRadius  = NULL;              // core radius of each soliton
static double (*Soliton_Center)[3]  = NULL;              // center coordinates of each soliton
static double  *Soliton_ScaleL      = NULL;              // L/D: length/density scale factors of each soliton
                                                         //      (defined as the ratio between the core radii/peak
                                                         //      density of the target and reference soliton profiles)
static double  *Soliton_ScaleD      = NULL;
static bool     first_run_flag      = true;              // flag suggesting first run (for determining whether write header in log file or not )

static double   Fake_Soliton_CM_MaxR;                          // maximum radius for determining System CM
static double   Fake_Soliton_CM_TolErrR;                       // maximum allowed errors for determining System CM

#ifdef PARTICLE
static long     NPar_AllRank_Check;
static char     Particle_Data_Filename[MAX_STRING];      // filename of the particles mass, initial position, initial velocity data
static char     Particle_Log_Filename[MAX_STRING];       // filename for recording particle data
static double  *Particle_Data_Table = NULL;              // particle data table [mass/position/velocity]
#endif
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
void Validate(void)
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Validating test problem %d ...\n", TESTPROB_ID );


// errors
#  if ( MODEL != ELBDM )
   Aux_Error( ERROR_INFO, "MODEL != ELBDM !!\n" );
#  endif

#  ifndef GRAVITY
   Aux_Error( ERROR_INFO, "GRAVITY must be enabled !!\n" );
#  endif

#  ifdef COMOVING
   Aux_Error( ERROR_INFO, "COMOVING must be disabled !!\n" );
#  endif

#  ifdef GRAVITY
   if ( OPT__BC_POT != BC_POT_ISOLATED )
      Aux_Error( ERROR_INFO, "must adopt isolated BC for gravity --> reset OPT__BC_POT !!\n" );
#  endif


// warnings
   if ( MPI_Rank == 0 )
   {
      if ( !OPT__INIT_RESTRICT )
         Aux_Message( stderr, "WARNING : it's recommended to enable OPT__INIT_RESTRICT !!\n" );
   }


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Validating test problem %d ... done\n", TESTPROB_ID );

} // FUNCTION : Validate



#if ( MODEL == ELBDM  &&  defined GRAVITY )
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
// ReadPara->Add( "KEY_IN_THE_FILE",           &VARIABLE,                   DEFAULT,       MIN,              MAX               );
// ********************************************************************************************************************************
   ReadPara->Add( "Soliton_N",                 &Soliton_N,                 -1,             1,                NoMax_int         );
   ReadPara->Add( "Soliton_RSeed",             &Soliton_RSeed,              0,             NoMin_int,        NoMax_int         );
   ReadPara->Add( "Soliton_Distortion_Flag",   &Soliton_Distortion,         0,             NoMin_int,        NoMax_int         );
   ReadPara->Add( "Soliton_CoreRadiusAll",     &Soliton_CoreRadiusAll,      NoDef_double,  NoMin_double,     NoMax_double      );
   ReadPara->Add( "Soliton_EmptyRegion",       &Soliton_EmptyRegion,        0.0,           NoMin_double,     NoMax_double      );
   ReadPara->Add( "Soliton_CM_MaxR",           &Fake_Soliton_CM_MaxR,       NoMax_double,  Eps_double,       NoMax_double      );
   ReadPara->Add( "Soliton_CM_TolErrR",        &Fake_Soliton_CM_TolErrR,    0.0,           Eps_double,       NoMax_double       );
   ReadPara->Add( "Soliton_DensProf_Filename",  Soliton_DensProf_Filename,  Useless_str,   Useless_str,      Useless_str       );
   ReadPara->Add( "Soliton_Distortion_Filename", Soliton_Distortion_Filename, Useless_str,   Useless_str,      Useless_str       );
#ifdef PARTICLE
   if (amr->Par->NPar_Active_AllRank>0)
   {
      ReadPara->Add( "Particle_Data_Filename",     Particle_Data_Filename,     Useless_str,   Useless_str,      Useless_str       );
      ReadPara->Add( "Particle_Log_Filename",      Particle_Log_Filename,      Useless_str,   Useless_str,      Useless_str       );
   }
#endif

   ReadPara->Read( FileName );

   delete ReadPara;

// (1-2) set the default values

// (1-3) check the runtime parameters
   if ( Soliton_RSeed >= 0  &&  Soliton_EmptyRegion < 0.0 )
      Aux_Error( ERROR_INFO, "Soliton_EmptyRegion (%14.7e) < 0.0 !!\n", Soliton_EmptyRegion );

   if ( Soliton_CoreRadiusAll == NoDef_double )
      Aux_Error( ERROR_INFO, "Runtime parameter \"Soliton_CoreRadiusAll\" is not set !!\n" );
   
// (2) set the problem-specific derived parameters
// (2-1) allocate memory
   Soliton_CoreRadius       = new double [Soliton_N];
   Soliton_Center           = new double [Soliton_N][3];
   Soliton_ScaleL           = new double [Soliton_N];
   Soliton_ScaleD           = new double [Soliton_N];

// (2-2) soliton core radii
   if ( Soliton_CoreRadiusAll > 0.0 )
   {
      for (int t=0; t<Soliton_N; t++)  Soliton_CoreRadius[t] = Soliton_CoreRadiusAll;
   }

   else
   {
//    for Soliton_CoreRadiusAll <= 0.0, comment out the following line and hard code the core radius of each soliton
      Aux_Error( ERROR_INFO, "for Soliton_CoreRadiusAll <= 0.0, please comment out this error check and hard code "
                             "the core radius of each soliton !!\n" );
//    for (int t=0; t<Soliton_N; t++)  Soliton_CoreRadius[t] = XXX;
   }

// (2-3) soliton centers
   if ( Soliton_RSeed >= 0 )
   {
      const double Coord_Min[3] = { Soliton_EmptyRegion, Soliton_EmptyRegion, Soliton_EmptyRegion };
      const double Coord_Max[3] = { amr->BoxSize[0] - Soliton_EmptyRegion,
                                    amr->BoxSize[1] - Soliton_EmptyRegion,
                                    amr->BoxSize[2] - Soliton_EmptyRegion };
      srand( Soliton_RSeed );

      for (int t=0; t<Soliton_N; t++)
      for (int d=0; d<3; d++)
         Soliton_Center[t][d] = ( (double)rand()/RAND_MAX )*(Coord_Max[d]-Coord_Min[d]) + Coord_Min[d];
   }

   else
   {
      if ( Soliton_N == 1 )
      {
         for (int d=0; d<3; d++)    Soliton_Center[0][d] = 0.5*amr->BoxSize[d];
      }

      else
      {
//       for Soliton_RSeed<0, comment out the following line and hard code the center of each soliton
         Aux_Error( ERROR_INFO, "for Soliton_RSeed < 0 and Soliton_N > 1, please comment out this error check and hard code "
                                "the center of each soliton !!\n" );

         /*
         for (int t=0; t<Soliton_N; t++)
         for (int d=0; d<3; d++)          Soliton_Center[t][d] = XXX;
         */
      }
   } // if ( Soliton_RSeed >= 0 ) ... else ...


// (3) load the reference soliton density profile and evaluate the scale factors
   if ( OPT__INIT != INIT_BY_RESTART )
   {
//    load the reference profile
      const bool RowMajor_No_soliton_prof             = false;    // load data into the column-major order
      const bool AllocMem_Yes_soliton_prof            = true;     // allocate memory for Soliton_DensProf
      const int  NCol_soliton_prof                    = 2;        // total number of columns to load
      const int  Col_soliton_prof[NCol_soliton_prof]  = {0, 1};   // target columns: (radius, density)
      Soliton_DensProf_NBin = Aux_LoadTable( Soliton_DensProf, Soliton_DensProf_Filename, NCol_soliton_prof, Col_soliton_prof, RowMajor_No_soliton_prof, AllocMem_Yes_soliton_prof );
//    load the soliton distortion ratio from file
      if (Soliton_Distortion==1)
      {
         const bool RowMajor_No_soliton_ratio             = false;    // load data into the column-major order
         const bool AllocMem_Yes_soliton_ratio            = true;     // allocate memory for Soliton_DensProf
         const int  NCol_soliton_ratio                    = 3;        // total number of columns to load
         const int  Col_soliton_ratio[NCol_soliton_ratio] = {0, 1, 2};   // target columns: (radius, density)
          
         const int Soliton_Number_Check = Aux_LoadTable( Soliton_Distortion_Ratio, Soliton_Distortion_Filename, NCol_soliton_ratio, Col_soliton_ratio, RowMajor_No_soliton_ratio, AllocMem_Yes_soliton_ratio );
         if (Soliton_Number_Check!=Soliton_N)
            Aux_Error( ERROR_INFO, "Soliton_Number_Check (%d) != Soliton_N (%d) !! Exit!!", Soliton_Number_Check, Soliton_N);
      }
      
#ifdef PARTICLE
      if (amr->Par->NPar_Active_AllRank>0)
      {
         const bool RowMajor_No_particle_data            = false;                 // load data into the column-major order
         const bool AllocMem_Yes_particle_data           = true;                  // allocate memory for Soliton_DensProf
         const int NCol_particle_data                    = 7;                     // total number of columns to load for particle data
         const int Col_particle_data[NCol_particle_data] = {0, 1, 2, 3, 4, 5, 6}; // target columns: (mass, position_x, position_y, position_z, velocity_x, velocity_y, velocity_z)
         if (amr->Par->NPar_Active_AllRank>0) 
            NPar_AllRank_Check = Aux_LoadTable( Particle_Data_Table, Particle_Data_Filename, NCol_particle_data, Col_particle_data, RowMajor_No_particle_data, AllocMem_Yes_particle_data );
      }
#endif

//    get the core radius of the reference profile
      const double *RadiusRef = Soliton_DensProf + 0*Soliton_DensProf_NBin;
      const double *DensRef   = Soliton_DensProf + 1*Soliton_DensProf_NBin;
      const double  DensCore  = 0.5*DensRef[0];   // define core radius as the half-density radius

      double CoreRadiusRef = NULL_REAL;

      for (int b=1; b<Soliton_DensProf_NBin-1; b++)
      {
         if ( DensRef[b] >= DensCore  &&  DensRef[b+1] <= DensCore )
         {
            CoreRadiusRef = 0.5*( RadiusRef[b] + RadiusRef[b+1] );
            break;
         }
      }

      if ( CoreRadiusRef == NULL_REAL )
         Aux_Error( ERROR_INFO, "cannot determine the reference core radius !!\n" );


//    evaluate the scale factors of each soliton
      for (int t=0; t<Soliton_N; t++)
      {
         Soliton_ScaleL[t] = Soliton_CoreRadius[t] / CoreRadiusRef;
         Soliton_ScaleD[t] = 1.0 / ( 4.0*M_PI*NEWTON_G*SQR(ELBDM_ETA)*POW4(Soliton_ScaleL[t]) );
      }
   } // if ( OPT__INIT != INIT_BY_RESTART )


// (4) reset other general-purpose parameters
//     --> a helper macro PRINT_WARNING is defined in TestProb.h
   const long   End_Step_Default = __INT_MAX__;
   const double End_T_Default    = 1.0e3;

   if ( END_STEP < 0 ) {
      END_STEP = End_Step_Default;
      PRINT_WARNING( "END_STEP", END_STEP, FORMAT_LONG );
   }

   if ( END_T < 0.0 ) {
      END_T = End_T_Default;
      PRINT_WARNING( "END_T", END_T, FORMAT_REAL );
   }


// (5) make a note
   if ( MPI_Rank == 0 )
   {
      Aux_Message( stdout, "======================================================================================\n" );
      Aux_Message( stdout, "  test problem ID                           = %d\n",     TESTPROB_ID                      );
      Aux_Message( stdout, "  total number of solitons                  = %d\n",     Soliton_N                        );
      Aux_Message( stdout, "  random seed for setting the center coord. = %d\n",     Soliton_RSeed                    );
      if ( OPT__INIT != INIT_BY_RESTART )
          Aux_Message( stdout, "  number of bins of the density profile     = %d\n",     Soliton_DensProf_NBin        );
      Aux_Message( stdout, "  soliton distortion flag                   = %d\n",     Soliton_Distortion               );
      Aux_Message( stdout, "  size of the soliton-free zone             = %13.7e\n", Soliton_EmptyRegion              );
      Aux_Message( stdout, "  density profile filename                  = %s\n",     Soliton_DensProf_Filename        );
      Aux_Message( stdout, "  soliton distortion ratio filename         = %s\n",     Soliton_Distortion_Filename      );
#ifdef PARTICLE
      if (amr->Par->NPar_Active_AllRank>0)
      {
          Aux_Message( stdout, "  particle data filename                    = %s\n",     Particle_Data_Filename           );
          Aux_Message( stdout, "  particle log filename                     = %s\n",     Particle_Log_Filename            );
      }
#endif
      
      Aux_Message( stdout, "\n" );
      Aux_Message( stdout, "  Soliton info:\n" );
      if (Soliton_Distortion==1)
      {
         if ( OPT__INIT != INIT_BY_RESTART )
         {
             Aux_Message( stdout, "  %7s  %13s  %13s  %13s  %13s  %13s  %13s  %13s  %13s  %13s\n",
                          "ID", "CoreRadius", "ScaleL", "ScaleD", "Center_X", "Center_Y", "Center_Z" ,
                          "Distortion_Ratio_X", "Distortion_Ratio_y", "Distortion_z");
             for (int t=0; t<Soliton_N; t++)
             Aux_Message( stdout, "  %7d  %13.6e  %13.6e  %13.6e  %13.6e  %13.6e  %13.6e  %13.6e  %13.6e  %13.6e\n",
                          t, Soliton_CoreRadius[t], Soliton_ScaleL[t], Soliton_ScaleD[t],
                          Soliton_Center[t][0], Soliton_Center[t][1], Soliton_Center[t][2],
                          Soliton_Distortion_Ratio[0], Soliton_Distortion_Ratio[1], Soliton_Distortion_Ratio[2]);
         }
      }
      else
      {
         if ( OPT__INIT != INIT_BY_RESTART )
         {
             Aux_Message( stdout, "  %7s  %13s  %13s  %13s  %13s  %13s  %13s\n",
                          "ID", "CoreRadius", "ScaleL", "ScaleD", "Center_X", "Center_Y", "Center_Z");
             for (int t=0; t<Soliton_N; t++)
             Aux_Message( stdout, "  %7d  %13.6e  %13.6e  %13.6e  %13.6e  %13.6e  %13.6e\n",
                          t, Soliton_CoreRadius[t], Soliton_ScaleL[t], Soliton_ScaleD[t],
                          Soliton_Center[t][0], Soliton_Center[t][1], Soliton_Center[t][2]);
         }

      }
      Aux_Message( stdout, "======================================================================================\n" );
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

   const double *Table_Radius  = Soliton_DensProf + 0*Soliton_DensProf_NBin;  // radius
   const double *Table_Density = Soliton_DensProf + 1*Soliton_DensProf_NBin;  // density

   double r_tar, r_ref, dens_ref;


// initialize density as zero since there may be multiple solitons
   fluid[DENS] = 0.0;

// loop over all solitons to get the total density
   for (int t=0; t<Soliton_N; t++)
   {
      if (Soliton_Distortion==1)
         r_tar = sqrt( SQR((x-Soliton_Center[t][0])/Soliton_Distortion_Ratio[0]) + 
                       SQR((y-Soliton_Center[t][1])/Soliton_Distortion_Ratio[1]) + 
                       SQR((z-Soliton_Center[t][2])/Soliton_Distortion_Ratio[2]) );
      else
         r_tar = sqrt( SQR(x-Soliton_Center[t][0]) + SQR(y-Soliton_Center[t][1]) + SQR(z-Soliton_Center[t][2]) );

//    rescale radius (target radius --> reference radius)
      r_ref = r_tar / Soliton_ScaleL[t];

//    linear interpolation
      dens_ref = Mis_InterpolateFromTable( Soliton_DensProf_NBin, Table_Radius, Table_Density, r_ref );

      if ( dens_ref == NULL_REAL )
      {
         if      ( r_ref <  Table_Radius[0] )
            dens_ref = Table_Density[0];

         else if ( r_ref >= Table_Radius[Soliton_DensProf_NBin-1] )
            dens_ref = Table_Density[Soliton_DensProf_NBin-1];

         else
            Aux_Error( ERROR_INFO, "interpolation failed at radius %13.7e (min/max radius = %13.7e/%13.7e) !!\n",
                       r_ref, Table_Radius[0], Table_Radius[Soliton_DensProf_NBin-1] );
      }

//    rescale density (reference density --> target density) and add to the fluid array
      fluid[DENS] += dens_ref*Soliton_ScaleD[t];
   } // for (int t=0; t<Soliton_N; t++)


// set the real and imaginary parts
   fluid[REAL] = sqrt( fluid[DENS] );
   fluid[IMAG] = 0.0;                  // imaginary part is always zero --> no initial velocity

} // FUNCTION : SetGridIC


#ifdef PARTICLE
//-------------------------------------------------------------------------------------------------------
// Function    :  Par_Init_ByFunction_Toy_Model
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
//                ParPosX/Y/Z   : Particle position array with the size of NPar_ThisRank
//                ParVelX/Y/Z   : Particle velocity array with the size of NPar_ThisRank
//                ParTime       : Particle time     array with the size of NPar_ThisRank
//                AllAttribute  : Pointer array for all particle attributes
//                                --> Dimension = [PAR_NATT_TOTAL][NPar_ThisRank]
//                                --> Use the attribute indices defined in Field.h (e.g., Idx_ParCreTime)
//                                    to access the data
//
// Return      :  ParMass, ParPosX/Y/Z, ParVelX/Y/Z, ParTime, AllAttribute
//-------------------------------------------------------------------------------------------------------
void Par_Init_ByFunction_Toy_Model( const long NPar_ThisRank, const long NPar_AllRank,
                                    real *ParMass, real *ParPosX, real *ParPosY, real *ParPosZ,
                                    real *ParVelX, real *ParVelY, real *ParVelZ, real *ParTime,
                                    real *AllAttribute[PAR_NATT_TOTAL] )
{

   if (amr->Par->NPar_Active_AllRank<=0)
      return;

   if ( MPI_Rank == 0 )
   {
       Aux_Message( stdout, "%s ...\n", __FUNCTION__ );
       if ( NPar_AllRank_Check!=NPar_AllRank )
          Aux_Error( ERROR_INFO, "NPar_AllRank_check(%ld) != NPar_AllRank(%ld) !!\n" );
   }


   real *Mass_AllRank   = NULL;
   real *Pos_AllRank[3] = { NULL, NULL, NULL };
   real *Vel_AllRank[3] = { NULL, NULL, NULL };

// only the master rank will construct the initial condition
   if ( MPI_Rank == 0 )
   {
      double  TotM;
      const double *Mass_table        = Particle_Data_Table;
      const double *Position_table[3];
      const double *Velocity_table[3];

      Mass_AllRank = new real [NPar_AllRank];
      for (int d=0; d<3; d++)
      {
         Position_table[d] = Particle_Data_Table+(1+d)*NPar_AllRank;
         Velocity_table[d] = Particle_Data_Table+(4+d)*NPar_AllRank;
         Pos_AllRank[d] = new real [NPar_AllRank];
         Vel_AllRank[d] = new real [NPar_AllRank];
      }


//    set particle attributes
      for (long p=0; p<NPar_AllRank; p++)
      {
//       mass
         Mass_AllRank[p] = Mass_table[p];
         TotM += Mass_AllRank[p];

//       position
         for (int d=0; d<3; d++)    Pos_AllRank[d][p] = Position_table[d][p];

//       check periodicity
         for (int d=0; d<3; d++)
         {
            if ( OPT__BC_FLU[d*2] == BC_FLU_PERIODIC )
               Pos_AllRank[d][p] = FMOD( Pos_AllRank[d][p]+(real)amr->BoxSize[d], (real)amr->BoxSize[d] );
         }

//       velocity
         for (int d=0; d<3; d++)    Vel_AllRank[d][p] = Velocity_table[d][p];

      } // for (long p=0; p<NPar_AllRank; p++)

      Aux_Message( stdout, "=====================================================================================================\n" );
      Aux_Message( stdout, "Total mass = %13.7e\n",  TotM );
      for (long p=0; p<NPar_AllRank; p++)
           Aux_Message( stdout, "Pisition for particle #%ld is ( %13.7e,%13.7e,%13.7e )\n", p, Pos_AllRank[0][p], Pos_AllRank[1][p], Pos_AllRank[2][p] );
      for (long p=0; p<NPar_AllRank; p++)
           Aux_Message( stdout, "Velocity for particle #%ld is ( %13.7e,%13.7e,%13.7e )\n", p, Vel_AllRank[0][p], Vel_AllRank[1][p], Vel_AllRank[2][p] );
      Aux_Message( stdout, "=====================================================================================================\n" );

   } // if ( MPI_Rank == 0 )


// synchronize all particles to the physical time on the base level
   for (long p=0; p<NPar_ThisRank; p++)   ParTime[p] = Time[0];


// get the number of particles in each rank and set the corresponding offsets
   if ( NPar_AllRank > (long)__INT_MAX__ )
      Aux_Error( ERROR_INFO, "NPar_Active_AllRank (%ld) exceeds the maximum integer (%ld) --> MPI will likely fail !!\n",
                 NPar_AllRank, (long)__INT_MAX__ );

   int NSend[MPI_NRank], SendDisp[MPI_NRank];
   int NPar_ThisRank_int = NPar_ThisRank;    // (i) convert to "int" and (ii) remove the "const" declaration
                                             // --> (ii) is necessary for OpenMPI version < 1.7

   MPI_Gather( &NPar_ThisRank_int, 1, MPI_INT, NSend, 1, MPI_INT, 0, MPI_COMM_WORLD );

   if ( MPI_Rank == 0 )
   {
      SendDisp[0] = 0;
      for (int r=1; r<MPI_NRank; r++)  SendDisp[r] = SendDisp[r-1] + NSend[r-1];
   }


// send particle attributes from the master rank to all ranks
   real *Mass   =   ParMass;
   real *Pos[3] = { ParPosX, ParPosY, ParPosZ };
   real *Vel[3] = { ParVelX, ParVelY, ParVelZ };

#  ifdef FLOAT8
   MPI_Scatterv( Mass_AllRank, NSend, SendDisp, MPI_DOUBLE, Mass, NPar_ThisRank, MPI_DOUBLE, 0, MPI_COMM_WORLD );

   for (int d=0; d<3; d++)
   {
      MPI_Scatterv( Pos_AllRank[d], NSend, SendDisp, MPI_DOUBLE, Pos[d], NPar_ThisRank, MPI_DOUBLE, 0, MPI_COMM_WORLD );
      MPI_Scatterv( Vel_AllRank[d], NSend, SendDisp, MPI_DOUBLE, Vel[d], NPar_ThisRank, MPI_DOUBLE, 0, MPI_COMM_WORLD );
   }

#  else
   MPI_Scatterv( Mass_AllRank, NSend, SendDisp, MPI_FLOAT,  Mass, NPar_ThisRank, MPI_FLOAT,  0, MPI_COMM_WORLD );

   for (int d=0; d<3; d++)
   {
      MPI_Scatterv( Pos_AllRank[d], NSend, SendDisp, MPI_FLOAT,  Pos[d], NPar_ThisRank, MPI_FLOAT,  0, MPI_COMM_WORLD );
      MPI_Scatterv( Vel_AllRank[d], NSend, SendDisp, MPI_FLOAT,  Vel[d], NPar_ThisRank, MPI_FLOAT,  0, MPI_COMM_WORLD );
   }
#  endif


   if ( MPI_Rank == 0 )
   {
      delete [] Mass_AllRank;

      for (int d=0; d<3; d++)
      {
         delete [] Pos_AllRank[d];
         delete [] Vel_AllRank[d];
      }
   }


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : Par_Init_ByFunction_Toy_Model


//-------------------------------------------------------------------------------------------------------
// Function    :  Record_Particle_Data
// Description :  Output the particle position and velocity
//
// Parameter   :  FileName : Output file name
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
static void Record_Particle_Data( char *FileName )
{

// check
//   if ( MPI_Rank == 0  &&  Aux_CheckFileExist(FileName) )
//      Aux_Message( stderr, "WARNING : file \"%s\" already exists and will be overwritten !!\n", FileName );

   FILE *File;
// header
   if ( MPI_Rank == 0 )
   {
      if ( first_run_flag )
      {
          if ( !Aux_CheckFileExist(FileName) )
              File = fopen( FileName, "w" );
          else
              File = fopen( FileName, "a" );

          fprintf( File, "#Time                    Step                    Active Particles   ");

          for (int v=0; v<PAR_NATT_TOTAL; v++)
              fprintf( File, "  %*s", (v==0)?20:21, ParAttLabel[v] );
          fprintf( File, "\n" );
          first_run_flag = false;
      }
      else
          File = fopen( FileName, "a" );
      fprintf( File, "%20.14e    %13ld    %13ld          ", Time[0], Step, amr->Par->NPar_Active_AllRank );
      fclose( File );
   }

// data
   for (int TargetMPIRank=0; TargetMPIRank<MPI_NRank; TargetMPIRank++)
   {
      if ( MPI_Rank == TargetMPIRank )
      {
         File = fopen( FileName, "a" );

         for (long p=0; p<amr->Par->NPar_AcPlusInac; p++)
         {
//          skip inactive particles
            if ( amr->Par->Mass[p] < 0.0 )   continue;

            for (int v=0; v<PAR_NATT_TOTAL; v++)   fprintf( File, "  %21.14e", amr->Par->Attribute[v][p] );

            fprintf( File, "\n" );
         }

         fclose( File );
      } // if ( MPI_Rank == TargetMPIRank )

      MPI_Barrier( MPI_COMM_WORLD );
   } // for (int TargetMPIRank=0; TargetMPIRank<MPI_NRank; TargetMPIRank++)


} // FUNCTION : Record_Particle_Data


#endif // #ifdef PARTICLE


//-------------------------------------------------------------------------------------------------------
// Function    :  End_Soliton_Toy_Model
// Description :  Free memory before terminating the program
//
// Note        :  1. Linked to the function pointer "End_User_Ptr" to replace "End_User()"
//
// Parameter   :  None
//-------------------------------------------------------------------------------------------------------
void End_Soliton_Toy_Model()
{

   delete [] Soliton_DensProf;
   delete [] Soliton_CoreRadius;
   delete [] Soliton_Center;
   delete [] Soliton_ScaleL;
   delete [] Soliton_ScaleD;
#ifdef PARTICLE
   delete [] Particle_Data_Table;
#endif

} // FUNCTION : End_Soliton_Toy_model


//-------------------------------------------------------------------------------------------------------
// Function    :  BC_Toy_Model
// Description :  Set the extenral boundary condition
//
// Note        :  1. Linked to the function pointer "BC_User_Ptr"
//
// Parameter   :  fluid    : Fluid field to be set
//                x/y/z    : Physical coordinates
//                Time     : Physical time
//                lv       : Refinement level
//                AuxArray : Auxiliary array
//
// Return      :  fluid
//-------------------------------------------------------------------------------------------------------
void BC_Toy_Model( real fluid[], const double x, const double y, const double z, const double Time,
         const int lv, double AuxArray[] )
{

   fluid[REAL] = (real)0.0;
   fluid[IMAG] = (real)0.0;
   fluid[DENS] = (real)0.0;

} // FUNCTION : BC_Toy_Model


//-------------------------------------------------------------------------------------------------------
// Function    :  GetCenterOfMass
// Description :  Record the center of mass (CM)
//
// Note        :  1. Invoked by Record_EridanusII() recursively
//                2. Only include cells within CM_MaxR from CM_Old[] when updating CM
//
// Parameter   :  CM_Old[] : Previous CM
//                CM_New[] : New CM to be returned
//                CM_MaxR  : Maximum radius to compute CM
//
// Return      :  CM_New[]
//-------------------------------------------------------------------------------------------------------
void GetCenterOfMass( const double CM_Old[], double CM_New[], const double CM_MaxR )
{

   const double CM_MaxR2          = SQR( CM_MaxR );
   const double HalfBox[3]        = { 0.5*amr->BoxSize[0], 0.5*amr->BoxSize[1], 0.5*amr->BoxSize[2] };
   const bool   Periodic          = ( OPT__BC_FLU[0] == BC_FLU_PERIODIC );
   const bool   IntPhase_No       = false;
   const real   MinDens_No        = -1.0;
   const real   MinPres_No        = -1.0;
   const real   MinTemp_No        = -1.0;
   const bool   DE_Consistency_No = false;
#  ifdef PARTICLE
#  ifdef LOAD_BALANCE
   const bool   SibBufPatch       = true;
   const bool   FaSibBufPatch     = true;
#  else
   const bool   SibBufPatch       = NULL_BOOL;
   const bool   FaSibBufPatch     = NULL_BOOL;
#  endif
#  endif // #ifdef PARTICLE

   int   *PID0List = NULL;
   double M_ThisRank, MR_ThisRank[3], M_AllRank, MR_AllRank[3];
   real (*TotalDens)[PS1][PS1][PS1];

   M_ThisRank = 0.0;
   for (int d=0; d<3; d++)    MR_ThisRank[d] = 0.0;


   for (int lv=0; lv<NLEVEL; lv++)
   {
//    initialize the particle density array (rho_ext) and collect particles to the target level

//    get the total density on grids
      TotalDens = new real [ amr->NPatchComma[lv][1] ][PS1][PS1][PS1];
      PID0List  = new int  [ amr->NPatchComma[lv][1]/8 ];

      for (int PID0=0, t=0; PID0<amr->NPatchComma[lv][1]; PID0+=8, t++)    PID0List[t] = PID0;

      Prepare_PatchData( lv, Time[lv], TotalDens[0][0][0], NULL, 0, amr->NPatchComma[lv][1]/8, PID0List, _TOTAL_DENS, _NONE,
                         OPT__RHO_INT_SCHEME, INT_NONE, UNIT_PATCH, NSIDE_00, IntPhase_No, OPT__BC_FLU, BC_POT_NONE,
                         MinDens_No, MinPres_No, MinTemp_No, DE_Consistency_No );

      delete [] PID0List;


//    free memory for collecting particles from other ranks and levels, and free density arrays with ghost zones (rho_ext)
#     ifdef PARTICLE
      Par_CollectParticle2OneLevel_FreeMemory( lv, SibBufPatch, FaSibBufPatch );

      Prepare_PatchData_FreeParticleDensityArray( lv );
#     endif


//    calculate the center of mass
      const double dh = amr->dh[lv];
      const double dv = CUBE( dh );

      for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
      {
//       skip non-leaf patches
         if ( amr->patch[0][lv][PID]->son != -1 )  continue;

         const double x0 = amr->patch[0][lv][PID]->EdgeL[0] + 0.5*dh;
         const double y0 = amr->patch[0][lv][PID]->EdgeL[1] + 0.5*dh;
         const double z0 = amr->patch[0][lv][PID]->EdgeL[2] + 0.5*dh;

         double x, y, z, dx, dy, dz;

         for (int k=0; k<PS1; k++)  {  z = z0 + k*dh;  dz = z - CM_Old[2];
                                       if ( Periodic ) {
                                          if      ( dz > +HalfBox[2] )  {  z -= amr->BoxSize[2];  dz -= amr->BoxSize[2];  }
                                          else if ( dz < -HalfBox[2] )  {  z += amr->BoxSize[2];  dz += amr->BoxSize[2];  }
                                       }
         for (int j=0; j<PS1; j++)  {  y = y0 + j*dh;  dy = y - CM_Old[1];
                                       if ( Periodic ) {
                                          if      ( dy > +HalfBox[1] )  {  y -= amr->BoxSize[1];  dy -= amr->BoxSize[1];  }
                                          else if ( dy < -HalfBox[1] )  {  y += amr->BoxSize[1];  dy += amr->BoxSize[1];  }
                                       }
         for (int i=0; i<PS1; i++)  {  x = x0 + i*dh;  dx = x - CM_Old[0];
                                       if ( Periodic ) {
                                          if      ( dx > +HalfBox[0] )  {  x -= amr->BoxSize[0];  dx -= amr->BoxSize[0];  }
                                          else if ( dx < -HalfBox[0] )  {  x += amr->BoxSize[0];  dx += amr->BoxSize[0];  }
                                       }

//          only include cells within CM_MaxR
            const double R2 = SQR(dx) + SQR(dy) + SQR(dz);
            if ( R2 < CM_MaxR2 )
            {
               const double dm = TotalDens[PID][k][j][i]*dv;

               M_ThisRank     += dm;
               MR_ThisRank[0] += dm*x;
               MR_ThisRank[1] += dm*y;
               MR_ThisRank[2] += dm*z;
            }
         }}}
      } // for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)

      delete [] TotalDens;
   } // for (int lv=0; lv<NLEVEL; lv++)


// collect data from all ranks to calculate the CM
// --> note that all ranks will get CM_New[]
   MPI_Allreduce( &M_ThisRank, &M_AllRank, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
   MPI_Allreduce( MR_ThisRank, MR_AllRank, 3, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );

   for (int d=0; d<3; d++)    CM_New[d] = MR_AllRank[d] / M_AllRank;

// map the new CM back to the simulation domain
   if ( Periodic )
   for (int d=0; d<3; d++)
   {
      if      ( CM_New[d] >= amr->BoxSize[d] )  CM_New[d] -= amr->BoxSize[d];
      else if ( CM_New[d] < 0.0              )  CM_New[d] += amr->BoxSize[d];

   }

   for (int d=0; d<3; d++)
      if ( CM_New[d] >= amr->BoxSize[d]  ||  CM_New[d] < 0.0 )
         Aux_Error( ERROR_INFO, "CM_New[%d] = %14.7e lies outside the domain !!\n", d, CM_New[d] );

} // FUNCTION : GetCenterOfMass


//-------------------------------------------------------------------------------------------------------
// Function    :  Record_Fake_Soliton_CoM
// Description :  Record the maximum density and center coordinates
//
// Note        :  1. It will also record the real and imaginary parts associated with the maximum density
//                2. For the center coordinates, it will record the position of maximum density, minimum potential,
//                   and center-of-mass
//                3. Output filename is fixed to "Record__Center"
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Record_Fake_Soliton_CoM(void)
{

   const char filename_center  [] = "Record__Center";
   const int  CountMPI            = 10;

   double dens, max_dens_loc=-__DBL_MAX__, max_dens_pos_loc[3], real_loc, imag_loc;
   double pote, min_pote_loc=+__DBL_MAX__, min_pote_pos_loc[3];
   double send[CountMPI], (*recv)[CountMPI]=new double [MPI_NRank][CountMPI];

   const long   DensMode          = _DENS;
   const bool   IntPhase_No       = false;
   const real   MinDens_No        = -1.0;
   const real   MinPres_No        = -1.0;
   const real   MinTemp_No        = -1.0;
   const bool   DE_Consistency_No = false;
#  ifdef PRATICLE
#  ifdef LOAD_BALANCE
   const bool   SibBufPatch       = true;
   const bool   FaSibBufPatch     = true;
#  else
   const bool   SibBufPatch       = NULL_BOOL;
   const bool   FaSibBufPatch     = NULL_BOOL;
#  endif
#  endif // #ifdef PARTICLE


// collect local data
   for (int lv=0; lv<NLEVEL; lv++)
   {
//    initialize the particle density array (rho_ext) and collect particles to the target level

//    get the total density on grids
      real (*TotalDens)[PS1][PS1][PS1] = new real [ amr->NPatchComma[lv][1] ][PS1][PS1][PS1];
      int   *PID0List                  = new int  [ amr->NPatchComma[lv][1]/8 ];

      for (int PID0=0, t=0; PID0<amr->NPatchComma[lv][1]; PID0+=8, t++)    PID0List[t] = PID0;

      Prepare_PatchData( lv, Time[lv], TotalDens[0][0][0], NULL, 0, amr->NPatchComma[lv][1]/8, PID0List, DensMode, _NONE,
                         OPT__RHO_INT_SCHEME, INT_NONE, UNIT_PATCH, NSIDE_00, IntPhase_No, OPT__BC_FLU, BC_POT_NONE,
                         MinDens_No, MinPres_No, MinTemp_No, DE_Consistency_No );

      delete [] PID0List;


//    free memory for collecting particles from other ranks and levels, and free density arrays with ghost zones (rho_ext)

      for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
      {
//       skip non-leaf patches
         if ( amr->patch[0][lv][PID]->son != -1 )  continue;

         for (int k=0; k<PS1; k++)  {  const double z = amr->patch[0][lv][PID]->EdgeL[2] + (k+0.5)*amr->dh[lv];
         for (int j=0; j<PS1; j++)  {  const double y = amr->patch[0][lv][PID]->EdgeL[1] + (j+0.5)*amr->dh[lv];
         for (int i=0; i<PS1; i++)  {  const double x = amr->patch[0][lv][PID]->EdgeL[0] + (i+0.5)*amr->dh[lv];

//          dens = amr->patch[ amr->FluSg[lv] ][lv][PID]->fluid[DENS][k][j][i];
            dens = TotalDens[PID][k][j][i];
            pote = amr->patch[ amr->PotSg[lv] ][lv][PID]->pot[k][j][i];

            if ( dens > max_dens_loc )
            {
               max_dens_loc        = dens;
               real_loc            = amr->patch[ amr->FluSg[lv] ][lv][PID]->fluid[REAL][k][j][i];
               imag_loc            = amr->patch[ amr->FluSg[lv] ][lv][PID]->fluid[IMAG][k][j][i];
               max_dens_pos_loc[0] = x;
               max_dens_pos_loc[1] = y;
               max_dens_pos_loc[2] = z;
            }

            if ( pote < min_pote_loc )
            {
               min_pote_loc        = pote;
               min_pote_pos_loc[0] = x;
               min_pote_pos_loc[1] = y;
               min_pote_pos_loc[2] = z;
            }
         }}}
      } // for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)

      delete [] TotalDens;
   } // for (int lv=0; lv<NLEVEL; lv++)


// gather data to the root rank
   send[0] = max_dens_loc;
   send[1] = real_loc;
   send[2] = imag_loc;
   send[3] = max_dens_pos_loc[0];
   send[4] = max_dens_pos_loc[1];
   send[5] = max_dens_pos_loc[2];
   send[6] = min_pote_loc;
   send[7] = min_pote_pos_loc[0];
   send[8] = min_pote_pos_loc[1];
   send[9] = min_pote_pos_loc[2];

   MPI_Gather( send, CountMPI, MPI_DOUBLE, recv[0], CountMPI, MPI_DOUBLE, 0, MPI_COMM_WORLD );


// record the maximum density and center coordinates
   double max_dens      = -__DBL_MAX__;
   double min_pote      = +__DBL_MAX__;
   int    max_dens_rank = -1;
   int    min_pote_rank = -1;

   if ( MPI_Rank == 0 )
   {
      for (int r=0; r<MPI_NRank; r++)
      {
         if ( recv[r][0] > max_dens )
         {
            max_dens      = recv[r][0];
            max_dens_rank = r;
         }

         if ( recv[r][6] < min_pote )
         {
            min_pote      = recv[r][6];
            min_pote_rank = r;
         }
      }

      if ( max_dens_rank < 0  ||  max_dens_rank >= MPI_NRank )
         Aux_Error( ERROR_INFO, "incorrect max_dens_rank (%d) !!\n", max_dens_rank );

      if ( min_pote_rank < 0  ||  min_pote_rank >= MPI_NRank )
         Aux_Error( ERROR_INFO, "incorrect min_pote_rank (%d) !!\n", min_pote_rank );

      static bool FirstTime = true;

      if ( FirstTime )
      {
         if ( Aux_CheckFileExist(filename_center) )
            Aux_Message( stderr, "WARNING : file \"%s\" already exists !!\n", filename_center );
         else
         {
            FILE *file_center = fopen( filename_center, "w" );
            fprintf( file_center, "#%19s  %10s  %14s  %14s  %14s  %14s  %14s  %14s  %14s  %14s  %14s  %14s  %10s  %14s  %14s  %14s\n",
                     "Time", "Step", "Dens", "Real", "Imag", "Dens_x", "Dens_y", "Dens_z", "Pote", "Pote_x", "Pote_y", "Pote_z",
                     "NIter", "CM_h_x", "CM_h_y", "CM_h_z" );
            fclose( file_center );
         }

         FirstTime = false;
      }

      FILE *file_center = fopen( filename_center, "a" );
      fprintf( file_center, "%20.14e  %10ld  %14.7e  %14.7e  %14.7e  %14.7e  %14.7e  %14.7e  %14.7e  %14.7e  %14.7e  %14.7e",
               Time[0], Step, recv[max_dens_rank][0], recv[max_dens_rank][1], recv[max_dens_rank][2], recv[max_dens_rank][3],
                              recv[max_dens_rank][4], recv[max_dens_rank][5], recv[min_pote_rank][6], recv[min_pote_rank][7],
                              recv[min_pote_rank][8], recv[min_pote_rank][9] );
      fclose( file_center );
   } // if ( MPI_Rank == 0 )


// compute the center of mass until convergence
   const double TolErrR2 = SQR( Fake_Soliton_CM_TolErrR );
   const int    NIterMax = 10;

   double dR2, CM_Old[3], CM_New[3];
   int NIter = 0;

// set an initial guess by the peak density position
   if ( MPI_Rank == 0 )
      for (int d=0; d<3; d++)    CM_Old[d] = recv[max_dens_rank][3+d];

   MPI_Bcast( CM_Old, 3, MPI_DOUBLE, 0, MPI_COMM_WORLD );

   while ( true )
   {
      GetCenterOfMass( CM_Old, CM_New, Fake_Soliton_CM_MaxR );

      dR2 = SQR( CM_Old[0] - CM_New[0] )
          + SQR( CM_Old[1] - CM_New[1] )
          + SQR( CM_Old[2] - CM_New[2] );
      NIter ++;

      if ( dR2 <= TolErrR2  ||  NIter >= NIterMax )
         break;
      else
         memcpy( CM_Old, CM_New, sizeof(double)*3 );
   }

   if ( MPI_Rank == 0 )
   {
      if ( dR2 > TolErrR2 )
         Aux_Message( stderr, "WARNING : dR (%13.7e) > Soliton_CM_TolErrR (%13.7e) !!\n", sqrt(dR2), Fake_Soliton_CM_TolErrR );

      FILE *file_center = fopen( filename_center, "a" );
      fprintf( file_center, "  %10d  %14.7e  %14.7e  %14.7e\n", NIter, CM_New[0], CM_New[1], CM_New[2] );
      fclose( file_center );
   }

   delete [] recv;

} // FUNCTION : Record_Fake_Soliton_CoM


//-------------------------------------------------------------------------------------------------------
// Function    :  Record_COM_and_Particle_Pos
// Description :  Record center of mass and position of particles
//
// Note        :  1. It will call center of mass routine 
//                2. For the center coordinates, it will record the position of maximum density, minimum potential,
//                   and center-of-mass
//                3. Output filename is fixed to "Record__Center"
//                4. It will call record particle position routine
//                5. Output filename is 
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
//
void Record_COM_and_Particle_Data(void)
{
   Record_Fake_Soliton_CoM();
#ifdef PARTICLE
   if (amr->Par->NPar_Active_AllRank>0)
      Record_Particle_Data(Particle_Log_Filename);
#endif
}


#endif // #if ( MODEL == ELBDM  &&  defined GRAVITY )


//-------------------------------------------------------------------------------------------------------
// Function    :  Init_TestProb_ELBDM_Soliton
// Description :  Test problem initializer
//
// Note        :  None
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Init_TestProb_ELBDM_Soliton_Toy_Model()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );


// validate the compilation flags and runtime parameters
   Validate();


#  if ( MODEL == ELBDM  &&  defined GRAVITY )
// set the problem-specific runtime parameters
   SetParameter();


   Init_Function_User_Ptr  = SetGridIC;
# ifdef PARTICLE
   Par_Init_ByFunction_Ptr = Par_Init_ByFunction_Toy_Model;
#endif
   Aux_Record_User_Ptr     = Record_COM_and_Particle_Data;
   BC_User_Ptr             = BC_Toy_Model;
   End_User_Ptr            = End_Soliton_Toy_Model;
#  endif // #if ( MODEL == ELBDM  &&  defined GRAVITY )


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : Init_TestProb_ELBDM_Soliton
