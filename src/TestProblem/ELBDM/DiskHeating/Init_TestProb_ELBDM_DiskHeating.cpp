#include "GAMER.h"
#include "TestProb.h"

#ifdef SUPPORT_GSL
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#endif


static void SetGridIC( real fluid[], const double x, const double y, const double z, const double Time,
                       const int lv, double AuxArray[] );
static void BC( real Array[], const int ArraySize[], real fluid[], const int NVar_Flu,
                const int GhostSize, const int idx[], const double pos[], const double Time,
                const int lv, const int TFluVarIdxList[], double AuxArray[] );
static void End_DiskHeating();
static double Get_Dispersion( double r );
static double Halo_Density( double r );

// problem-specific global variables
// =======================================================================================
static RandomNumber_t *RNG = NULL;
       double  Cen[3];                    // center
static bool    AddFixedHalo;              // add a fixed halo, must enable OPT__FREEZE_FLUID
static bool    HaloUseTable;              // 0 = from analytical profile, 1 = from table
       double  m_22;                      // ELBDM particle mass, used for soliton profile of the fixed halo
       double  CoreRadius;                // soliton radius of the fixed halo (in kpc)
static double  Rho_0;                     // halo rho_0 (in 1.0e+10 Msun*kpc^-3)
static double  Rs;                        // halo Rs (in kpc)
static double  Alpha;                     // dimensionless, used for alpha-beta-gamma density profile of the fixed halo
static double  Beta;                      // dimensionless, used for alpha-beta-gamma density profile of the fixed halo
static double  Gamma;                     // dimensionless, used for alpha-beta-gamma density profile of the fixed halo
static char    DensTableFile[MAX_STRING]; // fixed halo density profile filename
static double *DensTable = NULL;          // density table, radius should be in kpc and density should be in g/cm^3
static int     DensTable_Nbin;            // number of bins of density table
static bool    AddParWhenRestart;         // add a new disk to an existing snapshot, must enable OPT__RESTART_RESET
static bool    AddParWhenRestartByFile;   // add a new disk via DiskHeatingParticleIC
static long    AddParWhenRestartNPar;     // particle number of the new disk
static int     NewDisk_RSeed;             // random seed for setting new disk particle position and velocity
static double  Disk_Mass;                 // thin disk total mass (code unit)
static double  Disk_R;                    // thin disk scale radius (code unit)
static char    DispTableFile[MAX_STRING]; // velocity dispersion table filename used for thin disk
static double *DispTable = NULL;          // velocity dispersion table, radius should be in kpc and dispersion should be in km/s
static int     DispTable_Nbin;            // number of bins of velocity dispersion table
// =======================================================================================

#ifdef PARTICLE
void Par_Init_ByFunction_DiskHeating( const long NPar_ThisRank, const long NPar_AllRank,
                                      real_par *ParMass, real_par *ParPosX, real_par *ParPosY, real_par *ParPosZ,
                                      real_par *ParVelX, real_par *ParVelY, real_par *ParVelZ, real_par *ParTime,
                                      long_par *ParType, real_par *AllAttributeFlt[PAR_NATT_FLT_TOTAL],
                                      long_par *AllAttributeInt[PAR_NATT_INT_TOTAL] );
static void Init_NewDiskRestart();
static void Init_NewDiskVelocity();
#endif
#ifdef GRAVITY
void Init_ExtPot_Soliton();
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
#  if ( MODEL != ELBDM )
   Aux_Error( ERROR_INFO, "MODEL != ELBDM !!\n" );
#  endif

#  ifndef GRAVITY
   Aux_Error( ERROR_INFO, "GRAVITY must be enabled !!\n" );
#  endif

#  ifndef PARTICLE
   Aux_Error( ERROR_INFO, "PARTICLE must be enabled !!\n" );
#  endif

#  ifdef COMOVING
   Aux_Error( ERROR_INFO, "COMOVING must be disabled !!\n" );
#  endif


// warnings
   if ( MPI_Rank == 0 )
   {
      if ( FLAG_BUFFER_SIZE < 5 )
         Aux_Message( stderr, "WARNING : it's recommended to set FLAG_BUFFER_SIZE >= 5 for this test !!\n" );

      if ( !OPT__RECORD_CENTER )
         Aux_Message( stderr, "WARNING : it's recommended to set OPT__RECORD_CENTER = 1 for this test !!\n" );
   } // if ( MPI_Rank == 0 )


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Validating test problem %d ... done\n", TESTPROB_ID );

} // FUNCTION : Validate



#if ( MODEL == ELBDM )
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
//                3. Must NOT call any EoS routine here since it hasn't been initialized at this point
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
// ReadPara->Add( "KEY_IN_THE_FILE",         &VARIABLE,                DEFAULT,       MIN,              MAX               );
// ********************************************************************************************************************************
   ReadPara->Add( "CenX",                    &Cen[0],                  NoDef_double,  NoMin_double,     NoMax_double      );
   ReadPara->Add( "CenY",                    &Cen[1],                  NoDef_double,  NoMin_double,     NoMax_double      );
   ReadPara->Add( "CenZ",                    &Cen[2],                  NoDef_double,  NoMin_double,     NoMax_double      );
   ReadPara->Add( "AddFixedHalo",            &AddFixedHalo,            false,         Useless_bool,     Useless_bool      );
   ReadPara->Add( "HaloUseTable",            &HaloUseTable,            false,         Useless_bool,     Useless_bool      );
   ReadPara->Add( "m_22",                    &m_22,                    0.4,           Eps_double,       NoMax_double      );
   ReadPara->Add( "CoreRadius",              &CoreRadius,              1.0,           Eps_double,       NoMax_double      );
   ReadPara->Add( "Rho_0",                   &Rho_0,                   1.0,           Eps_double,       NoMax_double      );
   ReadPara->Add( "Rs",                      &Rs,                      1.0,           Eps_double,       NoMax_double      );
   ReadPara->Add( "Alpha",                   &Alpha,                   1.0,           Eps_double,       NoMax_double      );
   ReadPara->Add( "Beta",                    &Beta,                    1.0,           Eps_double,       NoMax_double      );
   ReadPara->Add( "Gamma",                   &Gamma,                   1.0,           Eps_double,       NoMax_double      );
   ReadPara->Add( "DensTableFile",           DensTableFile,            NoDef_str,     Useless_str,      Useless_str       );
   ReadPara->Add( "AddParWhenRestart",       &AddParWhenRestart,       false,         Useless_bool,     Useless_bool      );
   ReadPara->Add( "AddParWhenRestartByFile", &AddParWhenRestartByFile, true,          Useless_bool,     Useless_bool      );
   ReadPara->Add( "AddParWhenRestartNPar",   &AddParWhenRestartNPar,   (long)0,       (long)0,          NoMax_long        );
   ReadPara->Add( "NewDisk_RSeed",           &NewDisk_RSeed,           1002,          0,                NoMax_int         );
   ReadPara->Add( "Disk_Mass",               &Disk_Mass,               1.0,           Eps_double,       NoMax_double      );
   ReadPara->Add( "Disk_R",                  &Disk_R,                  1.0,           Eps_double,       NoMax_double      );
   ReadPara->Add( "DispTableFile",           DispTableFile,            NoDef_str,     Useless_str,      Useless_str       );


   ReadPara->Read( FileName );

   delete ReadPara;

// (1-2) set the default values

// (1-3) check the runtime parameters


// (2) set the problem-specific derived parameters
// use density table as background fixed halo profile
   if ( HaloUseTable == 1 ) {

//    read the density table
      const bool RowMajor_No  =  false;    // load data into the column-major order
      const bool AllocMem_Yes =  true;     // allocate memory for DensTable
      const int  NCol         =  2;        // total number of columns to load
      const int  Col[NCol]    = {0, 1};    // target columns: (radius, density)

      DensTable_Nbin = Aux_LoadTable( DensTable, DensTableFile, NCol, Col, RowMajor_No, AllocMem_Yes );

      double *DensTable_r = DensTable + 0*DensTable_Nbin;
      double *DensTable_d = DensTable + 1*DensTable_Nbin;

//    convert to code unit and log-log scale (input units of radius and density should be fixed to kpc and g/cm^3 respectively)
      for (int b=0; b<DensTable_Nbin; b++)
      {
         DensTable_r[b] = log( DensTable_r[b]*Const_kpc/UNIT_L );
         DensTable_d[b] = log( DensTable_d[b]*CUBE(UNIT_L)/UNIT_M );
      }
   }


// (3) reset other general-purpose parameters
//     --> a helper macro PRINT_WARNING is defined in TestProb.h
   const long   End_Step_Default = 2500;
   const double End_T_Default    = 2.5e-1;

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
      Aux_Message( stdout, "  test problem ID         = %d\n",     TESTPROB_ID             );
      Aux_Message( stdout, "  CenX                    = %13.7e\n", Cen[0]                  );
      Aux_Message( stdout, "  CenY                    = %13.7e\n", Cen[1]                  );
      Aux_Message( stdout, "  Cenz                    = %13.7e\n", Cen[2]                  );
      Aux_Message( stdout, "  AddFixedHalo            = %d\n",     AddFixedHalo            );
      Aux_Message( stdout, "  HaloUseTable            = %d\n",     HaloUseTable            );
      Aux_Message( stdout, "  m_22                    = %13.7e\n", m_22                    );
      Aux_Message( stdout, "  CoreRadius              = %13.7e\n", CoreRadius              );
      Aux_Message( stdout, "  Rho_0                   = %13.7e\n", Rho_0                   );
      Aux_Message( stdout, "  Rs                      = %13.7e\n", Rs                      );
      Aux_Message( stdout, "  Alpha                   = %13.7e\n", Alpha                   );
      Aux_Message( stdout, "  Beta                    = %13.7e\n", Beta                    );
      Aux_Message( stdout, "  Gamma                   = %13.7e\n", Gamma                   );
      Aux_Message( stdout, "  DensTableFile           = %s\n",     DensTableFile           );
      Aux_Message( stdout, "  AddParWhenRestart       = %d\n",     AddParWhenRestart       );
      Aux_Message( stdout, "  AddParWhenRestartByFile = %d\n",     AddParWhenRestartByFile );
      Aux_Message( stdout, "  AddParWhenRestartNPar   = %d\n",     AddParWhenRestartNPar   );
      Aux_Message( stdout, "  NewDisk_RSeed           = %d\n",     NewDisk_RSeed           );
      Aux_Message( stdout, "  Disk_Mass               = %13.7e\n", Disk_Mass               );
      Aux_Message( stdout, "  Disk_R                  = %13.7e\n", Disk_R                  );
      Aux_Message( stdout, "  DispTableFile           = %s\n",     DispTableFile           );
      Aux_Message( stdout, "=============================================================================\n" );
   }


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Setting runtime parameters ... done\n" );

} // FUNCTION : SetParameter



//-------------------------------------------------------------------------------------------------------
// Function    :  Halo_Density
// Description :  alpha-beta-gamma density profile for fixed background halo with soliton as function of radius
//
// Note        :  1. r should have unit in kpc
//                2. Returned density is in 1.0e10*Msun/kpc^3
//-------------------------------------------------------------------------------------------------------
double Halo_Density( double r )
{

   double rho_halo = Rho_0 / pow( r/Rs, Alpha ) / pow( 1.0 + pow(r/Rs,Beta), (Gamma-Alpha)/Beta );
   if ( fabs(rho_halo) <  __FLT_MIN__ )   rho_halo = 0.0;

   double rho_soliton = 0.0019 / m_22 / m_22 * pow( 1.0/CoreRadius/pow(1.0 + 0.091*pow(r/CoreRadius,2.0), 2.0), 4.0 );
   if ( fabs(rho_soliton) <  __FLT_MIN__ )   rho_soliton = 0.0;

   double rho_max = 0.0019 / m_22 / m_22 * pow( 1.0 / CoreRadius, 4.0 );

   if ( (rho_halo + rho_soliton) > rho_max )    return rho_max;
   else return (rho_halo + rho_soliton);

} // FUNCTION : Halo_Density



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

// set the output array
   double dens = 0.0;

   if ( AddFixedHalo && OPT__FREEZE_FLUID )  // add a fixed halo at the background
   {
      const double dx = x - Cen[0];
      const double dy = y - Cen[1];
      const double dz = z - Cen[2];
      const double r  = SQRT( dx*dx + dy*dy + dz*dz );

      if ( HaloUseTable )
      {
         const double *DensTable_r = DensTable + 0*DensTable_Nbin;
         const double *DensTable_d = DensTable + 1*DensTable_Nbin;

         if      (  r < exp( DensTable_r[0] )  )
            dens = exp( DensTable_d[0] );

         else if (  r > exp( DensTable_r[DensTable_Nbin-1] )  )
            dens = 0.0;

         else
            dens = exp( Mis_InterpolateFromTable(DensTable_Nbin, DensTable_r, DensTable_d, log(r)) );
      }
      else
      {
          const double Unit_D_GALIC = 1.0e10*Const_Msun/CUBE(Const_kpc);
          const double r_in_kpc     = r*UNIT_L/Const_kpc;
          dens = Halo_Density( r_in_kpc )*Unit_D_GALIC/UNIT_D;
      }
   }

   fluid[DENS] = (real)dens;
   fluid[REAL] = sqrt( fluid[DENS] );
   fluid[IMAG] = 0.0;

} // FUNCTION : SetGridIC



//-------------------------------------------------------------------------------------------------------
// Function    :  BC
// Description :  Boundary conditions
//
// Note        :
//-------------------------------------------------------------------------------------------------------
void BC( real Array[], const int ArraySize[], real fluid[], const int NVar_Flu,
         const int GhostSize, const int idx[], const double pos[], const double Time,
         const int lv, const int TFluVarIdxList[], double AuxArray[] )
{

// simply call the IC function
   SetGridIC( fluid, pos[0], pos[1], pos[2], Time, lv, AuxArray );

} // FUNCTION : BC



//-------------------------------------------------------------------------------------------------------
// Function    :  End_DiskHeating
// Description :  Free memory before terminating the program
//
// Note        :  1. Linked to the function pointer "End_User_Ptr" to replace "End_User()"
//-------------------------------------------------------------------------------------------------------
void End_DiskHeating()
{

   delete [] DensTable;
   delete [] DispTable;

} // FUNCTION : End_DiskHeating



#ifdef PARTICLE
#ifdef GRAVITY
//-------------------------------------------------------------------------------------------------------
// Function    :  Init_NewDiskRestart()
// Description :  Add a new disk from an existing snapshot
//
// Note        :  Must enable OPT__RESTART_RESET and AddParWhenRestart
//-------------------------------------------------------------------------------------------------------
void Init_NewDiskRestart()
{

   if ( AddParWhenRestart && !OPT__RESTART_RESET )
      Aux_Error( ERROR_INFO, "OPT__RESTART_RESET should be turned on to enable AddParWhenRestart !!\n" );

   if ( amr->Par->Init != PAR_INIT_BY_RESTART  ||  !OPT__RESTART_RESET  ||  !AddParWhenRestart )   return;


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );


   const long NNewPar      = ( MPI_Rank == 0 ) ? AddParWhenRestartNPar : 0;
   const long NPar_AllRank = NNewPar;
   real_par *NewParAttFlt[PAR_NATT_FLT_TOTAL];
   long_par *NewParAttInt[PAR_NATT_INT_TOTAL];

   for (int v=0; v<PAR_NATT_FLT_TOTAL; v++)   NewParAttFlt[v] = new real_par [NNewPar];
   for (int v=0; v<PAR_NATT_INT_TOTAL; v++)   NewParAttInt[v] = new long_par [NNewPar];

// set particle attributes
   real_par *Time_AllRank   =   NewParAttFlt[PAR_TIME];
   real_par *Mass_AllRank   =   NewParAttFlt[PAR_MASS];
   real_par *Pos_AllRank[3] = { NewParAttFlt[PAR_POSX], NewParAttFlt[PAR_POSY], NewParAttFlt[PAR_POSZ] };
   real_par *Vel_AllRank[3] = { NewParAttFlt[PAR_VELX], NewParAttFlt[PAR_VELY], NewParAttFlt[PAR_VELZ] };
   long_par *Type_AllRank   =   NewParAttInt[PAR_TYPE];

   if ( AddParWhenRestartByFile ) // add new disk via DiskHeatingParticleIC
   {
//    load data
      if ( MPI_Rank == 0 ) {
//       NOTE: DiskHeatingParticleIC uses the floating-point type for particle type and assumes single precision
         const char FileName[] = "DiskHeatingParticleIC";
         const int  NParAtt    = 8;    // mass, pos*3, vel*3, type

//       check
         if ( !Aux_CheckFileExist(FileName) )
            Aux_Error( ERROR_INFO, "file \"%s\" does not exist !!\n", FileName );

         FILE *FileTemp = fopen( FileName, "rb" );

         fseek( FileTemp, 0, SEEK_END );

         const long ExpectSize = long(NParAtt)*NPar_AllRank*sizeof(real_par);
         const long FileSize   = ftell( FileTemp );
         if ( FileSize != ExpectSize )
            Aux_Error( ERROR_INFO, "size of the file <%s> = %ld != expect = %ld !!\n",
                       FileName, FileSize, ExpectSize );
         fclose( FileTemp );

         Aux_Message( stdout, "   Loading data ... " );

         real_par *ParData_AllRank = new real_par [ NPar_AllRank*NParAtt ];

//       note that fread() may fail for large files if sizeof(size_t) == 4 instead of 8
         FILE *File = fopen( FileName, "rb" );

         for (int v=0; v<NParAtt; v++)
         {
            fseek( File, v*NPar_AllRank*sizeof(real_par), SEEK_SET );
            fread( ParData_AllRank+v*NPar_AllRank, sizeof(real_par), NPar_AllRank, File );
         }

         fclose( File );
         Aux_Message( stdout, "done\n" );
         Aux_Message( stdout, "   Storing data into particle repository ... " );

         real *ParData1 = new real [NParAtt];

         for (long p=0; p<NPar_AllRank; p++)
         {
//          collect data for the target particle
//          [att][id]
            for (int v=0; v<NParAtt; v++)
               ParData1[v] = ParData_AllRank[ v*NPar_AllRank + p ];

            Time_AllRank[p] = Time[0];
//          mass
            Mass_AllRank[p] = ParData1[0];
//          label
            Type_AllRank[p] = (long_par)ParData1[7]; // 1=CDM halo, 2=disk

//          position
            Pos_AllRank[0][p] = ParData1[1];
            Pos_AllRank[1][p] = ParData1[2];
            Pos_AllRank[2][p] = ParData1[3];
//          velocity
            Vel_AllRank[0][p] = ParData1[4];
            Vel_AllRank[1][p] = ParData1[5];
            Vel_AllRank[2][p] = ParData1[6];
         } // for (long p=0; p<NPar_AllRank; p++)

         Aux_Message( stdout, "done\n" );

//       free memory
         delete [] ParData_AllRank;
         delete [] ParData1;
      } // if ( MPI_Rank == 0 )
   } // if ( AddParWhenRestartByFile )

   else // add a thin disk
   {
#     ifndef SUPPORT_GSL
      Aux_Error( ERROR_INFO, "SUPPORT_GSL must be enabled when AddParWhenRestart=1 and AddParWhenRestartByFile=0 !!\n" );
#     endif

//    read velocity dispersion table
      const bool RowMajor_No  =  false;   // load data into the column-major order
      const bool AllocMem_Yes =  true;    // allocate memory for DispTable
      const int  NCol         =  2;       // total number of columns to load
      const int  Col[NCol]    = {0, 1};   // target columns: (radius, dispersion)

      DispTable_Nbin = Aux_LoadTable( DispTable, DispTableFile, NCol, Col, RowMajor_No, AllocMem_Yes );

      double *DispTable_r = DispTable + 0*DispTable_Nbin;
      double *DispTable_d = DispTable + 1*DispTable_Nbin;

//    convert to code unit (input units of radius and velocity dispersion should be fixed to kpc and km/s respectively)
      for ( int b = 0; b < DispTable_Nbin; b++ )
      {
         DispTable_r[b] = DispTable_r[b]*Const_kpc/UNIT_L;
         DispTable_d[b] = DispTable_d[b]*1e5/UNIT_V;
      }

      if ( MPI_Rank == 0 )
      {
         const double ParM = Disk_Mass / NPar_AllRank;
         double Ran, RanR, RanV, R, Rold, f, f_, phi, RanVec[3];

         Aux_Message( stdout, " Particle Mass = %13.7e\n", ParM );

//       initialize the RNG
         RNG = new RandomNumber_t( 1 );
         RNG->SetSeed( 0, NewDisk_RSeed );

         for (long p=0; p<NPar_AllRank; p++)
         {
            Time_AllRank[p] = Time[0];
//          mass
            Mass_AllRank[p] = ParM;
//          label
            Type_AllRank[p] = (long_par)3;      // use 3 to represent thin disk particles

//          position: statisfying surface density Sigma=Disk_Mass/(2*pi*Disk_R**2)*exp(-R/Disk_R)
            Ran  = RNG->GetValue( 0, 0.0, 1.0 );
            R = 1.0;

            do
            {
               f    = (1.0 + R) * exp(-R) + Ran - 1.0;
               f_   = -R * exp(-R);
               Rold = R;
               R    = R - f / f_;
            }
            while( fabs(R - Rold) / R > 1e-7 );

            RanR      = Disk_R*R;
            phi       = 2*M_PI*RNG->GetValue( 0, 0.0, 1.0 );
            RanVec[0] = RanR*cos(phi);
            RanVec[1] = RanR*sin(phi);
            RanVec[2] = 0;
            for (int d=0; d<3; d++)    Pos_AllRank[d][p] = RanVec[d] + Cen[d];

         } // for ( long p = 0; p < NPar_AllRank; p++ )

         delete RNG;
      } // if ( MPI_Rank == 0 )
   } // if ( !AddParWhenRestartByFile )


// add particles here
   Par_AddParticleAfterInit( NNewPar, NewParAttFlt, NewParAttInt );

// free memory
   for (int v=0; v<PAR_NATT_FLT_TOTAL; v++)   delete [] NewParAttFlt[v];
   for (int v=0; v<PAR_NATT_INT_TOTAL; v++)   delete [] NewParAttInt[v];

#  if ( defined PARTICLE  &&  defined LOAD_BALANCE )
   const double Par_Weight = amr->LB->Par_Weight;
#  else
   const double Par_Weight = 0.0;
#  endif
#  ifdef LOAD_BALANCE
   const UseLBFunc_t UseLB = USELB_YES;
#  else
   const UseLBFunc_t UseLB = USELB_NO;
#  endif

   for (int lv=0; lv<MAX_LEVEL; lv++)
   {
      if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Refining level %d ...\n", lv );

      Flag_Real( lv, UseLB );

      Refine( lv, UseLB );

#     ifdef LOAD_BALANCE
      LB_Init_LoadBalance( true, true, Par_Weight, true, false, lv+1 );
#     endif

      if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Refining level %d ... done\n", lv );
   } // for (int lv=OPT__UM_IC_LEVEL; lv<MAX_LEVEL; lv++)


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : Init_NewDiskRestart



//-------------------------------------------------------------------------------------------------------
// Function    :  Init_NewDiskVelocity()
// Description :  Initialize velocities of the new disk after Poisson solver
//
// Note        :  Must enable OPT__RESTART_RESET and AddParWhenRestart and disable AddParWhenRestartByFile
//-------------------------------------------------------------------------------------------------------
void Init_NewDiskVelocity()
{

   if ( amr->Par->Init != PAR_INIT_BY_RESTART  ||  !OPT__RESTART_RESET  ||
        !AddParWhenRestart  ||  AddParWhenRestartByFile )   return;


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );

#  ifndef SUPPORT_GSL
   Aux_Error( ERROR_INFO, "SUPPORT_GSL must be enabled when AddParWhenRestart=1 and AddParWhenRestartByFile=0 !!\n" );
#  endif

   const long NPar_ThisRank = amr->Par->NPar_AcPlusInac;
   long_par *ParType   =   amr->Par->Type;
   real_par *ParPos[3] = { amr->Par->PosX, amr->Par->PosY, amr->Par->PosZ };
   real_par *ParVel[3] = { amr->Par->VelX, amr->Par->VelY, amr->Par->VelZ };
   real_par *ParAcc[3] = { amr->Par->AccX, amr->Par->AccY, amr->Par->AccZ };

   real ParRadius[2];
   real NormParRadius[2];
   double V_acc, RanV[3], sigma;
   double ParR;

#  ifdef SUPPORT_GSL
// initialize the RNG
   gsl_rng *random_generator;
   random_generator = gsl_rng_alloc(gsl_rng_ranlxd1);
   gsl_rng_set(random_generator, NewDisk_RSeed + MPI_Rank);

   for (long p=0; p<NPar_ThisRank; p++)
   {
      if ( ParType[p] == 3 ) // applies for thin disk particles
      {
         ParRadius[0]     = ParPos[0][p] - Cen[0];
         ParRadius[1]     = ParPos[1][p] - Cen[1];
         ParR             = sqrt( SQR(ParRadius[0]) + SQR(ParRadius[1]) );
         NormParRadius[0] = ParRadius[0]/ ParR;
         NormParRadius[1] = ParRadius[1]/ ParR;

//       compute radial acceleration
         V_acc = sqrt(  fabs( ParRadius[0]*ParAcc[0][p] + ParRadius[1]*ParAcc[1][p] )  );

//       add velocity dispersion
         sigma   = Get_Dispersion( ParR );
         RanV[0] = gsl_ran_gaussian( random_generator, sigma );
         RanV[1] = gsl_ran_gaussian( random_generator, sigma );
         RanV[2] = gsl_ran_gaussian( random_generator, sigma );

//       compute particle velocities using acceleration fields
         ParVel[0][p] = - V_acc*NormParRadius[1]+RanV[0];
         ParVel[1][p] =   V_acc*NormParRadius[0]+RanV[1];
         ParVel[2][p] =   RanV[2];
      }
   } // for (long p=0; p<NPar_ThisRank; p++)

   gsl_rng_free( random_generator );
#  endif // #ifdef SUPPORT_GSL


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : Init_NewDiskVelocity



//-------------------------------------------------------------------------------------------------------
// Function    :  Get_Dispersion
// Description :  Get thin disk velocity dispersion from table
//
//-------------------------------------------------------------------------------------------------------
double Get_Dispersion( double r )
{

   double disp = 0.0;
   const double *DispTable_r = DispTable + 0*DispTable_Nbin;
   const double *DispTable_d = DispTable + 1*DispTable_Nbin;

   if      ( r < DispTable_r[0] )
      disp = DispTable_d[0];

   else if ( r > DispTable_r[DispTable_Nbin-1] )
      disp = DispTable_d[DispTable_Nbin-1];

   else
      disp = Mis_InterpolateFromTable( DispTable_Nbin, DispTable_r, DispTable_d, r );

   return disp;

} // FUNCTION : Get_Dispersion
#endif // #ifdef GRAVITY
#endif // #ifdef PARTICLE
#endif // #if ( MODEL == ELBDM )



//-------------------------------------------------------------------------------------------------------
// Function    :  Init_TestProb_ELBDM_DiskHeating
// Description :  Test problem initializer
//
// Note        :  None
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Init_TestProb_ELBDM_DiskHeating()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );


// validate the compilation flags and runtime parameters
   Validate();


#  if ( MODEL == ELBDM )
// set the problem-specific runtime parameters
   SetParameter();


// set the function pointers
   Init_Function_User_Ptr     = SetGridIC;
   BC_User_Ptr                = BC;
   End_User_Ptr               = End_DiskHeating;
#  ifdef GRAVITY
   Init_ExtPot_Ptr            = Init_ExtPot_Soliton;
#  endif
#  ifdef PARTICLE
   Par_Init_ByFunction_Ptr    = Par_Init_ByFunction_DiskHeating;
   Init_User_Ptr              = Init_NewDiskRestart;
   Init_User_AfterPoisson_Ptr = Init_NewDiskVelocity;
#  endif
#  endif // #if ( MODEL == ELBDM )


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : Init_TestProb_ELBDM_DiskHeating
