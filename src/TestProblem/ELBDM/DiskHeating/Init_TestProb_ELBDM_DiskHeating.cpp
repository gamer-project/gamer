#include "GAMER.h"
#include "TestProb.h"

#ifdef SUPPORT_GSL
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#endif // #ifdef SUPPORT_GSL


static void SetGridIC( real fluid[], const double x, const double y, const double z, const double Time,
                const int lv, double AuxArray[] );

static void BC( real fluid[], const double x, const double y, const double z, const double Time,
                const int lv, double AuxArray[] );
static void End_DiskHeating();
static void AddNewParticleAttribute();
static double get_value_from_interpolation(double r, double *array_radius, double *array_value);
static double get_value_from_interpolation_log(double r, double *array_radius, double *array_value);
static double get_dispersion(double r);
static void GetCenterOfMass( const double CM_Old[], double CM_New[], const double CM_MaxR );
static void Record_Max();

void Init_ExtPot_Soliton();
// problem-specific global variables
// =======================================================================================
FieldIdx_t ParTypeIdx_DiskHeating = Idx_Undefined;
static RandomNumber_t *RNG = NULL;
static double pi = 3.1415926;
static double G = 6.67408E-8;

static double CM_MaxR;      // maximum radius for determining CM
static double CM_TolErrR;   // maximum allowed errors for determining CM
       double Cen[3];       // center
static bool AddFixedHalo;  // add a fixed halo, must enable OPT__FREEZE_FLUID
static bool HaloUseTable;   // 0 = from analytical profile, 1 = from table
       double m_22;         // ELBDM particle mass, used for soliton profile of the fixed halo
       double CoreRadius;   // soliton radius of the fixed halo (in kpc)
static double Rho_0;        // halo rho_0 (in 1.0e+10 Msun*kpc^-3)
static double Rs;           // halo Rs (in kpc)
static double Alpha;        // dimensionless
static double Beta;         // dimensionless
static double Gamma;        // dimensionless
static char   DensTableFile[1000];   // fixed halo density profile filename
static int    DensNbin;              // Nbins of density table
static double *DensTable_r = NULL;   // radius of density table
static double *DensTable_d = NULL;   // density of density table
static bool   AddParWhenRestart;        // add a new disk to an existing snapshot, must enable OPT__RESTART_RESET
static bool   AddParWhenRestartByfile;  // add a new disk via PAR_IC
static int    AddParWhenRestartNPar;    // particle number of the new disk
static double Disk_Mass;                // total mass of the new disk
static double Disk_R;                   // scale radius of the new disk
static char   DispTableFile[1000];   // velocity dispersion filename
static int    DispNbin;              // Nbins of velocity dispersion table
static double *DispTable_r = NULL;   // radius of velocity dispersion table
static double *DispTable_d = NULL;   // velocity dispersion of velocity dispersion table

// =======================================================================================
#ifdef PARTICLE
void Par_Init_ByFunction_DiskHeating( const long NPar_ThisRank, const long NPar_AllRank,
                                real *ParMass, real *ParPosX, real *ParPosY, real *ParPosZ,
                                real *ParVelX, real *ParVelY, real *ParVelZ, real *ParTime,
                                real *ParType, real *AllAttribute[PAR_NATT_TOTAL]);
static void Init_NewDiskRestart();
static void Par_AfterAcceleration( const long NPar_ThisRank, const long NPar_AllRank,
                                 real *ParMass, real *ParPosX, real *ParPosY, real *ParPosZ,
                                 real *ParVelX, real *ParVelY, real *ParVelZ,
                                 real *ParAccX, real *ParAccY, real *ParAccZ,
                                 real *ParTime, real *ParType, real *AllAttribute[PAR_NATT_TOTAL] );


void (*Par_AfterAcceleration_Ptr)( const long NPar_ThisRank, const long NPar_AllRank,
                                 real *ParMass, real *ParPosX, real *ParPosY, real *ParPosZ,
                                 real *ParVelX, real *ParVelY, real *ParVelZ,
                                 real *ParAccX, real *ParAccY, real *ParAccZ,
                                 real *ParTime, real *ParType, real *AllAttribute[PAR_NATT_TOTAL] ) = Par_AfterAcceleration;

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

#  ifdef PARTICLE
//   if ( OPT__INIT == INIT_BY_FUNCTION  &&  amr->Par->Init != PAR_INIT_BY_FUNCTION )
//      Aux_Error( ERROR_INFO, "please set PAR_INIT = 1 (by FUNCTION) !!\n" );

   if ( PAR_NATT_USER != 1 )
      Aux_Error( ERROR_INFO, "please set PAR_NATT_USER = 1 in the Makefile !!\n" );
#  endif


// warnings
   if ( MPI_Rank == 0 )
   {
#     ifndef DUAL_ENERGY
         Aux_Message( stderr, "WARNING : it's recommended to enable DUAL_ENERGY for this test !!\n" );
#     endif

      if ( FLAG_BUFFER_SIZE < 5 )
         Aux_Message( stderr, "WARNING : it's recommended to set FLAG_BUFFER_SIZE >= 5 for this test !!\n" );
   } // if ( MPI_Rank == 0 )



   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Validating test problem %d ... done\n", TESTPROB_ID );

} // FUNCTION : Validate



// replace HYDRO by the target model (e.g., MHD/ELBDM) and also check other compilation flags if necessary (e.g., GRAVITY/PARTICLE)
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
   ReadPara->Add( "CM_MaxR",                 &CM_MaxR,                -1.0,           Eps_double,       NoMax_double      );
   ReadPara->Add( "CM_TolErrR",              &CM_TolErrR,             -1.0,           NoMin_double,     NoMax_double      );
   ReadPara->Add( "CenX",                    &Cen[0],                  NoDef_double,  NoMin_double,     NoMax_double      );
   ReadPara->Add( "CenY",                    &Cen[1],                  NoDef_double,  NoMin_double,     NoMax_double      );
   ReadPara->Add( "CenZ",                    &Cen[2],                  NoDef_double,  NoMin_double,     NoMax_double      );
   ReadPara->Add( "AddFixedHalo",            &AddFixedHalo,            false,         Useless_bool,     Useless_bool      );
   ReadPara->Add( "HaloUseTable",            &HaloUseTable,            false,         Useless_bool,     Useless_bool      );
   ReadPara->Add( "m_22",                    &m_22,                    0.81,          Eps_double,       NoMax_double      );
   ReadPara->Add( "CoreRadius",              &CoreRadius,              1.0,           Eps_double,       NoMax_double      );
   ReadPara->Add( "Rho_0",                   &Rho_0,                   1.0,           Eps_double,       NoMax_double      );
   ReadPara->Add( "Rs",                      &Rs,                      1.0,           Eps_double,       NoMax_double      );
   ReadPara->Add( "Alpha",                   &Alpha,                   1.0,           Eps_double,       NoMax_double      );
   ReadPara->Add( "Beta",                    &Beta,                    1.0,           Eps_double,       NoMax_double      );
   ReadPara->Add( "Gamma",                   &Gamma,                   1.0,           Eps_double,       NoMax_double      );
   ReadPara->Add( "DensTableFile",           DensTableFile,            NoDef_str,     Useless_str,      Useless_str       );
   ReadPara->Add( "DensNbin",                &DensNbin,                0,             0,                NoMax_int         );
   ReadPara->Add( "AddParWhenRestart",       &AddParWhenRestart,       false,         Useless_bool,     Useless_bool      );
   ReadPara->Add( "AddParWhenRestartByfile", &AddParWhenRestartByfile, true,          Useless_bool,     Useless_bool      );
   ReadPara->Add( "AddParWhenRestartNPar",   &AddParWhenRestartNPar,   0,             NoMin_int,        NoMax_int         );
   ReadPara->Add( "Disk_Mass",               &Disk_Mass,               1.0,           Eps_double,       NoMax_double      );
   ReadPara->Add( "Disk_R",                  &Disk_R,                  1.0,           Eps_double,       NoMax_double      );
   ReadPara->Add( "DispTableFile",           DispTableFile,            NoDef_str,     Useless_str,      Useless_str       );
   ReadPara->Add( "DispNbin",                &DispNbin,                0,             0,                NoMax_int         );


   ReadPara->Read( FileName );

   delete ReadPara;

// (1-2) set the default values

// (1-3) check the runtime parameters


// (2) set the problem-specific derived parameters
   // use density table as background fixed halo profile
   if ( HaloUseTable == 1 ) {
      if ( DensNbin == 0){
         if ( MPI_Rank == 0 )    Aux_Error( ERROR_INFO, "DensNbin should be set!\n" );
      }
      int STR_SIZE = 1024;
      int counter, k;
      FILE *input_table;
      char buff[STR_SIZE];
      double radius, value;
      DensTable_r = new double [DensNbin];
      DensTable_d = new double [DensNbin];

      // read the density table
      if ( MPI_Rank == 0 )
      {
         counter = 0;
         input_table = fopen(DensTableFile,"r");
         if (!input_table)
            Aux_Error( ERROR_INFO, "File %s cannot be opened! Exit!\n", DensTableFile);
         else
            Aux_Message( stdout, "Loading %s for reading density table succeed!\n", DensTableFile);
         while (!feof(input_table))
         {
            fgets(buff, STR_SIZE, input_table);
            if (buff==NULL||buff[0]=='\0')
               break;
            else if (buff[0]!='#')
            {
               sscanf(buff, "%lf %lf", &radius, &value);
               DensTable_r[counter] = radius;
               DensTable_d[counter] = value*CUBE(UNIT_L)/UNIT_M;
               counter ++;
            }
            memset(buff, '\0', STR_SIZE);
         }
         fclose(input_table);
      }
      MPI_Bcast(DensTable_r, DensNbin, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      MPI_Bcast(DensTable_d, DensNbin, MPI_DOUBLE, 0, MPI_COMM_WORLD);
   }
// (3) reset other general-purpose parameters
//     --> a helper macro PRINT_WARNING is defined in TestProb.h
   const long   End_Step_Default = __INT_MAX__;
   const double End_T_Default    = __FLT_MAX__; //2.5

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
      Aux_Message( stdout, "  test problem ID           = %d\n",         TESTPROB_ID             );
      Aux_Message( stdout, "  CM_MaxR                   = %13.7e\n",     CM_MaxR                 );
      Aux_Message( stdout, "  CM_TolErrR                = %13.7e\n",     CM_TolErrR              );
      Aux_Message( stdout, "  CenX                      = %13.7e\n",     Cen[0]                  );
      Aux_Message( stdout, "  CenY                      = %13.7e\n",     Cen[1]                  );
      Aux_Message( stdout, "  Cenz                      = %13.7e\n",     Cen[2]                  );
      Aux_Message( stdout, "  AddFixedHalo              = %d\n",         AddFixedHalo           );
      Aux_Message( stdout, "  HaloUseTable              = %d\n",         HaloUseTable            );
      Aux_Message( stdout, "  m_22                      = %13.7e\n",     m_22                    );
      Aux_Message( stdout, "  CoreRadius                = %13.7e\n",     CoreRadius              );
      Aux_Message( stdout, "  Rho_0                     = %13.7e\n",     Rho_0                   );
      Aux_Message( stdout, "  Rs                        = %13.7e\n",     Rs                      );
      Aux_Message( stdout, "  Alpha                     = %13.7e\n",     Alpha                   );
      Aux_Message( stdout, "  Beta                      = %13.7e\n",     Beta                    );
      Aux_Message( stdout, "  Gamma                     = %13.7e\n",     Gamma                   );
      Aux_Message( stdout, "  DensTableFile             = %s\n",         DensTableFile           );
      Aux_Message( stdout, "  DensNbin                  = %d\n",         DensNbin                );
      Aux_Message( stdout, "  AddParWhenRestart         = %d\n",         AddParWhenRestart       );
      Aux_Message( stdout, "  AddParWhenRestartByfile   = %d\n",         AddParWhenRestartByfile );
      Aux_Message( stdout, "  AddParWhenRestartNPar     = %d\n",         AddParWhenRestartNPar   );
      Aux_Message( stdout, "  Disk_Mass                 = %13.7e\n",     Disk_Mass               );
      Aux_Message( stdout, "  Disk_R                    = %13.7e\n",     Disk_R                  );
      Aux_Message( stdout, "  DispTableFile             = %s\n",         DispTableFile           );
      Aux_Message( stdout, "  DispNbin                  = %d\n",         DispNbin                );
      Aux_Message( stdout, "=============================================================================\n" );
   }


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Setting runtime parameters ... done\n" );

} // FUNCTION : SetParameter

// analytical density profile for fixed background halo
double halo_density(double r) {

   double rho_halo = Rho_0/pow(r/Rs, Alpha)/pow(1 + pow(r/Rs,Beta),(Gamma-Alpha)/Beta);
   if ( fabs(rho_halo) <  __FLT_MIN__) rho_halo = 0;

   double rho_soliton =  0.0019 / m_22 / m_22 * pow(1.0 / CoreRadius/ pow(1 + 0.091 * pow(r / CoreRadius, 2), 2), 4);
   if ( fabs(rho_soliton) <  __FLT_MIN__) rho_soliton = 0;

   double rho_max =  0.0019 / m_22 / m_22 * pow(1.0 / CoreRadius, 4);
   if ((rho_halo + rho_soliton) > rho_max) return rho_max;
   else return(rho_halo + rho_soliton);

}

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

   if ( AddFixedHalo && OPT__FREEZE_FLUID)  // add a fixed halo at the background
   {
      const double dx     = (x - Cen[0]);
      const double dy     = (y - Cen[1]);
      const double dz     = (z - Cen[2]);
      const double r      = SQRT( dx*dx + dy*dy + dz*dz );
      const double Unit_L_GALIC = 3.085678e21; //  1.0 kpc
      const double Unit_M_GALIC = 1.989e43;    //  1.0e10 solar masses
      const double Unit_D_GALIC = Unit_M_GALIC/CUBE(Unit_L_GALIC);
      const double r_in_kpc = r*UNIT_L/Unit_L_GALIC;
      if (HaloUseTable)
      {
         if (r_in_kpc < DensTable_r[0]) dens = DensTable_d[0];
         else if (r_in_kpc > DensTable_r[DensNbin-1]) dens = 0;
         else dens = get_value_from_interpolation_log(r_in_kpc, DensTable_r, DensTable_d);
      }
      else dens = halo_density(r_in_kpc)*Unit_D_GALIC/UNIT_D;
   }

//ELBDM example
   fluid[DENS] = (real)dens;
   fluid[REAL] = sqrt( fluid[DENS] );
   fluid[IMAG] = 0.0;

} // FUNCTION : SetGridIC

void BC( real fluid[], const double x, const double y, const double z, const double Time,
         const int lv, double AuxArray[] )
{

   SetGridIC( fluid, x, y, z, Time, lv, AuxArray );

} // FUNCTION : BC

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
   const real   MinEntr_No        = -1.0;
   const bool   DE_Consistency_No = false;
#  ifdef PARTICLE
   const bool   TimingSendPar_No  = false;
   const bool   PredictParPos_No  = false;
   const bool   JustCountNPar_No  = false;
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
/*
#     ifdef PARTICLE
      Prepare_PatchData_InitParticleDensityArray( lv );

      Par_CollectParticle2OneLevel( lv, _PAR_MASS|_PAR_POSX|_PAR_POSY|_PAR_POSZ|_PAR_TYPE ,PredictParPos_No, NULL_REAL, SibBufPatch, FaSibBufPatch, JustCountNPar_No,
                                    TimingSendPar_No );
#     endif
*/
//    get the total density on grids
      TotalDens = new real [ amr->NPatchComma[lv][1] ][PS1][PS1][PS1];
      PID0List  = new int  [ amr->NPatchComma[lv][1]/8 ];

      for (int PID0=0, t=0; PID0<amr->NPatchComma[lv][1]; PID0+=8, t++)    PID0List[t] = PID0;

      Prepare_PatchData( lv, Time[lv], TotalDens[0][0][0], NULL, 0, amr->NPatchComma[lv][1]/8, PID0List, _DENS, _NONE,
                         OPT__RHO_INT_SCHEME, INT_NONE, UNIT_PATCH, NSIDE_00, IntPhase_No, OPT__BC_FLU, BC_POT_NONE,
                         MinDens_No, MinPres_No, MinTemp_No, MinEntr_No, DE_Consistency_No );

      delete [] PID0List;


//    free memory for collecting particles from other ranks and levels, and free density arrays with ghost zones (rho_ext)
/*
#     ifdef PARTICLE
      Par_CollectParticle2OneLevel_FreeMemory( lv, SibBufPatch, FaSibBufPatch );

      Prepare_PatchData_FreeParticleDensityArray( lv );
#     endif
*/

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
// Function    :  Record_Max
// Description :  Record the maximum density and center coordinates
//
// Note        :  1. It will also record the real and imaginary parts associated with the maximum density
//                2. For the center coordinates, it will record the position of maximum density, minimum potential,
//                   and center-of-mass
//                3. Output filenames are fixed to "Record__MaxDens" and "Record__Center"
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Record_Max()
{

   const char filename_max_dens[] = "Record__MaxDens";
   const char filename_center  [] = "Record__Center";
   const int  CountMPI            = 10;

   double dens, max_dens_loc=-__DBL_MAX__, max_dens_pos_loc[3], real_loc, imag_loc;
   double pote, min_pote_loc=+__DBL_MAX__, min_pote_pos_loc[3];
   double send[CountMPI], (*recv)[CountMPI]=new double [MPI_NRank][CountMPI];


// collect local data
   for (int lv=0; lv<NLEVEL; lv++)
   for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
   {
//    skip non-leaf patches
      if ( amr->patch[0][lv][PID]->son != -1 )  continue;

      for (int k=0; k<PS1; k++)  {  const double z = amr->patch[0][lv][PID]->EdgeL[2] + (k+0.5)*amr->dh[lv];
      for (int j=0; j<PS1; j++)  {  const double y = amr->patch[0][lv][PID]->EdgeL[1] + (j+0.5)*amr->dh[lv];
      for (int i=0; i<PS1; i++)  {  const double x = amr->patch[0][lv][PID]->EdgeL[0] + (i+0.5)*amr->dh[lv];

         dens = amr->patch[ amr->FluSg[lv] ][lv][PID]->fluid[DENS][k][j][i];
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
   }


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
         if ( Aux_CheckFileExist(filename_max_dens) )
            Aux_Message( stderr, "WARNING : file \"%s\" already exists !!\n", filename_max_dens );
         else
         {
            FILE *file_max_dens = fopen( filename_max_dens, "w" );
            fprintf( file_max_dens, "#%19s   %10s   %14s   %14s   %14s\n", "Time", "Step", "Dens", "Real", "Imag" );
            fclose( file_max_dens );
         }

         if ( Aux_CheckFileExist(filename_center) )
            Aux_Message( stderr, "WARNING : file \"%s\" already exists !!\n", filename_center );
         else
         {
            FILE *file_center = fopen( filename_center, "w" );
            fprintf( file_center, "#%19s  %10s  %14s  %14s  %14s  %14s  %14s  %14s  %14s  %14s  %10s  %14s  %14s  %14s\n",
                     "Time", "Step", "Dens", "Dens_x", "Dens_y", "Dens_z", "Pote", "Pote_x", "Pote_y", "Pote_z",
                     "NIter", "CM_x", "CM_y", "CM_z" );
            fclose( file_center );
         }

         FirstTime = false;
      }

      FILE *file_max_dens = fopen( filename_max_dens, "a" );
      fprintf( file_max_dens, "%20.14e   %10ld   %14.7e   %14.7e   %14.7e\n",
               Time[0], Step, recv[max_dens_rank][0], recv[max_dens_rank][1], recv[max_dens_rank][2] );
      fclose( file_max_dens );

      FILE *file_center = fopen( filename_center, "a" );
      fprintf( file_center, "%20.14e  %10ld  %14.7e  %14.7e  %14.7e  %14.7e  %14.7e  %14.7e  %14.7e  %14.7e",
               Time[0], Step, recv[max_dens_rank][0], recv[max_dens_rank][3], recv[max_dens_rank][4], recv[max_dens_rank][5],
                              recv[min_pote_rank][6], recv[min_pote_rank][7], recv[min_pote_rank][8], recv[min_pote_rank][9] );
      fclose( file_center );
   } // if ( MPI_Rank == 0 )


// compute the center of mass until convergence
   const double TolErrR2 = SQR( CM_TolErrR );
   const int    NIterMax = 10;

   double dR2, CM_Old[3], CM_New[3];
   int NIter = 0;

// set an initial guess by the peak density position
   if ( MPI_Rank == 0 )
      for (int d=0; d<3; d++)    CM_Old[d] = recv[max_dens_rank][3+d];

   MPI_Bcast( CM_Old, 3, MPI_DOUBLE, 0, MPI_COMM_WORLD );

   while ( true )
   {
      GetCenterOfMass( CM_Old, CM_New, CM_MaxR );

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
         Aux_Message( stderr, "WARNING : dR (%13.7e) > CM_TolErrR (%13.7e) !!\n", sqrt(dR2), CM_TolErrR );

      FILE *file_center = fopen( filename_center, "a" );
      fprintf( file_center, "  %10d  %14.7e  %14.7e  %14.7e\n", NIter, CM_New[0], CM_New[1], CM_New[2] );
      fclose( file_center );
   }

   delete [] recv;

} // FUNCTION : Record_Max


double get_value_from_interpolation(double r, double *array_radius, double *array_value)
{
  int i = 0;
  while(r > array_radius[i])
  {
    i ++;
  }
  return array_value[i-1]+(array_value[i]-array_value[i-1])*(r-array_radius[i-1])/(array_radius[i]-array_radius[i-1]);
} // FUNCTION : get_value_from_interpolation


double get_value_from_interpolation_log(double r, double *array_radius, double *array_value)
{
  int i = 0;
  while(r > array_radius[i])
  {
    i ++;
  }
  return exp(log(array_value[i-1])+(log(array_value[i])-log(array_value[i-1]))*(log(r)-log(array_radius[i-1]))/(log(array_radius[i])-log(array_radius[i-1])));
} // FUNCTION : get_value_from_interpolation_log

void End_DiskHeating()
{
   delete [] DensTable_r;
   delete [] DensTable_d;
   delete [] DispTable_r;
   delete [] DispTable_d;
} // FUNCTION : End_DiskHeating


# ifdef PARTICLE

//-------------------------------------------------------------------------------------------------------
//// Function    :  AddNewParticleAttribute
//// Description :  Add the problem-specific particle attributes
////
//// Note        :  1. Ref: https://github.com/gamer-project/gamer/wiki/Adding-New-Simulations#v-add-problem-specific-grid-fields-and-particle-attributes
////                2. Invoke AddParticleField() for each of the problem-specific particle attribute:
////                   --> Attribute label sent to AddParticleField() will be used as the output name of the attribute
////                   --> Attribute index returned by AddParticleField() can be used to access the particle attribute data
////                3. Pre-declared attribute indices are put in Field.h
////
//// Parameter   :  None
////
//// Return      :  None
////-------------------------------------------------------------------------------------------------------
void AddNewParticleAttribute()
{

   // "Idx_ParLabel" has been predefined in Field.h
   if ( ParTypeIdx_DiskHeating == Idx_Undefined )  ParTypeIdx_DiskHeating = AddParticleAttribute( "ParNumber" );

} // FUNCTION : AddNewParticleAttribute
#ifdef GRAVITY

//-------------------------------------------------------------------------------------------------------
//// Function    :  Init_NewDiskRestart()
//// Description :  Add a new disk from an existing snapshot
//// Note        :  Must enable OPT__RESTART_RESET and AddParWhenRestart
////-------------------------------------------------------------------------------------------------------
void Init_NewDiskRestart()
{

   if ( amr->Par->Init != PAR_INIT_BY_RESTART  || !OPT__RESTART_RESET || !AddParWhenRestart )   return;

   const long   NNewPar        = ( MPI_Rank == 0 ) ? AddParWhenRestartNPar : 0;
   const long   NPar_AllRank   = NNewPar;
   real *NewParAtt[PAR_NATT_TOTAL];

   for (int v=0; v<PAR_NATT_TOTAL; v++)   NewParAtt[v] = new real [NNewPar];

// set particle attributes
   real *Time_AllRank   = NewParAtt[PAR_TIME];
   real *Mass_AllRank   = NewParAtt[PAR_MASS];
   real *Pos_AllRank[3] = { NewParAtt[PAR_POSX], NewParAtt[PAR_POSY], NewParAtt[PAR_POSZ] };
   real *Vel_AllRank[3] = { NewParAtt[PAR_VELX], NewParAtt[PAR_VELY], NewParAtt[PAR_VELZ] };
   real *Type_AllRank   = NewParAtt[PAR_TYPE];
   real *Label_AllRank  = NewParAtt[ParTypeIdx_DiskHeating];

   if ( AddParWhenRestartByfile ) // add new disk via PAR_IC
   {
//    load data
      if ( MPI_Rank == 0 ){
         const char FileName[]     = "PAR_IC";
         const int  NParAtt        = 8;             // mass, pos*3, vel*3, type
         const int  NParAttPerLoad = 1;

//       check
         if ( !Aux_CheckFileExist(FileName) )
            Aux_Error( ERROR_INFO, "file \"%s\" does not exist !!\n", FileName );

         FILE *FileTemp = fopen( FileName, "rb" );

         fseek( FileTemp, 0, SEEK_END );

         const long ExpectSize = long(NParAtt)*NPar_AllRank*sizeof(real);
         const long FileSize   = ftell( FileTemp );
         if ( FileSize != ExpectSize )
            Aux_Error( ERROR_INFO, "size of the file <%s> = %ld != expect = %ld !!\n",
                       FileName, FileSize, ExpectSize );
         fclose( FileTemp );

         Aux_Message( stdout, "   Loading data ... " );

         real *ParData_AllRank = new real [ NPar_AllRank*NParAtt ];

//       note that fread() may fail for large files if sizeof(size_t) == 4 instead of 8
         FILE *File = fopen( FileName, "rb" );

         for (int v=0; v<NParAtt; v+=NParAttPerLoad)
         {
            fseek( File, v*NPar_AllRank*sizeof(real), SEEK_SET );
            fread( ParData_AllRank+v*NPar_AllRank, sizeof(real), long(NParAttPerLoad)*NPar_AllRank, File );
         }

         fclose( File );
         Aux_Message( stdout, "done\n" );
         Aux_Message( stdout, "   Storing data into particle repository ... " );

         real *ParData1 = new real [NParAtt];

         for ( long p = 0; p < NPar_AllRank; p++)
         {
//          collect data for the target particle
//          [att][id]
            for (int v=0; v<NParAtt; v++)
               ParData1[v] = ParData_AllRank[ v*NPar_AllRank + p ];

            Time_AllRank[p] = Time[0];
//          mass
            Mass_AllRank[p] = ParData1[0];
//          label
            Type_AllRank[p] = ParData1[7];
            Label_AllRank[p] = ParData1[7];
//          position
            Pos_AllRank[0][p] = ParData1[1];
            Pos_AllRank[1][p] = ParData1[2];
            Pos_AllRank[2][p] = ParData1[3];
//          velocity
            Vel_AllRank[0][p] = ParData1[4];
            Vel_AllRank[1][p] = ParData1[5];
            Vel_AllRank[2][p] = ParData1[6];

         } // for ( long p = 0; p < NPar_AllRank; p++)
         Aux_Message( stdout, "done\n" );

      } // if ( MPI_Rank == 0 )
   }// if ( AddParWhenRestartByfile )

   else // add a thin disk
   {

#     ifndef SUPPORT_GSL
      Aux_Error( ERROR_INFO, "SUPPORT_GSL must be enabled when AddParWhenRestart=1 and AddParWhenRestartByfile=0 !!\n" );
#     endif

      // read velocity dispersion table
      if ( DispNbin == 0){
         if ( MPI_Rank == 0 )    Aux_Error( ERROR_INFO, "DispNbin should be set!\n" );
      }
      int STR_SIZE = 1024;
      int counter, k;
      FILE *input_table;
      char buff[STR_SIZE];
      double radius, value;
      DispTable_r = new double [DispNbin];
      DispTable_d = new double [DispNbin];

      // velocity dispersion table for thin disk
      if ( MPI_Rank == 0 )
      {
         counter = 0;
         input_table = fopen(DispTableFile,"r");
         if (!input_table)
         Aux_Error( ERROR_INFO, "File %s cannot be opened! Exit!\n", DispTableFile);
         else
            Aux_Message( stdout, "Loading %s for reading density table succeed!\n", DispTableFile);
         while (!feof(input_table))
         {
            fgets(buff, STR_SIZE, input_table);
            if (buff==NULL||buff[0]=='\0')
               break;
            else if (buff[0]!='#')
            {
               sscanf(buff, "%lf %lf", &radius, &value);
               DispTable_r[counter] = radius;
               DispTable_d[counter] = value*1e5/UNIT_V;
               counter ++;
            }
            memset(buff, '\0', STR_SIZE);
         }
         fclose(input_table);
      }
      MPI_Bcast(DispTable_r, DispNbin, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      MPI_Bcast(DispTable_d, DispNbin, MPI_DOUBLE, 0, MPI_COMM_WORLD);

      if ( MPI_Rank == 0 )
      {
         Aux_Message( stdout, "%s ...\n", __FUNCTION__ );

         const double ParM = Disk_Mass / NPar_AllRank;
         double Ran, RanR, RanV, R, Rold, f, f_, phi, RanVec[3];
         Aux_Message(stdout, " Particle Mass = %13.7e\n", ParM) ;

//       initialize the RNG
         RNG = new RandomNumber_t( 1 );
         RNG->SetSeed( 0, 1002 );

         for ( long p = 0; p < NPar_AllRank; p++)
         {
            Time_AllRank[p] = Time[0];
//          mass
            Mass_AllRank[p] = ParM;
//          label
            Type_AllRank[p] = 3;
            Label_AllRank[p] = 3;

//          position
            Ran  = RNG->GetValue( 0, 0.0, 1.0);
            R = 1.0;
            do
            {
               f = (1 + R) * exp(-R) + Ran - 1;
               f_ = -R * exp(-R);
               Rold = R;
               R = R - f / f_;
            }
            while(fabs(R - Rold) / R > 1e-7);
            RanR = Disk_R*R;
            phi = 2*pi*RNG->GetValue( 0, 0.0, 1.0 );
            RanVec[0] = RanR*cos(phi);
            RanVec[1] = RanR*sin(phi);
            RanVec[2] = 0;
            for (int d = 0; d < 3; d++) Pos_AllRank[d][p] = RanVec[d] + Cen[d];

//          velocity
            double V_acc, sigma, r;
            double RanV[3];
            V_acc = sqrt(Ran*Disk_Mass/RanR);
            Vel_AllRank[0][p] = - V_acc*sin(phi);
            Vel_AllRank[1][p] =   V_acc*cos(phi);
            Vel_AllRank[2][p] =   0;

         } // for ( long p = 0; p < NPar_AllRank; p++)

      } //if ( MPI_Rank == 0 )
   } //if ( !AddParWhenRestartByfile )

// add particles here
   Par_AddParticleAfterInit( NNewPar, NewParAtt );
// free memory
   for (int v=0; v<PAR_NATT_TOTAL; v++)   delete [] NewParAtt[v];

// initialize the k-space Green's function for the isolated BC.
   if ( OPT__BC_POT == BC_POT_ISOLATED )  Init_GreenFuncK();

#  if ( defined PARTICLE  &&  defined LOAD_BALANCE )
   const double Par_Weight    = amr->LB->Par_Weight;
#  else
   const double Par_Weight    = 0.0;
#  endif
#  ifdef LOAD_BALANCE
   const UseLBFunc_t UseLB    = USELB_YES;
#  else
   const UseLBFunc_t UseLB    = USELB_NO;
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

    if ( !AddParWhenRestartByfile )
    {
//    evaluate the initial average density if it is not set yet (may already be set in Init_ByRestart)
      if ( AveDensity_Init <= 0.0 )    Poi_GetAverageDensity();

//    evaluate the gravitational potential
      if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", "Calculating gravitational potential" );

      for (int lv=0; lv<NLEVEL; lv++)
      {
         if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Lv %2d ... ", lv );

         Buf_GetBufferData( lv, amr->FluSg[lv], NULL_INT, NULL_INT, DATA_GENERAL, _DENS, _NONE, Rho_ParaBuf, USELB_YES );

         Gra_AdvanceDt( lv, Time[lv], NULL_REAL, NULL_REAL, NULL_INT, amr->PotSg[lv], true, false, false, false, false );

         if ( lv > 0 )

         Buf_GetBufferData( lv, NULL_INT, NULL_INT, amr->PotSg[lv], POT_FOR_POISSON, _POTE, _NONE, Pot_ParaBuf, USELB_YES );

         if ( MPI_Rank == 0 )    Aux_Message( stdout, "done\n" );
      } // for (int lv=0; lv<NLEVEL; lv++)

      if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", "Calculating gravitational potential" );
//    initialize particle acceleration
      if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", "Calculating particle acceleration" );

      const bool StoreAcc_Yes    = true;
      const bool UseStoredAcc_No = false;

      for (int lv=0; lv<NLEVEL; lv++)
      Par_UpdateParticle( lv, amr->PotSgTime[lv][ amr->PotSg[lv] ], NULL_REAL, PAR_UPSTEP_ACC_ONLY, StoreAcc_Yes, UseStoredAcc_No );

      if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", "Calculating particle acceleration" );

      // compute particle velocities using acceleration fields
      Par_AfterAcceleration_Ptr( amr->Par->NPar_AcPlusInac, amr->Par->NPar_Active_AllRank,
                                 amr->Par->Mass, amr->Par->PosX, amr->Par->PosY, amr->Par->PosZ,
                                 amr->Par->VelX, amr->Par->VelY, amr->Par->VelZ,
                                 amr->Par->AccX, amr->Par->AccY, amr->Par->AccZ,
                                 amr->Par->Time, amr->Par->Type, amr->Par->Attribute );

   }// if ( !AddParWhenRestartByfile )

} // FUNCTION : Init_NewDiskRestart

////-------------------------------------------------------------------------------------------------------
//// Function    :  Par_AfterAcceleration()
//// Description :  Compute particle veolcities using acceleration fields
////-------------------------------------------------------------------------------------------------------
void Par_AfterAcceleration(const long NPar_ThisRank, const long NPar_AllRank, real *ParMass,
                                            real *ParPosX, real *ParPosY, real *ParPosZ,
                                            real *ParVelX, real *ParVelY, real *ParVelZ,
                                            real *ParAccX, real *ParAccY, real *ParAccZ,
                                            real *ParTime, real *ParType,  real *AllAttribute[PAR_NATT_TOTAL])
{
#  ifdef SUPPORT_GSL

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );

   real *ParPos[3] = { ParPosX, ParPosY, ParPosZ };
   real *ParVel[3] = { ParVelX, ParVelY, ParVelZ };
   real *ParAcc[3] = { ParAccX, ParAccY, ParAccZ };

   real ParRadius[2];
   real NormParRadius[2];
   double V_acc, RanV[3], sigma;
   double ParR;

   // initialize the RNG
   gsl_rng *random_generator;
   random_generator = gsl_rng_alloc(gsl_rng_ranlxd1);
   gsl_rng_set(random_generator, 42 + MPI_Rank);

   for (long p=0; p<NPar_ThisRank; p++)
   {
      if ( ParMass[p] < 0.0 )  continue;
      if ( ParType[p] == 3){
         ParRadius[0] = ParPos[0][p] - Cen[0];
         ParRadius[1] = ParPos[1][p] - Cen[1];
         ParR = sqrt( SQR(ParRadius[0]) + SQR(ParRadius[1]) );
         NormParRadius[0] = ParRadius[0]/ ParR;
         NormParRadius[1] = ParRadius[1]/ ParR;

         // compute radial acceleration
         V_acc = sqrt(fabs(ParRadius[0]*ParAcc[0][p]+ParRadius[1]*ParAcc[1][p]));

         // add velocity dispersion
         sigma = get_dispersion(ParR);
         RanV[0] = gsl_ran_gaussian(random_generator, sigma);
         RanV[1] = gsl_ran_gaussian(random_generator, sigma);
         RanV[2] = gsl_ran_gaussian(random_generator, sigma);

         ParVel[0][p] = - V_acc*NormParRadius[1]+RanV[0];
         ParVel[1][p] =   V_acc*NormParRadius[0]+RanV[1];
         ParVel[2][p] =   RanV[2];
      }

   }

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );
#  endif //ifdef SUPPORT_GSL

} // FUNCTION : Par_AfterAcceleration


double get_dispersion(double r){
   double disp = 0.0;
   const double r_in_kpc = r*UNIT_L/3.085678e21;
   if (r_in_kpc < DispTable_r[0]) disp = DispTable_d[0];
   else if (r_in_kpc > DispTable_r[DispNbin-1]) disp = DispTable_d[DispNbin-1];
   else disp = get_value_from_interpolation(r_in_kpc, DispTable_r, DispTable_d);

   return disp;
} // FUNCTION : get_dispersion




#endif //ifdef GRAVITY
#endif //ifdef PARTICLE
#endif // #if ( MODEL == ELBDM )



//-------------------------------------------------------------------------------------------------------
// Function    :  Init_TestProb_Template
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


// replace HYDRO by the target model (e.g., MHD/ELBDM) and also check other compilation flags if necessary (e.g., GRAVITY/PARTICLE)
#  if ( MODEL == ELBDM )
// set the problem-specific runtime parameters
   SetParameter();


// procedure to enable a problem-// 1. define a user-specified function (example functions are given below)
// 1. define a user-specified function (example functions are given below)
// 2. declare its function prototype on the top of this file
// 3. set the corresponding function pointer below to the new problem-specific function
// 4. enable the corresponding runtime option in "Input__Parameter"
//    --> for instance, enable OPT__OUTPUT_USER for Output_User_Ptr
   Init_Function_User_Ptr         = SetGridIC;
#  ifdef MHD
//   Init_Function_BField_User_Ptr  = SetBFieldIC;
#  endif
// comment out Init_ByFile_User_Ptr to use the default
//   Init_ByFile_User_Ptr           = NULL; // option: OPT__INIT=3;             example: Init/Init_ByFile.cpp -> Init_ByFile_Default()
//   Init_Field_User_Ptr            = NULL; // set NCOMP_PASSIVE_USER;          example: TestProblem/Hydro/Plummer/Init_TestProb_Hydro_Plummer.cpp --> AddNewField()
//   Flag_User_Ptr                  = NULL; // option: OPT__FLAG_USER;          example: Refine/Flag_User.cpp
//   Mis_GetTimeStep_User_Ptr       = NULL; // option: OPT__DT_USER;            example: Miscellaneous/Mis_GetTimeStep_User.cpp
   BC_User_Ptr                    = BC; // option: OPT__BC_FLU_*=4;         example: TestProblem/ELBDM/ExtPot/Init_TestProb_ELBDM_ExtPot.cpp --> BC()
#  ifdef MHD
//   BC_BField_User_Ptr             = NULL; // option: OPT__BC_FLU_*=4;
#  endif
//   Flu_ResetByUser_Func_Ptr       = NULL; // option: OPT__RESET_FLUID;        example: Fluid/Flu_ResetByUser.cpp
//   Init_DerivedField_User_Ptr     = NULL; // option: OPT__OUTPUT_USER_FIELD;  example: Fluid/Flu_DerivedField_User.cpp
//   Output_User_Ptr                = NULL; // option: OPT__OUTPUT_USER;        example: TestProblem/Hydro/AcousticWave/Init_TestProb_Hydro_AcousticWave.cpp --> OutputError()
   Aux_Record_User_Ptr            = Record_Max; // option: OPT__RECORD_USER;        example: Auxiliary/Aux_Record_User.cpp
//   Init_User_Ptr                  = NULL; // option: none;                    example: none
   End_User_Ptr                   = End_DiskHeating; // option: none;                    example: TestProblem/Hydro/ClusterMerger_vs_Flash/Init_TestProb_ClusterMerger_vs_Flash.cpp --> End_ClusterMerger()
#  ifdef GRAVITY
//   Init_ExtAcc_Ptr                = NULL; // option: OPT__EXT_ACC;            example: SelfGravity/CPU_Gravity/CPU_ExtAcc_PointMass.cpp
//   End_ExtAcc_Ptr                 = NULL;
   Init_ExtPot_Ptr                = Init_ExtPot_Soliton; // option: OPT__EXT_POT;            example: SelfGravity/CPU_Poisson/CPU_ExtPot_PointMass.cpp
//   End_ExtPot_Ptr                 = NULL;
//   Poi_AddExtraMassForGravity_Ptr = NULL; // option: OPT__GRAVITY_EXTRA_MASS; example: none
//   Poi_UserWorkBeforePoisson_Ptr  = NULL; // option: none;                    example: SelfGravity/Poi_UserWorkBeforePoisson.cpp
#  endif
#  ifdef PARTICLE
   Par_Init_ByFunction_Ptr        = Par_Init_ByFunction_DiskHeating; // option: PAR_INIT=1;              example: Particle/Par_Init_ByFunction.cpp
   Par_Init_Attribute_User_Ptr    = AddNewParticleAttribute; // set PAR_NATT_USER;               example: TestProblem/Hydro/AGORA_IsolatedGalaxy/Init_TestProb_Hydro_AGORA_IsolatedGalaxy.cpp --> AddNewParticleAttribute()
   Init_User_Ptr                  = Init_NewDiskRestart;
#  endif
#  if ( EOS == EOS_USER )
//   EoS_Init_Ptr                   = NULL; // option: EOS in the Makefile;     example: EoS/User_Template/CPU_EoS_User_Template.cpp
//   EoS_End_Ptr                    = NULL;
#  endif
#  endif // #if ( MODEL == HYDRO )
//   Src_Init_User_Ptr              = NULL; // option: SRC_USER;                example: SourceTerms/User_Template/CPU_Src_User_Template.cpp



   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : Init_TestProb_ELBDM_DiskHeating
