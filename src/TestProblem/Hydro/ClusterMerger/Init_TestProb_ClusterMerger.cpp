#include "GAMER.h"
#include <string>



// problem-specific global variables
// =======================================================================================
// (1) read from Input__TestProb
       int      Merger_Coll_NumHalos;     // number of clusters
static char   (*Merger_File_Prof)[1000];  // profile table of clusters
       char   (*Merger_File_Par) [1000];  // particle file of clusters
static bool    *Merger_Coll_IsGas;        // (true/false) --> does cluster have gas
       double (*Merger_Coll_Pos)[3];      // position of clusters
       double (*Merger_Coll_Vel)[3];      // velocity of clusters
       double  *CM_BH_Mass;               // black hole mass of clusters
       double  *Jet_HalfHeight;           // half height of the cylinder-shape jet source of clusters
       double  *Jet_Radius;               // radius of the cylinder-shape jet source of clusters

       bool     AGN_feedback;             // turn on/off (1/0) AGN feedback
       int      Accretion_Mode;           // 1: hot mode; 2: code mode; 3: combine (hot + cold)
       double   eta;                      // mass loading factor in jet feedback
       double   eps_f;                    // the radiative efficiency in jet feedback
       double   eps_m;                    // the fraction of total energy that goes into the thermal energy in jet feedback
       double   R_acc;                    // accretion radius: compute the accretion rate
       double   R_dep;                    // radius to deplete the accreted gas
       int      JetDirection_case;        // methods for choosing the jet direction
                                          //    1: fixed at x-axis;
                                          //    2: import from table (generate JetDirection.txt)
                                          //    3: align with angular momentum

static char     JetDirection_file[1000];  // jet direction file

       bool     Merger_Coll_UseMetals;    // (true/false) --> do the clusters have a metal field
       bool     Merger_Coll_LabelCenter;  // (true/false) --> label the particle closest to the center of each cluster
       bool     AdjustBHPos;              // (true/false) --> Adjust the BH position
       bool     AdjustBHVel;              // (true/false) --> Adjust the BH velocity
       double   AdjustPeriod;             // the time interval of adjustment
       bool     fixBH;                    // fix the BH at the simulation box center and set its velocity to be zero (1 cluster only)
// ---------------------------------------------------------------------------------------
// (2) read from files
static int     *Merger_NBin;              // number of radial bins of clusters
static double **Table_R = NULL;           // radius      of clusters
static double **Table_D = NULL;           // density     of clusters
static double **Table_P = NULL;           // pressure    of clusters
static double **Table_M = NULL;           // metallicity of clusters

       long    *NPar_EachCluster;         // number of particles in each cluster
       long     NPar_AllCluster = 0L;     // number of particles in all  clusters

       int      JetDirection_NBin;        // number of bins of the jet direction table
static double  *JetDirection = NULL;      // jet direction[time/theta_1/phi_1/theta_2/phi_2/theta_3/phi_3]
       double  *CM_Jet_Time_table;        // the time  table of jet direction
       double **CM_Jet_Theta_table;       // the theta table of jet direction for clusters
       double **CM_Jet_Phi_table;         // the phi   table of jet direction for clusters
// ---------------------------------------------------------------------------------------
// (3) variables to record
       double (*CM_ClusterCen)[3];        // the center  of each cluster
       double (*CM_BH_Pos)[3];            // BH position of each cluster
       double (*CM_BH_Vel)[3];            // BH velocity of each cluster
       double  *CM_BH_Mdot_tot;           // the total accretion rate of BHs
       double  *CM_BH_Mdot_hot;           // the hot   accretion rate of BHs
       double  *CM_BH_Mdot_cold;          // the cold  accretion rate of BHs
       double  *CM_Jet_Mdot;              // the feedback injeciton rate of mass
       double  *CM_Jet_Pdot;              // the feedback injeciton rate of momentum
       double  *CM_Jet_Edot;              // the feedback injeciton rate of total energy
       double (*CM_Jet_Vec)[3];           // jet direction
       double (*CM_RAcc_GasVel)[3];       // average gas velocity inside the accretion radius
       double  *CM_RAcc_SoundSpeed;       // average sound speed  inside the accreiton radius
       double  *CM_RAcc_GasDens;          // average gas density  inside the accreiton radius
       double  *CM_RAcc_RelativeVel;      // relative velocity between BH and gas for each cluster inside the accretion radius
       double  *CM_RAcc_ColdGasMass;      // cold gas mass        inside the accretion radius
       double  *CM_RAcc_GasMass;          // total gas mass       inside the accretion radius
       double  *CM_RAcc_ParMass;          // total DM mass        inside the accretion radius
       int     *CM_Cluster_NPar_close;    // total number of particles inside the 10 times the accrection radius of each cluster
       double  *CM_Bondi_SinkMass;        // total mass change             in the feedback region in one global time-step
       double  *CM_Bondi_SinkMomX;        // total x-momentum change       ...
       double  *CM_Bondi_SinkMomY;        // total y-momentum change       ...
       double  *CM_Bondi_SinkMomZ;        // total z-momentum change       ...
       double  *CM_Bondi_SinkMomXAbs;     // total |x-momentum| change     ...
       double  *CM_Bondi_SinkMomYAbs;     // total |y-momentum| change     ...
       double  *CM_Bondi_SinkMomZAbs;     // total |z-momentum| change     ...
       double  *CM_Bondi_SinkE;           // total injected energy         ...
       double  *CM_Bondi_SinkEk;          // total injected kinetic energy ...
       double  *CM_Bondi_SinkEt;          // total injected thermal energy ...
       int     *CM_Bondi_SinkNCell;       // total number of finest cells within the feedback region
// ---------------------------------------------------------------------------------------
// (4) other variables
       int      AdjustCount = 0;          // count the number of adjustments
       int      Merger_Coll_NumBHs;       // number of BHs in the simulation

       double  *E_inj_exp;                // the expected amount of injected energy
       double  *M_inj_exp;                // the expected amount of injected gas mass
       double (*ang_mom_sum)[3];

       double  *Jet_WaveK;                // jet wavenumber used in the sin() function to have smooth bidirectional jets
       double  *V_cyl;                    // the volume of jet source
       double  *M_inj;                    // the injected density
       double  *P_inj;                    // the injected momentum
       double  *E_inj;                    // the injected energy
       double  *normalize_const;          // the exact normalization constant

#ifdef MASSIVE_PARTICLES
       long_par *CM_ClusterIdx_Cur;       // the current cluster index
#endif

static FieldIdx_t *ColorFieldsIdx;
#ifdef MASSIVE_PARTICLES
       FieldIdx_t  Idx_ParHalo = Idx_Undefined;
#endif
// =======================================================================================


// problem-specific function prototypes
// =======================================================================================
// (1) external functions
#ifdef MASSIVE_PARTICLES
void AddNewParticleAttribute_ClusterMerger();
void Par_Init_ByFunction_ClusterMerger( const long NPar_ThisRank, const long NPar_AllRank,
                                        real_par *ParMass, real_par *ParPosX, real_par *ParPosY, real_par *ParPosZ,
                                        real_par *ParVelX, real_par *ParVelY, real_par *ParVelZ, real_par *ParTime,
                                        long_par *ParType, real_par *AllAttributeFlt[PAR_NATT_FLT_TOTAL],
                                        long_par *AllAttributeInt[PAR_NATT_INT_TOTAL] );
#endif

void Aux_Record_ClusterMerger();
bool Flag_ClusterMerger( const int i, const int j, const int k, const int lv, const int PID, const double *Threshold );
void AddNewField_ClusterMerger();
void Init_User_ClusterMerger();
int  Flu_ResetByUser_Func_ClusterMerger( real fluid[], const double Emag, const double x, const double y, const double z,
                                         const double Time, const double dt, const int lv, double AuxArray[] );
void Flu_ResetByUser_API_ClusterMerger( const int lv, const int FluSg, const int MagSg, const double TimeNew, const double dt );

extern void (*Flu_ResetByUser_API_Ptr)( const int lv, const int FluSg, const int MagSg, const double TimeNew, const double dt );
// ---------------------------------------------------------------------------------------
// (2) internal functions
static void   AllocateBHVarArray();

#ifdef SUPPORT_HDF5
static int    Read_Num_Points_ClusterMerger( std::string filename );
static void   Read_Profile_ClusterMerger( std::string filename, std::string fieldname, double field[] );
#ifdef MASSIVE_PARTICLES
static long   Read_Particle_Number_ClusterMerger( std::string filename );
#endif

static herr_t LoadField( const char *FieldName, void *FieldPtr, const hid_t H5_SetID_Target,
                         const hid_t H5_TypeID_Target );
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
void Validate()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Validating test problem %d ...\n", TESTPROB_ID );

#  if ( MODEL != HYDRO )
   Aux_Error( ERROR_INFO, "MODEL != HYDRO !!\n" );
#  endif

#  ifndef GRAVITY
   Aux_Error( ERROR_INFO, "GRAVITY must be enabled !!\n" );
#  endif

#  ifndef MASSIVE_PARTICLES
   Aux_Error( ERROR_INFO, "MASSIVE_PARTICLES must be enabled !!\n" );
#  endif

#  ifndef SUPPORT_HDF5
   Aux_Error( ERROR_INFO, "SUPPORT_HDF5 must be enabled !!\n" );
#  endif

#  ifdef COMOVING
   Aux_Error( ERROR_INFO, "COMOVING must be disabled !!\n" );
#  endif

#  ifdef PARTICLE
   if ( PAR_NATT_INT_USER != 1 )
      Aux_Error( ERROR_INFO, "PAR_NATT_INT_USER must be set to 1 in the Makefile !!");
#  endif

   if ( !OPT__UNIT )
      Aux_Error( ERROR_INFO, "OPT__UNIT must be enabled !!\n" );

   for (int f=0; f<6; f++)
     if ( OPT__BC_FLU[f] == BC_FLU_PERIODIC )
        Aux_Error( ERROR_INFO, "do not use periodic BC (OPT__BC_FLU* = 1) for this test !!\n" );

#  ifdef GRAVITY
   if ( OPT__BC_POT == BC_POT_PERIODIC )
      Aux_Error( ERROR_INFO, "do not use periodic BC (OPT__BC_POT = 1) for this test !!\n" );
#  endif

#  ifdef MASSIVE_PARTICLES
   if ( OPT__INIT == INIT_BY_FUNCTION  &&  amr->Par->Init != PAR_INIT_BY_FUNCTION )
      Aux_Error( ERROR_INFO, "please set PAR_INIT = 1 (by FUNCTION) !!\n" );
#  endif

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Validating test problem %d ... done\n", TESTPROB_ID );

} // FUNCTION : Validate



#if ( MODEL == HYDRO  &&  defined MASSIVE_PARTICLES )
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
// LOAD_PARA( load_mode, "KEY_IN_THE_FILE",         &VARIABLE,                 DEFAULT,            MIN,           MAX            );
// ********************************************************************************************************************************
   LOAD_PARA( load_mode, "Merger_Coll_NumHalos",    &Merger_Coll_NumHalos,     2,                  1,             3              );
   for (int c=0; c<Merger_Coll_NumHalos; c++)
   {
      char Merger_File_Prof_name [MAX_STRING];
      char Merger_File_Par_name  [MAX_STRING];
      char Merger_Coll_IsGas_name[MAX_STRING];
      char Merger_Coll_PosX_name [MAX_STRING];
      char Merger_Coll_PosY_name [MAX_STRING];
      char Merger_Coll_VelX_name [MAX_STRING];
      char Merger_Coll_VelY_name [MAX_STRING];
      char CM_BH_Mass_name       [MAX_STRING];
      char Jet_HalfHeight_name   [MAX_STRING];
      char Jet_Radius_name       [MAX_STRING];

      sprintf( Merger_File_Prof_name,  "Merger_File_Prof%d",  c+1 );
      sprintf( Merger_File_Par_name,   "Merger_File_Par%d",   c+1 );
      sprintf( Merger_Coll_IsGas_name, "Merger_Coll_IsGas%d", c+1 );
      sprintf( Merger_Coll_PosX_name,  "Merger_Coll_PosX%d",  c+1 );
      sprintf( Merger_Coll_PosY_name,  "Merger_Coll_PosY%d",  c+1 );
      sprintf( Merger_Coll_VelX_name,  "Merger_Coll_VelX%d",  c+1 );
      sprintf( Merger_Coll_VelY_name,  "Merger_Coll_VelY%d",  c+1 );
      sprintf( CM_BH_Mass_name,        "Bondi_MassBH%d",      c+1 );
      sprintf( Jet_HalfHeight_name,    "Jet_HalfHeight%d",    c+1 );
      sprintf( Jet_Radius_name,        "Jet_Radius%d",        c+1 );

      LOAD_PARA( load_mode, Merger_File_Prof_name,   Merger_File_Prof[c],      NoDef_str,          Useless_str,   Useless_str    );
      LOAD_PARA( load_mode, Merger_File_Par_name,    Merger_File_Par[c],       NoDef_str,          Useless_str,   Useless_str    );
      LOAD_PARA( load_mode, Merger_Coll_IsGas_name, &Merger_Coll_IsGas[c],     true,               Useless_bool,  Useless_bool   );
      LOAD_PARA( load_mode, Merger_Coll_PosX_name,  &Merger_Coll_Pos[c][0],   -1.0,                NoMin_double,  NoMax_double   );
      LOAD_PARA( load_mode, Merger_Coll_PosY_name,  &Merger_Coll_Pos[c][1],   -1.0,                NoMin_double,  NoMax_double   );
      LOAD_PARA( load_mode, Merger_Coll_VelX_name,  &Merger_Coll_Vel[c][0],   -1.0,                NoMin_double,  NoMax_double   );
      LOAD_PARA( load_mode, Merger_Coll_VelY_name,  &Merger_Coll_Vel[c][1],   -1.0,                NoMin_double,  NoMax_double   );
      LOAD_PARA( load_mode, CM_BH_Mass_name,        &CM_BH_Mass[c],           -1.0,                Eps_double,    NoMax_double   );
      LOAD_PARA( load_mode, Jet_HalfHeight_name,    &Jet_HalfHeight[c],       -1.0,                Eps_double,    NoMax_double   );
      LOAD_PARA( load_mode, Jet_Radius_name,        &Jet_Radius[c],           -1.0,                Eps_double,    NoMax_double   );
   }
   LOAD_PARA( load_mode, "Merger_Coll_UseMetals",   &Merger_Coll_UseMetals,    true,               Useless_bool,  Useless_bool   );
   LOAD_PARA( load_mode, "Merger_Coll_LabelCenter", &Merger_Coll_LabelCenter,  true,               Useless_bool,  Useless_bool   );
   LOAD_PARA( load_mode, "R_acc",                   &R_acc,                   -1.0,                NoMin_double,  NoMax_double   );
   LOAD_PARA( load_mode, "R_dep",                   &R_dep,                   -1.0,                NoMin_double,  NoMax_double   );
   LOAD_PARA( load_mode, "AGN_feedback",            &AGN_feedback,             false,              Useless_bool,  Useless_bool   );
   LOAD_PARA( load_mode, "Accretion_Mode",          &Accretion_Mode,           1,                  1,             3              );
   LOAD_PARA( load_mode, "eta",                     &eta,                     -1.0,                NoMin_double,  NoMax_double   );
   LOAD_PARA( load_mode, "eps_f",                   &eps_f,                   -1.0,                NoMin_double,  NoMax_double   );
   LOAD_PARA( load_mode, "eps_m",                   &eps_m,                   -1.0,                NoMin_double,  NoMax_double   );
   LOAD_PARA( load_mode, "AdjustBHPos",             &AdjustBHPos,              false,              Useless_bool,  Useless_bool   );
   LOAD_PARA( load_mode, "AdjustBHVel",             &AdjustBHVel,              false,              Useless_bool,  Useless_bool   );
   LOAD_PARA( load_mode, "AdjustPeriod",            &AdjustPeriod,            -1.0,                NoMin_double,  NoMax_double   );
   LOAD_PARA( load_mode, "JetDirection_case",       &JetDirection_case,        1,                  1,             3              );
   LOAD_PARA( load_mode, "JetDirection_file",        JetDirection_file,        "JetDirection.txt", Useless_str,   Useless_str    );
   LOAD_PARA( load_mode, "fixBH",                   &fixBH,                    false,              Useless_bool,  Useless_bool   );

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
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void SetParameter()
{

#  ifdef SUPPORT_HDF5

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Setting runtime parameters ...\n" );

// (1) load the problem-specific runtime parameters
// (1-1) read parameters from Input__TestProb
   const char FileName[] = "Input__TestProb";
   ReadPara_t *ReadPara  = new ReadPara_t;

// (1-1-1) get the number of halo first
   ReadPara->Add( "Merger_Coll_NumHalos", &Merger_Coll_NumHalos, 2, 1, 3 );
   ReadPara->Read( FileName );
   delete ReadPara;

// (1-1-2) allocate runtime parameter memories
   Merger_File_Prof  = new char   [ Merger_Coll_NumHalos ][ 1000 ];
   Merger_File_Par   = new char   [ Merger_Coll_NumHalos ][ 1000 ];
   Merger_Coll_IsGas = new bool   [ Merger_Coll_NumHalos ];
   Merger_Coll_Pos   = new double [ Merger_Coll_NumHalos ][ 3 ];
   Merger_Coll_Vel   = new double [ Merger_Coll_NumHalos ][ 3 ];
   Jet_HalfHeight    = new double [ Merger_Coll_NumHalos ];
   Jet_Radius        = new double [ Merger_Coll_NumHalos ];

   CM_BH_Mass        = new double [ Merger_Coll_NumHalos ];

// (1-1-3) load the rest cluster parameters
   ReadPara = new ReadPara_t;

   LoadInputTestProb( LOAD_READPARA, ReadPara, NULL );

   ReadPara->Read( FileName );

   delete ReadPara;

// validate that we have the correct number of passive scalars
   if ( Merger_Coll_NumHalos + (int)Merger_Coll_UseMetals != NCOMP_PASSIVE_USER )
      Aux_Error( ERROR_INFO, "please set NCOMP_PASSIVE_USER (currently %d) == Merger_Coll_NumHalos + Merger_Coll_UseMetals (currently %d) in the Makefile !!\n",
                 NCOMP_PASSIVE_USER, Merger_Coll_NumHalos + (int)Merger_Coll_UseMetals );

// set the correct parameters when fixing the BH
   if ( Merger_Coll_NumHalos != 1  &&  fixBH )
   {
      fixBH = false;
      Aux_Message( stdout, "WARNING! Reset fixBH to be false for multiple clusters!\n" );
   }

   if ( fixBH )
   {
      AdjustBHPos = false;
      AdjustBHVel = false;
   }

// convert to code units
   R_acc             *= Const_kpc / UNIT_L;
   R_dep             *= Const_kpc / UNIT_L;
   for (int c=0; c<Merger_Coll_NumHalos; c++)
   {
      Merger_Coll_Pos[c][0] *= Const_kpc / UNIT_L;
      Merger_Coll_Pos[c][1] *= Const_kpc / UNIT_L;
      Merger_Coll_Vel[c][0] *= (Const_km/Const_s) / UNIT_V;
      Merger_Coll_Vel[c][1] *= (Const_km/Const_s) / UNIT_V;
      CM_BH_Mass     [c]    *= Const_Msun / UNIT_M;
      Jet_HalfHeight [c]    *= Const_kpc / UNIT_L;
      Jet_Radius     [c]    *= Const_kpc / UNIT_L;
   }
   AdjustPeriod      *= Const_Myr / UNIT_T;

   ColorFieldsIdx = new FieldIdx_t [ Merger_Coll_NumHalos ];
   for (int c=0; c<Merger_Coll_NumHalos; c++)   ColorFieldsIdx[c] = Idx_Undefined;

   CM_ClusterIdx_Cur = new long_par [ Merger_Coll_NumHalos ];
   for (int c=0; c<Merger_Coll_NumHalos; c++)   CM_ClusterIdx_Cur[c] = c;

   if ( OPT__INIT != INIT_BY_RESTART )
   {
      Table_R     = new double* [ Merger_Coll_NumHalos ];
      Table_D     = new double* [ Merger_Coll_NumHalos ];
      Table_P     = new double* [ Merger_Coll_NumHalos ];
      Table_M     = new double* [ Merger_Coll_NumHalos ];
      Merger_NBin = new int     [ Merger_Coll_NumHalos ];

//    (2) load the radial profiles
      for (int c=0; c<Merger_Coll_NumHalos; c++)
      {
         if ( ! Merger_Coll_IsGas[c] )   continue;

         const std::string filename( Merger_File_Prof[c] );

         if ( !Aux_CheckFileExist( Merger_File_Prof[c] ) )
            Aux_Error( ERROR_INFO, "file \"%s\" does not exist !!\n", Merger_File_Prof[c] );

         if ( MPI_Rank == 0 )
         {
            Merger_NBin[c] = Read_Num_Points_ClusterMerger( filename );
            Aux_Message( stdout, "num_points%d = %d\n", c+1, Merger_NBin[c] );
         }

         MPI_Bcast( &Merger_NBin[c], 1, MPI_INT, 0, MPI_COMM_WORLD );

         Table_R[c] = new double [ Merger_NBin[c] ];
         Table_D[c] = new double [ Merger_NBin[c] ];
         Table_P[c] = new double [ Merger_NBin[c] ];
         Table_M[c] = new double [ Merger_NBin[c] ];

         if ( MPI_Rank == 0 )
         {
            Read_Profile_ClusterMerger( filename, "/fields/radius",   Table_R[c] );
            Read_Profile_ClusterMerger( filename, "/fields/density",  Table_D[c] );
            Read_Profile_ClusterMerger( filename, "/fields/pressure", Table_P[c] );
            if ( Merger_Coll_UseMetals )
               Read_Profile_ClusterMerger( filename, "/fields/metallicity", Table_M[c] );
            else
               for (int i=0; i<Merger_NBin[c]; i++)   Table_M[c][i] = 0.0;
//          convert to code units (assuming the input units are cgs)
            for (int b=0; b<Merger_NBin[c]; b++)
            {
               Table_R[c][b] /= UNIT_L;
               Table_D[c][b] /= UNIT_D;
               Table_P[c][b] /= UNIT_P;
            }
         } // if ( MPI_Rank == 0 )

         MPI_Bcast( Table_R[c], Merger_NBin[c], MPI_DOUBLE, 0, MPI_COMM_WORLD );
         MPI_Bcast( Table_D[c], Merger_NBin[c], MPI_DOUBLE, 0, MPI_COMM_WORLD );
         MPI_Bcast( Table_P[c], Merger_NBin[c], MPI_DOUBLE, 0, MPI_COMM_WORLD );
         MPI_Bcast( Table_M[c], Merger_NBin[c], MPI_DOUBLE, 0, MPI_COMM_WORLD );
      } // for (int c=0; c<Merger_Coll_NumHalos; c++)

//    (2-2) initialize the BH/Halo position and velocity
      if ( fixBH )
      {
         Merger_Coll_Pos[0][0] = amr->BoxCenter[0];
         Merger_Coll_Pos[0][1] = amr->BoxCenter[1];
         Merger_Coll_Vel[0][0] = 0.0;
         Merger_Coll_Vel[0][1] = 0.0;
      }
      for (int c=0; c<Merger_Coll_NumHalos; c++)
      {
         Merger_Coll_Pos[c][2] = amr->BoxCenter[2];
         Merger_Coll_Vel[c][2] = 0.0;
      }

//    set the number of black holes to be the same as the number of clusters initially
      Merger_Coll_NumBHs = Merger_Coll_NumHalos;

//    allocate BH related memories
      AllocateBHVarArray();

//    set initial values
      for (int c=0; c<Merger_Coll_NumBHs; c++)
      {
         CM_BH_Mdot_tot [c] = 0.0;
         CM_BH_Mdot_hot [c] = 0.0;
         CM_BH_Mdot_cold[c] = 0.0;

         CM_Cluster_NPar_close[c] = 0;

         for (int d=0; d<3; d++)   CM_ClusterCen[c][d] = Merger_Coll_Pos[c][d];
         for (int d=0; d<3; d++)   CM_BH_Pos    [c][d] = CM_ClusterCen  [c][d];
         for (int d=0; d<3; d++)   CM_BH_Vel    [c][d] = Merger_Coll_Vel[c][d];

         E_inj_exp[c] = 0.0;
         M_inj_exp[c] = 0.0;
         ang_mom_sum[c][0] = 1.0;
         ang_mom_sum[c][1] = 0.0;
         ang_mom_sum[c][2] = 0.0;
      }

//    (3) determine particle number
      NPar_EachCluster = new long [ Merger_Coll_NumHalos ];
      for (int c=0; c<Merger_Coll_NumHalos; c++)
      {
//       check file existence
         if ( !Aux_CheckFileExist( Merger_File_Par[c] ) )
            Aux_Error( ERROR_INFO, "file \"%s\" does not exist !!\n", Merger_File_Par[c] );

         const std::string pfilename( Merger_File_Par[c] );

         if ( MPI_Rank == 0 )
         {
            NPar_EachCluster[c] = Read_Particle_Number_ClusterMerger( pfilename );

            Aux_Message( stdout, "   Number of particles in cluster %d = %ld\n",
                         c+1, NPar_EachCluster[c] );
         } // if ( MPI_Rank == 0 )

         MPI_Bcast( &NPar_EachCluster[c], 1, MPI_LONG, 0, MPI_COMM_WORLD );

         NPar_AllCluster += NPar_EachCluster[c];
      } // for (int c=0; c<Merger_Coll_NumHalos; c++)

//    overwrite the total number of particles
      amr->Par->NPar_Active_AllRank = NPar_AllCluster;
      PRINT_RESET_PARA( amr->Par->NPar_Active_AllRank, FORMAT_LONG, "(PAR_NPAR in Input__Parameter)" );
   } // if ( OPT__INIT != INIT_BY_RESTART )


// (4) load the jet direction table
   if ( AGN_feedback  &&  JetDirection_case == 2 )
   {
//    allocate memories
      CM_Jet_Theta_table = new double* [ Merger_Coll_NumBHs ];
      CM_Jet_Phi_table   = new double* [ Merger_Coll_NumBHs ];

      const bool RowMajor_No  = false;                    // load data into the column-major order
      const bool AllocMem_Yes = true;                     // allocate memory for JetDirection
      const int  NCol         = 1 + 2*Merger_Coll_NumBHs; // total number of columns to load
                                                          // target columns: (time, theta_1, phi_1, ...)
      int Col[NCol];
      for (int c=0; c<NCol; c++)   Col[c] = c;

      if ( !Aux_CheckFileExist( JetDirection_file ) )
         Aux_Error( ERROR_INFO, "file \"%s\" does not exist !!\n", JetDirection_file );

      JetDirection_NBin = Aux_LoadTable( JetDirection, JetDirection_file, NCol, Col, RowMajor_No, AllocMem_Yes );
      CM_Jet_Time_table = JetDirection + 0*JetDirection_NBin;
      for (int c=0; c<Merger_Coll_NumBHs; c++)
      {
         CM_Jet_Theta_table[c] = JetDirection+(1+2*c)*JetDirection_NBin;
         CM_Jet_Phi_table[c]   = JetDirection+(2+2*c)*JetDirection_NBin;
      }

      for (int b=0; b<JetDirection_NBin; b++)   CM_Jet_Time_table[b] *= Const_Myr/UNIT_T;
   }


// (5) reset other general-purpose parameters
//     --> a helper macro PRINT_WARNING is defined in TestProb.h
   const long   End_Step_Default = __INT_MAX__;
   const double End_T_Default    = 10.0*Const_Gyr/UNIT_T;

   if ( END_STEP < 0 )
   {
      END_STEP = End_Step_Default;
      PRINT_RESET_PARA( END_STEP, FORMAT_LONG, "" );
   }

   if ( END_T < 0.0 )
   {
      END_T = End_T_Default;
      PRINT_RESET_PARA( END_T, FORMAT_REAL, "" );
   }

   if ( Merger_Coll_LabelCenter  &&  !OPT__RECORD_USER )
   {
      OPT__RECORD_USER = true;
      PRINT_RESET_PARA( OPT__RECORD_USER, FORMAT_BOOL, "" );
   }


// (6) make a note
   if ( MPI_Rank == 0 )
   {
      Aux_Message( stdout, "=============================================================================\n" );
      Aux_Message( stdout, "  test problem ID           = %d\n",           TESTPROB_ID );
      Aux_Message( stdout, "  number of clusters        = %d\n",           Merger_Coll_NumHalos );
      Aux_Message( stdout, "  turn on AGN feedback      = %s\n",          (AGN_feedback)? "yes":"no" );
      for (int c=0; c<Merger_Coll_NumHalos; c++) {
      Aux_Message( stdout, "  profile file %d            = %s\n",          c+1,  Merger_File_Prof[c] );
      Aux_Message( stdout, "  particle file %d           = %s\n",          c+1,  Merger_File_Par[c] );
      Aux_Message( stdout, "  cluster %d w/ gas          = %s\n",          c+1, (Merger_Coll_IsGas[c])? "yes":"no" );
      Aux_Message( stdout, "  cluster %d x-position      = %g\n",          c+1,  Merger_Coll_Pos[c][0] );
      Aux_Message( stdout, "  cluster %d y-position      = %g\n",          c+1,  Merger_Coll_Pos[c][1] );
      Aux_Message( stdout, "  cluster %d x-velocity      = %g\n",          c+1,  Merger_Coll_Vel[c][0] );
      Aux_Message( stdout, "  cluster %d y-velocity      = %g\n",          c+1,  Merger_Coll_Vel[c][1] );
      if ( AGN_feedback ) {
      Aux_Message( stdout, "  cluster %d BH mass         = %g\n",          c+1,  CM_BH_Mass[c]     );
      Aux_Message( stdout, "  cluster %d jet half-height = %g\n",          c+1,  Jet_HalfHeight[c] );
      Aux_Message( stdout, "  cluster %d jet radius      = %g\n",          c+1,  Jet_Radius[c]     ); }
      } // for (int c=0; c<Merger_Coll_NumHalos; c++)

      Aux_Message( stdout, "  use metals                = %s\n",          (Merger_Coll_UseMetals)? "yes":"no" );
      Aux_Message( stdout, "  label cluster centers     = %s\n",          (Merger_Coll_LabelCenter)? "yes":"no" );
      if ( AGN_feedback ) {
      Aux_Message( stdout, "  BH fixed                  = %s\n",          (fixBH)? "yes":"no" );
      Aux_Message( stdout, "  accretion mode            = %d\n",          Accretion_Mode      );
      Aux_Message( stdout, "  eta                       = %g\n",          eta                 );
      Aux_Message( stdout, "  eps_f                     = %g\n",          eps_f               );
      Aux_Message( stdout, "  eps_m                     = %g\n",          eps_m               );
      Aux_Message( stdout, "  accretion radius          = %g\n",          R_acc               );
      Aux_Message( stdout, "  depletion radius          = %g\n",          R_dep               );
      Aux_Message( stdout, "  jet direction case        = %d\n",          JetDirection_case   );
      if ( JetDirection_case == 2 ) {
      Aux_Message( stdout, "  jet direction file        = %s\n",          JetDirection_file   );
      }
      Aux_Message( stdout, "  adjust BH position        = %s\n",          (AdjustBHPos)? "yes":"no" );
      Aux_Message( stdout, "  adjust BH velocity        = %s\n",          (AdjustBHVel)? "yes":"no" );
      Aux_Message( stdout, "  adjust period             = %g\n",          AdjustPeriod        );
      } // if ( AGN_feedback )

      Aux_Message( stdout, "=============================================================================\n" );

//    check if the accretion region is larger than the jet cylinder
      for (int c=0; c<Merger_Coll_NumHalos; c++)
      {
         if ( R_acc < Jet_HalfHeight[c] )   Aux_Message( stderr, "WARNING : R_acc (%14.8e) is less than Jet_HalfHeight%d (%14.8e) !!\n", R_acc, c+1, Jet_HalfHeight[c] );
         if ( R_acc < Jet_Radius    [c] )   Aux_Message( stderr, "WARNING : R_acc (%14.8e) is less than Jet_Radius%d (%14.8e) !!\n",     R_acc, c+1, Jet_Radius[c] );
      }
   } // if ( MPI_Rank == 0 )

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Setting runtime parameters ... done\n" );

#  endif // #ifdef SUPPORT_HDF5

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

   bool AnyGasHalo = false;
   for (int c=0; c<Merger_Coll_NumHalos; c++)   AnyGasHalo |= Merger_Coll_IsGas[c];
   if ( ! AnyGasHalo )   return;

   const double pos_in[3] = { x, y, z };

   real Dens = 0.0, MomX = 0.0, MomY = 0.0, MomZ = 0.0, Pres = 0.0, Eint = 0.0, Etot = 0.0, Metl = 0.0;

   for (int c=0; c<Merger_Coll_NumHalos; c++)
   {
      real   dens, pres, metl;
      double r, rr;
      double rmax = Table_R[c][Merger_NBin[c]-1];

//    for each cell, we sum up the density and pressure from each halo and then calculate the weighted velocity
      if ( Merger_Coll_IsGas[c] )
      {
         r    = DIST_3D_DBL( pos_in, Merger_Coll_Pos[c] );
         rr   = r < rmax ? r : rmax;
         dens = Mis_InterpolateFromTable( Merger_NBin[c], Table_R[c], Table_D[c], rr );
         pres = Mis_InterpolateFromTable( Merger_NBin[c], Table_R[c], Table_P[c], rr );
         metl = Mis_InterpolateFromTable( Merger_NBin[c], Table_R[c], Table_M[c], rr );
      }
      else
      {
         r    = HUGE_NUMBER;
         dens = 0.0;
         pres = 0.0;
         metl = 0.0;
      }

      if ( dens == NULL_REAL )   dens = Table_D[c][Merger_NBin[c]-1];
      if ( pres == NULL_REAL )   pres = Table_P[c][Merger_NBin[c]-1];
      if ( metl == NULL_REAL )   metl = Table_M[c][Merger_NBin[c]-1];

      Dens += dens;
      Pres += pres;

      if ( r <= rmax )
      {
         MomX += Merger_Coll_Vel[c][0]*dens;
         MomY += Merger_Coll_Vel[c][1]*dens;
      }

      Metl += metl*dens;

      if ( Merger_Coll_IsGas[c] )
      {
         if ( r < rmax )   fluid[ColorFieldsIdx[c]] = dens;
         else              fluid[ColorFieldsIdx[c]] = 0.0;
      }
   } // for (int c=0; c<Merger_Coll_NumHalos; c++)

// compute the total gas energy
   Eint = EoS_DensPres2Eint_CPUPtr( Dens, Pres, NULL, EoS_AuxArray_Flt, EoS_AuxArray_Int, h_EoS_Table ); // assuming EoS requires no passive scalars
   Etot = Hydro_ConEint2Etot( Dens, MomX, MomY, MomZ, Eint, 0.0 ); // do NOT include magnetic energy here

   fluid[DENS] = Dens;
   fluid[MOMX] = MomX;
   fluid[MOMY] = MomY;
   fluid[MOMZ] = MomZ;
   fluid[ENGY] = Etot;
   if ( Merger_Coll_UseMetals )   fluid[Idx_Metal] = Metl;

} // FUNCTION : SetGridIC
#endif // #if ( MODEL == HYDRO  &&  defined MASSIVE_PARTICLES )



#ifdef SUPPORT_HDF5
//-------------------------------------------------------------------------------------------------------
// Function    :  Output_HDF5_User_ClusterMerger
// Description :  Store the problem specific parameter in HDF5 outputs (Data_*) under User group
//
// Note         : 1. This function only works in MPI_RANK == 0
//                2. We supports int, uint, long, ulong, bool, float, double, and string datatype
//                3. There MUST be more than one parameter to be stored
//                4. The pointer of the data MUST still exist outside the function, e.g. global variables
//
// Parameter   :  HDF5_OutUser : the structure storing the parameters
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Output_HDF5_User_ClusterMerger( HDF5_Output_t *HDF5_OutUser )
{

   HDF5_OutUser->Add( "Merger_Coll_NumBHs", &Merger_Coll_NumBHs );
   for (int c=0; c<Merger_Coll_NumBHs; c++)
   {
      for (int d=0; d<3; d++)
      {
         char BH_Pos_name[50], ClusterCen_name[50], BH_Vel_name[50];
         sprintf( BH_Pos_name,     "BH_Pos_%d_%d",     c, d );
         sprintf( ClusterCen_name, "ClusterCen_%d_%d", c, d );
         sprintf( BH_Vel_name,     "BH_Vel_%d_%d",     c, d );
         HDF5_OutUser->Add( BH_Pos_name,     &CM_BH_Pos[c][d]     );
         HDF5_OutUser->Add( ClusterCen_name, &CM_ClusterCen[c][d] );
         HDF5_OutUser->Add( BH_Vel_name,     &CM_BH_Vel[c][d]     );
      }
      char BH_Mass_name[50], BH_Mdot_tot_name[50], BH_Mdot_hot_name[50], BH_Mdot_cold_name[50];
      sprintf( BH_Mass_name,      "BH_Mass_%d",      c );
      sprintf( BH_Mdot_tot_name,  "BH_Mdot_tot_%d",  c );
      sprintf( BH_Mdot_hot_name,  "BH_Mdot_hot_%d",  c );
      sprintf( BH_Mdot_cold_name, "BH_Mdot_cold_%d", c );
      HDF5_OutUser->Add( BH_Mass_name,      &CM_BH_Mass  [c] );
      HDF5_OutUser->Add( BH_Mdot_tot_name,  &CM_BH_Mdot_tot [c] );
      HDF5_OutUser->Add( BH_Mdot_hot_name,  &CM_BH_Mdot_hot [c] );
      HDF5_OutUser->Add( BH_Mdot_cold_name, &CM_BH_Mdot_cold[c] );
   }
   HDF5_OutUser->Add( "AdjustCount", &AdjustCount );

#  ifdef MASSIVE_PARTICLES
   for (int c=0; c<Merger_Coll_NumHalos; c++)
   {
      char CM_ClusterIdx_Cur_name[50];
      sprintf( CM_ClusterIdx_Cur_name, "CM_ClusterIdx_Cur_%d", c );
      HDF5_OutUser->Add( CM_ClusterIdx_Cur_name, &CM_ClusterIdx_Cur[c] );
   }
#  endif

} // FUNCTION : Output_HDF5_User_ClusterMerger
#endif // #ifdef SUPPORT_HDF5



//-------------------------------------------------------------------------------------------------------
// Function    :  End_ClusterMerger
// Description :  Free memory before terminating the program
//
// Note        :  1. Linked to the function pointer "End_User_Ptr" to replace "End_User()"
//
// Parameter   :  None
//-------------------------------------------------------------------------------------------------------
void End_ClusterMerger()
{

#  ifdef SUPPORT_HDF5
   if ( OPT__INIT != INIT_BY_RESTART )
   {
      delete [] NPar_EachCluster;
      delete [] Merger_NBin;
      for (int c=0; c<Merger_Coll_NumHalos; c++)
      {
         if ( ! Merger_Coll_IsGas[c] )   continue;
         delete [] Table_R[c];
         delete [] Table_D[c];
         delete [] Table_P[c];
         delete [] Table_M[c];
      }
      delete [] Table_R;
      delete [] Table_D;
      delete [] Table_P;
      delete [] Table_M;
   }
#  endif

   delete [] Merger_File_Prof;
   delete [] Merger_File_Par;
   delete [] Merger_Coll_IsGas;
   delete [] Merger_Coll_Pos;
   delete [] Merger_Coll_Vel;

   delete [] Jet_HalfHeight;
   delete [] Jet_Radius;

#  ifdef MASSIVE_PARTICLES
   delete [] CM_ClusterIdx_Cur;
#  endif
   delete [] CM_Cluster_NPar_close;
   delete [] CM_ClusterCen;
   delete [] CM_BH_Pos;
   delete [] CM_BH_Vel;
   delete [] CM_BH_Mass;
   delete [] CM_BH_Mdot_tot;
   delete [] CM_BH_Mdot_hot;
   delete [] CM_BH_Mdot_cold;
   delete [] CM_Jet_Mdot;
   delete [] CM_Jet_Pdot;
   delete [] CM_Jet_Edot;
   delete [] CM_Jet_Vec;
   delete [] CM_RAcc_GasVel;
   delete [] CM_RAcc_SoundSpeed;
   delete [] CM_RAcc_GasDens;
   delete [] CM_RAcc_RelativeVel;
   delete [] CM_RAcc_ColdGasMass;
   delete [] CM_RAcc_GasMass;
   delete [] CM_RAcc_ParMass;
   delete [] CM_Bondi_SinkMass;
   delete [] CM_Bondi_SinkMomX;
   delete [] CM_Bondi_SinkMomY;
   delete [] CM_Bondi_SinkMomZ;
   delete [] CM_Bondi_SinkMomXAbs;
   delete [] CM_Bondi_SinkMomYAbs;
   delete [] CM_Bondi_SinkMomZAbs;
   delete [] CM_Bondi_SinkE;
   delete [] CM_Bondi_SinkEk;
   delete [] CM_Bondi_SinkEt;
   delete [] CM_Bondi_SinkNCell;

   delete [] ColorFieldsIdx;

   if ( AGN_feedback  &&  JetDirection_case == 2 )
   {
      delete [] CM_Jet_Theta_table;
      delete [] CM_Jet_Phi_table;
   }

   delete [] Jet_WaveK;
   delete [] V_cyl;
   delete [] M_inj;
   delete [] P_inj;
   delete [] E_inj;
   delete [] normalize_const;

   delete [] E_inj_exp;
   delete [] M_inj_exp;
   delete [] ang_mom_sum;

} // FUNCTION : End_ClusterMerger



#ifdef MHD
void SetBFieldIC( real magnetic[], const double x, const double y, const double z, const double Time,
                  const int lv, double AuxArray[] )
{

// setting the B-field in this problem is always handled by reading from a uniform grid in the file B_IC

   return;

} // FUNCTION : SetBFieldIC
#endif // #ifdef MHD



//-------------------------------------------------------------------------------------------------------
// Function    :  Init_TestProb_Hydro_ClusterMerger
// Description :  Test problem initializer
//
// Note        :  None
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Init_TestProb_Hydro_ClusterMerger()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );


// validate the compilation flags and runtime parameters
   Validate();


#  if ( MODEL == HYDRO  &&  defined MASSIVE_PARTICLES )
// set the problem-specific runtime parameters
   SetParameter();


// set the function pointers of various problem-specific routines
   Init_Function_User_Ptr         = SetGridIC;
   Flag_User_Ptr                  = Flag_ClusterMerger;
   End_User_Ptr                   = End_ClusterMerger;
   Aux_Record_User_Ptr            = Aux_Record_ClusterMerger;
   Par_Init_ByFunction_Ptr        = Par_Init_ByFunction_ClusterMerger;
   Init_Field_User_Ptr            = AddNewField_ClusterMerger;
   Par_Init_Attribute_User_Ptr    = AddNewParticleAttribute_ClusterMerger;
   Flu_ResetByUser_Func_Ptr       = Flu_ResetByUser_Func_ClusterMerger;
   Flu_ResetByUser_API_Ptr        = Flu_ResetByUser_API_ClusterMerger;
   Init_User_Ptr                  = Init_User_ClusterMerger;

#  ifdef MHD
   Init_Function_BField_User_Ptr = SetBFieldIC;
#  endif
#  ifdef SUPPORT_HDF5
   Output_HDF5_UserPara_Ptr       = Output_HDF5_User_ClusterMerger;
   Output_HDF5_InputTest_Ptr      = LoadInputTestProb;
#  endif
#  endif // if ( MODEL == HYDRO  &&  defined MASSIVE_PARTICLES )


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : Init_TestProb_Hydro_ClusterMerger



#ifdef SUPPORT_HDF5
//-------------------------------------------------------------------------------------------------------
// Function    :  Read_Num_Points_ClusterMerger
// Description :  TBF
//
// Note        :  TBF
//
// Parameter   :  TBF
//
// Return      :  TBF
//-------------------------------------------------------------------------------------------------------
int Read_Num_Points_ClusterMerger( std::string filename )
{

   hid_t   file_id, dataset, dataspace;
   herr_t  status;
   hsize_t dims[1], maxdims[1];

   int rank;

   file_id   = H5Fopen( filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT );
   dataset   = H5Dopen( file_id, "/fields/radius", H5P_DEFAULT );
   dataspace = H5Dget_space( dataset );
   rank      = H5Sget_simple_extent_dims( dataspace, dims, maxdims );

   H5Fclose( file_id );

   return (int)dims[0];

} // FUNCTION : Read_Num_Points_ClusterMerger



//-------------------------------------------------------------------------------------------------------
// Function    :  Read_Profile_ClusterMerger
// Description :  TBF
//
// Note        :  TBF
//
// Parameter   :  TBF
//
// Return      :  TBF
//-------------------------------------------------------------------------------------------------------
void Read_Profile_ClusterMerger( std::string filename, std::string fieldname, double field[] )
{

   hid_t  file_id, dataset;
   herr_t status;

   file_id = H5Fopen( filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT );
   dataset = H5Dopen( file_id, fieldname.c_str(), H5P_DEFAULT );
   status  = H5Dread( dataset, H5T_NATIVE_DOUBLE, H5S_ALL,
                      H5S_ALL, H5P_DEFAULT, field );
   H5Dclose( dataset );

   H5Fclose( file_id );

   return;

} // FUNCTION : Read_Profile_ClusterMerger



//-------------------------------------------------------------------------------------------------------
// Function    :  Read_Particle_Number_ClusterMerger
// Description :  TBF
//
// Note        :  TBF
//
// Parameter   :  TBF
//
// Return      :  TBF
//-------------------------------------------------------------------------------------------------------
long Read_Particle_Number_ClusterMerger( std::string filename )
{

   hid_t   file_id, dataset, dataspace;
   herr_t  status;
   hsize_t dims[1], maxdims[1];

   int rank;

   file_id   = H5Fopen( filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT );
   dataset   = H5Dopen( file_id, "particle_mass", H5P_DEFAULT );
   dataspace = H5Dget_space( dataset );
   rank      = H5Sget_simple_extent_dims( dataspace, dims, maxdims );

   H5Sclose( dataspace );
   H5Dclose( dataset );
   H5Fclose( file_id );

   return (long)dims[0];

} // FUNCTION : Read_Particle_Number_ClusterMerger
#endif // #ifdef SUPPORT_HDF5


#if ( MODEL == HYDRO )
//-------------------------------------------------------------------------------------------------------
// Function    :  AddNewField_ClusterMerger
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
void AddNewField_ClusterMerger()
{

   if ( Merger_Coll_UseMetals )
      Idx_Metal = AddField( "Metal", FIXUP_FLUX_YES, FIXUP_REST_YES, FLOOR_YES, NORMALIZE_NO, INTERP_FRAC_NO );

   for (int c=0; c<Merger_Coll_NumHalos; c++)
   {
      if ( ColorFieldsIdx[c] == Idx_Undefined )
      {
         char ColorField_name[MAX_STRING];
         sprintf( ColorField_name, "ColorField%d", c );
         ColorFieldsIdx[c] = AddField( ColorField_name, FIXUP_FLUX_YES, FIXUP_REST_YES, FLOOR_YES, NORMALIZE_NO, INTERP_FRAC_NO );
      }
   }

} // FUNCTION : AddNewField_ClusterMerger
#endif // #if ( MODEL == HYDRO )



#ifdef MASSIVE_PARTICLES
//-------------------------------------------------------------------------------------------------------
// Function    :  AddNewParticleAttribute_ClusterMerger
// Description :  Adding user-defined particle attributes
//
// Note        :  1. Invoked by Par_Init_Attribute() using the function pointer "Par_Init_Attribute_User_Ptr",
//                   which must be set by a test problem initializer
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void AddNewParticleAttribute_ClusterMerger()
{

   if ( Idx_ParHalo == Idx_Undefined )
      Idx_ParHalo = AddParticleAttributeInt( "ParHalo" );

} // FUNCTION : AddNewParticleAttribute_ClusterMerger
#endif // #ifdef MASSIVE_PARTICLES



//-------------------------------------------------------------------------------------------------------
// Function    :  Init_User_ClusterMerger
// Description :  Restart only: initialize the BH variables (position, velocity, mass) from particles labelled with PTYPE_CEN
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
void Init_User_ClusterMerger()
{

   if ( OPT__INIT != INIT_BY_RESTART )   return;

   if ( OPT__RESTART_RESET )
      Aux_Error( ERROR_INFO, "OPT__RESTART_RESET should be disabled !!\n" );

#  ifdef SUPPORT_HDF5
   const char FileName[] = "RESTART";

   hid_t  H5_FileID, H5_SetID_OutputUser, H5_TypeID_OutputUser;
   herr_t H5_Status;

   H5_FileID = H5Fopen( FileName, H5F_ACC_RDONLY, H5P_DEFAULT );
   if ( H5_FileID < 0 )
      Aux_Error( ERROR_INFO, "failed to open the restart HDF5 file \"%s\" !!\n", FileName );

   H5_SetID_OutputUser  = H5Dopen( H5_FileID, "User/UserPara", H5P_DEFAULT );
   if ( H5_SetID_OutputUser < 0 )
      Aux_Error( ERROR_INFO, "failed to open the dataset \"%s\" !!\n", "User/OutputUser" );

   H5_TypeID_OutputUser = H5Dget_type( H5_SetID_OutputUser );
   if ( H5_TypeID_OutputUser < 0 )
      Aux_Error( ERROR_INFO, "failed to open the datatype of \"%s\" !!\n", "User/OutputUser" );

   LoadField( "Merger_Coll_NumBHs", &Merger_Coll_NumBHs, H5_SetID_OutputUser, H5_TypeID_OutputUser );

// allocate BH related memories
   AllocateBHVarArray();

   for (int c=0; c<Merger_Coll_NumBHs; c++)
   {
      for (int d=0; d<3; d++)
      {
         char BH_Pos_name[50], ClusterCen_name[50], BH_Vel_name[50];
         sprintf( BH_Pos_name,     "BH_Pos_%d_%d",     c, d );
         sprintf( ClusterCen_name, "ClusterCen_%d_%d", c, d );
         sprintf( BH_Vel_name,     "BH_Vel_%d_%d",     c, d );
         LoadField( BH_Pos_name,     &CM_BH_Pos[c][d],     H5_SetID_OutputUser, H5_TypeID_OutputUser );
         LoadField( ClusterCen_name, &CM_ClusterCen[c][d], H5_SetID_OutputUser, H5_TypeID_OutputUser );
         LoadField( BH_Vel_name,     &CM_BH_Vel[c][d],     H5_SetID_OutputUser, H5_TypeID_OutputUser );
      }
      char BH_Mass_name[50], BH_Mdot_tot_name[50], BH_Mdot_hot_name[50], BH_Mdot_cold_name[50];
      sprintf( BH_Mass_name,      "BH_Mass_%d",      c );
      sprintf( BH_Mdot_tot_name,  "BH_Mdot_tot_%d",  c );
      sprintf( BH_Mdot_hot_name,  "BH_Mdot_hot_%d",  c );
      sprintf( BH_Mdot_cold_name, "BH_Mdot_cold_%d", c );
      LoadField( BH_Mass_name,      &CM_BH_Mass     [c], H5_SetID_OutputUser, H5_TypeID_OutputUser );
      LoadField( BH_Mdot_tot_name,  &CM_BH_Mdot_tot [c], H5_SetID_OutputUser, H5_TypeID_OutputUser );
      LoadField( BH_Mdot_hot_name,  &CM_BH_Mdot_hot [c], H5_SetID_OutputUser, H5_TypeID_OutputUser );
      LoadField( BH_Mdot_cold_name, &CM_BH_Mdot_cold[c], H5_SetID_OutputUser, H5_TypeID_OutputUser );
   }
   LoadField( "AdjustCount", &AdjustCount, H5_SetID_OutputUser, H5_TypeID_OutputUser );
#  ifdef MASSIVE_PARTICLES
   for (int c=0; c<Merger_Coll_NumHalos; c++)
   {
      char CM_ClusterIdx_Cur_name[50];
      sprintf( CM_ClusterIdx_Cur_name, "CM_ClusterIdx_Cur_%d", c );
      LoadField( CM_ClusterIdx_Cur_name, &CM_ClusterIdx_Cur[c], H5_SetID_OutputUser, H5_TypeID_OutputUser );
   }
#  endif

   H5_Status = H5Tclose( H5_TypeID_OutputUser );
   H5_Status = H5Dclose( H5_SetID_OutputUser );
   H5_Status = H5Fclose( H5_FileID );

#  endif // #ifdef SUPPORT_HDF5

} // FUNCTION : Init_User_ClusterMerger



//-------------------------------------------------------------------------------------------------------
// Function    :  AllocateBHVarArray
// Description :  Allocate the array of BH variables
//
// Note        :  1. Merger_Coll_NumBHs must be assigned before calling this function
//                2. Memories are freed in End_ClusterMerger()
//-------------------------------------------------------------------------------------------------------
void AllocateBHVarArray()
{

   if ( Merger_Coll_NumBHs < 1 )   Aux_Error( ERROR_INFO, "Merger_Coll_NumBHs < 1 !!\n" );

   CM_Cluster_NPar_close = new int    [ Merger_Coll_NumBHs ];
   CM_ClusterCen         = new double [ Merger_Coll_NumBHs ][ 3 ];
   CM_BH_Pos             = new double [ Merger_Coll_NumBHs ][ 3 ];
   CM_BH_Vel             = new double [ Merger_Coll_NumBHs ][ 3 ];
   CM_BH_Mdot_tot        = new double [ Merger_Coll_NumBHs ];
   CM_BH_Mdot_hot        = new double [ Merger_Coll_NumBHs ];
   CM_BH_Mdot_cold       = new double [ Merger_Coll_NumBHs ];
   CM_Jet_Mdot           = new double [ Merger_Coll_NumBHs ];
   CM_Jet_Pdot           = new double [ Merger_Coll_NumBHs ];
   CM_Jet_Edot           = new double [ Merger_Coll_NumBHs ];
   CM_Jet_Vec            = new double [ Merger_Coll_NumBHs ][ 3 ];
   CM_RAcc_GasVel        = new double [ Merger_Coll_NumBHs ][ 3 ];
   CM_RAcc_SoundSpeed    = new double [ Merger_Coll_NumBHs ];
   CM_RAcc_GasDens       = new double [ Merger_Coll_NumBHs ];
   CM_RAcc_RelativeVel   = new double [ Merger_Coll_NumBHs ];
   CM_RAcc_ColdGasMass   = new double [ Merger_Coll_NumBHs ];
   CM_RAcc_GasMass       = new double [ Merger_Coll_NumBHs ];
   CM_RAcc_ParMass       = new double [ Merger_Coll_NumBHs ];
   CM_Bondi_SinkMass     = new double [ Merger_Coll_NumBHs ];
   CM_Bondi_SinkMomX     = new double [ Merger_Coll_NumBHs ];
   CM_Bondi_SinkMomY     = new double [ Merger_Coll_NumBHs ];
   CM_Bondi_SinkMomZ     = new double [ Merger_Coll_NumBHs ];
   CM_Bondi_SinkMomXAbs  = new double [ Merger_Coll_NumBHs ];
   CM_Bondi_SinkMomYAbs  = new double [ Merger_Coll_NumBHs ];
   CM_Bondi_SinkMomZAbs  = new double [ Merger_Coll_NumBHs ];
   CM_Bondi_SinkE        = new double [ Merger_Coll_NumBHs ];
   CM_Bondi_SinkEk       = new double [ Merger_Coll_NumBHs ];
   CM_Bondi_SinkEt       = new double [ Merger_Coll_NumBHs ];
   CM_Bondi_SinkNCell    = new int    [ Merger_Coll_NumBHs ];

   Jet_WaveK             = new double [ Merger_Coll_NumBHs ];
   V_cyl                 = new double [ Merger_Coll_NumBHs ];
   M_inj                 = new double [ Merger_Coll_NumBHs ];
   P_inj                 = new double [ Merger_Coll_NumBHs ];
   E_inj                 = new double [ Merger_Coll_NumBHs ];
   normalize_const       = new double [ Merger_Coll_NumBHs ];

   E_inj_exp             = new double [ Merger_Coll_NumBHs ];
   M_inj_exp             = new double [ Merger_Coll_NumBHs ];
   ang_mom_sum           = new double [ Merger_Coll_NumBHs ][ 3 ];

   for (int c=0; c<Merger_Coll_NumBHs; c++)
   {
      CM_Cluster_NPar_close[c] = 0;
      for (int d=0; d<3; d++)   CM_ClusterCen [c][d] = 0.0;
      for (int d=0; d<3; d++)   CM_BH_Pos     [c][d] = 0.0;
      for (int d=0; d<3; d++)   CM_BH_Vel     [c][d] = 0.0;
      CM_BH_Mdot_tot       [c] = 0.0;
      CM_BH_Mdot_hot       [c] = 0.0;
      CM_BH_Mdot_cold      [c] = 0.0;
      CM_Jet_Mdot          [c] = 0.0;
      CM_Jet_Pdot          [c] = 0.0;
      CM_Jet_Edot          [c] = 0.0;
      for (int d=0; d<3; d++)   CM_Jet_Vec    [c][d] = 0.0;
      for (int d=0; d<3; d++)   CM_RAcc_GasVel[c][d] = 0.0;
      CM_RAcc_SoundSpeed   [c] = 0.0;
      CM_RAcc_GasDens      [c] = 0.0;
      CM_RAcc_RelativeVel  [c] = 0.0;
      CM_RAcc_ColdGasMass  [c] = 0.0;
      CM_RAcc_GasMass      [c] = 0.0;
      CM_RAcc_ParMass      [c] = 0.0;
      CM_Bondi_SinkMass    [c] = 0.0;
      CM_Bondi_SinkMomX    [c] = 0.0;
      CM_Bondi_SinkMomY    [c] = 0.0;
      CM_Bondi_SinkMomZ    [c] = 0.0;
      CM_Bondi_SinkMomXAbs [c] = 0.0;
      CM_Bondi_SinkMomYAbs [c] = 0.0;
      CM_Bondi_SinkMomZAbs [c] = 0.0;
      CM_Bondi_SinkE       [c] = 0.0;
      CM_Bondi_SinkEk      [c] = 0.0;
      CM_Bondi_SinkEt      [c] = 0.0;
      CM_Bondi_SinkNCell   [c] = 0;

      Jet_WaveK            [c] = 0.0;
      V_cyl                [c] = 0.0;
      M_inj                [c] = 0.0;
      P_inj                [c] = 0.0;
      E_inj                [c] = 0.0;
      normalize_const      [c] = 0.0;

      E_inj_exp            [c] = 0.0;
      M_inj_exp            [c] = 0.0;
      for (int d=0; d<3; d++)   ang_mom_sum   [c][d] = 0.0;
   }

} // FUNCITON : AllocateBHVarArray



#ifdef SUPPORT_HDF5
//-------------------------------------------------------------------------------------------------------
// Function    :  LoadField
// Description :  Load a single field from the input compound dataset
//
// Note        :  1. This function works for arbitary datatype (int, float, char, 1D array ...)
//                2. Memory must be allocated for FieldPtr in advance with sufficent size (except for "char *")
//                3. For loading a string, which has (type(FieldPtr) = (char *)), the memory must be freed
//                   manually by calling free()
//
// Parameter   :  FieldName        : Name of the target field
//                FieldPtr         : Pointer to store the retrieved data
//                H5_SetID_Target  : HDF5 dataset  ID of the target compound variable
//                H5_TypeID_Target : HDF5 datatype ID of the target compound variable
//
// Return      :  Success/fail <-> 0/<0
//-------------------------------------------------------------------------------------------------------
herr_t LoadField( const char *FieldName, void *FieldPtr, const hid_t H5_SetID_Target,
                  const hid_t H5_TypeID_Target )
{

   int    H5_FieldIdx;
   size_t H5_FieldSize;
   hid_t  H5_TypeID_Field;    // datatype ID of the target field in the compound variable
   hid_t  H5_TypeID_Load;     // datatype ID for loading the target field
   herr_t H5_Status;


// load
   H5_FieldIdx = H5Tget_member_index( H5_TypeID_Target, FieldName );

   if ( H5_FieldIdx >= 0 )
   {
      H5_TypeID_Field  = H5Tget_member_type( H5_TypeID_Target, H5_FieldIdx );
      H5_FieldSize     = H5Tget_size( H5_TypeID_Field );

      H5_TypeID_Load   = H5Tcreate( H5T_COMPOUND, H5_FieldSize );
      H5_Status        = H5Tinsert( H5_TypeID_Load, FieldName, 0, H5_TypeID_Field );

      H5_Status        = H5Dread( H5_SetID_Target, H5_TypeID_Load, H5S_ALL, H5S_ALL, H5P_DEFAULT, FieldPtr );
      if ( H5_Status < 0 )    Aux_Error( ERROR_INFO, "failed to load the field \"%s\" !!\n", FieldName );

      H5_Status        = H5Tclose( H5_TypeID_Field );
      H5_Status        = H5Tclose( H5_TypeID_Load  );
   } // if ( H5_FieldIdx >= 0 )
   else
   {
      Aux_Error( ERROR_INFO, "target field \"%s\" does not exist in the restart file !!\n", FieldName );
      return -1;
   } // if ( H5_FieldIdx >= 0 ) ... else ...

   return 0;

} // FUNCTION : LoadField
#endif // #ifdef SUPPORT_HDF5
