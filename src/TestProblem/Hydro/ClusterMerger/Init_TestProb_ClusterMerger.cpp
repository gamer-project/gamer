#include "GAMER.h"

#include <string>

#ifdef SUPPORT_HDF5
#include "hdf5.h"
#endif

// problem-specific global variables
// =======================================================================================
       int     Merger_Coll_NumHalos;      // number of clusters
       bool    AGN_feedback;              // turn on/off (1/0) AGN feedback  
static char    Merger_File_Prof1[1000];   // profile table of cluster 1
static char    Merger_File_Prof2[1000];   // profile table of cluster 2
static char    Merger_File_Prof3[1000];   // profile table of cluster 3
       char    Merger_File_Par1 [1000];   // particle file of cluster 1
       char    Merger_File_Par2 [1000];   // particle file of cluster 2
       char    Merger_File_Par3 [1000];   // particle file of cluster 3
       double  Unit_R1;                   // unit of length in cluster 1's IC file
       double  Unit_D1;                   // unit of density in cluster 1's IC file
       double  Unit_P1;                   // unit of pressure in cluster 1's IC file
       double  Unit_R2;                   // unit of length in cluster 2's IC file
       double  Unit_D2;                   // unit of density in cluster 2's IC file
       double  Unit_P2;                   // unit of pressure in cluster 2's IC file
       double  Unit_R3;                   // unit of length in cluster 3's IC file
       double  Unit_D3;                   // unit of density in cluster 3's IC file
       double  Unit_P3;                   // unit of pressure in cluster 3's IC file
       bool    Merger_Coll_IsGas1;        // (true/false) --> does cluster 1 have gas
       bool    Merger_Coll_IsGas2;        // (true/false) --> does cluster 2 have gas
       bool    Merger_Coll_IsGas3;        // (true/false) --> does cluster 3 have gas
       bool    Merger_Coll_UseMetals;     // (true/false) --> do the clusters have a metal field
       bool    Merger_Coll_LabelCenter;   // (true/false) --> label the particle closest to the center of each cluster
       double  Merger_Coll_PosX1;         // x-position of the first cluster
       double  Merger_Coll_PosY1;         // y-position of the first cluster
       double  Merger_Coll_PosX2;         // x-position of the second cluster
       double  Merger_Coll_PosY2;         // y-position of the second cluster
       double  Merger_Coll_PosX3;         // x-position of the third cluster
       double  Merger_Coll_PosY3;         // y-position of the third cluster
       double  Merger_Coll_VelX1;         // x-velocity of the first cluster
       double  Merger_Coll_VelY1;         // y-velocity of the first cluster
       double  Merger_Coll_VelX2;         // x-velocity of the second cluster
       double  Merger_Coll_VelY2;         // y-velocity of the second cluster
       double  Merger_Coll_VelX3;         // x-velocity of the third cluster
       double  Merger_Coll_VelY3;         // y-velocity of the third cluster
       long    NPar_EachCluster[3];       // Number of particles in each cluster
       long    NPar_AllCluster;           // Number of particles in all clusters

static double *Table_R1 = NULL;           // radius of cluster 1
static double *Table_D1 = NULL;           // density of cluster 1
static double *Table_P1 = NULL;           // pressure of cluster 1
static double *Table_M1 = NULL;           // metallicity of cluster 1

static double *Table_R2 = NULL;           // radius of cluster 2
static double *Table_D2 = NULL;           // density of cluster 2
static double *Table_P2 = NULL;           // pressure of cluster 2
static double *Table_M2 = NULL;           // metallicity of cluster 2

static double *Table_R3 = NULL;           // radius of cluster 3
static double *Table_D3 = NULL;           // density of cluster 3
static double *Table_P3 = NULL;           // pressure of cluster 3
static double *Table_M3 = NULL;           // metallicity of cluster 3

static int     Merger_NBin1;              // number of radial bins of cluster 1
static int     Merger_NBin2;              // number of radial bins of cluster 2
static int     Merger_NBin3;              // number of radial bins of cluster 3

static FieldIdx_t ColorField1Idx = Idx_Undefined;
static FieldIdx_t ColorField2Idx = Idx_Undefined;
static FieldIdx_t ColorField3Idx = Idx_Undefined;

       int    Accretion_Mode;       // 1: hot mode; 2: code mode; 3: combine (hot + cold)   
       double eta;                  // mass loading factor in jet feedback
       double eps_f;                // the radiative efficiency in jet feedback
       double eps_m;                // the fraction of total energy that goes into the thermal energy in jet feedback
       double R_acc;                // accretion radius: compute the accretion rate
       double R_dep;                // radius to deplete the accreted gas 

       double CM_Bondi_SinkMass[3];       // total mass change            in the feedback region in one global time-step
       double CM_Bondi_SinkMomX[3];       // total x-momentum change      ...
       double CM_Bondi_SinkMomY[3];       // total y-momentum change      ...
       double CM_Bondi_SinkMomZ[3];       // total z-momentum change      ...
       double CM_Bondi_SinkMomXAbs[3];    // total |x-momentum| change    ...
       double CM_Bondi_SinkMomYAbs[3];    // total |y-momentum| change    ...
       double CM_Bondi_SinkMomZAbs[3];    // total |z-momentum| change    ...
       double CM_Bondi_SinkE[3];          // total injected energy ...
       double CM_Bondi_SinkEk[3];         // total injected kinetic energy ...
       double CM_Bondi_SinkEt[3];         // total injected thermal energy ...
       int    CM_Bondi_SinkNCell[3];      // total number of finest cells within the feedback region

       double Bondi_MassBH1;        // black hole mass of cluster 1
       double Bondi_MassBH2;        // black hole mass of cluster 2
       double Bondi_MassBH3;        // black hole mass of cluster 3
       double Mdot_BH1;             // the accretion rate of BH 1
       double Mdot_BH2;             // the accretion rate of BH 2
       double Mdot_BH3;             // the accretion rate of BH 3
       double Jet_HalfHeight1;      // half height of the cylinder-shape jet source of cluster 1    
       double Jet_HalfHeight2;      // half height of the cylinder-shape jet source of cluster 2       
       double Jet_HalfHeight3;      // half height of the cylinder-shape jet source of cluster 3
       double Jet_Radius1;          // radius of the cylinder-shape jet source of cluster 1
       double Jet_Radius2;          // radius of the cylinder-shape jet source of cluster 2
       double Jet_Radius3;          // radius of the cylinder-shape jet source of cluster 3
       double Mdot[3];              // the feedback injeciton rate of mass
       double Pdot[3];              // the feedback injeciton rate of momentum
       double Edot[3];              // the feedback injeciton rate of total energy
       double Jet_Vec[3][3];        // jet direction   
       double GasVel[3][3];         // average gas velocity inside the accretion radius
       double SoundSpeed[3];        // average sound speed inside the accreiton radius
       double GasDens[3];           // average gas density inside the accreiton radius
       double RelativeVel[3];       // relative velocity between BH and gas for each cluster
       double ColdGasMass[3];       // cold gas mass inside the accretion radius
       double GasMass[3];           // total gas mass inside the accretion radius
       double ParMass[3];           // total DM mass inside the accretion radius
       double ClusterCen[3][3];     // the center of each cluster
       double BH_Pos[3][3];         // BH position of each cluster
       double BH_Vel[3][3];         // BH velocity of each cluster

       int     JetDirection_NBin;     // number of bins of the jet direction table 
static double *JetDirection = NULL;   // jet direction[time/theta_1/phi_1/theta_2/phi_2/theta_3/phi_3]
       double *Time_table;            // the time table of jet direction 
       double *Theta_table[3];        // the theta table of jet direction for 3 clusters
       double *Phi_table[3];          // the phi table of jet direction for 3 clusters

       bool   AdjustBHPos;    // (true/false) --> Adjust the BH position 
       bool   AdjustBHVel;    // (true/false) --> Adjust the BH velocity 
       double AdjustPeriod;   // the time interval of adjustment
       int    AdjustCount = 0;  // count the number of adjustments
       int    JetDirection_case;  // Methods for choosing the jet direction: 1. Fixed at x-axis; 2. Import from table (generate JetDirection.txt); 3. Align with angular momentum
       bool   fixBH;   // fix the BH at the simulation box center and set its velocity to be zero (1 cluster only)
// =======================================================================================


// problem-specific function prototypes
#ifdef MASSIVE_PARTICLES
long Read_Particle_Number_ClusterMerger(std::string filename);
void Par_Init_ByFunction_ClusterMerger(const long NPar_ThisRank,
                                       const long NPar_AllRank,
                                       real_par *ParMass, real_par *ParPosX, real_par *ParPosY, real_par *ParPosZ,
                                       real_par *ParVelX, real_par *ParVelY, real_par *ParVelZ, real_par *ParTime,
                                       real_par *ParType, real_par *AllAttribute[PAR_NATT_TOTAL]);
void Aux_Record_ClusterMerger();
#endif

bool Flag_ClusterMerger( const int i, const int j, const int k, const int lv, const int PID, const double *Threshold );
int Read_Num_Points_ClusterMerger(std::string filename);
void Read_Profile_ClusterMerger(std::string filename, std::string fieldname,
                                double field[]);
void AddNewField_ClusterMerger();

int Flu_ResetByUser_Func_ClusterMerger( real fluid[], const double Emag, const double x, const double y, const double z,
                                        const double Time, const double dt, const int lv, double AuxArray[] );
void Flu_ResetByUser_API_ClusterMerger( const int lv, const int FluSg, const int MagSg, const double TimeNew, const double dt );

extern void (*Flu_ResetByUser_API_Ptr)( const int lv, const int FluSg, const int MagSg, const double TimeNew, const double dt );
void Output_ClusterMerger();
void Init_User_ClusterMerger();

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

#ifdef SUPPORT_HDF5

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Setting runtime parameters ...\n" );

// (1) load the problem-specific runtime parameters
   const char FileName[] = "Input__TestProb";
   ReadPara_t *ReadPara  = new ReadPara_t;

// add parameters in the following format:
// --> note that VARIABLE, DEFAULT, MIN, and MAX must have the same data type
// --> some handy constants (e.g., NoMin_int, Eps_float, ...) are defined in "include/ReadPara.h"
// ********************************************************************************************************************************
// ReadPara->Add( "KEY_IN_THE_FILE",        &VARIABLE,               DEFAULT,          MIN,           MAX            );
// ********************************************************************************************************************************
   ReadPara->Add( "Merger_Coll_NumHalos",   &Merger_Coll_NumHalos,   2,                1,             3              );
   ReadPara->Add( "AGN_feedback",           &AGN_feedback,           false,            Useless_bool,  Useless_bool   );
   ReadPara->Add( "Merger_Coll_IsGas1",     &Merger_Coll_IsGas1,     true,             Useless_bool,  Useless_bool   );
   ReadPara->Add( "Merger_Coll_IsGas2",     &Merger_Coll_IsGas2,     true,             Useless_bool,  Useless_bool   );
   ReadPara->Add( "Merger_Coll_IsGas3",     &Merger_Coll_IsGas3,     true,             Useless_bool,  Useless_bool   );
   ReadPara->Add( "Merger_File_Prof1",       Merger_File_Prof1,      NoDef_str,        Useless_str,   Useless_str    );
   ReadPara->Add( "Merger_File_Par1",        Merger_File_Par1,       NoDef_str,        Useless_str,   Useless_str    );
   ReadPara->Add( "Merger_File_Prof2",       Merger_File_Prof2,      NoDef_str,        Useless_str,   Useless_str    );
   ReadPara->Add( "Merger_File_Par2",        Merger_File_Par2,       NoDef_str,        Useless_str,   Useless_str    );
   ReadPara->Add( "Merger_File_Prof3",       Merger_File_Prof3,      NoDef_str,        Useless_str,   Useless_str    );
   ReadPara->Add( "Merger_File_Par3",        Merger_File_Par3,       NoDef_str,        Useless_str,   Useless_str    );
   ReadPara->Add( "Unit_R1",                &Unit_R1,                1.0,              NoMin_double,  NoMax_double   );
   ReadPara->Add( "Unit_D1",                &Unit_D1,                1.0,              NoMin_double,  NoMax_double   );
   ReadPara->Add( "Unit_P1",                &Unit_P1,                1.0,              NoMin_double,  NoMax_double   );
   ReadPara->Add( "Unit_R2",                &Unit_R2,                1.0,              NoMin_double,  NoMax_double   );
   ReadPara->Add( "Unit_D2",                &Unit_D2,                1.0,              NoMin_double,  NoMax_double   );
   ReadPara->Add( "Unit_P2",                &Unit_P2,                1.0,              NoMin_double,  NoMax_double   );
   ReadPara->Add( "Unit_R3",                &Unit_R3,                1.0,              NoMin_double,  NoMax_double   );
   ReadPara->Add( "Unit_D3",                &Unit_D3,                1.0,              NoMin_double,  NoMax_double   );
   ReadPara->Add( "Unit_P3",                &Unit_P3,                1.0,              NoMin_double,  NoMax_double   );
   ReadPara->Add( "Merger_Coll_PosX1",      &Merger_Coll_PosX1,     -1.0,              NoMin_double,  NoMax_double   );
   ReadPara->Add( "Merger_Coll_PosY1",      &Merger_Coll_PosY1,     -1.0,              NoMin_double,  NoMax_double   );
   ReadPara->Add( "Merger_Coll_PosX2",      &Merger_Coll_PosX2,     -1.0,              NoMin_double,  NoMax_double   );
   ReadPara->Add( "Merger_Coll_PosY2",      &Merger_Coll_PosY2,     -1.0,              NoMin_double,  NoMax_double   );
   ReadPara->Add( "Merger_Coll_PosX3",      &Merger_Coll_PosX3,     -1.0,              NoMin_double,  NoMax_double   );
   ReadPara->Add( "Merger_Coll_PosY3",      &Merger_Coll_PosY3,     -1.0,              NoMin_double,  NoMax_double   );
   ReadPara->Add( "Merger_Coll_VelX1",      &Merger_Coll_VelX1,     -1.0,              NoMin_double,  NoMax_double   );
   ReadPara->Add( "Merger_Coll_VelY1",      &Merger_Coll_VelY1,     -1.0,              NoMin_double,  NoMax_double   );
   ReadPara->Add( "Merger_Coll_VelX2",      &Merger_Coll_VelX2,     -1.0,              NoMin_double,  NoMax_double   );
   ReadPara->Add( "Merger_Coll_VelY2",      &Merger_Coll_VelY2,     -1.0,              NoMin_double,  NoMax_double   );
   ReadPara->Add( "Merger_Coll_VelX3",      &Merger_Coll_VelX3,     -1.0,              NoMin_double,  NoMax_double   );
   ReadPara->Add( "Merger_Coll_VelY3",      &Merger_Coll_VelY3,     -1.0,              NoMin_double,  NoMax_double   );
   ReadPara->Add( "Merger_Coll_UseMetals",  &Merger_Coll_UseMetals,  true,             Useless_bool,  Useless_bool   );
   ReadPara->Add( "Merger_Coll_LabelCenter",&Merger_Coll_LabelCenter,true,             Useless_bool,  Useless_bool   );
   ReadPara->Add( "Bondi_MassBH1",          &Bondi_MassBH1,          -1.0,             Eps_double,    NoMax_double   );
   ReadPara->Add( "Bondi_MassBH2",          &Bondi_MassBH2,          -1.0,             Eps_double,    NoMax_double   );
   ReadPara->Add( "Bondi_MassBH3",          &Bondi_MassBH3,          -1.0,             Eps_double,    NoMax_double   );
   ReadPara->Add( "R_acc",                  &R_acc,                  -1.0,             NoMin_double,  NoMax_double   );
   ReadPara->Add( "R_dep",                  &R_dep,                  -1.0,             NoMin_double,  NoMax_double   );
   ReadPara->Add( "Mdot_BH1",               &Mdot_BH1,               -1.0,             NoMin_double,  NoMax_double   );
   ReadPara->Add( "Mdot_BH2",               &Mdot_BH2,               -1.0,             NoMin_double,  NoMax_double   );
   ReadPara->Add( "Mdot_BH3",               &Mdot_BH3,               -1.0,             NoMin_double,  NoMax_double   );
   ReadPara->Add( "Jet_HalfHeight1",        &Jet_HalfHeight1,        -1.0,             Eps_double,    NoMax_double   );
   ReadPara->Add( "Jet_HalfHeight2",        &Jet_HalfHeight2,        -1.0,             Eps_double,    NoMax_double   );
   ReadPara->Add( "Jet_HalfHeight3",        &Jet_HalfHeight3,        -1.0,             Eps_double,    NoMax_double   );
   ReadPara->Add( "Jet_Radius1",            &Jet_Radius1,            -1.0,             Eps_double,    NoMax_double   );
   ReadPara->Add( "Jet_Radius2",            &Jet_Radius2,            -1.0,             Eps_double,    NoMax_double   );
   ReadPara->Add( "Jet_Radius3",            &Jet_Radius3,            -1.0,             Eps_double,    NoMax_double   );
   ReadPara->Add( "Accretion_Mode",         &Accretion_Mode,          1,               1,             3              );
   ReadPara->Add( "eta",                    &eta,                    -1.0,             NoMin_double,  NoMax_double   );
   ReadPara->Add( "eps_f",                  &eps_f,                  -1.0,             NoMin_double,  NoMax_double   );
   ReadPara->Add( "eps_m",                  &eps_m,                  -1.0,             NoMin_double,  NoMax_double   );
   ReadPara->Add( "AdjustBHPos",           &AdjustBHPos,            false,             Useless_bool,  Useless_bool   );
   ReadPara->Add( "AdjustBHVel",           &AdjustBHVel,            false,             Useless_bool,  Useless_bool   );
   ReadPara->Add( "AdjustPeriod",          &AdjustPeriod,            -1.0,             NoMin_double,  NoMax_double   );
   ReadPara->Add( "JetDirection_case",     &JetDirection_case,       1,                1,             3              );
   ReadPara->Add( "fixBH",                 &fixBH,                  false,             Useless_bool,  Useless_bool   );

   ReadPara->Read( FileName );

   delete ReadPara;

// Validate that we have the correct number of passive scalars
   if ( Merger_Coll_NumHalos + (int)Merger_Coll_UseMetals != NCOMP_PASSIVE_USER )
      Aux_Error( ERROR_INFO,
                 "please set NCOMP_PASSIVE_USER (currently %d) == Merger_Coll_NumHalos + Merger_Coll_UseMetals (currently %d) in the Makefile !!\n", NCOMP_PASSIVE_USER, Merger_Coll_NumHalos + (int)Merger_Coll_UseMetals );

// Set the correct parameters when fixing the BH
   if ( Merger_Coll_NumHalos != 1 && fixBH == true ){
      fixBH = false;
      Aux_Message( stdout, "WARNING! Reset fixBH to be false for multiple clusters!\n" );
   }
   if ( fixBH == true ){
      AdjustBHPos = false;
      AdjustBHVel = false;
   }

// convert to code units
   Merger_Coll_PosX1 *= Const_kpc / UNIT_L;
   Merger_Coll_PosY1 *= Const_kpc / UNIT_L;
   Merger_Coll_PosX2 *= Const_kpc / UNIT_L;
   Merger_Coll_PosY2 *= Const_kpc / UNIT_L;
   Merger_Coll_PosX3 *= Const_kpc / UNIT_L;
   Merger_Coll_PosY3 *= Const_kpc / UNIT_L;
   Merger_Coll_VelX1 *= (Const_km/Const_s) / UNIT_V;
   Merger_Coll_VelY1 *= (Const_km/Const_s) / UNIT_V;
   Merger_Coll_VelX2 *= (Const_km/Const_s) / UNIT_V;
   Merger_Coll_VelY2 *= (Const_km/Const_s) / UNIT_V;
   Merger_Coll_VelX3 *= (Const_km/Const_s) / UNIT_V;
   Merger_Coll_VelY3 *= (Const_km/Const_s) / UNIT_V;
   Bondi_MassBH1     *= Const_Msun/UNIT_M;
   Bondi_MassBH2     *= Const_Msun/UNIT_M;
   Bondi_MassBH3     *= Const_Msun/UNIT_M;
   R_acc             *= Const_kpc / UNIT_L;
   R_dep             *= Const_kpc / UNIT_L;
   Mdot_BH1          *= (Const_Msun/Const_yr) / (UNIT_M / UNIT_T);
   Mdot_BH2          *= (Const_Msun/Const_yr) / (UNIT_M / UNIT_T);
   Mdot_BH3          *= (Const_Msun/Const_yr) / (UNIT_M / UNIT_T);
   Jet_HalfHeight1   *= Const_kpc / UNIT_L;
   Jet_HalfHeight2   *= Const_kpc / UNIT_L;
   Jet_HalfHeight3   *= Const_kpc / UNIT_L;
   Jet_Radius1       *= Const_kpc / UNIT_L;
   Jet_Radius2       *= Const_kpc / UNIT_L;
   Jet_Radius3       *= Const_kpc / UNIT_L;
   AdjustPeriod      *= Const_Myr / UNIT_T;

   if ( OPT__INIT != INIT_BY_RESTART ) {

//    (2) load the radial profiles
      const std::string filename1(Merger_File_Prof1);
      const std::string filename2(Merger_File_Prof2);
      const std::string filename3(Merger_File_Prof3);

      // cluster 1
      if ( Merger_Coll_IsGas1 ) {

         if ( MPI_Rank == 0 ) {
            Merger_NBin1 = Read_Num_Points_ClusterMerger(filename1);
            Aux_Message(stdout, "num_points1 = %d\n", Merger_NBin1);
         }

         MPI_Bcast(&Merger_NBin1, 1, MPI_INT, 0, MPI_COMM_WORLD);

         Table_R1 = new double [Merger_NBin1];
         Table_D1 = new double [Merger_NBin1];
         Table_P1 = new double [Merger_NBin1];
         Table_M1 = new double [Merger_NBin1];

         if ( MPI_Rank == 0 ) {
            Read_Profile_ClusterMerger(filename1, "/fields/radius", Table_R1);
            Read_Profile_ClusterMerger(filename1, "/fields/density", Table_D1);
            Read_Profile_ClusterMerger(filename1, "/fields/pressure", Table_P1);
            if ( Merger_Coll_UseMetals )
               Read_Profile_ClusterMerger(filename1, "/fields/metallicity", Table_M1);
            else
               for ( int i; i < Merger_NBin1; i++ ) Table_M1[i] = 0.0;
            // convert to code units (assuming the input units are cgs)
            for ( int b=0; b<Merger_NBin1; b++ ) {
               Table_R1[b] *= Unit_R1/UNIT_L;
               Table_D1[b] *= Unit_D1/UNIT_D;
               Table_P1[b] *= Unit_P1/UNIT_P;
            }

         }

         MPI_Bcast(Table_R1, Merger_NBin1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
         MPI_Bcast(Table_D1, Merger_NBin1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
         MPI_Bcast(Table_P1, Merger_NBin1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
         MPI_Bcast(Table_M1, Merger_NBin1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

      } // if ( Merger_Coll_IsGas1 )

      // cluster 2
      if ( Merger_Coll_NumHalos > 1 && Merger_Coll_IsGas2) {

         if (MPI_Rank == 0) {
            Merger_NBin2 = Read_Num_Points_ClusterMerger(filename2);
            Aux_Message(stdout, "num_points2 = %d\n", Merger_NBin2);
         }

         MPI_Bcast(&Merger_NBin2, 1, MPI_INT, 0, MPI_COMM_WORLD);

         Table_R2 = new double [Merger_NBin2];
         Table_D2 = new double [Merger_NBin2];
         Table_P2 = new double [Merger_NBin2];
         Table_M2 = new double [Merger_NBin2];

         if ( MPI_Rank == 0 ) {
            Read_Profile_ClusterMerger(filename2, "/fields/radius", Table_R2);
            Read_Profile_ClusterMerger(filename2, "/fields/density", Table_D2);
            Read_Profile_ClusterMerger(filename2, "/fields/pressure", Table_P2);
            if ( Merger_Coll_UseMetals )
               Read_Profile_ClusterMerger(filename2, "/fields/metallicity", Table_M2);
            else
               for ( int i; i < Merger_NBin2; i++ ) Table_M2[i] = 0.0;
            // convert to code units (assuming the input units are cgs)
            for ( int b=0; b<Merger_NBin2; b++ ) {
               Table_R2[b] *= Unit_R2/UNIT_L;
               Table_D2[b] *= Unit_D2/UNIT_D;
               Table_P2[b] *= Unit_P2/UNIT_P;
            }

         }

         MPI_Bcast(Table_R2, Merger_NBin2, MPI_DOUBLE, 0, MPI_COMM_WORLD);
         MPI_Bcast(Table_D2, Merger_NBin2, MPI_DOUBLE, 0, MPI_COMM_WORLD);
         MPI_Bcast(Table_P2, Merger_NBin2, MPI_DOUBLE, 0, MPI_COMM_WORLD);
         MPI_Bcast(Table_M2, Merger_NBin2, MPI_DOUBLE, 0, MPI_COMM_WORLD);

      } // if ( Merger_Coll_NumHalos > 1 && Merger_Coll_IsGas2 )

      // cluster 3
      if ( Merger_Coll_NumHalos > 2 && Merger_Coll_IsGas3) {

         if ( MPI_Rank == 0 ) {
            Merger_NBin3 = Read_Num_Points_ClusterMerger(filename3);
            Aux_Message(stdout, "num_points3 = %d\n", Merger_NBin3);
         }

         MPI_Bcast(&Merger_NBin3, 1, MPI_INT, 0, MPI_COMM_WORLD);

         Table_R3 = new double [Merger_NBin3];
         Table_D3 = new double [Merger_NBin3];
         Table_P3 = new double [Merger_NBin3];
         Table_M3 = new double [Merger_NBin3];

         if ( MPI_Rank == 0 ) {
            Read_Profile_ClusterMerger(filename3, "/fields/radius", Table_R3);
            Read_Profile_ClusterMerger(filename3, "/fields/density", Table_D3);
            Read_Profile_ClusterMerger(filename3, "/fields/pressure", Table_P3);
            if ( Merger_Coll_UseMetals )
               Read_Profile_ClusterMerger(filename3, "/fields/metallicity", Table_M3);
            else
               for ( int i; i < Merger_NBin3; i++ ) Table_M3[i] = 0.0;
            // convert to code units (assuming the input units are cgs)
            for (int b=0; b<Merger_NBin3; b++) {
               Table_R3[b] *= Unit_R3/UNIT_L;
               Table_D3[b] *= Unit_D3/UNIT_D;
               Table_P3[b] *= Unit_P3/UNIT_P;
            }

         }

         MPI_Bcast(Table_R3, Merger_NBin3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
         MPI_Bcast(Table_D3, Merger_NBin3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
         MPI_Bcast(Table_P3, Merger_NBin3, MPI_DOUBLE, 0, MPI_COMM_WORLD);
         MPI_Bcast(Table_M3, Merger_NBin3, MPI_DOUBLE, 0, MPI_COMM_WORLD);

      } // if ( Merger_Coll_NumHalos > 2 && Merger_Coll_IsGas3 )

//    (2-2) Initialize the BH position and velocity
      double ClusterCenter[3][3] = {{ Merger_Coll_PosX1, Merger_Coll_PosY1, amr->BoxCenter[2] },
                                    { Merger_Coll_PosX2, Merger_Coll_PosY2, amr->BoxCenter[2] },
                                    { Merger_Coll_PosX3, Merger_Coll_PosY3, amr->BoxCenter[2] }};
      double CenterVel[3][3]     = {{ Merger_Coll_VelX1, Merger_Coll_VelY1, 0.0 },
                                    { Merger_Coll_VelX2, Merger_Coll_VelY2, 0.0 },
                                    { Merger_Coll_VelX3, Merger_Coll_VelY3, 0.0 }}; 

      for (int c=0; c<Merger_Coll_NumHalos; c++) {
         for (int d=0; d<3; d++)   ClusterCen[c][d] = ClusterCenter[c][d];
         for (int d=0; d<3; d++)   BH_Pos[c][d] = ClusterCen[c][d];
         for (int d=0; d<3; d++)   BH_Vel[c][d] = CenterVel[c][d];
      }

//    (3) Determine particle number

      // check file existence
      if ( !Aux_CheckFileExist(Merger_File_Par1) )
         Aux_Error( ERROR_INFO, "file \"%s\" does not exist !!\n", Merger_File_Par1 );

      if ( Merger_Coll_NumHalos > 1  &&  !Aux_CheckFileExist(Merger_File_Par2) )
         Aux_Error( ERROR_INFO, "file \"%s\" does not exist !!\n", Merger_File_Par2 );

      if ( Merger_Coll_NumHalos > 2  &&  !Aux_CheckFileExist(Merger_File_Par3) )
         Aux_Error( ERROR_INFO, "file \"%s\" does not exist !!\n", Merger_File_Par3 );

      const std::string pfilename1(Merger_File_Par1);
      const std::string pfilename2(Merger_File_Par2);
      const std::string pfilename3(Merger_File_Par3);

      if ( MPI_Rank == 0 ) {

         NPar_EachCluster[0] = Read_Particle_Number_ClusterMerger(pfilename1);

         Aux_Message( stdout, "   Number of particles in cluster 1 = %ld\n",
                      NPar_EachCluster[0] );

         if ( Merger_Coll_NumHalos > 1 ) {
            NPar_EachCluster[1] = Read_Particle_Number_ClusterMerger(pfilename2);
            Aux_Message( stdout, "   Number of particles in cluster 2 = %ld\n", NPar_EachCluster[1] );
         } else {
            NPar_EachCluster[1] = 0;
         }

         if ( Merger_Coll_NumHalos > 2 ) {
            NPar_EachCluster[2] = Read_Particle_Number_ClusterMerger(pfilename3);
            Aux_Message( stdout, "   Number of particles in cluster 3 = %ld\n", NPar_EachCluster[2] );
         } else {
            NPar_EachCluster[2] = 0;
         }

      }

      MPI_Bcast(NPar_EachCluster, 3, MPI_LONG, 0, MPI_COMM_WORLD);

      NPar_AllCluster = NPar_EachCluster[0] + NPar_EachCluster[1] + NPar_EachCluster[2];

      // overwrite the total number of particles
      amr->Par->NPar_Active_AllRank = NPar_AllCluster;
      PRINT_RESET_PARA( amr->Par->NPar_Active_AllRank, FORMAT_LONG, "(PAR_NPAR in Input__Parameter)" );

   } // if ( OPT__INIT != INIT_BY_RESTART )


// (4) Load the jet direction table
   if ( AGN_feedback  &&  JetDirection_case == 2 ) {
      const bool RowMajor_No  = false;    // load data into the column-major order
      const bool AllocMem_Yes = true;     // allocate memory for JetDirection 
      const int  NCol         = 7;        // total number of columns to load
      const int  Col[NCol] = {0, 1, 2, 3, 4, 5, 6};  // target columns: (time, theta_1, phi_1, theta_2, phi_2, theta_3, phi_3)
 
      JetDirection_NBin = Aux_LoadTable( JetDirection, "JetDirection.txt", NCol, Col, RowMajor_No, AllocMem_Yes );
      Time_table = JetDirection+0*JetDirection_NBin;
      for (int d=0; d<3; d++) {
         Theta_table[d] = JetDirection+(1+2*d)*JetDirection_NBin;
         Phi_table[d]   = JetDirection+(2+2*d)*JetDirection_NBin;
      }       
      for (int b=0; b<JetDirection_NBin; b++)   Time_table[b] *= Const_Myr/UNIT_T;
   }       


// (5) reset other general-purpose parameters
//     --> a helper macro PRINT_WARNING is defined in TestProb.h
   const long   End_Step_Default = __INT_MAX__;
   const double End_T_Default    = 10.0*Const_Gyr/UNIT_T;

   if ( END_STEP < 0 ) {
      END_STEP = End_Step_Default;
      PRINT_RESET_PARA( END_STEP, FORMAT_LONG, "" );
   }

   if ( END_T < 0.0 ) {
      END_T = End_T_Default;
      PRINT_RESET_PARA( END_T, FORMAT_REAL, "" );
   }

   if ( Merger_Coll_LabelCenter  &&  !OPT__RECORD_USER ) {
      OPT__RECORD_USER = true;
      PRINT_RESET_PARA( OPT__RECORD_USER, FORMAT_BOOL, "" );
   }


// (6) make a note
   if ( MPI_Rank == 0 )
   {
      Aux_Message( stdout, "=============================================================================\n" );
      Aux_Message( stdout, "  test problem ID        = %d\n",           TESTPROB_ID );
      Aux_Message( stdout, "  number of clusters     = %d\n",           Merger_Coll_NumHalos );
      Aux_Message( stdout, "  turn on AGN feedback   = %s\n",          (AGN_feedback)? "yes":"no" );
      if ( Merger_Coll_IsGas1 )
      Aux_Message( stdout, "  profile file 1         = %s\n",           Merger_File_Prof1 );
      Aux_Message( stdout, "  particle file 1        = %s\n",           Merger_File_Par1 );
      Aux_Message( stdout, "  cluster 1 w/ gas       = %s\n",          (Merger_Coll_IsGas1)? "yes":"no" );
      if ( Merger_Coll_IsGas1 )
      Aux_Message( stdout, "  cluster 1 x-position   = %g\n",           Merger_Coll_PosX1 );
      Aux_Message( stdout, "  cluster 1 y-position   = %g\n",           Merger_Coll_PosY1 );
      Aux_Message( stdout, "  cluster 1 x-velocity   = %g\n",           Merger_Coll_VelX1 );
      Aux_Message( stdout, "  cluster 1 y-velocity   = %g\n",           Merger_Coll_VelY1 );
      if ( Merger_Coll_NumHalos > 1 ) {
      if ( Merger_Coll_IsGas2 )
      Aux_Message( stdout, "  profile file 2         = %s\n",           Merger_File_Prof2 );
      Aux_Message( stdout, "  particle file 2        = %s\n",           Merger_File_Par2 );
      Aux_Message( stdout, "  cluster 2 w/ gas       = %s\n",          (Merger_Coll_IsGas2)? "yes":"no" );
      if ( Merger_Coll_IsGas2 )
      Aux_Message( stdout, "  cluster 2 x-position   = %g\n",           Merger_Coll_PosX2 );
      Aux_Message( stdout, "  cluster 2 y-position   = %g\n",           Merger_Coll_PosY2 );
      Aux_Message( stdout, "  cluster 2 x-velocity   = %g\n",           Merger_Coll_VelX2 );
      Aux_Message( stdout, "  cluster 2 y-velocity   = %g\n",           Merger_Coll_VelY2 );
      }
      if ( Merger_Coll_NumHalos > 2 ) {
      if ( Merger_Coll_IsGas3 )
      Aux_Message( stdout, "  profile file 3         = %s\n",           Merger_File_Prof3 );
      Aux_Message( stdout, "  particle file 3        = %s\n",           Merger_File_Par3 );
      Aux_Message( stdout, "  cluster 3 w/ gas       = %s\n",          (Merger_Coll_IsGas3)? "yes":"no" );
      if ( Merger_Coll_IsGas3 )
      Aux_Message( stdout, "  cluster 3 x-position   = %g\n",           Merger_Coll_PosX3 );
      Aux_Message( stdout, "  cluster 3 y-position   = %g\n",           Merger_Coll_PosY3 );
      Aux_Message( stdout, "  cluster 3 x-velocity   = %g\n",           Merger_Coll_VelX3 );
      Aux_Message( stdout, "  cluster 3 y-velocity   = %g\n",           Merger_Coll_VelY3 );
      }
      Aux_Message( stdout, "  use metals             = %s\n",          (Merger_Coll_UseMetals)? "yes":"no" );
      Aux_Message( stdout, "  label cluster centers  = %s\n",          (Merger_Coll_LabelCenter)? "yes":"no" );
      Aux_Message( stdout, "=============================================================================\n" );

//    Check if the accretion region is larger than the jet cylinder
      if ( R_acc < Jet_HalfHeight1 )  Aux_Message( stderr, "WARNING : R_acc is less than Jet_HalfHeight1!!\n");
      if ( R_acc < Jet_Radius1 )      Aux_Message( stderr, "WARNING : R_acc is less than Jet_Radius1!!\n");
      if ( Merger_Coll_NumHalos > 1 ) {
      if ( R_acc < Jet_HalfHeight2 )  Aux_Message( stderr, "WARNING : R_acc is less than Jet_HalfHeight2!!\n");
      if ( R_acc < Jet_Radius2 )      Aux_Message( stderr, "WARNING : R_acc is less than Jet_Radius2!!\n");
      }
      if ( Merger_Coll_NumHalos > 2 ) {
      if ( R_acc < Jet_HalfHeight3 )  Aux_Message( stderr, "WARNING : R_acc is less than Jet_HalfHeight3!!\n");
      if ( R_acc < Jet_Radius3 )      Aux_Message( stderr, "WARNING : R_acc is less than Jet_Radius3!!\n");
      }
   }

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "   Setting runtime parameters ... done\n" );

#endif // #ifdef SUPPORT_HDF5

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

   if ( !( Merger_Coll_IsGas1 || Merger_Coll_IsGas2 || Merger_Coll_IsGas3 ))
      return;

   const double ClusterCenter1[3] = { Merger_Coll_PosX1, Merger_Coll_PosY1, amr->BoxCenter[2] };
   const double ClusterCenter2[3] = { Merger_Coll_PosX2, Merger_Coll_PosY2, amr->BoxCenter[2] };
   const double ClusterCenter3[3] = { Merger_Coll_PosX3, Merger_Coll_PosY3, amr->BoxCenter[2] };

   double r1, r2, r3, Dens1, Dens2, Dens3, Pres1, Pres2, Pres3;
   double Metl1, Metl2, Metl3, rmax1, rmax2, rmax3;
   double Dens, MomX, MomY, MomZ, Pres, Eint, Etot, Metl;

   rmax1 = Table_R1[Merger_NBin1-1];
   rmax2 = Merger_Coll_NumHalos > 1 ? Table_R2[Merger_NBin2-1] : -1.0;
   rmax3 = Merger_Coll_NumHalos > 2 ? Table_R3[Merger_NBin3-1] : -1.0;

   //    for each cell, we sum up the density and pressure from each halo and then calculate the weighted velocity
   if ( Merger_Coll_IsGas1 ) {
      r1    = sqrt( SQR(x-ClusterCenter1[0]) + SQR(y-ClusterCenter1[1]) + SQR(z-ClusterCenter1[2]) );
      double rr1 = r1 < rmax1 ? r1 : rmax1;
      Dens1 = Mis_InterpolateFromTable( Merger_NBin1, Table_R1, Table_D1, rr1 );
      Pres1 = Mis_InterpolateFromTable( Merger_NBin1, Table_R1, Table_P1, rr1 );
      Metl1 = Mis_InterpolateFromTable( Merger_NBin1, Table_R1, Table_M1, rr1 );
   } else {
      r1    = HUGE_NUMBER;
      Dens1 = 0.0;
      Pres1 = 0.0;
      Metl1 = 0.0;
   }

   if ( Merger_Coll_NumHalos > 1 && Merger_Coll_IsGas2 ) {
      r2    = sqrt( SQR(x-ClusterCenter2[0]) + SQR(y-ClusterCenter2[1]) + SQR(z-ClusterCenter2[2]) );
      double rr2 = r2 < rmax2 ? r2 : rmax2;
      Dens2 = Mis_InterpolateFromTable( Merger_NBin2, Table_R2, Table_D2, rr2 );
      Pres2 = Mis_InterpolateFromTable( Merger_NBin2, Table_R2, Table_P2, rr2 );
      Metl2 = Mis_InterpolateFromTable( Merger_NBin2, Table_R2, Table_M2, rr2 );
   } else {
      r2    = HUGE_NUMBER;
      Dens2 = 0.0;
      Pres2 = 0.0;
      Metl2 = 0.0;
   }

   if ( Merger_Coll_NumHalos > 2 && Merger_Coll_IsGas3 ) {
      r3    = sqrt( SQR(x-ClusterCenter3[0]) + SQR(y-ClusterCenter3[1]) + SQR(z-ClusterCenter3[2]) );
      double rr3 = r3 < rmax3 ? r3 : rmax3;
      Dens3 = Mis_InterpolateFromTable( Merger_NBin3, Table_R3, Table_D3, rr3 );
      Pres3 = Mis_InterpolateFromTable( Merger_NBin3, Table_R3, Table_P3, rr3 );
      Metl3 = Mis_InterpolateFromTable( Merger_NBin3, Table_R3, Table_M3, rr3 );
   } else {
      r3    = HUGE_NUMBER;
      Dens3 = 0.0;
      Pres3 = 0.0;
      Metl3 = 0.0;
   }

   if ( Dens1 == NULL_REAL )
      Dens1 = Table_D1[Merger_NBin1-1];
   if ( Pres1 == NULL_REAL )
      Pres1 = Table_P1[Merger_NBin1-1];
   if ( Metl1 == NULL_REAL )
      Metl1 = Table_M1[Merger_NBin1-1];

   if ( Dens2 == NULL_REAL )
      Dens2 = Table_D2[Merger_NBin2-1];
   if ( Pres2 == NULL_REAL )
      Pres2 = Table_P2[Merger_NBin2-1];
   if ( Metl2 == NULL_REAL )
      Metl2 = Table_M2[Merger_NBin2-1];

   if ( Dens3 == NULL_REAL )
      Dens3 = Table_D3[Merger_NBin3-1];
   if ( Pres3 == NULL_REAL )
      Pres3 = Table_P3[Merger_NBin3-1];
   if ( Metl3 == NULL_REAL )
      Metl3 = Table_M3[Merger_NBin3-1];

   Dens = Dens1 + Dens2 + Dens3;
   Pres = Pres1 + Pres2 + Pres3;

   MomX = 0.0;
   MomY = 0.0;
   if ( r1 <= rmax1 ) {
      MomX += Merger_Coll_VelX1*Dens1;
      MomY += Merger_Coll_VelY1*Dens1;
   }
   if ( r2 <= rmax2 ) {
     MomX += Merger_Coll_VelX2*Dens2;
     MomY += Merger_Coll_VelY2*Dens2;
   }
   if ( r3 <= rmax3 ) {
     MomX += Merger_Coll_VelX3*Dens3;
     MomY += Merger_Coll_VelY3*Dens3;
   }
   MomZ = 0.0;

   // compute the total gas energy
   Eint = EoS_DensPres2Eint_CPUPtr( Dens, Pres, NULL, EoS_AuxArray_Flt,
                                    EoS_AuxArray_Int, h_EoS_Table );    // assuming EoS requires no passive scalars

   Etot = Hydro_ConEint2Etot( Dens, MomX, MomY, MomZ, Eint, 0.0 );      // do NOT include magnetic energy here

   fluid[DENS] = Dens;
   fluid[MOMX] = MomX;
   fluid[MOMY] = MomY;
   fluid[MOMZ] = MomZ;
   fluid[ENGY] = Etot;

   if ( Merger_Coll_UseMetals ) {
      Metl = Metl1*Dens1 + Metl2*Dens2 + Metl3*Dens3;
      fluid[Idx_Metal] = Metl;
   }

   if ( Merger_Coll_IsGas1 ) {
      if ( r1 < rmax1 )
         fluid[ColorField1Idx] = Dens1;
      else
         fluid[ColorField1Idx] = 0.0;
   }

   if ( Merger_Coll_NumHalos > 1 && Merger_Coll_IsGas2 ) {
      if ( r2 < rmax2 )
         fluid[ColorField2Idx] = Dens2;
      else
         fluid[ColorField2Idx] = 0.0;
   }

   if ( Merger_Coll_NumHalos > 2 && Merger_Coll_IsGas3 ) {
      if ( r3 < rmax3 )
         fluid[ColorField3Idx] = Dens3;
      else
         fluid[ColorField3Idx] = 0.0;
   }

} // FUNCTION : SetGridIC

#endif // #if ( MODEL == HYDRO  &&  defined MASSIVE_PARTICLES )


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

   delete [] Table_R1;
   delete [] Table_D1;
   delete [] Table_P1;
   delete [] Table_M1;

   delete [] Table_R2;
   delete [] Table_D2;
   delete [] Table_P2;
   delete [] Table_M2;

   delete [] Table_R3;
   delete [] Table_D3;
   delete [] Table_P3;
   delete [] Table_M3;

} // FUNCTION : End_ClusterMerger

#ifdef MHD

void SetBFieldIC( real magnetic[], const double x, const double y, const double z, const double Time,
                  const int lv, double AuxArray[] )
{

   // Setting the B-field in this problem is always handled by reading from a uniform
   // grid in the file B_IC

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

   Flu_ResetByUser_Func_Ptr       = Flu_ResetByUser_Func_ClusterMerger;
   Flu_ResetByUser_API_Ptr        = Flu_ResetByUser_API_ClusterMerger;
   Output_User_Ptr                = Output_ClusterMerger;
   Init_User_Ptr                  = Init_User_ClusterMerger;

#  ifdef MHD
   Init_Function_BField_User_Ptr  = SetBFieldIC;
#  endif
#  endif // if ( MODEL == HYDRO  &&  defined MASSIVE_PARTICLES )

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : Init_TestProb_Hydro_ClusterMerger

#ifdef SUPPORT_HDF5

int Read_Num_Points_ClusterMerger(std::string filename)
{

   hid_t   file_id, dataset, dataspace;
   herr_t  status;
   hsize_t dims[1], maxdims[1];

   int rank;

   file_id = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);

   dataset   = H5Dopen(file_id, "/fields/radius", H5P_DEFAULT);
   dataspace = H5Dget_space(dataset);
   rank      = H5Sget_simple_extent_dims(dataspace, dims, maxdims);

   H5Fclose(file_id);

   return (int)dims[0];

} // FUNCTION : Read_Num_Points_ClusterMerger

void Read_Profile_ClusterMerger(std::string filename, std::string fieldname,
                                double field[])
{

   hid_t   file_id, dataset;
   herr_t  status;

   file_id = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);

   dataset = H5Dopen(file_id, fieldname.c_str(), H5P_DEFAULT);
   status = H5Dread(dataset, H5T_NATIVE_DOUBLE, H5S_ALL,
                     H5S_ALL, H5P_DEFAULT, field);
   H5Dclose(dataset);

   H5Fclose(file_id);

   return;

} // FUNCTION : Read_Profile_ClusterMerger

long Read_Particle_Number_ClusterMerger(std::string filename)
{

   hid_t   file_id, dataset, dataspace;
   herr_t  status;
   hsize_t dims[1], maxdims[1];

   int rank;

   file_id = H5Fopen(filename.c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);

   dataset   = H5Dopen(file_id, "particle_mass", H5P_DEFAULT);
   dataspace = H5Dget_space(dataset);
   rank      = H5Sget_simple_extent_dims(dataspace, dims, maxdims);

   H5Sclose(dataspace);
   H5Dclose(dataset);
   H5Fclose(file_id);

   return (long)dims[0];

} // FUNCTION : Read_Particle_Number_ClusterMerger

#endif // #ifdef SUPPORT_HDF5

#if ( MODEL == HYDRO )
void AddNewField_ClusterMerger()
{

   if ( Merger_Coll_UseMetals )
      Idx_Metal = AddField( "Metal", FIXUP_FLUX_YES, FIXUP_REST_YES, NORMALIZE_NO, INTERP_FRAC_NO );
   if ( ColorField1Idx == Idx_Undefined )
      ColorField1Idx = AddField( "ColorField1", FIXUP_FLUX_YES, FIXUP_REST_YES, NORMALIZE_NO, INTERP_FRAC_NO );
   if ( Merger_Coll_NumHalos > 1 && ColorField2Idx == Idx_Undefined )
      ColorField2Idx = AddField( "ColorField2", FIXUP_FLUX_YES, FIXUP_REST_YES, NORMALIZE_NO, INTERP_FRAC_NO );
   if ( Merger_Coll_NumHalos > 2 && ColorField3Idx == Idx_Undefined )
      ColorField3Idx = AddField( "ColorField3", FIXUP_FLUX_YES, FIXUP_REST_YES, NORMALIZE_NO, INTERP_FRAC_NO );

}

#endif


// Restart only: initialize the BH variables (position, velocity, mass) from particles labelled with PTYPE_CEN 
void Init_User_ClusterMerger()
{
   if ( OPT__INIT == INIT_BY_RESTART ){

      if ( OPT__RESTART_RESET ) {
         printf("Error! OPT__RESTART_RESET should be disabled.\n");
         return;
      }    

      const char FileName[] = "BH_variable.bin";
      int TargetDumpID = DumpID-1;  //INIT_DUMPID;

      FILE* File_User = fopen(FileName, "rb");
      if ( File_User == NULL ) {
         Aux_Error( ERROR_INFO, "Error opening the file \"%s\"\n", FileName ); 
         return;
      }    

      if ( MPI_Rank == 0 ) {
   
         double BH_Mass[3] = {0.0, 0.0, 0.0};
   
         double Time;
         int dumpID;
         while (fread(&Time, sizeof(double), 1, File_User) == 1) {
            fread(&dumpID, sizeof(int), 1, File_User);
            printf("dumpID = %d, TargetDumpID = %d\n", dumpID, TargetDumpID);  
            if (dumpID == TargetDumpID) {
               fread(&Merger_Coll_NumHalos, sizeof(int), 1, File_User);
               for (int c=0; c<Merger_Coll_NumHalos; c++) {
                  for (int d=0; d<3; d++)   fread(&BH_Pos[c][d], sizeof(double), 1, File_User);
                  for (int d=0; d<3; d++)   fread(&ClusterCen[c][d], sizeof(double), 1, File_User);
                  for (int d=0; d<3; d++)   fread(&BH_Vel[c][d], sizeof(double), 1, File_User);
                  fread(&BH_Mass[c], sizeof(double), 1, File_User);
               }
               fread(&AdjustCount, sizeof(int), 1, File_User);
               break;
            }
            else {
               fseek(File_User, sizeof(int), SEEK_CUR);
               for (int c=0; c<Merger_Coll_NumHalos; c++) {
                  fseek(File_User, sizeof(double) * 3, SEEK_CUR);
                  fseek(File_User, sizeof(double) * 3, SEEK_CUR);
                  fseek(File_User, sizeof(double) * 3, SEEK_CUR);
                  fseek(File_User, sizeof(double), SEEK_CUR);
               }
               fseek(File_User, sizeof(int), SEEK_CUR);
            }
         }
   
         Bondi_MassBH1 = BH_Mass[0];
         Bondi_MassBH2 = BH_Mass[1];
         Bondi_MassBH3 = BH_Mass[2];
   
         fclose(File_User);

         printf("Restarting! BH_Pos[0][0] = %23.17e, BH_Pos[0][1] = %23.17e\n", BH_Pos[0][0], BH_Pos[0][1]);  
      }  // if ( MPI_Rank == 0 )
      MPI_Bcast( &Merger_Coll_NumHalos, 1, MPI_INT, 0, MPI_COMM_WORLD );
      for (int c=0; c<Merger_Coll_NumHalos; c++){ 
         MPI_Bcast( BH_Pos[c], 3, MPI_DOUBLE, 0, MPI_COMM_WORLD ); 
         MPI_Bcast( BH_Vel[c], 3, MPI_DOUBLE, 0, MPI_COMM_WORLD );
         MPI_Bcast( ClusterCen[c], 3, MPI_DOUBLE, 0, MPI_COMM_WORLD );
      }
      MPI_Bcast( &Bondi_MassBH1, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD ); 
      MPI_Bcast( &Bondi_MassBH2, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD ); 
      MPI_Bcast( &Bondi_MassBH3, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD ); 
      MPI_Bcast( &AdjustCount, 1, MPI_INT, 0, MPI_COMM_WORLD ); 

   } // if ( OPT__INIT == INIT_BY_RESTART )

   else   return;
}



//-------------------------------------------------------------------------------------------------------
// Function    :  OutputClusterMerger
// Description :  Output the BH variables (pos, vel, mass, AdjustCount) in the ClusterMerger problem
//
// Note        :  1. Enabled by the runtime option "OPT__OUTPUT_USER"
//                2. Ensure restart runs can get identical values 
//
// Parameter   :  None
//
// Return      :  None
//-------------------------------------------------------------------------------------------------------
void Output_ClusterMerger()
{ 
   const char FileName[] = "BH_variable.bin";
   static bool FirstTime = true;

   if ( FirstTime ) {
      if ( MPI_Rank == 0 ) {
         if ( Aux_CheckFileExist(FileName) )
            Aux_Message( stderr, "WARNING: file \"%s\" already exists !!\n", FileName );
      }
      FirstTime = false;             
   } // if ( FirstTime )     

   double BH_Mass[3] = { Bondi_MassBH1, Bondi_MassBH2, Bondi_MassBH3 };

   if ( MPI_Rank == 0 ) {
      FILE* File_User = fopen(FileName, "ab");
      
      double time = Time[0] * UNIT_T / Const_Myr;
      fwrite(&time, sizeof(double), 1, File_User);
      fwrite(&DumpID, sizeof(int), 1, File_User);
      fwrite(&Merger_Coll_NumHalos, sizeof(int), 1, File_User);

      for (int c=0; c<Merger_Coll_NumHalos; c++) {      
         for (int d=0; d<3; d++)   fwrite(&BH_Pos[c][d], sizeof(double), 1, File_User);
         for (int d=0; d<3; d++)   fwrite(&ClusterCen[c][d], sizeof(double), 1, File_User);
         for (int d=0; d<3; d++)   fwrite(&BH_Vel[c][d], sizeof(double), 1, File_User); 
         fwrite(&BH_Mass[c], sizeof(double), 1, File_User);
      }
      fwrite(&AdjustCount, sizeof(int), 1, File_User);
      fclose(File_User);

      printf("Saving BH_variables! DumpID = %d, BH_Pos[0][0] = %23.17e, BH_Pos[0][1] = %23.17e\n", DumpID, BH_Pos[0][0], BH_Pos[0][1]);
   }

} // FUNCTION : Output_ClusterMerger
