#ifndef __HDF5_TYPEDEF_H__
#define __HDF5_TYPEDEF_H__


/*===========================================================================
Data structures defined here are mainly used for creating the compound
datatypes in the HDF5 format
===========================================================================*/

#include "hdf5.h"
#include "Macro.h"
#include "CUFLU.h"
#ifdef GRAVITY
#include "CUPOT.h"
#endif

#ifdef FLOAT8
#  define H5T_GAMER_REAL H5T_NATIVE_DOUBLE
#else
#  define H5T_GAMER_REAL H5T_NATIVE_FLOAT
#endif

#ifdef GAMER_DEBUG
#  define DEBUG_HDF5
#endif

#ifdef PARTICLE
# ifdef FLOAT8_PAR
#  define H5T_GAMER_REAL_PAR H5T_NATIVE_DOUBLE
# else
#  define H5T_GAMER_REAL_PAR H5T_NATIVE_FLOAT
# endif
#endif // #ifdef PARTICLE




//-------------------------------------------------------------------------------------------------------
// Structure   :  KeyInfo_t
// Description :  Data structure for outputting the key simulation information
//
// Note        :  1. Some variables have different names from their internal names in the code
//                   --> Make the output file more readable
//                2. In general, this data structure stores simulation information that cannot be determined
//                   from the input parameter files (i.e., Input__* files)
//                   --> Therefore, this data structure is crucial for the restart process
//-------------------------------------------------------------------------------------------------------
struct KeyInfo_t
{

   int    FormatVersion;
   int    Model;
   int    Float8;
   int    Gravity;
   int    Particle;
   int    NLevel;
   int    NCompFluid;               // NCOMP_FLUID
   int    NCompPassive;             // NCOMP_PASSIVE
   int    PatchSize;
   int    DumpID;
   int    NX0     [3];
   int    BoxScale[3];
   int    NPatch   [NLEVEL];
   int    CellScale[NLEVEL];        // amr->scale[lv]
#  if ( MODEL == HYDRO )
   int    Magnetohydrodynamics;
   int    SRHydrodynamics;
   int    CosmicRay;
#  endif

   long   Step;
   long   AdvanceCounter[NLEVEL];
   int    NFieldStored;             // number of grid fields to be stored (excluding B field)
   int    NMagStored;               // NCOMP_MAG (declare it even when MHD is off)
#  ifdef PARTICLE
   long   Par_NPar;                 // amr->Par->NPar_Active_AllRank
   int    Par_NAttStored;           // PAR_NATT_STORED
   int    Float8_Par;
#  endif
#  ifdef COSMIC_RAY
   int    CR_Diffusion;
#  endif

   double BoxSize[3];
   double Time       [NLEVEL];
   double CellSize   [NLEVEL];      // amr->dh[lv]
   double dTime_AllLv[NLEVEL];
#  ifdef GRAVITY
   double AveDens_Init;             // AveDensity_Init
#  endif

   char  *CodeVersion;
   char  *DumpWallTime;
   char  *GitBranch;
   char  *GitCommit;
   long   UniqueDataID;

}; // struct KeyInfo_t



//-------------------------------------------------------------------------------------------------------
// Structure   :  Makefile_t
// Description :  Data structure for outputting the Makefile options
//
// Note        :  1. Some options are output only when necessary
//-------------------------------------------------------------------------------------------------------
struct Makefile_t
{

   int Model;
   int Gravity;
   int Comoving;
   int Particle;
   int NLevel;
   int MaxPatch;

   int UseGPU;
   int GAMER_Debug;
   int BitwiseReproducibility;
   int Timing;
   int TimingSolver;
   int Float8;
   int Serial;
   int LoadBalance;
   int OverlapMPI;
   int OpenMP;
   int GPU_Arch;
   int Laohu;
   int SupportHDF5;
   int SupportGSL;
   int SupportFFTW;
   int SupportLibYT;
#  ifdef SUPPORT_LIBYT
   int LibYTUsePatchGroup;
   int LibYTInteractive;
   int LibYTReload;
   int LibYTJupyter;
#  endif
   int SupportGrackle;
   int RandomNumber;

#  ifdef GRAVITY
   int PotScheme;
   int StorePotGhost;
   int UnsplitGravity;
#  endif

#  if   ( MODEL == HYDRO )
   int FluScheme;
#  ifdef LR_SCHEME
   int LRScheme;
#  endif
#  ifdef RSOLVER
   int RSolver;
#  endif
   int DualEnergy;
   int Magnetohydrodynamics;
   int SRHydrodynamics;
   int CosmicRay;
   int EoS;
   int BarotropicEoS;

#  elif ( MODEL == ELBDM )
   int ConserveMass;
   int Laplacian4th;
   int SelfInteraction4;

#  else
#  error : unsupported MODEL !!
#  endif // MODEL

#  ifdef PARTICLE
   int MassiveParticles;
   int Tracer;
   int StoreParAcc;
   int StarFormation;
   int Feedback;
   int Par_NAttUser;
   int Float8_Par;
#  endif

#  ifdef COSMIC_RAY
   int CR_Diffusion;
#  endif

}; // struct Makefile_t



//-------------------------------------------------------------------------------------------------------
// Structure   :  SymConst_t
// Description :  Data structure for outputting the symbolic constants
//
// Note        :  1. Symbolic constants are defined in "Macro.h, CUFLU.h, CUPOT.h, Particle.h"
//-------------------------------------------------------------------------------------------------------
struct SymConst_t
{

   int    NCompFluid;
   int    NCompPassive;
   int    PatchSize;
   int    Flu_NIn;
   int    Flu_NOut;
   int    Flu_NIn_T;
   int    Flu_NIn_S;
   int    Flu_NOut_S;
   int    NFluxFluid;
   int    NFluxPassive;
   int    Flu_GhostSize;
   int    Flu_Nxt;
   int    Debug_HDF5;
   int    SibOffsetNonperiodic;

#  ifdef LOAD_BALANCE
   int    SonOffsetLB;
#  endif

   double TinyNumber;
   double HugeNumber;
   double MaxError;


#  ifdef GRAVITY
   int    Gra_NIn;
   int    Pot_GhostSize;
   int    Gra_GhostSize;
   int    Rho_GhostSize;
   int    Pot_Nxt;
   int    Gra_Nxt;
   int    Rho_Nxt;

#  ifdef UNSPLIT_GRAVITY
   int    USG_GhostSizeF;
   int    USG_GhostSizeG;
   int    USG_NxtF;
   int    USG_NxtG;
#  endif

   int    ExtPot_BlockSize;
   int    Gra_BlockSize;
   int    ExtPotNAuxMax;
   int    ExtAccNAuxMax;
   int    ExtPotNGeneMax;

#  if   ( POT_SCHEME == SOR )
   int    Pot_BlockSize_z;
   int    SOR_RhoShared;
   int    SOR_CPotShared;
   int    SOR_UseShuffle;
   int    SOR_UsePadding;
   int    SOR_ModReduction;
#  elif ( POT_SCHEME == MG  )
   int    Pot_BlockSize_x;
#  endif

#  endif // #ifdef GRAVITY


#  ifdef PARTICLE
   int    Par_NAttStored;
   int    Par_NType;
#  ifdef GRAVITY
   int    RhoExt_GhostSize;
#  endif
   int    Debug_Particle;

   double ParList_GrowthFactor;
   double ParList_ReduceFactor;
#  endif


   int    BitRep_Flux;
#  ifdef MHD
   int    BitRep_Electric;
#  endif

   int    InterpMask;
   int    FB_SepFluOut;


#  if   ( MODEL == HYDRO )
   int    Flu_BlockSize_x;
   int    Flu_BlockSize_y;
   int    CheckUnphyInFluid;
   int    CharReconstruction;
   int    LR_Eint;
   int    CheckIntermediate;
   int    RSolverRescue;
   int    HLL_NoRefState;
   int    HLL_IncludeAllWaves;
   int    HLLC_WaveSpeed;
   int    HLLE_WaveSpeed;
#  ifdef MHD
   int    HLLD_WaveSpeed;
#  endif
#  ifdef N_FC_VAR
   int    N_FC_Var;
#  endif
#  ifdef N_SLOPE_PPM
   int    N_Slope_PPM;
#  endif
#  ifdef MHD
   int    EulerY;
#  endif
   int    MHM_CheckPredict;
   int    EoSNAuxMax;
   int    EoSNTableMax;

#  elif  ( MODEL == ELBDM )
   int    Flu_BlockSize_x;
   int    Flu_BlockSize_y;

#  else
#  error : ERROR : unsupported MODEL !!
#  endif // MODEL


   int    dt_Flu_BlockSize;
   int    dt_Flu_UseShuffle;
#  ifdef GRAVITY
   int    dt_Gra_BlockSize;
   int    dt_Gra_UseShuffle;
#  endif

   int    Src_BlockSize;
   int    Src_GhostSize;
   int    Src_Nxt;
   int    Src_NAuxDlep;
   int    Src_DlepProfNVar;
   int    Src_DlepProfNBinMax;
   int    Src_NAuxUser;

   int    Der_GhostSize;
   int    Der_Nxt;
   int    Der_NOut_Max;

#  ifdef FEEDBACK
   int    FB_GhostSize;
   int    FB_Nxt;
#  endif

   int    NFieldStoredMax;

}; // struct SymConst_t




//-------------------------------------------------------------------------------------------------------
// Structure   :  InputPara_t
// Description :  Data structure for outputting the run-time parameters
//
// Note        :  1. All run-time parameters are loaded from the files "Input__XXX"
//-------------------------------------------------------------------------------------------------------
struct InputPara_t
{

// simulation scale
   double BoxSize;
   int    NX0_Tot[3];
   int    MPI_NRank;
   int    MPI_NRank_X[3];
   int    OMP_NThread;
   double EndT;
   long   EndStep;

// test problems
   int   TestProb_ID;

// code units
   int    Opt__Unit;
   double Unit_L;
   double Unit_M;
   double Unit_T;
   double Unit_V;
   double Unit_D;
   double Unit_E;
   double Unit_P;
#  ifdef MHD
   double Unit_B;
#  endif

// boundary condition
   int    Opt__BC_Flu[6];
#  ifdef GRAVITY
   int    Opt__BC_Pot;
   double GFunc_Coeff0;
#  endif

// particle
#  ifdef PARTICLE
   int    Par_Init;
   int    Par_ICFormat;
   double Par_ICMass;
   int    Par_ICType;
   int    Par_ICFloat8;
   int    Par_Interp;
   int    Par_InterpTracer;
   int    Par_Integ;
   int    Par_IntegTracer;
   int    Par_ImproveAcc;
   int    Par_PredictPos;
   double Par_RemoveCell;
   int    Opt__FreezePar;
   int    Par_GhostSize;
   int    Par_GhostSizeTracer;
   int    Par_TracerVelCorr;
   char  *ParAttLabel[PAR_NATT_TOTAL];
#  endif

// cosmology
#  ifdef COMOVING
   double A_Init;
   double OmegaM0;
   double Hubble0;
#  endif

// time-step determination
   double Dt__Max;
   double Dt__Fluid;
   double Dt__FluidInit;
#  ifdef GRAVITY
   double Dt__Gravity;
#  endif
#  if ( MODEL == ELBDM )
   double Dt__Phase;
#  endif
#  ifdef PARTICLE
   double Dt__ParVel;
   double Dt__ParVelMax;
   double Dt__ParAcc;
#  endif
#  ifdef CR_DIFFUSION
   double Dt__CR_Diffusion;
#  endif
#  ifdef COMOVING
   double Dt__MaxDeltaA;
#  endif
#  ifdef SRHD
   int    Dt__SpeedOfLight;
#  endif
   double Dt__SyncParentLv;
   double Dt__SyncChildrenLv;
   int    Opt__DtUser;
   int    Opt__DtLevel;
   int    Opt__RecordDt;
   int    AutoReduceDt;
   double AutoReduceDtFactor;
   double AutoReduceDtFactorMin;
#  if ( MODEL == HYDRO )
   double AutoReduceMinModFactor;
   double AutoReduceMinModMin;
#  endif
   double AutoReduceIntMonoFactor;
   double AutoReduceIntMonoMin;

// domain refinement
   int    RegridCount;
   int    RefineNLevel;
   int    FlagBufferSize;
   int    FlagBufferSizeMaxM1Lv;
   int    FlagBufferSizeMaxM2Lv;
   int    MaxLevel;
   int    Opt__Flag_Rho;
   int    Opt__Flag_RhoGradient;
#  if ( MODEL == HYDRO )
   int    Opt__Flag_PresGradient;
   int    Opt__Flag_Vorticity;
   int    Opt__Flag_Jeans;
#  ifdef MHD
   int    Opt__Flag_Current;
#  endif
#  ifdef SRHD
   int    Opt__Flag_LrtzGradient;
#  endif
#  ifdef COSMIC_RAY
   int    Opt__Flag_CRay;
#  endif
#  endif
#  if ( MODEL == ELBDM )
   int    Opt__Flag_EngyDensity;
#  endif
   int    Opt__Flag_LohnerDens;
#  if ( MODEL == HYDRO )
   int    Opt__Flag_LohnerEngy;
   int    Opt__Flag_LohnerPres;
   int    Opt__Flag_LohnerTemp;
   int    Opt__Flag_LohnerEntr;
#  ifdef COSMIC_RAY
   int    Opt__Flag_LohnerCRay;
#  endif
#  endif
   int    Opt__Flag_LohnerForm;
   int    Opt__Flag_User;
   int    Opt__Flag_User_Num;
   int    Opt__Flag_Region;
#  ifdef PARTICLE
   int    Opt__Flag_NParPatch;
   int    Opt__Flag_NParCell;
   int    Opt__Flag_ParMassCell;
#  endif
   int    Opt__NoFlagNearBoundary;
   int    Opt__PatchCount;
#  ifdef PARTICLE
   int    Opt__ParticleCount;
#  endif
   int    Opt__ReuseMemory;
   int    Opt__MemoryPool;

// load balance
#  ifdef LOAD_BALANCE
   double LB_WLI_Max;
#  ifdef PARTICLE
   double LB_Par_Weight;
#  endif
   int    Opt__RecordLoadBalance;
#  endif
   int    Opt__MinimizeMPIBarrier;

// fluid solvers in HYDRO
#  if ( MODEL == HYDRO )
   double Gamma;
   double MolecularWeight;
   double MuNorm;
   double IsoTemp;
   double MinMod_Coeff;
   int    MinMod_MaxIter;
   int    Opt__LR_Limiter;
   int    Opt__1stFluxCorr;
   int    Opt__1stFluxCorrScheme;
#  ifdef DUAL_ENERGY
   double DualEnergySwitch;
#  endif
#  ifdef MHD
   int    Opt__SameInterfaceB;
#  endif
#  endif // HYDRO

// ELBDM solvers
#  if ( MODEL == ELBDM )
   double ELBDM_Mass;
   double ELBDM_PlanckConst;
#  ifdef QUARTIC_SELF_INTERACTION
   double ELBDM_Lambda;
#  endif
   double ELBDM_Taylor3_Coeff;
   int    ELBDM_Taylor3_Auto;
#  endif // ELBDM

// fluid solvers in different models
   int    Flu_GPU_NPGroup;
   int    GPU_NStream;
   int    Opt__FixUp_Flux;
   long   FixUpFlux_Var;
#  ifdef MHD
   int    Opt__FixUp_Electric;
#  endif
   int    Opt__FixUp_Restrict;
   long   FixUpRestrict_Var;
   int    Opt__CorrAfterAllSync;
   int    Opt__NormalizePassive;
   int    NormalizePassive_NVar;
   int    NormalizePassive_VarIdx[NCOMP_PASSIVE];
   int    Opt__IntFracPassive_LR;
   int    IntFracPassive_NVar;
   int    IntFracPassive_VarIdx[NCOMP_PASSIVE];
   char  *FieldLabel[NFIELD_STORED_MAX];
#  ifdef MHD
   char  *MagLabel[NCOMP_MAG];
#  endif
   int    Opt__OverlapMPI;
   int    Opt__ResetFluid;
   int    Opt__ResetFluidInit;
   int    Opt__FreezeFluid;
#  if ( MODEL == HYDRO  ||  MODEL == ELBDM )
   double MinDens;
#  endif
#  if ( MODEL == HYDRO )
   double MinPres;
   double MinEint;
   double MinTemp;
   double MinEntr;
   int    Opt__CheckPresAfterFlu;
   int    Opt__LastResortFloor;
   int    JeansMinPres;
   int    JeansMinPres_Level;
   int    JeansMinPres_NCell;
#  endif

// gravity
#  ifdef GRAVITY
   double NewtonG;
#  if   ( POT_SCHEME == SOR )
   double SOR_Omega;
   int    SOR_MaxIter;
   int    SOR_MinIter;
#  elif ( POT_SCHEME == MG )
   int    MG_MaxIter;
   int    MG_NPreSmooth;
   int    MG_NPostSmooth;
   double MG_ToleratedError;
#  endif
   int    Pot_GPU_NPGroup;
   int    Opt__GraP5Gradient;
   int    Opt__SelfGravity;
   int    Opt__ExtAcc;
   int    Opt__ExtPot;
   char  *ExtPotTable_Name;
   int    ExtPotTable_NPoint[3];
   double ExtPotTable_dh[3];
   double ExtPotTable_EdgeL[3];
   int    ExtPotTable_Float8;
   int    Opt__GravityExtraMass;
#  endif // #ifdef GRAVITY

// source terms
   int    Src_Deleptonization;
   int    Src_User;
   int    Src_GPU_NPGroup;

// Grackle
#  ifdef SUPPORT_GRACKLE
   int    Grackle_Activate;
   int    Grackle_Verbose;
   int    Grackle_Cooling;
   int    Grackle_Primordial;
   int    Grackle_Metal;
   int    Grackle_UV;
   int    Grackle_CMB_Floor;
   int    Grackle_PE_Heating;
   double Grackle_PE_HeatingRate;
   char  *Grackle_CloudyTable;
   int    Grackle_ThreeBodyRate;
   int    Grackle_CIE_Cooling;
   int    Grackle_H2_OpaApprox;
   int    Che_GPU_NPGroup;
#  endif

// star formation
#  ifdef STAR_FORMATION
   int    SF_CreateStar_Scheme;
   int    SF_CreateStar_RSeed;
   int    SF_CreateStar_DetRandom;
   int    SF_CreateStar_MinLevel;
   double SF_CreateStar_MinGasDens;
   double SF_CreateStar_MassEff;
   double SF_CreateStar_MinStarMass;
   double SF_CreateStar_MaxStarMFrac;
#  endif

// feedback
#  ifdef FEEDBACK
   int   FB_Level;
   int   FB_RSeed;
   int   FB_SNe;
   int   FB_User;
#  endif

// cosmic ray
#  ifdef COSMIC_RAY
   double CR_Gamma;
#  endif

// microphysics
#  ifdef CR_DIFFUSION
   double CR_Diffusion_ParaCoeff;
   double CR_Diffusion_PerpCoeff;
   double CR_Diffusion_MinB;
#  endif

// initialization
   int    Opt__Init;
   int    RestartLoadNRank;
   int    Opt__RestartReset;
   int    Opt__UM_IC_Level;
   int    Opt__UM_IC_NLevel;
   int    Opt__UM_IC_NVar;
   int    Opt__UM_IC_Format;
   int    Opt__UM_IC_Float8;
   int    Opt__UM_IC_Downgrade;
   int    Opt__UM_IC_Refine;
   int    Opt__UM_IC_LoadNRank;
   int    UM_IC_RefineRegion[NLEVEL-1][6];
   int    Opt__InitRestrict;
   int    Opt__InitGridWithOMP;
   int    Opt__GPUID_Select;
   int    Init_Subsampling_NCell;
#  ifdef MHD
   int    Opt__InitBFieldByVecPot;
#  endif
#  ifdef SUPPORT_FFTW
   int    Opt__FFTW_Startup;
#  endif

// interpolation schemes
   int    Opt__Int_Time;
#  if ( MODEL == HYDRO )
   int    Opt__Int_Prim;
#  endif
#  if ( MODEL == ELBDM )
   int    Opt__Int_Phase;
#  endif
   int    Opt__Flu_IntScheme;
   int    Opt__RefFlu_IntScheme;
#  ifdef MHD
   int    Opt__Mag_IntScheme;
   int    Opt__RefMag_IntScheme;
#  endif
#  ifdef GRAVITY
   int    Opt__Pot_IntScheme;
   int    Opt__Rho_IntScheme;
   int    Opt__Gra_IntScheme;
   int    Opt__RefPot_IntScheme;
#  endif
   double IntMonoCoeff;
#  ifdef MHD
   double IntMonoCoeffB;
#  endif
   int    Mono_MaxIter;
   int    IntOppSign0thOrder;

// data dump
   int    Opt__Output_Total;
   int    Opt__Output_Part;
   int    Opt__Output_User;
#  ifdef PARTICLE
   int    Opt__Output_Par_Mode;
#  endif
   int    Opt__Output_BasePS;
   int    Opt__Output_Base;
#  ifdef MHD
   int    Opt__Output_CC_Mag;
#  endif
#  ifdef GRAVITY
   int    Opt__Output_Pot;
#  endif
#  ifdef PARTICLE
   int    Opt__Output_ParDens;
#  endif
#  if ( MODEL == HYDRO )
   int    Opt__Output_Pres;
   int    Opt__Output_Temp;
   int    Opt__Output_Entr;
   int    Opt__Output_Cs;
   int    Opt__Output_DivVel;
   int    Opt__Output_Mach;
#  ifdef MHD
   int    Opt__Output_DivMag;
#  endif
#  ifdef SRHD
   int    Opt__Output_Lorentz;
   int    Opt__Output_3Velocity;
   int    Opt__Output_Enthalpy;
#  endif
#  endif // #if ( MODEL == HYDRO )
   int    Opt__Output_UserField;
   int    Opt__Output_Mode;
   int    Opt__Output_Restart;
   int    Opt__Output_Step;
   double Opt__Output_Dt;
   char  *Opt__Output_Text_Format_Flt;
   double Output_PartX;
   double Output_PartY;
   double Output_PartZ;
   int    InitDumpID;

// libyt jupyter interface
#  if ( defined(SUPPORT_LIBYT) && defined(LIBYT_JUPYTER) )
   int    Yt_JupyterUseConnectionFile;
#  endif

// miscellaneous
   int    Opt__Verbose;
   int    Opt__TimingBarrier;
   int    Opt__TimingBalance;
   int    Opt__TimingMPI;
   int    Opt__RecordNote;
   int    Opt__RecordUnphy;
   int    Opt__RecordMemory;
   int    Opt__RecordPerformance;
   int    Opt__ManualControl;
   int    Opt__RecordCenter;
   double COM_CenX;
   double COM_CenY;
   double COM_CenZ;
   double COM_MaxR;
   double COM_MinRho;
   double COM_TolErrR;
   int    COM_MaxIter;
   int    Opt__RecordUser;
   int    Opt__OptimizeAggressive;
   int    Opt__SortPatchByLBIdx;

// simulation checks
   int    Opt__Ck_Refine;
   int    Opt__Ck_ProperNesting;
   int    Opt__Ck_Conservation;
   double AngMom_OriginX;
   double AngMom_OriginY;
   double AngMom_OriginZ;
   int    Opt__Ck_NormPassive;
   int    Opt__Ck_Restrict;
   int    Opt__Ck_Finite;
   int    Opt__Ck_PatchAllocate;
   int    Opt__Ck_FluxAllocate;
#  if ( MODEL == HYDRO )
   int    Opt__Ck_Negative;
#  endif
   double Opt__Ck_MemFree;
#  ifdef PARTICLE
   int    Opt__Ck_Particle;
#  endif
#  ifdef MHD
   int    Opt__Ck_InterfaceB;
   int    Opt__Ck_DivergenceB;
#  endif
   int    Opt__Ck_InputFluid;

// flag tables
   double FlagTable_Rho         [NLEVEL-1];
   double FlagTable_RhoGradient [NLEVEL-1];
   double FlagTable_Lohner      [NLEVEL-1][5];
   hvl_t  FlagTable_User        [NLEVEL-1];
#  if   ( MODEL == HYDRO )
   double FlagTable_PresGradient[NLEVEL-1];
   double FlagTable_Vorticity   [NLEVEL-1];
   double FlagTable_Jeans       [NLEVEL-1];
#  ifdef MHD
   double FlagTable_Current     [NLEVEL-1];
#  endif
#  ifdef SRHD
   double FlagTable_LrtzGradient[NLEVEL-1];
#  endif
#  ifdef COSMIC_RAY
   double FlagTable_CRay        [NLEVEL-1];
#  endif
#  elif ( MODEL == ELBDM )
   double FlagTable_EngyDensity [NLEVEL-1][2];
#  endif
#  ifdef PARTICLE
   int    FlagTable_NParPatch   [NLEVEL-1];
   int    FlagTable_NParCell    [NLEVEL-1];
   double FlagTable_ParMassCell [NLEVEL-1];
#  endif

// user-defined derived fields
// --> always allocate DER_NOUT_MAX labels and units but only record UserDerField_Num of them in HDF5
// --> more convenient since storing dynamic arrays such as (*UserDerField_Label)[MAX_STRING] in HDF5 can be tricky
   int    UserDerField_Num;
   char  *UserDerField_Label[DER_NOUT_MAX];
   char  *UserDerField_Unit [DER_NOUT_MAX];

}; // struct InputPara_t



//-------------------------------------------------------------------------------------------------------
// Function    :  SyncHDF5File
// Description :  Force NFS to synchronize (dump data to the disk before return)
//
// Note        :  1. Synchronization is accomplished by first opening the target file with the "appending" mode
//                   and then close it immediately
//                2. It's an inline function
//
// Parameter   :  FileName : Target file name
//-------------------------------------------------------------------------------------------------------
inline void SyncHDF5File( const char *FileName )
{

   FILE *File      = fopen( FileName, "ab" );
   int CloseStatus = fclose( File );

   if ( CloseStatus != 0 )
      Aux_Message( stderr, "WARNING : failed to close the file \"%s\" for synchronization !!\n", FileName );

} // FUNCTION : SyncHDF5File



#endif // #ifndef __HDF5_TYPEDEF_H__
