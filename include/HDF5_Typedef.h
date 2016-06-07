#ifndef __HDF5_TYPEDEF_H__
#define __HDF5_TYPEDEF_H__


/*===========================================================================
Data structures defined here are mainly used for creating the compound
datatypes in the HDF5 format
===========================================================================*/

#include "hdf5.h"

#ifdef FLOAT8
#  define H5T_GAMER_REAL H5T_NATIVE_DOUBLE
#else
#  define H5T_GAMER_REAL H5T_NATIVE_FLOAT
#endif

#ifdef GAMER_DEBUG
#  define DEBUG_HDF5
#endif



//-------------------------------------------------------------------------------------------------------
// Structure   :  KeyInfo_t
// Description :  Data structure for outputting the key simulation information
//
// Note        :  1. Some variables have different names from their internal names in the code
//                   --> Make the output file more readable
//-------------------------------------------------------------------------------------------------------
struct KeyInfo_t
{

   int    FormatVersion;
   int    Model;
   int    Float8;
   int    Gravity;
   int    Particle;
   int    NLevel;
   int    PatchSize;
   int    DumpID;
#  ifdef GRAVITY
   int    OutputPot;
#  endif
   int    NX0     [3];
   int    BoxScale[3];
   int    NPatch   [NLEVEL];
   int    CellScale[NLEVEL];        // amr->scale[lv]

   long   Step;
   long   AdvanceCounter[NLEVEL];

   double BoxSize [3];
   double Time    [NLEVEL];
   double CellSize[NLEVEL];         // amr->dh[lv]
#  ifdef GRAVITY
   double AveDens;                  // AveDensity
#  endif

   char  *CodeVersion;
   char  *DumpWallTime;

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
   int IndividualDt;
   int Comoving;
   int Particle;

   int UseGPU;
   int GAMER_Optimization;
   int GAMER_Debug;
   int Timing;
   int TimingSolver;
   int Intel;
   int Float8;
   int Serial;
   int LoadBalance;
   int OverlapMPI;
   int OpenMP;
   int GPU_Arch;
   int Laohu;
   int SupportHDF5;

   int NLevel;
   int MaxPatch;

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
   int NPassive;

#  elif ( MODEL == MHD )
#  warning : WAIT MHD !!!

#  elif ( MODEL == ELBDM )
   int ConserveMass;
   int Laplacian4th;
   int SelfInteraction4;

#  else
#  error : unsupported MODEL !!
#  endif // MODEL

#  ifdef PARTICLE
   int StoreParAcc;
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

   int    NComp;
   int    PatchSize;
   int    Flu_NIn;
   int    Flu_NOut;
   int    NFlux;
   int    Flu_GhostSize;
   int    Flu_Nxt;
   int    Debug_HDF5;
   int    SibOffsetNonperiodic;

#  ifdef LOAD_BALANCE
   int    SonOffsetLB;
#  endif

   double TinyValue;


#  ifdef GRAVITY
   int    Gra_NIn;
   int    Pot_GhostSize;
   int    Gra_GhostSize;
   int    Rho_GhostSize;
   int    Pot_Nxt;
   int    Gra_Nxt;
   int    Rho_Nxt;

#  ifdef UNSPLIT_GRAVITY
   int    USG_GhostSize;
   int    USG_NxtF;
   int    USG_NxtG;
#  endif

   int    Gra_BlockSize_z;
   int    ExtPotNAuxMax;
   int    ExtAccNAuxMax;

#  if   ( POT_SCHEME == SOR )
   int    Pot_BlockSize_z;
   int    UsePSolver_10to14;
#  elif ( POT_SCHEME == MG  )
   int    Pot_BlockSize_x;
#  endif

#  endif // #ifdef GRAVITY


#  ifdef PARTICLE
   int    NPar_Var;
   int    NPar_Passive;
   int    PM_GhostSize;
   int    PM_Nxt;

   int    Debug_Particle;

   double ParList_GrowthFactor;
   double ParList_ReduceFactor;
#  endif


#  if   ( MODEL == HYDRO )
   int    Flu_BlockSize_x;
   int    Flu_BlockSize_y;
   int    CheckNegativeInFluid;
   int    CharReconstruction;
   int    CheckIntermediate;
   int    HLL_NoRefState;
   int    HLL_IncludeAllWaves;
   int    WAF_Dissipate;
   int    PositiveDensInFixUp;

#  ifdef N_FC_VAR
   int    N_FC_Var;
#  endif

#  ifdef N_SLOPE_PPM
   int    N_Slope_PPM;
#  endif

#  ifdef MIN_PRES_DENS
   double Min_Pres_Dens;
#  endif

#  ifdef MIN_PRES
   double Min_Pres;
#  endif

#  ifdef MAX_ERROR
   double MaxError;
#  endif

#  elif ( MODEL == MHD )
   int    Flu_BlockSize_x;
   int    Flu_BlockSize_y;
#  warning : WAIT MHD !!!


#  elif  ( MODEL == ELBDM )
   int    Flu_BlockSize_x;
   int    Flu_BlockSize_y;

#  else
#  error : ERROR : unsupported MODEL !!
#  endif // MODEL

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

// boundary condition
   int    Opt__BC_Flu[6];
#  ifdef GRAVITY
   int    Opt__BC_Pot;
   double GFunc_Coeff0;
#  endif

// particle
#  ifdef PARTICLE
   long   Par_NPar;
   int    Par_Init;
   int    Par_Interp;
   int    Par_Integ;
   int    Par_SyncDump;
   int    Par_ImproveAcc;
   int    Par_PredictPos;
   double Par_RemoveCell;
#  endif

// cosmology
#  ifdef COMOVING
   double A_Init;
   double OmegaM0;
#  endif

// time-step determination
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
#  ifdef COMOVING
   double Dt__MaxDeltaA;
#  endif
   int    Opt__AdaptiveDt;
   int    Opt__RecordDt;
   int    Opt__DtUser;
   
// domain refinement
   int    RegridCount;
   int    FlagBufferSize;
   int    MaxLevel;
   int    Opt__Flag_Rho;
   int    Opt__Flag_RhoGradient;
#  if ( MODEL == HYDRO ) 
   int    Opt__Flag_PresGradient;
#  endif
#  if ( MODEL == ELBDM ) 
   int    Opt__Flag_EngyDensity;
#  endif
   int    Opt__Flag_LohnerDens;
#  if ( MODEL == HYDRO ) 
   int    Opt__Flag_LohnerEngy;
   int    Opt__Flag_LohnerPres;
#  endif
   int    Opt__Flag_LohnerForm;
   int    Opt__Flag_User;
   int    Opt__Flag_Region;
   int    Opt__PatchCount;
#  ifdef PARTICLE
   int    Opt__ParLevel;
#  endif

// load balance
#  ifdef LOAD_BALANCE
   double LB_Input__WLI_Max;
#  endif

// fluid solvers in HYDRO
#  if ( MODEL == HYDRO )
   double Gamma;   
   double MinMod_Coeff;
   double EP_Coeff;
   int    Opt__LR_Limiter;
   int    Opt__WAF_Limiter;
   int    Opt__CorrUnphyScheme;
#  endif

// ELBDM solvers
#  if ( MODEL == ELBDM )
   double ELBDM_Mass;
   double ELBDM_PlanckConst;
#  ifdef QUARTIC_SELF_INTERACTION
   double ELBDM_Lambda;
#  endif
   double ELBDM_Taylor3_Coeff;
   int    ELBDM_Taylor3_Auto;
#  endif

// fluid solvers in both HYDRO/MHD/ELBDM
   int    Flu_GPU_NPGroup;    
   int    GPU_NStream;
   int    Opt__FixUp_Flux;
   int    Opt__FixUp_Restrict;
   int    Opt__OverlapMPI;
   int    Opt__ResetFluid;
   int    Opt__CorrUnphy;

// self-gravity
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
   int    Opt__GravityType;
   int    Opt__ExternalPot;
#  endif

// initialization
   int    Opt__Init;
   int    Opt__RestartHeader;   
   int    Opt__UM_Start_Level;
   int    Opt__UM_Start_NVar;
   int    Opt__UM_Start_Downgrade;
   int    Opt__UM_Start_Refine;
   int    Opt__UM_Factor_5over3;
   int    Opt__InitRestrict;
   int    Opt__GPUID_Select;
   int    Init_Subsampling_NCell;

// interpolation schemes
   int    Opt__Int_Time;   
#  if ( MODEL == ELBDM ) 
   int    Opt__Int_Phase;
#  endif
   int    Opt__Flu_IntScheme;
#  ifdef GRAVITY
   int    Opt__Pot_IntScheme;   
   int    Opt__Rho_IntScheme;
   int    Opt__Gra_IntScheme;
#  endif
   int    Opt__RefFlu_IntScheme;
#  ifdef GRAVITY
   int    Opt__RefPot_IntScheme;
#  endif
   double IntMonoCoeff;

// data dump
   int    Opt__Output_Total;
   int    Opt__Output_Part;
   int    Opt__Output_TestError; 
#  ifdef PARTICLE
   int    Opt__Output_Particle;
#  endif
   int    Opt__Output_BasePS; 
   int    Opt__Output_Base;
#  ifdef GRAVITY
   int    Opt__Output_Pot;
#  endif
   int    Opt__Output_Mode; 
   int    Opt__Output_Step;
   double Opt__Output_Dt;
   double Output_PartX;
   double Output_PartY;
   double Output_PartZ;
   int    InitDumpID;

// miscellaneous
   int    Opt__Verbose;
   int    Opt__TimingBalance;
   int    Opt__TimingMPI;   
   int    Opt__RecordMemory;
   int    Opt__RecordPerformance;
   int    Opt__ManualControl;
   int    Opt__RecordUser;

// simulation checks
   int    Opt__Ck_Refine;
   int    Opt__Ck_ProperNesting;
   int    Opt__Ck_Conservation;
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

// flag tables
   double FlagTable_Rho         [NLEVEL-1]; 
   double FlagTable_RhoGradient [NLEVEL-1]; 
   double FlagTable_Lohner      [NLEVEL-1][3];
   double FlagTable_User        [NLEVEL-1];
#  if   ( MODEL == HYDRO )
   double FlagTable_PresGradient[NLEVEL-1];
#  elif ( MODEL == ELBDM )
   double FlagTable_EngyDensity [NLEVEL-1][2];
#  endif
   
}; // struct InputPara_t



/*
//-------------------------------------------------------------------------------------------------------
// Structure   :  PatchInfo_t 
// Description :  Data structure for outputting the patch information
//
// Note        :  1. Father, Son and Sibling store the "global identification (GID)" instead of the "local
//                   patch identification (PID)"
//                   --> Patches at the same level at all MPI ranks will have different GID but can have
//                       the same PID
//                2. LB_Idx will be stored in all parallelization modes (SERIAL, NO_LOAD_BALANCE, LOAD_BALANCE)
//-------------------------------------------------------------------------------------------------------
struct PatchInfo_t
{

   int    Lv;
   int    GID;
   int    Father;
   int    Son;
   int    Sibling[26];
   int    Corner[3];

   ulong  PaddedCr1D;
   long   LBIdx;

   double EdgeL[3];
   double EdgeR[3];

#  ifdef PARTICLE
// add variables related to particles
#  endif

}; // struct PatchInfo_t



//-------------------------------------------------------------------------------------------------------
// Structure   :  PatchData_t 
// Description :  Data structure for outputting the patch data
//
// Note        :  1. "Pote" array MUST be the LAST member
//-------------------------------------------------------------------------------------------------------
struct PatchData_t 
{

#  if   ( MODEL == HYDRO )
   real Dens[PS1][PS1][PS1];
   real MomX[PS1][PS1][PS1];
   real MomY[PS1][PS1][PS1];
   real MomZ[PS1][PS1][PS1];
   real Engy[PS1][PS1][PS1];

#  elif ( MODEL == MHD )
#  warning : WAIT MHD !!!

#  elif ( MODEL == ELBDM )
   real Dens[PS1][PS1][PS1];
   real Real[PS1][PS1][PS1];
   real Imag[PS1][PS1][PS1];

#  else
#  error : unsupported MODEL !!
#  endif // MODEL

// add passive variables

// Pote MUST be the LAST member of the PatchData_t structure (since it's output only if OPT__OUTPUT_POT == true)
#  ifdef GRAVITY
   real Pote[PS1][PS1][PS1];
#  endif

}; // struct PatchData_t
*/



//-------------------------------------------------------------------------------------------------------
// Function    :  SyncHDF5File 
// Description :  Force NFS to synchronize (dump data to the disk before return)
//
// Note        :  1. Synchronization is accomplished by first opening the target file with the "appending" mode
//                   and then close it immediately
//                3. It's an inline function
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
