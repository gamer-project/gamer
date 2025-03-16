#ifndef __TYPEDEF_H__
#define __TYPEDEF_H__



// *****************************************************************************
// ** This header defines the different data types in GAMER.                 **
// ** Please DO NOT modify the number assigned to any enumerator constants!! **
// **                                                                        **
// ** We no longer use "enum" to make loading runtime parameters easier.     **
// ** --> Please use "typedef int" instead.                                  **
// *****************************************************************************


#include "Macro.h"


// single/double precision
#ifdef FLOAT8
typedef double real;
#else
typedef float  real;
#endif

#ifdef FLOAT8_PAR
typedef double real_par;
#else
typedef float  real_par;
#endif

#ifdef INT8_PAR
typedef long long_par;
#else
typedef int  long_par;
#endif

#ifdef SUPPORT_GRACKLE
#include <grackle_float.h>
#if   defined GRACKLE_FLOAT_8
typedef double real_che;
#elif defined GRACKLE_FLOAT_4
typedef float  real_che;
#else
#error : ERROR : GRACKLE_FLOAT_8 and GRACKLE_FLOAT_4 are not defined in Grackle library !!
#endif
#endif // #ifdef SUPPORT_GRACKLE

#if ( GRAMFE_SCHEME == GRAMFE_FFT )
#ifdef GRAMFE_FFT_FLOAT8
typedef double gramfe_fft_float;
#else
typedef float  gramfe_fft_float;
#endif
#endif // #if ( GRAMFE_SCHEME == GRAMFE_FFT )

#ifdef GRAMFE_MATMUL_FLOAT8
typedef double gramfe_matmul_float;
#else
typedef float  gramfe_matmul_float;
#endif


// short names for unsigned type
typedef unsigned short     ushort;
typedef unsigned int       uint;
typedef unsigned long int  ulong;


// test problem IDs
typedef int TestProbID_t;
const TestProbID_t
   TESTPROB_NONE                               =    0,

   TESTPROB_HYDRO_BLAST_WAVE                   =    1,
   TESTPROB_HYDRO_ACOUSTIC_WAVE                =    2,
   TESTPROB_HYDRO_BONDI                        =    3,
   TESTPROB_HYDRO_CLUSTER_MERGER               =    4,
   TESTPROB_HYDRO_AGORA_ISOLATED_GALAXY        =    5,
   TESTPROB_HYDRO_CAUSTIC                      =    6,
   TESTPROB_HYDRO_SPHERICAL_COLLAPSE           =    7,
   TESTPROB_HYDRO_KELVIN_HELMHOLTZ_INSTABILITY =    8,
   TESTPROB_HYDRO_RIEMANN                      =    9,
   TESTPROB_HYDRO_JET                          =   10,
   TESTPROB_HYDRO_PLUMMER                      =   11,
   TESTPROB_HYDRO_GRAVITY                      =   12,
   TESTPROB_HYDRO_MHD_ABC                      =   13,
   TESTPROB_HYDRO_MHD_ORSZAG_TANG_VORTEX       =   14,
   TESTPROB_HYDRO_MHD_LINEAR_WAVE              =   15,
   TESTPROB_HYDRO_JEANS_INSTABILITY            =   16,
   TESTPROB_HYDRO_PARTICLE_EQUILIBRIUM_IC      =   17,
   TESTPROB_HYDRO_PARTICLE_TEST                =   18,
   TESTPROB_HYDRO_ENERGY_POWER_SPECTRUM        =   19,
   TESTPROB_HYDRO_CR_SOUNDWAVE                 =   20,
   TESTPROB_HYDRO_CR_SHOCKTUBE                 =   21,
   TESTPROB_HYDRO_CR_DIFFUSION                 =   23,
   TESTPROB_HYDRO_BARRED_POT                   =   51,
   TESTPROB_HYDRO_JET_ICM_WALL                 =   52,
   TESTPROB_HYDRO_CDM_LSS                      =  100,
   TESTPROB_HYDRO_ZELDOVICH                    =  101,
   TESTPROB_HYDRO_GRACKLE_COMOVING             =  102,
   TESTPROB_ELBDM_EXTPOT                       = 1000,
   TESTPROB_ELBDM_JEANS_INSTABILITY_COMOVING   = 1001,
   TESTPROB_ELBDM_JEANS_INSTABILITY_PHYSICAL   = 1002,
   TESTPROB_ELBDM_SOLITON                      = 1003,
   TESTPROB_ELBDM_SELF_SIMILAR_HALO            = 1004,
   TESTPROB_ELBDM_VORTEX_PAIR_ROTATING         = 1005,
   TESTPROB_ELBDM_VORTEX_PAIR_LINEAR           = 1006,
   TESTPROB_ELBDM_ISOLATED_HALO                = 1007,
   TESTPROB_ELBDM_GAUSSIAN_WAVE_PACKET         = 1008,
   TESTPROB_ELBDM_LSS                          = 1009,
   TESTPROB_ELBDM_PLANE_WAVE                   = 1010,
   TESTPROB_ELBDM_PERTURBATION                 = 1011,
   TESTPROB_ELBDM_HALO_MERGER                  = 1012,
   TESTPROB_ELBDM_DISK_HEATING                 = 1013;

// program initialization options
typedef int OptInit_t;
const OptInit_t
   INIT_BY_FUNCTION = 1,
   INIT_BY_RESTART  = 2,
   INIT_BY_FILE     = 3;


// program initialization options for the magnetic field by vector potential
typedef int OptInitMagByVecPot_t;
const OptInitMagByVecPot_t
   INIT_MAG_BYVECPOT_NONE = 0,
   INIT_MAG_BYVECPOT_FILE = 1,
   INIT_MAG_BYVECPOT_FUNC = 2;


// data format for OPT__INIT=INIT_BY_FILE
typedef int UM_IC_Format_t;
const UM_IC_Format_t
   UM_IC_FORMAT_NONE = 0,
   UM_IC_FORMAT_VZYX = 1,
   UM_IC_FORMAT_ZYXV = 2;


// data format for PAR_INIT=PAR_INIT_BY_FILE
typedef int ParICFormat_t;
const ParICFormat_t
   PAR_IC_FORMAT_NONE   = 0,
   PAR_IC_FORMAT_ATT_ID = 1,
   PAR_IC_FORMAT_ID_ATT = 2;


// FFTW startup options
typedef int FFTWStartup_t;
const FFTWStartup_t
   FFTW_STARTUP_DEFAULT  = -1,
   FFTW_STARTUP_ESTIMATE = 0,
   FFTW_STARTUP_MEASURE  = 1,
   FFTW_STARTUP_PATIENT  = 2;


// program restart options
typedef int OptRestartH_t;
const OptRestartH_t
   RESTART_HEADER_SKIP  = 0,
   RESTART_HEADER_CHECK = 1;


// interpolation schemes
typedef int IntScheme_t;
const IntScheme_t
   INT_DEFAULT  = -1,
   INT_NONE     = 0,
   INT_MINMOD3D = 1,
   INT_MINMOD1D = 2,
   INT_VANLEER  = 3,
   INT_CQUAD    = 4,
   INT_QUAD     = 5,
   INT_CQUAR    = 6,
   INT_QUAR     = 7,
   INT_SPECTRAL = 8;


// data reconstruction TVD limiters
typedef int LR_Limiter_t;
const LR_Limiter_t
   LR_LIMITER_DEFAULT    = -1,
   LR_LIMITER_NONE       = 0,
   LR_LIMITER_VANLEER    = 1,
   LR_LIMITER_GMINMOD    = 2,
   LR_LIMITER_ALBADA     = 3,
   LR_LIMITER_VL_GMINMOD = 4,
   LR_LIMITER_EXTPRE     = 5,
   LR_LIMITER_CENTRAL    = 6,
   LR_LIMITER_ATHENA     = 7;


// data output formats
typedef int OptOutputFormat_t;
const OptOutputFormat_t
   OUTPUT_TOTAL_NONE     = 0,
   OUTPUT_FORMAT_HDF5    = 1,
   OUTPUT_FORMAT_CBINARY = 2;


// data output criteria
typedef int OptOutputMode_t;
const OptOutputMode_t
   OUTPUT_CONST_STEP = 1,
   OUTPUT_CONST_DT   = 2,
   OUTPUT_USE_TABLE  = 3;


// OPT__OUTPUT_PART options
typedef int OptOutputPart_t;
const OptOutputPart_t
   OUTPUT_PART_NONE = 0,
   OUTPUT_XY        = 1,
   OUTPUT_YZ        = 2,
   OUTPUT_XZ        = 3,
   OUTPUT_X         = 4,
   OUTPUT_Y         = 5,
   OUTPUT_Z         = 6,
   OUTPUT_DIAG      = 7,
   OUTPUT_BOX       = 8;


// OPT_OUTPUT_PAR_MODE options
typedef int OptOutputParMode_t;
const OptOutputParMode_t
   OUTPUT_PAR_NONE = 0,
   OUTPUT_PAR_TEXT = 1,
   OUTPUT_PAR_CBIN = 2;


// options in Prepare_PatchData()
typedef int PrepUnit_t;
const PrepUnit_t
   UNIT_PATCH      = 1,
   UNIT_PATCHGROUP = 2;

typedef int NSide_t;
const NSide_t
   NSIDE_00 = 0,
   NSIDE_06 = 6,
   NSIDE_26 = 26;


// use the load-balance alternative functions
typedef int UseLBFunc_t;
const UseLBFunc_t
   USELB_NO  = 0,
   USELB_YES = 1;


// enable check or not
typedef int Check_t;
const Check_t
   CHECK_OFF = 0,
   CHECK_ON  = 1;


// modes of Hydro_IsUnphysical()
typedef int IsUnphyMode_t;
const IsUnphyMode_t
   UNPHY_MODE_SING         = 0,  // check single field
   UNPHY_MODE_CONS         = 1,  // check conserved variables, including passive scalars
   UNPHY_MODE_PRIM         = 2,  // check primitive variables, including passive scalars
   UNPHY_MODE_PASSIVE_ONLY = 3;  // only check passive scalars


// verbosity levels of Hydro_IsUnphysical()
typedef int IsUnphVerb_t;
const IsUnphVerb_t
   UNPHY_SILENCE = 0,   // print nothing
   UNPHY_VERBOSE = 1;   // print out unphysical values


// whether the interpolated fields include all conserved variables in hydrodynamics
typedef bool AllCons_t;
const AllCons_t
   ALL_CONS_NO  = false,
   ALL_CONS_YES = true;


// locally reduce the monotonic coefficient or not
typedef bool ReduceOrFixMonoCoeff_t;
const ReduceOrFixMonoCoeff_t
   INT_FIX_MONO_COEFF    = false,   // fix the coefficient
   INT_REDUCE_MONO_COEFF = true;    // locally reduce the coefficient


// whether switch from conserved to primitive variables when interpolation fails
typedef bool IntPrim_t;
const IntPrim_t
   INT_PRIM_NO  = false,
   INT_PRIM_YES = true;


// target solver in InvokeSolver()
// --> must start from 0 because of the current TIMING_SOLVER implementation
// --> when adding new solvers, please modify the NSOLVER constant accordingly
const int NSOLVER = 7;

typedef int Solver_t;
const Solver_t
   FLUID_SOLVER               = 0
#ifdef GRAVITY
  ,POISSON_SOLVER             = 1
  ,GRAVITY_SOLVER             = 2
  ,POISSON_AND_GRAVITY_SOLVER = 3
#endif
#ifdef SUPPORT_GRACKLE
  ,GRACKLE_SOLVER             = 4
#endif
  ,DT_FLU_SOLVER              = 5
#ifdef GRAVITY
  ,DT_GRA_SOLVER              = 6
#endif
  ,SRC_SOLVER                 = 7
  ;


// target mode in Buf_GetBufferData() and LB_GetBufferData()
typedef int GetBufMode_t;
const GetBufMode_t
   DATA_GENERAL         = 1
  ,DATA_AFTER_FIXUP     = 2
  ,DATA_AFTER_REFINE    = 3
  ,DATA_RESTRICT        = 4
  ,COARSE_FINE_FLUX     = 5
#ifdef GRAVITY
  ,POT_FOR_POISSON      = 6
  ,POT_AFTER_REFINE     = 7
#endif
#ifdef MHD
  ,COARSE_FINE_ELECTRIC = 8
#endif
  ;


// fluid boundary conditions
typedef int OptFluBC_t;
const OptFluBC_t
   BC_FLU_NONE       = 0,
   BC_FLU_PERIODIC   = 1,
   BC_FLU_OUTFLOW    = 2,
   BC_FLU_REFLECTING = 3,
   BC_FLU_USER       = 4,
   BC_FLU_DIODE      = 5;


// gravity boundary conditions
typedef int OptPotBC_t;
const OptPotBC_t
#ifdef GRAVITY
   BC_POT_NONE     = 0,
   BC_POT_PERIODIC = 1,
   BC_POT_ISOLATED = 2;
#else
   BC_POT_NONE     = 0;
#endif


// particle schemes
#ifdef PARTICLE
typedef int ParInit_t;
const ParInit_t
   PAR_INIT_NONE        = 0,
   PAR_INIT_BY_FUNCTION = 1,
   PAR_INIT_BY_RESTART  = 2,
   PAR_INIT_BY_FILE     = 3;

typedef int ParInterp_t;
const ParInterp_t
   PAR_INTERP_NONE = 0,
   PAR_INTERP_NGP  = 1,
   PAR_INTERP_CIC  = 2,
   PAR_INTERP_TSC  = 3;

typedef int ParInteg_t;
const ParInteg_t
   PAR_INTEG_NONE  = 0,
   PAR_INTEG_EULER = 1,
   PAR_INTEG_KDK   = 2;

typedef int TracerInteg_t;
const TracerInteg_t
   TRACER_INTEG_NONE  = 0,
   TRACER_INTEG_EULER = 1,
   TRACER_INTEG_RK2   = 2;

typedef int ParUpStep_t;
const ParUpStep_t
   PAR_UPSTEP_PRED     = 1,
   PAR_UPSTEP_CORR     = 2,
   PAR_UPSTEP_ACC_ONLY = 3;

typedef int ParSync_t;
const ParSync_t
   PAR_SYNC_NONE  = 0,
   PAR_SYNC_TEMP  = 1,
   PAR_SYNC_FORCE = 2;

typedef int ParOutputDens_t;
const ParOutputDens_t
   PAR_OUTPUT_DENS_NONE     = 0,
   PAR_OUTPUT_DENS_PAR_ONLY = 1,
   PAR_OUTPUT_DENS_TOTAL    = 2;

typedef int ParPass2Son_t;
const ParPass2Son_t
   PAR_PASS2SON_GENERAL = 1,
   PAR_PASS2SON_EVOLVE  = 2;
#endif // #ifdef PARTICLE


// external acceleration (must be defined for the fluid solver even when GRAVITY is off)
typedef int OptExtAcc_t;
const OptExtAcc_t
   EXT_ACC_NONE  = 0,
   EXT_ACC_FUNC  = 1,
   EXT_ACC_TABLE = 2;


// external potential (must be defined for the fluid solver even when GRAVITY is off)
typedef int OptExtPot_t;
const OptExtPot_t
   EXT_POT_NONE  = 0,
   EXT_POT_FUNC  = 1,
   EXT_POT_TABLE = 2;


// different usages of external potential when computing total potential on level Lv
// --> ADD     : add external potential on Lv
//     SUB     : subtract external potential for preparing self-gravity potential on Lv-1
//     SUB_TINT: like SUB but for temporal interpolation
typedef int ExtPotUsage_t;
const ExtPotUsage_t
   EXT_POT_USAGE_ADD      = 0,
   EXT_POT_USAGE_SUB      = 1,
   EXT_POT_USAGE_SUB_TINT = 2;


// forms of the Lohner's error estimator
typedef int OptLohnerForm_t;
const OptLohnerForm_t
   LOHNER_FLASH1    = 1,
   LOHNER_FLASH2    = 2,
   LOHNER_FORM_INV1 = 3,
   LOHNER_FORM_INV2 = 4;


// OPT__1ST_FLUX_CORR and OPT__1ST_FLUX_CORR_SCHEME options
#if ( MODEL == HYDRO )
typedef int Opt1stFluxCorr_t;
const Opt1stFluxCorr_t
   FIRST_FLUX_CORR_NONE = 0,
   FIRST_FLUX_CORR_3D   = 1,
   FIRST_FLUX_CORR_3D1D = 2;

typedef int OptRSolver1st_t;
const OptRSolver1st_t
   RSOLVER_1ST_DEFAULT = -1,
   RSOLVER_1ST_NONE    = 0,
   RSOLVER_1ST_ROE     = 1,
   RSOLVER_1ST_HLLC    = 2,
   RSOLVER_1ST_HLLE    = 3,
   RSOLVER_1ST_HLLD    = 4;
#endif // #if ( MODEL == HYDRO )


// OPT__CORR_AFTER_ALL_SYNC options
typedef int OptCorrAfterSync_t;
const OptCorrAfterSync_t
   CORR_AFTER_SYNC_DEFAULT     = -1,
   CORR_AFTER_SYNC_NONE        = 0,
   CORR_AFTER_SYNC_EVERY_STEP  = 1,
   CORR_AFTER_SYNC_BEFORE_DUMP = 2;


// OPT__DT_LEVEL options
typedef int OptTimeStepLevel_t;
const OptTimeStepLevel_t
   DT_LEVEL_SHARED    = 1,
   DT_LEVEL_DIFF_BY_2 = 2,
   DT_LEVEL_FLEXIBLE  = 3;


// AddField() options
typedef int FixUpFlux_t;
const FixUpFlux_t
   FIXUP_FLUX_NO  = 0,
   FIXUP_FLUX_YES = 1;

typedef int FixUpRestrict_t;
const FixUpRestrict_t
   FIXUP_REST_NO  = 0,
   FIXUP_REST_YES = 1;

typedef int NormPassive_t;
const NormPassive_t
   NORMALIZE_NO  = 0,
   NORMALIZE_YES = 1;

typedef int IntFracPassive_t;
const IntFracPassive_t
   INTERP_FRAC_NO  = 0,
   INTERP_FRAC_YES = 1;


// field types
typedef int FieldIdx_t;


// Grackle
#ifdef SUPPORT_GRACKLE
// map to the "primordial_chemistry" option of Grackle
typedef int GracklePriChe_t;
const GracklePriChe_t
   GRACKLE_PRI_CHE_CLOUDY = 0,
   GRACKLE_PRI_CHE_NSPE6  = 1,
   GRACKLE_PRI_CHE_NSPE9  = 2,
   GRACKLE_PRI_CHE_NSPE12 = 3;
#endif


// star formation
#ifdef STAR_FORMATION
// schemes of creating new star particles
typedef int SF_CreateStarScheme_t;
const SF_CreateStarScheme_t
   SF_CREATE_STAR_SCHEME_NONE  = 0,
   SF_CREATE_STAR_SCHEME_AGORA = 1;
#endif


// ELBDM_REMOVE_MOTION_CM options
#if ( MODEL == ELBDM )
typedef int ELBDMRemoveMotionCM_t;
const ELBDMRemoveMotionCM_t
   ELBDM_REMOVE_MOTION_CM_NONE       = 0,
   ELBDM_REMOVE_MOTION_CM_INIT       = 1,
   ELBDM_REMOVE_MOTION_CM_EVERY_STEP = 2;
#endif


// options in Aux_ComputeProfile() and Aux_FindExtrema()
typedef int PatchType_t;
const PatchType_t
   PATCH_LEAF                 = 0,
   PATCH_NONLEAF              = 1,
   PATCH_BOTH                 = 2,
   PATCH_LEAF_PLUS_MAXNONLEAF = 3;


// options in Aux_FindExtrema()
typedef int ExtremaMode_t;
const ExtremaMode_t
   EXTREMA_MIN = 1,
   EXTREMA_MAX = 2;


// options in LoadInputTestProb()
typedef int LoadParaMode_t;
const LoadParaMode_t
   LOAD_READPARA    = 1,
   LOAD_HDF5_OUTPUT = 2;


// function pointers
typedef real (*EoS_GUESS_t)    ( const real Con[], real* const Constant, const double AuxArray_Flt[],
                                 const int AuxArray_Int[], const real *const Table[EOS_NTABLE_MAX] );
typedef void (*EoS_H2TEM_t)    ( const real HTilde, real* const Temp, real* const DiffTemp,
                                 const real Passive[], const double AuxArray_Flt[],
                                 const int AuxArray_Int[], const real *const Table[EOS_NTABLE_MAX] );
typedef real (*EoS_TEM2H_t)    ( const real Temp, const real Passive[], const double AuxArray_Flt[],
                                 const int AuxArray_Int[], const real *const Table[EOS_NTABLE_MAX] );
typedef real (*EoS_DE2P_t)     ( const real Dens, const real Eint, const real Passive[],
                                 const double AuxArray_Flt[], const int AuxArray_Int[],
                                 const real *const Table[EOS_NTABLE_MAX] );
typedef real (*EoS_DP2E_t)     ( const real Dens, const real Pres, const real Passive[],
                                 const double AuxArray_Flt[], const int AuxArray_Int[],
                                 const real *const Table[EOS_NTABLE_MAX] );
typedef real (*EoS_DP2C_t)     ( const real Dens, const real Pres, const real Passive[],
                                 const double AuxArray_Flt[], const int AuxArray_Int[],
                                 const real *const Table[EOS_NTABLE_MAX] );
typedef real (*EoS_DE2T_t)     ( const real Dens, const real Eint, const real Passive[],
                                 const double AuxArray_Flt[], const int AuxArray_Int[],
                                 const real *const Table[EOS_NTABLE_MAX] );
typedef real (*EoS_DT2P_t)     ( const real Dens, const real Temp, const real Passive[],
                                 const double AuxArray_Flt[], const int AuxArray_Int[],
                                 const real *const Table[EOS_NTABLE_MAX] );
typedef real (*EoS_DE2S_t)     ( const real Dens, const real Eint, const real Passive[],
                                 const double AuxArray_Flt[], const int AuxArray_Int[],
                                 const real *const Table[EOS_NTABLE_MAX] );
typedef void (*EoS_GENE_t)     ( const int Mode, real Out[], const real In_Flt[], const int In_Int[],
                                 const double AuxArray_Flt[], const int AuxArray_Int[],
                                 const real *const Table[EOS_NTABLE_MAX] );
#ifdef COSMIC_RAY
typedef real (*EoS_CRE2CRP_t)  ( const real E_CR,
                                 const double AuxArray_Flt[], const int AuxArray_Int[],
                                 const real *const Table[EOS_NTABLE_MAX] );
#endif
typedef void (*ExtAcc_t)       ( real Acc[], const double x, const double y, const double z, const double Time,
                                 const double UserArray[] );
typedef real (*ExtPot_t)       ( const double x, const double y, const double z, const double Time,
                                 const double UserArray_Flt[], const int UserArray_Int[],
                                 const ExtPotUsage_t Usage, const real PotTable[], void **GenePtr );
typedef void (*IntSchemeFunc_t)( real CData[], const int CSize[3], const int CStart[3], const int CRange[3],
                                 real FData[], const int FSize[3], const int FStart[3], const int NComp,
                                 const bool UnwrapPhase, const bool Monotonic[], const real MonoCoeff, const bool OppSign0thOrder );



#endif  // #ifndef __TYPEDEF_H__
