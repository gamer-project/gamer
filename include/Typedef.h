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
   TESTPROB_HYDRO_CLUSTER_MERGER_VS_FLASH      =    4,
   TESTPROB_HYDRO_AGORA_ISOLATED_GALAXY        =    5,
   TESTPROB_HYDRO_CAUSTIC                      =    6,
   TESTPROB_HYDRO_SPHERICAL_COLLAPSE           =    7,
   TESTPROB_HYDRO_KELVIN_HELMHOLTZ_INSTABILITY =    8,
   TESTPROB_HYDRO_RIEMANN                      =    9,

   TESTPROB_ELBDM_EXTPOT                       = 1000;


// program initialization options
typedef int OptInit_t;
const OptInit_t
   INIT_BY_FUNCTION = 1,
   INIT_BY_RESTART  = 2,
   INIT_BY_FILE     = 3;


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
   INT_QUAR     = 7;


// data reconstruction TVD limiters
typedef int LR_Limiter_t;
const LR_Limiter_t
   LR_LIMITER_NONE = 0,
   VANLEER         = 1,
   GMINMOD         = 2,
   ALBADA          = 3,
   VL_GMINMOD      = 4,
   EXTPRE          = 5;


// TVD limiters for the WAF scheme
typedef int WAF_Limiter_t;
const WAF_Limiter_t
   WAF_LIMITER_NONE = 0,
   WAF_SUPERBEE     = 1,
   WAF_VANLEER      = 2,
   WAF_ALBADA       = 3,
   WAF_MINBEE       = 4;


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
   OUTPUT_DIAG      = 7;


// options in "Prepare_PatchData"
typedef int PrepUnit_t;
const PrepUnit_t
   UNIT_PATCH      = 1,
   UNIT_PATCHGROUP = 2;

typedef int NSide_t;
const NSide_t
   NSIDE_00 = 0,
   NSIDE_06 = 6,
   NSIDE_26 = 26;


// use the load-balance alternative function in "Buf_GetBufferData" and "Flag_Real"
typedef int UseLBFunc_t;
const UseLBFunc_t
   USELB_NO  = 0,
   USELB_YES = 1;


// enable check or not
typedef int Check_t;
const Check_t
   CHECK_OFF = 0,
   CHECK_ON  = 1;


// target solver in "InvokeSolvers"
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
  ;


// target mode in "Buf_GetBufferData and LB_GetBufferData"
typedef int GetBufMode_t;
const GetBufMode_t
#ifdef GRAVITY
   DATA_GENERAL      = 1,
   DATA_AFTER_FIXUP  = 2,
   DATA_AFTER_REFINE = 3,
   DATA_RESTRICT     = 4,
   COARSE_FINE_FLUX  = 5,
   POT_FOR_POISSON   = 6,
   POT_AFTER_REFINE  = 7;
#else
   DATA_GENERAL      = 1,
   DATA_AFTER_FIXUP  = 2,
   DATA_AFTER_REFINE = 3,
   DATA_RESTRICT     = 4,
   COARSE_FINE_FLUX  = 5;
#endif // #ifdef GRAVITY ... else ...


// fluid boundary conditions
typedef int OptFluBC_t;
const OptFluBC_t
   BC_FLU_NONE       = 0,
   BC_FLU_PERIODIC   = 1,
   BC_FLU_OUTFLOW    = 2,
   BC_FLU_REFLECTING = 3,
   BC_FLU_USER       = 4;


// the gravity boundary conditions
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
   PAR_INTERP_NONE    = 0,
   PAR_INTERP_NGP     = 1,
   PAR_INTERP_CIC     = 2,
   PAR_INTERP_TSC     = 3;

typedef int ParInteg_t;
const ParInteg_t
   PAR_INTEG_NONE    = 0,
   PAR_INTEG_EULER   = 1,
   PAR_INTEG_KDK     = 2;

typedef int ParUpStep_t;
const ParUpStep_t
   PAR_UPSTEP_PRED     = 1,
   PAR_UPSTEP_CORR     = 2,
   PAR_UPSTEP_ACC_ONLY = 3;

typedef int ParSync_t;
const ParSync_t
   PAR_SYNC_NONE    = 0,
   PAR_SYNC_TEMP    = 1,
   PAR_SYNC_FORCE   = 2;

typedef int ParOutputDens_t;
const ParOutputDens_t
   PAR_OUTPUT_DENS_NONE     = 0,
   PAR_OUTPUT_DENS_PAR_ONLY = 1,
   PAR_OUTPUT_DENS_TOTAL    = 2;
#endif // #ifdef PARTICLE


// the gravity types (this type needs to be defined for the Fluid solver even when GRAVITY is off)
typedef int OptGravityType_t;
const OptGravityType_t
   GRAVITY_NONE     = 0,
   GRAVITY_SELF     = 1,
   GRAVITY_EXTERNAL = 2,
   GRAVITY_BOTH     = 3;


// forms of the Lohner's error estimator
typedef int OptLohnerForm_t;
const OptLohnerForm_t
   LOHNER_FLASH1    = 1,
   LOHNER_FLASH2    = 2,
   LOHNER_FORM_INV1 = 3,
   LOHNER_FORM_INV2 = 4;


// OPT__1ST_FLUX_CORR and OPT__1ST_FLUX_CORR_SCHEME options
#if ( MODEL == HYDRO || MODEL == MHD )
typedef int Opt1stFluxCorr_t;
const Opt1stFluxCorr_t
   FIRST_FLUX_CORR_NONE    = 0,
   FIRST_FLUX_CORR_3D      = 1,
   FIRST_FLUX_CORR_3D1D    = 2;

typedef int OptRSolver1st_t;
const OptRSolver1st_t
   RSOLVER_1ST_NONE    = 0,
   RSOLVER_1ST_ROE     = 1,
   RSOLVER_1ST_HLLC    = 2,
   RSOLVER_1ST_HLLE    = 3;
#endif // #if ( MODEL == HYDRO || MODEL == MHD )


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


// Grackle
#ifdef SUPPORT_GRACKLE
// original Grackle or the reduced CPU/GPU implementation in GAMER
typedef int GrackleMode_t;
const GrackleMode_t
   GRACKLE_MODE_NONE  = 0,
   GRACKLE_MODE_ORI   = 1,
   GRACKLE_MODE_GAMER = 2;

// primordial chemistry
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



#endif  // #ifndef __TYPEDEF_H__
