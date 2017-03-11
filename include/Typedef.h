#ifndef __TYPEDEF_H__
#define __TYPEDEF_H__



// ****************************************************************************
// ** This header defines the "typedef and enum" in GAMER.                   **
// ** Please DO NOT modify the number assigned to any enumerator constant !! **
// ****************************************************************************

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


// options of the program initialization
enum OptInit_t { INIT_STARTOVER=1, INIT_RESTART=2, INIT_UM=3 };


// options of the program restart
enum OptRestartH_t { RESTART_HEADER_SKIP=0, RESTART_HEADER_CHECK=1 };


// interpolation scheme
enum IntScheme_t { INT_DEFAULT=-1, INT_NONE=0, INT_MINMOD3D=1, INT_MINMOD1D=2, INT_VANLEER=3, INT_CQUAD=4, INT_QUAD=5,
                   INT_CQUAR=6, INT_QUAR=7 };


// data reconstruction TVD limiters
enum LR_Limiter_t  { LR_LIMITER_NONE=0, VANLEER=1, GMINMOD=2, ALBADA=3, VL_GMINMOD=4, EXTPRE=5 };
enum WAF_Limiter_t { WAF_LIMITER_NONE=0, WAF_SUPERBEE=1, WAF_VANLEER=2, WAF_ALBADA=3, WAF_MINBEE=4 };


// options of the data output format
enum OptOutputFormat_t{ OUTPUT_TOTAL_NONE=0, OUTPUT_FORMAT_HDF5=1, OUTPUT_FORMAT_CBINARY=2 };


// options of the data output criteria
enum OptOutputMode_t { OUTPUT_CONST_STEP=1, OUTPUT_CONST_DT=2, OUTPUT_USE_TABLE=3 };


// options of the output part
enum OptOutputPart_t { OUTPUT_PART_NONE=0, OUTPUT_XY=1, OUTPUT_YZ=2, OUTPUT_XZ=3, OUTPUT_X=4, OUTPUT_Y=5, OUTPUT_Z=6,
                       OUTPUT_DIAG=7 };


// options in "Prepare_PatchData"
enum PrepUnit_t { UNIT_PATCH=1, UNIT_PATCHGROUP=2 };
enum NSide_t    { NSIDE_00=0, NSIDE_06=6, NSIDE_26=26 };


// use the load-balance alternative function in "Buf_GetBufferData" and "Flag_Real"
enum UseLBFunc_t { USELB_NO=0, USELB_YES=1 };


// enable check or not
enum Check_t { CHECK_OFF=0, CHECK_ON=1 };


// targeted solver in "InvokeSolvers"
#ifdef GRAVITY
enum Solver_t { FLUID_SOLVER=0, POISSON_SOLVER=1, GRAVITY_SOLVER=2, POISSON_AND_GRAVITY_SOLVER=3 };
#else
enum Solver_t { FLUID_SOLVER=0 };
#endif


// targeted mode in "Buf_GetBufferData and LB_GetBufferData"
#ifdef GRAVITY
enum GetBufMode_t { DATA_GENERAL=1, DATA_AFTER_FIXUP=2, DATA_AFTER_REFINE=3, DATA_RESTRICT=4, COARSE_FINE_FLUX=5,
                    POT_FOR_POISSON=6, POT_AFTER_REFINE=7 };
#else
enum GetBufMode_t { DATA_GENERAL=1, DATA_AFTER_FIXUP=2, DATA_AFTER_REFINE=3, DATA_RESTRICT=4, COARSE_FINE_FLUX=5};
#endif // #ifdef GRAVITY ... else ...


// options of the fluid boundary condition
enum OptFluBC_t { BC_FLU_NONE=0, BC_FLU_PERIODIC=1, BC_FLU_OUTFLOW=2, BC_FLU_REFLECTING=3, BC_FLU_USER=4 };


// options of the gravity boundary condition
#ifdef GRAVITY
enum OptPotBC_t { BC_POT_NONE=0, BC_POT_PERIODIC=1, BC_POT_ISOLATED=2 };
#else
enum OptPotBC_t { BC_POT_NONE=0 };
#endif


// options of particles
#ifdef PARTICLE
enum ParInit_t       { PAR_INIT_NONE=0, PAR_INIT_BY_FUNCTION=1, PAR_INIT_BY_RESTART=2, PAR_INIT_BY_FILE=3 };
enum ParInterp_t     { PAR_INTERP_DEFAULT=-1, PAR_INTERP_NONE=0, PAR_INTERP_NGP=1, PAR_INTERP_CIC=2, PAR_INTERP_TSC=3 };
enum ParInteg_t      { PAR_INTEG_DEFAULT=-1, PAR_INTEG_NONE=0, PAR_INTEG_EULER=1, PAR_INTEG_KDK=2 };
enum ParUpStep_t     { PAR_UPSTEP_PRED=1, PAR_UPSTEP_CORR=2, PAR_UPSTEP_ACC_ONLY };
enum ParSync_t       { PAR_SYNC_DEFAULT=-1, PAR_SYNC_NONE=0, PAR_SYNC_TEMP=1, PAR_SYNC_FORCE=2 };
enum ParOutputDens_t { PAR_OUTPUT_DENS_NONE=0, PAR_OUTPUT_DENS_PAR_ONLY=1, PAR_OUTPUT_DENS_TOTAL=2 };
#endif // #ifdef PARTICLE


// options of the gravity types (this type needs to be defined for the Fluid solver even when GRAVITY is off)
//#ifdef GRAVITY
enum OptGravityType_t { GRAVITY_NONE=0, GRAVITY_SELF=1, GRAVITY_EXTERNAL=2, GRAVITY_BOTH=3 };
//#endif


// options of the form of the Lohner's error estimator
enum OptLohnerForm_t { LOHNER_DEFAULT=-1, LOHNER_FLASH1=0, LOHNER_FLASH2=1, LOHNER_FORM_INV1=2, LOHNER_FORM_INV2=3 };


// options for OPT__1ST_FLUX_CORR_SCHEME
#if ( MODEL == HYDRO || MODEL == MHD )
enum OptRSolver_t { RSOLVER_DEFAULT=-1, RSOLVER_NONE=0, RSOLVER_ROE=1, RSOLVER_HLLC=2, RSOLVER_HLLE=3 };
#endif

// options for OPT__CORR_AFTER_ALL_SYNC
enum OptCorrAfterSync_t { CORR_DEFAULT=-1, CORR_NONE=0, CORR_EVERY_STEP=1, CORR_BEFORE_DUMP=2 };



#endif  // #ifndef __TYPEDEF_H__
