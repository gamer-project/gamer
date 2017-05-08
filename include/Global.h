#ifndef __GLOBAL_H__
#define __GLOBAL_H__



#include "Macro.h"
#include "AMR.h"


// **********************************************************************************************************
// ** Variables in CAPITAL letters are loaded from the parameter file "Input__Parameter". Please refer to  **
// ** that file for detailed descritpions of each variable.                                                **
// **********************************************************************************************************


// 1. common global variables
// ============================================================================================================
extern AMR_t     *amr;                                // AMR tree structure
extern double     Time[NLEVEL];                       // "present"  physical time at each level
extern double     Time_Prev[NLEVEL];                  // "previous" physical time at each level
extern long       AdvanceCounter[NLEVEL];             // number of sub-steps that each level has been evolved
extern long       NCorrUnphy[NLEVEL];                 // number of cells corrected by either OPT__1ST_FLUX_CORR or MIN_DENS/PRES
extern long       Step;                               // number of main steps
extern double     dTime_Base;                         // time interval to update physical time at the base level

extern real       MinDtInfo_Fluid[NLEVEL];            // maximum CFL speed at each level
extern double     FlagTable_Rho        [NLEVEL-1];    // refinement criterion of density
extern double     FlagTable_RhoGradient[NLEVEL-1];    // refinement criterion of density gradient
extern double     FlagTable_Lohner     [NLEVEL-1][4]; // refinement criterion based on Lohner's error estimator
extern double     FlagTable_User       [NLEVEL-1];    // user-defined refinement criterion
extern double    *DumpTable;                          // dump table recording the physical times to output data
extern int        DumpTable_NDump;                    // number of data dumps in the dump table
extern int        DumpID;                             // index of the current output file
extern double     DumpTime;                           // time of the next dump (for OPT__OUTPUT_MODE=1)

extern int        MPI_Rank;                           // MPI rank ID in the MPI_COMM_WORLD
extern int        MPI_Rank_X[3];                      // order of this MPI process in x/y/z directions
extern int        MPI_SibRank[26];                    // sibling MPI rank ID: same order as sibling[26]
extern int        NX0[3];                             //  of base-level cells per process in x/y/z directions
extern int        NPatchTotal[NLEVEL];                // total number of patches in all ranks
extern int       *BaseP;                              // table recording the IDs of the base-level patches
extern int        Flu_ParaBuf;                        // number of parallel buffers to exchange all fluid
                                                      // variables for the fluid solver and fluid refinement

extern char      *PassiveFieldName_Grid[NCOMP_PASSIVE];
extern int        PassiveNorm_NVar;
extern int        PassiveNorm_VarIdx[NCOMP_PASSIVE];

extern double     BOX_SIZE, DT__FLUID, DT__FLUID_INIT, END_T, OUTPUT_DT;
extern long int   END_STEP;
extern int        NX0_TOT[3], OUTPUT_STEP, REGRID_COUNT, FLU_GPU_NPGROUP, OMP_NTHREAD;
extern int        MPI_NRank, MPI_NRank_X[3], GPU_NSTREAM, FLAG_BUFFER_SIZE, MAX_LEVEL;

extern int        OPT__UM_START_LEVEL, OPT__UM_START_NVAR, OPT__GPUID_SELECT, OPT__PATCH_COUNT;
extern int        INIT_DUMPID, INIT_SUBSAMPLING_NCELL, OPT__TIMING_BARRIER, OPT__REUSE_MEMORY, RESTART_LOAD_NRANK;
extern double     OUTPUT_PART_X, OUTPUT_PART_Y, OUTPUT_PART_Z;
extern double     OPT__CK_MEMFREE, INT_MONO_COEFF, UNIT_L, UNIT_M, UNIT_T, UNIT_V, UNIT_D, UNIT_E, UNIT_P;
extern bool       OPT__FLAG_RHO, OPT__FLAG_RHO_GRADIENT, OPT__FLAG_USER, OPT__FLAG_LOHNER_DENS, OPT__FLAG_REGION;
extern bool       OPT__DT_USER, OPT__RECORD_DT, OPT__RECORD_MEMORY, OPT__ADAPTIVE_DT, OPT__MEMORY_POOL;
extern bool       OPT__FIXUP_RESTRICT, OPT__INIT_RESTRICT, OPT__VERBOSE, OPT__MANUAL_CONTROL, OPT__UNIT;
extern bool       OPT__INT_TIME, OPT__OUTPUT_TEST_ERROR, OPT__OUTPUT_BASE, OPT__OVERLAP_MPI, OPT__TIMING_BALANCE;
extern bool       OPT__OUTPUT_BASEPS, OPT__CK_REFINE, OPT__CK_PROPER_NESTING, OPT__CK_FINITE, OPT__RECORD_PERFORMANCE;
extern bool       OPT__CK_RESTRICT, OPT__CK_PATCH_ALLOCATE, OPT__FIXUP_FLUX, OPT__CK_FLUX_ALLOCATE, OPT__CK_NORMALIZE_PASSIVE;
extern bool       OPT__UM_START_DOWNGRADE, OPT__UM_START_REFINE, OPT__UM_FACTOR_5OVER3, OPT__TIMING_MPI;
extern bool       OPT__CK_CONSERVATION, OPT__RESET_FLUID, OPT__RECORD_USER, OPT__NORMALIZE_PASSIVE;

extern OptInit_t          OPT__INIT;
extern OptRestartH_t      OPT__RESTART_HEADER;
extern IntScheme_t        OPT__FLU_INT_SCHEME, OPT__REF_FLU_INT_SCHEME;
extern OptOutputFormat_t  OPT__OUTPUT_TOTAL;
extern OptOutputPart_t    OPT__OUTPUT_PART;
extern OptOutputMode_t    OPT__OUTPUT_MODE;
extern OptFluBC_t         OPT__BC_FLU[6];          // boundary conditions of fluid at (-x,+x,-y,+y,-z,+z) faces
extern OptLohnerForm_t    OPT__FLAG_LOHNER_FORM;
extern OptCorrAfterSync_t OPT__CORR_AFTER_ALL_SYNC;



// 2. global variables for different applications
// ============================================================================================================
// (2-1) fluid solver in different models
#if   ( MODEL == HYDRO )
extern double           FlagTable_PresGradient[NLEVEL-1], FlagTable_Vorticity[NLEVEL-1], FlagTable_Jeans[NLEVEL-1];
extern double           GAMMA, MINMOD_COEFF, EP_COEFF, MOLECULAR_WEIGHT;
extern LR_Limiter_t     OPT__LR_LIMITER;
extern WAF_Limiter_t    OPT__WAF_LIMITER;
extern Opt1stFluxCorr_t OPT__1ST_FLUX_CORR;
extern OptRSolver_t     OPT__1ST_FLUX_CORR_SCHEME;
extern bool             OPT__FLAG_PRES_GRADIENT, OPT__FLAG_LOHNER_ENGY, OPT__FLAG_LOHNER_PRES, OPT__FLAG_LOHNER_TEMP;
extern bool             OPT__FLAG_VORTICITY, OPT__FLAG_JEANS;
extern int              OPT__CK_NEGATIVE;
extern double           MIN_DENS, MIN_PRES;
#ifdef DUAL_ENERGY
extern double           DUAL_ENERGY_SWITCH;
#endif

#elif ( MODEL == MHD )
#warning WAIT MHD !!!
extern double           MIN_DENS, MIN_PRES;
#ifdef DUAL_ENERGY
extern double           DUAL_ENERGY_SWITCH;
#endif

#elif ( MODEL == ELBDM )
extern double           DT__PHASE, FlagTable_EngyDensity[NLEVEL-1][2];
extern bool             OPT__FLAG_ENGY_DENSITY, OPT__INT_PHASE, ELBDM_TAYLOR3_AUTO;
extern double           ELBDM_TAYLOR3_COEFF, ELBDM_MASS, ELBDM_PLANCK_CONST, ELBDM_ETA, MIN_DENS;
extern real             MinDtInfo_Phase[NLEVEL];            // maximum time derivative of phase at each level
#ifdef QUARTIC_SELF_INTERACTION
extern double           ELBDM_LAMBDA;
#endif

#else
#  error : ERROR : unsupported MODEL !!
#endif // MODEL


// (2-2) self-gravity
// ============================================================================================================
#ifdef GRAVITY
extern double     AveDensity_Init;                    // initial average mass density (in all levels)
extern real       MinDtInfo_Gravity[NLEVEL];          // maximum gravitational acceleration at each level
extern int        Pot_ParaBuf;                        // number of parallel buffers to exchange potential for the
                                                      // Poisson/Gravity solvers and the potential refinement
extern int        Rho_ParaBuf;                        // number of parallel buffers to exchange density for the Poisson solver

extern real      *GreenFuncK;
extern double     GFUNC_COEFF0;
extern double     DT__GRAVITY;
extern double     NEWTON_G;
extern int        POT_GPU_NPGROUP;
extern bool       OPT__OUTPUT_POT, OPT__GRA_P5_GRADIENT, OPT__EXTERNAL_POT;
extern double     SOR_OMEGA;
extern int        SOR_MAX_ITER, SOR_MIN_ITER;
extern double     MG_TOLERATED_ERROR;
extern int        MG_MAX_ITER, MG_NPRE_SMOOTH, MG_NPOST_SMOOTH;

extern IntScheme_t      OPT__POT_INT_SCHEME, OPT__RHO_INT_SCHEME, OPT__GRA_INT_SCHEME, OPT__REF_POT_INT_SCHEME;
extern OptPotBC_t       OPT__BC_POT;
extern OptGravityType_t OPT__GRAVITY_TYPE;
#endif


// (2-3) cosmology simulations
// ============================================================================================================
#ifdef COMOVING
extern double     A_INIT, OMEGA_M0, DT__MAX_DELTA_A, HUBBLE0;
#endif


// (2-4) load balance
// ============================================================================================================
#ifdef LOAD_BALANCE
extern double     LB_INPUT__WLI_MAX;                  // LB->WLI_Max loaded from "Input__Parameter"
#ifdef PARTICLE
extern double     LB_INPUT__PAR_WEIGHT;               // LB->Par_Weight loaded from "Input__Parameter"
#endif
extern bool       OPT__RECORD_LOAD_BALANCE;
#endif


// (2-5) particle
// ============================================================================================================
#ifdef PARTICLE
extern double          DT__PARVEL, DT__PARVEL_MAX, DT__PARACC;
extern real            MinDtInfo_ParVelAcc[2][NLEVEL];
extern bool            OPT__OUTPUT_PAR_TEXT, OPT__CK_PARTICLE, OPT__FLAG_NPAR_CELL, OPT__FLAG_PAR_MASS_CELL;
extern int             OPT__PARTICLE_COUNT, OPT__FLAG_NPAR_PATCH, FlagTable_NParPatch[NLEVEL-1], FlagTable_NParCell[NLEVEL-1];
extern double          FlagTable_ParMassCell[NLEVEL-1];
extern ParOutputDens_t OPT__OUTPUT_PAR_DENS;
extern char           *PassiveFieldName_Par[PAR_NPASSIVE];
#endif


// (2-6) yt inline analysis
// ============================================================================================================
#ifdef SUPPORT_LIBYT
extern char            YT_SCRIPT[MAX_STRING];
extern yt_verbose      YT_VERBOSE;
#endif



// 3. CPU (host) arrays for transferring data between CPU and GPU
// ============================================================================================================
extern real       (*h_Flu_Array_F_In [2])[FLU_NIN ][  FLU_NXT   *FLU_NXT   *FLU_NXT   ];
extern real       (*h_Flu_Array_F_Out[2])[FLU_NOUT][8*PATCH_SIZE*PATCH_SIZE*PATCH_SIZE];
extern real       (*h_Flux_Array[2])[9][NFLUX_TOTAL][4*PATCH_SIZE*PATCH_SIZE];
extern double     (*h_Corner_Array_F [2])[3];
extern real        *h_MinDtInfo_Fluid_Array[2];
#ifdef DUAL_ENERGY
extern char       (*h_DE_Array_F_Out[2])[8*PATCH_SIZE*PATCH_SIZE*PATCH_SIZE];
#endif

#ifdef GRAVITY
extern real       (*h_Rho_Array_P    [2])[RHO_NXT][RHO_NXT][RHO_NXT];
extern real       (*h_Pot_Array_P_In [2])[POT_NXT][POT_NXT][POT_NXT];
extern real       (*h_Pot_Array_P_Out[2])[GRA_NXT][GRA_NXT][GRA_NXT];
extern real       (*h_Flu_Array_G    [2])[GRA_NIN][PATCH_SIZE][PATCH_SIZE][PATCH_SIZE];
extern double     (*h_Corner_Array_G [2])[3];
#ifdef DUAL_ENERGY
extern char       (*h_DE_Array_G     [2])[PATCH_SIZE][PATCH_SIZE][PATCH_SIZE];
#endif

#ifdef UNSPLIT_GRAVITY
extern real       (*h_Pot_Array_USG_F[2])[USG_NXT_F ][USG_NXT_F ][USG_NXT_F ];
extern real       (*h_Pot_Array_USG_G[2])[USG_NXT_G ][USG_NXT_G ][USG_NXT_G ];
extern real       (*h_Flu_Array_USG_G[2])[GRA_NIN-1][PS1][PS1][PS1];
#endif
#endif



// 4. GPU (device) global memory arrays and timers
// ============================================================================================================
/*** These global variables are NOT included here. Instead, they are included by individual files
     only if necessary. ***/



#endif // #ifndef __GLOBAL_H__
