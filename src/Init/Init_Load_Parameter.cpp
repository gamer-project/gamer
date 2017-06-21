#include "Copyright.h"
#include "GAMER.h"
#include "ReadPara.h"




//-------------------------------------------------------------------------------------------------------
// Function    :  Init_Load_Parameter
// Description :  Load the runtime parameters
//
// Note        :  Currently the filename is fixed to "Input__Parameter"
//-------------------------------------------------------------------------------------------------------
void Init_Load_Parameter()
{

   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ...\n", __FUNCTION__ );


   const char  FileName[] = "Input__Parameter";
   ReadPara_t *ReadPara   = new ReadPara_t;

// add parameters in the following format:
// --> note that VARIABLE, DEFAULT, MIN, and MAX must have the same data type
// --> some handy constants (e.g., NoMin_int, Eps_float, ...) are defined in "ReadPara.h"
// ********************************************************************************************************************************
// ReadPara->Add( "KEY_IN_THE_INPUT_TABLE",     &VARIABLE,                       DEFAULT,          MIN,           MAX            );
// ********************************************************************************************************************************

// simulation scale
   ReadPara->Add( "BOX_SIZE",                   &BOX_SIZE,                       -1.0,             Eps_double,    NoMax_double   );
   ReadPara->Add( "NX0_TOT_X",                  &NX0_TOT[0],                     -1,               PS2,           NoMax_int      );
   ReadPara->Add( "NX0_TOT_Y",                  &NX0_TOT[1],                     -1,               PS2,           NoMax_int      );
   ReadPara->Add( "NX0_TOT_Z",                  &NX0_TOT[2],                     -1,               PS2,           NoMax_int      );
   ReadPara->Add( "MPI_NRANK",                  &MPI_NRank,                      -1,               1,             NoMax_int      );
// do not check MPI_NRANK_X/Y/Z since they can be negative during restart and will be deprecated in the future
   ReadPara->Add( "MPI_NRANK_X",                &MPI_NRank_X[0],                 -1,               NoMin_int,     NoMax_int      );
   ReadPara->Add( "MPI_NRANK_Y",                &MPI_NRank_X[1],                 -1,               NoMin_int,     NoMax_int      );
   ReadPara->Add( "MPI_NRANK_Z",                &MPI_NRank_X[2],                 -1,               NoMin_int,     NoMax_int      );
// do not check OMP_NTHREAD since it may be reset by Init_ResetDefaultParameter()
   ReadPara->Add( "OMP_NTHREAD",                &OMP_NTHREAD,                    -1,               NoMin_int,     NoMax_int      );
// do not check END_T and END_STEP since they may be reset by test problems or restart
   ReadPara->Add( "END_T",                      &END_T,                          -1.0,             NoMin_double,  NoMax_double   );
   ReadPara->Add( "END_STEP",                   &END_STEP,                       -1L,              NoMin_long,    NoMax_long     );


// test problems
   ReadPara->Add( "TESTPROB_ID",                &TESTPROB_ID,                     0,               0,             NoMax_int      );


// code units
   ReadPara->Add( "OPT__UNIT",                  &OPT__UNIT,                       false,           Useless_bool,  Useless_bool   );
// do not check units since they are validated by Init_Unit()
   ReadPara->Add( "UNIT_L",                     &UNIT_L,                         -1.0,             NoMin_double,  NoMax_double   );
   ReadPara->Add( "UNIT_M",                     &UNIT_M,                         -1.0,             NoMin_double,  NoMax_double   );
   ReadPara->Add( "UNIT_T",                     &UNIT_T,                         -1.0,             NoMin_double,  NoMax_double   );
   ReadPara->Add( "UNIT_V",                     &UNIT_V,                         -1.0,             NoMin_double,  NoMax_double   );
   ReadPara->Add( "UNIT_D",                     &UNIT_D,                         -1.0,             NoMin_double,  NoMax_double   );


// boundary conditions
   ReadPara->Add( "OPT__BC_FLU_XM",             &OPT__BC_FLU[0],                 -1,               1,             4              );
   ReadPara->Add( "OPT__BC_FLU_XP",             &OPT__BC_FLU[1],                 -1,               1,             4              );
   ReadPara->Add( "OPT__BC_FLU_YM",             &OPT__BC_FLU[2],                 -1,               1,             4              );
   ReadPara->Add( "OPT__BC_FLU_YP",             &OPT__BC_FLU[3],                 -1,               1,             4              );
   ReadPara->Add( "OPT__BC_FLU_ZM",             &OPT__BC_FLU[4],                 -1,               1,             4              );
   ReadPara->Add( "OPT__BC_FLU_ZP",             &OPT__BC_FLU[5],                 -1,               1,             4              );
#  ifdef GRAVITY
   ReadPara->Add( "OPT__BC_POT",                &OPT__BC_POT,                    -1,               1,             2              );
// do not check GFUNC_COEFF0 since it may be reset by Init_ResetDefaultParameter()
   ReadPara->Add( "GFUNC_COEFF0",               &GFUNC_COEFF0,                   -1.0,             NoMin_double,  NoMax_double   );
#  endif


// particle
#  ifdef PARTICLE
// do no check PAR_NPAR since it may be reset by restart
   ReadPara->Add( "PAR_NPAR",                   &amr->Par->NPar_Active_AllRank,  -1L,              NoMin_long,    NoMax_long     );
   ReadPara->Add( "PAR_INIT",                   &amr->Par->Init,                 -1,               1,             3              );
   ReadPara->Add( "PAR_INTERP",                 &amr->Par->Interp,                PAR_INTERP_TSC,  1,             3              );
   ReadPara->Add( "PAR_INTEG",                  &amr->Par->Integ,                 PAR_INTEG_KDK,   1,             2              );
   ReadPara->Add( "PAR_IMPROVE_ACC",            &amr->Par->ImproveAcc,            true,            Useless_bool,  Useless_bool   );
   ReadPara->Add( "PAR_PREDICT_POS",            &amr->Par->PredictPos,            true,            Useless_bool,  Useless_bool   );
// do not check PAR_REMOVE_CELL since it may be reset by Init_ResetDefaultParameter()
   ReadPara->Add( "PAR_REMOVE_CELL",            &amr->Par->RemoveCell,           -1.0,             NoMin_double,  NoMax_double   );
#  endif // #ifdef PARTICLE


// cosmology
#  ifdef COMOVING
   ReadPara->Add( "A_INIT",                     &A_INIT,                         -1.0,             Eps_double,    NoMax_double   );
   ReadPara->Add( "OMEGA_M0",                   &OMEGA_M0,                       -1.0,             0.0,           1.0            );
   ReadPara->Add( "HUBBLE0",                    &HUBBLE0,                        -1.0,             Eps_double,    1.0            );
#  endif


// time-step
// do not check DT__FLUID/FLUID_INIT/GRAVITY/PARVEL_MAX since they may be reset by Init_ResetDefaultParameter()
   ReadPara->Add( "DT__FLUID",                  &DT__FLUID,                      -1.0,             NoMin_double,  NoMax_double   );
   ReadPara->Add( "DT__FLUID_INIT",             &DT__FLUID_INIT,                 -1.0,             NoMin_double,  NoMax_double   );
#  ifdef GRAVITY
   ReadPara->Add( "DT__GRAVITY",                &DT__GRAVITY,                    -1.0,             NoMin_double,  NoMax_double   );
#  endif
#  if ( MODEL == ELBDM )
   ReadPara->Add( "DT__PHASE",                  &DT__PHASE,                       0.0,             0.0,           NoMax_double   );
#  endif
#  ifdef PARTICLE
   ReadPara->Add( "DT__PARVEL",                 &DT__PARVEL,                      0.5,             0.0,           NoMax_double   );
   ReadPara->Add( "DT__PARVEL_MAX",             &DT__PARVEL_MAX,                 -1.0,             NoMin_double,  NoMax_double   );
   ReadPara->Add( "DT__PARACC",                 &DT__PARACC,                      0.5,             0.0,           NoMax_double   );
#  endif
#  ifdef COMOVING
   ReadPara->Add( "DT__MAX_DELTA_A",            &DT__MAX_DELTA_A,                 0.01,            0.0,           NoMax_double   );
#  endif
   ReadPara->Add( "OPT__DT_USER",               &OPT__DT_USER,                    false,           Useless_bool,  Useless_bool   );
   ReadPara->Add( "OPT__ADAPTIVE_DT",           &OPT__ADAPTIVE_DT,                false,           Useless_bool,  Useless_bool   );
   ReadPara->Add( "OPT__RECORD_DT",             &OPT__RECORD_DT,                  true,            Useless_bool,  Useless_bool   );


// grid refinement
   ReadPara->Add( "REGRID_COUNT",               &REGRID_COUNT,                    4,               1,             NoMax_int      );
   ReadPara->Add( "FLAG_BUFFER_SIZE",           &FLAG_BUFFER_SIZE,                PS1,             0,             PS1            );
   ReadPara->Add( "MAX_LEVEL",                  &MAX_LEVEL,                       TOP_LEVEL,       0,             TOP_LEVEL      );
   ReadPara->Add( "OPT__FLAG_RHO",              &OPT__FLAG_RHO,                   false,           Useless_bool,  Useless_bool   );
   ReadPara->Add( "OPT__FLAG_RHO_GRADIENT",     &OPT__FLAG_RHO_GRADIENT,          false,           Useless_bool,  Useless_bool   );
   ReadPara->Add( "OPT__FLAG_LOHNER_DENS",      &OPT__FLAG_LOHNER_DENS,           false,           Useless_bool,  Useless_bool   );
#  if ( MODEL == HYDRO   ||  MODEL == MHD )
   ReadPara->Add( "OPT__FLAG_PRES_GRADIENT",    &OPT__FLAG_PRES_GRADIENT,         false,           Useless_bool,  Useless_bool   );
   ReadPara->Add( "OPT__FLAG_VORTICITY",        &OPT__FLAG_VORTICITY,             false,           Useless_bool,  Useless_bool   );
   ReadPara->Add( "OPT__FLAG_JEANS",            &OPT__FLAG_JEANS,                 false,           Useless_bool,  Useless_bool   );
   ReadPara->Add( "OPT__FLAG_LOHNER_ENGY",      &OPT__FLAG_LOHNER_ENGY,           false,           Useless_bool,  Useless_bool   );
   ReadPara->Add( "OPT__FLAG_LOHNER_PRES",      &OPT__FLAG_LOHNER_PRES,           false,           Useless_bool,  Useless_bool   );
   ReadPara->Add( "OPT__FLAG_LOHNER_TEMP",      &OPT__FLAG_LOHNER_TEMP,           false,           Useless_bool,  Useless_bool   );
#  endif
#  if ( MODEL == ELBDM )
   ReadPara->Add( "OPT__FLAG_ENGY_DENSITY",     &OPT__FLAG_ENGY_DENSITY,          false,           Useless_bool,  Useless_bool   );
#  endif
   ReadPara->Add( "OPT__FLAG_LOHNER_FORM",      &OPT__FLAG_LOHNER_FORM,           LOHNER_FLASH2,   1,             4              );
   ReadPara->Add( "OPT__FLAG_USER",             &OPT__FLAG_USER,                  false,           Useless_bool,  Useless_bool   );
   ReadPara->Add( "OPT__FLAG_REGION",           &OPT__FLAG_REGION,                false,           Useless_bool,  Useless_bool   );
#  ifdef PARTICLE
   ReadPara->Add( "OPT__FLAG_NPAR_PATCH",       &OPT__FLAG_NPAR_PATCH,            0,               0,             2              );
   ReadPara->Add( "OPT__FLAG_NPAR_CELL",        &OPT__FLAG_NPAR_CELL,             false,           Useless_bool,  Useless_bool   );
   ReadPara->Add( "OPT__FLAG_PAR_MASS_CELL",    &OPT__FLAG_PAR_MASS_CELL,         false,           Useless_bool,  Useless_bool   );
#  endif
   ReadPara->Add( "OPT__PATCH_COUNT",           &OPT__PATCH_COUNT,                1,               0,             2              );
#  ifdef PARTICLE
   ReadPara->Add( "OPT__PARTICLE_COUNT",        &OPT__PARTICLE_COUNT,             1,               0,             2              );
#  endif
   ReadPara->Add( "OPT__REUSE_MEMORY",          &OPT__REUSE_MEMORY,               2,               0,             2              );
   ReadPara->Add( "OPT__MEMORY_POOL",           &OPT__MEMORY_POOL,                false,           Useless_bool,  Useless_bool   );


// load balance
#  ifdef LOAD_BALANCE
   ReadPara->Add( "LB_INPUT__WLI_MAX",          &LB_INPUT__WLI_MAX,               0.1,             0.0,           NoMax_double   );
#  ifdef PARTICLE
   ReadPara->Add( "LB_INPUT__PAR_WEIGHT",       &LB_INPUT__PAR_WEIGHT,            0.0,             0.0,           NoMax_double   );
#  endif
   ReadPara->Add( "OPT__RECORD_LOAD_BALANCE",   &OPT__RECORD_LOAD_BALANCE,        true,            Useless_bool,  Useless_bool   );
#  endif


// Grackle
#  ifdef SUPPORT_GRACKLE
   ReadPara->Add( "GRACKLE_ENABLE",             &GRACKLE_ENABLE,                  true,            Useless_bool,  Useless_bool   );
   ReadPara->Add( "GRACKLE_VERBOSE",            &GRACKLE_VERBOSE,                 true,            Useless_bool,  Useless_bool   );
   ReadPara->Add( "GRACKLE_COOLING",            &GRACKLE_COOLING,                 true,            Useless_bool,  Useless_bool   );
   ReadPara->Add( "GRACKLE_PRIMORDIAL",         &GRACKLE_PRIMORDIAL,              0,               0,             3              );
   ReadPara->Add( "GRACKLE_METAL",              &GRACKLE_METAL,                   false,           Useless_bool,  Useless_bool   );
   ReadPara->Add( "GRACKLE_UV",                 &GRACKLE_UV,                      false,           Useless_bool,  Useless_bool   );
   ReadPara->Add( "GRACKLE_CLOUDY_TABLE",        GRACKLE_CLOUDY_TABLE,            Useless_str,     Useless_str,   Useless_str    );
#  endif


// fluid solvers in HYDRO and MHD
#  if ( MODEL == HYDRO )
   ReadPara->Add( "GAMMA",                      &GAMMA,                           5.0/3.0,         Eps_double,    NoMax_double   );
   ReadPara->Add( "MOLECULAR_WEIGHT",           &MOLECULAR_WEIGHT,                0.6,             Eps_double,    NoMax_double   );
   ReadPara->Add( "MINMOD_COEFF",               &MINMOD_COEFF,                    2.0,             1.0,           2.0            );
   ReadPara->Add( "EP_COEFF",                   &EP_COEFF,                        1.25,            1.0,           NoMax_double   );
   ReadPara->Add( "OPT__LR_LIMITER",            &OPT__LR_LIMITER,                 VL_GMINMOD,      0,             5              );
   ReadPara->Add( "OPT__WAF_LIMITER",           &OPT__WAF_LIMITER,                WAF_VANLEER,     0,             4              );
   ReadPara->Add( "OPT__1ST_FLUX_CORR",         &OPT__1ST_FLUX_CORR,              FIRST_FLUX_CORR_3D1D, 0,        2              );
   ReadPara->Add( "OPT__1ST_FLUX_CORR_SCHEME",  &OPT__1ST_FLUX_CORR_SCHEME,       RSOLVER_1ST_ROE, 0,             3              );
#  ifdef DUAL_ENERGY
   ReadPara->Add( "DUAL_ENERGY_SWITCH",         &DUAL_ENERGY_SWITCH,              2.0e-2,          0.0,           NoMax_double   );
#  endif

#  elif ( MODEL == MHD )
#  warning : WAIT MHD !!!
#  endif // #if ( MODEL == HYDRO/MHD )


// fluid solver in ELBDM
#  if ( MODEL == ELBDM )
   ReadPara->Add( "ELBDM_MASS",                 &ELBDM_MASS,                     -1.0,             Eps_double,    NoMax_double   );
// do not check ELBDM_PLANCK_CONST since it may be reset by Init_Unit()
   ReadPara->Add( "ELBDM_PLANCK_CONST",         &ELBDM_PLANCK_CONST,             -1.0,             NoMin_double,  NoMax_double   );
#  ifdef QUARTIC_SELF_INTERACTION
   ReadPara->Add( "ELBDM_LAMBDA",               &ELBDM_LAMBDA,                    1.0,             NoMin_double,  NoMax_double   );
#  endif
   ReadPara->Add( "ELBDM_TAYLOR3_COEFF",        &ELBDM_TAYLOR3_COEFF,             1.0/6.0,         NoMin_double,  NoMax_double   );
   ReadPara->Add( "ELBDM_TAYLOR3_AUTO",         &ELBDM_TAYLOR3_AUTO,              true,            Useless_bool,  Useless_bool   );
#  endif // #if ( MODEL == ELBDM )


// fluid solvers in all models
// do not check FLU_GPU_NPGROUP and GPU_NSTREAM since they may be reset by either Init_Set_Default_SOR_Parameter or CUAPI_Set_Default_GPU_Parameter()
   ReadPara->Add( "FLU_GPU_NPGROUP",            &FLU_GPU_NPGROUP,                -1,               NoMin_int,     NoMax_int      );
   ReadPara->Add( "GPU_NSTREAM",                &GPU_NSTREAM,                    -1,               NoMin_int,     NoMax_int      );
   ReadPara->Add( "OPT__FIXUP_FLUX",            &OPT__FIXUP_FLUX,                 true,            Useless_bool,  Useless_bool   );
   ReadPara->Add( "OPT__FIXUP_RESTRICT",        &OPT__FIXUP_RESTRICT,             true,            Useless_bool,  Useless_bool   );
// do not check OPT__CORR_AFTER_ALL_SYNC since it may be reset by Init_ResetDefaultParameter()
   ReadPara->Add( "OPT__CORR_AFTER_ALL_SYNC",   &OPT__CORR_AFTER_ALL_SYNC,       -1,               NoMin_int,     NoMax_int      );
   ReadPara->Add( "OPT__NORMALIZE_PASSIVE",     &OPT__NORMALIZE_PASSIVE,          true,            Useless_bool,  Useless_bool   );
   ReadPara->Add( "OPT__OVERLAP_MPI",           &OPT__OVERLAP_MPI,                false,           Useless_bool,  Useless_bool   );
   ReadPara->Add( "OPT__RESET_FLUID",           &OPT__RESET_FLUID,                false,           Useless_bool,  Useless_bool   );
#  if ( MODEL == HYDRO  ||  MODEL == MHD  ||  MODEL == ELBDM )
   ReadPara->Add( "MIN_DENS",                   &MIN_DENS,                        0.0,             0.0,           NoMax_double   );
#  endif
#  if ( MODEL == HYDRO  ||  MODEL == MHD )
   ReadPara->Add( "MIN_PRES",                   &MIN_PRES,                        0.0,             0.0,           NoMax_double   );
#  endif


// self-gravity
#  ifdef GRAVITY
// do not check NEWTON_G since it may be reset by Init_Unit()
   ReadPara->Add( "NEWTON_G",                   &NEWTON_G,                       -1.0,             NoMin_double,  NoMax_double   );
// do not check SOR_XXX since they may be reset by Init_Set_Default_SOR_Parameter()
   ReadPara->Add( "SOR_OMEGA",                  &SOR_OMEGA,                      -1.0,             NoMin_double,  NoMax_double   );
   ReadPara->Add( "SOR_MAX_ITER",               &SOR_MAX_ITER,                   -1,               NoMin_int,     NoMax_int      );
   ReadPara->Add( "SOR_MIN_ITER",               &SOR_MIN_ITER,                   -1,               NoMin_int,     NoMax_int      );
// do not check MG_XXX since they may be reset by Init_Set_Default_MG_Parameter()
   ReadPara->Add( "MG_MAX_ITER",                &MG_MAX_ITER,                    -1,               NoMin_int,     NoMax_int      );
   ReadPara->Add( "MG_NPRE_SMOOTH",             &MG_NPRE_SMOOTH,                 -1,               NoMin_int,     NoMax_int      );
   ReadPara->Add( "MG_NPOST_SMOOTH",            &MG_NPOST_SMOOTH,                -1,               NoMin_int,     NoMax_int      );
   ReadPara->Add( "MG_TOLERATED_ERROR",         &MG_TOLERATED_ERROR,             -1.0,             NoMin_double,  NoMax_double   );
// do not check POT_GPU_NPGROUP since it may be reset by either Init_Set_Default_SOR_Parameter or CUAPI_Set_Default_GPU_Parameter()
   ReadPara->Add( "POT_GPU_NPGROUP",            &POT_GPU_NPGROUP,                -1,               NoMin_int,     NoMax_int      );
   ReadPara->Add( "OPT__GRA_P5_GRADIENT",       &OPT__GRA_P5_GRADIENT,            false,           Useless_bool,  Useless_bool   );
   ReadPara->Add( "OPT__GRAVITY_TYPE",          &OPT__GRAVITY_TYPE,              -1,               1,             3              );
   ReadPara->Add( "OPT__EXTERNAL_POT",          &OPT__EXTERNAL_POT,               false,           Useless_bool,  Useless_bool   );
#  endif // #ifdef GRAVITY


// initialization
   ReadPara->Add( "OPT__INIT",                  &OPT__INIT,                      -1,               1,             3              );
   ReadPara->Add( "RESTART_LOAD_NRANK",         &RESTART_LOAD_NRANK,              1,               1,             NoMax_int      );
   ReadPara->Add( "OPT__RESTART_HEADER",        &OPT__RESTART_HEADER,             RESTART_HEADER_CHECK, 0,        1              );
   ReadPara->Add( "OPT__UM_START_LEVEL",        &OPT__UM_START_LEVEL,             0,               0,             TOP_LEVEL      );
   ReadPara->Add( "OPT__UM_START_NVAR",         &OPT__UM_START_NVAR,              1,               1,             NCOMP_TOTAL    );
   ReadPara->Add( "OPT__UM_START_DOWNGRADE",    &OPT__UM_START_DOWNGRADE,         true,            Useless_bool,  Useless_bool   );
   ReadPara->Add( "OPT__UM_START_REFINE",       &OPT__UM_START_REFINE,            true,            Useless_bool,  Useless_bool   );
   ReadPara->Add( "OPT__UM_FACTOR_5OVER3",      &OPT__UM_FACTOR_5OVER3,           false,           Useless_bool,  Useless_bool   );
   ReadPara->Add( "OPT__INIT_RESTRICT",         &OPT__INIT_RESTRICT,              true,            Useless_bool,  Useless_bool   );
   ReadPara->Add( "OPT__GPUID_SELECT",          &OPT__GPUID_SELECT,              -2,              -3,             NoMax_int      );
   ReadPara->Add( "INIT_SUBSAMPLING_NCELL",     &INIT_SUBSAMPLING_NCELL,          0,               0,             NoMax_int      );


// interpolation schemes
   ReadPara->Add( "OPT__INT_TIME",              &OPT__INT_TIME,                   true,            Useless_bool,  Useless_bool   );
#  if ( MODEL == ELBDM )
   ReadPara->Add( "OPT__INT_PHASE",             &OPT__INT_PHASE,                  true,            Useless_bool,  Useless_bool   );
#  endif
// do not check OPT__FLU_INT_SCHEME and OPT__REF_FLU_INT_SCHEME since they may be reset by Init_ResetDefaultParameter()
   ReadPara->Add( "OPT__FLU_INT_SCHEME",        &OPT__FLU_INT_SCHEME,             INT_DEFAULT,     NoMin_int,     NoMax_int      );
   ReadPara->Add( "OPT__REF_FLU_INT_SCHEME",    &OPT__REF_FLU_INT_SCHEME,         INT_DEFAULT,     NoMin_int,     NoMax_int      );
#  ifdef GRAVITY
   ReadPara->Add( "OPT__POT_INT_SCHEME",        &OPT__POT_INT_SCHEME,             INT_QUAD,        4,             5              );
   ReadPara->Add( "OPT__RHO_INT_SCHEME",        &OPT__RHO_INT_SCHEME,             INT_CQUAD,       1,             7              );
   ReadPara->Add( "OPT__GRA_INT_SCHEME",        &OPT__GRA_INT_SCHEME,             INT_QUAD,        1,             7              );
   ReadPara->Add( "OPT__REF_POT_INT_SCHEME",    &OPT__REF_POT_INT_SCHEME,         INT_QUAD,        1,             7              );
#  endif
   ReadPara->Add( "INT_MONO_COEFF",             &INT_MONO_COEFF,                  2.0,             1.0,           4.0            );


// data dump
   ReadPara->Add( "OPT__OUTPUT_TOTAL",          &OPT__OUTPUT_TOTAL,               1,               0,             2              );
   ReadPara->Add( "OPT__OUTPUT_PART",           &OPT__OUTPUT_PART,                0,               0,             7              );
   ReadPara->Add( "OPT__OUTPUT_USER",           &OPT__OUTPUT_USER,                false,           Useless_bool,  Useless_bool   );
#  ifdef PARTICLE
   ReadPara->Add( "OPT__OUTPUT_PAR_TEXT",       &OPT__OUTPUT_PAR_TEXT,            false,           Useless_bool,  Useless_bool   );
#  endif
   ReadPara->Add( "OPT__OUTPUT_BASEPS",         &OPT__OUTPUT_BASEPS,              false,           Useless_bool,  Useless_bool   );
   ReadPara->Add( "OPT__OUTPUT_BASE",           &OPT__OUTPUT_BASE,                false,           Useless_bool,  Useless_bool   );
#  ifdef GRAVITY
   ReadPara->Add( "OPT__OUTPUT_POT",            &OPT__OUTPUT_POT,                 false,           Useless_bool,  Useless_bool   );
#  endif
#  ifdef PARTICLE
   ReadPara->Add( "OPT__OUTPUT_PAR_DENS",       &OPT__OUTPUT_PAR_DENS,            PAR_OUTPUT_DENS_PAR_ONLY, 0,    2              );
#  endif
   ReadPara->Add( "OPT__OUTPUT_MODE",           &OPT__OUTPUT_MODE,               -1,               1,             3              );
// do not check OUTPUT_STEP and OUTPUT_DT since they depend on OPT__OUTPUT_MODE
   ReadPara->Add( "OUTPUT_STEP",                &OUTPUT_STEP,                    -1,               NoMin_int,     NoMax_int      );
   ReadPara->Add( "OUTPUT_DT",                  &OUTPUT_DT,                      -1.0,             NoMin_double,  NoMax_double   );
// do not check OUTPUT_PART_X/Y/Z since they depend on OPT__OUTPUT_PART
   ReadPara->Add( "OUTPUT_PART_X",              &OUTPUT_PART_X,                  -1.0,             NoMin_double,  NoMax_double   );
   ReadPara->Add( "OUTPUT_PART_Y",              &OUTPUT_PART_Y,                  -1.0,             NoMin_double,  NoMax_double   );
   ReadPara->Add( "OUTPUT_PART_Z",              &OUTPUT_PART_Z,                  -1.0,             NoMin_double,  NoMax_double   );
   ReadPara->Add( "INIT_DUMPID",                &INIT_DUMPID,                    -1,               NoMin_int,     NoMax_int      );


// yt inline analysis
#  ifdef SUPPORT_LIBYT
   ReadPara->Add( "YT_SCRIPT",                   YT_SCRIPT,                       Useless_str,     Useless_str,   Useless_str    );
   ReadPara->Add( "YT_VERBOSE",                 &YT_VERBOSE,                      1,               0,             3              );
#  endif


// miscellaneous
   ReadPara->Add( "OPT__VERBOSE",               &OPT__VERBOSE,                    false,           Useless_bool,  Useless_bool   );
// do not check OPT__TIMING_BARRIER since it depends on other options
   ReadPara->Add( "OPT__TIMING_BARRIER",        &OPT__TIMING_BARRIER,            -1,               NoMin_int,     NoMax_int      );
   ReadPara->Add( "OPT__TIMING_BALANCE",        &OPT__TIMING_BALANCE,             false,           Useless_bool,  Useless_bool   );
   ReadPara->Add( "OPT__TIMING_MPI",            &OPT__TIMING_MPI,                 false,           Useless_bool,  Useless_bool   );
   ReadPara->Add( "OPT__RECORD_MEMORY",         &OPT__RECORD_MEMORY,              true,            Useless_bool,  Useless_bool   );
   ReadPara->Add( "OPT__RECORD_PERFORMANCE",    &OPT__RECORD_PERFORMANCE,         true,            Useless_bool,  Useless_bool   );
   ReadPara->Add( "OPT__MANUAL_CONTROL",        &OPT__MANUAL_CONTROL,             true,            Useless_bool,  Useless_bool   );
   ReadPara->Add( "OPT__RECORD_USER",           &OPT__RECORD_USER,                false,           Useless_bool,  Useless_bool   );


// simulation checks
   ReadPara->Add( "OPT__CK_REFINE",             &OPT__CK_REFINE,                  false,           Useless_bool,  Useless_bool   );
   ReadPara->Add( "OPT__CK_PROPER_NESTING",     &OPT__CK_PROPER_NESTING,          false,           Useless_bool,  Useless_bool   );
   ReadPara->Add( "OPT__CK_CONSERVATION",       &OPT__CK_CONSERVATION,            false,           Useless_bool,  Useless_bool   );
   ReadPara->Add( "OPT__CK_NORMALIZE_PASSIVE",  &OPT__CK_NORMALIZE_PASSIVE,       false,           Useless_bool,  Useless_bool   );
   ReadPara->Add( "OPT__CK_RESTRICT",           &OPT__CK_RESTRICT,                false,           Useless_bool,  Useless_bool   );
   ReadPara->Add( "OPT__CK_FINITE",             &OPT__CK_FINITE,                  false,           Useless_bool,  Useless_bool   );
   ReadPara->Add( "OPT__CK_PATCH_ALLOCATE",     &OPT__CK_PATCH_ALLOCATE,          false,           Useless_bool,  Useless_bool   );
   ReadPara->Add( "OPT__CK_FLUX_ALLOCATE",      &OPT__CK_FLUX_ALLOCATE,           false,           Useless_bool,  Useless_bool   );
#  if ( MODEL == HYDRO  ||  MODEL == MHD )
   ReadPara->Add( "OPT__CK_NEGATIVE",           &OPT__CK_NEGATIVE,                0,               0,             3              );
#  endif
   ReadPara->Add( "OPT__CK_MEMFREE",            &OPT__CK_MEMFREE,                 1.0,             0.0,           NoMax_double   );
#  ifdef PARTICLE
   ReadPara->Add( "OPT__CK_PARTICLE",           &OPT__CK_PARTICLE,                false,           Useless_bool,  Useless_bool   );
#  endif



// load parameters
// ********************************************************************************************************************************
   ReadPara->Read( FileName );
// ********************************************************************************************************************************


// free memory
   delete ReadPara;


   if ( MPI_Rank == 0 )    Aux_Message( stdout, "%s ... done\n", __FUNCTION__ );

} // FUNCTION : Init_Load_Parameter

