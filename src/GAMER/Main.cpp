#include "GAMER.h"
#include "CUFLU.h"

#ifdef GRAVITY
#include "CUPOT.h"
#endif



// ***********************
// **  GLOBAL VARIABLES **
// ***********************

// 1. common global variables
// =======================================================================================================
AMR_t               *amr = NULL;

double               Time[NLEVEL]           = { 0.0 };
double               dTime_AllLv[NLEVEL]    = { 0.0 };
long                 AdvanceCounter[NLEVEL] = { 0 };
long                 NCorrUnphy[NLEVEL]     = { 0 };
long                 Step                   = 0;
int                  DumpID                 = 0;
double               DumpTime               = 0.0;

double               dTime_Base;
double               Time_Prev            [NLEVEL];
double               FlagTable_Rho        [NLEVEL-1];
double               FlagTable_RhoGradient[NLEVEL-1];
double               FlagTable_Lohner     [NLEVEL-1][4];
double               FlagTable_User       [NLEVEL-1];
double              *DumpTable = NULL;
int                  DumpTable_NDump;
int                  PassiveNorm_NVar;
int                  PassiveNorm_VarIdx[NCOMP_PASSIVE];

int                  MPI_Rank, MPI_Rank_X[3], MPI_SibRank[26], NX0[3], NPatchTotal[NLEVEL];
int                 *BaseP = NULL;
int                  Flu_ParaBuf;

char                *PassiveFieldName_Grid[NCOMP_PASSIVE];

double               BOX_SIZE, DT__FLUID, DT__FLUID_INIT, END_T, OUTPUT_DT, DT__SYNC_PARENT_LV, DT__SYNC_CHILDREN_LV;
long                 END_STEP;
int                  NX0_TOT[3], OUTPUT_STEP, REGRID_COUNT, FLU_GPU_NPGROUP, OMP_NTHREAD;
int                  MPI_NRank, MPI_NRank_X[3];
int                  GPU_NSTREAM, FLAG_BUFFER_SIZE, FLAG_BUFFER_SIZE_MAXM1_LV, FLAG_BUFFER_SIZE_MAXM2_LV, MAX_LEVEL;

IntScheme_t          OPT__FLU_INT_SCHEME, OPT__REF_FLU_INT_SCHEME;
double               OUTPUT_PART_X, OUTPUT_PART_Y, OUTPUT_PART_Z, AUTO_REDUCE_DT_FACTOR, AUTO_REDUCE_DT_FACTOR_MIN;
double               OPT__CK_MEMFREE, INT_MONO_COEFF, UNIT_L, UNIT_M, UNIT_T, UNIT_V, UNIT_D, UNIT_E, UNIT_P;
int                  OPT__UM_IC_LEVEL, OPT__UM_IC_NVAR, OPT__GPUID_SELECT, OPT__PATCH_COUNT;
int                  INIT_DUMPID, INIT_SUBSAMPLING_NCELL, OPT__TIMING_BARRIER, OPT__REUSE_MEMORY, RESTART_LOAD_NRANK;
bool                 OPT__FLAG_RHO, OPT__FLAG_RHO_GRADIENT, OPT__FLAG_USER, OPT__FLAG_LOHNER_DENS, OPT__FLAG_REGION;
bool                 OPT__DT_USER, OPT__RECORD_DT, OPT__RECORD_MEMORY, OPT__MEMORY_POOL, OPT__RESTART_RESET;
bool                 OPT__FIXUP_RESTRICT, OPT__INIT_RESTRICT, OPT__VERBOSE, OPT__MANUAL_CONTROL, OPT__UNIT;
bool                 OPT__INT_TIME, OPT__OUTPUT_USER, OPT__OUTPUT_BASE, OPT__OVERLAP_MPI, OPT__TIMING_BALANCE;
bool                 OPT__OUTPUT_BASEPS, OPT__CK_REFINE, OPT__CK_PROPER_NESTING, OPT__CK_FINITE, OPT__RECORD_PERFORMANCE;
bool                 OPT__CK_RESTRICT, OPT__CK_PATCH_ALLOCATE, OPT__FIXUP_FLUX, OPT__CK_FLUX_ALLOCATE, OPT__CK_NORMALIZE_PASSIVE;
bool                 OPT__UM_IC_DOWNGRADE, OPT__UM_IC_REFINE, OPT__TIMING_MPI;
bool                 OPT__CK_CONSERVATION, OPT__RESET_FLUID, OPT__RECORD_USER, OPT__NORMALIZE_PASSIVE, AUTO_REDUCE_DT;
bool                 OPT__OPTIMIZE_AGGRESSIVE, OPT__INIT_GRID_WITH_OMP;
TestProbID_t         TESTPROB_ID;
OptInit_t            OPT__INIT;
OptOutputFormat_t    OPT__OUTPUT_TOTAL;
OptOutputPart_t      OPT__OUTPUT_PART;
OptOutputMode_t      OPT__OUTPUT_MODE;
OptFluBC_t           OPT__BC_FLU[6];
OptLohnerForm_t      OPT__FLAG_LOHNER_FORM;
OptCorrAfterSync_t   OPT__CORR_AFTER_ALL_SYNC;
OptTimeStepLevel_t   OPT__DT_LEVEL;


// 2. global variables for different applications
// =======================================================================================================
// (2-1) fluid solver in different models
#if   ( MODEL == HYDRO )
double               FlagTable_PresGradient[NLEVEL-1], FlagTable_Vorticity[NLEVEL-1], FlagTable_Jeans[NLEVEL-1];
double               GAMMA, MINMOD_COEFF, EP_COEFF, MOLECULAR_WEIGHT;
LR_Limiter_t         OPT__LR_LIMITER;
WAF_Limiter_t        OPT__WAF_LIMITER;
Opt1stFluxCorr_t     OPT__1ST_FLUX_CORR;
OptRSolver1st_t      OPT__1ST_FLUX_CORR_SCHEME;
bool                 OPT__FLAG_PRES_GRADIENT, OPT__FLAG_LOHNER_ENGY, OPT__FLAG_LOHNER_PRES, OPT__FLAG_LOHNER_TEMP;
bool                 OPT__FLAG_VORTICITY, OPT__FLAG_JEANS, JEANS_MIN_PRES;
int                  OPT__CK_NEGATIVE, JEANS_MIN_PRES_LEVEL, JEANS_MIN_PRES_NCELL;
double               MIN_DENS, MIN_PRES;
#ifdef DUAL_ENERGY
double               DUAL_ENERGY_SWITCH;
#endif

#elif ( MODEL == MHD )
#warning : WAIT MHD !!!
double               MIN_DENS, MIN_PRES;
#ifdef DUAL_ENERGY
double               DUAL_ENERGY_SWITCH;
#endif

#elif ( MODEL == ELBDM )
double               DT__PHASE, FlagTable_EngyDensity[NLEVEL-1][2];
bool                 OPT__FLAG_ENGY_DENSITY, OPT__INT_PHASE;
bool                 ELBDM_TAYLOR3_AUTO;
double               ELBDM_TAYLOR3_COEFF;
double               ELBDM_MASS, ELBDM_PLANCK_CONST, ELBDM_ETA, MIN_DENS;
#ifdef QUARTIC_SELF_INTERACTION
double               ELBDM_LAMBDA;
#endif

#else
#error : unsupported MODEL !!
#endif // MODEL

// (2-2) self-gravity
#ifdef GRAVITY
double               AveDensity_Init = -1.0;    // initialize it as <= 0 to check if it is properly set later
int                  Pot_ParaBuf, Rho_ParaBuf;

real                *GreenFuncK      = NULL;
double               GFUNC_COEFF0;
double               DT__GRAVITY;
double               NEWTON_G;
int                  POT_GPU_NPGROUP;
bool                 OPT__OUTPUT_POT, OPT__GRA_P5_GRADIENT, OPT__EXTERNAL_POT;
double               SOR_OMEGA;
int                  SOR_MAX_ITER, SOR_MIN_ITER;
double               MG_TOLERATED_ERROR;
int                  MG_MAX_ITER, MG_NPRE_SMOOTH, MG_NPOST_SMOOTH;
IntScheme_t          OPT__POT_INT_SCHEME, OPT__RHO_INT_SCHEME, OPT__GRA_INT_SCHEME, OPT__REF_POT_INT_SCHEME;
OptPotBC_t           OPT__BC_POT;
OptGravityType_t     OPT__GRAVITY_TYPE;
#endif

// (2-3) cosmological simulations
#ifdef COMOVING
double               A_INIT, OMEGA_M0, DT__MAX_DELTA_A, HUBBLE0;
#endif

// (2-4) load balance
#ifdef LOAD_BALANCE
double               LB_INPUT__WLI_MAX;
#ifdef PARTICLE
double               LB_INPUT__PAR_WEIGHT;
#endif
bool                 OPT__RECORD_LOAD_BALANCE;
#endif
bool                 OPT__MINIMIZE_MPI_BARRIER;

// (2-5) particle
#ifdef PARTICLE
double               DT__PARVEL, DT__PARVEL_MAX, DT__PARACC;
bool                 OPT__OUTPUT_PAR_TEXT, OPT__CK_PARTICLE, OPT__FLAG_NPAR_CELL, OPT__FLAG_PAR_MASS_CELL;
int                  OPT__PARTICLE_COUNT, OPT__FLAG_NPAR_PATCH, FlagTable_NParPatch[NLEVEL-1], FlagTable_NParCell[NLEVEL-1];
double               FlagTable_ParMassCell[NLEVEL-1];
ParOutputDens_t      OPT__OUTPUT_PAR_DENS;
char                *PassiveFieldName_Par[PAR_NPASSIVE];
#endif

// (2-6) yt inline analysis
#ifdef SUPPORT_LIBYT
char                 YT_SCRIPT[MAX_STRING];
yt_verbose           YT_VERBOSE;
#endif

// (2-7) Grackle
#ifdef SUPPORT_GRACKLE
GrackleMode_t        GRACKLE_MODE;
bool                 GRACKLE_VERBOSE;
bool                 GRACKLE_COOLING;
GracklePriChe_t      GRACKLE_PRIMORDIAL;
bool                 GRACKLE_METAL;
bool                 GRACKLE_UV;
bool                 GRACKLE_CMB_FLOOR;
bool                 GRACKLE_PE_HEATING;
double               GRACKLE_PE_HEATING_RATE;
char                 GRACKLE_CLOUDY_TABLE[MAX_STRING];
int                  CHE_GPU_NPGROUP;
#endif

// (2-8) star formation
#ifdef STAR_FORMATION
SF_CreateStarScheme_t SF_CREATE_STAR_SCHEME;
int                   SF_CREATE_STAR_RSEED;
bool                  SF_CREATE_STAR_DET_RANDOM;
int                   SF_CREATE_STAR_MIN_LEVEL;
double                SF_CREATE_STAR_MIN_GAS_DENS;
double                SF_CREATE_STAR_MASS_EFF;
double                SF_CREATE_STAR_MIN_STAR_MASS;
double                SF_CREATE_STAR_MAX_STAR_MFRAC;
#endif


// 3. CPU (host) arrays for transferring data between CPU and GPU
// =======================================================================================================
// (3-1) fluid solver
real (*h_Flu_Array_F_In [2])[FLU_NIN ][  FLU_NXT   *FLU_NXT   *FLU_NXT   ] = { NULL, NULL };
real (*h_Flu_Array_F_Out[2])[FLU_NOUT][8*PATCH_SIZE*PATCH_SIZE*PATCH_SIZE] = { NULL, NULL };
real (*h_Flux_Array[2])[9][NFLUX_TOTAL][4*PATCH_SIZE*PATCH_SIZE]           = { NULL, NULL };
double (*h_Corner_Array_F[2])[3]                                           = { NULL, NULL };
#ifdef DUAL_ENERGY
char (*h_DE_Array_F_Out[2])[8*PATCH_SIZE*PATCH_SIZE*PATCH_SIZE]            = { NULL, NULL };
#endif

#ifdef GRAVITY
// (3-2) gravity solver
real (*h_Rho_Array_P    [2])[RHO_NXT][RHO_NXT][RHO_NXT]                    = { NULL, NULL };
real (*h_Pot_Array_P_In [2])[POT_NXT][POT_NXT][POT_NXT]                    = { NULL, NULL };
real (*h_Pot_Array_P_Out[2])[GRA_NXT][GRA_NXT][GRA_NXT]                    = { NULL, NULL };
real (*h_Flu_Array_G    [2])[GRA_NIN][PATCH_SIZE][PATCH_SIZE][PATCH_SIZE]  = { NULL, NULL };
double (*h_Corner_Array_G [2])[3]                                          = { NULL, NULL };
#ifdef DUAL_ENERGY
char (*h_DE_Array_G    [2])[PATCH_SIZE][PATCH_SIZE][PATCH_SIZE]            = { NULL, NULL };
#endif

// (3-3) unsplit gravity correction
#ifdef UNSPLIT_GRAVITY
real (*h_Pot_Array_USG_F[2])[USG_NXT_F][USG_NXT_F][USG_NXT_F]              = { NULL, NULL };
real (*h_Pot_Array_USG_G[2])[USG_NXT_G][USG_NXT_G][USG_NXT_G]              = { NULL, NULL };
real (*h_Flu_Array_USG_G[2])[GRA_NIN-1][PS1][PS1][PS1]                     = { NULL, NULL };
#endif
#endif

// (3-4) Grackle chemistry
#ifdef SUPPORT_GRACKLE
real (*h_Che_Array[2])                                                     = { NULL, NULL };
grackle_field_data *Che_FieldData                                          = NULL;
code_units Che_Units;
#endif

// (3-5) dt solver
real  *h_dt_Array_T[2]                                                     = { NULL, NULL };
real (*h_Flu_Array_T[2])[NCOMP_FLUID][ CUBE(PS1) ]                         = { NULL, NULL };
#ifdef GRAVITY
real (*h_Pot_Array_T[2])[ CUBE(GRA_NXT) ]                                  = { NULL, NULL };
#endif


// 4. GPU (device) global memory arrays
// =======================================================================================================
#ifdef GPU
// (4-1) fluid solver
real (*d_Flu_Array_F_In )[FLU_NIN ][ FLU_NXT*FLU_NXT*FLU_NXT ]             = NULL;
real (*d_Flu_Array_F_Out)[FLU_NOUT][ PS2*PS2*PS2 ]                         = NULL;
real (*d_Flux_Array)[9][NFLUX_TOTAL][ PS2*PS2 ]                            = NULL;
double (*d_Corner_Array_F)[3]                                              = NULL;
#ifdef DUAL_ENERGY
char (*d_DE_Array_F_Out)[ PS2*PS2*PS2 ]                                    = NULL;
#endif
#if ( MODEL == HYDRO )
#if ( FLU_SCHEME == MHM  ||  FLU_SCHEME == MHM_RP  ||  FLU_SCHEME == CTU )
real (*d_PriVar)     [NCOMP_TOTAL][ FLU_NXT*FLU_NXT*FLU_NXT ]              = NULL;
real (*d_Slope_PPM_x)[NCOMP_TOTAL][ N_SLOPE_PPM*N_SLOPE_PPM*N_SLOPE_PPM ]  = NULL;
real (*d_Slope_PPM_y)[NCOMP_TOTAL][ N_SLOPE_PPM*N_SLOPE_PPM*N_SLOPE_PPM ]  = NULL;
real (*d_Slope_PPM_z)[NCOMP_TOTAL][ N_SLOPE_PPM*N_SLOPE_PPM*N_SLOPE_PPM ]  = NULL;
real (*d_FC_Var_xL)  [NCOMP_TOTAL][ N_FC_VAR*N_FC_VAR*N_FC_VAR ]           = NULL;
real (*d_FC_Var_xR)  [NCOMP_TOTAL][ N_FC_VAR*N_FC_VAR*N_FC_VAR ]           = NULL;
real (*d_FC_Var_yL)  [NCOMP_TOTAL][ N_FC_VAR*N_FC_VAR*N_FC_VAR ]           = NULL;
real (*d_FC_Var_yR)  [NCOMP_TOTAL][ N_FC_VAR*N_FC_VAR*N_FC_VAR ]           = NULL;
real (*d_FC_Var_zL)  [NCOMP_TOTAL][ N_FC_VAR*N_FC_VAR*N_FC_VAR ]           = NULL;
real (*d_FC_Var_zR)  [NCOMP_TOTAL][ N_FC_VAR*N_FC_VAR*N_FC_VAR ]           = NULL;
real (*d_FC_Flux_x)  [NCOMP_TOTAL][ N_FC_FLUX*N_FC_FLUX*N_FC_FLUX ]        = NULL;
real (*d_FC_Flux_y)  [NCOMP_TOTAL][ N_FC_FLUX*N_FC_FLUX*N_FC_FLUX ]        = NULL;
real (*d_FC_Flux_z)  [NCOMP_TOTAL][ N_FC_FLUX*N_FC_FLUX*N_FC_FLUX ]        = NULL;
#endif // FLU_SCHEME
#elif ( MODEL == MHD )
#warning : WAIT MHD !!!
#endif // MODEL

#ifdef GRAVITY
// (4-2) gravity solver
real (*d_Rho_Array_P    )[ RHO_NXT*RHO_NXT*RHO_NXT ]                       = NULL;
real (*d_Pot_Array_P_In )[ POT_NXT*POT_NXT*POT_NXT ]                       = NULL;
real (*d_Pot_Array_P_Out)[ GRA_NXT*GRA_NXT*GRA_NXT ]                       = NULL;
real (*d_Flu_Array_G    )[GRA_NIN][ PS1*PS1*PS1 ]                          = NULL;
double (*d_Corner_Array_G )[3]                                             = NULL;
#ifdef DUAL_ENERGY
char (*d_DE_Array_G     )[ PS1*PS1*PS1 ]                                   = NULL;
#endif

// (4-3) unsplit gravity correction
#ifdef UNSPLIT_GRAVITY
real (*d_Pot_Array_USG_F)[ USG_NXT_F*USG_NXT_F*USG_NXT_F ]                 = NULL;
real (*d_Pot_Array_USG_G)[ USG_NXT_G*USG_NXT_G*USG_NXT_G ]                 = NULL;
real (*d_Flu_Array_USG_G)[GRA_NIN-1][ PS1*PS1*PS1        ]                 = NULL;
#endif
#endif

// (4-4) Grackle chemistry

// (4-5) dt solver
real *d_dt_Array_T                                                         = NULL;
real (*d_Flu_Array_T)[NCOMP_FLUID][ CUBE(PS1) ]                            = NULL;
#ifdef GRAVITY
real (*d_Pot_Array_T)[ CUBE(GRA_NXT) ]                                     = NULL;
#endif
#endif // #ifdef GPU


// 5. timers
// =======================================================================================================
#ifdef TIMING
Timer_t *Timer_Main[7];
Timer_t *Timer_MPI[3];
Timer_t *Timer_dt         [NLEVEL];
Timer_t *Timer_Flu_Advance[NLEVEL];
Timer_t *Timer_Gra_Advance[NLEVEL];
Timer_t *Timer_Che_Advance[NLEVEL];
Timer_t *Timer_SF         [NLEVEL];
Timer_t *Timer_FixUp      [NLEVEL];
Timer_t *Timer_Flag       [NLEVEL];
Timer_t *Timer_Refine     [NLEVEL];
Timer_t *Timer_GetBuf     [NLEVEL][9];
Timer_t *Timer_Lv         [NLEVEL];
Timer_t *Timer_Par_Update [NLEVEL][3];
Timer_t *Timer_Par_2Sib   [NLEVEL];
Timer_t *Timer_Par_2Son   [NLEVEL];
Timer_t *Timer_Par_Collect[NLEVEL];
Timer_t *Timer_Par_MPI    [NLEVEL][6];
#endif

#ifdef TIMING_SOLVER
Timer_t *Timer_Pre         [NLEVEL][NSOLVER];
Timer_t *Timer_Sol         [NLEVEL][NSOLVER];
Timer_t *Timer_Clo         [NLEVEL][NSOLVER];
Timer_t *Timer_Poi_PreRho  [NLEVEL];
Timer_t *Timer_Poi_PreFlu  [NLEVEL];
Timer_t *Timer_Poi_PrePot_C[NLEVEL];
Timer_t *Timer_Poi_PrePot_F[NLEVEL];
#endif


// function pointer for recording the user-specified info
extern void (*Aux_Record_User_Ptr)();




//-------------------------------------------------------------------------------------------------------
// Function    :  main
// Description :  GAMER main function
//-------------------------------------------------------------------------------------------------------
int main( int argc, char *argv[] )
{

// initialization
// ======================================================================================================
   Timer_t Timer_Total;
   Timer_Total.Start();

#  ifdef TIMING
   Timer_t  Timer_Init, Timer_Other;
   Timer_Init.Start();
#  endif


   Init_GAMER( &argc, &argv );

   Aux_TakeNote();

#  ifdef GPU
   CUAPI_DiagnoseDevice();
#  endif

   Output_DumpData( 0 );

   if ( OPT__PATCH_COUNT > 0 )            Aux_Record_PatchCount();
   if ( OPT__RECORD_MEMORY )              Aux_GetMemInfo();
   if ( OPT__RECORD_USER  &&
        Aux_Record_User_Ptr != NULL )     Aux_Record_User_Ptr();
#  ifdef PARTICLE
   if ( OPT__PARTICLE_COUNT > 0 )         Par_Aux_Record_ParticleCount();
#  endif

   Aux_Check();

#  ifdef TIMING
   Aux_ResetTimer();
#  endif


#  ifdef TIMING
   Timer_Init.Stop();
#  endif
// ======================================================================================================



// main loop
// ======================================================================================================
   MPI_Barrier( MPI_COMM_WORLD );


   while ( (Time[0]-END_T < -1.e-10)  &&  (Step < END_STEP) )
   {

#     ifdef TIMING
      MPI_Barrier( MPI_COMM_WORLD );
      Timer_Main[0]->Start();    // timer for one iteration
#     endif


//    1. advance all physical attributes by one global time-step
//    ---------------------------------------------------------------------------------------------------
      TIMING_FUNC(   EvolveLevel( 0, NULL_REAL ),   Timer_Main[2]   );

      Step ++;
//    ---------------------------------------------------------------------------------------------------


//    2. apply various corrections
//       --> synchronize particles, restrict data, recalculate potential and particle acceleration, ...
//    ---------------------------------------------------------------------------------------------------
      if ( OPT__CORR_AFTER_ALL_SYNC == CORR_AFTER_SYNC_EVERY_STEP )
      TIMING_FUNC(   Flu_CorrAfterAllSync(),     Timer_Main[6]   );
//    ---------------------------------------------------------------------------------------------------


//    3. output data and execute auxiliary functions
//    ---------------------------------------------------------------------------------------------------
      TIMING_FUNC(   Output_DumpData( 1 ),            Timer_Main[3]   );

      if ( OPT__PATCH_COUNT == 1 )
      TIMING_FUNC(   Aux_Record_PatchCount(),         Timer_Main[4]   );

      if ( OPT__RECORD_MEMORY )
      TIMING_FUNC(   Aux_GetMemInfo(),                Timer_Main[4]   );

      if ( OPT__RECORD_USER  &&  Aux_Record_User_Ptr != NULL )
      TIMING_FUNC(   Aux_Record_User_Ptr(),           Timer_Main[4]   );

      TIMING_FUNC(   Aux_Record_CorrUnphy(),          Timer_Main[4]   );

#     ifdef PARTICLE
      if ( OPT__PARTICLE_COUNT == 1 )
      TIMING_FUNC(   Par_Aux_Record_ParticleCount(),  Timer_Main[4]   );
#     endif

      TIMING_FUNC(   Aux_Check(),                     Timer_Main[4]   );
//    ---------------------------------------------------------------------------------------------------


//    4. perform yt inline analysis
//    ---------------------------------------------------------------------------------------------------
#     ifdef SUPPORT_LIBYT
      YT_Inline();
#     endif
//    ---------------------------------------------------------------------------------------------------


//    5. check whether to manually terminate the run
//    ---------------------------------------------------------------------------------------------------
      int Terminate = false;

//    enable this functionality only if OPT__MANUAL_CONTROL is on
      if ( OPT__MANUAL_CONTROL )
      TIMING_FUNC(   End_StopManually( Terminate ),   Timer_Main[4]   );
//    ---------------------------------------------------------------------------------------------------


//    6. check whether to redistribute all patches for LOAD_BALANCE
//    ---------------------------------------------------------------------------------------------------
#     ifdef LOAD_BALANCE
      if ( OPT__TIMING_BARRIER ) MPI_Barrier( MPI_COMM_WORLD );
#     ifdef TIMING
      Timer_Main[5]->Start();    // timer for load balance
#     endif

      if ( LB_EstimateLoadImbalance() > amr->LB->WLI_Max )
      {
         if ( MPI_Rank == 0 )
         {
            Aux_Message( stdout, "Weighted load-imbalance factor (%13.7e) > threshold (%13.7e) ",
                         amr->LB->WLI, amr->LB->WLI_Max );
            Aux_Message( stdout, "--> redistributing all patches ...\n" );
         }

         const bool   Redistribute_Yes = true;
         const bool   ResetLB_Yes      = true;
#        ifdef PARTICLE
         const double ParWeight        = amr->LB->Par_Weight;
#        else
         const double ParWeight        = 0.0;
#        endif
         LB_Init_LoadBalance( Redistribute_Yes, ParWeight, ResetLB_Yes );

         if ( OPT__PATCH_COUNT > 0 )         Aux_Record_PatchCount();

#        ifdef PARTICLE
         if ( OPT__PARTICLE_COUNT > 0 )      Par_Aux_Record_ParticleCount();
#        endif
      } // if ( LB_EstimateLoadImbalance() > amr->LB->WLI_Max )

#     ifdef TIMING
      Timer_Main[5]->Stop();
#     endif
#     endif // #ifdef LOAD_BALANCE
//    ---------------------------------------------------------------------------------------------------


//    7. record timing
//    ---------------------------------------------------------------------------------------------------
#     ifdef TIMING
      MPI_Barrier( MPI_COMM_WORLD );
      Timer_Main[0]->Stop();

      Timer_Other.Start();

      if ( OPT__RECORD_PERFORMANCE )
      Aux_Record_Performance( Timer_Main[0]->GetValue() );

      Aux_Record_Timing();

      Aux_ResetTimer();

      Timer_Other.Stop();
#     endif
//    ---------------------------------------------------------------------------------------------------


      if ( Terminate )  break;

   } // while ( (Time[0]-END_T < -1.e-10)  &&  (Step < END_STEP) )


   MPI_Barrier( MPI_COMM_WORLD );
   Timer_Total.Stop();
// ======================================================================================================



// termination
// ======================================================================================================
// output the final result
   Output_DumpData( 2 );


// record the total simulation time
#  ifdef TIMING
   Aux_AccumulatedTiming( Timer_Total.GetValue(), Timer_Init.GetValue(), Timer_Other.GetValue() );
#  endif

   if ( MPI_Rank == 0 )
   {
      FILE *Note = fopen( "Record__Note", "a" );
      fprintf( Note, "\n" );
      fprintf( Note, "Total Processing Time : %lf s\n", Timer_Total.GetValue() );
      fprintf( Note, "\n" );
      fclose( Note );
   }


   End_GAMER();
// ======================================================================================================

   return 0;

} // FUNCTION : Main

