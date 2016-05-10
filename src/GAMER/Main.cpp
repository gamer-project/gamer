#include "Copyright.h"
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
AMR_t            *amr = NULL;

double            Time[NLEVEL]           = { 0.0 };
long              AdvanceCounter[NLEVEL] = { 0 }; 
long              NCorrUnphy[NLEVEL]     = { 0 };
long              Step                   = 0;
int               DumpID                 = 0;
double            DumpTime               = 0.0;

double            dTime_Base;
double            Time_Prev            [NLEVEL];
real              MinDtInfo_Fluid      [NLEVEL];
double            FlagTable_Rho        [NLEVEL-1]; 
double            FlagTable_RhoGradient[NLEVEL-1]; 
double            FlagTable_Lohner     [NLEVEL-1][3];
double            FlagTable_User       [NLEVEL-1];
double           *DumpTable = NULL;
int               DumpTable_NDump;

int               MPI_Rank, MPI_Rank_X[3], MPI_SibRank[26], NX0[3], NPatchTotal[NLEVEL];
int              *BaseP = NULL; 
int               Flu_ParaBuf;

double            BOX_SIZE, DT__FLUID, DT__FLUID_INIT, END_T, OUTPUT_DT;
long              END_STEP;
int               NX0_TOT[3], OUTPUT_STEP, REGRID_COUNT, FLU_GPU_NPGROUP, OMP_NTHREAD;
int               MPI_NRank, MPI_NRank_X[3], GPU_NSTREAM, FLAG_BUFFER_SIZE, MAX_LEVEL;

IntScheme_t       OPT__FLU_INT_SCHEME, OPT__REF_FLU_INT_SCHEME;
double            OUTPUT_PART_X, OUTPUT_PART_Y, OUTPUT_PART_Z;
double            OPT__CK_MEMFREE, INT_MONO_COEFF;
int               OPT__UM_START_LEVEL, OPT__UM_START_NVAR, OPT__GPUID_SELECT, OPT__PATCH_COUNT;
int               INIT_DUMPID, INIT_SUBSAMPLING_NCELL;
bool              OPT__FLAG_RHO, OPT__FLAG_RHO_GRADIENT, OPT__FLAG_USER, OPT__FLAG_LOHNER_DENS, OPT__FLAG_REGION;
bool              OPT__DT_USER, OPT__RECORD_DT, OPT__RECORD_MEMORY, OPT__ADAPTIVE_DT;
bool              OPT__FIXUP_RESTRICT, OPT__INIT_RESTRICT, OPT__VERBOSE, OPT__MANUAL_CONTROL;
bool              OPT__INT_TIME, OPT__OUTPUT_TEST_ERROR, OPT__OUTPUT_BASE, OPT__OVERLAP_MPI, OPT__TIMING_BALANCE;
bool              OPT__OUTPUT_BASEPS, OPT__CK_REFINE, OPT__CK_PROPER_NESTING, OPT__CK_FINITE, OPT__RECORD_PERFORMANCE;
bool              OPT__CK_RESTRICT, OPT__CK_PATCH_ALLOCATE, OPT__FIXUP_FLUX, OPT__CK_FLUX_ALLOCATE;
bool              OPT__UM_START_DOWNGRADE, OPT__UM_START_REFINE, OPT__UM_FACTOR_5OVER3, OPT__TIMING_MPI;
bool              OPT__CK_CONSERVATION, OPT__RESET_FLUID, OPT__RECORD_USER, OPT__CORR_UNPHY;
OptInit_t         OPT__INIT;
OptRestartH_t     OPT__RESTART_HEADER;
OptOutputFormat_t OPT__OUTPUT_TOTAL;
OptOutputPart_t   OPT__OUTPUT_PART;
OptOutputMode_t   OPT__OUTPUT_MODE;
OptFluBC_t        OPT__BC_FLU[6];
OptLohnerForm_t   OPT__FLAG_LOHNER_FORM;


// 2. global variables for different applications
// =======================================================================================================
// (2-1) fluid solver in different models
#if   ( MODEL == HYDRO )
double         FlagTable_PresGradient[NLEVEL-1];
double         GAMMA, MINMOD_COEFF, EP_COEFF;
LR_Limiter_t   OPT__LR_LIMITER;
WAF_Limiter_t  OPT__WAF_LIMITER;
OptRSolver_t   OPT__CORR_UNPHY_SCHEME;
bool           OPT__FLAG_PRES_GRADIENT, OPT__FLAG_LOHNER_ENGY, OPT__FLAG_LOHNER_PRES;
int            OPT__CK_NEGATIVE;

#elif ( MODEL == MHD )
#warning : WAIT MHD !!!

#elif ( MODEL == ELBDM )
double         DT__PHASE, FlagTable_EngyDensity[NLEVEL-1][2];
bool           OPT__FLAG_ENGY_DENSITY, OPT__INT_PHASE;
bool           ELBDM_TAYLOR3_AUTO;
double         ELBDM_TAYLOR3_COEFF;
double         ELBDM_MASS, ELBDM_PLANCK_CONST, ELBDM_ETA;
#ifdef QUARTIC_SELF_INTERACTION
double         ELBDM_LAMBDA;
#endif
real           MinDtInfo_Phase[NLEVEL];

#else
#error : unsupported MODEL !!
#endif // MODEL

// (2-2) self-gravity
#ifdef GRAVITY
double           AveDensity = -1.0;    // initialize it as <= 0 to check if it is properly set later
real             MinDtInfo_Gravity[NLEVEL];
int              Pot_ParaBuf, Rho_ParaBuf;

real            *GreenFuncK      = NULL;
double           GFUNC_COEFF0;
double           DT__GRAVITY;
double           NEWTON_G;
int              POT_GPU_NPGROUP;
bool             OPT__OUTPUT_POT, OPT__GRA_P5_GRADIENT, OPT__EXTERNAL_POT;
double           SOR_OMEGA;
int              SOR_MAX_ITER, SOR_MIN_ITER;
double           MG_TOLERATED_ERROR;
int              MG_MAX_ITER, MG_NPRE_SMOOTH, MG_NPOST_SMOOTH;
IntScheme_t      OPT__POT_INT_SCHEME, OPT__RHO_INT_SCHEME, OPT__GRA_INT_SCHEME, OPT__REF_POT_INT_SCHEME;
OptPotBC_t       OPT__BC_POT;
OptGravityType_t OPT__GRAVITY_TYPE;
#endif

// (2-3) cosmological simulations
#ifdef COMOVING
double         A_INIT, OMEGA_M0, DT__MAX_DELTA_A;
#endif

// (2-4) load balance
#ifdef LOAD_BALANCE
double         LB_INPUT__WLI_MAX;
#endif

// (2-5) particle
#ifdef PARTICLE 
double         DT__PARVEL;
real           MinDtInfo_ParVel[NLEVEL];
bool           OPT__OUTPUT_PARTICLE, OPT__CK_PARTICLE, OPT__PAR_LEVEL;
#endif


// 3. CPU (host) arrays for transferring data bewteen CPU and GPU
// =======================================================================================================
// (3-1) fluid solver
real (*h_Flu_Array_F_In [2])[FLU_NIN ][  FLU_NXT   *FLU_NXT   *FLU_NXT   ] = { NULL, NULL };
real (*h_Flu_Array_F_Out[2])[FLU_NOUT][8*PATCH_SIZE*PATCH_SIZE*PATCH_SIZE] = { NULL, NULL };
real (*h_Flux_Array[2])[9][NFLUX][4*PATCH_SIZE*PATCH_SIZE]                 = { NULL, NULL };
double (*h_Corner_Array_F [2])[3]                                          = { NULL, NULL };
real *h_MinDtInfo_Fluid_Array[2]                                           = { NULL, NULL };

#ifdef GRAVITY
// (3-2) gravity solver
real (*h_Rho_Array_P    [2])[RHO_NXT][RHO_NXT][RHO_NXT]                    = { NULL, NULL };
real (*h_Pot_Array_P_In [2])[POT_NXT][POT_NXT][POT_NXT]                    = { NULL, NULL };
real (*h_Pot_Array_P_Out[2])[GRA_NXT][GRA_NXT][GRA_NXT]                    = { NULL, NULL };
real (*h_Flu_Array_G    [2])[GRA_NIN][PATCH_SIZE][PATCH_SIZE][PATCH_SIZE]  = { NULL, NULL };
double (*h_Corner_Array_G [2])[3]                                          = { NULL, NULL };

// (3-3) unsplit gravity correction
#ifdef UNSPLIT_GRAVITY
real (*h_Pot_Array_USG_F[2])[USG_NXT_F][USG_NXT_F][USG_NXT_F]              = { NULL, NULL };
real (*h_Pot_Array_USG_G[2])[USG_NXT_G][USG_NXT_G][USG_NXT_G]              = { NULL, NULL };
real (*h_Flu_Array_USG_G[2])[GRA_NIN-1][PS1][PS1][PS1]                     = { NULL, NULL };
#endif
#endif


// 4. GPU (device) global memory arrays
// =======================================================================================================
#ifdef GPU
// (4-1) fluid solver
real (*d_Flu_Array_F_In )[FLU_NIN ][ FLU_NXT*FLU_NXT*FLU_NXT ]    = NULL;
real (*d_Flu_Array_F_Out)[FLU_NOUT][ PS2*PS2*PS2 ]                = NULL;
real (*d_Flux_Array)[9][NFLUX][ PS2*PS2 ]                         = NULL;
double (*d_Corner_Array_F)[3]                                     = NULL;
real  *d_MinDtInfo_Fluid_Array                                    = NULL;
#if ( MODEL == HYDRO )
#if ( FLU_SCHEME == MHM  ||  FLU_SCHEME == MHM_RP  ||  FLU_SCHEME == CTU )
real (*d_PriVar)     [5][ FLU_NXT*FLU_NXT*FLU_NXT ]               = NULL;
real (*d_Slope_PPM_x)[5][ N_SLOPE_PPM*N_SLOPE_PPM*N_SLOPE_PPM ]   = NULL;
real (*d_Slope_PPM_y)[5][ N_SLOPE_PPM*N_SLOPE_PPM*N_SLOPE_PPM ]   = NULL;
real (*d_Slope_PPM_z)[5][ N_SLOPE_PPM*N_SLOPE_PPM*N_SLOPE_PPM ]   = NULL;
real (*d_FC_Var_xL)  [5][ N_FC_VAR*N_FC_VAR*N_FC_VAR ]            = NULL;
real (*d_FC_Var_xR)  [5][ N_FC_VAR*N_FC_VAR*N_FC_VAR ]            = NULL;
real (*d_FC_Var_yL)  [5][ N_FC_VAR*N_FC_VAR*N_FC_VAR ]            = NULL;
real (*d_FC_Var_yR)  [5][ N_FC_VAR*N_FC_VAR*N_FC_VAR ]            = NULL;
real (*d_FC_Var_zL)  [5][ N_FC_VAR*N_FC_VAR*N_FC_VAR ]            = NULL;
real (*d_FC_Var_zR)  [5][ N_FC_VAR*N_FC_VAR*N_FC_VAR ]            = NULL;
real (*d_FC_Flux_x)  [5][ N_FC_FLUX*N_FC_FLUX*N_FC_FLUX ]         = NULL;
real (*d_FC_Flux_y)  [5][ N_FC_FLUX*N_FC_FLUX*N_FC_FLUX ]         = NULL;
real (*d_FC_Flux_z)  [5][ N_FC_FLUX*N_FC_FLUX*N_FC_FLUX ]         = NULL;
#endif // FLU_SCHEME
#elif ( MODEL == MHD )
#warning : WAIT MHD !!!
#endif // MODEL

#ifdef GRAVITY
// (4-2) gravity solver
real (*d_Rho_Array_P    )[ RHO_NXT*RHO_NXT*RHO_NXT ]              = NULL;
real (*d_Pot_Array_P_In )[ POT_NXT*POT_NXT*POT_NXT ]              = NULL;
real (*d_Pot_Array_P_Out)[ GRA_NXT*GRA_NXT*GRA_NXT ]              = NULL;
real (*d_Flu_Array_G    )[GRA_NIN][ PS1*PS1*PS1 ]                 = NULL;
double (*d_Corner_Array_G )[3]                                    = NULL;

// (3-3) unsplit gravity correction
#ifdef UNSPLIT_GRAVITY
real (*d_Pot_Array_USG_F)[ USG_NXT_F*USG_NXT_F*USG_NXT_F ]        = NULL;
real (*d_Pot_Array_USG_G)[ USG_NXT_G*USG_NXT_G*USG_NXT_G ]        = NULL;
real (*d_Flu_Array_USG_G)[GRA_NIN-1][ PS1*PS1*PS1        ]        = NULL;
#endif
#endif

// (4-3) CUDA stream (put in CUAPI_MemAllocate_Fluid.cu)
//cudaStream_t *Stream                                              = NULL;
#endif // #ifdef GPU


// 5. timers
// =======================================================================================================
#ifdef TIMING
Timer_t *Timer_Main[6];
Timer_t *Timer_MPI[3];
Timer_t *Timer_Flu_Advance[NLEVEL];
Timer_t *Timer_Gra_Advance[NLEVEL];
Timer_t *Timer_FixUp      [NLEVEL];
Timer_t *Timer_Flag       [NLEVEL];
Timer_t *Timer_Refine     [NLEVEL];
Timer_t *Timer_GetBuf     [NLEVEL][8];
Timer_t *Timer_Lv         [NLEVEL];
#endif

#ifdef TIMING_SOLVER
Timer_t *Timer_Pre         [NLEVEL][4];
Timer_t *Timer_Sol         [NLEVEL][4];
Timer_t *Timer_Clo         [NLEVEL][4];
Timer_t *Timer_Poi_PreRho  [NLEVEL];
Timer_t *Timer_Poi_PreFlu  [NLEVEL];
Timer_t *Timer_Poi_PrePot_C[NLEVEL];
Timer_t *Timer_Poi_PrePot_F[NLEVEL];
#endif




//-------------------------------------------------------------------------------------------------------
// Function    :  main
// Description :  GAMER main function 
//-------------------------------------------------------------------------------------------------------
int main( int argc, char *argv[] )
{

// initialization 
// ======================================================================================================
   Timer_t  Timer_Total( 1 );
   Timer_Total.Start();

#  ifdef TIMING
   Timer_t  Timer_Init( 1 ), Timer_Other( 1 );
   Timer_Init.Start();
#  endif


   Init_GAMER( &argc, &argv );

   Aux_TakeNote();

#  ifdef GPU
   CUAPI_DiagnoseDevice();
#  endif

   Output_DumpData( 0 );

   if ( OPT__PATCH_COUNT > 0 )   Aux_PatchCount();
   if ( OPT__RECORD_MEMORY   )   Aux_GetMemInfo();
   if ( OPT__RECORD_USER     )   Aux_RecordUser();
#  ifdef PARTICLE
   if ( OPT__PAR_LEVEL       )   Par_Aux_GetParticleLevel();
#  endif

   Aux_Check();

#  ifdef TIMING
   Aux_ResetTimer();
#  endif


#  ifdef TIMING
   Timer_Init.Stop( false );
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


//    a. determine the time step
//    ---------------------------------------------------------------------------------------------------
      TIMING_FUNC(   Mis_GetTimeStep(),   Timer_Main[1],   false   );

      if ( MPI_Rank == 0  &&  Step%1 == 0 )   
      {
         Aux_Message( stdout, "Time : %13.7e -> %13.7e,    Step : %7ld -> %7ld,    dt = %14.7e\n", 
                      Time[0], Time[0]+dTime_Base, Step, Step+1, dTime_Base );
      }

//    ---------------------------------------------------------------------------------------------------


//    b. advance all physical attributes by one step
//    ---------------------------------------------------------------------------------------------------
      TIMING_FUNC(   EvolveLevel( 0, dTime_Base ),   Timer_Main[2],   false   );

      Step ++;
//    ---------------------------------------------------------------------------------------------------


//    c. restrict all data and re-calculate potential in the debug mode (in order to check the RESTART process)
//    ---------------------------------------------------------------------------------------------------
#     ifdef GAMER_DEBUG
      for (int lv=MAX_LEVEL-1; lv>=0; lv--)
      {
         if ( NPatchTotal[lv+1] == 0 )    continue;

         if ( OPT__VERBOSE  &&  MPI_Rank == 0 )
            Aux_Message( stdout, "   DEBUG: restrict data at Lv %2d            ... ", lv );

//       to ensure the same round-off errors for the father patches with newly refined son patches in the RESTART process
         Flu_Restrict( lv, amr->FluSg[lv+1], amr->FluSg[lv], NULL_INT, NULL_INT, _FLU );

#        ifdef LOAD_BALANCE 
         LB_GetBufferData( lv, amr->FluSg[lv], NULL_INT, DATA_RESTRICT, _FLU, NULL_INT );
#        endif

         Buf_GetBufferData( lv, amr->FluSg[lv], NULL_INT, DATA_GENERAL, _FLU, Flu_ParaBuf, USELB_YES );

         if ( OPT__VERBOSE  &&  MPI_Rank == 0 )    Aux_Message( stdout, "done\n" );
      }

#     ifdef GRAVITY
      if ( OPT__GRAVITY_TYPE == GRAVITY_SELF  ||  OPT__GRAVITY_TYPE == GRAVITY_BOTH )
      for (int lv=0; lv<=MAX_LEVEL; lv++)
      {
         if ( NPatchTotal[lv] == 0 )   break;

         if ( OPT__VERBOSE  &&  MPI_Rank == 0 )
            Aux_Message( stdout, "   DEBUG: recalculate potential at Lv %2d    ... ", lv );

#        ifdef COMOVING
         const double Poi_Coeff = 4.0*M_PI*NEWTON_G*Time[lv];
#        else
         const double Poi_Coeff = 4.0*M_PI*NEWTON_G;
#        endif

//       collect particles from all descendant patches
#        ifdef PARTICLE
         Par_CollectParticleForPoisson( lv );
#        endif

         if ( lv == 0 )    
            CPU_PoissonSolver_FFT( Poi_Coeff, amr->PotSg[lv], Time[lv] );

         else              
         {
            Buf_GetBufferData( lv, amr->FluSg[lv], NULL_INT, DATA_GENERAL, _DENS, Rho_ParaBuf, USELB_YES );

            InvokeSolver( POISSON_SOLVER, lv, Time[lv], NULL_REAL, NULL_REAL, Poi_Coeff, NULL_INT, amr->PotSg[lv], false, false );
         }

         Buf_GetBufferData( lv, NULL_INT, amr->PotSg[lv], POT_FOR_POISSON, _POTE, Pot_ParaBuf, USELB_YES );

//       free variables of descendant particles
#        ifdef PARTICLE
         Par_ReleaseParticleForPoisson( lv );
#        endif

         if ( OPT__VERBOSE  &&  MPI_Rank == 0 )    Aux_Message( stdout, "done\n" );

      } // for (int lv=0; lv<=MAX_LEVEL; lv++)
#     endif // #ifdef GRAVITY
#     endif // #ifdef GAMER_DEBUG


//    d. output data and execute auxiliary functions
//    ---------------------------------------------------------------------------------------------------
      TIMING_FUNC(   Output_DumpData( 1 ),       Timer_Main[3],   false   );

      if ( OPT__PATCH_COUNT == 1 )
      TIMING_FUNC(   Aux_PatchCount(),           Timer_Main[4],   false   );

      if ( OPT__RECORD_MEMORY )   
      TIMING_FUNC(   Aux_GetMemInfo(),           Timer_Main[4],   false   );

      if ( OPT__RECORD_USER     )
      TIMING_FUNC(   Aux_RecordUser(),           Timer_Main[4],   false   );

      if ( OPT__CORR_UNPHY )
      TIMING_FUNC(   Aux_RecordCorrUnphy(),      Timer_Main[4],   false   );

#     ifdef PARTICLE
      if ( OPT__PAR_LEVEL )
      TIMING_FUNC(   Par_Aux_GetParticleLevel(), Timer_Main[4],   false   );
#     endif

      TIMING_FUNC(   Aux_Check(),                Timer_Main[4],   false   );
//    ---------------------------------------------------------------------------------------------------


//    e. check whether to manually terminate the run 
//    ---------------------------------------------------------------------------------------------------
      int Terminate = false;

//    enable this functionality only if OPT__MANUAL_CONTROL is on
      if ( OPT__MANUAL_CONTROL )
      TIMING_FUNC(   End_StopManually( Terminate ),   Timer_Main[4],   false   );
//    ---------------------------------------------------------------------------------------------------


//    f. check whether to redistribute all patches for LOAD_BALANCE
//    ---------------------------------------------------------------------------------------------------
#     ifdef LOAD_BALANCE
      const bool DuringRestart_No = false;

      if ( amr->LB->WLI > amr->LB->WLI_Max )
      {
         if ( MPI_Rank == 0 )    
         {
            Aux_Message( stdout, "Weighted load-imbalance factor (%13.7e) > threshold (%13.7e) ", 
                         amr->LB->WLI, amr->LB->WLI_Max );
            Aux_Message( stdout, "--> redistributing all patches ...\n" );
         }

         TIMING_FUNC(   amr->LB->reset(),                        Timer_Main[5],   false   );
         TIMING_FUNC(   LB_Init_LoadBalance( DuringRestart_No ), Timer_Main[5],   false   );
         TIMING_FUNC(   Aux_PatchCount(),                        Timer_Main[5],   false   );
      }
#     endif
//    ---------------------------------------------------------------------------------------------------


//    g. record timing
//    ---------------------------------------------------------------------------------------------------
#     ifdef TIMING
      MPI_Barrier( MPI_COMM_WORLD );
      Timer_Main[0]->Stop( false );

      Timer_Other.Start();

      if ( OPT__RECORD_PERFORMANCE )
      Aux_RecordPerformance( Timer_Main[0]->GetValue(0) );

      Aux_RecordTiming();

      Aux_ResetTimer();

      Timer_Other.Stop( false );
#     endif
//    ---------------------------------------------------------------------------------------------------


      if ( Terminate )  break;

   } // while ( (Time[0]-END_T < -1.e-10)  &&  (Step < END_STEP) )


   MPI_Barrier( MPI_COMM_WORLD );
   Timer_Total.Stop( false );
// ======================================================================================================



// termination
// ======================================================================================================
// output the final result
   Output_DumpData( 2 );


// record the total simulation time
#  ifdef TIMING
   Aux_AccumulatedTiming( Timer_Total.GetValue(0), Timer_Init.GetValue(0), Timer_Other.GetValue(0) );
#  endif

   if ( MPI_Rank == 0 )
   {
      FILE *Note = fopen( "Record__Note", "a" );
      fprintf( Note, "\n" );
      fprintf( Note, "Total Processing Time : %lf s\n", Timer_Total.GetValue( 0 ) );
      fprintf( Note, "\n" );
      fclose( Note );
   }


   End_GAMER();
// ======================================================================================================

   return 0;

} // FUNCTION : Main

