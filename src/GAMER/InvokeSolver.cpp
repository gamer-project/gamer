#include "Copyright.h"
#include "GAMER.h"

static void Preparation_Step( const Solver_t TSolver, const int lv, const double TimeNew, const double TimeOld, const int NPG, 
                              const int *PID0_List, const int ArrayID );
static void Solver( const Solver_t TSolver, const int lv, const double TimeNew, const double TimeOld,
                    const int NPG, const int ArrayID, const double dt, const double Poi_Coeff );
static void Closing_Step( const Solver_t TSolver, const int lv, const int SaveSg_Flu, const int SaveSg_Pot, const int NPG,
                          const int *PID0_List, const int ArrayID, const double dt );

extern Timer_t *Timer_Pre         [NLEVEL][4];
extern Timer_t *Timer_Sol         [NLEVEL][4];
extern Timer_t *Timer_Clo         [NLEVEL][4];
#ifdef GRAVITY
extern Timer_t *Timer_Poi_PreRho  [NLEVEL];
extern Timer_t *Timer_Poi_PreFlu  [NLEVEL];
extern Timer_t *Timer_Poi_PrePot_C[NLEVEL];
extern Timer_t *Timer_Poi_PrePot_F[NLEVEL];
#endif




//-------------------------------------------------------------------------------------------------------
// Function    :  InvokeSolver
// Description :  Invoke the GPU (or CPU) solvers and enable the concurrent execution between CPU and GPU
//
// Note        :  a. Use the input parameter "TSolver" to control the targeted solver
//                b. Each solver involves three steps
//                   --> 1. preparation step : prepare the input data 
//                       2. execution   step : invoke the solvers --> advance solutions or evaluate potentials
//                       3. closing     step : store the updated data
//                c. Currently the fluid solver can only store the updated data in the different sandglass from 
//                   the input data
//                d. For LOAD_BALANCE, one can turn on the option "LB_INPUT__OVERLAP_MPI" to enable the 
//                   overlapping between MPI communication and CPU/GPU computation
//
// Parameter   :  TSolver        : Targeted solver
//                                 --> FLUID_SOLVER               : Fluid / ELBDM solver
//                                     POISSON_SOLVER             : Poisson solver
//                                     GRAVITY_SOLVER             : Gravity solver
//                                     POISSON_AND_GRAVITY_SOLVER : Poisson + Gravity solvers
//                lv             : Targeted refinement level 
//                TimeNew        : Targeted physical time to reach
//                TimeOld        : Physical time before update
//                                 --> For Fluid and Gravity solver, this function updates physical time from TimeOld to TimeNew
//                                     For Poisson solver, this function calculates potential at **TimeNew**
//                dt             : Time interval to advance solution for the fluid and gravity solvers
//                                 (can be different from TimeNew-TimeOld if COMOVING is on)
//                Poi_Coeff      : Coefficient in front of the RHS in the Poisson eq.
//                SaveSg_Flu     : Sandglass to store the updated fluid data (for both the fluid and gravity solvers)
//                SaveSg_Pot     : Sandglass to store the updated potential data (for the Poisson solver)
//                OverlapMPI     : true --> Overlap MPI time with CPU/GPU computation
//                Overlap_Sync   : true  --> Advance the patches which cannot be overlapped with MPI communication
//                                 false --> Advance the patches which can    be overlapped with MPI communication
//                                 (useful only if "OverlapMPI == true")
//-------------------------------------------------------------------------------------------------------
void InvokeSolver( const Solver_t TSolver, const int lv, const double TimeNew, const double TimeOld, const double dt,
                   const double Poi_Coeff, const int SaveSg_Flu, const int SaveSg_Pot, 
                   const bool OverlapMPI, const bool Overlap_Sync )
{

// check
// currently the fluid solver can only store the updated data in the different sandglass from the input data
   if ( TSolver == FLUID_SOLVER  &&  SaveSg_Flu == amr->FluSg[lv] )
      Aux_Error( ERROR_INFO, "SaveSg_Flu (%d) == amr->FluSg (%d) in the fluid solver at level %d !!\n", 
                 SaveSg_Flu, amr->FluSg[lv], lv );

#  ifndef GRAVITY
   if ( TSolver != FLUID_SOLVER )
      Aux_Error( ERROR_INFO, "incorrect parameter %s = %d !!\n", "TSolver", TSolver );
#  endif

   if ( TSolver == FLUID_SOLVER  &&  ( SaveSg_Flu != 0 &&  SaveSg_Flu != 1 )  )
      Aux_Error( ERROR_INFO, "incorrect SaveSg_Flu (%d) !!\n", SaveSg_Flu );

#  ifdef GRAVITY
   if (  ( TSolver == GRAVITY_SOLVER || TSolver == POISSON_AND_GRAVITY_SOLVER ) && ( SaveSg_Flu != 0 &&  SaveSg_Flu != 1 )  )
      Aux_Error( ERROR_INFO, "incorrect SaveSg_Flu (%d) !!\n", SaveSg_Flu );

   if (  ( TSolver == POISSON_SOLVER || TSolver == POISSON_AND_GRAVITY_SOLVER ) && ( SaveSg_Pot != 0 &&  SaveSg_Pot != 1 )  )
      Aux_Error( ERROR_INFO, "incorrect SaveSg_Pot (%d) !!\n", SaveSg_Pot );
#  endif


// maximum number of patch groups to be updated at a time
#  ifdef GRAVITY
   const int NPG_Max = ( TSolver == FLUID_SOLVER ) ? FLU_GPU_NPGROUP : POT_GPU_NPGROUP;
#  else
   const int NPG_Max = FLU_GPU_NPGROUP;
#  endif

   int *PID0_List    = NULL;  // list recording the patch indicies with LocalID==0 to be udpated
   bool AllocateList = false; // whether to allocate PID0_List or not
   int  ArrayID      = 0;     // array index to load and store data ( 0 or 1 )
   int  NPG[2];               // number of patch groups to be updated at a time 
   int  NTotal;               // total number of patch groups to be updated 
   int  Disp;                 // index displacement in PID0_List
 
   if ( OverlapMPI )
   {
#     ifdef LOAD_BALANCE
      if ( TSolver == FLUID_SOLVER )  
      {
         if ( Overlap_Sync )
         {
            NTotal    = amr->LB->OverlapMPI_FluSyncN   [lv];
            PID0_List = amr->LB->OverlapMPI_FluSyncPID0[lv]; 
         }
         else
         {
            NTotal    = amr->LB->OverlapMPI_FluAsyncN   [lv];
            PID0_List = amr->LB->OverlapMPI_FluAsyncPID0[lv]; 
         }
      }

#     ifdef GRAVITY
      else              
      {
         if ( Overlap_Sync )
         {
            NTotal    = amr->LB->OverlapMPI_PotSyncN   [lv];
            PID0_List = amr->LB->OverlapMPI_PotSyncPID0[lv]; 
         }
         else
         {
            NTotal    = amr->LB->OverlapMPI_PotAsyncN   [lv];
            PID0_List = amr->LB->OverlapMPI_PotAsyncPID0[lv]; 
         }
      }
#     endif 

#     else // #ifdef LOAD_BALANCE ... else ...
      Aux_Error( ERROR_INFO, "MPI overlapping is NOT supported if LOAD_BALANCE is off !!\n" );
#     endif // #ifdef LOAD_BALANCE ... else ...

   } // if ( OverlapMPI )

   else
   {
      AllocateList = true;
      NTotal       = amr->NPatchComma[lv][1] / 8;
      PID0_List    = new int [NTotal];

      for (int t=0; t<NTotal; t++)  PID0_List[t] = 8*t;
   } // if ( OverlapMPI ) ... else ...

   NPG[ArrayID] = ( NPG_Max < NTotal ) ? NPG_Max : NTotal;


//-------------------------------------------------------------------------------------------------------------
   TIMING_SYNC(   Preparation_Step( TSolver, lv, TimeNew, TimeOld, NPG[ArrayID], PID0_List, ArrayID ),
                  Timer_Pre[lv][TSolver]  );
//-------------------------------------------------------------------------------------------------------------


//-------------------------------------------------------------------------------------------------------------
   TIMING_SYNC(   Solver( TSolver, lv, TimeNew, TimeOld, NPG[ArrayID], ArrayID, dt, Poi_Coeff ),
                  Timer_Sol[lv][TSolver]  );
//-------------------------------------------------------------------------------------------------------------


   for (Disp=NPG_Max; Disp<NTotal; Disp+=NPG_Max)
   {

      ArrayID      = 1 - ArrayID;
      NPG[ArrayID] = ( NPG_Max < NTotal-Disp ) ? NPG_Max : NTotal-Disp;


//-------------------------------------------------------------------------------------------------------------
      TIMING_SYNC(   Preparation_Step( TSolver, lv, TimeNew, TimeOld, NPG[ArrayID], PID0_List+Disp, ArrayID ),  
                     Timer_Pre[lv][TSolver]  );
//-------------------------------------------------------------------------------------------------------------


//-------------------------------------------------------------------------------------------------------------
#     ifdef GPU
      CUAPI_Synchronize();
#     endif
//-------------------------------------------------------------------------------------------------------------


//-------------------------------------------------------------------------------------------------------------
      TIMING_SYNC(   Solver( TSolver, lv, TimeNew, TimeOld, NPG[ArrayID], ArrayID, dt, Poi_Coeff ), 
                     Timer_Sol[lv][TSolver]  );
//-------------------------------------------------------------------------------------------------------------


//-------------------------------------------------------------------------------------------------------------
      TIMING_SYNC(   Closing_Step( TSolver, lv, SaveSg_Flu, SaveSg_Pot, NPG[1-ArrayID], PID0_List+Disp-NPG_Max, 1-ArrayID, dt ),
                     Timer_Clo[lv][TSolver]  ); 
//-------------------------------------------------------------------------------------------------------------

   } // for (int Disp=NPG_Max; Disp<NTotal; Disp+=NPG_Max)


//-------------------------------------------------------------------------------------------------------------
#  ifdef GPU
   CUAPI_Synchronize();
#  endif
//-------------------------------------------------------------------------------------------------------------


//-------------------------------------------------------------------------------------------------------------
   TIMING_SYNC(   Closing_Step( TSolver, lv, SaveSg_Flu, SaveSg_Pot, NPG[ArrayID], PID0_List+Disp-NPG_Max, ArrayID, dt ),
                  Timer_Clo[lv][TSolver]  ); 
//-------------------------------------------------------------------------------------------------------------
     

   if ( AllocateList )  delete [] PID0_List;

} // FUNCTION : InvokeSolver



//-------------------------------------------------------------------------------------------------------
// Function    :  Preparation_Step
// Description :  Prepare the input data for CPU/GPU solvers 
//
// Note        :  Use the input parameter "TSolver" to control the targeted solver 
//
// Parameter   :  TSolver     : Targeted solver
//                              --> FLUID_SOLVER               : Fluid / ELBDM solver
//                                  POISSON_SOLVER             : Poisson solver
//                                  GRAVITY_SOLVER             : Gravity solver
//                                  POISSON_AND_GRAVITY_SOLVER : Poisson + Gravity solvers
//                lv          : Targeted refinement level 
//                TimeNew     : Targeted physical time to reach
//                TimeOld     : Physical time before update
//                              --> For Fluid   solver, it prepares data at TimeOld
//                              --> For Gravity solver, it prepares data at TimeNew
//                                  (TimeOld data will also prepared for UNSPLIT_GRAVITY)
//                                  For Poisson solver, it prepares data at TimeNew
//                NPG         : Number of patch groups to be prepared at a time
//                PID0_List   : List recording the patch indicies with LocalID==0 to be udpated
//                ArrayID     : Array index to load and store data ( 0 or 1 )
//-------------------------------------------------------------------------------------------------------
void Preparation_Step( const Solver_t TSolver, const int lv, const double TimeNew, const double TimeOld, const int NPG, 
                       const int *PID0_List, const int ArrayID )
{

#  ifndef UNSPLIT_GRAVITY
   real (*h_Pot_Array_USG_F[2])[USG_NXT_F][USG_NXT_F][USG_NXT_F] = { NULL, NULL };
#  endif


   switch ( TSolver )
   {
      case FLUID_SOLVER :  
         Flu_Prepare( lv, TimeOld, h_Flu_Array_F_In[ArrayID][0][0], h_Pot_Array_USG_F[ArrayID][0][0][0],
                      h_Corner_Array_F[ArrayID], NPG, PID0_List );
         break;

#     ifdef GRAVITY
      case POISSON_SOLVER :  
         TIMING_SYNC(   Poi_Prepare_Rho( lv, TimeNew, h_Rho_Array_P   [ArrayID], NPG, PID0_List ), 
                        Timer_Poi_PreRho[lv]   );

         TIMING_SYNC(   Poi_Prepare_Pot( lv, TimeNew, h_Pot_Array_P_In[ArrayID], NPG, PID0_List ),
                        Timer_Poi_PrePot_C[lv]   );
         break;

      case GRAVITY_SOLVER :  
         TIMING_SYNC(   Gra_Prepare_Flu( lv,          h_Flu_Array_G   [ArrayID], NPG, PID0_List ),
                        Timer_Poi_PreFlu[lv]   );

         if ( OPT__GRAVITY_TYPE == GRAVITY_SELF  ||  OPT__GRAVITY_TYPE == GRAVITY_BOTH )
         TIMING_SYNC(   Gra_Prepare_Pot( lv, TimeNew, h_Pot_Array_P_Out[ArrayID], NPG, PID0_List ),
                        Timer_Poi_PrePot_F[lv]   );

         if ( OPT__GRAVITY_TYPE == GRAVITY_EXTERNAL  ||  OPT__GRAVITY_TYPE == GRAVITY_BOTH  ||  OPT__EXTERNAL_POT )
         TIMING_SYNC(   Gra_Prepare_Corner( lv, h_Corner_Array_G[ArrayID], NPG, PID0_List ),
                        Timer_Poi_PreFlu[lv]   );

#        ifdef UNSPLIT_GRAVITY
//       use the same timer "Timer_Poi_PreFlu" as Gra_Prepare_Flu and Gra_Prepare_Corner
         TIMING_SYNC(   Gra_Prepare_USG( lv, TimeOld, h_Pot_Array_USG_G[ArrayID], h_Flu_Array_USG_G[ArrayID],
                        NPG, PID0_List ),
                        Timer_Poi_PreFlu[lv]   );
#        endif
         break;

      case POISSON_AND_GRAVITY_SOLVER :  
         TIMING_SYNC(   Poi_Prepare_Rho( lv, TimeNew, h_Rho_Array_P   [ArrayID], NPG, PID0_List ), 
                        Timer_Poi_PreRho[lv]   );

         TIMING_SYNC(   Poi_Prepare_Pot( lv, TimeNew, h_Pot_Array_P_In[ArrayID], NPG, PID0_List ),
                        Timer_Poi_PrePot_C[lv]   );

         TIMING_SYNC(   Gra_Prepare_Flu( lv,          h_Flu_Array_G   [ArrayID], NPG, PID0_List ),
                        Timer_Poi_PreFlu[lv]   );

         if ( OPT__GRAVITY_TYPE == GRAVITY_EXTERNAL  ||  OPT__GRAVITY_TYPE == GRAVITY_BOTH  ||  OPT__EXTERNAL_POT )
         TIMING_SYNC(   Gra_Prepare_Corner( lv, h_Corner_Array_G[ArrayID], NPG, PID0_List ),
                        Timer_Poi_PreFlu[lv]   );

#        ifdef UNSPLIT_GRAVITY
//       use the same timer "Timer_Poi_PreFlu" as Gra_Prepare_Flu and Gra_Prepare_Corner
         TIMING_SYNC(   Gra_Prepare_USG( lv, TimeOld, h_Pot_Array_USG_G[ArrayID], h_Flu_Array_USG_G[ArrayID],
                        NPG, PID0_List ),
                        Timer_Poi_PreFlu[lv]   );
#        endif
         break;
#     endif // #ifdef GARVITY

      default :
         Aux_Error( ERROR_INFO, "incorrect parameter %s = %d !!\n", "TSolver", TSolver );

   } // switch ( TSolver )

} // FUNCTION : Preparation_Step



//-------------------------------------------------------------------------------------------------------
// Function    :  Solver
// Description :  Invoke the CPU/GPU solvers
//
// Note        :  Use the input parameter "TSolver" to control the targeted solver 
//
// Parameter   :  TSolver     : Targeted solver
//                              --> FLUID_SOLVER               : Fluid / ELBDM solver
//                                  POISSON_SOLVER             : Poisson solver
//                                  GRAVITY_SOLVER             : Gravity solver
//                                  POISSON_AND_GRAVITY_SOLVER : Poisson + Gravity solvers
//                lv          : Targeted refinement level 
//                TimeNew     : Targeted physical time to reach (only useful for adding external potential)
//                TimeOld     : Physical time before update     (only useful for adding external potential with UNSPLIT_GRAVITY)
//                NPG         : Number of patch groups to be upcated at a time
//                ArrayID     : Array index to load and store data ( 0 or 1 )
//                dt          : Time interval to advance solution (for the fluid and gravity solvers)
//                Poi_Coeff   : Coefficient in front of the RHS in the Poisson eq.
//-------------------------------------------------------------------------------------------------------
void Solver( const Solver_t TSolver, const int lv, const double TimeNew, const double TimeOld,
             const int NPG, const int ArrayID, const double dt, const double Poi_Coeff )
{

   const double dh = amr->dh[lv];

#  ifdef GRAVITY
   const bool POISSON_ON  = true;
   const bool GRAVITY_ON  = true;
   const bool POISSON_OFF = false;
   const bool GRAVITY_OFF = false;
#  else
   const OptGravityType_t OPT__GRAVITY_TYPE = GRAVITY_NONE;
#  endif // #ifdef GRAVITY ... else ...

// define useless variables in different models
#  if ( MODEL != ELBDM )
   const double ELBDM_ETA               = NULL_REAL;
   const double ELBDM_TAYLOR3_COEFF     = NULL_REAL;
   const bool ELBDM_TAYLOR3_AUTO        = NULL_BOOL;
#  endif

#  if ( MODEL != HYDRO )
   const bool   Flu_XYZ                 = true;
   const double GAMMA                   = NULL_REAL;
   const double MINMOD_COEFF            = NULL_REAL;
   const double EP_COEFF                = NULL_REAL;
   const LR_Limiter_t  OPT__LR_LIMITER  = LR_LIMITER_NONE;
   const WAF_Limiter_t OPT__WAF_LIMITER = WAF_LIMITER_NONE;
#  else
   const bool Flu_XYZ                   = 1 - ( AdvanceCounter[lv]%2 );    // forward/backward sweep
#  endif

#  ifndef QUARTIC_SELF_INTERACTION
   const double ELBDM_LAMBDA            = NULL_REAL;
#  endif

#  ifndef UNSPLIT_GRAVITY
   real (*h_Pot_Array_USG_F[2])[USG_NXT_F ][USG_NXT_F ][USG_NXT_F ] = { NULL, NULL };
#  ifdef GRAVITY
   real (*h_Pot_Array_USG_G[2])[USG_NXT_G ][USG_NXT_G ][USG_NXT_G ] = { NULL, NULL };
   real (*h_Flu_Array_USG_G[2])[GRA_NIN-1][PS1][PS1][PS1]           = { NULL, NULL };
#  endif
#  endif

#  if ( MODEL != HYDRO  &&  MODEL != ELBDM )
#  error : ERROR : ADD THE MODEL-DEPENDENT USELESS VARIABLES FOR THE NEW MODELS HERE
#  endif


   switch ( TSolver )
   {
      case FLUID_SOLVER :  

#        ifdef GPU
         CUAPI_Asyn_FluidSolver( h_Flu_Array_F_In[ArrayID], h_Flu_Array_F_Out[ArrayID], h_Flux_Array[ArrayID], 
                                 h_Corner_Array_F[ArrayID], h_MinDtInfo_Fluid_Array[ArrayID], h_Pot_Array_USG_F[ArrayID],
                                 NPG, dt, dh, GAMMA, OPT__FIXUP_FLUX, Flu_XYZ, 
                                 OPT__LR_LIMITER, MINMOD_COEFF, EP_COEFF, OPT__WAF_LIMITER, ELBDM_ETA, 
                                 ELBDM_TAYLOR3_COEFF, ELBDM_TAYLOR3_AUTO, OPT__ADAPTIVE_DT, TimeOld, OPT__GRAVITY_TYPE,
                                 GPU_NSTREAM );
#        else
         CPU_FluidSolver       ( h_Flu_Array_F_In[ArrayID], h_Flu_Array_F_Out[ArrayID], h_Flux_Array[ArrayID], 
                                 h_Corner_Array_F[ArrayID], h_MinDtInfo_Fluid_Array[ArrayID], h_Pot_Array_USG_F[ArrayID],
                                 NPG, dt, dh, GAMMA, OPT__FIXUP_FLUX, Flu_XYZ, 
                                 OPT__LR_LIMITER, MINMOD_COEFF, EP_COEFF, OPT__WAF_LIMITER, ELBDM_ETA, 
                                 ELBDM_TAYLOR3_COEFF, ELBDM_TAYLOR3_AUTO, OPT__ADAPTIVE_DT, TimeOld, OPT__GRAVITY_TYPE );
#        endif
         break;


#     ifdef GRAVITY
      case POISSON_SOLVER :  

#        ifdef GPU     
         CUAPI_Asyn_PoissonGravitySolver( h_Rho_Array_P[ArrayID], h_Pot_Array_P_In[ArrayID], 
                                          h_Pot_Array_P_Out[ArrayID], NULL, NULL,
                                          NULL, NULL,
                                          NPG, dt, dh, SOR_MIN_ITER, SOR_MAX_ITER, 
                                          SOR_OMEGA, MG_MAX_ITER, MG_NPRE_SMOOTH, MG_NPOST_SMOOTH, 
                                          MG_TOLERATED_ERROR, Poi_Coeff, OPT__POT_INT_SCHEME, 
                                          NULL_BOOL, ELBDM_ETA, NULL_REAL, POISSON_ON, GRAVITY_OFF, GPU_NSTREAM,
                                          GRAVITY_NONE, NULL_REAL, NULL_REAL, NULL_BOOL );
#        else
         CPU_PoissonGravitySolver       ( h_Rho_Array_P[ArrayID], h_Pot_Array_P_In[ArrayID], 
                                          h_Pot_Array_P_Out[ArrayID], NULL, NULL,
                                          NULL, NULL,
                                          NPG, dt, dh, SOR_MIN_ITER, SOR_MAX_ITER, 
                                          SOR_OMEGA, MG_MAX_ITER, MG_NPRE_SMOOTH, MG_NPOST_SMOOTH, 
                                          MG_TOLERATED_ERROR, Poi_Coeff, OPT__POT_INT_SCHEME, 
                                          NULL_BOOL, ELBDM_ETA, NULL_REAL, POISSON_ON, GRAVITY_OFF,
                                          GRAVITY_NONE, NULL_REAL, NULL_REAL, NULL_BOOL );
#        endif
         break;


      case GRAVITY_SOLVER :  
              
#        ifdef GPU     
         CUAPI_Asyn_PoissonGravitySolver( NULL, NULL, 
                                          h_Pot_Array_P_Out[ArrayID], h_Flu_Array_G[ArrayID], h_Corner_Array_G[ArrayID],
                                          h_Pot_Array_USG_G[ArrayID], h_Flu_Array_USG_G[ArrayID],
                                          NPG, dt, dh, NULL_INT, NULL_INT, 
                                          NULL_REAL, NULL_INT, NULL_INT, NULL_INT, 
                                          NULL_REAL, NULL_REAL, (IntScheme_t)NULL_INT, 
                                          OPT__GRA_P5_GRADIENT, ELBDM_ETA, ELBDM_LAMBDA, POISSON_OFF, GRAVITY_ON, GPU_NSTREAM,
                                          OPT__GRAVITY_TYPE, TimeNew, TimeOld, OPT__EXTERNAL_POT );
#        else
         CPU_PoissonGravitySolver       ( NULL, NULL, 
                                          h_Pot_Array_P_Out[ArrayID], h_Flu_Array_G[ArrayID], h_Corner_Array_G[ArrayID],
                                          h_Pot_Array_USG_G[ArrayID], h_Flu_Array_USG_G[ArrayID],
                                          NPG, dt, dh, NULL_INT, NULL_INT, 
                                          NULL_REAL, NULL_INT, NULL_INT, NULL_INT, 
                                          NULL_REAL, NULL_REAL, (IntScheme_t)NULL_INT, 
                                          OPT__GRA_P5_GRADIENT, ELBDM_ETA, ELBDM_LAMBDA, POISSON_OFF, GRAVITY_ON, 
                                          OPT__GRAVITY_TYPE, TimeNew, TimeOld, OPT__EXTERNAL_POT );
#        endif
         break;


      case POISSON_AND_GRAVITY_SOLVER :  

#        ifdef GPU     
         CUAPI_Asyn_PoissonGravitySolver( h_Rho_Array_P[ArrayID], h_Pot_Array_P_In[ArrayID], 
                                          h_Pot_Array_P_Out[ArrayID], h_Flu_Array_G[ArrayID], h_Corner_Array_G[ArrayID],
                                          h_Pot_Array_USG_G[ArrayID], h_Flu_Array_USG_G[ArrayID],
                                          NPG, dt, dh, SOR_MIN_ITER, SOR_MAX_ITER, 
                                          SOR_OMEGA, MG_MAX_ITER, MG_NPRE_SMOOTH, MG_NPOST_SMOOTH, 
                                          MG_TOLERATED_ERROR, Poi_Coeff, OPT__POT_INT_SCHEME, 
                                          OPT__GRA_P5_GRADIENT, ELBDM_ETA, ELBDM_LAMBDA, POISSON_ON, GRAVITY_ON, GPU_NSTREAM,
                                          OPT__GRAVITY_TYPE, TimeNew, TimeOld, OPT__EXTERNAL_POT );
#        else
         CPU_PoissonGravitySolver       ( h_Rho_Array_P[ArrayID], h_Pot_Array_P_In[ArrayID], 
                                          h_Pot_Array_P_Out[ArrayID], h_Flu_Array_G[ArrayID], h_Corner_Array_G[ArrayID],
                                          h_Pot_Array_USG_G[ArrayID], h_Flu_Array_USG_G[ArrayID],
                                          NPG, dt, dh, SOR_MIN_ITER, SOR_MAX_ITER, 
                                          SOR_OMEGA, MG_MAX_ITER, MG_NPRE_SMOOTH, MG_NPOST_SMOOTH, 
                                          MG_TOLERATED_ERROR, Poi_Coeff, OPT__POT_INT_SCHEME, 
                                          OPT__GRA_P5_GRADIENT, ELBDM_ETA, ELBDM_LAMBDA, POISSON_ON, GRAVITY_ON,
                                          OPT__GRAVITY_TYPE, TimeNew, TimeOld, OPT__EXTERNAL_POT );
#        endif
         break;

#     endif // #ifdef GRAVITY


      default :
         Aux_Error( ERROR_INFO, "incorrect parameter %s = %d !!\n", "TSolver", TSolver );

   } // switch ( TSolver )

} // FUNCTION : Solver



//-------------------------------------------------------------------------------------------------------
// Function    :  Closing_Step
// Description :  Store the updated solutions back to the patch pointers 
//
// Note        :  Use the input parameter "TSolver" to control the targeted solver 
//
// Parameter   :  TSolver     : Targeted solver
//                              --> FLUID_SOLVER               : Fluid / ELBDM solver
//                                  POISSON_SOLVER             : Poisson solver
//                                  GRAVITY_SOLVER             : Gravity solver
//                                  POISSON_AND_GRAVITY_SOLVER : Poisson + Gravity solvers
//                lv          : Targeted refinement level 
//                SaveSg_Flu  : Sandglass to store the updated fluid data (for both the fluid and gravity solvers)
//                SaveSg_Pot  : Sandglass to store the updated potential data (for the Poisson solver)
//                NPG         : Number of patch groups to be evaluated at a time
//                PID0_List   : List recording the patch indicies with LocalID==0 to be udpated
//                ArrayID     : Array index to load and store data ( 0 or 1 )
//                dt          : Time interval to advance solution (for OPT__CORR_UNPHY in Flu_Close)
//-------------------------------------------------------------------------------------------------------
void Closing_Step( const Solver_t TSolver, const int lv, const int SaveSg_Flu, const int SaveSg_Pot, const int NPG, 
                   const int *PID0_List, const int ArrayID, const double dt )
{

   switch ( TSolver )
   {
      case FLUID_SOLVER :   
         Flu_Close( lv, SaveSg_Flu, h_Flux_Array[ArrayID], h_Flu_Array_F_Out[ArrayID], 
                    h_MinDtInfo_Fluid_Array[ArrayID], NPG, PID0_List, OPT__ADAPTIVE_DT,
                    h_Flu_Array_F_In[ArrayID], dt );
         break;

#     ifdef GRAVITY
      case POISSON_SOLVER : 
         Poi_Close( lv, SaveSg_Pot, h_Pot_Array_P_Out[ArrayID], NPG, PID0_List ); 
         break;

      case GRAVITY_SOLVER :  
         Gra_Close( lv, SaveSg_Flu, h_Flu_Array_G    [ArrayID], NPG, PID0_List );
         break;

      case POISSON_AND_GRAVITY_SOLVER :  
         Poi_Close( lv, SaveSg_Pot, h_Pot_Array_P_Out[ArrayID], NPG, PID0_List ); 
         Gra_Close( lv, SaveSg_Flu, h_Flu_Array_G    [ArrayID], NPG, PID0_List );
         break;
#     endif

      default:
         Aux_Error( ERROR_INFO, "incorrect parameter %s = %d !!\n", "TSolver", TSolver );

   } // switch ( TSolver )

} // FUNCTION : Closing_Step


