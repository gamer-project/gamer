#include "GAMER.h"

#ifdef GRAVITY

#ifdef TIMING
extern Timer_t *Timer_Gra_Advance[NLEVEL];
extern Timer_t *Timer_GetBuf     [NLEVEL][8];
extern Timer_t *Timer_Par_Collect[NLEVEL];
#endif

//   Aux_Message( stdout, "DEBUG: Rank: %d, %s, lv: %d, TF: %d, %ld <-> %ld\n", MPI_Rank, __FUNCTION__, lv, FullRefinedLv, (long)NPatchTotal[lv], NX0_TOT[0]*NX0_TOT[1]*NX0_TOT[2]/512*(1L<<lv)*(1L<<lv)*(1L<<lv) );
   // Aux_Message( stdout, "\n" );
   // Aux_Message( stdout, "Track flu prepare %d %24.16e %d %24.16e %s %d\n", amr->FluSg[lv], amr->FluSgTime[lv][amr->FluSg[lv]], 1-amr->FluSg[lv], amr->FluSgTime[lv][1-amr->FluSg[lv]], __FUNCTION__, __LINE__ );
   // Aux_Message( stdout, "Track pot prepare %d %24.16e %d %24.16e %s %d\n", amr->PotSg[lv], amr->PotSgTime[lv][amr->PotSg[lv]], 1-amr->PotSg[lv], amr->PotSgTime[lv][1-amr->PotSg[lv]], __FUNCTION__, __LINE__ );
//      Aux_Message( stdout, "DEBUG: %s %d\n", __FILE__, __LINE__ );



//-------------------------------------------------------------------------------------------------------
// Function    :  Gra_AdvanceDt
// Description :  Solve the Poisson equation and advance the fluid variables by the gravitational acceleration
//
// Note        :  1. Poisson solver : lv = 0 : invoke CPU_PoissonSolver_FFT()
//                                    lv > 0 : invoke InvokeSolver()
//                2. Gravity solver : invoke InvokeSolver()
//                3. The updated potential and fluid variables will be stored in the same sandglass
//                4. PotSg at lv=0 will be updated here, but PotSg at at lv>0 and FluSg at lv>=0 will NOT be updated
//                   (they will be updated in EvolveLevel instead)
//                   --> It is because the lv-0 Poisson and Gravity solvers are invoked separately, and Gravity solver
//                       needs to call Prepare_PatchData to get the updated potential
//
// Parameter   :  lv           : Target refinement level
//                TimeNew      : Target physical time to reach
//                TimeOld      : Physical time before update
//                               --> For Gravity solver, this function updates physical time from TimeOld to TimeNew
//                                   For Poisson solver, this function calculates potential at **TimeNew**
//                dt           : Time interval to advance solution (can be different from TimeNew-TimeOld if COMOVING is on)
//                SaveSg_Flu   : Sandglass to store the updated fluid data (for the gravity solver)
//                SaveSg_Pot   : Sandglass to store the updated potential data  (for the Poisson solver)
//                Poisson      : true --> invoke the Poisson solver to evaluate the gravitational potential
//                                        (including both self-gravity potential and external potential)
//                Gravity      : true --> invoke the Gravity solver to evolve fluid by the gravitational acceleration
//                                        (including self-gravity, external potentional, and external acceleration)
//                OverlapMPI   : true --> Overlap MPI time with CPU/GPU computation
//                Overlap_Sync : true  --> Advance the patches which cannot be overlapped with MPI communication
//                               false --> Advance the patches which can    be overlapped with MPI communication
//                               (useful only if "OverlapMPI == true")
//                Timing       : enable timing --> disable it in Flu_CorrAfterAllSync()
//-------------------------------------------------------------------------------------------------------
void Gra_AdvanceDt( const int lv, const double TimeNew, const double TimeOld, const double dt,
                    const int SaveSg_Flu, const int SaveSg_Pot, const bool Poisson, const bool Gravity,
                    const bool OverlapMPI, const bool Overlap_Sync, const bool Timing )
{

// check
   if ( Poisson  &&  ( SaveSg_Pot != 0 &&  SaveSg_Pot != 1 )  )
      Aux_Error( ERROR_INFO, "incorrect SaveSg_Pot (%d) !!\n", SaveSg_Pot );

   if ( Gravity  &&  ( SaveSg_Flu != 0 &&  SaveSg_Flu != 1 )  )
      Aux_Error( ERROR_INFO, "incorrect SaveSg_Flu (%d) !!\n", SaveSg_Flu );

   if ( !Poisson  &&  !Gravity )
   {
      Aux_Message( stderr, "WARNING : nothing to do in %s !!\n", __FUNCTION__ );
      return;
   }

   if (  !Poisson  &&  Gravity  &&  ( OPT__SELF_GRAVITY || OPT__EXT_POT )  )
      Aux_Message( stderr, "WARNING : Poisson=off but Gravity=on --> ARE YOU SURE ?!\n" );


// whether we actually need potential
   const bool UsePot = (  Poisson  &&  ( OPT__SELF_GRAVITY || OPT__EXT_POT )  );


// coefficient in front of the RHS in the Poisson eq.
#  ifdef COMOVING
   const double Poi_Coeff = 4.0*M_PI*NEWTON_G*TimeNew;   // use TimeNew for calculating potential
#  else
   const double Poi_Coeff = 4.0*M_PI*NEWTON_G;
#  endif


// initialize the particle density array (rho_ext) and collect particles to the target level
#  ifdef PARTICLE
   const bool TimingSendPar_Yes = Timing;
   const bool JustCountNPar_No  = false;
#  ifdef LOAD_BALANCE
   const bool PredictPos        = amr->Par->PredictPos;
   const bool SibBufPatch       = true;
   const bool FaSibBufPatch     = true;
#  else
   const bool PredictPos        = false;
   const bool SibBufPatch       = NULL_BOOL;
   const bool FaSibBufPatch     = NULL_BOOL;
#  endif

   if ( UsePot )
   {
      TIMING_FUNC(   Par_CollectParticle2OneLevel( lv, _PAR_MASS|_PAR_POSX|_PAR_POSY|_PAR_POSZ, _PAR_TYPE, PredictPos,
                                                   TimeNew, SibBufPatch, FaSibBufPatch, JustCountNPar_No, TimingSendPar_Yes ),
                     Timer_Par_Collect[lv],   Timing   );

      TIMING_FUNC(   Prepare_PatchData_InitParticleDensityArray( lv, TimeNew ),
                     Timer_Par_Collect[lv],   Timing   );
   }
#  endif // #ifdef PARTICLE


// user-specified work before invoking the Poisson/gravity solvers
// --> call it even when UsePot==false in order to support external acceleration (OPT__EXT_ACC)
   if ( Poi_UserWorkBeforePoisson_Ptr != NULL )
      TIMING_FUNC(   Poi_UserWorkBeforePoisson_Ptr( TimeNew, lv ),
                     Timer_Gra_Advance[lv],   ( Timing && lv == 0 )   );


   bool FullRefinedLv = false;
   if ( (long)NPatchTotal[lv] == NX0_TOT[0]*NX0_TOT[1]*NX0_TOT[2]/512*(1L<<lv)*(1L<<lv)*(1L<<lv) )   FullRefinedLv = true;


// the base-level Poisson solver is implemented using the FFTW library (with CPUs only)
   if ( FullRefinedLv )
   {
      if ( ! FFTW_Inited[lv] )   Init_FFTW( lv );

//    initialize the k-space Green's function for the isolated BC
      if ( ! GreenFuncK_Inited[lv] )
      if ( OPT__SELF_GRAVITY  &&  OPT__BC_POT == BC_POT_ISOLATED )    Init_GreenFuncK( lv );

//    do not need to calculate the gravitational potential if self-gravity is disabled
      if ( UsePot )
      {
#        ifdef SUPPORT_FFTW
         if ( OPT__SELF_GRAVITY )
         TIMING_FUNC(   CPU_PoissonSolver_FFT( Poi_Coeff, SaveSg_Pot, TimeNew, lv ),
                        Timer_Gra_Advance[lv],   Timing   );
#        endif

         if ( OPT__EXT_POT )
         TIMING_FUNC(   CPU_ExtPotSolver_FullyRefinedLevel( CPUExtPot_Ptr, ExtPot_AuxArray_Flt, ExtPot_AuxArray_Int,
                                                            h_ExtPotTable, h_ExtPotGenePtr,
                                                            TimeNew, OPT__SELF_GRAVITY, SaveSg_Pot, lv ),
                        Timer_Gra_Advance[lv],   Timing   );

         amr->PotSg    [lv]             = SaveSg_Pot;
         amr->PotSgTime[lv][SaveSg_Pot] = TimeNew;

//       note that the MPI bandwidth achieved in the following command may be much lower than normal
//       because of switching back from the MPI buffer used by FFTW
         TIMING_FUNC(   Buf_GetBufferData( lv, NULL_INT, NULL_INT, SaveSg_Pot, POT_FOR_POISSON, _POTE, _NONE, Pot_ParaBuf, USELB_YES ),
                        Timer_GetBuf[lv][1],   Timing  );

//       must call Poi_StorePotWithGhostZone AFTER collecting potential for buffer patches
#        ifdef STORE_POT_GHOST
         TIMING_FUNC(   Poi_StorePotWithGhostZone( lv, SaveSg_Pot, true ),
                        Timer_Gra_Advance[lv],   Timing   );
#        endif
      }

      if ( Gravity )
      {
//       TIMING_FUNC(   InvokeSolver( GRAVITY_SOLVER, lv, TimeNew, TimeOld, dt, NULL_REAL, SaveSg_Flu, NULL_INT, NULL_INT,
//                                    OverlapMPI, Overlap_Sync ),
//                      Timer_Gra_Advance[lv],   Timing  );

         TIMING_FUNC(   InvokeSolver( GRAVITY_SOLVER, lv, TimeNew, TimeOld, dt, NULL_REAL, SaveSg_Flu, NULL_INT, NULL_INT,
                                      false, false ),
                        Timer_Gra_Advance[lv],   Timing   );

         amr->FluSg[lv] = SaveSg_Flu;
      } // if ( Gravity )
   } // if ( FullRefinedLv )


   else // if ( FullRefinedLv )
   {
      if      (  Poisson  &&  !Gravity )
      {
#        if ( POT_SCHEME == HYPRE_POI )
         // Buf_GetBufferData( lv-1, NULL_INT, NULL_INT,   amr->PotSg[lv-1], POT_FOR_POISSON, _POTE, _NONE, Pot_ParaBuf, USELB_YES );
         // Buf_GetBufferData( lv-1, NULL_INT, NULL_INT, 1-amr->PotSg[lv-1], POT_FOR_POISSON, _POTE, _NONE, Pot_ParaBuf, USELB_YES );
         Hypre_SolvePoisson( SaveSg_Pot, lv, TimeNew, Poi_Coeff );
         TIMING_FUNC(   Buf_GetBufferData( lv, NULL_INT, NULL_INT, SaveSg_Pot, POT_FOR_POISSON, _POTE, _NONE, Pot_ParaBuf, USELB_YES ),
                        Timer_GetBuf[lv][1],   Timing  );
//       must call Poi_StorePotWithGhostZone AFTER collecting potential for buffer patches
#        ifdef STORE_POT_GHOST
         TIMING_FUNC(   Poi_StorePotWithGhostZone( lv, SaveSg_Pot, true ),
                        Timer_Gra_Advance[lv],   Timing   );
#        endif
#        else
         InvokeSolver( POISSON_SOLVER,             lv, TimeNew, TimeOld, NULL_REAL, Poi_Coeff, NULL_INT,   NULL_INT, SaveSg_Pot,
                       OverlapMPI, Overlap_Sync );
#        endif
      }

      else if ( !Poisson  &&   Gravity )
      {
         InvokeSolver( GRAVITY_SOLVER,             lv, TimeNew, TimeOld, dt,        NULL_REAL, SaveSg_Flu, NULL_INT, NULL_INT,
                       OverlapMPI, Overlap_Sync );
      }
      else if (  Poisson  &&   Gravity )
      {
#        if ( POT_SCHEME == HYPRE_POI )
         Hypre_SolvePoisson( SaveSg_Pot, lv, TimeNew, Poi_Coeff );
         amr->PotSg    [lv]             = SaveSg_Pot;
         amr->PotSgTime[lv][SaveSg_Pot] = TimeNew;
         TIMING_FUNC(   Buf_GetBufferData( lv, NULL_INT, NULL_INT, SaveSg_Pot, POT_FOR_POISSON, _POTE, _NONE, Pot_ParaBuf, USELB_YES ),
                        Timer_GetBuf[lv][1],   Timing  );
//       must call Poi_StorePotWithGhostZone AFTER collecting potential for buffer patches
#        ifdef STORE_POT_GHOST
         TIMING_FUNC(   Poi_StorePotWithGhostZone( lv, SaveSg_Pot, true ),
                        Timer_Gra_Advance[lv],   Timing   );
#        endif
         // Aux_Message( stdout, "%s update pot SG\n", __FILE__ );
         // Buf_GetBufferData( lv, NULL_INT, NULL_INT,   amr->PotSg[lv], POT_FOR_POISSON, _POTE, _NONE, Pot_ParaBuf, USELB_YES );
         // Buf_GetBufferData( lv, amr->FluSg[lv], NULL_INT, amr->PotSg[lv], DATA_GENERAL, _DENS|_POTE, _NONE, Rho_ParaBuf, USELB_YES );
         // Buf_GetBufferData( lv, NULL_INT, NULL_INT, amr->PotSg[lv], DATA_GENERAL, _POTE, _NONE, Rho_ParaBuf, USELB_YES );
         // Buf_GetBufferData( lv, amr->FluSg[lv], NULL_INT, amr->PotSg[lv], DATA_GENERAL, _DENS|_MOMX|_MOMY|_MOMZ|_ENGY|_POTE, _NONE, Rho_ParaBuf, USELB_YES );


         InvokeSolver( GRAVITY_SOLVER,             lv, TimeNew, TimeOld, dt,        NULL_REAL, SaveSg_Flu, NULL_INT, NULL_INT,
                       OverlapMPI, Overlap_Sync );
#        else
         InvokeSolver( POISSON_AND_GRAVITY_SOLVER, lv, TimeNew, TimeOld, dt,        Poi_Coeff, SaveSg_Flu, NULL_INT, SaveSg_Pot,
                       OverlapMPI, Overlap_Sync );
#        endif
      }
   }


// free memory for collecting particles from other ranks and levels, and free density arrays with ghost zones (rho_ext)
#  ifdef PARTICLE
   if ( UsePot )
   {
//    don't use the TIMING_FUNC macro since we don't want to call MPI_Barrier here even when OPT__TIMING_BARRIER is on
//    --> otherwise OPT__TIMING_BALANCE will fail because all ranks are synchronized before and after Gra_AdvanceDt
#     ifdef TIMING
      if ( Timing )  Timer_Par_Collect[lv]->Start();
#     endif

      Par_CollectParticle2OneLevel_FreeMemory( lv, SibBufPatch, FaSibBufPatch );
      Prepare_PatchData_FreeParticleDensityArray( lv );

#     ifdef TIMING
      if ( Timing )  Timer_Par_Collect[lv]->Stop();
#     endif
   }
#  endif

} // FUNCTION : Gra_AdvanceDt



#endif // #ifdef GRAVITY
