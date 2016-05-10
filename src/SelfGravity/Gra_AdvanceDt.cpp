#include "Copyright.h"
#include "GAMER.h"

#ifdef GRAVITY



extern Timer_t *Timer_Gra_Advance[NLEVEL];
extern Timer_t *Timer_GetBuf     [NLEVEL][8];




//-------------------------------------------------------------------------------------------------------
// Function    :  Gra_AdvanceDt
// Description :  Solve the Poisson equation and advance the fluid variables by the gravitational acceleration
//
// Note        :  1. Poisson solver : lv = 0 : invoke the function "CPU_PoissonSolver_FFT"
//                                    lv > 0 : invoke the function "InvokeSolver"
//                2. Gravity solver : invoke the function "InvokeSolver"
//                3. The updated potential and fluid variables will be stored in the same sandglass 
//                4. PotSg at lv=0 will be updated here, but PotSg at at lv>0 and FluSg at lv>=0 will NOT be updated
//                   (they will be updated in EvolveLevel instead)
//                   --> It is because the lv-0 Poisson and Gravity solvers are invoked separately, and Gravity solver
//                       needs to call PreparePatchData to get the updated potential
//
// Parameter   :  lv             : Targeted refinement level 
//                TimeNew        : Targeted physical time to reach
//                TimeOld        : Physical time before update
//                                 --> For Gravity solver, this function updates physical time from TimeOld to TimeNew
//                                     For Poisson solver, this function calculates potential at **TimeNew**
//                dt             : Time interval to advance solution (can be different from TimeNew-TimeOld if COMOVING is on)
//                SaveSg_Flu     : Sandglass to store the updated fluid data (for the gravity solver)
//                SaveSg_Pot     : Sandglass to store the updated potential data  (for the Poisson solver)
//                Poisson        : true --> invoke the Poisson solver to evaluate the gravitational potential
//                Gravity        : true --> invoke the Gravity solver to evolve fluid by the gravitational acceleration
//                OverlapMPI     : true --> Overlap MPI time with CPU/GPU computation
//                Overlap_Sync   : true  --> Advance the patches which cannot be overlapped with MPI communication
//                                 false --> Advance the patches which can    be overlapped with MPI communication
//                                 (useful only if "OverlapMPI == true")
//-------------------------------------------------------------------------------------------------------
void Gra_AdvanceDt( const int lv, const double TimeNew, const double TimeOld, const double dt,
                    const int SaveSg_Flu, const int SaveSg_Pot, const bool Poisson, const bool Gravity,
                    const bool OverlapMPI, const bool Overlap_Sync )
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

   if ( !Poisson  &&  Gravity  &&  OPT__GRAVITY_TYPE != GRAVITY_EXTERNAL )
      Aux_Message( stderr, "WARNING : Poisson=off, Gravity=on, GravityType!=external --> ARE YOU SURE ?!\n" );


// coefficient in front of the RHS in the Poisson eq.
#  ifdef COMOVING
   const double Poi_Coeff = 4.0*M_PI*NEWTON_G*TimeNew;   // use TimeNew for calculating potential
#  else
   const double Poi_Coeff = 4.0*M_PI*NEWTON_G;
#  endif


// collect particles from all descendant patches
#  ifdef PARTICLE
   if ( Poisson )    Par_CollectParticleForPoisson( lv );
#  endif


// the base-level Poisson solver is implemented using the FFTW library (with CPUs only)
   if ( lv == 0 )    
   {
//    do not need to calculate the gravitational potential if self-gravity is disabled
      if ( Poisson )
      {
         TIMING_FUNC(   CPU_PoissonSolver_FFT( Poi_Coeff, SaveSg_Pot, TimeNew ),
                        Timer_Gra_Advance[lv],   false   );   

         amr->PotSg    [0]             = SaveSg_Pot;
         amr->PotSgTime[0][SaveSg_Pot] = TimeNew;

         TIMING_FUNC(   Buf_GetBufferData( lv, NULL_INT, SaveSg_Pot, POT_FOR_POISSON, _POTE, Pot_ParaBuf, USELB_YES ),
                        Timer_GetBuf[lv][1],   true   );
      }

      if ( Gravity )  
      {
//       TIMING_FUNC(   InvokeSolver( GRAVITY_SOLVER, lv, TimeNew, TimeOld, dt, NULL_REAL, SaveSg_Flu, NULL_INT, 
//                                    OverlapMPI, Overlap_Sync ),
//                      Timer_Gra_Advance[lv],   true   );
  
         TIMING_FUNC(   InvokeSolver( GRAVITY_SOLVER, lv, TimeNew, TimeOld, dt, NULL_REAL, SaveSg_Flu, NULL_INT, false, false ),
                        Timer_Gra_Advance[lv],   true   );

         if ( OPT__RESET_FLUID )    Flu_ResetByUser( lv, SaveSg_Flu, TimeNew );

         TIMING_FUNC(   Buf_GetBufferData( lv, SaveSg_Flu, NULL_INT, DATA_GENERAL, _FLU, Flu_ParaBuf, USELB_YES ),
                        Timer_GetBuf[lv][2],   true   );

         amr    ->FluSg[0] = SaveSg_Flu;
      } // if ( Gravity )
   } // if ( lv == 0 )


   else // lv > 0 
   {
      if      (  Poisson  &&  !Gravity )
      InvokeSolver( POISSON_SOLVER,             lv, TimeNew, TimeOld, NULL_REAL, Poi_Coeff, NULL_INT,   SaveSg_Pot,
                    OverlapMPI, Overlap_Sync );

      else if ( !Poisson  &&   Gravity )  
      InvokeSolver( GRAVITY_SOLVER,             lv, TimeNew, TimeOld, dt,        NULL_REAL, SaveSg_Flu, NULL_INT,
                    OverlapMPI, Overlap_Sync );

      else if (  Poisson  &&   Gravity )  
      InvokeSolver( POISSON_AND_GRAVITY_SOLVER, lv, TimeNew, TimeOld, dt,        Poi_Coeff, SaveSg_Flu, SaveSg_Pot,
                    OverlapMPI, Overlap_Sync );
      
      if ( Gravity  &&  OPT__RESET_FLUID )   Flu_ResetByUser( lv, SaveSg_Flu, TimeNew );
   }


// free variables of descendant particles
#  ifdef PARTICLE
   if ( Poisson )    Par_ReleaseParticleForPoisson( lv );
#  endif

} // FUNCTION : Gra_AdvanceDt



#ifdef PARTICLE
//-------------------------------------------------------------------------------------------------------
// Function    :  Par_CollectParticleForPoisson
// Description :  Collect particles from all descendants (sons, grandsons, ...) for the Poisson solver
//
// Note        :  1. Invoke Par_CollectParticleFromDescendant for all Lv=lv patches and Lv=lv-1 patches
//                   adjacent to the lv <-> lv-1 boundaries
//                2. One must call Par_ReleaseParticleForPoisson later to properly free the array "ParList_Desc"
//                3. This function will also be invoded by Main when DEBUG is on
//
// Parameter   :  lv : Targeted refinement level 
//-------------------------------------------------------------------------------------------------------
void Par_CollectParticleForPoisson( const int lv )
{

// collect particles for all lv patches
   for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)   Par_CollectParticleFromDescendant( lv, PID );

// collect particles for all lv-1 patches adjacent to the lv <-> lv-1 boundaries
   if ( lv > 0 )
   {
      const int FaLv = lv - 1;
      int SibPID;

      for (int PID=0; PID<amr->NPatchComma[FaLv][1]; PID++)
      {
         if ( amr->patch[0][FaLv][PID]->son == -1 )
         {
            for (int s=0; s<26; s++)
            {
               SibPID = amr->patch[0][FaLv][PID]->sibling[s];

               if ( SibPID >= 0  &&  amr->patch[0][FaLv][SibPID]->son != -1 )
               {
                  Par_CollectParticleFromDescendant( FaLv, PID );
                  break;
               }
            }
         }
      } // for (int PID=0; PID<amr->NPatchComma[FaLv][1]; PID++)
   } // if ( lv > 0 )

} // FUNCTION : Par_CollectParticleForPoisson



//-------------------------------------------------------------------------------------------------------
// Function    :  Par_ReleaseParticleForPoisson
// Description :  Release the memory allocated by Par_CollectParticleForPoisson
//
// Note        :  1. This function will also be invoded by Main when DEBUG is on
//
// Parameter   :  lv : Targeted refinement level 
//-------------------------------------------------------------------------------------------------------
void Par_ReleaseParticleForPoisson( const int lv )
{

   for (int TLv=lv; TLv>=0 && TLv>=lv-1; TLv--)
   for (int PID=0; PID<amr->NPatchComma[TLv][1]; PID++)
   {
      if ( amr->patch[0][TLv][PID]->ParList_Desc != NULL )
      {
         delete [] amr->patch[0][TLv][PID]->ParList_Desc;
         amr->patch[0][TLv][PID]->ParList_Desc = NULL;
      }

      amr->patch[0][TLv][PID]->NPar_Desc = -1;
   }

} // FUNCTION : Par_ReleaseParticleForPoisson
#endif // #ifdef PARTICLE



#endif // #ifdef GRAVITY
