#include "GAMER.h"

#ifdef PARTICLE




//-------------------------------------------------------------------------------------------------------
// Function    :  Par_GetTimeStep_VelAcc
// Description :  Estimate the evolution time-step from the maximum particle velocity and acceleration
//                --> dt_vel = DT__PARVEL*dh/v_max, where v_max = max(vx,vy,vz,all_particles)
//                    dt_acc = DT__PARACC*(dh/a_max)^0.5, where a_max = max(accx,accy,accz,all_particles)
//
// Note        :  1. This function should be applied to both physical and comoving coordinates and always
//                   return the evolution time-step (dt) actually used in various solvers
//                   --> Physical coordinates : dt = physical time interval
//                       Comoving coordinates : dt = delta(scale_factor) / ( Hubble_parameter*scale_factor^3 )
//                   --> We convert dt back to the physical time interval, which equals "delta(scale_factor)"
//                       in the comoving coordinates, in Mis_GetTimeStep()
//                2. Particle acceleration criterion is used only when DT__PARACC > 0.0
//                3. Particles in non-leaf patches are also included by default
//                   --> Controlled by "IncNonleaf"
//
// Parameter   :  dt_vel : Evolution time-step estimated from the particle velocity
//                         --> Call-by-reference
//                dt_acc : Evolution time-step estimated from the particle acceleration
//                         --> Call-by-reference
//                lv     : Target refinement level
//
// Return      :  dt_vel, dt_acc
//-------------------------------------------------------------------------------------------------------
void Par_GetTimeStep_VelAcc( double &dt_vel, double &dt_acc, const int lv )
{

   const bool  IncNonleaf     = true;
   const real_par *Vel[3]     = { amr->Par->VelX, amr->Par->VelY, amr->Par->VelZ };
#  ifdef STORE_PAR_ACC
   const real_par *Acc[3]     = { amr->Par->AccX, amr->Par->AccY, amr->Par->AccZ };
#  else
   const real_par *Acc[3]     = { NULL, NULL, NULL };
#  endif
#  ifdef MASSIVE_PARTICLES
   const bool  UseAcc         = ( DT__PARACC > 0.0 );
#  else
   const bool  UseAcc         = false;
#  endif
#  ifdef OPENMP
   const int   NT             = OMP_NTHREAD;   // number of OpenMP threads
#  else
   const int   NT             = 1;
#  endif
   const long_par *ParType    = amr->Par->Type;

#  ifndef STORE_PAR_ACC
   if ( UseAcc )
      Aux_Error( ERROR_INFO, "DT__PARACC (%14.7e) > 0.0 when STORE_PAR_ACC is off !!\n", DT__PARACC );
#  endif


// collect particles for non-leaf patches
   const bool PredictPos_No    = false;
   const bool SibBufPatch_No   = false;
   const bool FaSibBufPatch_No = false;
   const bool JustCountNPar_No = false;
   const bool TimingSendPar_No = false;
#  ifdef STORE_PAR_ACC
   const int  ParAccBIdx       = _PAR_ACC;
#  else
   const int  ParAccBIdx       = 0;
#  endif

   if ( IncNonleaf )
      Par_CollectParticle2OneLevel( lv, _PAR_VEL|((UseAcc)?ParAccBIdx:0), _PAR_TYPE, PredictPos_No,
                                    NULL_REAL, SibBufPatch_No, FaSibBufPatch_No, JustCountNPar_No,
                                    TimingSendPar_No );


// get the maximum particle velocity and acceleration on the target level
   real_par MaxVel, MaxAcc;
   long NParVel=0, NParAcc=0;

   real_par *MaxVel_OMP = new real_par [NT];
   real_par *MaxAcc_OMP = new real_par [NT];

   MaxVel = (real_par)0.0;  // don't assign negative values since we assume it to be positive-definite
   MaxAcc = (real_par)0.0;

   for (int t=0; t<NT; t++)
   {
      MaxVel_OMP[t] = (real_par)0.0;
      MaxAcc_OMP[t] = (real_par)0.0;
   }

//###NOTE: OpenMP may not improve performance here
#  pragma omp parallel
   {
#     ifdef OPENMP
      const int TID = omp_get_thread_num();
#     else
      const int TID = 0;
#     endif

#     pragma omp for reduction( +:NParVel, NParAcc ) schedule( PAR_OMP_SCHED, PAR_OMP_SCHED_CHUNK )
      for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
      {
         int   NParThisPatch;
         long *ParList = NULL;
         bool  UseCopy;

//       leaf patches
         if ( amr->patch[0][lv][PID]->son == -1 )
         {
            NParThisPatch = amr->patch[0][lv][PID]->NPar;
            ParList       = amr->patch[0][lv][PID]->ParList;
            UseCopy       = false;
         }

//       include non-leaf patches
         else if ( IncNonleaf )
         {
            NParThisPatch = amr->patch[0][lv][PID]->NPar_Copy;
#           ifdef LOAD_BALANCE
            ParList       = NULL;
            UseCopy       = true;
#           else
            ParList       = amr->patch[0][lv][PID]->ParList_Copy;
            UseCopy       = false;
#           endif
         } // if ( amr->patch[0][lv][PID]->son == -1 ) ... elif ...

//       exclude non-leaf patches
         else
         {
            NParThisPatch = 0;
            ParList       = NULL;
            UseCopy       = false;
         } // if ( amr->patch[0][lv][PID]->son == -1 ) ... elif ... else ...


         if ( UseCopy )
         {
//          ParAttFlt/Int_Copy[] is only defined in LOAD_BALANCE
            real_par *Vel_Copy[3] = { NULL, NULL, NULL };
            real_par *Acc_Copy[3] = { NULL, NULL, NULL };
            long_par *Typ_Copy    = NULL;
#           ifdef LOAD_BALANCE
            for (int d=0; d<3; d++)
            {
               Vel_Copy[d] = amr->patch[0][lv][PID]->ParAttFlt_Copy[ PAR_VELX + d ];
#              ifdef STORE_PAR_ACC
               Acc_Copy[d] = amr->patch[0][lv][PID]->ParAttFlt_Copy[ PAR_ACCX + d ];
#              endif
            }
            Typ_Copy = amr->patch[0][lv][PID]->ParAttInt_Copy[ PAR_TYPE ];
#           endif

            for (int p=0; p<NParThisPatch; p++)
            {
//             don't check velocity for tracer particles unless enabling OPT__FREEZE_FLUID
//             since the fluid CFL condition should be sufficient
//             --> this ensures that tracer particles do not change the simulation results
               if ( Typ_Copy[p] != PTYPE_TRACER  ||  OPT__FREEZE_FLUID )
               {
                  for (int d=0; d<3; d++)
                  MaxVel_OMP[TID] = MAX( MaxVel_OMP[TID], FABS(Vel_Copy[d][p]) );

                  NParVel ++;
               }

//             don't check acceleration for tracer particles
               if ( UseAcc  &&  Typ_Copy[p] != PTYPE_TRACER )
               {
                  for (int d=0; d<3; d++)
                  MaxAcc_OMP[TID] = MAX( MaxAcc_OMP[TID], FABS(Acc_Copy[d][p]) );

                  NParAcc ++;
               }
            }
         } // if ( UseCopy )

         else
         {
            for (int p=0; p<NParThisPatch; p++)
            {
               const long ParID = ParList[p];

//             don't check velocity for tracer particles unless enabling OPT__FREEZE_FLUID
//             since the fluid CFL condition should be sufficient
//             --> this ensures that tracer particles do not change the simulation results
               if ( ParType[ParID] != PTYPE_TRACER  ||  OPT__FREEZE_FLUID )
               {
                  for (int d=0; d<3; d++)
                  MaxVel_OMP[TID] = MAX( MaxVel_OMP[TID], FABS(Vel[d][ParID]) );

                  NParVel ++;
               }

//             don't check acceleration for tracer particles
               if ( UseAcc  &&  ParType[ParID] != PTYPE_TRACER )
               {
                  for (int d=0; d<3; d++)
                  MaxAcc_OMP[TID] = MAX( MaxAcc_OMP[TID], FABS(Acc[d][ParID]) );

                  NParAcc ++;
               }
            }
         } // if ( UseCopy ) ... else ...
      } // for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)
   } // end of OpenMP parallel region


// compare the maximum velocity and acceleration evaluated by different OMP threads
   for (int t=0; t<NT; t++)
   {
      MaxVel = MAX( MaxVel_OMP[t], MaxVel );

      if ( UseAcc )
      MaxAcc = MAX( MaxAcc_OMP[t], MaxAcc );
   }

   delete [] MaxVel_OMP;
   delete [] MaxAcc_OMP;


// get the time-step in this rank
   double dt_vel_local, dt_acc_local;
   dt_vel_local =       amr->dh[lv] / (double)MaxVel;

   if ( UseAcc )
   dt_acc_local = sqrt( amr->dh[lv] / (double)MaxAcc );


// get the minimum time-step in all ranks
   MPI_Allreduce( &dt_vel_local, &dt_vel, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD );
#  ifndef SERIAL
   MPI_Reduce( (MPI_Rank==0)?MPI_IN_PLACE:&NParVel, &NParVel, 1, MPI_LONG, MPI_SUM, 0, MPI_COMM_WORLD );
#  endif

   if ( UseAcc ) {
   MPI_Allreduce( &dt_acc_local, &dt_acc, 1, MPI_DOUBLE, MPI_MIN, MPI_COMM_WORLD );
#  ifndef SERIAL
   MPI_Reduce( (MPI_Rank==0)?MPI_IN_PLACE:&NParAcc, &NParAcc, 1, MPI_LONG, MPI_SUM, 0, MPI_COMM_WORLD );
#  endif
   }


// verify the minimum time-step
   if ( !Aux_IsFinite(dt_vel)  &&  NParVel > 0  &&  MPI_Rank == 0 )
   {
      Aux_Message( stderr, "WARNING : time-step estimation by particle velocity is incorrect (dt_vel = %13.7e, NParVel = %ld) !!\n", dt_vel, NParVel );
      Aux_Message( stderr, "          --> Likely all particles have zero velocity\n" );

      if ( DT__PARVEL_MAX < 0.0 )
      Aux_Message( stderr, "          --> You might want to set DT__PARVEL_MAX properly\n" );
   }

   if ( UseAcc  &&  !Aux_IsFinite(dt_acc)  &&  NParAcc > 0  &&  MPI_Rank == 0 )
      Aux_Error( ERROR_INFO, "time-step estimation by particle acceleration is incorrect (dt_acc = %13.7e, NParAcc = %ld) !!\n", dt_acc, NParAcc );


// multiply by the safety factor
   dt_vel *= DT__PARVEL;
   if ( DT__PARVEL_MAX >= 0.0 )  dt_vel = MIN( dt_vel, DT__PARVEL_MAX );

   if ( UseAcc )
   dt_acc *= DT__PARACC;


// free memory allocated by Par_CollectParticle2OneLevel()
   if ( IncNonleaf )
      Par_CollectParticle2OneLevel_FreeMemory( lv, SibBufPatch_No, FaSibBufPatch_No );

} // FUNCTION : Par_GetTimeStep_VelAcc



#endif // #ifdef PARTICLE
