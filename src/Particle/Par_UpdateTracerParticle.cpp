#include "GAMER.h"

#if ( defined PARTICLE  &&  defined TRACER )




//-------------------------------------------------------------------------------------------------------
// Function    :  Par_UpdateTracerParticle
// Description :  Update tracer particle position and velocity at the target level
//
// Note        :  1. Does not take into account the "periodic B.C." when updating particles
//                   --> After update, the particle position may lie outside the simulation box
//                   --> It will be corrected by Par_PassParticle2Sibling()
//                2. Particle time may not be synchronized before invoking this function
//                   --> For example, particles just cross from coarse (lv) to fine (lv+1) grids may have time greater than
//                       other particles at lv+1. Also, particles just cross from fine (lv) to coarse (lv-1) grids may have
//                       time less than other particles at lv-1.
//                   --> This routine updates *all* tracer particles (even for those with time > TimeNew)
//                       --> Different from Par_UpdateParticle()
//                       --> For particles with time > TimeNew, they will be evolved backward in time
//                       --> So after invoking this function, all tracer particles on this level will
//                           be synchronized
//                3. The MapOnly mode maps the velocity to the particles, but does not advance
//                   their positions
//
// Parameter   :  lv      : Target refinement level
//                TimeNew : Target physical time to reach
//                TimeOld : Physical time before update
//                MapOnly : If true, only map the velocity to the particles
//                          and do not advance them
//-------------------------------------------------------------------------------------------------------
void Par_UpdateTracerParticle( const int lv, const double TimeNew, const double TimeOld,
                               const bool MapOnly )
{

   const bool     IntPhase_No       = false;
   const bool     DE_Consistency_No = false;
   const real     MinDens_No        = -1.0;
   const real     MinPres_No        = -1.0;
   const real     MinTemp_No        = -1.0;
   const real     MinEntr_No        = -1.0;
   const double   dh                = amr->dh[lv];
   const double   _dh               = 1.0/dh;
   const int      ParGhost          = amr->Par->GhostSizeTracer;
   const int      VelSize           = PS1 + 2*ParGhost;
   const bool     UseTracers_Yes    = true;
#  ifdef COMOVING
   const real_par dt_com            = (real_par)Mis_dTime2dt( TimeOld, TimeNew-TimeOld );
#  endif

         real_par *ParPos[3] = { amr->Par->PosX, amr->Par->PosY, amr->Par->PosZ };
         real_par *ParVel[3] = { amr->Par->VelX, amr->Par->VelY, amr->Par->VelZ };
         real_par *ParTime   = amr->Par->Time;
   const long_par *ParType   = amr->Par->Type;


// get the maximum number of particles in a single patch
// --> must use "NPar" instead of "NParType[(int)PTYPE_TRACER]" since currently
//     both Vel_Temp[] and InterpParPos[] still allocate memory for non-tracer particles
   int NParMax = 0;
   for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)   NParMax = MAX( NParMax, amr->patch[0][lv][PID]->NPar );


// nothing to do if there is no particle
   if ( NParMax <= 0 )  return;


// OpenMP parallel region
#  pragma omp parallel
   {

// per-thread variables
   real *VelX = new real [ 8*CUBE(VelSize) ];   // 8: number of patches per patch group
   real *VelY = new real [ 8*CUBE(VelSize) ];
   real *VelZ = new real [ 8*CUBE(VelSize) ];

   real_par **Vel_Temp     = NULL;
   real_par **InterpParPos = NULL;
   Aux_AllocateArray2D( Vel_Temp,     3, NParMax );
   Aux_AllocateArray2D( InterpParPos, 3, NParMax );

   bool     GotYou;
   long     ParID;
   real_par dt;

// loop over all **real** patch groups
#  pragma omp for schedule( PAR_OMP_SCHED, PAR_OMP_SCHED_CHUNK )
   for (int PID0=0; PID0<amr->NPatchComma[lv][1]; PID0+=8)
   {
//    1. find the patch groups with target tracer particles
//    --> use patch group as the calculation unit since Prepare_PatchData() only works with patch group
//    --> disadvantage: some patches may not have tracer particles ... (they will be skipped later)
      GotYou = false;

      for (int PID=PID0; PID<PID0+8; PID++)
      {
         if ( amr->patch[0][lv][PID]->NParType[(int)PTYPE_TRACER] > 0 )    GotYou = true;

         if ( GotYou )  break;
      }

//    nothing to do if there are no target tracer particles in the target patch group
      if ( !GotYou )    continue;


//    2. prepare the velocity data for the patch group with particles (need NSIDE_26 for ParGhost>0)
      Prepare_PatchData( lv, TimeNew, VelX, NULL, ParGhost, 1, &PID0, _VELX, _NONE,
                         OPT__FLU_INT_SCHEME, INT_NONE, UNIT_PATCH, NSIDE_26, IntPhase_No,
                         OPT__BC_FLU, BC_POT_NONE, MinDens_No, MinPres_No, MinTemp_No, MinEntr_No, DE_Consistency_No );
      Prepare_PatchData( lv, TimeNew, VelY, NULL, ParGhost, 1, &PID0, _VELY, _NONE,
                         OPT__FLU_INT_SCHEME, INT_NONE, UNIT_PATCH, NSIDE_26, IntPhase_No,
                         OPT__BC_FLU, BC_POT_NONE, MinDens_No, MinPres_No, MinTemp_No, MinEntr_No, DE_Consistency_No );
      Prepare_PatchData( lv, TimeNew, VelZ, NULL, ParGhost, 1, &PID0, _VELZ, _NONE,
                         OPT__FLU_INT_SCHEME, INT_NONE, UNIT_PATCH, NSIDE_26, IntPhase_No,
                         OPT__BC_FLU, BC_POT_NONE, MinDens_No, MinPres_No, MinTemp_No, MinEntr_No, DE_Consistency_No );


      for (int PID=PID0, P=0; PID<PID0+8; PID++, P++)
      {
//       3. compute the particle velocity
//       skip patches with no tracer particles
         if ( amr->patch[0][lv][PID]->NParType[(int)PTYPE_TRACER] == 0 )   continue;

         double EdgeL[3], EdgeR[3];

         for (int d=0; d<3; d++) {
            EdgeL[d] = amr->patch[0][lv][PID]->EdgeL[d] - dh*ParGhost;
            EdgeR[d] = amr->patch[0][lv][PID]->EdgeR[d] + dh*ParGhost;
         }

         for (int p=0; p<amr->patch[0][lv][PID]->NPar; p++)
         {
            ParID = amr->patch[0][lv][PID]->ParList[p];

//          skip massive particles
            if ( ParType[ParID] != PTYPE_TRACER )
               continue;

            if ( MapOnly )
               for (int d=0; d<3; d++)    InterpParPos[d][p] = ParPos[d][ParID];

            else
            {
//             determine time-step
               dt = (real_par)TimeNew - ParTime[ParID];

//             convert time-step for comoving
#              ifdef COMOVING
               if ( ParTime[ParID] == (real_par)TimeOld )    dt = dt_com;   // avoid redundant calculations
               else                                          dt = (real_par)Mis_dTime2dt( (double)ParTime[ParID], (double)dt );
#              endif

//             predict the positions at TimeNew
               for (int d=0; d<3; d++)    InterpParPos[d][p] = ParPos[d][ParID] + dt*ParVel[d][ParID];
            }
         } // for (int p=0; p<amr->patch[0][lv][PID]->NPar; p++)

         Par_MapMesh2Particles( EdgeL, EdgeR, _dh, VelSize, VelX+P*CUBE(VelSize),
                                amr->patch[0][lv][PID]->NPar, InterpParPos, ParType,
                                amr->patch[0][lv][PID]->ParList, UseTracers_Yes, Vel_Temp[0],
                                amr->Par->TracerVelCorr );
         Par_MapMesh2Particles( EdgeL, EdgeR, _dh, VelSize, VelY+P*CUBE(VelSize),
                                amr->patch[0][lv][PID]->NPar, InterpParPos, ParType,
                                amr->patch[0][lv][PID]->ParList, UseTracers_Yes, Vel_Temp[1],
                                amr->Par->TracerVelCorr );
         Par_MapMesh2Particles( EdgeL, EdgeR, _dh, VelSize, VelZ+P*CUBE(VelSize),
                                amr->patch[0][lv][PID]->NPar, InterpParPos, ParType,
                                amr->patch[0][lv][PID]->ParList, UseTracers_Yes, Vel_Temp[2],
                                amr->Par->TracerVelCorr );

//       4. update particles
         for (int p=0; p<amr->patch[0][lv][PID]->NPar; p++)
         {
            ParID = amr->patch[0][lv][PID]->ParList[p];

//          skip massive particles
            if ( ParType[ParID] != PTYPE_TRACER )
               continue;

//          4.1 MapOnly mode
            if ( MapOnly )
            {
               for (int d=0; d<3; d++)
                  ParVel[d][ParID] = Vel_Temp[d][p];
            }

//          4.2 Euler method
            else if ( amr->Par->IntegTracer == TRACER_INTEG_EULER )
            {
               for (int d=0; d<3; d++) {
                  ParPos[d][ParID] = InterpParPos[d][p];
                  ParVel[d][ParID] = Vel_Temp    [d][p];
               }

               ParTime[ParID] = (real_par)TimeNew;
            }

//          4.3 RK2 scheme (position only)
            else if ( amr->Par->IntegTracer == TRACER_INTEG_RK2 )
            {
//             determine time-step
               dt = (real_par)TimeNew - ParTime[ParID];

//             convert time-step for comoving
#              ifdef COMOVING
               if ( ParTime[ParID] == (real_par)TimeOld )    dt = dt_com;   // avoid redundant calculations
               else                                          dt = (real_par)Mis_dTime2dt( (double)ParTime[ParID], (double)dt );
#              endif

               for (int d=0; d<3; d++)
                  InterpParPos[d][p] = ParPos[d][ParID] +
                     (real_par)0.5*dt*( Vel_Temp[d][p] + ParVel[d][ParID] );
            } // amr->Par->IntegTracer

         } // for (int p=0; p<amr->patch[0][lv][PID]->NPar; p++)`


//       5. RK2 scheme (velocity)
         if ( !MapOnly  &&  amr->Par->IntegTracer == TRACER_INTEG_RK2 )
         {
            Par_MapMesh2Particles( EdgeL, EdgeR, _dh, VelSize, VelX+P*CUBE(VelSize),
                                   amr->patch[0][lv][PID]->NPar, InterpParPos, ParType,
                                   amr->patch[0][lv][PID]->ParList, UseTracers_Yes, Vel_Temp[0],
                                   amr->Par->TracerVelCorr );
            Par_MapMesh2Particles( EdgeL, EdgeR, _dh, VelSize, VelY+P*CUBE(VelSize),
                                   amr->patch[0][lv][PID]->NPar, InterpParPos, ParType,
                                   amr->patch[0][lv][PID]->ParList, UseTracers_Yes, Vel_Temp[1],
                                   amr->Par->TracerVelCorr );
            Par_MapMesh2Particles( EdgeL, EdgeR, _dh, VelSize, VelZ+P*CUBE(VelSize),
                                   amr->patch[0][lv][PID]->NPar, InterpParPos, ParType,
                                   amr->patch[0][lv][PID]->ParList, UseTracers_Yes, Vel_Temp[2],
                                   amr->Par->TracerVelCorr );

            for (int p=0; p<amr->patch[0][lv][PID]->NPar; p++) {

               ParID = amr->patch[0][lv][PID]->ParList[p];

//             skip massive particles
               if ( ParType[ParID] != PTYPE_TRACER )
                  continue;

               for (int d=0; d<3; d++) {
                  ParPos[d][ParID] = InterpParPos[d][p];
                  ParVel[d][ParID] = Vel_Temp    [d][p];
               }

               ParTime[ParID] = (real_par)TimeNew;
            } // for (int p=0; p<amr->patch[0][lv][PID]->NPar; p++)
         } // if ( !MapOnly  &&  amr->Par->IntegTracer == TRACER_INTEG_RK2 )
      } // for (int PID=PID0, P=0; PID<PID0+8; PID++, P++)
   } // for (int PID0=0; PID0<amr->NPatchComma[lv][1]; PID0+=8)


// 6. free memory
   delete [] VelX;
   delete [] VelY;
   delete [] VelZ;

   Aux_DeallocateArray2D( Vel_Temp     );
   Aux_DeallocateArray2D( InterpParPos );

   } // end of OpenMP parallel region

} // FUNCTION : Par_UpdateTracerParticle



#endif // #if ( defined PARTICLE  &&  defined TRACER )
