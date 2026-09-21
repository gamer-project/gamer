#include "GAMER.h"

#if ( defined PARTICLE  &&  defined GRAVITY  &&  defined STORE_PAR_POT )




//-------------------------------------------------------------------------------------------------------
// Function    :  Par_UpdateParticlePotential
// Description :  Interpolate the gravitational potential onto all massive particles at the target level
//                and store the result in Par->Pot[]
//
// Note        :  1. Enabled by the runtime option OPT__OUTPUT_PAR_POT (in addition to compiling with STORE_PAR_POT)
//                   --> Returns immediately if OPT__OUTPUT_PAR_POT is false
//                2. Reuses the generic mesh-to-particle interpolation routine Par_MapMesh2Particles(),
//                   called here with UseTracers=false so that only massive particles are mapped
//                   --> The interpolation order follows amr->Par->Interp (the massive-particle scheme),
//                       and the ghost zones are sized by amr->Par->GhostSize accordingly
//                3. Purely a diagnostic snapshot of the potential at each particle's current position
//                   --> Unlike Par_UpdateParticle(), it does not feed back into the particle integration
//                4. Skips patches/particle groups with no massive particles
//                5. Tracer particles are skipped and keep whatever value Par->Pot[] already holds
//                   (initialized to 0 for every particle; see Particle::InitRepo/AddOneParticle)
//
// Parameter   :  lv       : Target refinement level
//                PrepTime : Target physical time for preparing the potential data
//
// Return      :  amr->Par->Pot[]
//-------------------------------------------------------------------------------------------------------
void Par_UpdateParticlePotential( const int lv, const double PrepTime )
{

   if ( !OPT__OUTPUT_PAR_POT )   return;

   const bool     IntPhase_No        = false;
   const bool     DE_Consistency_No  = false;
   const real     MinDens_No         = -1.0;
   const real     MinPres_No         = -1.0;
   const real     MinTemp_No         = -1.0;
   const real     MinEntr_No         = -1.0;
   const double   dh                 = amr->dh[lv];
   const double   _dh                = 1.0/dh;
   const int      ParGhost           = amr->Par->GhostSize;
   const int      PotSize            = PS1 + 2*ParGhost;
   const bool     UseTracers_No      = false;
   const bool     CorrectVelocity_No = false;

         real_par *ParPos[3] = { amr->Par->PosX, amr->Par->PosY, amr->Par->PosZ };
         real_par *ParPot    = amr->Par->Pot;
   const long_par *ParType   = amr->Par->Type;


// get the maximum number of particles in a single patch
   int NParMax = 0;
   for (int PID=0; PID<amr->NPatchComma[lv][1]; PID++)   NParMax = MAX( NParMax, amr->patch[0][lv][PID]->NPar );


// nothing to do if there is no particle
   if ( NParMax <= 0 )  return;


// OpenMP parallel region
#  pragma omp parallel
   {

// per-thread variables
   real *Pot = new real [ 8*CUBE(PotSize) ];   // 8: number of patches per patch group

   real_par **ParPotTemp  = NULL;
   real_par **InterpParPos = NULL;
   Aux_AllocateArray2D( ParPotTemp,   1, NParMax );
   Aux_AllocateArray2D( InterpParPos, 3, NParMax );

   bool GotYou;
   long ParID;

// loop over all **real** patch groups
#  pragma omp for schedule( PAR_OMP_SCHED, PAR_OMP_SCHED_CHUNK )
   for (int PID0=0; PID0<amr->NPatchComma[lv][1]; PID0+=8)
   {
//    1. find the patch groups with target massive particles
//    --> use patch group as the calculation unit since Prepare_PatchData() only works with patch group
      GotYou = false;

      for (int PID=PID0; PID<PID0+8; PID++)
      {
         if ( amr->patch[0][lv][PID]->NPar - amr->patch[0][lv][PID]->NParType[(int)PTYPE_TRACER] > 0 )   GotYou = true;

         if ( GotYou )  break;
      }

//    nothing to do if there are no target massive particles in the target patch group
      if ( !GotYou )    continue;


//    2. prepare the potential data for the patch group with particles (need NSIDE_26 for ParGhost>0)
      Prepare_PatchData( lv, PrepTime, Pot, NULL, ParGhost, 1, &PID0, _POTE, _NONE,
                         OPT__GRA_INT_SCHEME, INT_NONE, UNIT_PATCH, NSIDE_26, IntPhase_No,
                         OPT__BC_FLU, OPT__BC_POT, MinDens_No, MinPres_No, MinTemp_No, MinEntr_No, DE_Consistency_No );

      for (int PID=PID0, P=0; PID<PID0+8; PID++, P++)
      {
//       skip patches with no massive particles
         if ( amr->patch[0][lv][PID]->NPar - amr->patch[0][lv][PID]->NParType[(int)PTYPE_TRACER] == 0 )   continue;

         double EdgeL[3], EdgeR[3];

         for (int d=0; d<3; d++) {
            EdgeL[d] = amr->patch[0][lv][PID]->EdgeL[d] - dh*ParGhost;
            EdgeR[d] = amr->patch[0][lv][PID]->EdgeR[d] + dh*ParGhost;
         }

//       3. collect the current position of every particle in this patch
//       --> Par_MapMesh2Particles() will skip the tracer ones internally (UseTracers_No)
         for (int p=0; p<amr->patch[0][lv][PID]->NPar; p++)
         {
            ParID = amr->patch[0][lv][PID]->ParList[p];

            for (int d=0; d<3; d++)    InterpParPos[d][p] = ParPos[d][ParID];
         }

         Par_MapMesh2Particles( EdgeL, EdgeR, _dh, PotSize, Pot+P*CUBE(PotSize),
                                amr->patch[0][lv][PID]->NPar, InterpParPos, ParType,
                                amr->patch[0][lv][PID]->ParList, UseTracers_No, ParPotTemp[0],
                                CorrectVelocity_No );

//       4. store the interpolated potential
         for (int p=0; p<amr->patch[0][lv][PID]->NPar; p++)
         {
            ParID = amr->patch[0][lv][PID]->ParList[p];

//          skip tracer particles (Par_MapMesh2Particles() did not fill their ParPotTemp[0][p] slot)
            if ( ParType[ParID] == PTYPE_TRACER )    continue;

            ParPot[ParID] = ParPotTemp[0][p];
         } // for (int p=0; p<amr->patch[0][lv][PID]->NPar; p++)
      } // for (int PID=PID0, P=0; PID<PID0+8; PID++, P++)
   } // for (int PID0=0; PID0<amr->NPatchComma[lv][1]; PID0+=8)


// 5. free memory
   delete [] Pot;

   Aux_DeallocateArray2D( ParPotTemp   );
   Aux_DeallocateArray2D( InterpParPos );

   } // end of OpenMP parallel region

} // FUNCTION : Par_UpdateParticlePotential



#endif // #if ( defined PARTICLE  &&  defined GRAVITY  &&  defined STORE_PAR_POT )
