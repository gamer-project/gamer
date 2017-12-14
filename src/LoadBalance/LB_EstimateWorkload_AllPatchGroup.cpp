#include "GAMER.h"

#ifdef LOAD_BALANCE




//-------------------------------------------------------------------------------------------------------
// Function    :  LB_EstimateWorkload_AllPatchGroup
// Description :  Estimate the workload of all patch groups in this rank
//
// Note        :  1. Workload of each patch is estimated as "1.0 + NParThisPatch*ParWeight/PATCH_SIZE^3",
//                   where "NParThisPatch" is the number of particles in this patch and "ParWeight"
//                   is the relative weighting of particles (see below)
//                   --> Workload of a single patch (without particles) is normalized to 1.0
//                2. Workload of each patch **includes particles in the children patches"
//                   --> For non-leaf patches, this function will collect particles from the leaf patches
//                3. This function assumes that "NPatchTotal[lv]" has already been set by invoking the
//                   function "Mis_GetTotalPatchNumber( lv )"
//
// Parameter   :  lv        : Target refinement level
//                ParWeight : Relative workload weighting of particles
//                            --> <= 0.0 : do not consider particle weighting
//                Load_PG   : Estimated workload of all patch groups in this rank
//                            --> Must be preallocated with th size "amr->NPatchComma[lv][1]/8"
//
// Return      :  Load_PG
//-------------------------------------------------------------------------------------------------------
void LB_EstimateWorkload_AllPatchGroup( const int lv, const double ParWeight, double *Load_PG )
{

// check
   if ( Load_PG == NULL )  Aux_Error( ERROR_INFO, "Load_PG == NULL !!\n" );


// 1. workload of cells --> assuming the weighting of each patch == 1.0
   const int NPG_ThisRank = amr->NPatchComma[lv][1] / 8;

   for (int t=0; t<NPG_ThisRank; t++)  Load_PG[t] = 8.0; // 8 patches per patch group


// 2. workload of particles
#  ifdef PARTICLE
   if ( ParWeight > 0.0 )
   {
//    renormalize the load-balance weighting of one particle so that the weighting of one patch is 1.0
      const double ParWeight_Norm = ParWeight / (double)CUBE(PS1);

//    get the number of particles in each patch
      const bool PredictPos_No     = false;
      const bool SibBufPatch_No    = false;
      const bool FaSibBufPatch_No  = false;
      const bool JustCountNPar_Yes = true;
      const bool TimingSendPar_No  = false;

      Par_CollectParticle2OneLevel( lv, PredictPos_No, NULL_REAL, SibBufPatch_No, FaSibBufPatch_No, JustCountNPar_Yes,
                                    TimingSendPar_No );

      for (int t=0; t<NPG_ThisRank; t++)
      for (int PID=t*8; PID<(t+1)*8; PID++)
      {
         int NParThisPatch;

         if ( amr->patch[0][lv][PID]->son == -1 )  NParThisPatch = amr->patch[0][lv][PID]->NPar;
         else                                      NParThisPatch = amr->patch[0][lv][PID]->NPar_Copy;

#        ifdef DEBUG_PARTICLE
         if ( NParThisPatch < 0 )
            Aux_Error( ERROR_INFO, "NPar (%d) has not been calculated (lv %d, PID %d) !!\n",
                       NParThisPatch, lv, PID );
#        endif

//       add the load-balance weighting of all particles in this patch
         Load_PG[t] += NParThisPatch*ParWeight_Norm;
      } // for t ... PID ...

//    free memory allocated by Par_CollectParticle2OneLevel
      Par_CollectParticle2OneLevel_FreeMemory( lv, SibBufPatch_No, FaSibBufPatch_No );

   } // if ( ParWeight_Norm > 0.0 )
#  endif // #ifdef PARTICLE

} // FUNCTION : LB_EstimateWorkload_AllPatchGroup



#endif // #ifdef LOAD_BALANCE
