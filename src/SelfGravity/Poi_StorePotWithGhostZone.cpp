#include "GAMER.h"

#if ( defined GRAVITY  &&  defined STORE_POT_GHOST )




//-------------------------------------------------------------------------------------------------------
// Function    :  Poi_StorePotWithGhostZone
// Description :  Fill up the potential array pot_ext[] including ghost zones for each target patch
//
// Note        :  1. Called by Gra_AdvancedDt() after the base-level FFT solver, EvolveLevel() after grid
//                   refinement, and Flu_CorrAfterAllSync() when OPT__CORR_AFTER_ALL_SYNC is enabled
//                2. For potential at lv>0, pot_ext[] is filled by Poi_Close() directly (except after grid
//                   refinement), and thus NO need to call this function for that.
//                3. After grid refinement, only newly-allocated patches need to set pot_ext[]
//                   --> We set pot_ext[0][0][0] == POT_EXT_NEED_INIT for newly-allocated patches
//                       so as to distinguish them from other existing patches
//                4. Currently pot_ext[] is only used for Par->ImproveAcc and star formation routines
//                5. Prepare_PatchData() will use extrapolation to fill the ghost zones outside the
//                   non-periodic boundaries
//                6. Currently the number of ghost zones is set by GRA_GHOST_SIZE
//                7. Only apply to real patches
//
// Parameter   :  lv       : Target refinement level
//                PotSg    : Target potential sandglass
//                AllPatch : true  --> work on all real patches at lv
//                                     (used after the root-level FFT solver and Flu_CorrAfterAllSync())
//                           false --> only work on real patches with pot_ext[0][0][0] == POT_EXT_NEED_INIT
//                                     (used after grid refinement)
//-------------------------------------------------------------------------------------------------------
void Poi_StorePotWithGhostZone( const int lv, const int PotSg, const bool AllPatch )
{

   const bool   IntPhase_No       = false;
   const bool   DE_Consistency_No = false;
   const real   MinDens_No        = -1.0;
   const real   MinPres_No        = -1.0;
   const double PrepPotTime       = amr->PotSgTime[lv][PotSg];
   const int    PotGhost          = GRA_GHOST_SIZE;
   const int    PotSize           = PS1 + 2*PotGhost;
   const int    PotSizeCube       = CUBE(PotSize);

// OpenMP parallel region
#  pragma omp parallel
   {
//    per-thread variables
      real *Pot = new real [ 8*PotSizeCube ];   // 8: number of patches per patch group

#     pragma omp for schedule( runtime )
      for (int PID0=0; PID0<amr->NPatchComma[lv][1]; PID0+=8)
      {
         if ( AllPatch  ||  amr->patch[PotSg][lv][PID0]->pot_ext[0][0][0] == POT_EXT_NEED_INIT )
         {
            Prepare_PatchData( lv, PrepPotTime, Pot, NULL, PotGhost, 1, &PID0, _POTE, _NONE,
                               OPT__REF_POT_INT_SCHEME, INT_NONE, UNIT_PATCH, NSIDE_26, IntPhase_No,
                               OPT__BC_FLU, OPT__BC_POT, MinDens_No, MinPres_No, DE_Consistency_No );

            for (int PID=PID0, P=0; PID<PID0+8; PID++, P++)
               memcpy( amr->patch[PotSg][lv][PID]->pot_ext, Pot+P*PotSizeCube, PotSizeCube*sizeof(real) );
         }
      }

      delete [] Pot;
   } // end of OpenMP parallel region

} // FUNCTION : Poi_StorePotWithGhostZone



#endif // #if ( defined GRAVITY  &&  defined STORE_POT_GHOST )
