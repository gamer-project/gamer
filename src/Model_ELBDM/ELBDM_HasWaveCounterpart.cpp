#include "GAMER.h"

#if ( MODEL == ELBDM && ELBDM_SCHEME == ELBDM_HYBRID )

//-------------------------------------------------------------------------------------------------------
// Function    :  ELBDM_HasWaveCounterpart
// Description :  Check whether cell [I, J, K] in patch indexed by GID GID0 has wave counterpart on refined levels by traversing the global AMR tree
//
// Note        :  1. This function requires LB_GlobalPatch* tree to be initialised beforehand
//
// Parameter   :  I   : x-index of patch GID0
//             :  J   : y-index of patch GID0
//             :  K   : z-index of patch GID0
//             : GID0 : global ID of patch
//             : tree : pointer to array of LB_GlobalPatch objects
//                      needs to initialised beforehand by calling LB_GlobalPatch* tree = LB_GatherTree(pc, MPI_Node);
//
// Return      :  "true"  if cell [I, J, K] in patch GID0 has    wave counterpart
//                "false" if cell [I, J, K] in patch GID0 has NO wave counterpart
//-------------------------------------------------------------------------------------------------------
bool ELBDM_HasWaveCounterpart(int I, int J, int K, long GID0, LB_GlobalPatch* tree)
{
   int lv = tree[GID0].level;

// convert Coordinates of cell to global integer coordinate system
   int Coordinates[3] = {I, J, K};
   for ( int l = 0; l < 3; ++l )
      Coordinates[l] = tree[GID0].corner[l] + Coordinates[l] * amr->scale[lv];

// during first iteration, we can cover the eight patches in the patch group
// imagine that patches in patch group are children of a common father even at level 0
// we set currentLv = 0 in loop and iterate over patch group
// once we find the patch that the coordinate is in, we iterate over that patch's children

   int  currentLv  = lv - 1;
   bool isInside   = true;

   long currentGID = GID0;
   long sonGID     = GID0;

// traverse the tree until we get to leave nodes
   while ( sonGID != -1 ) {
      currentLv += 1;

//    loop over patches in patch block or sons and check which patch the cell {K, J, I} belongs to
      for (int GID=sonGID; GID<sonGID+8; GID++)
      {
         isInside = true;
//       check whether cell {I, J, K} is within g,iven patch
         for ( int l = 0; l < 3; ++l ) {
            if (Coordinates[l] < tree[GID].corner[l] || Coordinates[l] >=  tree[GID].corner[l] + PS1 * amr->scale[currentLv])
               isInside = false;
         }

//       if it is within the given patch, go to the patch's son
         if ( isInside ) {
            currentGID = GID;
            sonGID     = tree[GID].son;
            break;
         }
#        ifdef GAMER_DEBUG
         if ( !isInside && (GID == (sonGID + 7)) ) {
            Aux_Error(ERROR_INFO, "Coordinates in ELBDM_HasWaveCounterpart not in correct patch block!! \n");
         }
#        endif
      }
   }

   return amr->use_wave_flag[currentLv];
} // FUNCTION : ELBDM_HasWaveCounterpart

#endif // #if ( MODEL == ELBDM && ELBDM_SCHEME == ELBDM_HYBRID )