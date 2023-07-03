#include "GAMER.h"

#if ( MODEL == ELBDM && ELBDM_SCHEME == ELBDM_HYBRID )


// check whether cell {I, J, K} is within given patch
bool IsInsidePatch(int Coordinates[3], long GID, LB_GlobalPatch* tree)
{

   bool isInside = true;

   for ( int l = 0; l < 3; ++l ) {
      if (Coordinates[l] < tree[GID].corner[l] || Coordinates[l] >=  tree[GID].corner[l] + PS1 * amr->scale[tree[GID].level])
         isInside = false;
   }

   return isInside;
}

//-------------------------------------------------------------------------------------------------------
// Function    :  ELBDM_HasWaveCounterpart
// Description :  Check whether cell [I, J, K] in patch indexed by GID GID has wave counterpart on refined levels by traversing the global AMR tree
//
// Note        :  1. This function requires LB_GlobalPatch* tree to be initialised beforehand
//
// Parameter   :  I   : x-index of patch GID
//             :  J   : y-index of patch GID
//             :  K   : z-index of patch GID
//             : GID  : global ID of patch
//             : GID  : global ID of patch
//             : tree : pointer to array of LB_GlobalPatch objects
//                      needs to initialised beforehand by calling LB_GlobalPatch* tree = LB_GatherTree(pc, MPI_Node);
//
// Return      :  "true"  if cell [I, J, K] in patch GID has    wave counterpart
//                "false" if cell [I, J, K] in patch GID has NO wave counterpart
//-------------------------------------------------------------------------------------------------------
bool ELBDM_HasWaveCounterpart(int I, int J, int K, long GID0, long GID, LB_GlobalPatch* tree)
{
   int lv = tree[GID].level;

// convert coordinates of cell to global integer coordinate system
   int Coordinates[3] = {I, J, K};
   for ( int l = 0; l < 3; ++l )
      Coordinates[l] = tree[GID0].corner[l] + Coordinates[l] * amr->scale[lv];

// skip calculation
   if ( !IsInsidePatch(Coordinates, GID, tree) ) {
      return false;
   } else {
      if ( amr->use_wave_flag[lv] ) return true;
   }

   int  currentLv  = lv;
   long currentGID = GID;
   long sonGID     = tree[currentGID].son;

// traverse the tree up until leave nodes
   while ( sonGID != -1 ) {
      currentLv += 1;

//    loop over patches in patch block or sons and check if cell {K, J, I} belongs to patch
      for (long currentGID=sonGID; currentGID<sonGID+8; currentGID++)
      {

//       if it is within the given patch,
//       if not iterate over son's sons
         if ( IsInsidePatch(Coordinates, currentGID, tree) ) {
//          first check if son is on wave level
            if (amr->use_wave_flag[currentLv]) return true;
            else
            {
               sonGID     = tree[currentGID].son;
               break;
            }
         } else if ( currentGID == (sonGID + 7) )
         {
            return false;
         }
      }
   }

   return false;
} // FUNCTION : ELBDM_HasWaveCounterpart

#endif // #if ( MODEL == ELBDM && ELBDM_SCHEME == ELBDM_HYBRID )