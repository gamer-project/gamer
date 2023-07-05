#include "GAMER.h"

#if ( MODEL == ELBDM && ELBDM_SCHEME == ELBDM_HYBRID )


//-------------------------------------------------------------------------------------------------------
// Function    :  ELBDM_HasWaveCounterpart
// Description :  Check whether cell [I, J, K] in patch indexed by GID GID has wave counterpart on refined levels by traversing the global AMR Tree
//
// Note        :  1. This function requires LB_GlobalPatch* Tree to be initialised beforehand
//
// Parameter   :  I   : x-index of patch GID
//             :  J   : y-index of patch GID
//             :  K   : z-index of patch GID
//             : GID0 : global ID of patch group
//             : GID  : global ID of patch
//             : Tree : pointer to array of LB_GlobalPatch objects
//                      needs to initialised beforehand by calling LB_GlobalPatch* Tree = LB_GatherTree(pc, MPI_Node);
//
// Return      :  "true"  if cell [I, J, K] in patch GID has    wave counterpart
//                "false" if cell [I, J, K] in patch GID has NO wave counterpart
//-------------------------------------------------------------------------------------------------------
bool ELBDM_HasWaveCounterpart(int I, int J, int K, long GID0, long GID, LB_GlobalPatch* Tree)
{
   long ChildGID = ELBDM_FindRefinedGID(I, J, K, GID0, GID, Tree);
   if ( ChildGID == -1 ) return false;

   bool HasWaveCounterpart = amr->use_wave_flag[Tree[ChildGID].level];

   return HasWaveCounterpart;
} // FUNCTION : ELBDM_HasWaveCounterpart

#endif // #if ( MODEL == ELBDM && ELBDM_SCHEME == ELBDM_HYBRID )