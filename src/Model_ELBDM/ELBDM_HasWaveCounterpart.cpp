#include "GAMER.h"

#if ( ELBDM_SCHEME == ELBDM_HYBRID )




//-------------------------------------------------------------------------------------------------------
// Function    :  ELBDM_HasWaveCounterpart
// Description :  Check whether cell [I, J, K] in patch indexed by GID has wave counterpart on refined levels by traversing the global AMR Tree
//
// Note        :  1. This function requires LB_GlobalPatch* Tree to be initialised beforehand
//
// Parameter   :  I   : x-index relative to patch GID0
//             :  J   : y-index relative to patch GID0
//             :  K   : z-index relative to patch GID0
//             : GID0 : global ID of patch group
//             : GID  : global ID of patch
//             : Tree : pointer to array of LB_GlobalPatch objects
//                      needs to be initialised beforehand by calling LB_GlobalPatch* Tree = LB_GatherTree(pc, MPI_Node);
//
// Return      :  "true"  if cell [I, J, K] in patch GID has    wave counterpart
//                "false" if cell [I, J, K] in patch GID has NO wave counterpart
//-------------------------------------------------------------------------------------------------------
bool ELBDM_HasWaveCounterpart( const int I, const int J, const int K, const long GID0, const long GID, const LB_GlobalTree& GlobalTree )
{

// convert to global coordinates
   const int X = GlobalTree.Local2Global( I, 0, GID0 );
   const int Y = GlobalTree.Local2Global( J, 1, GID0 );
   const int Z = GlobalTree.Local2Global( K, 2, GID0 );

// find GID of child
   const long ChildGID = GlobalTree.FindRefinedCounterpart( X, Y, Z, GID );
   if ( ChildGID == -1 ) return false;

   const bool HasWaveCounterpart = amr->use_wave_flag[ GlobalTree[ChildGID].level ];

   return HasWaveCounterpart;

} // FUNCTION : ELBDM_HasWaveCounterpart



#endif // #if ( ELBDM_SCHEME == ELBDM_HYBRID )
