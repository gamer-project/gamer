#include "GAMER.h"

#if ( MODEL == ELBDM && ELBDM_SCHEME == ELBDM_HYBRID )


// check whether cell {I, J, K} is within given patch
bool IsInsidePatch(int Coordinates[3], long GID, LB_GlobalPatch* Tree)
{

   bool isInside = true;

   for ( int l = 0; l < 3; ++l ) {
      if (Coordinates[l] < Tree[GID].corner[l] || Coordinates[l] >=  Tree[GID].corner[l] + PS1 * amr->scale[Tree[GID].level])
         isInside = false;
   }

   return isInside;
}

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
//             : GID  : global ID of patch within patch group
//             : Tree : pointer to array of LB_GlobalPatch objects
//                      needs to initialised beforehand by calling LB_GlobalPatch* Tree = LB_GatherTree(pc, MPI_Node);
//
// Return      :  GID
//-------------------------------------------------------------------------------------------------------
long ELBDM_FindRefinedGID(int I, int J, int K, long GID0, long GID, LB_GlobalPatch* Tree)
{


#  ifdef GAMER_DEBUG
// sanity check
   if (GID0 == -1)
   {
      Aux_Error(ERROR_INFO, "GID0 == -1 !!\n");
   }
   if (GID == -1)
   {
      Aux_Error(ERROR_INFO, "GID  == -1!\n");
   }
#  endif

// convert coordinates of cell to global integer coordinate system
   int Coordinates[3] = {I, J, K};
   for ( int l = 0; l < 3; ++l )
      Coordinates[l] = Tree[GID0].corner[l] + Coordinates[l] * amr->scale[Tree[GID].level];

// skip calculation if coordinates are not inside patch GID
   if ( !IsInsidePatch(Coordinates, GID, Tree) ) {
      return -1;
   }

   long FaGID  = GID;
   long SonGID = Tree[FaGID].son;

   printf("I %d J %d K %d %ld %ld\n", I, J, K, FaGID, SonGID);

// traverse the tree up until leave nodes
   while ( SonGID != -1 ) {
//    loop over patches in patch block and check if cell {K, J, I} belongs to patch
      for (int LocalPID = 0; LocalPID < 8; LocalPID++)
      {

         if ( IsInsidePatch(Coordinates, SonGID + LocalPID, Tree) )
         {
            FaGID   = SonGID + LocalPID;
            SonGID  = Tree[FaGID].son;
            break;

         }
#        ifdef GAMER_DEBUG
         else if ( LocalPID == 7 )
         {
            Aux_Error(ERROR_INFO, "Global coordinates {%d, %d, %d} in father patch (GID = %ld, lv = %d), but not in any son patch (GID = %ld, lv = %d)!!\n",
            Coordinates[0], Coordinates[1], Coordinates[2], FaGID, Tree[FaGID].level, SonGID, Tree[SonGID].level);
         }
#        endif
      }
   }

#  ifdef GAMER_DEBUG
// sanity check
   if (SonGID != -1)
   {
      Aux_Error(ERROR_INFO, "SonGID == -1!!\n");
   }
   if (FaGID == -1)
   {
      Aux_Error(ERROR_INFO, "FaGID == -1!!\n");
   }

// check whether GID is ancestor of FaGID
// number of generations between GID and FaGID
   int  NGenerations = Tree[FaGID].level - Tree[GID].level;

// iterate back through ancestors
   long AncestorGID = FaGID;
   for (int i = 0; i < NGenerations; ++i)
   {
      AncestorGID = Tree[AncestorGID].father;
   }

   if (Tree[AncestorGID].level != Tree[GID].level)
   {
      Aux_Error(ERROR_INFO, "GID (%ld) and Ancestor (%ld) on different levels!!\n", GID, AncestorGID);
   }

   if (AncestorGID != GID)
   {
      Aux_Error(ERROR_INFO, "GID (GID = %ld, lv = %d) and Ancestor (%ld) are not related!!\n", GID, Tree[GID].level, AncestorGID);
   }

#  endif

   return FaGID;
} // FUNCTION : ELBDM_FindLevel

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