#ifndef __GATHER_TREE_H__
#define __GATHER_TREE_H__



#include <stdio.h>
#include "Macro.h"



class NonCopyable
{
  protected:
    NonCopyable() {}
    ~NonCopyable() {}
  private:
    NonCopyable( const NonCopyable & );
    NonCopyable& operator = ( const NonCopyable & );
}; // class NonCopyable




struct LB_GlobalPatch
{
   int    corner[3];
   int    sibling[26];
   int    father;
   int    son;
   long   LB_Idx;
   int    level;
   double EdgeL[3];
   double EdgeR[3];
   ulong  PaddedCr1D;
   int    MPI_Rank;

#  ifdef PARTICLE
   int    NPar;
#  endif
}; // struct LB_GlobalPatch



struct LB_PatchCount : private NonCopyable
{
   LB_PatchCount();
   ~LB_PatchCount();

   long NPatchAllLv;
   int  NPatchLocalAllLv;
   int  NPatchLocal[NLEVEL];      // number of patches per level on MPI_Rank
   int  (*NPatchAllRank)[NLEVEL]; // number of patches in [MPI rank][level]
   int  GID_Offset[NLEVEL];       // offsets that can be used to convert local PID at level lv to GID via GID = PID + GID_Offset[lv]
   int  GID_LvStart[NLEVEL];      // global patch index at which level starts

   bool isInitialised;
}; // struct LB_PatchCount



// allocate memory and store pointers to lists with local patch information
struct LB_LocalPatchExchangeList : private NonCopyable
{
   LB_LocalPatchExchangeList();
   ~LB_LocalPatchExchangeList();

   long    *LBIdxList_Local     [NLEVEL];      // load balance ids
   int    (*CrList_Local        [NLEVEL])[3];  // patch corners
   int     *FaList_Local        [NLEVEL];      // father GIDs
   int     *SonList_Local       [NLEVEL];      // son GIDs
   int    (*SibList_Local       [NLEVEL])[26]; // sibling GIDs
   double (*EdgeLList_Local     [NLEVEL])[3];  // left edge of the patch
   double (*EdgeRList_Local     [NLEVEL])[3];  // right edge of the patch
   ulong   *PaddedCr1DList_Local[NLEVEL];      // 1D corner coordinate padded with two base-level patches on each side
   int     *MPI_RankList_Local  [NLEVEL];      // MPI rank where patch resides
#  ifdef PARTICLE
   int     *NParList_Local      [NLEVEL];      // particle GIDs
#  endif

   bool  isInitialised;
   long *LBIdxList_Sort         [NLEVEL];
   int  *LBIdxList_Sort_IdxTable[NLEVEL];
   bool  LBIdxisInitialised;
}; // struct LB_LocalPatchExchangeList



// allocate memory and store pointers to lists with global patch information
struct LB_GlobalPatchExchangeList : private NonCopyable
{
   LB_GlobalPatchExchangeList( LB_PatchCount& pc, int root );
   ~LB_GlobalPatchExchangeList();

   long    *LBIdxList_AllLv;            // load balance ids
   int    (*CrList_AllLv         )[3];  // patch corners
   int     *FaList_AllLv;               // father GIDs
   int     *SonList_AllLv;              // son GIDs
   int    (*SibList_AllLv        )[26]; // sibling GIDs
   double (*EdgeLList_AllLv      )[3];  // left edge of the patch
   double (*EdgeRList_AllLv      )[3];  // right edge of the patch
   ulong   *PaddedCr1DList_AllLv;       // 1D corner coordinate padded with two base-level patches on each side
   int     *MPI_RankList_AllLv;         // MPI rank where patch resides
#  ifdef PARTICLE
   int     *NParList_AllLv;             // particle numbers
#  endif
   bool isAllocated;
   bool isInitialised;
}; // struct LB_GlobalPatchExchangeList



#endif // #ifndef __GATHER_TREE_H__
