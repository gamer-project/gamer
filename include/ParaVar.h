#ifndef __PARAVAR_H__
#define __PARAVAR_H__



#include "Macro.h"




//-------------------------------------------------------------------------------------------------------
// Structure   :  ParaVar_t
// Description :  Data structure collecting variables related to the GAMER parallelization
//                (mainly for the load-imbalance implementation)
//
// Data Member :  BounP_NList       : Number of boundary patches in 26 sibling directions
//                BounP_IDList      : IDs of boundary patches in 26 sibling directions
//                BounP_PosList     : Positions of boundary patches recorded in "BounP_IDList"
//                SendP_NList       : Number of boundary patches to be sent to 26 sibling ranks
//                SendP_IDList      : IDs of boundary patches to be sent to 26 sibling ranks
//                RecvP_NList       : Number of buffer patches to receive data from 26 sibling ranks
//                RecvP_IDList      : IDs of buffer patches to receive data from 26 sibling ranks
//                BounFlag_NList    : Number of flag information to be sent to 26 sibling ranks
//                BounFlag_PostList : Position of each flag recorded in the array "BounFlag_NList"
//                BuffFlag_NList    : Number of flag information to receive from 26 sibling ranks
//                BuffFlag_PostList : Position of each flag recorded in the array "BuffFlag_NList"
//                SendF_NList       : Number of patches to send flux data to 6 sibling ranks
//                SendF_IDList      : IDs of patches to send flux data to 6 sibling ranks
//                RecvF_NList       : Number of patches to receive flux data from 6 sibling ranks
//                RecvF_IDList      : IDs of patches to receive fluix data from 6 sibling ranks
//                SubDomain_EdgeL   : Left  edge of the sub-domain
//                SubDomain_EdgeR   : Right edge of the sub-domain
//
// Method      :  ParaVar_t   : Constructor
//               ~ParaVar_t   : Destructor
//                Lvdelete    : Delete and reset all variables in the given level
//-------------------------------------------------------------------------------------------------------
struct ParaVar_t
{

// data members
// ===================================================================================
   int  BounP_NList      [NLEVEL  ][26];
   int *BounP_IDList     [NLEVEL  ][26];
   int *BounP_PosList    [NLEVEL  ][26];

   int  SendP_NList      [NLEVEL  ][26];
   int *SendP_IDList     [NLEVEL  ][26];
   int  RecvP_NList      [NLEVEL  ][26];
   int *RecvP_IDList     [NLEVEL  ][26];

   int  BounFlag_NList   [NLEVEL  ][26];
   int *BounFlag_PosList [NLEVEL  ][26];
   int  BuffFlag_NList   [NLEVEL  ][26];
   int *BuffFlag_PosList [NLEVEL  ][26];

   int  SendF_NList      [NLEVEL-1][ 6];
   int *SendF_IDList     [NLEVEL-1][ 6];
   int  RecvF_NList      [NLEVEL-1][ 6];
   int *RecvF_IDList     [NLEVEL-1][ 6];

   double SubDomain_EdgeL[3];
   double SubDomain_EdgeR[3];



   //===================================================================================
   // Constructor :  ParaVar_t
   // Description :  Constructor of the structure "ParaVar_t"
   //
   // Note        :  Initialize all pointers as NULL and all counters as 0
   //===================================================================================
   ParaVar_t()
   {
      for (int lv=0; lv<NLEVEL; lv++)
      for (int s=0; s<26; s++)
      {
         BounP_IDList      [lv][s] = NULL;
         BounP_PosList     [lv][s] = NULL;
         SendP_IDList      [lv][s] = NULL;
         RecvP_IDList      [lv][s] = NULL;
         BounFlag_PosList  [lv][s] = NULL;
         BuffFlag_PosList  [lv][s] = NULL;

         BounP_NList       [lv][s] = 0;
         SendP_NList       [lv][s] = 0;
         RecvP_NList       [lv][s] = 0;
         BounFlag_NList    [lv][s] = 0;
         BuffFlag_NList    [lv][s] = 0;
      }

      for (int lv=0; lv<NLEVEL-1; lv++)
      for (int s=0; s<6; s++)
      {
         SendF_IDList      [lv][s] = NULL;
         RecvF_IDList      [lv][s] = NULL;

         RecvF_NList       [lv][s] = 0;
         SendF_NList       [lv][s] = 0;
      }

      for (int d=0; d<3; d++)
      {
         SubDomain_EdgeL[d] = NULL_REAL;
         SubDomain_EdgeR[d] = NULL_REAL;
      }
   } // METHOD : ParaVar_t



   //===================================================================================
   // Constructor :  ~ParaVar_t
   // Description :  Destructor of the structure "ParaVar_t"
   //
   // Note        :  Deallocate memory previously allocated and reset all XXX_NList
   //===================================================================================
   ~ParaVar_t()
   {
      for (int lv=0; lv<NLEVEL; lv++)     Lvdelete( lv );
   } // METHOD : ~ParaVar_t



   //===================================================================================
   // Constructor :  Lvdelete
   // Description :  Deallocate memory previously allocated and reset all XXX_NList
   //                for the target refinement level
   //===================================================================================
   void Lvdelete( const int lv )
   {
      for (int s=0; s<26; s++)
      {
         if ( BounP_IDList    [lv][s] != NULL )    delete [] BounP_IDList    [lv][s];
         if ( BounP_PosList   [lv][s] != NULL )    delete [] BounP_PosList   [lv][s];
         if ( SendP_IDList    [lv][s] != NULL )    delete [] SendP_IDList    [lv][s];
         if ( RecvP_IDList    [lv][s] != NULL )    delete [] RecvP_IDList    [lv][s];
         if ( BounFlag_PosList[lv][s] != NULL )    delete [] BounFlag_PosList[lv][s];
         if ( BuffFlag_PosList[lv][s] != NULL )    delete [] BuffFlag_PosList[lv][s];

         BounP_IDList    [lv][s] = NULL;
         BounP_PosList   [lv][s] = NULL;
         SendP_IDList    [lv][s] = NULL;
         RecvP_IDList    [lv][s] = NULL;
         BounFlag_PosList[lv][s] = NULL;
         BuffFlag_PosList[lv][s] = NULL;

         BounP_NList     [lv][s] = 0;
         SendP_NList     [lv][s] = 0;
         RecvP_NList     [lv][s] = 0;
         BounFlag_NList  [lv][s] = 0;
         BuffFlag_NList  [lv][s] = 0;
      }

      if ( lv != NLEVEL-1 )
      {
         for (int s=0; s<6; s++)
         {
            if ( SendF_IDList[lv][s] != NULL )     delete [] SendF_IDList[lv][s];
            if ( RecvF_IDList[lv][s] != NULL )     delete [] RecvF_IDList[lv][s];

            SendF_IDList[lv][s] = NULL;
            RecvF_IDList[lv][s] = NULL;

            SendF_NList [lv][s] = 0;
            RecvF_NList [lv][s] = 0;
         }
      }
   } // METHOD : Lvdelete


}; // struct ParaVar_t



#endif // #ifndef __PARAVAR_H__
