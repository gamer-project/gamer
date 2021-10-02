#include "ExtractUniform.h"




//-------------------------------------------------------------------------------------------------------
// Function    :  Buf_SortBoundaryPatch
// Description :  Sort the boundary patches according to their positions recorded in the PosList array
//
// Note        :  Bubble sort
//
// Parameter   :  NPatch   : The number of patches to be sorted
//                IDList   : The patch indices list
//                PosList  : The patch positions list
//-------------------------------------------------------------------------------------------------------
void Buf_SortBoundaryPatch( const int NPatch, int *IDList, int *PosList )
{

   int TempID, TempPos;

   for (int p=0; p<NPatch-1;   p++)
   for (int t=0; t<NPatch-1-p; t++)
   {
      if ( PosList[t] > PosList[t+1] )
      {
         TempPos        = PosList[t];
         TempID         = IDList [t];

         PosList[t]     = PosList[t+1];
         IDList [t]     = IDList [t+1];

         PosList[t+1]   = TempPos;
         IDList [t+1]   = TempID;
      }
   }

}
