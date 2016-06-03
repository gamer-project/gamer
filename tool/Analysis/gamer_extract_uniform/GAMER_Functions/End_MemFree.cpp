#include "ExtractUniform.h"




//-------------------------------------------------------------------------------------------------------
// Function    :  End_MemFree
// Description :  Free memory allocated by the function "Init_MemoryAllocate"
//-------------------------------------------------------------------------------------------------------
void End_MemFree()
{

// a. deallocate parallelization variables
   delete [] BaseP;

   for (int lv=0; lv<NLEVEL; lv++)
   for (int s=0; s<26; s++)
   {
      if ( ParaVar.BounP_IDList     [lv][s] != NULL )    delete [] ParaVar.BounP_IDList      [lv][s];
      if ( ParaVar.BounP_PosList    [lv][s] != NULL )    delete [] ParaVar.BounP_PosList     [lv][s];
      if ( ParaVar.SendP_IDList     [lv][s] != NULL )    delete [] ParaVar.SendP_IDList      [lv][s];
      if ( ParaVar.RecvP_IDList     [lv][s] != NULL )    delete [] ParaVar.RecvP_IDList      [lv][s];
   }


// b. deallocate patches
   for (int lv=NLEVEL-1; lv>=0; lv--)
   {
      for (int PID=0; PID<amr.num[lv]; PID++)
      {
         delete amr.patch[lv][PID];

         amr.patch[lv][PID] = NULL;
      }

      amr.num[lv] = 0;
   }

}



