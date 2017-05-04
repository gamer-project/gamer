#include "GAMER.h"

extern double *AGORA_VcProf[2];




//-------------------------------------------------------------------------------------------------------
// Function    :  End_TestProb
// Description :  Release memory for test problem
//
// Note        :  None
//
// Parameter   :  None
//-------------------------------------------------------------------------------------------------------
void End_TestProb()
{

   for (int v=0; v<2; v++)
   {
      delete [] AGORA_VcProf[v];
      AGORA_VcProf[v] = NULL;
   }

} // FUNCTION : End_TestProb
