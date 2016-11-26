#include "GAMER.h"

#ifdef PARTICLE


extern double *ClusterMerger_Prof[3];




//-------------------------------------------------------------------------------------------------------
// Function    :  End_TestProb
// Description :  Release memory for the N particles force test
//
// Note        :  None
//
// Parameter   :  None 
//-------------------------------------------------------------------------------------------------------
void End_TestProb()
{  

   for (int v=0; v<3; v++)    delete [] ClusterMerger_Prof[v];

} // FUNCTION : End_TestProb



#endif // #ifdef PARTICLE
