#include "GAMER.h"

#if ( MODEL == ELBDM )


extern double *SelSimHalo_Prof[3];




//-------------------------------------------------------------------------------------------------------
// Function    :  End_TestProb
// Description :  Release memory for the self-similar halo test
//
// Note        :  None
//
// Parameter   :  None
//-------------------------------------------------------------------------------------------------------
void End_TestProb()
{

   for (int v=0; v<3; v++)
   {
      if ( SelSimHalo_Prof[v] != NULL )
      {
         delete [] SelSimHalo_Prof[v];
         SelSimHalo_Prof[v] = NULL;
      }
   }

} // FUNCTION : End_TestProb



#endif // #if ( MODEL == ELBDM )
