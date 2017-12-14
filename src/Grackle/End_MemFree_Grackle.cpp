#include "GAMER.h"

#if ( !defined GPU  &&  defined SUPPORT_GRACKLE )




//-------------------------------------------------------------------------------------------------------
// Function    :  End_MemFree_Grackle
// Description :  Free memory previously allocated by Init_MemAllocate_Grackle()
//
// Note        :  1. Only work when using CPUs only
//                2. Invoked by End_MemFree()
//
// Parameter   :  None
//-------------------------------------------------------------------------------------------------------
void End_MemFree_Grackle()
{

// nothing to do if Grackle is disabled
   if ( GRACKLE_MODE == GRACKLE_MODE_NONE )  return;


   for (int t=0; t<2; t++)
   {
      if ( h_Che_Array[t] != NULL )    delete [] h_Che_Array[t];
      h_Che_Array[t] = NULL;
   }

} // FUNCTION : End_MemFree_Grackle



#endif // #if ( !defined GPU  &&  defined SUPPORT_GRACKLE )
