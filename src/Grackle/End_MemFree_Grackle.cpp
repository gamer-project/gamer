#include "GAMER.h"

#ifdef SUPPORT_GRACKLE




//-------------------------------------------------------------------------------------------------------
// Function    :  End_MemFree_Grackle
// Description :  Free memory previously allocated by Init_MemAllocate_Grackle()
//
// Note        :  1. Work even when GPU is enabled
//                2. Invoked by End_MemFree()
//
// Parameter   :  None
//-------------------------------------------------------------------------------------------------------
void End_MemFree_Grackle()
{

// nothing to do if Grackle is disabled
   if ( !GRACKLE_ACTIVATE )   return;


   for (int t=0; t<2; t++)
   {
      if ( h_Che_Array[t] != NULL )    delete [] h_Che_Array[t];
      h_Che_Array[t] = NULL;
   }

} // FUNCTION : End_MemFree_Grackle



#endif // #ifdef SUPPORT_GRACKLE
