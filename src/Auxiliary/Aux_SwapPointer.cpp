#include "GAMER.h"




//-------------------------------------------------------------------------------------------------------
// Function    :  Aux_SwapPointer
// Description :  Swap the two input pointers
//
// Note        :  To swap two pointers ptr1 and ptr2, one should call "Aux_SwapPointer( (void**)&ptr1, (void**)&ptr2 )"
//
// Parameter   :  Ptr1/2   : Pointers to the pointers to be swapped
//-------------------------------------------------------------------------------------------------------
void Aux_SwapPointer( void **Ptr1, void **Ptr2 )
{

   void* Ptr_tmp;

   Ptr_tmp = *Ptr1;
   *Ptr1   = *Ptr2;
   *Ptr2   = Ptr_tmp;

} // FUNCTION : Aux_SwapPointer
