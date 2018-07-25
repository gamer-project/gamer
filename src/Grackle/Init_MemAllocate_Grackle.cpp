#include "GAMER.h"

#ifdef SUPPORT_GRACKLE


// global variables for accessing h_Che_Array[]
// --> also used by Grackle_Prepare.cpp and Grackle_Close.cpp
//int Che_NField = NULL_INT;
int Che_NField = 4;




//-------------------------------------------------------------------------------------------------------
// Function    :  Init_MemAllocate_Grackle
// Description :  Allocate the CPU memory for the Grackle solver
//
// Note        :  1. Work even when GPU is enabled
//                2. Invoked by Init_MemAllocate()
//
// Parameter   :  Che_NPG : Number of patch groups to be evaluated at a time
//-------------------------------------------------------------------------------------------------------
void Init_MemAllocate_Grackle( const int Che_NPG )
{

// nothing to do if Grackle is disabled
   if ( !GRACKLE_ACTIVATE )   return;


// set Che_NField ...


// allocate the input/output array for the Grackle solver
   for (int t=0; t<2; t++)
      h_Che_Array[t] = new real [ (long)Che_NField*(long)Che_NPG*(long)CUBE(PS2) ];

} // FUNCTION : Init_MemAllocate_Grackle



#endif // #ifdef SUPPORT_GRACKLE
