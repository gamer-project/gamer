#include "GAMER.h"

#if ( !defined GPU  &&  defined SUPPORT_GRACKLE )




//-------------------------------------------------------------------------------------------------------
// Function    :  Init_MemAllocate_Grackle
// Description :  Allocate the CPU memory for the Grackle solver
//
// Note        :  1. Only work when using CPUs only
//                2. Prepare CHE_NPREP variables
//                   --> CHE_NPREP = 4 currently
//                   --> [mass density, specific internal energy, kinematic energy density, metal density]
//                3. Invoked by Init_MemAllocate()
//                4. Use patches instead of patch groups as the allocation unit
//
// Parameter   :  Che_NPG : Number of patch groups to be evaluated at a time
//-------------------------------------------------------------------------------------------------------
void Init_MemAllocate_Grackle( const int Che_NPG )
{

// nothing to do if Grackle is disabled
   if ( GRACKLE_MODE == GRACKLE_MODE_NONE )  return;


   for (int t=0; t<2; t++)
   {
      h_Che_Array[t] = new real [ (long)CHE_NPREP*(long)Che_NPG*(long)CUBE(PS2) ];
   }

} // FUNCTION : Init_MemAllocate_Grackle



#endif // #if ( !defined GPU  &&  defined SUPPORT_GRACKLE )
