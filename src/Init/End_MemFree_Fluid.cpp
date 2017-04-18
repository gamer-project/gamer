#include "Copyright.h"
#ifndef GPU

#include "GAMER.h"




//-------------------------------------------------------------------------------------------------------
// Function    :  End_MemFree_Fluid
// Description :  Free memory previously allocated by the function "Init_MemAllocate_Fluid"
//
// Parameter   :  None
//-------------------------------------------------------------------------------------------------------
void End_MemFree_Fluid()
{

   for (int t=0; t<2; t++)
   {
      if ( h_Flu_Array_F_In       [t] != NULL )    delete [] h_Flu_Array_F_In       [t];
      if ( h_Flu_Array_F_Out      [t] != NULL )    delete [] h_Flu_Array_F_Out      [t];
      if ( h_Flux_Array           [t] != NULL )    delete [] h_Flux_Array           [t];
#     ifdef UNSPLIT_GRAVITY
      if ( h_Corner_Array_F       [t] != NULL )    delete [] h_Corner_Array_F       [t];
#     endif
      if ( h_MinDtInfo_Fluid_Array[t] != NULL )    delete [] h_MinDtInfo_Fluid_Array[t];
#     ifdef DUAL_ENERGY
      if ( h_DE_Array_F_Out       [t] != NULL )    delete [] h_DE_Array_F_Out       [t];
#     endif

      h_Flu_Array_F_In       [t] = NULL;
      h_Flu_Array_F_Out      [t] = NULL;
      h_Flux_Array           [t] = NULL;
#     ifdef UNSPLIT_GRAVITY
      h_Corner_Array_F       [t] = NULL;
#     endif
      h_MinDtInfo_Fluid_Array[t] = NULL;
#     ifdef DUAL_ENERGY
      h_DE_Array_F_Out       [t] = NULL;
#     endif
   }

} // FUNCTION : End_MemFree_Fluid



#endif // #ifndef GPU
