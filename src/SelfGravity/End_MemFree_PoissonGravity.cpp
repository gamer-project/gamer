#include "GAMER.h"

#if ( !defined GPU  &&  defined GRAVITY )




//-------------------------------------------------------------------------------------------------------
// Function    :  End_MemFree_PoissonGravity
// Description :  Free memory previously allocated by Init_MemAllocate_PoissonGravity()
//
// Parameter   :  None
//-------------------------------------------------------------------------------------------------------
void End_MemFree_PoissonGravity()
{

   for (int t=0; t<2; t++)
   {
      delete [] h_Rho_Array_P    [t];  h_Rho_Array_P    [t] = NULL;
      delete [] h_Pot_Array_P_In [t];  h_Pot_Array_P_In [t] = NULL;
      delete [] h_Pot_Array_P_Out[t];  h_Pot_Array_P_Out[t] = NULL;
#     ifdef UNSPLIT_GRAVITY
      delete [] h_Pot_Array_USG_G[t];  h_Pot_Array_USG_G[t] = NULL;
      delete [] h_Flu_Array_USG_G[t];  h_Flu_Array_USG_G[t] = NULL;
#     endif
      delete [] h_Flu_Array_G    [t];  h_Flu_Array_G    [t] = NULL;
      delete [] h_Corner_Array_G [t];  h_Corner_Array_G [t] = NULL;
#     ifdef DUAL_ENERGY
      delete [] h_DE_Array_G     [t];  h_DE_Array_G     [t] = NULL;
#     endif
#     ifdef MHD
      delete [] h_EngyB_Array_G  [t];  h_EngyB_Array_G  [t] = NULL;
#     endif
      delete [] h_Pot_Array_T    [t];  h_Pot_Array_T    [t] = NULL;
   }

   delete [] GreenFuncK;   GreenFuncK = NULL;

} // FUNCTION : End_MemFree_PoissonGravity



#endif // #if ( !defined GPU  &&  defined GRAVITY )
