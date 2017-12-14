#include "GAMER.h"

#if ( !defined GPU  &&  defined GRAVITY )




//-------------------------------------------------------------------------------------------------------
// Function    :  End_MemFree_PoissonGravity
// Description :  Free memory previously allocated by the function "Init_MemAllocate_PoissonGravity"
//
// Parameter   :  None
//-------------------------------------------------------------------------------------------------------
void End_MemFree_PoissonGravity()
{

   for (int t=0; t<2; t++)
   {
      if ( h_Rho_Array_P    [t] != NULL )    delete [] h_Rho_Array_P    [t];
      if ( h_Pot_Array_P_In [t] != NULL )    delete [] h_Pot_Array_P_In [t];
      if ( h_Pot_Array_P_Out[t] != NULL )    delete [] h_Pot_Array_P_Out[t];
#     ifdef UNSPLIT_GRAVITY
      if ( h_Pot_Array_USG_G[t] != NULL )    delete [] h_Pot_Array_USG_G[t];
      if ( h_Flu_Array_USG_G[t] != NULL )    delete [] h_Flu_Array_USG_G[t];
#     endif
      if ( h_Flu_Array_G    [t] != NULL )    delete [] h_Flu_Array_G    [t];
      if ( h_Corner_Array_G [t] != NULL )    delete [] h_Corner_Array_G [t];
#     ifdef DUAL_ENERGY
      if ( h_DE_Array_G     [t] != NULL )    delete [] h_DE_Array_G     [t];
#     endif
      if ( h_Pot_Array_T    [t] != NULL )    delete [] h_Pot_Array_T    [t];

      h_Rho_Array_P    [t] = NULL;
      h_Pot_Array_P_In [t] = NULL;
      h_Pot_Array_P_Out[t] = NULL;
#     ifdef UNSPLIT_GRAVITY
      h_Pot_Array_USG_G[t] = NULL;
      h_Flu_Array_USG_G[t] = NULL;
#     endif
      h_Flu_Array_G    [t] = NULL;
      h_Corner_Array_G [t] = NULL;
#     ifdef DUAL_ENERGY
      h_DE_Array_G     [t] = NULL;
#     endif
      h_Pot_Array_T    [t] = NULL;
   }

   if ( GreenFuncK != NULL )  delete [] GreenFuncK;

} // FUNCTION : End_MemFree_PoissonGravity



#endif // #if ( !defined GPU  &&  defined GRAVITY )
