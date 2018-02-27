#include "CUAPI.h"

#if ( defined GPU  &&  defined GRAVITY )



extern real (*d_Rho_Array_P    )[ RHO_NXT*RHO_NXT*RHO_NXT ];
extern real (*d_Pot_Array_P_In )[ POT_NXT*POT_NXT*POT_NXT ];
extern real (*d_Pot_Array_P_Out)[ GRA_NXT*GRA_NXT*GRA_NXT ];
#ifdef UNSPLIT_GRAVITY
extern real (*d_Pot_Array_USG_F)[ USG_NXT_F*USG_NXT_F*USG_NXT_F ];
extern real (*d_Pot_Array_USG_G)[ USG_NXT_G*USG_NXT_G*USG_NXT_G ];
extern real (*d_Flu_Array_USG_G)[GRA_NIN-1][ PS1*PS1*PS1 ];
#endif
extern real (*d_Flu_Array_G    )[GRA_NIN  ][ PS1*PS1*PS1 ];
extern double (*d_Corner_Array_G)[3];
#ifdef DUAL_ENERGY
extern char (*d_DE_Array_G     )[ PS1*PS1*PS1 ];
#endif
extern real (*d_Pot_Array_T    )[ CUBE(GRA_NXT) ];




//-------------------------------------------------------------------------------------------------------
// Function    :  CUAPI_MemFree_PoissonGravity
// Description :  Free the device and host memory previously allocated by CUAPI_MemAllocate_PoissonGravity()
//
// Parameter   :  None
//-------------------------------------------------------------------------------------------------------
void CUAPI_MemFree_PoissonGravity()
{

// free the device memory
   if ( d_Rho_Array_P     != NULL )    CUDA_CHECK_ERROR(  cudaFree( d_Rho_Array_P     )  );
   if ( d_Pot_Array_P_In  != NULL )    CUDA_CHECK_ERROR(  cudaFree( d_Pot_Array_P_In  )  );
   if ( d_Pot_Array_P_Out != NULL )    CUDA_CHECK_ERROR(  cudaFree( d_Pot_Array_P_Out )  );
#  ifdef UNSPLIT_GRAVITY
   if ( d_Pot_Array_USG_F != NULL )    CUDA_CHECK_ERROR(  cudaFree( d_Pot_Array_USG_F )  );
   if ( d_Pot_Array_USG_G != NULL )    CUDA_CHECK_ERROR(  cudaFree( d_Pot_Array_USG_G )  );
   if ( d_Flu_Array_USG_G != NULL )    CUDA_CHECK_ERROR(  cudaFree( d_Flu_Array_USG_G )  );
#  endif
   if ( d_Flu_Array_G     != NULL )    CUDA_CHECK_ERROR(  cudaFree( d_Flu_Array_G     )  );
   if ( d_Corner_Array_G  != NULL )    CUDA_CHECK_ERROR(  cudaFree( d_Corner_Array_G  )  );
#  ifdef DUAL_ENERGY
   if ( d_DE_Array_G      != NULL )    CUDA_CHECK_ERROR(  cudaFree( d_DE_Array_G      )  );
#  endif
   if ( d_Pot_Array_T     != NULL )    CUDA_CHECK_ERROR(  cudaFree( d_Pot_Array_T     )  );

   d_Rho_Array_P     = NULL;
   d_Pot_Array_P_In  = NULL;
   d_Pot_Array_P_Out = NULL;
#  ifdef UNSPLIT_GRAVITY
   d_Pot_Array_USG_F = NULL;
   d_Pot_Array_USG_G = NULL;
   d_Flu_Array_USG_G = NULL;
#  endif
   d_Flu_Array_G     = NULL;
   d_Corner_Array_G  = NULL;
#  ifdef DUAL_ENERGY
   d_DE_Array_G      = NULL;
#  endif
   d_Pot_Array_T     = NULL;


// free the host memory allocated by CUDA
   for (int t=0; t<2; t++)
   {
      if ( h_Rho_Array_P    [t] != NULL )    CUDA_CHECK_ERROR(  cudaFreeHost( h_Rho_Array_P    [t] )  );
      if ( h_Pot_Array_P_In [t] != NULL )    CUDA_CHECK_ERROR(  cudaFreeHost( h_Pot_Array_P_In [t] )  );
      if ( h_Pot_Array_P_Out[t] != NULL )    CUDA_CHECK_ERROR(  cudaFreeHost( h_Pot_Array_P_Out[t] )  );
#     ifdef UNSPLIT_GRAVITY
      if ( h_Pot_Array_USG_F[t] != NULL )    CUDA_CHECK_ERROR(  cudaFreeHost( h_Pot_Array_USG_F[t] )  );
      if ( h_Pot_Array_USG_G[t] != NULL )    CUDA_CHECK_ERROR(  cudaFreeHost( h_Pot_Array_USG_G[t] )  );
      if ( h_Flu_Array_USG_G[t] != NULL )    CUDA_CHECK_ERROR(  cudaFreeHost( h_Flu_Array_USG_G[t] )  );
#     endif
      if ( h_Flu_Array_G    [t] != NULL )    CUDA_CHECK_ERROR(  cudaFreeHost( h_Flu_Array_G    [t] )  );
      if ( h_Corner_Array_G [t] != NULL )    CUDA_CHECK_ERROR(  cudaFreeHost( h_Corner_Array_G [t] )  );
#     ifdef DUAL_ENERGY
      if ( h_DE_Array_G     [t] != NULL )    CUDA_CHECK_ERROR(  cudaFreeHost( h_DE_Array_G     [t] )  );
#     endif
      if ( h_Pot_Array_T    [t] != NULL )    CUDA_CHECK_ERROR(  cudaFreeHost( h_Pot_Array_T    [t] )  );

      h_Rho_Array_P    [t] = NULL;
      h_Pot_Array_P_In [t] = NULL;
      h_Pot_Array_P_Out[t] = NULL;
#     ifdef UNSPLIT_GRAVITY
      h_Pot_Array_USG_F[t] = NULL;
      h_Pot_Array_USG_G[t] = NULL;
      h_Flu_Array_USG_G[t] = NULL;
#     endif
      h_Flu_Array_G    [t] = NULL;
      h_Corner_Array_G [t] = NULL;
#     ifdef DUAL_ENERGY
      h_DE_Array_G     [t] = NULL;
#     endif
      h_Pot_Array_T    [t] = NULL;
   } // for (int t=0; t<2; t++)

} // FUNCTION : CUAPI_MemFree_PoissonGravity



#endif // #if ( defined GPU  &&  defined GRAVITY )
