#include "Copyright.h"
#include "GAMER.h"
#include "CUFLU.h"

#ifdef GPU



extern real (*d_Flu_Array_F_In )[FLU_NIN ][ FLU_NXT*FLU_NXT*FLU_NXT ];
extern real (*d_Flu_Array_F_Out)[FLU_NOUT][ PS2*PS2*PS2 ];
extern real (*d_Flux_Array)[9][NFLUX][ PS2*PS2 ];
#ifdef UNSPLIT_GRAVITY
extern double (*d_Corner_Array_F)[3];
#endif
extern real  *d_MinDtInfo_Fluid_Array;

// global memory arrays in different models
#if   ( MODEL == HYDRO )
#if ( FLU_SCHEME == MHM  ||  FLU_SCHEME == MHM_RP  ||  FLU_SCHEME == CTU )
extern real (*d_PriVar)     [5][ FLU_NXT*FLU_NXT*FLU_NXT ];
extern real (*d_Slope_PPM_x)[5][ N_SLOPE_PPM*N_SLOPE_PPM*N_SLOPE_PPM ];
extern real (*d_Slope_PPM_y)[5][ N_SLOPE_PPM*N_SLOPE_PPM*N_SLOPE_PPM ];
extern real (*d_Slope_PPM_z)[5][ N_SLOPE_PPM*N_SLOPE_PPM*N_SLOPE_PPM ];
extern real (*d_FC_Var_xL)  [5][ N_FC_VAR*N_FC_VAR*N_FC_VAR ];
extern real (*d_FC_Var_xR)  [5][ N_FC_VAR*N_FC_VAR*N_FC_VAR ];
extern real (*d_FC_Var_yL)  [5][ N_FC_VAR*N_FC_VAR*N_FC_VAR ];
extern real (*d_FC_Var_yR)  [5][ N_FC_VAR*N_FC_VAR*N_FC_VAR ];
extern real (*d_FC_Var_zL)  [5][ N_FC_VAR*N_FC_VAR*N_FC_VAR ];
extern real (*d_FC_Var_zR)  [5][ N_FC_VAR*N_FC_VAR*N_FC_VAR ];
extern real (*d_FC_Flux_x)  [5][ N_FC_FLUX*N_FC_FLUX*N_FC_FLUX ];
extern real (*d_FC_Flux_y)  [5][ N_FC_FLUX*N_FC_FLUX*N_FC_FLUX ];
extern real (*d_FC_Flux_z)  [5][ N_FC_FLUX*N_FC_FLUX*N_FC_FLUX ];
#endif // #if ( FLU_SCHEME == MHM  ||  FLU_SCHEME == MHM_RP  ||  FLU_SCHEME == CTU )

#elif ( MODEL == MHD )
#warning : WAIT MHD !!!

#elif ( MODEL != ELBDM )
#warning : DO YOU WANT TO ADD SOMETHING HERE FOR THE NEW MODEL ??
#endif // MODEL

extern cudaStream_t *Stream;




//-------------------------------------------------------------------------------------------------------
// Function    :  CUAPI_MemFree_Fluid
// Description :  Free GPU and CPU memory previously allocated by the function "CUAPI_MemAllocate_Fluid"
//
// Parameter   :  GPU_NStream : Number of CUDA streams for the asynchronous memory copy
//-------------------------------------------------------------------------------------------------------
void CUAPI_MemFree_Fluid( const int GPU_NStream )
{

// free the device memory (in all models)
   if ( d_Flu_Array_F_In        != NULL )    CUDA_CHECK_ERROR(  cudaFree( d_Flu_Array_F_In        )  );
   if ( d_Flu_Array_F_Out       != NULL )    CUDA_CHECK_ERROR(  cudaFree( d_Flu_Array_F_Out       )  );
   if ( d_Flux_Array            != NULL )    CUDA_CHECK_ERROR(  cudaFree( d_Flux_Array            )  );
   if ( d_MinDtInfo_Fluid_Array != NULL )    CUDA_CHECK_ERROR(  cudaFree( d_MinDtInfo_Fluid_Array )  );
#  ifdef UNSPLIT_GRAVITY
   if ( d_Corner_Array_F        != NULL )    CUDA_CHECK_ERROR(  cudaFree( d_Corner_Array_F        )  );
#  endif

   d_Flu_Array_F_In        = NULL;
   d_Flu_Array_F_Out       = NULL;
   d_Flux_Array            = NULL;
#  ifdef UNSPLIT_GRAVITY
   d_Corner_Array_F        = NULL;
#  endif
   d_MinDtInfo_Fluid_Array = NULL;


// free the device memory (in different models)
#  if   ( MODEL == HYDRO )
#  if ( FLU_SCHEME == MHM  ||  FLU_SCHEME == MHM_RP  ||  FLU_SCHEME == CTU )
   if ( d_PriVar      != NULL )    CUDA_CHECK_ERROR(  cudaFree( d_PriVar      )  );

   if ( d_Slope_PPM_x != NULL )    CUDA_CHECK_ERROR(  cudaFree( d_Slope_PPM_x )  );
   if ( d_Slope_PPM_y != NULL )    CUDA_CHECK_ERROR(  cudaFree( d_Slope_PPM_y )  );
   if ( d_Slope_PPM_z != NULL )    CUDA_CHECK_ERROR(  cudaFree( d_Slope_PPM_z )  );

   if ( d_FC_Var_xL   != NULL )    CUDA_CHECK_ERROR(  cudaFree( d_FC_Var_xL   )  );
   if ( d_FC_Var_xR   != NULL )    CUDA_CHECK_ERROR(  cudaFree( d_FC_Var_xR   )  );
   if ( d_FC_Var_yL   != NULL )    CUDA_CHECK_ERROR(  cudaFree( d_FC_Var_yL   )  );
   if ( d_FC_Var_yR   != NULL )    CUDA_CHECK_ERROR(  cudaFree( d_FC_Var_yR   )  );
   if ( d_FC_Var_zL   != NULL )    CUDA_CHECK_ERROR(  cudaFree( d_FC_Var_zL   )  );
   if ( d_FC_Var_zR   != NULL )    CUDA_CHECK_ERROR(  cudaFree( d_FC_Var_zR   )  );

   if ( d_FC_Flux_x   != NULL )    CUDA_CHECK_ERROR(  cudaFree( d_FC_Flux_x   )  );
   if ( d_FC_Flux_y   != NULL )    CUDA_CHECK_ERROR(  cudaFree( d_FC_Flux_y   )  );
   if ( d_FC_Flux_z   != NULL )    CUDA_CHECK_ERROR(  cudaFree( d_FC_Flux_z   )  );

   d_PriVar      = NULL;

   d_Slope_PPM_x = NULL;
   d_Slope_PPM_y = NULL;
   d_Slope_PPM_z = NULL;

   d_FC_Var_xL   = NULL;
   d_FC_Var_xR   = NULL;
   d_FC_Var_yL   = NULL;
   d_FC_Var_yR   = NULL;
   d_FC_Var_zL   = NULL;
   d_FC_Var_zR   = NULL;

   d_FC_Flux_x   = NULL;
   d_FC_Flux_y   = NULL;
   d_FC_Flux_z   = NULL;
#  endif // #if ( FLU_SCHEME == MHM  ||  FLU_SCHEME == MHM_RP  ||  FLU_SCHEME == CTU )

#  elif ( MODEL == MHD )
#  warning : WAIT MHD !!!

#  elif ( MODEL != ELBDM )
#  warning : DO YOU WANT TO ADD SOMETHING HERE FOR THE NEW MODEL ??
#  endif // MODEL


// free the host memory allocated by CUDA
   for (int t=0; t<2; t++)
   {
      if ( h_Flu_Array_F_In       [t] != NULL ) CUDA_CHECK_ERROR(  cudaFreeHost( h_Flu_Array_F_In       [t] )  );
      if ( h_Flu_Array_F_Out      [t] != NULL ) CUDA_CHECK_ERROR(  cudaFreeHost( h_Flu_Array_F_Out      [t] )  );
      if ( h_Flux_Array           [t] != NULL ) CUDA_CHECK_ERROR(  cudaFreeHost( h_Flux_Array           [t] )  );
#     ifdef UNSPLIT_GRAVITY
      if ( h_Corner_Array_F       [t] != NULL ) CUDA_CHECK_ERROR(  cudaFreeHost( h_Corner_Array_F       [t] )  );
#     endif
      if ( h_MinDtInfo_Fluid_Array[t] != NULL ) CUDA_CHECK_ERROR(  cudaFreeHost( h_MinDtInfo_Fluid_Array[t] )  );

      h_Flu_Array_F_In       [t] = NULL;  
      h_Flu_Array_F_Out      [t] = NULL;
      h_Flux_Array           [t] = NULL;
#     ifdef UNSPLIT_GRAVITY
      h_Corner_Array_F       [t] = NULL;
#     endif
      h_MinDtInfo_Fluid_Array[t] = NULL;
   }


// destroy streams
   if ( Stream != NULL )
   {
      for (int s=0; s<GPU_NStream; s++)   
      {
         CUDA_CHECK_ERROR(  cudaStreamDestroy( Stream[s] )  );
      }

      delete [] Stream;
      Stream = NULL;
   }

} // FUNCTION : CUAPI_MemFree_Fluid



#endif // #ifdef GPU
