#include "CUAPI.h"
#include "CUFLU.h"

#ifdef GPU



extern real (*d_Flu_Array_F_In )[FLU_NIN ][ CUBE(FLU_NXT) ];
extern real (*d_Flu_Array_F_Out)[FLU_NOUT][ CUBE(PS2) ];
extern real (*d_Flux_Array)[9][NFLUX_TOTAL][ SQR(PS2) ];
#ifdef UNSPLIT_GRAVITY
extern double (*d_Corner_Array_F)[3];
#endif
#ifdef DUAL_ENERGY
extern char (*d_DE_Array_F_Out)[ CUBE(PS2) ];
#endif
extern real *d_dt_Array_T;
extern real (*d_Flu_Array_T)[NCOMP_FLUID][ CUBE(PS1) ];

// global memory arrays in different models
#if   ( MODEL == HYDRO )
#if ( FLU_SCHEME == MHM  ||  FLU_SCHEME == MHM_RP  ||  FLU_SCHEME == CTU )
extern real (*d_PriVar)      [NCOMP_TOTAL][ CUBE(FLU_NXT)     ];
extern real (*d_Slope_PPM)[3][NCOMP_TOTAL][ CUBE(N_SLOPE_PPM) ];
extern real (*d_FC_Var)   [6][NCOMP_TOTAL][ CUBE(N_FC_VAR)    ];
extern real (*d_FC_Flux)  [3][NCOMP_TOTAL][ CUBE(N_FC_FLUX)   ];
#endif // #if ( FLU_SCHEME == MHM  ||  FLU_SCHEME == MHM_RP  ||  FLU_SCHEME == CTU )

#elif ( MODEL == MHD )
#warning : WAIT MHD !!!

#elif ( MODEL != ELBDM )
#warning : DO YOU WANT TO ADD SOMETHING HERE FOR THE NEW MODEL ??
#endif // MODEL

extern cudaStream_t *Stream;




//-------------------------------------------------------------------------------------------------------
// Function    :  CUAPI_MemFree_Fluid
// Description :  Free the GPU and CPU memory previously allocated by CUAPI_MemAllocate_Fluid()
//
// Parameter   :  GPU_NStream : Number of CUDA streams for the asynchronous memory copy
//-------------------------------------------------------------------------------------------------------
void CUAPI_MemFree_Fluid( const int GPU_NStream )
{

// free the device memory (used by all models)
   if ( d_Flu_Array_F_In  != NULL ) {  CUDA_CHECK_ERROR(  cudaFree( d_Flu_Array_F_In  )  );  d_Flu_Array_F_In  = NULL; }
   if ( d_Flu_Array_F_Out != NULL ) {  CUDA_CHECK_ERROR(  cudaFree( d_Flu_Array_F_Out )  );  d_Flu_Array_F_Out = NULL; }
   if ( d_Flux_Array      != NULL ) {  CUDA_CHECK_ERROR(  cudaFree( d_Flux_Array      )  );  d_Flux_Array      = NULL; }
#  ifdef UNSPLIT_GRAVITY
   if ( d_Corner_Array_F  != NULL ) {  CUDA_CHECK_ERROR(  cudaFree( d_Corner_Array_F  )  );  d_Corner_Array_F  = NULL; }
#  endif
#  ifdef DUAL_ENERGY
   if ( d_DE_Array_F_Out  != NULL ) {  CUDA_CHECK_ERROR(  cudaFree( d_DE_Array_F_Out  )  );  d_DE_Array_F_Out  = NULL; }
#  endif
   if ( d_dt_Array_T      != NULL ) {  CUDA_CHECK_ERROR(  cudaFree( d_dt_Array_T      )  );  d_dt_Array_T      = NULL; }
   if ( d_Flu_Array_T     != NULL ) {  CUDA_CHECK_ERROR(  cudaFree( d_Flu_Array_T     )  );  d_Flu_Array_T     = NULL; }


// free the device memory (used by different models)
#  if   ( MODEL == HYDRO )
#  if ( FLU_SCHEME == MHM  ||  FLU_SCHEME == MHM_RP  ||  FLU_SCHEME == CTU )
   if ( d_PriVar    != NULL ) {  CUDA_CHECK_ERROR(  cudaFree( d_PriVar    )  );  d_PriVar    = NULL; }
   if ( d_Slope_PPM != NULL ) {  CUDA_CHECK_ERROR(  cudaFree( d_Slope_PPM )  );  d_Slope_PPM = NULL; }
   if ( d_FC_Var    != NULL ) {  CUDA_CHECK_ERROR(  cudaFree( d_FC_Var    )  );  d_FC_Var    = NULL; }
   if ( d_FC_Flux   != NULL ) {  CUDA_CHECK_ERROR(  cudaFree( d_FC_Flux   )  );  d_FC_Flux   = NULL; }
#  endif // #if ( FLU_SCHEME == MHM  ||  FLU_SCHEME == MHM_RP  ||  FLU_SCHEME == CTU )

#  elif ( MODEL == MHD )
#  warning : WAIT MHD !!!

#  elif ( MODEL != ELBDM )
#  warning : DO YOU WANT TO ADD SOMETHING HERE FOR THE NEW MODEL ??
#  endif // MODEL


// free the host memory allocated by CUDA
   for (int t=0; t<2; t++)
   {
      if ( h_Flu_Array_F_In [t] != NULL ) {  CUDA_CHECK_ERROR(  cudaFreeHost( h_Flu_Array_F_In [t] )  );  h_Flu_Array_F_In [t] = NULL; }
      if ( h_Flu_Array_F_Out[t] != NULL ) {  CUDA_CHECK_ERROR(  cudaFreeHost( h_Flu_Array_F_Out[t] )  );  h_Flu_Array_F_Out[t] = NULL; }
      if ( h_Flux_Array     [t] != NULL ) {  CUDA_CHECK_ERROR(  cudaFreeHost( h_Flux_Array     [t] )  );  h_Flux_Array     [t] = NULL; }
#     ifdef UNSPLIT_GRAVITY
      if ( h_Corner_Array_F [t] != NULL ) {  CUDA_CHECK_ERROR(  cudaFreeHost( h_Corner_Array_F [t] )  );  h_Corner_Array_F [t] = NULL; }
#     endif
#     ifdef DUAL_ENERGY
      if ( h_DE_Array_F_Out [t] != NULL ) {  CUDA_CHECK_ERROR(  cudaFreeHost( h_DE_Array_F_Out [t] )  );  h_DE_Array_F_Out [t] = NULL; }
#     endif
      if ( h_dt_Array_T     [t] != NULL ) {  CUDA_CHECK_ERROR(  cudaFreeHost( h_dt_Array_T     [t] )  );  h_dt_Array_T     [t] = NULL; }
      if ( h_Flu_Array_T    [t] != NULL ) {  CUDA_CHECK_ERROR(  cudaFreeHost( h_Flu_Array_T    [t] )  );  h_Flu_Array_T    [t] = NULL; }
   } // for (int t=0; t<2; t++)


// destroy streams
   if ( Stream != NULL )
   {
      for (int s=0; s<GPU_NStream; s++)   CUDA_CHECK_ERROR(  cudaStreamDestroy( Stream[s] )  );

      delete [] Stream;
      Stream = NULL;
   }

} // FUNCTION : CUAPI_MemFree_Fluid



#endif // #ifdef GPU
