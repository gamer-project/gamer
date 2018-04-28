#include "CUAPI.h"
#include "CUFLU.h"

#ifdef GPU



// *******************************************
// ** CUDA stream objects are declared here **
cudaStream_t *Stream;
// *******************************************


extern real (*d_Flu_Array_F_In )[FLU_NIN ][ FLU_NXT*FLU_NXT*FLU_NXT ];
extern real (*d_Flu_Array_F_Out)[FLU_NOUT][ PS2*PS2*PS2 ];
extern real (*d_Flux_Array)[9][NFLUX_TOTAL][ PS2*PS2 ];
#ifdef UNSPLIT_GRAVITY
extern real (*d_Pot_Array_USG_F)[ USG_NXT_F*USG_NXT_F*USG_NXT_F ];
extern double (*d_Corner_Array_F)[3];
#endif
#ifdef DUAL_ENERGY
extern char (*d_DE_Array_F_Out)[ PS2*PS2*PS2 ];
#endif
extern real *d_dt_Array_T;
extern real (*d_Flu_Array_T)[NCOMP_FLUID][ CUBE(PS1) ];

// global memory arrays in different models
#if ( MODEL == HYDRO )
#if ( FLU_SCHEME == MHM  ||  FLU_SCHEME == MHM_RP  ||  FLU_SCHEME == CTU )
extern real (*d_PriVar)     [NCOMP_TOTAL][ FLU_NXT*FLU_NXT*FLU_NXT ];
extern real (*d_Slope_PPM_x)[NCOMP_TOTAL][ N_SLOPE_PPM*N_SLOPE_PPM*N_SLOPE_PPM ];
extern real (*d_Slope_PPM_y)[NCOMP_TOTAL][ N_SLOPE_PPM*N_SLOPE_PPM*N_SLOPE_PPM ];
extern real (*d_Slope_PPM_z)[NCOMP_TOTAL][ N_SLOPE_PPM*N_SLOPE_PPM*N_SLOPE_PPM ];
extern real (*d_FC_Var_xL)  [NCOMP_TOTAL][ N_FC_VAR*N_FC_VAR*N_FC_VAR ];
extern real (*d_FC_Var_xR)  [NCOMP_TOTAL][ N_FC_VAR*N_FC_VAR*N_FC_VAR ];
extern real (*d_FC_Var_yL)  [NCOMP_TOTAL][ N_FC_VAR*N_FC_VAR*N_FC_VAR ];
extern real (*d_FC_Var_yR)  [NCOMP_TOTAL][ N_FC_VAR*N_FC_VAR*N_FC_VAR ];
extern real (*d_FC_Var_zL)  [NCOMP_TOTAL][ N_FC_VAR*N_FC_VAR*N_FC_VAR ];
extern real (*d_FC_Var_zR)  [NCOMP_TOTAL][ N_FC_VAR*N_FC_VAR*N_FC_VAR ];
extern real (*d_FC_Flux_x)  [NCOMP_TOTAL][ N_FC_FLUX*N_FC_FLUX*N_FC_FLUX ];
extern real (*d_FC_Flux_y)  [NCOMP_TOTAL][ N_FC_FLUX*N_FC_FLUX*N_FC_FLUX ];
extern real (*d_FC_Flux_z)  [NCOMP_TOTAL][ N_FC_FLUX*N_FC_FLUX*N_FC_FLUX ];
#endif // #if ( FLU_SCHEME == MHM  ||  FLU_SCHEME == MHM_RP  ||  FLU_SCHEME == CTU )

#elif ( MODEL == MHD )
#warning : WAIT MHD !!!

#elif ( MODEL != ELBDM )
#warning : DO YOU WANT TO ADD SOMETHING HERE FOR THE NEW MODEL ??
#endif // MODEL




//-------------------------------------------------------------------------------------------------------
// Function    :  CUAPI_MemAllocate_Fluid
// Description :  Allocate GPU and CPU memory for the fluid solver
//
// Parameter   :  Flu_NPG     : Number of patch groups evaluated simultaneously by GPU for the fluid solver
//                Pot_NPG     : Number of patch groups evaluated simultaneously by GPU for the gravity solver
//                              --> Here it is used only for the dt solver
//                GPU_NStream : Number of CUDA stream objects
//-------------------------------------------------------------------------------------------------------
void CUAPI_MemAllocate_Fluid( const int Flu_NPG, const int Pot_NPG, const int GPU_NStream )
{

// size of the global memory arrays in all models
   const int  Flu_NP            = 8*Flu_NPG;
#  ifdef GRAVITY
   const int  Pot_NP            = 8*Pot_NPG;
#  endif
   const long Flu_MemSize_F_In  = sizeof(real  )*Flu_NPG*FLU_NIN *FLU_NXT*FLU_NXT*FLU_NXT;
   const long Flu_MemSize_F_Out = sizeof(real  )*Flu_NPG*FLU_NOUT*PS2*PS2*PS2;
   const long Flux_MemSize      = sizeof(real  )*Flu_NPG*9*NFLUX_TOTAL*PS2*PS2;
#  ifdef UNSPLIT_GRAVITY
   const long Pot_MemSize_USG_F = sizeof(real  )*Flu_NPG*USG_NXT_F*USG_NXT_F*USG_NXT_F;
   const long Corner_MemSize    = sizeof(double)*Flu_NPG*3;
#  endif
#  ifdef DUAL_ENERGY
   const long DE_MemSize_F_Out  = sizeof(char  )*Flu_NPG*PS2*PS2*PS2;
#  endif
#  ifdef GRAVITY
   const long dt_MemSize_T      = sizeof(real  )*MAX( Flu_NP, Pot_NP ); // dt_Array_T is used for both DT_FLU_SOLVER and DT_GRA_SOLVER
#  else
   const long dt_MemSize_T      = sizeof(real  )*Flu_NP;
#  endif
   const long Flu_MemSize_T     = sizeof(real  )*Flu_NP*NCOMP_FLUID*CUBE(PS1);

// the size of the global memory arrays in different models
#  if   ( MODEL == HYDRO )
#  if ( FLU_SCHEME == MHM  ||  FLU_SCHEME == MHM_RP  ||  FLU_SCHEME == CTU )
   const long PriVar_MemSize    = Flu_MemSize_F_In;
   const long FC_Var_MemSize    = sizeof(real)*Flu_NPG*NCOMP_TOTAL*N_FC_VAR*N_FC_VAR*N_FC_VAR;
   const long FC_Flux_MemSize   = sizeof(real)*Flu_NPG*NCOMP_TOTAL*N_FC_FLUX*N_FC_FLUX*N_FC_FLUX;

#  if ( LR_SCHEME == PPM )
   const long Slope_PPM_MemSize = sizeof(real)*Flu_NPG*NCOMP_TOTAL*N_SLOPE_PPM*N_SLOPE_PPM*N_SLOPE_PPM;
#  endif

#  endif // #if ( FLU_SCHEME == MHM  ||  FLU_SCHEME == MHM_RP  ||  FLU_SCHEME == CTU )

#  elif ( MODEL == MHD )
#  warning : WAIT MHD !!!

#  elif ( MODEL != ELBDM )
#  warning : DO YOU WANT TO ADD SOMETHING HERE FOR THE NEW MODEL ??
#  endif // MODEL


// output the total memory requirement
   long TotalSize = Flu_MemSize_F_In + Flu_MemSize_F_Out + dt_MemSize_T + Flu_MemSize_T;

   if ( amr->WithFlux )
   TotalSize += Flux_MemSize;

#  ifdef UNSPLIT_GRAVITY
   TotalSize += Pot_MemSize_USG_F;

   if ( OPT__GRAVITY_TYPE == GRAVITY_EXTERNAL  ||  OPT__GRAVITY_TYPE == GRAVITY_BOTH  ||  OPT__EXTERNAL_POT )
   TotalSize += Corner_MemSize;
#  endif

#  ifdef DUAL_ENERGY
   TotalSize += DE_MemSize_F_Out;
#  endif

#  if   ( MODEL == HYDRO )

#  if ( FLU_SCHEME == MHM  ||  FLU_SCHEME == MHM_RP  ||  FLU_SCHEME == CTU )
   TotalSize += PriVar_MemSize + 6*FC_Var_MemSize + 3*FC_Flux_MemSize;

#  if ( LR_SCHEME == PPM )
   TotalSize += 3*Slope_PPM_MemSize;
#  endif // PPM
#  endif // MHM/MHM_RP/CTU

#  elif ( MODEL == MHD )
#  warning : WAIT MHD !!!

#  elif ( MODEL != ELBDM )
#  warning : DO YOU WANT TO ADD SOMETHING HERE FOR THE NEW MODEL ??
#  endif // MODEL

   if ( MPI_Rank == 0 )
      Aux_Message( stdout, "NOTE : total memory requirement in GPU fluid solver = %ld MB\n", TotalSize/(1<<20) );


// allocate the device memory (in all models)
   CUDA_CHECK_ERROR(  cudaMalloc( (void**) &d_Flu_Array_F_In,        Flu_MemSize_F_In        )  );
   CUDA_CHECK_ERROR(  cudaMalloc( (void**) &d_Flu_Array_F_Out,       Flu_MemSize_F_Out       )  );

   if ( amr->WithFlux )
   CUDA_CHECK_ERROR(  cudaMalloc( (void**) &d_Flux_Array,            Flux_MemSize            )  );

#  ifdef UNSPLIT_GRAVITY
   CUDA_CHECK_ERROR(  cudaMalloc( (void**) &d_Pot_Array_USG_F,       Pot_MemSize_USG_F       )  );

   if ( OPT__GRAVITY_TYPE == GRAVITY_EXTERNAL  ||  OPT__GRAVITY_TYPE == GRAVITY_BOTH  ||  OPT__EXTERNAL_POT )
   CUDA_CHECK_ERROR(  cudaMalloc( (void**) &d_Corner_Array_F,        Corner_MemSize          )  );
#  endif

#  ifdef DUAL_ENERGY
   CUDA_CHECK_ERROR(  cudaMalloc( (void**) &d_DE_Array_F_Out,        DE_MemSize_F_Out        )  );
#  endif

   CUDA_CHECK_ERROR(  cudaMalloc( (void**) &d_dt_Array_T,            dt_MemSize_T            )  );
   CUDA_CHECK_ERROR(  cudaMalloc( (void**) &d_Flu_Array_T,           Flu_MemSize_T           )  );


// allocate the device memory (in different models)
#  if   ( MODEL == HYDRO )
#  if ( FLU_SCHEME == MHM  ||  FLU_SCHEME == MHM_RP  ||  FLU_SCHEME == CTU )
   CUDA_CHECK_ERROR(  cudaMalloc( (void**) &d_PriVar, PriVar_MemSize )  );

#  if ( LR_SCHEME == PPM )
   CUDA_CHECK_ERROR(  cudaMalloc( (void**) &d_Slope_PPM_x, Slope_PPM_MemSize )  );
   CUDA_CHECK_ERROR(  cudaMalloc( (void**) &d_Slope_PPM_y, Slope_PPM_MemSize )  );
   CUDA_CHECK_ERROR(  cudaMalloc( (void**) &d_Slope_PPM_z, Slope_PPM_MemSize )  );
#  endif

   CUDA_CHECK_ERROR(  cudaMalloc( (void**) &d_FC_Var_xL, FC_Var_MemSize )  );
   CUDA_CHECK_ERROR(  cudaMalloc( (void**) &d_FC_Var_xR, FC_Var_MemSize )  );
   CUDA_CHECK_ERROR(  cudaMalloc( (void**) &d_FC_Var_yL, FC_Var_MemSize )  );
   CUDA_CHECK_ERROR(  cudaMalloc( (void**) &d_FC_Var_yR, FC_Var_MemSize )  );
   CUDA_CHECK_ERROR(  cudaMalloc( (void**) &d_FC_Var_zL, FC_Var_MemSize )  );
   CUDA_CHECK_ERROR(  cudaMalloc( (void**) &d_FC_Var_zR, FC_Var_MemSize )  );

   CUDA_CHECK_ERROR(  cudaMalloc( (void**) &d_FC_Flux_x, FC_Flux_MemSize )  );
   CUDA_CHECK_ERROR(  cudaMalloc( (void**) &d_FC_Flux_y, FC_Flux_MemSize )  );
   CUDA_CHECK_ERROR(  cudaMalloc( (void**) &d_FC_Flux_z, FC_Flux_MemSize )  );
#  endif // #if ( FLU_SCHEME == MHM  ||  FLU_SCHEME == MHM_RP  ||  FLU_SCHEME == CTU )

#  elif ( MODEL == MHD )
#  warning : WAIT MHD !!!

#  elif ( MODEL != ELBDM )
#  warning : DO YOU WANT TO ADD SOMETHING HERE FOR THE NEW MODEL ??
#  endif // MODEL


// allocate the host memory by CUDA
   for (int t=0; t<2; t++)
   {
      CUDA_CHECK_ERROR(  cudaMallocHost( (void**) &h_Flu_Array_F_In [t], Flu_MemSize_F_In        )  );
      CUDA_CHECK_ERROR(  cudaMallocHost( (void**) &h_Flu_Array_F_Out[t], Flu_MemSize_F_Out       )  );

      if ( amr->WithFlux )
      CUDA_CHECK_ERROR(  cudaMallocHost( (void**) &h_Flux_Array     [t], Flux_MemSize            )  );

#     ifdef UNSPLIT_GRAVITY
      CUDA_CHECK_ERROR(  cudaMallocHost( (void**) &h_Pot_Array_USG_F[t], Pot_MemSize_USG_F       )  );

      if ( OPT__GRAVITY_TYPE == GRAVITY_EXTERNAL  ||  OPT__GRAVITY_TYPE == GRAVITY_BOTH  ||  OPT__EXTERNAL_POT )
      CUDA_CHECK_ERROR(  cudaMallocHost( (void**) &h_Corner_Array_F [t], Corner_MemSize          )  );
#     endif

#     ifdef DUAL_ENERGY
      CUDA_CHECK_ERROR(  cudaMallocHost( (void**) &h_DE_Array_F_Out [t], DE_MemSize_F_Out        )  );
#     endif

      CUDA_CHECK_ERROR(  cudaMallocHost( (void**) &h_dt_Array_T     [t], dt_MemSize_T            )  );
      CUDA_CHECK_ERROR(  cudaMallocHost( (void**) &h_Flu_Array_T    [t], Flu_MemSize_T           )  );
   } // for (int t=0; t<2; t++)


// create streams
   Stream = new cudaStream_t [GPU_NStream];
   for (int s=0; s<GPU_NStream; s++)      CUDA_CHECK_ERROR(  cudaStreamCreate( &Stream[s] )  );

} // FUNCTION : CUAPI_MemAllocate_Fluid



#endif // #ifdef GPU
