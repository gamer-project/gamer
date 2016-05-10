#include "Copyright.h"
#include "GAMER.h"
#include "CUFLU.h"

#ifdef GPU



// *******************************************
// ** CUDA stream objects are declared here **
cudaStream_t *Stream;
// *******************************************


extern real (*d_Flu_Array_F_In )[FLU_NIN ][ FLU_NXT*FLU_NXT*FLU_NXT ];
extern real (*d_Flu_Array_F_Out)[FLU_NOUT][ PS2*PS2*PS2 ];
extern real (*d_Flux_Array)[9][NFLUX][ PS2*PS2 ];
#ifdef UNSPLIT_GRAVITY
extern real (*d_Pot_Array_USG_F)[ USG_NXT_F*USG_NXT_F*USG_NXT_F ];
extern double (*d_Corner_Array_F)[3];
#endif
extern real  *d_MinDtInfo_Fluid_Array;

// global memory arrays in different models
#if ( MODEL == HYDRO )
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




//-------------------------------------------------------------------------------------------------------
// Function    :  CUAPI_MemAllocate_Fluid
// Description :  Allocate GPU and CPU memory for the fluid solver
//
// Parameter   :  Flu_NPatchGroup   : Number of patch groups evaluated simultaneously by GPU 
//                GPU_NStream       : Number of CUDA stream objects
//-------------------------------------------------------------------------------------------------------
void CUAPI_MemAllocate_Fluid( const int Flu_NPatchGroup, const int GPU_NStream )
{

// the size of the global memory arrays in all models
   const long Flu_MemSize_F_In        = sizeof(real  )*Flu_NPatchGroup*FLU_NIN *FLU_NXT*FLU_NXT*FLU_NXT;
   const long Flu_MemSize_F_Out       = sizeof(real  )*Flu_NPatchGroup*FLU_NOUT*PS2*PS2*PS2;
   const long Flux_MemSize            = sizeof(real  )*Flu_NPatchGroup*9*NFLUX*PS2*PS2;
   const long MinDtInfo_Fluid_MemSize = sizeof(real  )*Flu_NPatchGroup;
#  ifdef UNSPLIT_GRAVITY
   const long Pot_MemSize_USG_F       = sizeof(real  )*Flu_NPatchGroup*USG_NXT_F*USG_NXT_F*USG_NXT_F;
   const long Corner_MemSize          = sizeof(double)*Flu_NPatchGroup*3;
#  endif

// the size of the global memory arrays in different models
#  if   ( MODEL == HYDRO )
#  if ( FLU_SCHEME == MHM  ||  FLU_SCHEME == MHM_RP  ||  FLU_SCHEME == CTU )
   const long PriVar_MemSize    = Flu_MemSize_F_In;
   const long FC_Var_MemSize    = sizeof(real)*Flu_NPatchGroup*NCOMP*N_FC_VAR*N_FC_VAR*N_FC_VAR;
   const long FC_Flux_MemSize   = sizeof(real)*Flu_NPatchGroup*NCOMP*N_FC_FLUX*N_FC_FLUX*N_FC_FLUX;

#  if ( LR_SCHEME == PPM )    
   const long Slope_PPM_MemSize = sizeof(real)*Flu_NPatchGroup*NCOMP*N_SLOPE_PPM*N_SLOPE_PPM*N_SLOPE_PPM;
#  endif 

#  endif // #if ( FLU_SCHEME == MHM  ||  FLU_SCHEME == MHM_RP  ||  FLU_SCHEME == CTU )

#  elif ( MODEL == MHD )
#  warning : WAIT MHD !!!

#  elif ( MODEL != ELBDM )
#  warning : DO YOU WANT TO ADD SOMETHING HERE FOR THE NEW MODEL ??
#  endif // MODEL


// output the total memory requirement
   long TotalSize = Flu_MemSize_F_In + Flu_MemSize_F_Out;

   if ( amr->WithFlux )
   TotalSize += Flux_MemSize;

#  ifdef UNSPLIT_GRAVITY
   TotalSize += Pot_MemSize_USG_F;

   if ( OPT__GRAVITY_TYPE == GRAVITY_EXTERNAL  ||  OPT__GRAVITY_TYPE == GRAVITY_BOTH  ||  OPT__EXTERNAL_POT )
   TotalSize += Corner_MemSize;
#  endif
   
   if ( OPT__ADAPTIVE_DT )
   TotalSize += MinDtInfo_Fluid_MemSize;

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

   if ( OPT__ADAPTIVE_DT )
   CUDA_CHECK_ERROR(  cudaMalloc( (void**) &d_MinDtInfo_Fluid_Array, MinDtInfo_Fluid_MemSize )  );


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
      CUDA_CHECK_ERROR(  cudaMallocHost( (void**) &h_Flu_Array_F_In       [t], Flu_MemSize_F_In        )  );
      CUDA_CHECK_ERROR(  cudaMallocHost( (void**) &h_Flu_Array_F_Out      [t], Flu_MemSize_F_Out       )  );

      if ( amr->WithFlux )
      CUDA_CHECK_ERROR(  cudaMallocHost( (void**) &h_Flux_Array           [t], Flux_MemSize            )  );

#     ifdef UNSPLIT_GRAVITY
      CUDA_CHECK_ERROR(  cudaMallocHost( (void**) &h_Pot_Array_USG_F      [t], Pot_MemSize_USG_F       )  );

      if ( OPT__GRAVITY_TYPE == GRAVITY_EXTERNAL  ||  OPT__GRAVITY_TYPE == GRAVITY_BOTH  ||  OPT__EXTERNAL_POT )
      CUDA_CHECK_ERROR(  cudaMallocHost( (void**) &h_Corner_Array_F       [t], Corner_MemSize          )  );
#     endif

      if ( OPT__ADAPTIVE_DT )
      CUDA_CHECK_ERROR(  cudaMallocHost( (void**) &h_MinDtInfo_Fluid_Array[t], MinDtInfo_Fluid_MemSize )  );
   }


// create streams
   Stream = new cudaStream_t [GPU_NStream];
   for (int s=0; s<GPU_NStream; s++)      CUDA_CHECK_ERROR(  cudaStreamCreate( &Stream[s] )  );

} // FUNCTION : CUAPI_MemAllocate_Fluid



#endif // #ifdef GPU
