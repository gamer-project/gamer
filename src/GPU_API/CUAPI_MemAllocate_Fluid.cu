#include "CUAPI.h"
#include "CUFLU.h"

#ifdef GPU



// *******************************************
// ** CUDA stream objects are declared here **
cudaStream_t *Stream;
// *******************************************


extern real (*d_Flu_Array_F_In )[FLU_NIN ][ CUBE(FLU_NXT) ];
extern real (*d_Flu_Array_F_Out)[FLU_NOUT][ CUBE(PS2) ];
extern real (*d_Flux_Array)[9][NFLUX_TOTAL][ SQR(PS2) ];
#ifdef UNSPLIT_GRAVITY
extern real (*d_Pot_Array_USG_F)[ CUBE(USG_NXT_F) ];
extern double (*d_Corner_Array_F)[3];
#endif
#ifdef DUAL_ENERGY
extern char (*d_DE_Array_F_Out)[ CUBE(PS2) ];
#endif
#ifdef MHD
extern real (*d_Mag_Array_F_In )[NCOMP_MAG][ FLU_NXT_P1*SQR(FLU_NXT) ];
extern real (*d_Mag_Array_F_Out)[NCOMP_MAG][ PS2P1*SQR(PS2)          ];
extern real (*d_Ele_Array      )[9][NCOMP_ELE][ PS2P1*PS2 ];
extern real (*d_Mag_Array_T)[NCOMP_MAG][ PS1P1*SQR(PS1) ];
extern real (*d_Mag_Array_S_In)[NCOMP_MAG][ SRC_NXT_P1*SQR(SRC_NXT) ];
#endif
extern real *d_dt_Array_T;
extern real (*d_Flu_Array_T)[FLU_NIN_T][ CUBE(PS1) ];
extern real (*d_Flu_Array_S_In )[FLU_NIN_S ][ CUBE(SRC_NXT) ];
extern real (*d_Flu_Array_S_Out)[FLU_NOUT_S][ CUBE(PS1)     ];
extern double (*d_Corner_Array_S)[3];
#if ( FLU_SCHEME == MHM  ||  FLU_SCHEME == MHM_RP  ||  FLU_SCHEME == CTU )
extern real (*d_PriVar)      [NCOMP_LR            ][ CUBE(FLU_NXT)     ];
extern real (*d_Slope_PPM)[3][NCOMP_LR            ][ CUBE(N_SLOPE_PPM) ];
extern real (*d_FC_Var)   [6][NCOMP_TOTAL_PLUS_MAG][ CUBE(N_FC_VAR)    ];
extern real (*d_FC_Flux)  [3][NCOMP_TOTAL_PLUS_MAG][ CUBE(N_FC_FLUX)   ];
#ifdef MHD
extern real (*d_FC_Mag_Half)[NCOMP_MAG][ FLU_NXT_P1*SQR(FLU_NXT) ];
extern real (*d_EC_Ele     )[NCOMP_MAG][ CUBE(N_EC_ELE)          ];
#endif
#endif // FLU_SCHEME

#if ( MODEL == ELBDM )
extern bool (*d_IsCompletelyRefined);
#endif

#if ( ELBDM_SCHEME == ELBDM_HYBRID )
extern bool (*d_HasWaveCounterpart)[ CUBE(HYB_NXT) ];
#endif

#if ( GRAMFE_SCHEME == GRAMFE_MATMUL )
extern gramfe_matmul_float (*d_Flu_TimeEvo)[ 2*FLU_NXT ];
#endif

#if ( MODEL != HYDRO  &&  MODEL != ELBDM )
#  warning : DO YOU WANT TO ADD SOMETHING HERE FOR THE NEW MODEL ??
#endif




//-------------------------------------------------------------------------------------------------------
// Function    :  CUAPI_MemAllocate_Fluid
// Description :  Allocate GPU and CPU memory for the fluid solver
//
// Parameter   :  Flu_NPG     : Number of patch groups evaluated simultaneously by GPU for the fluid solver
//                Pot_NPG     : Number of patch groups evaluated simultaneously by GPU for the gravity solver
//                              --> Here it is used only for the dt solver
//                Src_NPG     : Number of patch groups evaluated simultaneously by GPU for the source-term solver
//                GPU_NStream : Number of CUDA stream objects
//
// Return      :  GAMER_SUCCESS / GAMER_FAILED
//-------------------------------------------------------------------------------------------------------
int CUAPI_MemAllocate_Fluid( const int Flu_NPG, const int Pot_NPG, const int Src_NPG, const int GPU_NStream )
{

// size of the global memory arrays in all models
   const int  Flu_NP              = 8*Flu_NPG;
#  ifdef GRAVITY
   const int  Pot_NP              = 8*Pot_NPG;
#  endif
   const int  Src_NP              = 8*Src_NPG;
   const long Flu_MemSize_F_In    = sizeof(real  )*Flu_NPG*FLU_NIN *CUBE(FLU_NXT);
   const long Flu_MemSize_F_Out   = sizeof(real  )*Flu_NPG*FLU_NOUT*CUBE(PS2);
   const long Flux_MemSize        = sizeof(real  )*Flu_NPG*9*NFLUX_TOTAL*SQR(PS2);
#  ifdef UNSPLIT_GRAVITY
   const long Pot_MemSize_USG_F   = sizeof(real  )*Flu_NPG*CUBE(USG_NXT_F);
   const long Corner_MemSize_F    = sizeof(double)*Flu_NPG*3;
#  endif
#  ifdef DUAL_ENERGY
   const long DE_MemSize_F_Out    = sizeof(char  )*Flu_NPG*CUBE(PS2);
#  endif
#  ifdef MHD
   const long Mag_MemSize_F_In    = sizeof(real  )*Flu_NPG*NCOMP_MAG*FLU_NXT_P1*SQR(FLU_NXT);
   const long Mag_MemSize_F_Out   = sizeof(real  )*Flu_NPG*NCOMP_MAG*PS2P1*SQR(PS2);
   const long Ele_MemSize         = sizeof(real  )*Flu_NPG*9*NCOMP_ELE*PS2P1*PS2;
   const long Mag_MemSize_T       = sizeof(real  )*Flu_NP*NCOMP_MAG*PS1P1*SQR(PS1);
   const long Mag_MemSize_S_In    = sizeof(real  )*Src_NP*NCOMP_MAG*SRC_NXT_P1*SQR(SRC_NXT);
#  endif
#  ifdef GRAVITY
   const long dt_MemSize_T        = sizeof(real  )*MAX( Flu_NP, Pot_NP ); // dt_Array_T is used for both DT_FLU_SOLVER and DT_GRA_SOLVER
#  else
   const long dt_MemSize_T        = sizeof(real  )*Flu_NP;
#  endif
   const long Flu_MemSize_T       = sizeof(real  )*Flu_NP*FLU_NIN_T *CUBE(PS1);
   const long Flu_MemSize_S_In    = sizeof(real  )*Src_NP*FLU_NIN_S *CUBE(SRC_NXT);
   const long Flu_MemSize_S_Out   = sizeof(real  )*Src_NP*FLU_NOUT_S*CUBE(PS1);
   const long Corner_MemSize_S    = sizeof(double)*Src_NP*3;

// the size of the global memory arrays in different models
#  if ( FLU_SCHEME == MHM  ||  FLU_SCHEME == MHM_RP  ||  FLU_SCHEME == CTU )
   const long PriVar_MemSize      = sizeof(real  )*Flu_NPG  *NCOMP_LR            *CUBE(FLU_NXT);
   const long FC_Var_MemSize      = sizeof(real  )*Flu_NPG*6*NCOMP_TOTAL_PLUS_MAG*CUBE(N_FC_VAR);
   const long FC_Flux_MemSize     = sizeof(real  )*Flu_NPG*3*NCOMP_TOTAL_PLUS_MAG*CUBE(N_FC_FLUX);
#  if ( LR_SCHEME == PPM )
   const long Slope_PPM_MemSize   = sizeof(real  )*Flu_NPG*3*NCOMP_LR            *CUBE(N_SLOPE_PPM);
#  endif
#  ifdef MHD
   const long FC_Mag_Half_MemSize = sizeof(real  )*Flu_NPG  *NCOMP_MAG*FLU_NXT_P1*SQR(FLU_NXT);
   const long EC_Ele_MemSize      = sizeof(real  )*Flu_NPG  *NCOMP_MAG*CUBE(N_EC_ELE);
#  endif
#  endif // FLU_SCHEME

#  if ( MODEL == ELBDM )
   const long Flu_MemSize_IsCompletelyRefined = sizeof(bool )*Flu_NPG;
#  endif

#  if ( ELBDM_SCHEME == ELBDM_HYBRID )
   const long Flu_MemSize_HasWaveCounterpart  = sizeof(bool )*Flu_NPG*CUBE(HYB_NXT);
#  endif

#  if ( GRAMFE_SCHEME == GRAMFE_MATMUL )
   const long GramFE_TimeEvo_MemSize          = sizeof(gramfe_matmul_float)*2*PS2*FLU_NXT;
#  endif

#  if ( MODEL != HYDRO  &&  MODEL != ELBDM )
#     warning : DO YOU WANT TO ADD SOMETHING HERE FOR THE NEW MODEL ??
#  endif


// output the total memory requirement
   long TotalSize = Flu_MemSize_F_In + Flu_MemSize_F_Out + dt_MemSize_T + Flu_MemSize_T;

   if ( amr->WithFlux )
   TotalSize += Flux_MemSize;

#  ifdef UNSPLIT_GRAVITY
   TotalSize += Pot_MemSize_USG_F;

   if ( OPT__EXT_ACC )
   TotalSize += Corner_MemSize_F;
#  endif

#  ifdef DUAL_ENERGY
   TotalSize += DE_MemSize_F_Out;
#  endif

#  ifdef MHD
   TotalSize += Mag_MemSize_F_In + Mag_MemSize_F_Out + Mag_MemSize_T;

   if ( amr->WithElectric )
   TotalSize += Ele_MemSize;
#  endif

#  if ( FLU_SCHEME == MHM  ||  FLU_SCHEME == MHM_RP  ||  FLU_SCHEME == CTU )
   TotalSize += PriVar_MemSize + FC_Var_MemSize + FC_Flux_MemSize;

#  if ( LR_SCHEME == PPM )
   TotalSize += Slope_PPM_MemSize;
#  endif

#  ifdef MHD
   TotalSize += FC_Mag_Half_MemSize + EC_Ele_MemSize;
#  endif
#  endif // MHM/MHM_RP/CTU

#  if ( MODEL == ELBDM )
   TotalSize += Flu_MemSize_IsCompletelyRefined;
#  endif

#  if ( ELBDM_SCHEME == ELBDM_HYBRID )
   TotalSize += Flu_MemSize_HasWaveCounterpart;
#  endif

#  if ( GRAMFE_SCHEME == GRAMFE_MATMUL )
   TotalSize += GramFE_TimeEvo_MemSize;
#  endif

#  if ( MODEL != HYDRO  &&  MODEL != ELBDM )
#     warning : DO YOU WANT TO ADD SOMETHING HERE FOR THE NEW MODEL ??
#  endif

   if ( SrcTerms.Any )
   {
      TotalSize += Flu_MemSize_S_In + Flu_MemSize_S_Out;
#     ifdef MHD
      TotalSize += Mag_MemSize_S_In;
#     endif
      TotalSize += Corner_MemSize_S;
   }

   if ( MPI_Rank == 0 )
      Aux_Message( stdout, "NOTE : total memory requirement in GPU fluid solver = %ld MB\n", TotalSize/(1<<20) );


// allocate the device memory
   CUDA_CHECK_MALLOC(  cudaMalloc( (void**) &d_Flu_Array_F_In,       Flu_MemSize_F_In     )  );
   CUDA_CHECK_MALLOC(  cudaMalloc( (void**) &d_Flu_Array_F_Out,      Flu_MemSize_F_Out    )  );

   if ( amr->WithFlux )
   CUDA_CHECK_MALLOC(  cudaMalloc( (void**) &d_Flux_Array,           Flux_MemSize         )  );

#  ifdef UNSPLIT_GRAVITY
   CUDA_CHECK_MALLOC(  cudaMalloc( (void**) &d_Pot_Array_USG_F,      Pot_MemSize_USG_F    )  );

   if ( OPT__EXT_ACC )
   CUDA_CHECK_MALLOC(  cudaMalloc( (void**) &d_Corner_Array_F,       Corner_MemSize_F     )  );
#  endif

#  ifdef DUAL_ENERGY
   CUDA_CHECK_MALLOC(  cudaMalloc( (void**) &d_DE_Array_F_Out,       DE_MemSize_F_Out     )  );
#  endif

#  ifdef MHD
   CUDA_CHECK_MALLOC(  cudaMalloc( (void**) &d_Mag_Array_F_In,       Mag_MemSize_F_In     )  );
   CUDA_CHECK_MALLOC(  cudaMalloc( (void**) &d_Mag_Array_F_Out,      Mag_MemSize_F_Out    )  );

   if ( amr->WithElectric )
   CUDA_CHECK_MALLOC(  cudaMalloc( (void**) &d_Ele_Array,            Ele_MemSize          )  );

   CUDA_CHECK_MALLOC(  cudaMalloc( (void**) &d_Mag_Array_T,          Mag_MemSize_T        )  );
#  endif

   CUDA_CHECK_MALLOC(  cudaMalloc( (void**) &d_dt_Array_T,           dt_MemSize_T         )  );
   CUDA_CHECK_MALLOC(  cudaMalloc( (void**) &d_Flu_Array_T,          Flu_MemSize_T        )  );

#  if ( FLU_SCHEME == MHM  ||  FLU_SCHEME == MHM_RP  ||  FLU_SCHEME == CTU )
   CUDA_CHECK_MALLOC(  cudaMalloc( (void**) &d_FC_Var,               FC_Var_MemSize       )  );

   CUDA_CHECK_MALLOC(  cudaMalloc( (void**) &d_FC_Flux,              FC_Flux_MemSize      )  );

   CUDA_CHECK_MALLOC(  cudaMalloc( (void**) &d_PriVar,               PriVar_MemSize       )  );

#  if ( LR_SCHEME == PPM )
   CUDA_CHECK_MALLOC(  cudaMalloc( (void**) &d_Slope_PPM,            Slope_PPM_MemSize    )  );
#  endif
#  ifdef MHD
   CUDA_CHECK_MALLOC(  cudaMalloc( (void**) &d_FC_Mag_Half,          FC_Mag_Half_MemSize  )  );
   CUDA_CHECK_MALLOC(  cudaMalloc( (void**) &d_EC_Ele,               EC_Ele_MemSize       )  );
#  endif
#  endif // #if ( FLU_SCHEME == MHM  ||  FLU_SCHEME == MHM_RP  ||  FLU_SCHEME == CTU )

   if ( SrcTerms.Any ) {
   CUDA_CHECK_MALLOC(  cudaMalloc( (void**) &d_Flu_Array_S_In,       Flu_MemSize_S_In     )  );
   CUDA_CHECK_MALLOC(  cudaMalloc( (void**) &d_Flu_Array_S_Out,      Flu_MemSize_S_Out    )  );
#  ifdef MHD
   CUDA_CHECK_MALLOC(  cudaMalloc( (void**) &d_Mag_Array_S_In,       Mag_MemSize_S_In     )  );
#  endif
   CUDA_CHECK_MALLOC(  cudaMalloc( (void**) &d_Corner_Array_S,       Corner_MemSize_S     )  );
   }


#  if ( MODEL == ELBDM )
   CUDA_CHECK_MALLOC(  cudaMalloc( (void**) &d_IsCompletelyRefined,  Flu_MemSize_IsCompletelyRefined )  );
#  endif

#  if ( ELBDM_SCHEME == ELBDM_HYBRID )
   CUDA_CHECK_MALLOC(  cudaMalloc( (void**) &d_HasWaveCounterpart,   Flu_MemSize_HasWaveCounterpart  )  );
#  endif

#  if ( GRAMFE_SCHEME == GRAMFE_MATMUL )
   CUDA_CHECK_MALLOC(  cudaMalloc( (void**) &d_Flu_TimeEvo,          GramFE_TimeEvo_MemSize          )  );
#  endif


#  if ( MODEL != HYDRO  &&  MODEL != ELBDM )
#     warning : DO YOU WANT TO ADD SOMETHING HERE FOR THE NEW MODEL ??
#  endif


// allocate the host memory by CUDA
   for (int t=0; t<2; t++)
   {
      CUDA_CHECK_MALLOC(  cudaMallocHost( (void**) &h_Flu_Array_F_In     [t],  Flu_MemSize_F_In     )  );
      CUDA_CHECK_MALLOC(  cudaMallocHost( (void**) &h_Flu_Array_F_Out    [t],  Flu_MemSize_F_Out    )  );

      if ( amr->WithFlux )
      CUDA_CHECK_MALLOC(  cudaMallocHost( (void**) &h_Flux_Array         [t],  Flux_MemSize         )  );

#     ifdef UNSPLIT_GRAVITY
      CUDA_CHECK_MALLOC(  cudaMallocHost( (void**) &h_Pot_Array_USG_F    [t],  Pot_MemSize_USG_F    )  );

      if ( OPT__EXT_ACC )
      CUDA_CHECK_MALLOC(  cudaMallocHost( (void**) &h_Corner_Array_F     [t],  Corner_MemSize_F     )  );
#     endif

#     ifdef DUAL_ENERGY
      CUDA_CHECK_MALLOC(  cudaMallocHost( (void**) &h_DE_Array_F_Out     [t],  DE_MemSize_F_Out     )  );
#     endif

#     ifdef MHD
      CUDA_CHECK_MALLOC(  cudaMallocHost( (void**) &h_Mag_Array_F_In     [t],  Mag_MemSize_F_In     )  );
      CUDA_CHECK_MALLOC(  cudaMallocHost( (void**) &h_Mag_Array_F_Out    [t],  Mag_MemSize_F_Out    )  );

      if ( amr->WithElectric )
      CUDA_CHECK_MALLOC(  cudaMallocHost( (void**) &h_Ele_Array          [t],  Ele_MemSize          )  );

      CUDA_CHECK_MALLOC(  cudaMallocHost( (void**) &h_Mag_Array_T        [t],  Mag_MemSize_T        )  );
#     endif

      CUDA_CHECK_MALLOC(  cudaMallocHost( (void**) &h_dt_Array_T         [t],  dt_MemSize_T         )  );
      CUDA_CHECK_MALLOC(  cudaMallocHost( (void**) &h_Flu_Array_T        [t],  Flu_MemSize_T        )  );

      if ( SrcTerms.Any ) {
      CUDA_CHECK_MALLOC(  cudaMallocHost( (void**) &h_Flu_Array_S_In     [t],  Flu_MemSize_S_In     )  );
      CUDA_CHECK_MALLOC(  cudaMallocHost( (void**) &h_Flu_Array_S_Out    [t],  Flu_MemSize_S_Out    )  );
#     ifdef MHD
      CUDA_CHECK_MALLOC(  cudaMallocHost( (void**) &h_Mag_Array_S_In     [t],  Mag_MemSize_S_In     )  );
#     endif
      CUDA_CHECK_MALLOC(  cudaMallocHost( (void**) &h_Corner_Array_S     [t],  Corner_MemSize_S     )  );
      }

#     if ( MODEL == ELBDM )
      CUDA_CHECK_MALLOC(  cudaMallocHost( (void**) &h_IsCompletelyRefined[t],  Flu_MemSize_IsCompletelyRefined  )  );
#     endif

#     if ( ELBDM_SCHEME == ELBDM_HYBRID )
      CUDA_CHECK_MALLOC(  cudaMallocHost( (void**) &h_HasWaveCounterpart [t],  Flu_MemSize_HasWaveCounterpart   )  );
#     endif
   } // for (int t=0; t<2; t++)

#  if ( GRAMFE_SCHEME == GRAMFE_MATMUL )
   CUDA_CHECK_MALLOC(  cudaMallocHost( (void**) &h_GramFE_TimeEvo,  GramFE_TimeEvo_MemSize )  );
#  endif

// create streams
   Stream = new cudaStream_t [GPU_NStream];
   for (int s=0; s<GPU_NStream; s++)      CUDA_CHECK_ERROR(  cudaStreamCreate( &Stream[s] )  );


   return GAMER_SUCCESS;

} // FUNCTION : CUAPI_MemAllocate_Fluid



#endif // #ifdef GPU
