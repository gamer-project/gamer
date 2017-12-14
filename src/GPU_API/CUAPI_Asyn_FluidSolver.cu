#include "GAMER.h"
#include "CUFLU.h"

#ifdef GPU



#if   ( MODEL == HYDRO )
#if   ( FLU_SCHEME == RTVD )
__global__ void CUFLU_FluidSolver_RTVD( real g_Fluid_In []   [NCOMP_TOTAL][ FLU_NXT*FLU_NXT*FLU_NXT ],
                                        real g_Fluid_Out[]   [NCOMP_TOTAL][ PS2*PS2*PS2 ],
                                        real g_Flux     [][9][NCOMP_TOTAL][ PS2*PS2 ],
                                        const double g_Corner[][3],
                                        const real g_Pot_USG[][ USG_NXT_F*USG_NXT_F*USG_NXT_F ],
                                        const real dt, const real _dh, const real Gamma, const bool StoreFlux,
                                        const bool XYZ, const real MinDens, const real MinPres );
#elif ( FLU_SCHEME == WAF )
__global__ void CUFLU_FluidSolver_WAF( real g_Fluid_In []   [NCOMP_TOTAL][ FLU_NXT*FLU_NXT*FLU_NXT ],
                                       real g_Fluid_Out[]   [NCOMP_TOTAL][ PS2*PS2*PS2 ],
                                       real g_Flux     [][9][NCOMP_TOTAL][ PS2*PS2 ],
                                       const double g_Corner[][3],
                                       const real g_Pot_USG[][ USG_NXT_F*USG_NXT_F*USG_NXT_F ],
                                       const real dt, const real _dh, const real Gamma, const bool StoreFlux,
                                       const bool XYZ, const WAF_Limiter_t WAF_Limiter, const real MinDens, const real MinPres );
#elif ( FLU_SCHEME == MHM  ||  FLU_SCHEME == MHM_RP )
__global__ void CUFLU_FluidSolver_MHM( const real g_Fluid_In[]   [NCOMP_TOTAL][ FLU_NXT*FLU_NXT*FLU_NXT ],
                                       real g_Fluid_Out     []   [NCOMP_TOTAL][ PS2*PS2*PS2 ],
                                       char DE_Array_Out    []                [ PS2*PS2*PS2 ],
                                       real g_Flux          [][9][NCOMP_TOTAL][ PS2*PS2 ],
                                       const double g_Corner[][3],
                                       const real g_Pot_USG[] [ USG_NXT_F*USG_NXT_F*USG_NXT_F ],
                                       real g_PriVar     [][NCOMP_TOTAL][ FLU_NXT*FLU_NXT*FLU_NXT ],
                                       real g_Slope_PPM_x[][NCOMP_TOTAL][ N_SLOPE_PPM*N_SLOPE_PPM*N_SLOPE_PPM],
                                       real g_Slope_PPM_y[][NCOMP_TOTAL][ N_SLOPE_PPM*N_SLOPE_PPM*N_SLOPE_PPM],
                                       real g_Slope_PPM_z[][NCOMP_TOTAL][ N_SLOPE_PPM*N_SLOPE_PPM*N_SLOPE_PPM],
                                       real g_FC_Var_xL  [][NCOMP_TOTAL][ N_FC_VAR*N_FC_VAR*N_FC_VAR ],
                                       real g_FC_Var_xR  [][NCOMP_TOTAL][ N_FC_VAR*N_FC_VAR*N_FC_VAR ],
                                       real g_FC_Var_yL  [][NCOMP_TOTAL][ N_FC_VAR*N_FC_VAR*N_FC_VAR ],
                                       real g_FC_Var_yR  [][NCOMP_TOTAL][ N_FC_VAR*N_FC_VAR*N_FC_VAR ],
                                       real g_FC_Var_zL  [][NCOMP_TOTAL][ N_FC_VAR*N_FC_VAR*N_FC_VAR ],
                                       real g_FC_Var_zR  [][NCOMP_TOTAL][ N_FC_VAR*N_FC_VAR*N_FC_VAR ],
                                       real g_FC_Flux_x  [][NCOMP_TOTAL][ N_FC_FLUX*N_FC_FLUX*N_FC_FLUX ],
                                       real g_FC_Flux_y  [][NCOMP_TOTAL][ N_FC_FLUX*N_FC_FLUX*N_FC_FLUX ],
                                       real g_FC_Flux_z  [][NCOMP_TOTAL][ N_FC_FLUX*N_FC_FLUX*N_FC_FLUX ],
                                       const real dt, const real _dh, const real Gamma, const bool StoreFlux,
                                       const LR_Limiter_t LR_Limiter, const real MinMod_Coeff,
                                       const real EP_Coeff, const double Time, const OptGravityType_t GravityType,
                                       const real MinDens, const real MinPres, const real DualEnergySwitch,
                                       const bool NormPassive, const int NNorm,
                                       const bool JeansMinPres, const real JeansMinPres_Coeff );
#elif ( FLU_SCHEME == CTU )
__global__ void CUFLU_FluidSolver_CTU( const real g_Fluid_In[]   [NCOMP_TOTAL][ FLU_NXT*FLU_NXT*FLU_NXT ],
                                       real g_Fluid_Out     []   [NCOMP_TOTAL][ PS2*PS2*PS2 ],
                                       char DE_Array_Out    []                [ PS2*PS2*PS2 ],
                                       real g_Flux          [][9][NCOMP_TOTAL][ PS2*PS2 ],
                                       const double g_Corner[][3],
                                       const real g_Pot_USG[] [ USG_NXT_F*USG_NXT_F*USG_NXT_F ],
                                       real g_PriVar     [][NCOMP_TOTAL][ FLU_NXT*FLU_NXT*FLU_NXT ],
                                       real g_Slope_PPM_x[][NCOMP_TOTAL][ N_SLOPE_PPM*N_SLOPE_PPM*N_SLOPE_PPM],
                                       real g_Slope_PPM_y[][NCOMP_TOTAL][ N_SLOPE_PPM*N_SLOPE_PPM*N_SLOPE_PPM],
                                       real g_Slope_PPM_z[][NCOMP_TOTAL][ N_SLOPE_PPM*N_SLOPE_PPM*N_SLOPE_PPM],
                                       real g_FC_Var_xL  [][NCOMP_TOTAL][ N_FC_VAR*N_FC_VAR*N_FC_VAR ],
                                       real g_FC_Var_xR  [][NCOMP_TOTAL][ N_FC_VAR*N_FC_VAR*N_FC_VAR ],
                                       real g_FC_Var_yL  [][NCOMP_TOTAL][ N_FC_VAR*N_FC_VAR*N_FC_VAR ],
                                       real g_FC_Var_yR  [][NCOMP_TOTAL][ N_FC_VAR*N_FC_VAR*N_FC_VAR ],
                                       real g_FC_Var_zL  [][NCOMP_TOTAL][ N_FC_VAR*N_FC_VAR*N_FC_VAR ],
                                       real g_FC_Var_zR  [][NCOMP_TOTAL][ N_FC_VAR*N_FC_VAR*N_FC_VAR ],
                                       real g_FC_Flux_x  [][NCOMP_TOTAL][ N_FC_FLUX*N_FC_FLUX*N_FC_FLUX ],
                                       real g_FC_Flux_y  [][NCOMP_TOTAL][ N_FC_FLUX*N_FC_FLUX*N_FC_FLUX ],
                                       real g_FC_Flux_z  [][NCOMP_TOTAL][ N_FC_FLUX*N_FC_FLUX*N_FC_FLUX ],
                                       const real dt, const real _dh, const real Gamma, const bool StoreFlux,
                                       const LR_Limiter_t LR_Limiter, const real MinMod_Coeff,
                                       const real EP_Coeff, const double Time, const OptGravityType_t GravityType,
                                       const real MinDens, const real MinPres, const real DualEnergySwitch,
                                       const bool NormPassive, const int NNorm,
                                       const bool JeansMinPres, const real JeansMinPres_Coeff );
#endif // FLU_SCHEME

#elif ( MODEL == MHD )
#warning : WAIT MHD !!!

#elif ( MODEL == ELBDM )
__global__ void CUFLU_ELBDMSolver( real g_Fluid_In [][FLU_NIN ][ FLU_NXT*FLU_NXT*FLU_NXT ],
                                   real g_Fluid_Out[][FLU_NOUT][ PS2*PS2*PS2 ],
                                   real g_Flux     [][9][NFLUX_TOTAL][ PS2*PS2 ],
                                   const real dt, const real _dh, const real Eta, const bool StoreFlux,
                                   const real Taylor3_Coeff, const bool XYZ, const real MinDens );

#else
#error : ERROR : unsupported MODEL !!
#endif // MODEL


// device pointers
extern real (*d_Flu_Array_F_In )[FLU_NIN ][ FLU_NXT*FLU_NXT*FLU_NXT ];
extern real (*d_Flu_Array_F_Out)[FLU_NOUT][ PS2*PS2*PS2 ];
extern real (*d_Flux_Array)[9][NFLUX_TOTAL][ PS2*PS2 ];
extern double (*d_Corner_Array_F)[3];
#if ( MODEL == HYDRO  ||  MODEL == MHD )
#ifdef DUAL_ENERGY
extern char (*d_DE_Array_F_Out)[ PS2*PS2*PS2 ];
#else
static char (*d_DE_Array_F_Out)[ PS2*PS2*PS2 ] = NULL;
#endif
#endif
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
#endif // FLU_SCHEME
#elif ( MODEL == MHD )
#warning : WAIT MHD !!!
#endif // MODEL

#ifdef UNSPLIT_GRAVITY
extern real (*d_Pot_Array_USG_F)[ USG_NXT_F*USG_NXT_F*USG_NXT_F ];
#else
#if ( MODEL == HYDRO  ||  MODEL == MHD )
static real (*d_Pot_Array_USG_F)[ USG_NXT_F*USG_NXT_F*USG_NXT_F ] = NULL;
#endif
#endif // #ifdef UNSPLIT_GRAVITY ... else ...

extern cudaStream_t *Stream;




//-------------------------------------------------------------------------------------------------------
// Function    :  CUAPI_Asyn_FluidSolver
// Description :  1. MODEL == HYDRO : use GPU to solve the Euler equations by different schemes
//                                    --> invoke the kernel "CUFLU_FluidSolver_XXX"
//                2. MODEL == ELBDM : use GPU to solve the kinematic operator in the Schrodinger's equations
//                                    --> invoke the kernel "CUFLU_ELBDMSolver"
//
//                ***********************************************************
//                **                Asynchronous Function                  **
//                **                                                       **
//                **  will return before the execution in GPU is complete  **
//                ***********************************************************
//
// Note        :  1. Use streams for the asychronous memory copy between device and host
//                2. Prefix "d" : for pointers pointing to the "Device" memory space
//                   Prefix "h" : for pointers pointing to the "Host"   memory space
//                3. Use the input pamameter "XYZ" to control the order of update for dimensional-splitting
//                   methods (RTVD/WAF)
//                4. Currently five hydro schemes are supported :
//                   1. Relaxing TVD scheme                            (RTVD  ) -->   split
//                   2. Weighted-Average-Flux scheme                   (WAF   ) -->   split
//                   3. MUSCL-Hancock scheme                           (MHM   ) --> unsplit
//                   4. MUSCL-Hancock scheme with Riemann prediction   (MHM_RP) --> unsplit
//                   5. Corner-Transport-Upwind scheme                 (CTU   ) --> unsplit
//
// Parameter   :  h_Flu_Array_In       : Host array to store the input fluid variables
//                h_Flu_Array_Out      : Host array to store the output fluid variables
//                h_DE_Array_Out       : Host array to store the dual-energy status
//                h_Flux_Array         : Host array to store the output fluxes
//                h_Corner_Array       : Host array storing the physical corner coordinates of each patch group
//                h_Pot_Array_USG      : Host array storing the input potential for UNSPLIT_GRAVITY
//                NPatchGroup          : Number of patch groups evaluated simultaneously by GPU
//                dt                   : Time interval to advance solution
//                dh                   : Grid size
//                Gamma                : Ratio of specific heats
//                StoreFlux            : true --> store the coarse-fine fluxes
//                XYZ                  : true  : x->y->z ( forward sweep)
//                                       false : z->y->x (backward sweep)
//                                       ~ useless in directionally unsplit schemes
//                LR_Limiter           : Slope limiter for the data reconstruction in the MHM/MHM_RP/CTU schemes
//                                       (0/1/2/3/4) = (vanLeer/generalized MinMod/vanAlbada/
//                                                      vanLeer + generalized MinMod/extrema-preserving) limiter
//                MinMod_Coeff         : Coefficient of the generalized MinMod limiter
//                EP_Coeff             : Coefficient of the extrema-preserving limiter
//                WAF_Limiter          : Flux limiter for the WAF scheme
//                                       (0/1/2/3) = (SuperBee/vanLeer/vanAlbada/MinBee)
//                ELBDM_Eta            : Particle mass / Planck constant
//                ELBDM_Taylor3_Coeff  : Coefficient in front of the third term in the Taylor expansion for ELBDM
//                ELBDM_Taylor3_Auto   : true --> Determine ELBDM_Taylor3_Coeff automatically by invoking the
//                                                function "ELBDM_SetTaylor3Coeff"
//                Time                 : Current physical time                                     (for UNSPLIT_GRAVITY only)
//                GravityType          : Types of gravity --> self-gravity, external gravity, both (for UNSPLIT_GRAVITY only)
//                GPU_NStream          : Number of CUDA streams for the asynchronous memory copy
//                MinDens/Pres         : Minimum allowed density and pressure
//                DualEnergySwitch     : Use the dual-energy formalism if E_int/E_kin < DualEnergySwitch
//                NormPassive          : true --> normalize passive scalars so that the sum of their mass density
//                                                is equal to the gas mass density
//                NNorm                : Number of passive scalars to be normalized
//                                       --> Should be set to the global variable "PassiveNorm_NVar"
//                JeansMinPres         : Apply minimum pressure estimated from the Jeans length
//                JeansMinPres_Coeff   : Coefficient used by JeansMinPres = G*(Jeans_NCell*Jeans_dh)^2/(Gamma*pi);
//
// Useless parameters in HYDRO : ELBDM_Eta
// Useless parameters in ELBDM : h_Flux_Array, Gamma, LR_Limiter, MinMod_Coeff, EP_Coeff, WAF_Limite, MinPres
//-------------------------------------------------------------------------------------------------------
void CUAPI_Asyn_FluidSolver( real h_Flu_Array_In [][FLU_NIN    ][ FLU_NXT*FLU_NXT*FLU_NXT ],
                             real h_Flu_Array_Out[][FLU_NOUT   ][ PS2*PS2*PS2 ],
                             char h_DE_Array_Out[][ PS2*PS2*PS2 ],
                             real h_Flux_Array[][9][NFLUX_TOTAL][ PS2*PS2 ],
                             const double h_Corner_Array[][3],
                             real h_Pot_Array_USG[][USG_NXT_F][USG_NXT_F][USG_NXT_F],
                             const int NPatchGroup, const real dt, const real dh, const real Gamma, const bool StoreFlux,
                             const bool XYZ, const LR_Limiter_t LR_Limiter, const real MinMod_Coeff, const real EP_Coeff,
                             const WAF_Limiter_t WAF_Limiter, const real ELBDM_Eta, real ELBDM_Taylor3_Coeff,
                             const bool ELBDM_Taylor3_Auto, const double Time, const OptGravityType_t GravityType,
                             const int GPU_NStream, const real MinDens, const real MinPres, const real DualEnergySwitch,
                             const bool NormPassive, const int NNorm,
                             const bool JeansMinPres, const real JeansMinPres_Coeff )
{

// check
#  ifdef GAMER_DEBUG
#  if   ( MODEL == HYDRO )
   if ( LR_Limiter != VANLEER  &&  LR_Limiter != GMINMOD  &&  LR_Limiter != ALBADA  &&  LR_Limiter != EXTPRE  &&
        LR_Limiter != VL_GMINMOD  &&  LR_Limiter != LR_LIMITER_NONE )
      Aux_Error( ERROR_INFO, "unsupported limiter (%d) !!\n", LR_Limiter );

#  ifdef UNSPLIT_GRAVITY
   if ( GravityType == GRAVITY_SELF  ||  GravityType == GRAVITY_BOTH )
   {
      if ( h_Pot_Array_USG   == NULL )   Aux_Error( ERROR_INFO, "h_Pot_Array_USG == NULL !!\n" );
      if ( d_Pot_Array_USG_F == NULL )   Aux_Error( ERROR_INFO, "d_Pot_Array_USG_F == NULL !!\n" );
   }

   if ( GravityType == GRAVITY_EXTERNAL  ||  GravityType == GRAVITY_BOTH )
   {
      if ( h_Corner_Array   == NULL )    Aux_Error( ERROR_INFO, "h_Corner_Array == NULL !!\n" );
      if ( d_Corner_Array_F == NULL )    Aux_Error( ERROR_INFO, "d_Corner_Array_F == NULL !!\n" );
   }
#  endif

#  elif ( MODEL == MHD )
#  warning : WAIT MHD !!!

#  elif ( MODEL == ELBDM )

#  else
#  warning : DO YOU WANT TO ADD SOMETHING HERE FOR THE NEW MODEL ??
#  endif

   if ( StoreFlux )
   {
      if ( d_Flux_Array == NULL )   Aux_Error( ERROR_INFO, "device Flux_Array is not allocated !!\n" );
      if ( h_Flux_Array == NULL )   Aux_Error( ERROR_INFO, "host Flux_Array is not allocated !!\n" );
   }
#  endif // #ifdef GAMER_DEBUG


   const real _dh = (real)1.0/dh;
   const dim3 BlockDim_FluidSolver ( FLU_BLOCK_SIZE_X, FLU_BLOCK_SIZE_Y, 1 ); // for the fluidsolvers

// model-dependent operations
#  if   ( MODEL == HYDRO )

#  elif ( MODEL == MHD )
#  warning : WAIT MHD !!!

#  elif ( MODEL == ELBDM )
// evaluate the optimized Taylor expansion coefficient
   if ( ELBDM_Taylor3_Auto )  ELBDM_Taylor3_Coeff = ELBDM_SetTaylor3Coeff( dt, dh, ELBDM_Eta );

#  else
#  error : ERROR : unsupported MODEL !!
#  endif // MODEL

   int *NPatch_per_Stream  = new int [GPU_NStream];
   int *UsedPatch          = new int [GPU_NStream];
   int *Flu_MemSize_In     = new int [GPU_NStream];
   int *Flu_MemSize_Out    = new int [GPU_NStream];
   int *Flux_MemSize       = new int [GPU_NStream];
#  ifdef UNSPLIT_GRAVITY
   int *USG_MemSize        = new int [GPU_NStream];
   int *Corner_MemSize     = new int [GPU_NStream];
#  endif
#  ifdef DUAL_ENERGY
   int *DE_MemSize_Out     = new int [GPU_NStream];
#  endif


// set the number of patches of each stream
   UsedPatch[0] = 0;

   if ( GPU_NStream == 1 )    NPatch_per_Stream[0] = NPatchGroup;
   else
   {
      for (int s=0; s<GPU_NStream-1; s++)
      {
         NPatch_per_Stream[s] = NPatchGroup / GPU_NStream;
         UsedPatch[s+1] = UsedPatch[s] + NPatch_per_Stream[s];
      }

      NPatch_per_Stream[GPU_NStream-1] = NPatchGroup - UsedPatch[GPU_NStream-1];
   }


// set the size of data to be transferred into GPU in each stream
   for (int s=0; s<GPU_NStream; s++)
   {
      Flu_MemSize_In [s] = sizeof(real  )*NPatch_per_Stream[s]*FLU_NIN *CUBE(FLU_NXT);
      Flu_MemSize_Out[s] = sizeof(real  )*NPatch_per_Stream[s]*FLU_NOUT*CUBE(PS2);
      Flux_MemSize   [s] = sizeof(real  )*NPatch_per_Stream[s]*NFLUX_TOTAL*9*PS2*PS2;
#     ifdef UNSPLIT_GRAVITY
      USG_MemSize    [s] = sizeof(real  )*NPatch_per_Stream[s]*CUBE(USG_NXT_F);
      Corner_MemSize [s] = sizeof(double)*NPatch_per_Stream[s]*3;
#     endif
#     ifdef DUAL_ENERGY
      DE_MemSize_Out [s] = sizeof(char  )*NPatch_per_Stream[s]*CUBE(PS2);
#     endif
   }


// a. copy data from host to device
//=========================================================================================
   for (int s=0; s<GPU_NStream; s++)
   {
      if ( NPatch_per_Stream[s] == 0 )    continue;

      CUDA_CHECK_ERROR(  cudaMemcpyAsync( d_Flu_Array_F_In  + UsedPatch[s], h_Flu_Array_In  + UsedPatch[s],
                         Flu_MemSize_In[s], cudaMemcpyHostToDevice, Stream[s] )  );

#     ifdef UNSPLIT_GRAVITY
      CUDA_CHECK_ERROR(  cudaMemcpyAsync( d_Pot_Array_USG_F + UsedPatch[s], h_Pot_Array_USG + UsedPatch[s],
                         USG_MemSize   [s], cudaMemcpyHostToDevice, Stream[s] )  );

      if ( GravityType == GRAVITY_EXTERNAL  ||  GravityType == GRAVITY_BOTH )
      CUDA_CHECK_ERROR(  cudaMemcpyAsync( d_Corner_Array_F  + UsedPatch[s], h_Corner_Array  + UsedPatch[s],
                         Corner_MemSize[s], cudaMemcpyHostToDevice, Stream[s] )  );
#     endif
   } // for (int s=0; s<GPU_NStream; s++)


// b. execute the kernel
//=========================================================================================
   for (int s=0; s<GPU_NStream; s++)
   {
      if ( NPatch_per_Stream[s] == 0 )    continue;

#     if   ( MODEL == HYDRO )

#        if   ( FLU_SCHEME == RTVD )

         CUFLU_FluidSolver_RTVD <<< NPatch_per_Stream[s], BlockDim_FluidSolver, 0, Stream[s] >>>
                                ( d_Flu_Array_F_In  + UsedPatch[s],
                                  d_Flu_Array_F_Out + UsedPatch[s],
                                  d_Flux_Array      + UsedPatch[s],
                                  d_Corner_Array_F  + UsedPatch[s],
                                  d_Pot_Array_USG_F + UsedPatch[s],
                                  dt, _dh, Gamma, StoreFlux, XYZ, MinDens, MinPres );

#        elif ( FLU_SCHEME == WAF )

         CUFLU_FluidSolver_WAF <<< NPatch_per_Stream[s], BlockDim_FluidSolver, 0, Stream[s] >>>
                               ( d_Flu_Array_F_In  + UsedPatch[s],
                                 d_Flu_Array_F_Out + UsedPatch[s],
                                 d_Flux_Array      + UsedPatch[s],
                                 d_Corner_Array_F  + UsedPatch[s],
                                 d_Pot_Array_USG_F + UsedPatch[s],
                                 dt, _dh, Gamma, StoreFlux, XYZ, WAF_Limiter, MinDens, MinPres );

#        elif ( FLU_SCHEME == MHM  ||  FLU_SCHEME == MHM_RP )

         CUFLU_FluidSolver_MHM <<< NPatch_per_Stream[s], BlockDim_FluidSolver, 0, Stream[s] >>>
                               ( d_Flu_Array_F_In  + UsedPatch[s],
                                 d_Flu_Array_F_Out + UsedPatch[s],
                                 d_DE_Array_F_Out  + UsedPatch[s],
                                 d_Flux_Array      + UsedPatch[s],
                                 d_Corner_Array_F  + UsedPatch[s],
                                 d_Pot_Array_USG_F + UsedPatch[s],
                                 d_PriVar          + UsedPatch[s],
                                 d_Slope_PPM_x     + UsedPatch[s],
                                 d_Slope_PPM_y     + UsedPatch[s],
                                 d_Slope_PPM_z     + UsedPatch[s],
                                 d_FC_Var_xL       + UsedPatch[s],
                                 d_FC_Var_xR       + UsedPatch[s],
                                 d_FC_Var_yL       + UsedPatch[s],
                                 d_FC_Var_yR       + UsedPatch[s],
                                 d_FC_Var_zL       + UsedPatch[s],
                                 d_FC_Var_zR       + UsedPatch[s],
                                 d_FC_Flux_x       + UsedPatch[s],
                                 d_FC_Flux_y       + UsedPatch[s],
                                 d_FC_Flux_z       + UsedPatch[s],
                                 dt, _dh, Gamma, StoreFlux, LR_Limiter, MinMod_Coeff, EP_Coeff,
                                 Time, GravityType, MinDens, MinPres, DualEnergySwitch, NormPassive, NNorm,
                                 JeansMinPres, JeansMinPres_Coeff );

#        elif ( FLU_SCHEME == CTU )

         CUFLU_FluidSolver_CTU <<< NPatch_per_Stream[s], BlockDim_FluidSolver, 0, Stream[s] >>>
                               ( d_Flu_Array_F_In  + UsedPatch[s],
                                 d_Flu_Array_F_Out + UsedPatch[s],
                                 d_DE_Array_F_Out  + UsedPatch[s],
                                 d_Flux_Array      + UsedPatch[s],
                                 d_Corner_Array_F  + UsedPatch[s],
                                 d_Pot_Array_USG_F + UsedPatch[s],
                                 d_PriVar          + UsedPatch[s],
                                 d_Slope_PPM_x     + UsedPatch[s],
                                 d_Slope_PPM_y     + UsedPatch[s],
                                 d_Slope_PPM_z     + UsedPatch[s],
                                 d_FC_Var_xL       + UsedPatch[s],
                                 d_FC_Var_xR       + UsedPatch[s],
                                 d_FC_Var_yL       + UsedPatch[s],
                                 d_FC_Var_yR       + UsedPatch[s],
                                 d_FC_Var_zL       + UsedPatch[s],
                                 d_FC_Var_zR       + UsedPatch[s],
                                 d_FC_Flux_x       + UsedPatch[s],
                                 d_FC_Flux_y       + UsedPatch[s],
                                 d_FC_Flux_z       + UsedPatch[s],
                                 dt, _dh, Gamma, StoreFlux, LR_Limiter, MinMod_Coeff, EP_Coeff,
                                 Time, GravityType, MinDens, MinPres, DualEnergySwitch, NormPassive, NNorm,
                                 JeansMinPres, JeansMinPres_Coeff );

#        else

#        error : unsupported GPU hydro scheme

#        endif // FLU_SCHEME

#     elif ( MODEL == MHD )
#     warning :: WAIT MHD !!!

#     elif ( MODEL == ELBDM )

         CUFLU_ELBDMSolver <<< NPatch_per_Stream[s], BlockDim_FluidSolver, 0, Stream[s] >>>
                               ( d_Flu_Array_F_In  + UsedPatch[s],
                                 d_Flu_Array_F_Out + UsedPatch[s],
                                 d_Flux_Array      + UsedPatch[s],
                                 dt, _dh, ELBDM_Eta, StoreFlux, ELBDM_Taylor3_Coeff, XYZ, MinDens );

#     else

#        error : unsupported MODEL !!

#     endif // MODEL


      CUDA_CHECK_ERROR( cudaGetLastError() );
   } // for (int s=0; s<GPU_NStream; s++)


// c. copy data from device to host
//=========================================================================================
   for (int s=0; s<GPU_NStream; s++)
   {
      if ( NPatch_per_Stream[s] == 0 )    continue;

      CUDA_CHECK_ERROR(  cudaMemcpyAsync( h_Flu_Array_Out + UsedPatch[s], d_Flu_Array_F_Out + UsedPatch[s],
                         Flu_MemSize_Out[s], cudaMemcpyDeviceToHost, Stream[s] )  );

      if ( StoreFlux )
      CUDA_CHECK_ERROR(  cudaMemcpyAsync( h_Flux_Array    + UsedPatch[s], d_Flux_Array      + UsedPatch[s],
                         Flux_MemSize[s],    cudaMemcpyDeviceToHost, Stream[s] )  );

#     if   ( MODEL == HYDRO )
#     ifdef DUAL_ENERGY
      CUDA_CHECK_ERROR(  cudaMemcpyAsync( h_DE_Array_Out  + UsedPatch[s], d_DE_Array_F_Out  + UsedPatch[s],
                         DE_MemSize_Out[s],  cudaMemcpyDeviceToHost, Stream[s] )  );
#     endif

#     elif ( MODEL == MHD )
#     warning : WAIT MHD !!!
#     endif // MODEL
   } // for (int s=0; s<GPU_NStream; s++)


   delete [] NPatch_per_Stream;
   delete [] UsedPatch;
   delete [] Flu_MemSize_In;
   delete [] Flu_MemSize_Out;
   delete [] Flux_MemSize;
#  ifdef UNSPLIT_GRAVITY
   delete [] USG_MemSize;
   delete [] Corner_MemSize;
#  endif
#  ifdef DUAL_ENERGY
   delete [] DE_MemSize_Out;
#  endif

} // FUNCTION : CUAPI_Asyn_FluidSolver



#endif // #ifdef GPU
