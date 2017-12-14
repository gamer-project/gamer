#include "Macro.h"
#include "CUFLU.h"

#if ( defined GPU  &&  MODEL == HYDRO  &&  FLU_SCHEME == CTU )



#include "CUFLU_Shared_FluUtility.cu"
#include "CUFLU_Shared_DataReconstruction.cu"
#include "CUFLU_Shared_ComputeFlux.cu"
#include "CUFLU_Shared_FullStepUpdate.cu"

static __device__ void CUFLU_TGradient_Correction( real g_FC_Var_xL[][NCOMP_TOTAL][ N_FC_VAR*N_FC_VAR*N_FC_VAR ],
                                                   real g_FC_Var_xR[][NCOMP_TOTAL][ N_FC_VAR*N_FC_VAR*N_FC_VAR ],
                                                   real g_FC_Var_yL[][NCOMP_TOTAL][ N_FC_VAR*N_FC_VAR*N_FC_VAR ],
                                                   real g_FC_Var_yR[][NCOMP_TOTAL][ N_FC_VAR*N_FC_VAR*N_FC_VAR ],
                                                   real g_FC_Var_zL[][NCOMP_TOTAL][ N_FC_VAR*N_FC_VAR*N_FC_VAR ],
                                                   real g_FC_Var_zR[][NCOMP_TOTAL][ N_FC_VAR*N_FC_VAR*N_FC_VAR ],
                                                   const real g_FC_Flux_x[][NCOMP_TOTAL][ N_FC_FLUX*N_FC_FLUX*N_FC_FLUX ],
                                                   const real g_FC_Flux_y[][NCOMP_TOTAL][ N_FC_FLUX*N_FC_FLUX*N_FC_FLUX ],
                                                   const real g_FC_Flux_z[][NCOMP_TOTAL][ N_FC_FLUX*N_FC_FLUX*N_FC_FLUX ],
                                                   const real dt, const real _dh, const real Gamma,
                                                   const real MinDens, const real MinPres );




#ifdef UNSPLIT_GRAVITY
#include "CUPOT.h"
__constant__ double ExtAcc_AuxArray_d_Flu[EXT_ACC_NAUX_MAX];

//-------------------------------------------------------------------------------------------------------
// Function    :  CUFLU_FluidSolver_SetConstMem_ExtAcc
// Description :  Set the constant memory of ExtAcc_AuxArray_d_Flu used by CUFLU_FluidSolver_CTU/MHM
//
// Note        :  Adopt the suggested approach for CUDA version >= 5.0
//
// Parameter   :  None
//
// Return      :  0/-1 : successful/failed
//---------------------------------------------------------------------------------------------------
int CUFLU_FluidSolver_SetConstMem_ExtAcc( double ExtAcc_AuxArray_h[] )
{

   if (  cudaSuccess != cudaMemcpyToSymbol( ExtAcc_AuxArray_d_Flu, ExtAcc_AuxArray_h, EXT_ACC_NAUX_MAX*sizeof(double),
                                            0, cudaMemcpyHostToDevice)  )
      return -1;

   else
      return 0;

} // FUNCTION : CUFLU_FluidSolver_SetConstMem_ExtAcc
#endif



#if ( NCOMP_PASSIVE > 0 )
__constant__ int NormIdx_d[NCOMP_PASSIVE];

//-------------------------------------------------------------------------------------------------------
// Function    :  CUFLU_FluidSolver_SetConstMem_NormIdx
// Description :  Set the constant memory of NormIdx_d used by CUFLU_FluidSolver_CTU/MHM
//
// Note        :  Adopt the suggested approach for CUDA version >= 5.0
//
// Parameter   :  None
//
// Return      :  0/-1 : successful/failed
//---------------------------------------------------------------------------------------------------
int CUFLU_FluidSolver_SetConstMem_NormIdx( int NormIdx_h[] )
{

   if (  cudaSuccess != cudaMemcpyToSymbol( NormIdx_d, NormIdx_h, NCOMP_PASSIVE*sizeof(int),
                                            0, cudaMemcpyHostToDevice)  )
      return -1;

   else
      return 0;

} // FUNCTION : CUFLU_FluidSolver_SetConstMem_NormIdx

#else
__constant__ int *NormIdx_d = NULL;

#endif // #if ( NCOMP_PASSIVE > 0 ) ... else ...




//-------------------------------------------------------------------------------------------------------
// Function    :  CUFLU_FluidSolver_CTU
// Description :  GPU fluid solver based on the Corner-Transport-Upwind (CTU) scheme
//
// Note        :  1. Prefix "g" for pointers pointing to the "Global" memory space
//                   Prefix "s" for pointers pointing to the "Shared" memory space
//                2. Ref : Stone et al., ApJS, 178, 137 (2008)
//                3. Each patch group requires about 2.9*10^7 flops with CTU + PPM + Roe solver
//                   --> 206 GFLOPS is achieved in one C2050 GPU
//
// Parameter   :  g_Fluid_In         : Global memory array storing the input fluid variables
//                g_Fluid_Out        : Global memory array to store the output fluid variables
//                g_DE_Out           : Global memory array to store the output dual-energy status
//                g_Flux             : Global memory array to store the output fluxes
//                g_Corner           : Global memory array storing the physical corner coordinates of each patch group (USELESS CURRENTLY)
//                g_Pot_USG          : Global memory array storing the input potential for UNSPLIT_GRAVITY (NOT SUPPORTED in RTVD)
//                g_PriVar           : Global memory array to store the primitive variables
//                g_Slope_PPM_x      : Global memory array to store the x-slope for the PPM reconstruction
//                g_Slope_PPM_y      : Global memory array to store the y-slope for the PPM reconstruction
//                g_Slope_PPM_z      : Global memory array to store the z-slope for the PPM reconstruction
//                g_FC_Var_xL        : Global memory array to store the half-step variables on the -x surface
//                g_FC_Var_xR        : Global memory array to store the half-step variables on the +x surface
//                g_FC_Var_yL        : Global memory array to store the half-step variables on the -y surface
//                g_FC_Var_yR        : Global memory array to store the half-step variables on the +y surface
//                g_FC_Var_zL        : Global memory array to store the half-step variables on the -z surface
//                g_FC_Var_zR        : Global memory array to store the half-step variables on the +z surface
//                g_FC_Flux_x        : Global memory array to store the face-centered fluxes in the x direction
//                g_FC_Flux_y        : Global memory array to store the face-centered fluxes in the y direction
//                g_FC_Flux_z        : Global memory array to store the face-centered fluxes in the z direction
//                dt                 : Time interval to advance solution
//                _dh                : 1 / grid size
//                Gamma              : Ratio of specific heats
//                StoreFlux          : true --> store the coarse-fine fluxes
//                LR_Limiter         : Slope limiter for the data reconstruction in the MHM/MHM_RP/CTU schemes
//                                     (0/1/2/3/4) = (vanLeer/generalized MinMod/vanAlbada/
//                                                    vanLeer + generalized MinMod/extrema-preserving) limiter
//                MinMod_Coeff       : Coefficient of the generalized MinMod limiter
//                EP_Coeff           : Coefficient of the extrema-preserving limiter
//                Time               : Current physical time                                     (for UNSPLIT_GRAVITY only)
//                GravityType        : Types of gravity --> self-gravity, external gravity, both (for UNSPLIT_GRAVITY only)
//                MinDens/Pres       : Minimum allowed density and pressure
//                DualEnergySwitch   : Use the dual-energy formalism if E_int/E_kin < DualEnergySwitch
//                NormPassive        : true --> normalize passive scalars so that the sum of their mass density
//                                              is equal to the gas mass density
//                NNorm              : Number of passive scalars to be normalized
//                                     --> Should be set to the global variable "PassiveNorm_NVar"
//                JeansMinPres       : Apply minimum pressure estimated from the Jeans length
//                JeansMinPres_Coeff : Coefficient used by JeansMinPres = G*(Jeans_NCell*Jeans_dh)^2/(Gamma*pi);
//-------------------------------------------------------------------------------------------------------
__global__ void CUFLU_FluidSolver_CTU( const real g_Fluid_In[]   [NCOMP_TOTAL][ FLU_NXT*FLU_NXT*FLU_NXT ],
                                       real g_Fluid_Out     []   [NCOMP_TOTAL][ PS2*PS2*PS2 ],
                                       char g_DE_Out        []                [ PS2*PS2*PS2 ],
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
                                       const bool JeansMinPres, const real JeansMinPres_Coeff )
{

#  ifdef UNSPLIT_GRAVITY
   const bool CorrHalfVel_Yes = true;
#  endif
   const bool CorrHalfVel_No  = false;

// 1. conserved variables --> primitive variables
   CUFLU_Con2Pri_AllGrids( g_Fluid_In, g_PriVar, Gamma, MinPres, NormPassive, NNorm, NormIdx_d,
                           JeansMinPres, JeansMinPres_Coeff );
   __syncthreads();


// 2. evaluate the half-step face-centered solution
   CUFLU_DataReconstruction( g_PriVar, g_Slope_PPM_x, g_Slope_PPM_y, g_Slope_PPM_z, g_FC_Var_xL, g_FC_Var_xR,
                             g_FC_Var_yL, g_FC_Var_yR, g_FC_Var_zL, g_FC_Var_zR, FLU_NXT, FLU_GHOST_SIZE-1,
                             Gamma, LR_Limiter, MinMod_Coeff, EP_Coeff, dt, _dh, MinDens, MinPres,
                             NormPassive, NNorm, NormIdx_d );
   __syncthreads();


// 3. evaluate the face-centered half-step fluxes by solving the Riemann problem
   CUFLU_ComputeFlux( g_FC_Var_xL, g_FC_Var_xR, g_FC_Var_yL, g_FC_Var_yR, g_FC_Var_zL, g_FC_Var_zR,
                      g_FC_Flux_x, g_FC_Flux_y, g_FC_Flux_z, NULL, false, 0, Gamma,
                      CorrHalfVel_No, NULL, NULL, NULL_REAL, NULL_REAL, NULL_REAL, GRAVITY_NONE, NULL, MinPres );
   __syncthreads();


// 4. correct the face-centered variables by the transverse flux gradients
   CUFLU_TGradient_Correction( g_FC_Var_xL, g_FC_Var_xR, g_FC_Var_yL, g_FC_Var_yR, g_FC_Var_zL, g_FC_Var_zR,
                               g_FC_Flux_x, g_FC_Flux_y, g_FC_Flux_z, dt, _dh, Gamma, MinDens, MinPres );
   __syncthreads();


// 5. evaluate the face-centered full-step fluxes by solving the Riemann problem with the corrected data
#  ifdef UNSPLIT_GRAVITY
   CUFLU_ComputeFlux( g_FC_Var_xL, g_FC_Var_xR, g_FC_Var_yL, g_FC_Var_yR, g_FC_Var_zL, g_FC_Var_zR,
                      g_FC_Flux_x, g_FC_Flux_y, g_FC_Flux_z, g_Flux, StoreFlux, 1, Gamma,
                      CorrHalfVel_Yes, g_Pot_USG, g_Corner, dt, _dh, Time, GravityType, ExtAcc_AuxArray_d_Flu, MinPres );
#  else
   CUFLU_ComputeFlux( g_FC_Var_xL, g_FC_Var_xR, g_FC_Var_yL, g_FC_Var_yR, g_FC_Var_zL, g_FC_Var_zR,
                      g_FC_Flux_x, g_FC_Flux_y, g_FC_Flux_z, g_Flux, StoreFlux, 1, Gamma,
                      CorrHalfVel_No, NULL, NULL, NULL_REAL, NULL_REAL, NULL_REAL, GRAVITY_NONE, NULL, MinPres );
#  endif

   __syncthreads();


// 6. evaluate the full-step solution
   CUFLU_FullStepUpdate( g_Fluid_In, g_Fluid_Out, g_DE_Out, g_FC_Flux_x, g_FC_Flux_y, g_FC_Flux_z,
                         dt, _dh, Gamma, MinDens, MinPres, DualEnergySwitch, NormPassive, NNorm, NormIdx_d );

} // FUNCTION : CUFLU_FluidSolver_CTU



//-------------------------------------------------------------------------------------------------------
// Function    :  CUFLU_TGradient_Correction
// Description :  Correct the face-centered variables by the transverse flux gradients
//
// Note        :  1. Prefix "g" for pointers pointing to the "Global" memory space
//                   Prefix "s" for pointers pointing to the "Shared" memory space
//                2. The function is asynchronous
//                   --> "__syncthreads" must be called before using the output data
//                3. The sizes of the arrays g_FC_Var_x/y/z and g_FC_Flux_XX in each direction are assumed
//                   to be "N_FC_VAR" (which is always PS2+2 in the currennt cases)
//
// Parameter   :  g_FC_Var_xL  : Global memory array to store the input and output -x face-centered conserved variables
//                g_FC_Var_xR  : Global memory array to store the input and output +x face-centered conserved variables
//                g_FC_Var_yL  : Global memory array to store the input and output -y face-centered conserved variables
//                g_FC_Var_yR  : Global memory array to store the input and output +y face-centered conserved variables
//                g_FC_Var_zL  : Global memory array to store the input and output -z face-centered conserved variables
//                g_FC_Var_zR  : Global memory array to store the input and output +z face-centered conserved variables
//                g_FC_Flux_x  : Global memory array storing the face-centered fluxes in the x direction
//                g_FC_Flux_y  : Global memory array storing the face-centered fluxes in the y direction
//                g_FC_Flux_z  : Global memory array storing the face-centered fluxes in the z direction
//                dt           : Time interval to advance solution
//                _dh          : 1 / grid size
//                Gamma        : Ratio of specific heats
//                MinDens/Pres : Minimum allowed density and pressure
//-------------------------------------------------------------------------------------------------------
__device__ void CUFLU_TGradient_Correction( real g_FC_Var_xL[][NCOMP_TOTAL][ N_FC_VAR*N_FC_VAR*N_FC_VAR ],
                                            real g_FC_Var_xR[][NCOMP_TOTAL][ N_FC_VAR*N_FC_VAR*N_FC_VAR ],
                                            real g_FC_Var_yL[][NCOMP_TOTAL][ N_FC_VAR*N_FC_VAR*N_FC_VAR ],
                                            real g_FC_Var_yR[][NCOMP_TOTAL][ N_FC_VAR*N_FC_VAR*N_FC_VAR ],
                                            real g_FC_Var_zL[][NCOMP_TOTAL][ N_FC_VAR*N_FC_VAR*N_FC_VAR ],
                                            real g_FC_Var_zR[][NCOMP_TOTAL][ N_FC_VAR*N_FC_VAR*N_FC_VAR ],
                                            const real g_FC_Flux_x[][NCOMP_TOTAL][ N_FC_FLUX*N_FC_FLUX*N_FC_FLUX ],
                                            const real g_FC_Flux_y[][NCOMP_TOTAL][ N_FC_FLUX*N_FC_FLUX*N_FC_FLUX ],
                                            const real g_FC_Flux_z[][NCOMP_TOTAL][ N_FC_FLUX*N_FC_FLUX*N_FC_FLUX ],
                                            const real dt, const real _dh, const real Gamma,
                                            const real MinDens, const real MinPres )
{

   const uint  bx         = blockIdx.x;
   const uint  tx         = threadIdx.x;
   const uint  dID_FC_Var = blockDim.x;
   const uint3 dID_In     = make_uint3( 1, N_FC_VAR, N_FC_VAR*N_FC_VAR );
   const real  dt_dh2     = (real)0.5*dt*_dh;
   const real  Gamma_m1   = Gamma - (real)1.0;
   const real _Gamma_m1   = (real)1.0/Gamma_m1;

   uint   ID_In_xL, ID_In_yL, ID_In_zL, ID_In_R, ID_FC_Var;
   uint3  ID3d;
   real   FC_Var_xL, FC_Var_xR, FC_Var_yL, FC_Var_yR, FC_Var_zL, FC_Var_zR;
   real   FC_Flux_xL, FC_Flux_xR, FC_Flux_yL, FC_Flux_yR, FC_Flux_zL, FC_Flux_zR;
   real   TGrad1, TGrad2, Corr;
   bool   Inner_x, Inner_y, Inner_z;
   bool   Inner_xy, Inner_yz, Inner_xz;


#  define Load( Input, Output, ID, v )    (  Output = Input[bx][v][ID]  )
#  define Dump( Input, Output, ID, v )    (  Output[bx][v][ID] = Input  )

#  define Correct( FC_Var_L, FC_Var_R, FC_Flux_T1_L, FC_Flux_T1_R, FC_Flux_T2_L, FC_Flux_T2_R )    \
   {                                                                                               \
      TGrad1 = FC_Flux_T1_R - FC_Flux_T1_L;                                                        \
      TGrad2 = FC_Flux_T2_R - FC_Flux_T2_L;                                                        \
      Corr   = -dt_dh2*( TGrad1 + TGrad2 );                                                        \
                                                                                                   \
      FC_Var_L += Corr;                                                                            \
      FC_Var_R += Corr;                                                                            \
   } // Correct

#  define Correct_1v( v )                                                                                   \
   {                                                                                                        \
      /* load the face-centered variables */                                                                \
      if ( Inner_yz )                                                                                       \
      {                                                                                                     \
         Load( g_FC_Var_xL, FC_Var_xL, ID_FC_Var, v );                                                      \
         Load( g_FC_Var_xR, FC_Var_xR, ID_FC_Var, v );                                                      \
      }                                                                                                     \
                                                                                                            \
      if ( Inner_xz )                                                                                       \
      {                                                                                                     \
         Load( g_FC_Var_yL, FC_Var_yL, ID_FC_Var, v );                                                      \
         Load( g_FC_Var_yR, FC_Var_yR, ID_FC_Var, v );                                                      \
      }                                                                                                     \
                                                                                                            \
      if ( Inner_xy )                                                                                       \
      {                                                                                                     \
         Load( g_FC_Var_zL, FC_Var_zL, ID_FC_Var, v );                                                      \
         Load( g_FC_Var_zR, FC_Var_zR, ID_FC_Var, v );                                                      \
      }                                                                                                     \
                                                                                                            \
                                                                                                            \
      /* load the face-ceneterd fluxes */                                                                   \
      if ( Inner_x )                                                                                        \
      {                                                                                                     \
         Load( g_FC_Flux_x, FC_Flux_xL, ID_In_xL, v );                                                      \
         Load( g_FC_Flux_x, FC_Flux_xR, ID_In_R,  v );                                                      \
      }                                                                                                     \
                                                                                                            \
      if ( Inner_y )                                                                                        \
      {                                                                                                     \
         Load( g_FC_Flux_y, FC_Flux_yL, ID_In_yL, v );                                                      \
         Load( g_FC_Flux_y, FC_Flux_yR, ID_In_R,  v );                                                      \
      }                                                                                                     \
                                                                                                            \
      if ( Inner_z )                                                                                        \
      {                                                                                                     \
         Load( g_FC_Flux_z, FC_Flux_zL, ID_In_zL, v );                                                      \
         Load( g_FC_Flux_z, FC_Flux_zR, ID_In_R,  v );                                                      \
      }                                                                                                     \
                                                                                                            \
                                                                                                            \
      /* compute the transverse gradient and correct the face-centered variables */                         \
      if ( Inner_yz )   Correct( FC_Var_xL, FC_Var_xR, FC_Flux_yL, FC_Flux_yR, FC_Flux_zL, FC_Flux_zR );    \
      if ( Inner_xz )   Correct( FC_Var_yL, FC_Var_yR, FC_Flux_xL, FC_Flux_xR, FC_Flux_zL, FC_Flux_zR );    \
      if ( Inner_xy )   Correct( FC_Var_zL, FC_Var_zR, FC_Flux_xL, FC_Flux_xR, FC_Flux_yL, FC_Flux_yR );    \
                                                                                                            \
                                                                                                            \
      /* store the corrected face-centered variables back to the global arrays */                           \
      if ( Inner_yz )                                                                                       \
      {                                                                                                     \
         Dump( FC_Var_xL, g_FC_Var_xL, ID_FC_Var, v );                                                      \
         Dump( FC_Var_xR, g_FC_Var_xR, ID_FC_Var, v );                                                      \
      }                                                                                                     \
                                                                                                            \
      if ( Inner_xz )                                                                                       \
      {                                                                                                     \
         Dump( FC_Var_yL, g_FC_Var_yL, ID_FC_Var, v );                                                      \
         Dump( FC_Var_yR, g_FC_Var_yR, ID_FC_Var, v );                                                      \
      }                                                                                                     \
                                                                                                            \
      if ( Inner_xy )                                                                                       \
      {                                                                                                     \
         Dump( FC_Var_zL, g_FC_Var_zL, ID_FC_Var, v );                                                      \
         Dump( FC_Var_zR, g_FC_Var_zR, ID_FC_Var, v );                                                      \
      }                                                                                                     \
   } // Correct_1v


   ID_FC_Var = tx;

// loop over all cells
   while ( ID_FC_Var < N_FC_VAR*N_FC_VAR*N_FC_VAR )
   {
//    calculate the array indices
      ID_In_R  = ID_FC_Var;
      ID_In_xL = ID_In_R - dID_In.x;
      ID_In_yL = ID_In_R - dID_In.y;
      ID_In_zL = ID_In_R - dID_In.z;

      ID3d.x   = ID_FC_Var%N_FC_VAR;
      ID3d.y   = ID_FC_Var%(N_FC_VAR*N_FC_VAR)/N_FC_VAR;
      ID3d.z   = ID_FC_Var/(N_FC_VAR*N_FC_VAR);

      Inner_x  = ( ID3d.x != 0U  &&  ID3d.x != N_FC_VAR-1 );
      Inner_y  = ( ID3d.y != 0U  &&  ID3d.y != N_FC_VAR-1 );
      Inner_z  = ( ID3d.z != 0U  &&  ID3d.z != N_FC_VAR-1 );
      Inner_xy = Inner_x && Inner_y;
      Inner_yz = Inner_y && Inner_z;
      Inner_xz = Inner_x && Inner_z;


//    apply transverse flux correction
      for (int v=0; v<NCOMP_TOTAL; v++)   Correct_1v( v );


//    ensure positive density and pressure
//    --> make sure that CUFLU_CheckMinPresInEngy() is applied to the **corrected** density
//###OPTIMIZATION: check before calling Dump() to reduce global memory access
      if ( Inner_yz )
      {
         FluVar ConVar;

         ConVar.Rho = g_FC_Var_xL[bx][0][ID_FC_Var];
         ConVar.Px  = g_FC_Var_xL[bx][1][ID_FC_Var];
         ConVar.Py  = g_FC_Var_xL[bx][2][ID_FC_Var];
         ConVar.Pz  = g_FC_Var_xL[bx][3][ID_FC_Var];
         ConVar.Egy = g_FC_Var_xL[bx][4][ID_FC_Var];

         ConVar.Rho                    = FMAX( ConVar.Rho, MinDens );
         g_FC_Var_xL[bx][0][ID_FC_Var] = ConVar.Rho;
         g_FC_Var_xL[bx][4][ID_FC_Var] = CUFLU_CheckMinPresInEngy( ConVar, Gamma_m1, _Gamma_m1, MinPres );


         ConVar.Rho = g_FC_Var_xR[bx][0][ID_FC_Var];
         ConVar.Px  = g_FC_Var_xR[bx][1][ID_FC_Var];
         ConVar.Py  = g_FC_Var_xR[bx][2][ID_FC_Var];
         ConVar.Pz  = g_FC_Var_xR[bx][3][ID_FC_Var];
         ConVar.Egy = g_FC_Var_xR[bx][4][ID_FC_Var];

         ConVar.Rho                    = FMAX( ConVar.Rho, MinDens );
         g_FC_Var_xR[bx][0][ID_FC_Var] = ConVar.Rho;
         g_FC_Var_xR[bx][4][ID_FC_Var] = CUFLU_CheckMinPresInEngy( ConVar, Gamma_m1, _Gamma_m1, MinPres );


//       floor passive scalars
#        if ( NCOMP_PASSIVE > 0 )
         for (int v=NCOMP_FLUID; v<NCOMP_TOTAL; v++)
         {
            g_FC_Var_xL[bx][v][ID_FC_Var] = FMAX( g_FC_Var_xL[bx][v][ID_FC_Var], TINY_NUMBER );
            g_FC_Var_xR[bx][v][ID_FC_Var] = FMAX( g_FC_Var_xR[bx][v][ID_FC_Var], TINY_NUMBER );
         }
#        endif
      } // if ( Inner_yz )

      if ( Inner_xz )
      {
         FluVar ConVar;

         ConVar.Rho = g_FC_Var_yL[bx][0][ID_FC_Var];
         ConVar.Px  = g_FC_Var_yL[bx][1][ID_FC_Var];
         ConVar.Py  = g_FC_Var_yL[bx][2][ID_FC_Var];
         ConVar.Pz  = g_FC_Var_yL[bx][3][ID_FC_Var];
         ConVar.Egy = g_FC_Var_yL[bx][4][ID_FC_Var];

         ConVar.Rho                    = FMAX( ConVar.Rho, MinDens );
         g_FC_Var_yL[bx][0][ID_FC_Var] = ConVar.Rho;
         g_FC_Var_yL[bx][4][ID_FC_Var] = CUFLU_CheckMinPresInEngy( ConVar, Gamma_m1, _Gamma_m1, MinPres );


         ConVar.Rho = g_FC_Var_yR[bx][0][ID_FC_Var];
         ConVar.Px  = g_FC_Var_yR[bx][1][ID_FC_Var];
         ConVar.Py  = g_FC_Var_yR[bx][2][ID_FC_Var];
         ConVar.Pz  = g_FC_Var_yR[bx][3][ID_FC_Var];
         ConVar.Egy = g_FC_Var_yR[bx][4][ID_FC_Var];

         ConVar.Rho                    = FMAX( ConVar.Rho, MinDens );
         g_FC_Var_yR[bx][0][ID_FC_Var] = ConVar.Rho;
         g_FC_Var_yR[bx][4][ID_FC_Var] = CUFLU_CheckMinPresInEngy( ConVar, Gamma_m1, _Gamma_m1, MinPres );


//       floor passive scalars
#        if ( NCOMP_PASSIVE > 0 )
         for (int v=NCOMP_FLUID; v<NCOMP_TOTAL; v++)
         {
            g_FC_Var_yL[bx][v][ID_FC_Var] = FMAX( g_FC_Var_yL[bx][v][ID_FC_Var], TINY_NUMBER );
            g_FC_Var_yR[bx][v][ID_FC_Var] = FMAX( g_FC_Var_yR[bx][v][ID_FC_Var], TINY_NUMBER );
         }
#        endif
      } // if ( Inner_xz )

      if ( Inner_xy )
      {
         FluVar ConVar;

         ConVar.Rho = g_FC_Var_zL[bx][0][ID_FC_Var];
         ConVar.Px  = g_FC_Var_zL[bx][1][ID_FC_Var];
         ConVar.Py  = g_FC_Var_zL[bx][2][ID_FC_Var];
         ConVar.Pz  = g_FC_Var_zL[bx][3][ID_FC_Var];
         ConVar.Egy = g_FC_Var_zL[bx][4][ID_FC_Var];

         ConVar.Rho                    = FMAX( ConVar.Rho, MinDens );
         g_FC_Var_zL[bx][0][ID_FC_Var] = ConVar.Rho;
         g_FC_Var_zL[bx][4][ID_FC_Var] = CUFLU_CheckMinPresInEngy( ConVar, Gamma_m1, _Gamma_m1, MinPres );


         ConVar.Rho = g_FC_Var_zR[bx][0][ID_FC_Var];
         ConVar.Px  = g_FC_Var_zR[bx][1][ID_FC_Var];
         ConVar.Py  = g_FC_Var_zR[bx][2][ID_FC_Var];
         ConVar.Pz  = g_FC_Var_zR[bx][3][ID_FC_Var];
         ConVar.Egy = g_FC_Var_zR[bx][4][ID_FC_Var];

         ConVar.Rho                    = FMAX( ConVar.Rho, MinDens );
         g_FC_Var_zR[bx][0][ID_FC_Var] = ConVar.Rho;
         g_FC_Var_zR[bx][4][ID_FC_Var] = CUFLU_CheckMinPresInEngy( ConVar, Gamma_m1, _Gamma_m1, MinPres );


//       floor passive scalars
#        if ( NCOMP_PASSIVE > 0 )
         for (int v=NCOMP_FLUID; v<NCOMP_TOTAL; v++)
         {
            g_FC_Var_zL[bx][v][ID_FC_Var] = FMAX( g_FC_Var_zL[bx][v][ID_FC_Var], TINY_NUMBER );
            g_FC_Var_zR[bx][v][ID_FC_Var] = FMAX( g_FC_Var_zR[bx][v][ID_FC_Var], TINY_NUMBER );
         }
#        endif
      } // if ( Inner_xy )


      ID_FC_Var += dID_FC_Var;

   } // while ( ID_FC_Var < N_FC_VAR*N_FC_VAR*N_FC_VAR )

#  undef Load
#  undef Dump
#  undef Correct
#  undef Correct_1v

} // FUNCTION : CUFLU_TGradient_Correction



#endif // #if ( defined GPU  &&  MODEL == HYDRO  &&  FLU_SCHEME == CTU )
