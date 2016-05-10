#include "Copyright.h"
#include "Macro.h"
#include "CUFLU.h"

#if (  defined GPU  &&  MODEL == HYDRO  &&  ( FLU_SCHEME == MHM || FLU_SCHEME == MHM_RP )  )



#include "CUFLU_Shared_FluUtility.cu"
#include "CUFLU_Shared_DataReconstruction.cu"
#include "CUFLU_Shared_ComputeFlux.cu"
#include "CUFLU_Shared_FullStepUpdate.cu"
#if   ( RSOLVER == EXACT )
#include "CUFLU_Shared_RiemannSolver_Exact.cu"
#elif ( RSOLVER == ROE )
#include "CUFLU_Shared_RiemannSolver_Roe.cu"
#elif ( RSOLVER == HLLE )
#include "CUFLU_Shared_RiemannSolver_HLLE.cu"
#elif ( RSOLVER == HLLC )
#include "CUFLU_Shared_RiemannSolver_HLLC.cu"
#endif

#if ( FLU_SCHEME == MHM_RP )
static __device__ void CUFLU_RiemannPredict_Flux( const real g_Fluid_In   [][5][ FLU_NXT*FLU_NXT*FLU_NXT ], 
                                                        real g_Half_Flux_x[][5][ N_FC_FLUX*N_FC_FLUX*N_FC_FLUX ],
                                                        real g_Half_Flux_y[][5][ N_FC_FLUX*N_FC_FLUX*N_FC_FLUX ],
                                                        real g_Half_Flux_z[][5][ N_FC_FLUX*N_FC_FLUX*N_FC_FLUX ],
                                                  const real Gamma );
static __device__ void CUFLU_RiemannPredict( const real g_Fluid_In   [][5][ FLU_NXT*FLU_NXT*FLU_NXT ], 
                                             const real g_Half_Flux_x[][5][ N_FC_FLUX*N_FC_FLUX*N_FC_FLUX ],
                                             const real g_Half_Flux_y[][5][ N_FC_FLUX*N_FC_FLUX*N_FC_FLUX ],
                                             const real g_Half_Flux_z[][5][ N_FC_FLUX*N_FC_FLUX*N_FC_FLUX ],
                                                   real g_Half_Var   [][5][ FLU_NXT*FLU_NXT*FLU_NXT ],
                                             const real dt, const real _dh, const real Gamma );
#endif

#ifdef UNSPLIT_GRAVITY
#include "CUPOT.h"
__constant__ double ExtAcc_AuxArray_d_Flu[EXT_ACC_NAUX_MAX];

//-------------------------------------------------------------------------------------------------------
// Function    :  CUFLU_FluidSolver_SetConstMem
// Description :  Set the constant memory used by CUFLU_FluidSolver_CTU
//
// Note        :  Adopt the suggested approach for CUDA version >= 5.0
//
// Parameter   :  None
//
// Return      :  0/-1 : successful/failed
//---------------------------------------------------------------------------------------------------
int CUFLU_FluidSolver_SetConstMem( double ExtAcc_AuxArray_h[] )
{

   if (  cudaSuccess != cudaMemcpyToSymbol( ExtAcc_AuxArray_d_Flu, ExtAcc_AuxArray_h, EXT_ACC_NAUX_MAX*sizeof(double),
                                            0, cudaMemcpyHostToDevice)  )
      return -1;

   else
      return 0;

} // FUNCTION : CUFLU_FluidSolver_SetConstMem
#endif




//-------------------------------------------------------------------------------------------------------
// Function    :  CUFLU_FluidSolver_MHM
// Description :  GPU fluid solver based on the MUSCL-Hancock scheme 
//
// Note        :  1. Prefix "g" for pointers pointing to the "Global" memory space
//                   Prefix "s" for pointers pointing to the "Shared" memory space
//                2. The three-dimensional evolution is achieved by using the unsplit method
//                3. Two half-step prediction schemes are supported, including "MHM" and "MHM_RP"
//                   MHM    : use interpolated face-centered values to calculate the half-step fluxes
//                   MHM_RP : use Riemann solver to calculate the half-step fluxes
//                4. Ref : 
//                   MHM    : "Riemann Solvers and Numerical Methods for Fluid Dynamics 
//                              - A Practical Introduction ~ by Eleuterio F. Toro"
//                   MHM_RP : Stone & Gardiner, NewA, 14, 139 (2009)
//
// Parameter   :  g_Fluid_In     : Global memory array storing the input fluid variables
//                g_Fluid_Out    : Global memory array to store the output fluid variables
//                g_Flux         : Global memory array to store the output fluxes
//                g_Corner       : Global memory array storing the physical corner coordinates of each patch group (USELESS CURRENTLY)
//                g_Pot_USG      : Global memory array storing the input potential for UNSPLIT_GRAVITY (NOT SUPPORTED in RTVD)
//                g_PriVar       : Global memory array to store the primitive variables
//                g_Slope_PPM_x  : Global memory array to store the x-slope for the PPM reconstruction
//                g_Slope_PPM_y  : Global memory array to store the y-slope for the PPM reconstruction
//                g_Slope_PPM_z  : Global memory array to store the z-slope for the PPM reconstruction
//                g_FC_Var_xL    : Global memory array to store the half-step variables on the -x surface
//                g_FC_Var_xR    : Global memory array to store the half-step variables on the +x surface
//                g_FC_Var_yL    : Global memory array to store the half-step variables on the -y surface
//                g_FC_Var_yR    : Global memory array to store the half-step variables on the +y surface
//                g_FC_Var_zL    : Global memory array to store the half-step variables on the -z surface
//                g_FC_Var_zR    : Global memory array to store the half-step variables on the +z surface
//                g_FC_Flux_x    : Global memory array to store the face-centered fluxes in the x direction
//                g_FC_Flux_y    : Global memory array to store the face-centered fluxes in the y direction
//                g_FC_Flux_z    : Global memory array to store the face-centered fluxes in the z direction
//                dt             : Time interval to advance solution
//                _dh            : 1 / grid size
//                Gamma          : Ratio of specific heats
//                StoreFlux      : true --> store the coarse-fine fluxes
//                LR_Limiter     : Slope limiter for the data reconstruction in the MHM/MHM_RP/CTU schemes
//                                 (0/1/2/3/4) = (vanLeer/generalized MinMod/vanAlbada/
//                                                vanLeer + generalized MinMod/extrema-preserving) limiter
//                MinMod_Coeff   : Coefficient of the generalized MinMod limiter
//                EP_Coeff       : Coefficient of the extrema-preserving limiter
//                Time           : Current physical time                                     (for UNSPLIT_GRAVITY only)
//                GravityType    : Types of gravity --> self-gravity, external gravity, both (for UNSPLIT_GRAVITY only)
//-------------------------------------------------------------------------------------------------------
__global__ void CUFLU_FluidSolver_MHM( const real g_Fluid_In[][5][ FLU_NXT*FLU_NXT*FLU_NXT ], 
                                       real g_Fluid_Out  [][5][ PS2*PS2*PS2 ],
                                       real g_Flux    [][9][5][ PS2*PS2 ], 
                                       const double g_Corner[][3],
                                       const real g_Pot_USG[] [ USG_NXT_F*USG_NXT_F*USG_NXT_F ],
                                       real g_PriVar     [][5][ FLU_NXT*FLU_NXT*FLU_NXT ],
                                       real g_Slope_PPM_x[][5][ N_SLOPE_PPM*N_SLOPE_PPM*N_SLOPE_PPM],
                                       real g_Slope_PPM_y[][5][ N_SLOPE_PPM*N_SLOPE_PPM*N_SLOPE_PPM],
                                       real g_Slope_PPM_z[][5][ N_SLOPE_PPM*N_SLOPE_PPM*N_SLOPE_PPM],
                                       real g_FC_Var_xL  [][5][ N_FC_VAR*N_FC_VAR*N_FC_VAR ], 
                                       real g_FC_Var_xR  [][5][ N_FC_VAR*N_FC_VAR*N_FC_VAR ], 
                                       real g_FC_Var_yL  [][5][ N_FC_VAR*N_FC_VAR*N_FC_VAR ], 
                                       real g_FC_Var_yR  [][5][ N_FC_VAR*N_FC_VAR*N_FC_VAR ],
                                       real g_FC_Var_zL  [][5][ N_FC_VAR*N_FC_VAR*N_FC_VAR ], 
                                       real g_FC_Var_zR  [][5][ N_FC_VAR*N_FC_VAR*N_FC_VAR ],
                                       real g_FC_Flux_x  [][5][ N_FC_FLUX*N_FC_FLUX*N_FC_FLUX ],
                                       real g_FC_Flux_y  [][5][ N_FC_FLUX*N_FC_FLUX*N_FC_FLUX ],
                                       real g_FC_Flux_z  [][5][ N_FC_FLUX*N_FC_FLUX*N_FC_FLUX ],
                                       const real dt, const real _dh, const real Gamma, const bool StoreFlux,
                                       const LR_Limiter_t LR_Limiter, const real MinMod_Coeff, 
                                       const real EP_Coeff, const double Time, const OptGravityType_t GravityType )
{

#  ifdef UNSPLIT_GRAVITY
   const bool CorrHalfVel_Yes = true;
#  else
   const bool CorrHalfVel_No  = false;
#  endif


// 1. half-step prediction
#  if   ( FLU_SCHEME == MHM_RP ) // a. use Riemann solver to calculate the half-step fluxes


// use pointers to avoid redundant memory consumption
   real (*const g_Half_Var)   [5][ FLU_NXT*FLU_NXT*FLU_NXT ]       = g_PriVar;
   real (*const g_Half_Flux_x)[5][ N_FC_FLUX*N_FC_FLUX*N_FC_FLUX ] = g_FC_Flux_x;
   real (*const g_Half_Flux_y)[5][ N_FC_FLUX*N_FC_FLUX*N_FC_FLUX ] = g_FC_Flux_y;
   real (*const g_Half_Flux_z)[5][ N_FC_FLUX*N_FC_FLUX*N_FC_FLUX ] = g_FC_Flux_z;


// (1.a-1) evaluate the face-centered half-step fluxes with the piecewise constant data reconstruction
   CUFLU_RiemannPredict_Flux( g_Fluid_In, g_Half_Flux_x, g_Half_Flux_y, g_Half_Flux_z, Gamma );
   __syncthreads();


// (1.a-2) evaluate the half-step cell-centered solution
   CUFLU_RiemannPredict( g_Fluid_In, g_Half_Flux_x, g_Half_Flux_y, g_Half_Flux_z, g_Half_Var, dt, _dh, Gamma );
   __syncthreads();


// (1.a-3) evaluate the half-step face-centered solution by data reconstruction
   CUFLU_DataReconstruction( g_Half_Var, g_Slope_PPM_x, g_Slope_PPM_y, g_Slope_PPM_z, g_FC_Var_xL, g_FC_Var_xR, 
                             g_FC_Var_yL, g_FC_Var_yR, g_FC_Var_zL, g_FC_Var_zR, N_HF_VAR, FLU_GHOST_SIZE-2, 
                             Gamma, LR_Limiter, MinMod_Coeff, EP_Coeff, dt, _dh );
   __syncthreads();


#  elif ( FLU_SCHEME == MHM ) // b. use interpolated face-centered values to calculate the half-step fluxes 


// (1.b-1) conserved variables --> primitive variables
   CUFLU_Con2Pri_AllGrids( g_Fluid_In, g_PriVar, Gamma );
   __syncthreads();


// (1.b-2) evaluate the half-step face-centered solution by data reconstruction
   CUFLU_DataReconstruction( g_PriVar, g_Slope_PPM_x, g_Slope_PPM_y, g_Slope_PPM_z, g_FC_Var_xL, g_FC_Var_xR, 
                             g_FC_Var_yL, g_FC_Var_yR, g_FC_Var_zL, g_FC_Var_zR, FLU_NXT, FLU_GHOST_SIZE-1, 
                             Gamma, LR_Limiter, MinMod_Coeff, EP_Coeff, dt, _dh );
   __syncthreads();


#  endif // #if ( FLU_SCHEME == MHM_RP ) ... else ...


// 2. evaluate the face-centered full-step fluxes by solving the Riemann problem
#  ifdef UNSPLIT_GRAVITY
   CUFLU_ComputeFlux( g_FC_Var_xL, g_FC_Var_xR, g_FC_Var_yL, g_FC_Var_yR, g_FC_Var_zL, g_FC_Var_zR, 
                      g_FC_Flux_x, g_FC_Flux_y, g_FC_Flux_z, g_Flux, StoreFlux, 1, Gamma,
                      CorrHalfVel_Yes, g_Pot_USG, g_Corner, dt, _dh, Time, GravityType, ExtAcc_AuxArray_d_Flu );
#  else
   CUFLU_ComputeFlux( g_FC_Var_xL, g_FC_Var_xR, g_FC_Var_yL, g_FC_Var_yR, g_FC_Var_zL, g_FC_Var_zR, 
                      g_FC_Flux_x, g_FC_Flux_y, g_FC_Flux_z, g_Flux, StoreFlux, 1, Gamma,
                      CorrHalfVel_No, NULL, NULL, NULL_REAL, NULL_REAL, NULL_REAL, GRAVITY_NONE, NULL );
#  endif
   __syncthreads();


// 3. evaluate the full-step solution
   CUFLU_FullStepUpdate( g_Fluid_In, g_Fluid_Out, g_FC_Flux_x, g_FC_Flux_y, g_FC_Flux_z, dt, _dh, Gamma );

} // FUNCTION : CUFLU_FluidSolver_MHM



#if ( FLU_SCHEME == MHM_RP )
//-------------------------------------------------------------------------------------------------------
// Function    :  CUFLU_RiemannPredict_Flux
// Description :  Evaluate the half-step face-centered fluxes by Riemann solver 
//
// Note        :  1. Work for the MUSCL-Hancock method + Riemann-prediction (MHM_RP)
//                2. Currently support the exact, Roe, HLLE, and HLLC solvers
//                3. Prefix "g" for pointers pointing to the "Global" memory space
//                   Prefix "s" for pointers pointing to the "Shared" memory space
//                4. The function is asynchronous
//                   --> "__syncthreads" must be called before using the output data
//                5. The size of the g_Half_Flux_x/y/z arrays are assumed to be "N_HF_FLUX"
//                   --> The fluxes at the left surface of the cell (i+1,j+1,k+1) in "g_Fluid_In" will
//                       be stored at "(k*N_HF_FLUX+j)*N_HF_FLUX+i" in g_Half_Flux_x/y/z
//
// Parameter   :  g_Fluid_In     : Global memory array storing the input fluid variables
//                g_Half_Flux_x  : Global memory array to store the face-centered fluxes in the x direction
//                g_Half_Flux_y  : Global memory array to store the face-centered fluxes in the y direction
//                g_Half_Flux_z  : Global memory array to store the face-centered fluxes in the z direction
//                Gamma          : Ratio of specific heats
//-------------------------------------------------------------------------------------------------------
__device__ void CUFLU_RiemannPredict_Flux( const real g_Fluid_In   [][5][ FLU_NXT*FLU_NXT*FLU_NXT ], 
                                                 real g_Half_Flux_x[][5][ N_FC_FLUX*N_FC_FLUX*N_FC_FLUX ],
                                                 real g_Half_Flux_y[][5][ N_FC_FLUX*N_FC_FLUX*N_FC_FLUX ],
                                                 real g_Half_Flux_z[][5][ N_FC_FLUX*N_FC_FLUX*N_FC_FLUX ],
                                           const real Gamma )
{

   const uint  bx       = blockIdx.x;
   const uint  tx       = threadIdx.x;
   const uint  dID      = blockDim.x;
   const uint3 dID_In   = make_uint3( 1, FLU_NXT, FLU_NXT*FLU_NXT );

#  if ( RSOLVER == EXACT )
   const real  Gamma_m1 = Gamma - (real)1.0;
#  endif

   uint   ID, ID_In, ID_Out, Nxy;
   uint3  ID3d;
   FluVar VarL, VarR, FC_Flux;

#  if ( RSOLVER == EXACT )
   FluVar *Useless = NULL;
#  endif


#  define Load( Input, Output, ID )    \
   {                                   \
      Output.Rho = Input[bx][0][ID];   \
      Output.Px  = Input[bx][1][ID];   \
      Output.Py  = Input[bx][2][ID];   \
      Output.Pz  = Input[bx][3][ID];   \
      Output.Egy = Input[bx][4][ID];   \
   } // Load

#  define Dump( Input, Output, ID )    \
   {                                   \
      Output[bx][0][ID] = Input.Rho;   \
      Output[bx][1][ID] = Input.Px;    \
      Output[bx][2][ID] = Input.Py;    \
      Output[bx][3][ID] = Input.Pz;    \
      Output[bx][4][ID] = Input.Egy;   \
   } // Dump

#  if   ( RSOLVER == EXACT )

      /* exact solver */
      #define RiemannSolver( Dir, VarL, VarR )                                                           \
      {                                                                                                  \
         VarL = CUFLU_Con2Pri( VarL, Gamma_m1 );                                                         \
         VarR = CUFLU_Con2Pri( VarR, Gamma_m1 );                                                         \
                                                                                                         \
         FC_Flux = CUFLU_RiemannSolver_Exact( Dir, *Useless, *Useless, *Useless, VarL, VarR, Gamma );    \
      } // RiemannSolver

#  elif ( RSOLVER == ROE )

      /* Roe solver */
      #define RiemannSolver( Dir, VarL, VarR )                                                           \
      {                                                                                                  \
         FC_Flux = CUFLU_RiemannSolver_Roe( Dir, VarL, VarR, Gamma );                                    \
      } // RiemannSolver

#  elif ( RSOLVER == HLLE )

      /* HLLE solver */
      #define RiemannSolver( Dir, VarL, VarR )                                                           \
      {                                                                                                  \
         FC_Flux = CUFLU_RiemannSolver_HLLE( Dir, VarL, VarR, Gamma );                                   \
      } // RiemannSolver

#  elif ( RSOLVER == HLLC )

      /* HLLC solver */
      #define RiemannSolver( Dir, VarL, VarR )                                                           \
      {                                                                                                  \
         FC_Flux = CUFLU_RiemannSolver_HLLC( Dir, VarL, VarR, Gamma );                                   \
      } // RiemannSolver

#  else

#     error : ERROR : unsupported Riemann solver (EXACT/ROE/HLLE/HLLC) !!

#  endif

#  define GetFlux( Dir, Nx, Ny, Gap_x, Gap_y, Gap_z, dID_In, g_Half_Flux )                                     \
   {                                                                                                           \
      ID  = tx;                                                                                                \
      Nxy = (Nx)*(Ny);                                                                                         \
                                                                                                               \
      /* loop over all cells */                                                                                \
      while ( ID < N_HF_FLUX*(N_HF_FLUX-1)*(N_HF_FLUX-1) )                                                     \
      {                                                                                                        \
         ID3d.x = ID%(Nx);                                                                                     \
         ID3d.y = ID%Nxy/(Nx);                                                                                 \
         ID3d.z = ID/Nxy;                                                                                      \
         ID_In  = __umul24( __umul24( ID3d.z+Gap_z, FLU_NXT   ) + ID3d.y+Gap_y, FLU_NXT    ) + ID3d.x+Gap_x;   \
         ID_Out = __umul24( __umul24( ID3d.z,       N_HF_FLUX ) + ID3d.y,       N_HF_FLUX  ) + ID3d.x;         \
                                                                                                               \
         Load( g_Fluid_In, VarL, ID_In        );                                                               \
         Load( g_Fluid_In, VarR, ID_In+dID_In );                                                               \
                                                                                                               \
         RiemannSolver( Dir, VarL, VarR );                                                                     \
                                                                                                               \
         Dump( FC_Flux, g_Half_Flux, ID_Out );                                                                 \
                                                                                                               \
         ID += dID;                                                                                            \
                                                                                                               \
      } /* while ( ID < N_HF_FLUX*(N_HF_FLUX-1)*(N_HF_FLUX-1) ) */                                             \
   } // GetFlux

   GetFlux( 0, N_HF_FLUX,   N_HF_FLUX-1, 0, 1, 1, dID_In.x, g_Half_Flux_x );
   GetFlux( 1, N_HF_FLUX-1, N_HF_FLUX,   1, 0, 1, dID_In.y, g_Half_Flux_y );
   GetFlux( 2, N_HF_FLUX-1, N_HF_FLUX-1, 1, 1, 0, dID_In.z, g_Half_Flux_z );

#  undef Load
#  undef Dump
#  undef RiemannSolver
#  undef GetFlux

} // FUNCTION : CUFLU_RiemannPredict_Flux



//-------------------------------------------------------------------------------------------------------
// Function    :  CUFLU_RiemannPredict
// Description :  Evolve the cell-centered variables by half time-step by using the Riemann solvers
//
// Note        :  1. Work for the MUSCL-Hancock method + Riemann-prediction (MHM_RP)
//                2. The input array "g_Fluid_In" should be conserved variables
//                3. For the performance consideration, the output data will be primitive variables
//                4. Prefix "g" for pointers pointing to the "Global" memory space
//                   Prefix "s" for pointers pointing to the "Shared" memory space
//                5. The function is asynchronous 
//                   --> "__syncthreads" must be called before using the output data
//                6. The size of the g_Half_Var array are assumed to be "N_HF_VAR"
//
// Parameter   :  g_Fluid_In     : Global memory array storing the input fluid variables
//                g_Half_Flux_x  : Global memory array storing the face-centered fluxes in the x direction
//                g_Half_Flux_y  : Global memory array storing the face-centered fluxes in the y direction
//                g_Half_Flux_z  : Global memory array storing the face-centered fluxes in the z direction
//                g_Half_Var     : Global memory array to store the half-step solution
//                dt             : Time interval to advance solution
//                _dh            : 1 / grid size
//                Gamma          : Ratio of specific heats
//-------------------------------------------------------------------------------------------------------
__device__ void CUFLU_RiemannPredict( const real g_Fluid_In   [][5][ FLU_NXT*FLU_NXT*FLU_NXT ], 
                                      const real g_Half_Flux_x[][5][ N_FC_FLUX*N_FC_FLUX*N_FC_FLUX ],
                                      const real g_Half_Flux_y[][5][ N_FC_FLUX*N_FC_FLUX*N_FC_FLUX ],
                                      const real g_Half_Flux_z[][5][ N_FC_FLUX*N_FC_FLUX*N_FC_FLUX ],
                                            real g_Half_Var   [][5][ FLU_NXT*FLU_NXT*FLU_NXT ],
                                      const real dt, const real _dh, const real Gamma )
{

   const uint  bx       = blockIdx.x;
   const uint  tx       = threadIdx.x;
   const real  dt_dh2   = (real)0.5*dt*_dh;
   const uint  dID_Out  = blockDim.x;
   const uint3 dID_Flux = make_uint3( 1, N_HF_FLUX, N_HF_FLUX*N_HF_FLUX );
   const real  Gamma_m1 = Gamma - (real)1.0;

   uint   ID_Out = tx;
   uint   ID_In, ID_Flux;
   uint3  ID3d;
   FluVar Var;
   real   FluxDiff;

#  if ( defined MIN_PRES_DENS  ||  defined MIN_PRES )
   const real _Gamma_m1 = (real)1.0 / Gamma_m1;
#  endif


// loop over all cells
   while ( ID_Out < N_HF_VAR*N_HF_VAR*N_HF_VAR )
   {
      ID3d.x  = ID_Out%N_HF_VAR;
      ID3d.y  = ID_Out%(N_HF_VAR*N_HF_VAR)/N_HF_VAR;
      ID3d.z  = ID_Out/(N_HF_VAR*N_HF_VAR);
      ID_In   = __umul24(  __umul24( ID3d.z+1, FLU_NXT   ) + ID3d.y+1, FLU_NXT    ) + ID3d.x+1;
      ID_Flux = __umul24(  __umul24( ID3d.z,   N_HF_FLUX ) + ID3d.y,   N_HF_FLUX  ) + ID3d.x;


//    half-step update
      Var.Rho = g_Fluid_In[bx][0][ID_In];
      Var.Px  = g_Fluid_In[bx][1][ID_In];
      Var.Py  = g_Fluid_In[bx][2][ID_In];
      Var.Pz  = g_Fluid_In[bx][3][ID_In];
      Var.Egy = g_Fluid_In[bx][4][ID_In];

#     define Update( comp, v )                                                                                 \
      {                                                                                                        \
         FluxDiff = dt_dh2 * (  g_Half_Flux_x[bx][v][ID_Flux+dID_Flux.x] - g_Half_Flux_x[bx][v][ID_Flux] +     \
                                g_Half_Flux_y[bx][v][ID_Flux+dID_Flux.y] - g_Half_Flux_y[bx][v][ID_Flux] +     \
                                g_Half_Flux_z[bx][v][ID_Flux+dID_Flux.z] - g_Half_Flux_z[bx][v][ID_Flux]  );   \
         Var.comp -= FluxDiff;                                                                                 \
      } // Update

      Update( Rho, 0 );
      Update( Px,  1 );
      Update( Py,  2 );
      Update( Pz,  3 );
      Update( Egy, 4 );

#     undef Update


//    enforce the pressure to be positive
#     if ( defined MIN_PRES_DENS  ||  defined MIN_PRES )
      Var.Egy = CUFLU_PositivePres_In_Engy( Var, Gamma_m1, _Gamma_m1 );
#     endif


//    conserved variables --> primitive variables
      Var = CUFLU_Con2Pri( Var, Gamma_m1 );


//    save the updated data back to the output global array
      g_Half_Var[bx][0][ID_Out] = Var.Rho;
      g_Half_Var[bx][1][ID_Out] = Var.Px;
      g_Half_Var[bx][2][ID_Out] = Var.Py;
      g_Half_Var[bx][3][ID_Out] = Var.Pz;
      g_Half_Var[bx][4][ID_Out] = Var.Egy;


      ID_Out += dID_Out;

   } // while ( ID_Out < N_HF_VAR*N_HF_VAR*N_HF_VAR )

} // FUNCTION : CUFLU_RiemannPredict
#endif // #if ( FLU_SCHEME == MHM_RP )



#endif // #if (  defined GPU  &&  MODEL == HYDRO  &&  ( FLU_SCHEME == MHM || FLU_SCHEME == MHM_RP )  )
