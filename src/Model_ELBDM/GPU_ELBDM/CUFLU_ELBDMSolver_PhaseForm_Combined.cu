#include "Macro.h"
#include "CUFLU.h"

#if ( MODEL == ELBDM && ELBDM_SCHEME == HYBRID )

#ifndef __umul24
#  define __umul24( a, b )   ( (a)*(b) )
#endif 
#ifndef __mul24
#  define  __mul24( a, b )   ( (a)*(b) )
#endif

// useful macros
# define to1D1(z,y,x) ( __umul24(z, FLU_NXT*FLU_NXT) + __umul24(y, FLU_NXT) + x )
# define to1D2(z,y,x) ( __umul24(z-FLU_GHOST_SIZE, PS2*PS2) + __umul24(y-FLU_GHOST_SIZE, PS2) + x-FLU_GHOST_SIZE )

#ifdef __CUDACC__
# define CGPU_FLU_BLOCK_SIZE_X FLU_BLOCK_SIZE_X
# define CGPU_FLU_BLOCK_SIZE_Y FLU_BLOCK_SIZE_Y
#else 
# define CGPU_FLU_BLOCK_SIZE_X 1
# define CGPU_FLU_BLOCK_SIZE_Y 1
#endif 

GPU_DEVICE
static uint get1D1(uint k, uint j, uint i, int XYZ) {
   switch ( XYZ )
   {
      case 0:  return to1D1( k, j, i );
      case 3:  return to1D1( k, i, j );
      case 6:  return to1D1( i, k, j );
   }
   return 0;
}

GPU_DEVICE
static uint get1D2(uint k, uint j, uint i, int XYZ) {
   switch ( XYZ )
   {
      case 0:  return to1D2( k, j, i );
      case 3:  return to1D2( k, i, j );
      case 6:  return to1D2( i, k, j );
   }
   return 0;
}

# define CELL_LOOP( NCell, leftGhost, rightGhost )    for ( (Idx   = tid, (NStep = (NCell) - (leftGhost) - (rightGhost), (si = Idx % NStep + (leftGhost), sj = Idx / NStep))); \
                                                             Idx   < NColumnOnce * NStep; \
                                                            (Idx  += NThread, (si = Idx % NStep + (leftGhost), sj = Idx / NStep)) )

# define GTR( a, b )     (  ( (a) > (b) ) ? (1) : (0)  )
# define LSS( a, b )     (  ( (a) < (b) ) ? (1) : (0)  )
# define SGN( a )        (  ( (a) > (0) ) ? (1) : ( (a) < (0) ) ? (-1) : (0) )

# ifdef SMOOTH_PHASE
# define  TWOPI          (  real(2 * M_PI)  )
# define _TWOPI          (  real(1.0) / TWOPI )
# define UNWRAP( ref, wrap ) (  round(((ref) - (wrap)) * _TWOPI) )
# define GRADIENT_RATIO(f, t) ( (f[t+1] - f[t  ]) / (f[t] - f[t-1] + (((f[t] - f[t-1]) == 0) ? 1e-8 : 0)))
# endif // # ifdef SMOOTH_PHASE

// First-order forward and backward gradient
# define GRADF1(In, t) ( In[t + 1] - In[t    ] )
# define GRADB1(In, t) ( In[t    ] - In[t - 1] )

// Second-order forward and backward gradient
# define GRADF2(In, t) ( real(1.0/2.0) * ( -3*In[t] + 4*In[t+1] - In[t+2]))
# define GRADB2(In, t) ( real(1.0/2.0) * (  3*In[t] - 4*In[t-1] + In[t-2]))

// Third-order forward and backward gradient
# define GRADF3(In, t) ( real(1.0/6.0) * ( -2*In[t-1] - 3*In[t] + 6*In[t+1] - In[t+2]))
# define GRADB3(In, t) ( real(1.0/6.0) * (  2*In[t+1] + 3*In[t] - 6*In[t-1] + In[t-2]))

// Second- and fourth-order laplacian
# define LAP2(In, t)   ( real(1.0/1.0 ) * (            + 1 *In[t-1] - 2 *In[t] + 1 *In[t+1]             ))
# define LAP4(In, t)   ( real(1.0/12.0) * ( -1*In[t-2] + 16*In[t-1] - 30*In[t] + 16*In[t+1] - 1*In[t+2] ))

// Second- and fourth-order centered gradient
# define GRADC2(In, t) ( real(1.0/2.0 ) * (            - 1*In[t-1] + 1*In[t+1]             ))
# define GRADC4(In, t) ( real(1.0/12.0) * (  1*In[t-2] - 8*In[t-1] + 8*In[t+1] - 1*In[t+2] ))

//First-order upwind flux reconstruction
# define UPWIND_FM(Rc, Vb, t)        (   FMAX(Vb, 0) * Rc[t-1] \
                                       + FMIN(Vb, 0) * Rc[t  ] )

//Second-order MUSCL flux reconstruction
# if ( HYBRID_SCHEME == HYBRID_MUSCL )

// VAN ALBADA LIMITER
# define LIMITER(In) ( (SQR(In) + In)/((real)1. + SQR(In)) )
// MC
//# define LIMITER(r) MAX(0, MIN(MIN((1. + r) / 2., 2.), 2. * r))
// VAN LEER
//# define LIMITER(r) ((r + FABS(r))/(1.0 + FABS(r)))
// SUPERBEE
//# define LIMITER(r) MAX(MAX(0, MIN(2*r, 1)), MIN(r, 2))

# define UPWIND_GRADIENT_RATIO(Rc, Vb, t) ( ((Rc[t-1] - Rc[t-2]) * GTR(Vb, 0)  \
                                           + (Rc[t+1] - Rc[t  ]) * LSS(Vb, 0)) \
                                           / (Rc[t] - Rc[t-1] + (((Rc[t] - Rc[t-1]) == 0) ? 1e-8 : 0)))

# define MUSCL_FM(Rc, Vb, t, _v) (   FMAX(Vb, 0) * Rc[t-1] \
                                   + FMIN(Vb, 0) * Rc[t  ] \
                                   + real(0.5) * FABS(Vb) * (1. - FABS(Vb * _v)) * LIMITER(UPWIND_GRADIENT_RATIO(Rc, Vb, t)) * (Rc[t] - Rc[t - 1]) )


//Second-order unlimited flux reconstruction
# elif  ( HYBRID_SCHEME == HYBRID_FROMM )

# define FROMM_FM(Rc, Vb, t, _v) (   FMAX(Vb, 0) * ( Rc[t-1] + real(0.25) * ( Rc[t  ] - Rc[t-2] ) * (1.0 - Vb * _v)) \
                                   + FMIN(Vb, 0) * ( Rc[t  ] - real(0.25) * ( Rc[t+1] - Rc[t-1] ) * (1.0 + Vb * _v)) )


//Third-order monotonic flux reconstruction
#elif ( HYBRID_SCHEME == HYBRID_PPM )
void PPM_INTERPOLATION(real* a_array, real* a_L_array, real* a_R_array, int i, int densityLimiter);
void PPM_LIMITER      (real* a_array, real* a_L_array, real* a_R_array, int i, int densityLimiter);
real PPM_FM           (real* a_array, real* a_L_array, real* a_R_array, real* v_L_array, int i, real dh, real dt);
# endif // # if ( HYBRID_SCHEME == HYBRID_PPM )

//Second- and fourth-order quantum pressure depending on the scheme used
# if (HYBRID_SCHEME == HYBRID_FROMM || HYBRID_SCHEME == HYBRID_MUSCL || HYBRID_SCHEME == HYBRID_UPWIND)
# define QUANTUM_PRESSURE(Sr, t)  ((real) 0.5 * (pow(GRADC2(Sr, t), 2) + LAP2(Sr, t)))
# else 
# define QUANTUM_PRESSURE(Sr, t)  ((real) 0.5 * (pow(GRADC4(Sr, t), 2) + LAP4(Sr, t)))
# endif // # if (HYBRID_SCHEME == HYBRID_FROMM || HYBRID_SCHEME == HYBRID_MUSCL || HYBRID_SCHEME == HYBRID_UPWIND)


#  if ( HYBRID_SCHEME == HYBRID_FROMM || HYBRID_SCHEME == HYBRID_MUSCL || HYBRID_SCHEME == HYBRID_UPWIND)
# define GHOST_ZONE_PER_STAGE     2                      // ghost zone per Runge-Kutta stage                       
#  else 
# define GHOST_ZONE_PER_STAGE     4                      // ghost zone per Runge-Kutta stage
#  endif 


//Osher-Sethian flux for Hamilton-Jacobi equation
# define OSHER_SETHIAN_FLUX(vp, vm) ((real) 0.5 * (pow(MIN(vp, 0), 2) + pow(MAX(vm, 0), 2)))


//#define N_TIME_LEVELS 1
//const real TIME_COEFFS[N_TIME_LEVELS]                = {1.0};
//const real RK_COEFFS  [N_TIME_LEVELS][N_TIME_LEVELS] = {{1.0}};
//const static real FLUX_COEFFS[N_TIME_LEVELS]         = {1.0};

//#define N_TIME_LEVELS 2
//GPU_DEVICE_VARIABLE
//const static real TIME_COEFFS[N_TIME_LEVELS]                = {1.0/2.0, 1.0};
//
//GPU_DEVICE_VARIABLE
//const static real RK_COEFFS  [N_TIME_LEVELS][N_TIME_LEVELS] = {{1.0, 0.0}, {1.0, 0.0}};
//
//#ifdef CONSERVE_MASS
//GPU_DEVICE_VARIABLE
//const static real FLUX_COEFFS[N_TIME_LEVELS]                = {0.0, 1.0};
//#endif 

#define N_TIME_LEVELS 3
GPU_DEVICE_VARIABLE
const static real TIME_COEFFS[N_TIME_LEVELS]                = {1.0, 1.0/4.0, 2.0/3.0};

GPU_DEVICE_VARIABLE
const static real RK_COEFFS  [N_TIME_LEVELS][N_TIME_LEVELS] = {{1.0, 0.0, 0.0}, {3.0/4.0, 1.0/4.0, 0.0}, {1.0/3.0, 0.0, 2.0/3.0}};

#ifdef CONSERVE_MASS
GPU_DEVICE_VARIABLE
const static real FLUX_COEFFS[N_TIME_LEVELS]                = {1.0/6.0, 1.0/6.0, 2.0/3.0};
#endif 

#define NO_LIMITER    0
#define SOME_LIMITER  1
#define FULL_LIMITER  2 

GPU_DEVICE
static void CUFLU_Advance( real g_Fluid_In [][FLU_NIN ][ CUBE(FLU_NXT) ],
                           #ifdef GAMER_DEBUG
                           real g_Fluid_Out[][FLU_NOUT][ CUBE(PS2) ],
                           #else
                           real g_Fluid_Out[][FLU_NIN ][ CUBE(PS2) ],
                           #endif 
                           real g_Flux     [][9][NFLUX_TOTAL][ SQR(PS2) ],
                           int NPatchGroup,
                           const real dt, const real _dh, const real Eta, const bool StoreFlux,
                           const uint j_gap, const uint k_gap,         
                           real s_In    [][N_TIME_LEVELS+1][FLU_NIN][FLU_NXT],
                           real s_LogRho[][FLU_NXT],
                           real s_Fm    [][2][FLU_NXT],
                           real s_Flux  [][FLU_NXT], 
                           bool s_RK1  [][FLU_NXT],
                           int  s_2PI   [][FLU_NXT],
                           const bool FinalOut, const int XYZ, const real MinDens );

//-------------------------------------------------------------------------------------------------------
// Function    :  CUFLU_ELBDMSolver_PhaseForm
// Description :  GPU solver for kinetic term in Hamilton Jacobi-Madelung equations
//
// Note        :  1. The three-dimensional evolution is achieved by applying x, y, and z operators successively.
//                   Since these operators commute, the order of applying them are irrelevant.
//                   --> Input parameter "XYZ" is actually useless#
//                2. Currently supports 3 modes: HYBRID_SCHEME == FIRST_ORDER, SECOND_ORDER, THIRD_ORDER set by global compile-time constants 
//                   They respectively use first-order upwinding, second-order PLM flux reconstructions and third-order flux reconstructions with higher-order finite differences
//                2. Prefix "g" for pointers pointing to the "Global" memory space
//                   Prefix "s" for pointers pointing to the "Shared" memory space
//
// Parameter   :  g_Fluid_In     : Global memory array storing the input variables
//                g_Fluid_Out    : Global memory array to store the output variables
//                g_Flux         : Global memory array to store the output fluxes (useful only if StoreFlux == true)
//                dt             : Time interval to advance solution
//                _dh            : 1 / grid size
//                Eta            : Particle mass / Planck constant
//                StoreFlux      : true --> store the coarse-fine fluxes
//                                      --> useful only if CONSERVE_MASS is defined
//                XYZ            : true  : x->y->z ( forward sweep)
//                                 false : z->y->x (backward sweep)
//                                 --> Meaningless if CONSERVE_MASS is off since the operators along different directions
//                                     commute
//                                 --> Meaningful if CONSERVE_MASS is on, in which the symmetry along different directions
//                                     are broken ...
//                MinDens        : Minimum allowed density
//-------------------------------------------------------------------------------------------------------



#ifdef __CUDACC__
__global__
void CUFLU_ELBDMSolver_PhaseForm( real g_Fluid_In [][FLU_NIN ][ CUBE(FLU_NXT) ],
                                  #ifdef GAMER_DEBUG
                                  real g_Fluid_Out[][FLU_NOUT ][ CUBE(PS2) ],
                                  #else
                                  real g_Fluid_Out[][FLU_NIN ][ CUBE(PS2) ],
                                  #endif 
                                  real g_Flux     [][9][NFLUX_TOTAL][ SQR(PS2) ],
                                  const real dt, const real _dh, const real Eta, const bool StoreFlux,
                                  const bool XYZ, const real MinDens )
#else
void CPU_ELBDMSolver_PhaseForm(   real g_Fluid_In [][FLU_NIN ][ CUBE(FLU_NXT)], 
                                  #ifdef GAMER_DEBUG
                                  real g_Fluid_Out[][FLU_NOUT ][ CUBE(PS2) ],
                                  #else
                                  real g_Fluid_Out[][FLU_NIN ][ CUBE(PS2) ],
                                  #endif 
                                  real g_Flux     [][9][NFLUX_TOTAL][ SQR(PS2) ],
                                  const int NPatchGroup, 
                                  const real dt, const real dh, const real Eta, const bool StoreFlux,
                                  const bool XYZ, const real MinDens )
#endif 
{


#  ifdef __CUDACC__
// create memories for columns of various intermediate fields in shared GPU memory
   __shared__ real s_In      [CGPU_FLU_BLOCK_SIZE_Y][N_TIME_LEVELS + 1][FLU_NIN][FLU_NXT];
   __shared__ real s_LogRho  [CGPU_FLU_BLOCK_SIZE_Y][FLU_NXT]; // one column of the 0.5 * log(rho) for every thread block
   __shared__ real s_Fm      [CGPU_FLU_BLOCK_SIZE_Y][2][FLU_NXT]; // one column of the fluxes for every thread block
   __shared__ bool s_RK1    [CGPU_FLU_BLOCK_SIZE_Y][FLU_NXT]; 

#  ifdef CONSERVE_MASS
   __shared__ real s_Flux    [CGPU_FLU_BLOCK_SIZE_Y][FLU_NXT];
#  else  // #  ifdef CONSERVE_MASS
              real (*s_Flux)[FLU_NXT] = NULL;  // useless if CONSERVE_MASS is off
#  endif // #  ifdef CONSERVE_MASS ... # else

#  ifdef SMOOTH_PHASE
   __shared__ int   s_2PI    [CGPU_FLU_BLOCK_SIZE_Y][FLU_NXT];
#  else // # ifdef SMOOTH_PHASE
              int (*s_2PI)   [FLU_NXT] = NULL; 
#  endif // # ifdef SMOOTH_PHASE


   const int NPatchGroup = 0; 

#  else // #  ifdef __CUDACC__
// allocate memory on stack within loop for CPU run
   real (*s_In)     [N_TIME_LEVELS + 1][FLU_NIN][FLU_NXT]      = NULL;
   real (*s_LogRho) [FLU_NXT]                              = NULL;
   real (*s_Fm)     [2][FLU_NXT]                           = NULL;
   bool (*s_RK1)   [FLU_NXT]                              = NULL;
   real (*s_Flux)   [FLU_NXT]                              = NULL;
   int  (*s_2PI)    [FLU_NXT]                              = NULL; 
   const real _dh                                          = real(1.0)/dh; 
#  endif

   if ( XYZ )
   {
      CUFLU_Advance( g_Fluid_In, g_Fluid_Out, g_Flux, NPatchGroup, dt, _dh, Eta, StoreFlux,
                                  0,              0, s_In, s_LogRho, s_Fm, s_Flux, s_RK1, s_2PI, false, 0, MinDens );
      CUFLU_Advance( g_Fluid_In, g_Fluid_Out, g_Flux, NPatchGroup, dt, _dh, Eta, StoreFlux,
                     FLU_GHOST_SIZE,              0, s_In, s_LogRho, s_Fm, s_Flux, s_RK1, s_2PI, false, 3, MinDens );
      CUFLU_Advance( g_Fluid_In, g_Fluid_Out, g_Flux, NPatchGroup, dt, _dh, Eta, StoreFlux,
                     FLU_GHOST_SIZE, FLU_GHOST_SIZE, s_In, s_LogRho, s_Fm, s_Flux, s_RK1, s_2PI,  true, 6, MinDens );
   }

   else
   {
      CUFLU_Advance( g_Fluid_In, g_Fluid_Out, g_Flux, NPatchGroup, dt, _dh, Eta, StoreFlux,
                                  0,              0, s_In, s_LogRho, s_Fm, s_Flux, s_RK1, s_2PI, false, 6, MinDens );
      CUFLU_Advance( g_Fluid_In, g_Fluid_Out, g_Flux, NPatchGroup, dt, _dh, Eta, StoreFlux,
                                  0, FLU_GHOST_SIZE, s_In, s_LogRho, s_Fm, s_Flux, s_RK1, s_2PI, false, 3, MinDens );
      CUFLU_Advance( g_Fluid_In, g_Fluid_Out, g_Flux, NPatchGroup, dt, _dh, Eta, StoreFlux,
                     FLU_GHOST_SIZE, FLU_GHOST_SIZE, s_In, s_LogRho, s_Fm, s_Flux, s_RK1, s_2PI,  true, 0, MinDens );
   }

} // FUNCTION : CUFLU_ELBDMSolver_PhaseForm



//-------------------------------------------------------------------------------------------------------
// Function    :  CUFLU_Advance
// Description :  Use CPU/GPU to advance solutions by one time-step
//
// Note        :  1. Based on solving the continuity equation via upwinding and higher-order flux reconstruction and the Hamilton-Jacobi equation via upwinding with higher-order central differences
//                2. Prefix "g" for pointers pointing to the "Global" memory space
//                   Prefix "s" for pointers pointing to the "Shared" memory space
//                3. The direction of the one dimensional sweep is determined by the input parameter "XYZ"
//
// Parameter   :  g_Fluid_In     : Global memory array storing the input variables
//                g_Fluid_Out    : Global memory array to store the output variables
//                g_Flux         : Global memory array to store the output fluxes (useful only if StoreFlux == true)
//                dt             : Time interval to advance solution
//                _dh            : 1 / grid size
//                Eta            : Particle mass / Planck constant
//                StoreFlux      : true --> store the coarse-fine fluxes
//                                   --> useful only if CONSERVE_MASS is defined
//                j_gap          : Number of useless grids on each side in the j direction (j may not be equal to y)
//                k_gap          : Number of useless grids on each side in the k direction (k mya not be equal to z)
//                s_In           : Shared memory array to store the input data and the solutions at different times
//                s_LogRho           : Shared memory array to store the density logarithms
//                s_Fm           : Shared memory array to store the boundary density fluxes during the computation
//                s_Flux         : Shared memory array to store the boundary fluxes
//                FinalOut       : true --> store the updated data to g_Fluid_Out
//                XYZ            : 0 : Update the solution in the x direction
//                                 3 : Update the solution in the y direction
//                                 6 : Update the solution in the z direction
//                                 --> This parameter is also used to determine the place to store the output fluxes
//                MinDens        : Minimum allowed density
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE
void CUFLU_Advance(  real g_Fluid_In [][FLU_NIN ][ CUBE(FLU_NXT) ],
                     #ifdef GAMER_DEBUG
                     real g_Fluid_Out[][FLU_NOUT ][ CUBE(PS2) ],
                     #else
                     real g_Fluid_Out[][FLU_NIN ][ CUBE(PS2) ],
                     #endif 
                     real g_Flux     [][9][NFLUX_TOTAL][ SQR(PS2) ],
                     int NPatchGroup,
                     const real dt, const real _dh, const real Eta, const bool StoreFlux,
                     const uint j_gap, const uint k_gap, 
                     real s_In     [][N_TIME_LEVELS + 1][FLU_NIN][FLU_NXT],
                     real s_LogRho [][FLU_NXT],
                     real s_Fm     [][2][FLU_NXT],
                     real s_Flux   [][FLU_NXT], 
                     bool s_RK1    [][FLU_NXT],
                     int  s_2PI    [][FLU_NXT],
                     const bool FinalOut,
                     const int XYZ, const real MinDens )
{

   const real dh           = real(1.0)/_dh;                      // grid spacing
   const real Coeff1       = real(1.0) * dt /(dh * Eta);           // coefficient for continuity equation
   const real Coeff2       = real(0.5) * dt /(dh * dh * Eta);      // coefficient for HJ-equation
   const real Coeff3       = dh * real(0.5) * Eta * real(3.0);     // coefficient for determining velocity timestep
   const real FluidMinDens = FMAX(real(1e-10), MinDens);          // minimum density while computing quantum pressure and when correcting negative density 

   const uint j_end        = FLU_NXT -  j_gap    ;          // last y-column to be updated

   const uint size_j       = FLU_NXT - (j_gap<<1);          // number of y-columns to be updated
   const uint size_k       = FLU_NXT - (k_gap<<1);          // number of z-columns to be updated
   const uint NColumnTotal = __umul24( size_j, size_k );    // total number of data columns to be updated

// openmp pragma for the CPU solver
#  ifndef __CUDACC__
#  pragma omp parallel
#  endif
   {
#     ifdef __CUDACC__
      const int bx = blockIdx.x;
#     else
//    in CPU mode, every thread works on one patch group at a time and corresponds to one block in the grid of the GPU solver
#     pragma omp for schedule( runtime ) private ( s_In, s_LogRho, s_Fm, s_Flux, s_RK1, s_2PI )
      for (int bx=0; bx<NPatchGroup; bx++)
#     endif
      {

         uint Column0 = 0;                // the total number of columns that have been updated
         uint Idx, Idx1, Idx2;            // temporary indices used for indexing column updates, writing data to g_Fluid_In, g_Fluid_Out
#        ifdef CONSERVE_MASS
         uint Idx3;                       // temporary index used for writing data to g_Flux
#        endif
#        ifdef SMOOTH_PHASE
         uint Idx4;
#        endif // # ifdef SMOOTH_PHASE
         uint delta_k;
         uint si, sj;                     // array indices used in the shared memory array
         uint g1, g2;                     // left and right ghost zones while updating each column

         uint NStep;                      // number of iterations for updating each column
         real De_New, Ph_New, v, qp, vp, vm;
         int l, l_min, l_max;
         uint time_level; 

#        ifdef GAMER_DEBUG
         real St_New;
#        endif 

#        ifdef __CUDACC__
//       use two-dimensional thread blocks in GPU mode
         const uint tx            = threadIdx.x;
         const uint ty            = threadIdx.y;
#        else  // # ifdef __CUDACC__
//       every block just has a single thread with temporary memory on the stack in CPU mode
         const uint tx            = 0;
         const uint ty            = 0;

//       create arrays for columns of various intermediate fields on the stack
         real s_In_1PG       [CGPU_FLU_BLOCK_SIZE_Y][N_TIME_LEVELS + 1][FLU_NIN][FLU_NXT];   // density and phase fields at all RK stages + RK 1 result
         real s_LogRho_1PG   [CGPU_FLU_BLOCK_SIZE_Y][FLU_NXT];                               // 1/2 * density logarithm 
         real s_Fm_1PG       [CGPU_FLU_BLOCK_SIZE_Y][2][FLU_NXT];                               // density flux

         bool s_RK1_1PG     [CGPU_FLU_BLOCK_SIZE_Y][FLU_NXT];

#        ifdef CONSERVE_MASS
         real s_Flux_1PG     [CGPU_FLU_BLOCK_SIZE_Y][FLU_NXT];
#        endif // #  ifdef CONSERVE_MASS

#        ifdef SMOOTH_PHASE
         int  s_2PI_1PG      [CGPU_FLU_BLOCK_SIZE_Y][FLU_NXT]; 
#        endif // # ifdef SMOOTH_PHASE

         s_In       = s_In_1PG;
         s_LogRho   = s_LogRho_1PG; 
         s_Fm       = s_Fm_1PG;
         s_RK1     = s_RK1_1PG;

#        ifdef CONSERVE_MASS
         s_Flux     = s_Flux_1PG;
#        endif // #  ifdef CONSERVE_MASS

#        ifdef SMOOTH_PHASE
         s_2PI      = s_2PI_1PG; 
#        endif 

#        endif // # ifdef __CUDACC__ ... # else 

         const uint tid          = __umul24( ty , CGPU_FLU_BLOCK_SIZE_X ) + tx;    // thread ID within block
               uint j            = j_gap + ty % size_j;                    // (i,j,k): array indices used in g_Fluid_In
               uint k            = k_gap + ty / size_j;                    // (i,j,k): array indices used in g_Fluid_In
               uint i            = tx + FLU_GHOST_SIZE;                    // (i,j,k): array indices used in g_Fluid_In
         uint NColumnOnce        = MIN( NColumnTotal, CGPU_FLU_BLOCK_SIZE_Y );     // number of columns updated per iteration
         const uint NThread      = CGPU_FLU_BLOCK_SIZE_X * CGPU_FLU_BLOCK_SIZE_Y;          // total number of threads within block


//       determine the array indices for loading the ghost-zone data for GPU solver
#        ifdef __CUDACC__
         bool LoadGhost = false;                // true --> load the ghost-zone data
         uint LoadGhost_i;
         int  LoadGhost_di, LoadGhost_dIdx1;


//       use the first 2*FLU_GHOST_SIZE threads to load the ghost zones in addition to their regular cells 
         if ( tx < 2*FLU_GHOST_SIZE )
         {
            LoadGhost = true;

            if ( tx < FLU_GHOST_SIZE )    LoadGhost_di = -FLU_GHOST_SIZE;
            else                          LoadGhost_di = -FLU_GHOST_SIZE + PS2;

            switch ( XYZ )
            {
               case 0:  LoadGhost_dIdx1 = LoadGhost_di;                                break;
               case 3:  LoadGhost_dIdx1 = __mul24( LoadGhost_di, FLU_NXT );            break;
               case 6:  LoadGhost_dIdx1 = __mul24( LoadGhost_di, FLU_NXT*FLU_NXT );    break;
            }

            LoadGhost_i = (int)i + LoadGhost_di;
         } // if ( tx < 2*FLU_GHOST_SIZE )
#        endif // #ifdef __CUDACC__

//       loop over all data columns
         while ( Column0 < NColumnTotal )
         {

//          1. load data into shared memory
            if ( tid < NColumnOnce*PS2 )
            {

               time_level = 0; 
            
#              ifdef __CUDACC__ 
//             1.1 determine the array indices for loading global memory data along different directions
               Idx1 = get1D1( k, j, i, XYZ );

//             1.2 load the interior data into shared memory at time_level 0
               s_In[ty][time_level][DENS][i] = g_Fluid_In[bx][DENS][Idx1];
               s_In[ty][time_level][PHAS][i] = g_Fluid_In[bx][PHAS][Idx1];
#              ifdef CONSERVE_MASS
               s_Flux[ty][i] = 0;
#              endif
               s_RK1[ty][i] = false;

//             1.3 load the ghost-zone data into shared memory
               if ( LoadGhost )
               {
                  s_In[ty][time_level][DENS][LoadGhost_i] = g_Fluid_In[bx][DENS][ (int)Idx1 + LoadGhost_dIdx1 ];
                  s_In[ty][time_level][PHAS][LoadGhost_i] = g_Fluid_In[bx][PHAS][ (int)Idx1 + LoadGhost_dIdx1 ];
#                 ifdef CONSERVE_MASS
                  s_Flux[ty][LoadGhost_i] = 0;
#                 endif 
                  s_RK1[ty][LoadGhost_i] = false;
               }
#              else // # ifdef __CUDACC__

               for (si = 0; si < FLU_NXT; ++si)
               {
                  Idx1 = get1D1( k, j, si, XYZ );
                  s_In[ty][time_level][DENS][si] = g_Fluid_In[bx][DENS][Idx1];
                  s_In[ty][time_level][PHAS][si] = g_Fluid_In[bx][PHAS][Idx1];
                  s_RK1[ty][si] = false;
#                 ifdef CONSERVE_MASS
                  s_Flux[ty][si] = 0; 
#                 endif         
               }
#              endif // # ifdef __CUDACC__  ... else
            } // if ( tid < NColumnOnce*PS2 )

//          1.4 sync data read into S_In
#           ifdef __CUDACC__
            __syncthreads();
#           endif // # ifdef __CUDACC_

//          1.5 smooth out 2 PI-discontinuitites in phase field
#           ifdef SMOOTH_PHASE



            bool hasChanged = false; 
            
//          1.5.1. compute magnitude of 2 pi-jumps wherever we detect discontinuity
            CELL_LOOP(FLU_NXT, 1, 1)
            {
               if ( GRADIENT_RATIO(s_In[sj][0][PHAS], si) < - real(0.0) ) {
                  s_2PI[sj][si] = UNWRAP(s_In[sj][0][PHAS][si - 1], s_In[sj][0][PHAS][si]);
                  if (si == FLU_NXT - 2) {
                     s_2PI[sj][FLU_NXT - 1] = UNWRAP(s_In[sj][0][PHAS][si], s_In[sj][0][PHAS][si + 1]);
                  }
               } else {
                  s_2PI[sj][si] = 0;
                  if (si == FLU_NXT - 2) {
                     s_2PI[sj][si + 1] = 0;
                  }
               }
            }

#           ifdef __CUDACC__
            __syncthreads();
#           endif // # ifdef __CUDACC_

//          1.5.2. add multiples of 2 pi to initial phase field to smoothen it
            CELL_LOOP(FLU_NXT, 1, 0)
            {
               for (Idx4 = 1; Idx4 <= si; ++Idx4) {
                  s_In[sj][time_level][PHAS][si] += s_2PI[sj][Idx4] * TWOPI;
               }
            }

#           ifdef __CUDACC__
            __syncthreads();
#           endif // # ifdef __CUDACC_
#           endif // # ifdef SMOOTH_PHASE

//          2. Runge-Kutta iterations
            for (time_level = 0; time_level < N_TIME_LEVELS; ++time_level) 
            {
               g1 = GHOST_ZONE_PER_STAGE *   time_level       ;
               g2 = GHOST_ZONE_PER_STAGE * ( time_level + 1 ) ;         

//             2.2 compute density logarithms
               CELL_LOOP(FLU_NXT, g1 + 1, g1 + 1)
               { 
                  s_LogRho[sj][si] = log(FMAX(s_In[sj][time_level][DENS][si], FluidMinDens));
               } 

//             2.3 sync s_Fm, s_LogRho and s_Flux
#              ifdef __CUDACC__ 
               __syncthreads();
#              endif 

//             2.1 check the velocity-dependent CFL-condition and switch to forward-Euler for updating the density wherever the CFL-condition is not met
               CELL_LOOP(FLU_NXT, g2, g2 - 1)
               {
                  v = _dh * GRADC2 (s_In[sj][time_level][PHAS], si);
//                dt = 1 / MaxdS_dx * 0.5 * ELBDM_ETA * DT__VELOCITY;
//                compute CFL condition timestep and quantum pressure term
                  qp     = real(1.0/2.0) * LAP2(s_LogRho[sj], si)  + real(1.0/4.0) * SQR(GRADC2(s_LogRho[sj], si));
                  
//                if the time step adopted in solver is larger than what velocity-dependent CFL condition allows, we switch to RK1 with a first-order upwind discretisation
                  if ( dt * FABS(v)  > Coeff3 || FABS(qp) > 0.15) {
//                   compute how far wrong information can propagate
                     l_min = si - (N_TIME_LEVELS - time_level) * 1;
                     l_max = si + (N_TIME_LEVELS - time_level) * 1 + 1;
                     if (l_min < 0)        l_min = 0;
                     if (l_max > FLU_NXT ) l_max = FLU_NXT; 
                     for (l = l_min; l < l_max; ++l) s_RK1[sj][l] = true;
                  }
               }


//             2.2 sync s_RK1
#              ifdef __CUDACC__ 
               __syncthreads();
#              endif 

//             2.3 compute backward density fluxes at all cell faces of real cells
               CELL_LOOP(FLU_NXT, g2, g2 - 1)
               {
                  v = _dh * GRADB1 (s_In[sj][time_level][PHAS], si);

                  if ( s_RK1[sj][si] ) {
                     s_Fm[sj][0][si] = UPWIND_FM(s_In[sj][time_level][DENS], v, si); 
                  } else {
#                 if ( HYBRID_SCHEME == HYBRID_UPWIND )
//                   access Rc[time_level][i, i-1], Pc[time_level][i, i-1]
                     s_Fm[sj][0][si] = UPWIND_FM(s_In[sj][time_level][DENS], v, si); 
#                 elif ( HYBRID_SCHEME == HYBRID_FROMM )
//                   access Rc[time_level][i, i-1, i-2], Pc[time_level][i, i-1]
                     s_Fm[sj][0][si] = FROMM_FM (s_In[sj][time_level][DENS], v, si, Coeff1); 
#                 elif ( HYBRID_SCHEME == HYBRID_MUSCL )
//                   access Rc[time_level][i, i-1, i-2], Pc[time_level][i, i-1]
                     s_Fm[sj][0][si] = MUSCL_FM (s_In[sj][time_level][DENS], v, si, Coeff1); 
#                 elif ( HYBRID_SCHEME == HYBRID_PPM ) 
//                   access rho_L[i, i-1], rho_R[i, i-1], v_L[i]
                     s_Fm[sj][0][si] = PPM_FM   (s_In[sj][time_level][DENS], rho_L, rho_R, v_L, si, dh, dt);
#                 endif
                  }

#                 ifdef CONSERVE_MASS
//                   2.2.1 update density fluxes
                     s_Flux[sj][si] += FLUX_COEFFS[time_level] * s_Fm[sj][0][si];
#                 endif
               }


//             2.3 sync s_Fm, s_LogRho and s_Flux
#              ifdef __CUDACC__ 
               __syncthreads();
#              endif 


//             3. update density and phase
               CELL_LOOP(FLU_NXT, g2, g2)
               {
//                3.1 evolve continutiy equation with density fluxes
//                  fp  = s_Fm[sj][si+1];
//                  fm  = s_Fm[sj][si  ];

//                3.2 evolve Hamilton-Jacobi equation with Osher-Sethian flux and quantum pressure discretisation
//                  qp  = QUANTUM_PRESSURE  (s_LogRho[sj], si);
//                  vp  = GRADF3(s_In[sj][time_level][PHAS], si);
//                  vm  = GRADB3(s_In[sj][time_level][PHAS], si);
//                  osf = OSHER_SETHIAN_FLUX(vp, vm); 
                  qp     = real(1.0/2.0) * LAP2(s_LogRho[sj], si)  + real(1.0/4.0) * SQR(GRADC2(s_LogRho[sj], si));

                  if ( s_RK1[sj][si] ) {
                     vp = GRADF1(s_In[sj][time_level][PHAS], si);
                     vm = GRADB1(s_In[sj][time_level][PHAS], si);
                  } else {
                     vp = GRADF3(s_In[sj][time_level][PHAS], si);
                     vm = GRADB3(s_In[sj][time_level][PHAS], si);
                  }


//                3.1 && 3.2
                  De_New = TIME_COEFFS[time_level] * Coeff1 * ( s_Fm[sj][0][si] - s_Fm[sj][0][si+1] );
                  Ph_New = TIME_COEFFS[time_level] * Coeff2 * ( - SQR(FMIN(vp, 0)) - SQR(FMAX(vm, 0)) + qp );

//                3.3 use N_TIME_LEVELS-stages RK-algorithm
                  for (uint tl = 0; tl < time_level + 1; ++tl) {
                     De_New += RK_COEFFS[time_level][tl] * s_In[sj][tl][DENS][si];
                     Ph_New += RK_COEFFS[time_level][tl] * s_In[sj][tl][PHAS][si];
                  }
                  
//                3.4 store RK1 results on time_level 0 
                  if ( time_level == 0 ) {
                     s_In[sj][N_TIME_LEVELS][DENS][si] = FMAX(De_New, FluidMinDens);
                     s_In[sj][N_TIME_LEVELS][PHAS][si] =      Ph_New;
                     s_Fm[sj][1][si]                   = s_Fm[sj][0][si];   
                  }

//                3.5 while computing the temporary results in RK algorithm, just write them to s_In
                  if ( time_level < N_TIME_LEVELS - 1 ) {
                     s_In[sj][time_level+1][DENS][si] = De_New;
                     s_In[sj][time_level+1][PHAS][si] = Ph_New;    
                  }

//                4. write back final results to g_Fluid_In[0] or g_Fluid_Out to save memory
                  else if ( time_level == N_TIME_LEVELS - 1 ) {

//                   4.1 handle the case that the velocity timestep criterion is not met -> no update
                     if ( s_RK1[sj][si] || De_New < 0 || De_New != De_New ) {
#                       ifdef GAMER_DEBUG
                        s_RK1[sj][si] = true; 
#                       endif // # ifdef GAMER_DEBUG
                        De_New        = s_In[sj][N_TIME_LEVELS][DENS][si];
                        Ph_New        = s_In[sj][N_TIME_LEVELS][PHAS][si];
#                       ifdef CONSERVE_MASS
                        s_Flux[sj][si] = s_Fm[sj][1][si];
#                       endif // #ifdef CONSERVE_MASS
                     }

                     
//                   4.2 data
                     if ( FinalOut )
                     {

//                      apply the the minimum density check
                        De_New = (De_New < MinDens) ? MinDens : De_New;

#                       ifndef __CUDACC__ 
                        Idx2 = get1D2( k, j, si, XYZ);
#                       else // #ifndef __CUDACC_
                        Idx2 = get1D2( k, j,  i, XYZ);
#                       endif // #ifndef __CUDACC_ ... else

                        g_Fluid_Out[bx][DENS][Idx2] = De_New;
                        g_Fluid_Out[bx][PHAS][Idx2] = Ph_New;

//                      write cells that have been updated with RK1 to stub field in debug mode
#                       ifdef GAMER_DEBUG
                        g_Fluid_Out[bx][STUB][Idx2] += (real) s_RK1[sj][si];
#                       endif 

                     } else { // if ( FinalOut )       
//                      do not store negative densities
                        De_New = (De_New < 0) ? MinDens : De_New;    

//                      restore original global phase field
#                       ifdef SMOOTH_PHASE
                        for (Idx4 = 1; Idx4 <= si; ++Idx4) {
                           Ph_New -= s_2PI[sj][Idx4] * TWOPI;
                        }
#                       endif // # ifdef SMOOTH_PHASE

#                       ifndef __CUDACC__ 
                        Idx1 = get1D1( k, j, si, XYZ );
#                       endif // #ifndef __CUDACC__ 

                        g_Fluid_In[bx][DENS][Idx1] = De_New;
                        g_Fluid_In[bx][PHAS][Idx1] = Ph_New;


//                      write cells that have not been updated to stub field in debug mode
#                       ifdef GAMER_DEBUG
                        if ( k >= FLU_GHOST_SIZE && k < PS2 + FLU_GHOST_SIZE &&  j >= FLU_GHOST_SIZE && j < PS2 + FLU_GHOST_SIZE ) {
#                          ifndef __CUDACC__ 
                           Idx2 = get1D2( k, j, si, XYZ);
                           St_New = (real) s_RK1[sj][si];
#                          else // #ifndef __CUDACC_
                           Idx2 = get1D2( k, j,  i, XYZ);
                           St_New = (real) s_RK1[sj][si];
#                          endif // #ifndef __CUDACC_ ... else
                              if (XYZ == 0) {
                                 g_Fluid_Out[bx][STUB][Idx2]  = St_New;
                              } else {
                                 g_Fluid_Out[bx][STUB][Idx2] += St_New;
                              }
                        }
#                       endif 
                     } // if ( FinalOut ) ... else

#                    ifdef CONSERVE_MASS
//                   4.3 fluxes (for the flux-correction operation)
                     if ( StoreFlux  &&  tx == 0 )
                     if ( k >= FLU_GHOST_SIZE  &&  k < FLU_NXT-FLU_GHOST_SIZE )
                     if ( j >= FLU_GHOST_SIZE  &&  j < FLU_NXT-FLU_GHOST_SIZE )
                     {
                        Idx3 = __umul24( k-FLU_GHOST_SIZE, PS2 ) + (j-FLU_GHOST_SIZE);

                        g_Flux[bx][XYZ+0][0][Idx3] = s_Flux[ty][  0 + FLU_GHOST_SIZE] / Eta;
                        g_Flux[bx][XYZ+1][0][Idx3] = s_Flux[ty][PS1 + FLU_GHOST_SIZE] / Eta;
                        g_Flux[bx][XYZ+2][0][Idx3] = s_Flux[ty][PS2 + FLU_GHOST_SIZE] / Eta;
                     }
#                    endif // # ifdef CONSERVE_MASS
                  } // if ( time_level < N_TIME_LEVELS - 1 ) {
               } // CELL_LOOP(FLU_NXT, g2, g2)
#              ifdef __CUDACC__ 
               __syncthreads();
#              endif 
            } // for (time_level = 0; time_level < N_TIME_LEVELS; ++time_level) 
            
//          5.3 reset the target array indices
            j += NColumnOnce;

            if ( j >= j_end )
            {
               delta_k  = ( j - j_end )/size_j + 1;
               k       += delta_k;
               j       -= __umul24( size_j, delta_k );
            }

//          5.4 update remaining number of columns
            Column0     += NColumnOnce;
            NColumnOnce  = MIN( NColumnTotal - Column0, CGPU_FLU_BLOCK_SIZE_Y );

         } // while ( Column0 < NColumnTotal )
      }
   }
} // FUNCTION : CUFLU_Advance


# if ( HYBRID_SCHEME == HYBRID_PPM )

//Accesses a_array[i-1, i, i+1, i+2] and fills a_R_array[i] and a_L_array[i+1]
void PPM_INTERPOLATION(real* a_array, real* a_L_array, real* a_R_array, int i, int densityLimiter) {
   real a, ap, app, am, amm, delta_a, delta_ap, delta_m, delta_mp, ap2, a_L, a_R; 

// -------------------------------------------------------------------------
//  interpolate the cell-centered data to the edges
// -------------------------------------------------------------------------
      a   = a_array[i];
      ap  = a_array[i+1];
      app = a_array[i+2];
      am  = a_array[i-1];

//    Average slope of parabola
      delta_a  = 0.5 * (ap - am);
      delta_ap = 0.5 * (app - a); 

      if ((ap - a) * ( a - am )   > 0.0)
         delta_m  = FMIN(FABS(delta_a ), 2.0 * FMIN(FABS(a - am), FABS(ap  - a ))) * SGN(delta_a);
      else
         delta_m = 0.0; 


      if ((app - ap) * ( ap - a ) > 0.0)
         delta_mp = FMIN(FABS(delta_ap), 2.0 * FMIN(FABS(ap - a), FABS(app - ap))) * SGN(delta_ap);
      else
         delta_mp = 0.0; 


      if ( densityLimiter == NO_LIMITER )
         ap2 = real(7.0/12.0) * ( a + ap ) - real(1.0/12.0) * ( app + am );
      else
         ap2 = a + real(1.0/2.0) * ( ap - a ) - real(1.0/6.0) * (delta_mp - delta_m);


      a_R_array[i  ] = ap2; 
      a_L_array[i+1] = ap2; 
}

//Accesses a_array[i, i+1, i-1], a_L_array[i], a_R_array[i]
void PPM_LIMITER(real* a_array, real* a_L_array, real* a_R_array, int i, int densityLimiter) {
   real a, ap, am, a_L, a_R;
   a   = a_array[i]; 
   ap  = a_array[i+1]; 
   am  = a_array[i-1]; 
   a_R = a_R_array[i];
   a_L = a_L_array[i]; 

// Set coefficients of the interpolating parabola such that it does not overshoot
   if (densityLimiter == FULL_LIMITER) {
//    1. If a is local extremum, set the interpolation function to be constant
      if ((a_R - a)*(a - a_L) <= (0.0)) {
         a_L = a;
         a_R = a;
      }
      
//    2. If a between a_R and a_L, but very close to them, the parabola might still overshoot
      if (+ (a_R - a_L)*(a_R - a_L) / real(6.0) < (a_R - a_L) * (a - real(1.0/2.0) * (a_L + a_R)))
         a_L = real(3.0) * a - real(2.0) * a_R;
      if (- (a_R - a_L)*(a_R - a_L) / real(6.0) > (a_R - a_L) * (a - real(1.0/2.0) * (a_L + a_R)))
         a_R = real(3.0) * a - real(2.0) * a_L;

//    make sure that we didn't over or undershoot -- this may not
//    be needed, but is discussed in Colella & Sekora (2008)
      a_R = FMAX(a_R, FMIN(a, ap));
      a_R = FMIN(a_R, FMAX(a, ap));
      a_L = FMAX(a_L, FMIN(am, a));
      a_L = FMIN(a_L, FMAX(am, a));
   }


   a_R_array[i] = a_R;
   a_L_array[i] = a_L;
}



//Accesses a_array[i, i-1], a_L_array[i, i-1], a_R_array[i, i-1]
real PPM_FM(real* a_array, real* a_L_array, real* a_R_array, real* v_L_array, int i, real dh, real dt) {
   real d_a, d_am, a_6, a_6m, x, y, fm_L, fm_R;
   
   d_a  = a_R_array[i  ] - a_L_array[i  ];
   d_am = a_R_array[i-1] - a_L_array[i-1];
   a_6  = real(6.0) * (a_array[i  ] - real(1.0/2.0) * ( a_R_array[i]   + a_L_array[i]  ));
   a_6m = real(6.0) * (a_array[i-1] - real(1.0/2.0) * ( a_R_array[i-1] + a_L_array[i-1]));

   x    = + v_L_array[i] * dt / dh;
   fm_L = a_R_array[i-1] - x/2.0 * (d_am  - ( 1.0 - 2.0/3.0 * x) * a_6m);
    
   x    = - v_L_array[i] * dt / dh;
   fm_R = a_L_array[i]   + x/2.0 * (d_a   + ( 1.0 - 2.0/3.0 * x) * a_6 );

   return fm_L * FMAX(  v_L_array[i], 0.0 ) + fm_R * FMIN(  v_L_array[i], 0.0 );
}
#endif // # if ( HYBRID_SCHEME == HYBRID_PPM )


#endif // #if ( defined GPU  &&  MODEL == ELBDM && ELBDM_SCHEME == HYBRID)
