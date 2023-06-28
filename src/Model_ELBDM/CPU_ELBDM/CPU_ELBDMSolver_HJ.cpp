#include "Macro.h"
#include "CUFLU.h"

#if ( MODEL == ELBDM && ELBDM_SCHEME == ELBDM_HYBRID )

// convert to 1D index with ghost boundary
# define to1D1(z,y,x) (  (z)                 * HYB_NXT * HYB_NXT +  (y)                 * HYB_NXT +  (x)                  )
// convert to 1D index without ghost boundary
# define to1D2(z,y,x) ( ((z)-HYB_GHOST_SIZE) * PS2     * PS2     + ((y)-HYB_GHOST_SIZE) * PS2     + ((x)-HYB_GHOST_SIZE)  )

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

# ifdef HYBRID_SMOOTH_PHASE
# define  TWOPI               (  real( 2 * M_PI )  )
# define _TWOPI               (  real( 1.0 ) / TWOPI )
# define UNWRAP( ref, wrap )  (  round(((ref) - (wrap)) * _TWOPI) )
# define GRADIENT_RATIO(f, t) ( (f[t+1] - f[t  ]) / (f[t] - f[t-1] + (((f[t] - f[t-1]) == 0) ? 1e-8 : 0)))
# endif // # ifdef HYBRID_SMOOTH_PHASE

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
# define UPWIND_FM(Rc, Vb, t)        (   MAX(Vb, 0) * Rc[t-1] \
                                       + MIN(Vb, 0) * Rc[t  ] )

//Second-order MUSCL flux reconstruction
# if ( HYBRID_SCHEME == HYBRID_MUSCL )

// VAN ALBADA LIMITER ( works best, VAN LEER also good, smooth functions are better then SUPERBEE for instance )
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

# define MUSCL_FM(Rc, Vb, t, _v) (   MAX(Vb, 0) * Rc[t-1] \
                                   + MIN(Vb, 0) * Rc[t  ] \
                                   + real(0.5) * FABS(Vb) * (1. - FABS(Vb * _v)) * LIMITER(UPWIND_GRADIENT_RATIO(Rc, Vb, t)) * (Rc[t] - Rc[t - 1]) )


//Second-order unlimited flux reconstruction
# elif  ( HYBRID_SCHEME == HYBRID_FROMM )

# define FROMM_FM(Rc, Vb, t, _v) (   MAX(Vb, 0) * ( Rc[t-1] + real(0.25) * ( Rc[t  ] - Rc[t-2] ) * (1.0 - Vb * _v)) \
                                   + MIN(Vb, 0) * ( Rc[t  ] - real(0.25) * ( Rc[t+1] - Rc[t-1] ) * (1.0 + Vb * _v)) )

# endif // # if ( HYBRID_SCHEME == HYBRID_MUSCL ) ... # else

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


//Options for different Runge-Kutta schemes

// First-order method (DT_HYBRID < 0.01 ):
//#define N_TIME_LEVELS 1
//const real TIME_COEFFS[N_TIME_LEVELS]                = {1.0};
//const real RK_COEFFS  [N_TIME_LEVELS][N_TIME_LEVELS] = {{1.0}};
//#ifdef CONSERVE_MASS
//const static real FLUX_COEFFS[N_TIME_LEVELS]         = {1.0};
//#endif

// Second-order method (DT_HYBRID = 0.05 instead of 0.4):
//#define N_TIME_LEVELS 2
//GPU_DEVICE_VARIABLE
//const static real TIME_COEFFS[N_TIME_LEVELS]                = {1.0/2.0, 1.0};#
//GPU_DEVICE_VARIABLE
//const static real RK_COEFFS  [N_TIME_LEVELS][N_TIME_LEVELS] = {{1.0, 0.0}, {1.0, 0.0}};
//
//#ifdef CONSERVE_MASS
//GPU_DEVICE_VARIABLE
//const static real FLUX_COEFFS[N_TIME_LEVELS]                = {0.0, 1.0};
//#endif

#define N_TIME_LEVELS 3
GPU_DEVICE_VARIABLE
const static double TIME_COEFFS[N_TIME_LEVELS]                = {1.0, 1.0/4.0, 2.0/3.0};

GPU_DEVICE_VARIABLE
const static double RK_COEFFS  [N_TIME_LEVELS][N_TIME_LEVELS] = {{1.0, 0.0, 0.0}, {3.0/4.0, 1.0/4.0, 0.0}, {1.0/3.0, 0.0, 2.0/3.0}};

#ifdef CONSERVE_MASS
GPU_DEVICE_VARIABLE
const static double FLUX_COEFFS[N_TIME_LEVELS]                = {1.0/6.0, 1.0/6.0, 2.0/3.0};
#endif

GPU_DEVICE
static void CUFLU_Advance( real g_Fluid_In [][FLU_NIN ][ CUBE(HYB_NXT) ],
                           #ifdef GAMER_DEBUG
                           real g_Fluid_Out[][FLU_NOUT][ CUBE(PS2) ],
                           #else
                           real g_Fluid_Out[][FLU_NIN ][ CUBE(PS2) ],
                           #endif
                           real g_Flux     [][9][NFLUX_TOTAL][ SQR(PS2) ],
                           const bool g_HasWaveCounterpart [] [ CUBE(HYB_NXT) ],
                           int NPatchGroup,
                           const real dt, const real _dh, const real Eta, const bool StoreFlux,
                           const uint j_gap, const uint k_gap,
                           real s_In    [][N_TIME_LEVELS+1][FLU_NIN][HYB_NXT],
                           real s_LogRho[]                          [HYB_NXT],
                           real s_QP    []                          [HYB_NXT],
                           real s_Fm    [][2]                       [HYB_NXT],
                           real s_Flux  []                          [HYB_NXT],
                           bool s_RK1   []                          [HYB_NXT],
                           int  s_2PI   []                          [HYB_NXT],
                           const bool FinalOut, const int XYZ, const real MinDens );

//-------------------------------------------------------------------------------------------------------
// Function    :  CUFLU_ELBDMSolver_HamiltonJacobi
// Description :  GPU solver for kinetic term in Hamilton Jacobi-Madelung equations
//
// Note        :  1. The three-dimensional evolution is achieved by applying x, y, and z operators successively.
//                   Since these operators commute, the order of applying them are irrelevant.
//                   --> Input parameter "XYZ" is actually useless#
//                2. Currently supports 3 modes: HYBRID_SCHEME == HYBRID_UPWIND, HYBRID_FROMM, HYBRID_MUSCL set by global compile-time constants
//                   They respectively use first-order upwinding, second-order PLM flux reconstructions without and with limiters
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
void CUFLU_ELBDMSolver_HamiltonJacobi( real g_Fluid_In [][FLU_NIN ][ CUBE(HYB_NXT) ],
                                       #ifdef GAMER_DEBUG
                                       real g_Fluid_Out[][FLU_NOUT][ CUBE(PS2) ],
                                       #else
                                       real g_Fluid_Out[][FLU_NIN ][ CUBE(PS2) ],
                                       #endif
                                       real g_Flux     [][9][NFLUX_TOTAL][ SQR(PS2) ],
                                       const bool g_HasWaveCounterpart [] [ CUBE(HYB_NXT) ],
                                       const real dt, const real _dh, const real Eta, const bool StoreFlux,
                                       const bool XYZ, const real MinDens )
#else
void CPU_ELBDMSolver_HamiltonJacobi(   real g_Fluid_In [][FLU_NIN ][ CUBE(HYB_NXT)],
                                       #ifdef GAMER_DEBUG
                                       real g_Fluid_Out[][FLU_NOUT][ CUBE(PS2) ],
                                       #else
                                       real g_Fluid_Out[][FLU_NIN ][ CUBE(PS2) ],
                                       #endif
                                       real g_Flux     [][9][NFLUX_TOTAL][ SQR(PS2) ],
                                       const bool g_HasWaveCounterpart [] [ CUBE(HYB_NXT) ],
                                       const int NPatchGroup,
                                       const real dt, const real dh, const real Eta, const bool StoreFlux,
                                       const bool XYZ, const real MinDens )
#endif
{
#  ifdef __CUDACC__
// create memories for columns of various intermediate fields in shared GPU memory
   __shared__ real s_In      [CGPU_FLU_BLOCK_SIZE_Y][N_TIME_LEVELS + 1][FLU_NIN][HYB_NXT];
   __shared__ real s_LogRho  [CGPU_FLU_BLOCK_SIZE_Y]                            [HYB_NXT];        // 0.5 * log(rho)
   __shared__ real s_QP      [CGPU_FLU_BLOCK_SIZE_Y]                            [HYB_NXT];        // the quantum pressure
   __shared__ real s_Fm      [CGPU_FLU_BLOCK_SIZE_Y][2]                         [HYB_NXT];        // the density fluxes for every thread block
   __shared__ bool s_RK1     [CGPU_FLU_BLOCK_SIZE_Y]                            [HYB_NXT];        // booleans indicating where to switch to first-order

#  ifdef CONSERVE_MASS
   __shared__ real s_Flux    [CGPU_FLU_BLOCK_SIZE_Y]                            [HYB_NXT];         // the average density fluxes
#  else  // #  ifdef CONSERVE_MASS
              real (*s_Flux)                                                    [HYB_NXT] = NULL;  // useless if CONSERVE_MASS is off
#  endif // #  ifdef CONSERVE_MASS ... # else

#  ifdef HYBRID_SMOOTH_PHASE
   __shared__ int   s_2PI    [CGPU_FLU_BLOCK_SIZE_Y][HYB_NXT];                                     // number of phase windings between neighbouring points
#  else // # ifdef HYBRID_SMOOTH_PHASE
              int (*s_2PI)   [HYB_NXT] = NULL;                                                     // useless if HYBRID_SMOOTH_PHASE is off
#  endif // # ifdef HYBRID_SMOOTH_PHASE


   const int NPatchGroup = 0;

#  else // #  ifdef __CUDACC__
// allocate memory on stack within loop for CPU run
   real (*s_In)     [N_TIME_LEVELS + 1][FLU_NIN][HYB_NXT] = NULL;
   real (*s_LogRho)                             [HYB_NXT] = NULL;
   real (*s_QP)                                 [HYB_NXT] = NULL;
   real (*s_Fm)     [2]                         [HYB_NXT] = NULL;
   bool (*s_RK1)                                [HYB_NXT] = NULL;
   real (*s_Flux)                               [HYB_NXT] = NULL;
   int  (*s_2PI)                                [HYB_NXT] = NULL;
   const real _dh                                         = real(1.0)/dh;
#  endif

   if ( XYZ )
   {
      CUFLU_Advance( g_Fluid_In, g_Fluid_Out, g_Flux, g_HasWaveCounterpart, NPatchGroup, dt, _dh, Eta, StoreFlux,
                                  0,              0, s_In, s_LogRho, s_QP, s_Fm, s_Flux, s_RK1, s_2PI, false, 0, MinDens );
      CUFLU_Advance( g_Fluid_In, g_Fluid_Out, g_Flux, g_HasWaveCounterpart, NPatchGroup, dt, _dh, Eta, StoreFlux,
                     HYB_GHOST_SIZE,              0, s_In, s_LogRho, s_QP, s_Fm, s_Flux, s_RK1, s_2PI, false, 3, MinDens );
      CUFLU_Advance( g_Fluid_In, g_Fluid_Out, g_Flux, g_HasWaveCounterpart, NPatchGroup, dt, _dh, Eta, StoreFlux,
                     HYB_GHOST_SIZE, HYB_GHOST_SIZE, s_In, s_LogRho, s_QP, s_Fm, s_Flux, s_RK1, s_2PI,  true, 6, MinDens );
   }

   else
   {
      CUFLU_Advance( g_Fluid_In, g_Fluid_Out, g_Flux, g_HasWaveCounterpart, NPatchGroup, dt, _dh, Eta, StoreFlux,
                                  0,              0, s_In, s_LogRho, s_QP, s_Fm, s_Flux, s_RK1, s_2PI, false, 6, MinDens );
      CUFLU_Advance( g_Fluid_In, g_Fluid_Out, g_Flux, g_HasWaveCounterpart, NPatchGroup, dt, _dh, Eta, StoreFlux,
                                  0, HYB_GHOST_SIZE, s_In, s_LogRho, s_QP, s_Fm, s_Flux, s_RK1, s_2PI, false, 3, MinDens );
      CUFLU_Advance( g_Fluid_In, g_Fluid_Out, g_Flux, g_HasWaveCounterpart, NPatchGroup, dt, _dh, Eta, StoreFlux,
                     HYB_GHOST_SIZE, HYB_GHOST_SIZE, s_In, s_LogRho, s_QP, s_Fm, s_Flux, s_RK1, s_2PI,  true, 0, MinDens );
   }

} // FUNCTION : CUFLU_ELBDMSolver_HamiltonJacobi



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
//                s_LogRho       : Shared memory array to store the density logarithms
//                s_QP           : Shared memory array to store the quantum pressure
//                s_Fm           : Shared memory array to store the boundary density fluxes during the computation
//                s_Flux         : Shared memory array to store the boundary fluxes
//                s_RK1          : Shared memory array to store where to use the RK1, first-order upwind scheme
//                s_2PI          : Shared memory array to store the 2 PI winding of the phase in case the option HYBRID_SMOOTH_PHASE is turned on
//                FinalOut       : true --> store the updated data to g_Fluid_Out
//                XYZ            : 0 : Update the solution in the x direction
//                                 3 : Update the solution in the y direction
//                                 6 : Update the solution in the z direction
//                                 --> This parameter is also used to determine the place to store the output fluxes
//                MinDens        : Minimum allowed density
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE
void CUFLU_Advance(  real g_Fluid_In [][FLU_NIN  ][ CUBE(HYB_NXT) ],
                     #ifdef GAMER_DEBUG
                     real g_Fluid_Out[][FLU_NOUT ][ CUBE(PS2) ],
                     #else
                     real g_Fluid_Out[][FLU_NIN  ][ CUBE(PS2) ],
                     #endif
                     real g_Flux     [][9][NFLUX_TOTAL][ SQR(PS2) ],
                     const bool g_HasWaveCounterpart [] [ CUBE(HYB_NXT) ],
                     int NPatchGroup,
                     const real dt, const real _dh, const real Eta, const bool StoreFlux,
                     const uint j_gap, const uint k_gap,
                     real s_In     [][N_TIME_LEVELS + 1][FLU_NIN][HYB_NXT],
                     real s_LogRho []                            [HYB_NXT],
                     real s_QP     []                            [HYB_NXT],
                     real s_Fm     [][2]                         [HYB_NXT],
                     real s_Flux   []                            [HYB_NXT],
                     bool s_RK1    []                            [HYB_NXT],
                     int  s_2PI    []                            [HYB_NXT],
                     const bool FinalOut,
                     const int XYZ, const real MinDens )
{

   const real dh           = real(1.0)/_dh;                                  // grid spacing
   const real Coeff1       = real(1.0) * dt /(dh * Eta);                     // coefficient for continuity equation
   const real Coeff2       = real(0.5) * dt /(dh * dh * Eta);                // coefficient for HJ-equation
#  ifndef HYBRID_IGNORE_FLUID_FAILURE
   const real Coeff3       = real(0.5) * real(3.0) * dh * dh * Eta / dt;     // coefficient for determining velocity timestep
#  endif // # ifndef HYBRID_IGNORE_FLUID_FAILURE
   const real FluidMinDens = MAX(real(1e-10), MinDens);                      // minimum density while computing quantum pressure

   const uint j_end        = HYB_NXT -  j_gap    ;          // last y-column to be updated
   const uint size_j       = HYB_NXT - (j_gap<<1);          // number of y-columns to be updated
   const uint size_k       = HYB_NXT - (k_gap<<1);          // number of z-columns to be updated
   const uint NColumnTotal = size_j * size_k;    // total number of data columns to be updated

// openmp pragma for the CPU solver
#  ifndef __CUDACC__
#  pragma omp parallel
#  endif
   {
#     ifdef __CUDACC__
      const int bx = blockIdx.x;
#     else
//    in CPU mode, every thread works on one patch group at a time and corresponds to one block in the grid of the GPU solver
#     pragma omp for schedule( runtime ) private ( s_In, s_LogRho, s_QP, s_Fm, s_Flux, s_RK1, s_2PI )
      for (int bx=0; bx<NPatchGroup; bx++)
#     endif
      {

         uint Column0 = 0;                // the total number of columns that have been updated
         uint Idx, Idx1, Idx2;            // temporary indices used for indexing column updates, writing data to g_Fluid_In, g_Fluid_Out
#        ifdef CONSERVE_MASS
         uint Idx3;                       // temporary index used for writing data to g_Flux
#        endif // # ifdef CONSERVE_MASS
#        ifdef HYBRID_SMOOTH_PHASE
         uint Idx4;
#        endif // # ifdef HYBRID_SMOOTH_PHASE
         uint delta_k;
         uint si, sj;                     // array indices used in the shared memory array
         uint g1, g2;                     // left and right ghost zones while updating each column

         uint NStep;                      // number of iterations for updating each column
         real De_New, Ph_New, v, vp, vm;
         uint time_level;

#        ifndef HYBRID_IGNORE_FLUID_FAILURE
         int l, l_min, l_max;
#        endif // # ifndef HYBRID_IGNORE_FLUID_FAILURE

#        ifdef GAMER_DEBUG
         real St_New;                     // for storing where we use RK1 in stub field when GAMER_DEBUG is on
#        endif // # ifdef GAMER_DEBUG

#        ifdef __CUDACC__
//       use two-dimensional thread blocks in GPU mode
         const uint tx            = threadIdx.x;
         const uint ty            = threadIdx.y;
#        else  // # ifdef __CUDACC__
//       every block just has a single thread with temporary memory on the stack in CPU mode
         const uint tx            = 0;
         const uint ty            = 0;

//       create arrays for columns of various intermediate fields on the stack
         real s_In_1PG       [CGPU_FLU_BLOCK_SIZE_Y][N_TIME_LEVELS + 1][FLU_NIN][HYB_NXT]; // density and phase fields at all RK stages + RK 1 result
         real s_LogRho_1PG   [CGPU_FLU_BLOCK_SIZE_Y]                            [HYB_NXT]; // 1/2 * density logarithm
         real s_QP_1PG       [CGPU_FLU_BLOCK_SIZE_Y]                            [HYB_NXT]; // quantum pressure
         real s_Fm_1PG       [CGPU_FLU_BLOCK_SIZE_Y][2]                         [HYB_NXT]; // density flux

         bool s_RK1_1PG      [CGPU_FLU_BLOCK_SIZE_Y]                            [HYB_NXT]; // use RK1

#        ifdef CONSERVE_MASS
         real s_Flux_1PG     [CGPU_FLU_BLOCK_SIZE_Y]                            [HYB_NXT]; // boundary fluxes for fixup flux
#        endif // #  ifdef CONSERVE_MASS

#        ifdef HYBRID_SMOOTH_PHASE
         int  s_2PI_1PG      [CGPU_FLU_BLOCK_SIZE_Y]                            [HYB_NXT]; // 2PI windings
#        endif // # ifdef HYBRID_SMOOTH_PHASE

         s_In       = s_In_1PG;
         s_LogRho   = s_LogRho_1PG;
         s_QP       = s_QP_1PG;
         s_Fm       = s_Fm_1PG;
         s_RK1      = s_RK1_1PG;

#        ifdef CONSERVE_MASS
         s_Flux     = s_Flux_1PG;
#        endif // #  ifdef CONSERVE_MASS

#        ifdef HYBRID_SMOOTH_PHASE
         s_2PI      = s_2PI_1PG;
#        endif

#        endif // # ifdef __CUDACC__ ... # else

         const uint tid          = ty * CGPU_FLU_BLOCK_SIZE_X + tx;                // thread ID within block
               uint j,k;                                                           // (j,k): array indices used in g_Fluid_In
         uint NColumnOnce        = MIN( NColumnTotal, CGPU_FLU_BLOCK_SIZE_Y );     // number of columns updated per iteration
         const uint NThread      = CGPU_FLU_BLOCK_SIZE_X * CGPU_FLU_BLOCK_SIZE_Y;  // total number of threads within block


//       loop over all data columns
         while ( Column0 < NColumnTotal )
         {

            time_level = 0;

//          1. load data into shared memory
            CELL_LOOP(HYB_NXT, 0, 0) {
               j = j_gap + ( sj + Column0 ) % size_j;
               k = k_gap + ( sj + Column0 ) / size_j;

//             1.1 determine the array indices for loading global memory data along different directions
               Idx1 = get1D1( k, j, si, XYZ );

//             1.2 load the interior data into shared memory at time_level 0
               s_In[sj][time_level][DENS][si] = g_Fluid_In[bx][DENS][Idx1];
               s_In[sj][time_level][PHAS][si] = g_Fluid_In[bx][PHAS][Idx1];

#              ifdef CONSERVE_MASS
               s_Flux[sj][si] = 0;
#              endif

               s_RK1[sj][si] = false;
            }

//          1.4 sync data read into S_In
#           ifdef __CUDACC__
            __syncthreads();
#           endif // # ifdef __CUDACC_

//          1.5 smooth out 2 PI-discontinuitites in phase field
#           ifdef HYBRID_SMOOTH_PHASE

            bool hasChanged = false;

//          1.5.1. compute magnitude of 2 pi-jumps wherever we detect discontinuity
            CELL_LOOP(HYB_NXT, 1, 1)
            {
               if ( GRADIENT_RATIO(s_In[sj][0][PHAS], si) < - real(0.0) ) {
                  s_2PI[sj][si] = UNWRAP(s_In[sj][0][PHAS][si - 1], s_In[sj][0][PHAS][si]);
                  if (si == HYB_NXT - 2) {
                     s_2PI[sj][HYB_NXT - 1] = UNWRAP(s_In[sj][0][PHAS][si], s_In[sj][0][PHAS][si + 1]);
                  }
               } else {
                  s_2PI[sj][si] = 0;
                  if (si == HYB_NXT - 2) {
                     s_2PI[sj][si + 1] = 0;
                  }
               }
            }

#           ifdef __CUDACC__
            __syncthreads();
#           endif // # ifdef __CUDACC_

//          1.5.2. add multiples of 2 pi to initial phase field to smoothen it
            CELL_LOOP(HYB_NXT, 1, 0)
            {
               for (Idx4 = 1; Idx4 <= si; ++Idx4) {
                  s_In[sj][time_level][PHAS][si] += s_2PI[sj][Idx4] * TWOPI;
               }
            }

#           ifdef __CUDACC__
            __syncthreads();
#           endif // # ifdef __CUDACC_
#           endif // # ifdef HYBRID_SMOOTH_PHASE

//          2. Runge-Kutta iterations
            for (time_level = 0; time_level < N_TIME_LEVELS; ++time_level)
            {
               g1 = GHOST_ZONE_PER_STAGE *   time_level       ;
               g2 = GHOST_ZONE_PER_STAGE * ( time_level + 1 ) ;

//             2.2 compute density logarithms
               CELL_LOOP(HYB_NXT, g1, g1)
               {
                  s_LogRho[sj][si] = LOG(MAX(s_In[sj][time_level][DENS][si], FluidMinDens));
               }

//             2.3 sync _LogRho
#              ifdef __CUDACC__
               __syncthreads();
#              endif

//             2.4 check the velocity-dependent CFL-condition and switch to forward-Euler for updating the density wherever the CFL-condition is not met
               CELL_LOOP(HYB_NXT, g1 + 1, g1 + 1)
               {
//                dt = 1 / MaxdS_dx * 0.5 * ELBDM_ETA * DT__HYBRID_VELOCITY;
//                compute CFL condition timestep and quantum pressure term
                  s_QP[sj][si] = real(1.0/2.0) * LAP2(s_LogRho[sj], si)  + real(1.0/4.0) * SQR(GRADC2(s_LogRho[sj], si));

#                 ifndef HYBRID_IGNORE_FLUID_FAILURE
//                if the time step adopted in solver is larger than what velocity-dependent CFL condition allows, we switch to RK1 with a first-order upwind discretisation
                  if (FABS(s_QP[sj][si]) > 0.15 || FABS(GRADC2 (s_In[sj][time_level][PHAS], si))  > Coeff3) {
//                   compute how far wrong information can propagate
                     l_min = si - (N_TIME_LEVELS - time_level) * 1;
                     l_max = si + (N_TIME_LEVELS - time_level) * 1 + 1;
                     if (l_min < 0)        l_min = 0;
                     if (l_max > HYB_NXT ) l_max = HYB_NXT;
                     for (l = l_min; l < l_max; ++l) s_RK1[sj][l] = true;
                  }
#                 endif // # ifndef HYBRID_IGNORE_FLUID_FAILURE
               }


//             2.5 sync s_RK1
#              ifdef __CUDACC__
               __syncthreads();
#              endif

//             2.6 compute backward density fluxes at all cell faces of real cells
               CELL_LOOP(HYB_NXT, g2, g2 - 1)
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
#                 endif
                  }

#                 ifdef CONSERVE_MASS
//                   2.2.1 update density fluxes
                     s_Flux[sj][si] += FLUX_COEFFS[time_level] * s_Fm[sj][0][si];
#                 endif
               }


//             2.7 sync s_Fm and s_Flux
#              ifdef __CUDACC__
               __syncthreads();
#              endif


//             3. update density and phase
               CELL_LOOP(HYB_NXT, g2, g2)
               {
//                3.1 evolve continutiy equation with density fluxes
//                  fp  = s_Fm[sj][si+1];
//                  fm  = s_Fm[sj][si  ];

//                3.2 evolve Hamilton-Jacobi equation with Osher-Sethian flux and quantum pressure discretisation
//                  qp  = QUANTUM_PRESSURE  (s_LogRho[sj], si);
//                  vp  = GRADF3(s_In[sj][time_level][PHAS], si);
//                  vm  = GRADB3(s_In[sj][time_level][PHAS], si);
//                  osf = OSHER_SETHIAN_FLUX(vp, vm);

                  if ( s_RK1[sj][si] ) {
                     vp = GRADF1(s_In[sj][time_level][PHAS], si);
                     vm = GRADB1(s_In[sj][time_level][PHAS], si);
                  } else {
                     vp = GRADF3(s_In[sj][time_level][PHAS], si);
                     vm = GRADB3(s_In[sj][time_level][PHAS], si);
                  }


//                3.1 && 3.2
                  De_New = TIME_COEFFS[time_level] * Coeff1 * ( s_Fm[sj][0][si] - s_Fm[sj][0][si+1] );
                  Ph_New = TIME_COEFFS[time_level] * Coeff2 * ( - SQR(MIN(vp, 0)) - SQR(MAX(vm, 0)) + s_QP[sj][si] );

//                3.3 use N_TIME_LEVELS-stages RK-algorithm
                  for (uint tl = 0; tl < time_level + 1; ++tl) {
                     De_New += RK_COEFFS[time_level][tl] * s_In[sj][tl][DENS][si];
                     Ph_New += RK_COEFFS[time_level][tl] * s_In[sj][tl][PHAS][si];
                  }

//                3.4 store RK1 results on time_level 0
                  if ( time_level == 0 ) {
                     s_In[sj][N_TIME_LEVELS][DENS][si] = MAX(De_New, FluidMinDens);
                     s_In[sj][N_TIME_LEVELS][PHAS][si] =     Ph_New;
                     s_Fm[sj][1][si]                   = s_Fm[sj][0][si];
                  }

//                3.5 while computing the temporary results in RK algorithm, just write them to s_In
                  if ( time_level < N_TIME_LEVELS - 1 ) {
                     s_In[sj][time_level+1][DENS][si] = De_New;
                     s_In[sj][time_level+1][PHAS][si] = Ph_New;
                  } else {
//                3.6 handle the case that the velocity timestep criterion is not met, we detect negative density or nan -> first-order update
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
                     s_In[sj][N_TIME_LEVELS][DENS][si] = De_New;
                     s_In[sj][N_TIME_LEVELS][PHAS][si] = Ph_New;
                  }
               }

#              ifdef  __CUDACC__
               __syncthreads();
#              endif

            } // if ( time_level < N_TIME_LEVELS - 1 ) {

//          4. write back final results to g_Fluid_In[0] or g_Fluid_Out to save memory

//          4.1 write FFT array back to output array
            if ( FinalOut )
            {
               CELL_LOOP(HYB_NXT, HYB_GHOST_SIZE, HYB_GHOST_SIZE)
               {
                  j = j_gap + ( sj + Column0 ) % size_j ;
                  k = k_gap + ( sj + Column0 ) / size_j;

                  Idx2 = get1D2( k, j, si, XYZ );

                  De_New = s_In[sj][N_TIME_LEVELS][DENS][si];
                  Ph_New = s_In[sj][N_TIME_LEVELS][PHAS][si];

//                apply the the minimum density check
                  De_New = (De_New < MinDens) ? MinDens : De_New;

                  g_Fluid_Out[bx][DENS][Idx2] = De_New;
                  g_Fluid_Out[bx][PHAS][Idx2] = Ph_New;

//                write cells that have been updated with RK1 to stub field in debug mode
#                 ifdef GAMER_DEBUG
                  g_Fluid_Out[bx][STUB][Idx2] += (real) s_RK1[sj][si];
#                 endif
               }

            } else { // if ( FinalOut )

               CELL_LOOP(HYB_NXT, HYB_GHOST_SIZE, HYB_GHOST_SIZE)
               {

                  j = j_gap + ( sj + Column0 ) % size_j ;
                  k = k_gap + ( sj + Column0 ) / size_j;

                  Idx1 = get1D1( k, j, si, XYZ );

                  De_New = s_In[sj][N_TIME_LEVELS][DENS][si];
                  Ph_New = s_In[sj][N_TIME_LEVELS][PHAS][si];

//                do not store negative densities
                  De_New = (De_New < 0) ? MinDens : De_New;

//                restore original global phase field
#                 ifdef HYBRID_SMOOTH_PHASE
                  for (Idx4 = 1; Idx4 <= si; ++Idx4)
                  {
                     Ph_New -= s_2PI[sj][Idx4] * TWOPI;
                  }
#                 endif // # ifdef HYBRID_SMOOTH_PHASE


                  g_Fluid_In[bx][DENS][Idx1] = De_New;
                  g_Fluid_In[bx][PHAS][Idx1] = Ph_New;


//                   write cells that have not been updated to stub field in debug mode
#                 ifdef GAMER_DEBUG
                  if ( k >= HYB_GHOST_SIZE && k < PS2 + HYB_GHOST_SIZE &&  j >= HYB_GHOST_SIZE && j < PS2 + HYB_GHOST_SIZE ) {
                     Idx2 = get1D2( k, j, si, XYZ );
                     if (XYZ == 0) {
                        g_Fluid_Out[bx][STUB][Idx2]  = St_New;
                     } else {
                        g_Fluid_Out[bx][STUB][Idx2] += St_New;
                     }
                  }
#                 endif
               }
            } // if ( FinalOut ) ... else

//          4.3 fluxes (for the flux-correction operation)
#           ifdef CONSERVE_MASS
            if ( StoreFlux  &&  tx == 0 )
            {
               j = j_gap + ty % size_j;
               k = k_gap + ty / size_j;

               if ( k >= HYB_GHOST_SIZE  &&  k < HYB_NXT-HYB_GHOST_SIZE )
               if ( j >= HYB_GHOST_SIZE  &&  j < HYB_NXT-HYB_GHOST_SIZE )
               {
                  Idx3 = ( k - HYB_GHOST_SIZE ) * PS2 + (j-HYB_GHOST_SIZE);

                  g_Flux[bx][XYZ+0][0][Idx3] = s_Flux[ty][  0 + HYB_GHOST_SIZE] / Eta;
                  g_Flux[bx][XYZ+1][0][Idx3] = s_Flux[ty][PS1 + HYB_GHOST_SIZE] / Eta;
                  g_Flux[bx][XYZ+2][0][Idx3] = s_Flux[ty][PS2 + HYB_GHOST_SIZE] / Eta;
               }
            }
#           endif // # ifdef CONSERVE_MASS


#           ifdef __CUDACC__
            __syncthreads();
#           endif

//          4.4 reset the target array indices
            j += NColumnOnce;

            if ( j >= j_end )
            {
               delta_k  = ( j - j_end )/size_j + 1;
               k       += delta_k;
               j       -= size_j * delta_k;
            }

//          4.5 update remaining number of columns
            Column0     += NColumnOnce;
            NColumnOnce  = MIN( NColumnTotal - Column0, CGPU_FLU_BLOCK_SIZE_Y );

         } // while ( Column0 < NColumnTotal )
      }
   }
} // FUNCTION : CUFLU_Advance

#endif // #if ( defined GPU  &&  MODEL == ELBDM && ELBDM_SCHEME == ELBDM_HYBRID)
