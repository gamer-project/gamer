#include "CUFLU.h"
#include "Macro.h"

#if ( MODEL == ELBDM && ELBDM_SCHEME == ELBDM_HYBRID )

// convert to 1D index with ghost boundary
# define to1D1(z,y,x) (  (z)                 * HYB_NXT * HYB_NXT +  (y)                 * HYB_NXT +  (x)                  )
// convert to 1D index without ghost boundary
# define to1D2(z,y,x) ( ((z)-HYB_GHOST_SIZE) * PS2     * PS2     + ((y)-HYB_GHOST_SIZE) * PS2     + ((x)-HYB_GHOST_SIZE)  )

# define to1D3(z,y,x) (  (z) * PS2     * PS2                     +  (y) * PS2                     +  (x)                  )

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

# define  TWOPI                (  real( 2 * M_PI )  )
# define _TWOPI                (  real( 1.0 ) / TWOPI )
# define UNWRAP( ref, wrap )   (  round(((ref) - (wrap)) * _TWOPI) )
# define EQUALISE( ref, wrap ) (  (wrap) + UNWRAP(ref, wrap) * TWOPI  )
# define GRADIENT_RATIO(f, t)  ( (f[t+1] - f[t  ]) / (f[t] - f[t-1] + (((f[t] - f[t-1]) == 0) ? 1e-8 : 0)))

// First-order forward and backward gradient
# define GRADF1(In, t) ( In[t + 1] - In[t    ] )
# define GRADB1(In, t) ( In[t    ] - In[t - 1] )

// First-order forward and backward gradient with phase unwrapping
# define UNWRAPGRADF1(In, t) ( In[t + 1] - EQUALISE(In[ t + 1], In[t    ] ) )
# define UNWRAPGRADB1(In, t) ( In[t    ] - EQUALISE(In[ t    ], In[t - 1] ) )

// Second-order forward and backward gradient
# define GRADF2(In, t) ( real(1.0/2.0) * ( -3*In[t] + 4*In[t+1] - In[t+2]))
# define GRADB2(In, t) ( real(1.0/2.0) * (  3*In[t] - 4*In[t-1] + In[t-2]))


// Second-order forward and backward gradient with phase unwrapping
# define UNWRAPGRADF2(In, t) ( real(1.0/2.0) * ( -3*In[t] + 4*EQUALISE(In[t], In[t+1]) - EQUALISE(In[t], In[t+2])))
# define UNWRAPGRADB2(In, t) ( real(1.0/2.0) * (  3*In[t] - 4*EQUALISE(In[t], In[t-1]) + EQUALISE(In[t], In[t-2])))

// Third-order forward and backward gradient
# define GRADF3(In, t) ( real(1.0/6.0) * ( -2*In[t-1] - 3*In[t] + 6*In[t+1] - In[t+2]))
# define GRADB3(In, t) ( real(1.0/6.0) * (  2*In[t+1] + 3*In[t] - 6*In[t-1] + In[t-2]))


// Third-order forward and backward gradient with phase unwrapping
# define UNWRAPGRADF3(In, t) ( real(1.0/6.0) * ( -2*EQUALISE(In[t], In[t-1]) - 3*In[t] + 6*EQUALISE(In[t], In[t+1]) - EQUALISE(In[t], In[t+2])))
# define UNWRAPGRADB3(In, t) ( real(1.0/6.0) * (  2*EQUALISE(In[t], In[t+1]) + 3*In[t] - 6*EQUALISE(In[t], In[t-1]) + EQUALISE(In[t], In[t-2])))

// Second- and fourth-order laplacian
# define LAP2(In, t)   ( real(1.0/1.0 ) * (            + 1 *In[t-1] - 2 *In[t] + 1 *In[t+1]             ))
# define LAP4(In, t)   ( real(1.0/12.0) * ( -1*In[t-2] + 16*In[t-1] - 30*In[t] + 16*In[t+1] - 1*In[t+2] ))

// Second- and fourth-order centered gradient
# define GRADC2(In, t) ( real(1.0/2.0 ) * (            - 1*In[t-1] + 1*In[t+1]             ))
# define GRADC4(In, t) ( real(1.0/12.0) * (  1*In[t-2] - 8*In[t-1] + 8*In[t+1] - 1*In[t+2] ))

//First-order upwind flux reconstruction
# define UPWIND_FM(Rc, Vb, t)        (   MAX(Vb, 0) * Rc[t-1] \
                                       + MIN(Vb, 0) * Rc[t  ] )

// VAN ALBADA LIMITER ( works best, VAN LEER also good, smooth functions are better then SUPERBEE for instance )
# define LIMITER(In) ( (SQR(In) + In)/((real)1. + SQR(In)) )
// MC:       # define LIMITER(r) MAX(0, MIN(MIN((1. + r) / 2., 2.), 2. * r))
// VAN LEER: # define LIMITER(r) ((r + FABS(r))/(1.0 + FABS(r)))
// SUPERBEE: # define LIMITER(r) MAX(MAX(0, MIN(2*r, 1)), MIN(r, 2))

# define UPWIND_GRADIENT_RATIO(Rc, Vb, t) ( ((Rc[t-1] - Rc[t-2]) * GTR(Vb, 0)  \
                                           + (Rc[t+1] - Rc[t  ]) * LSS(Vb, 0)) \
                                           / (Rc[t] - Rc[t-1] + (((Rc[t] - Rc[t-1]) == 0) ? 1e-8 : 0)))

# define MUSCL_FM(Rc, Vb, t, _v) (   MAX(Vb, 0) * Rc[t-1] \
                                   + MIN(Vb, 0) * Rc[t  ] \
                                   + real(0.5) * FABS(Vb) * (1. - FABS(Vb * _v)) * LIMITER(UPWIND_GRADIENT_RATIO(Rc, Vb, t)) * (Rc[t] - Rc[t - 1]) )

# define FROMM_FM(Rc, Vb, t, _v) (   MAX(Vb, 0) * ( Rc[t-1] + real(0.25) * ( Rc[t  ] - Rc[t-2] ) * (1.0 - Vb * _v)) \
                                   + MIN(Vb, 0) * ( Rc[t  ] - real(0.25) * ( Rc[t+1] - Rc[t-1] ) * (1.0 + Vb * _v)) )


# define GHOST_ZONE_PER_STAGE     2 // ghost zone per Runge-Kutta stage

//Osher-Sethian flux for Hamilton-Jacobi equation
# define OSHER_SETHIAN_FLUX(Vel_p, Vel_m) ((real) 0.5 * (pow(MIN(Vel_p, 0), 2) + pow(MAX(Vel_m, 0), 2)))


//Options for different Runge-Kutta schemes
/*
 First-order method (DT_HYBRID < 0.01 ):
#define N_TIME_LEVELS 1
const real TIME_COEFFS[N_TIME_LEVELS]                = {1.0};
const real RK_COEFFS  [N_TIME_LEVELS][N_TIME_LEVELS] = {{1.0}};
#ifdef CONSERVE_MASS
const static real FLUX_COEFFS[N_TIME_LEVELS]         = {1.0};
#endif
 Second-order method (DT_HYBRID = 0.05 instead of 0.4):
#define N_TIME_LEVELS 2
GPU_DEVICE_VARIABLE
const static real TIME_COEFFS[N_TIME_LEVELS]                = {1.0/2.0, 1.0};#
GPU_DEVICE_VARIABLE
const static real RK_COEFFS  [N_TIME_LEVELS][N_TIME_LEVELS] = {{1.0, 0.0}, {1.0, 0.0}};

#ifdef CONSERVE_MASS
GPU_DEVICE_VARIABLE
const static real FLUX_COEFFS[N_TIME_LEVELS]                = {0.0, 1.0};
#endif
*/

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
                           const bool g_IsCompletelyRefined [],
                           const bool g_HasWaveCounterpart[][ CUBE(HYB_NXT) ],
                           int NPatchGroup,
                           const real dt, const real _dh, const real Eta, const bool StoreFlux,
                           const uint j_gap, const uint k_gap,
                           real s_In    [][N_TIME_LEVELS + 1][FLU_NIN][HYB_NXT],
                           real s_Flux  [][HYB_NXT],
                           int  s_HasWaveCounterpart   [][HYB_NXT],
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
                                       const bool g_IsCompletelyRefined [],
                                       const bool g_HasWaveCounterpart[][ CUBE(HYB_NXT) ],
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
                                       const bool g_IsCompletelyRefined [],
                                       const bool g_HasWaveCounterpart[][ CUBE(HYB_NXT) ],
                                       const int NPatchGroup,
                                       const real dt, const real dh, const real Eta, const bool StoreFlux,
                                       const bool XYZ, const real MinDens )
#endif
{
#  ifdef __CUDACC__
// create memories for columns of various intermediate fields in shared GPU memory
   __shared__ real  s_In                 [CGPU_FLU_BLOCK_SIZE_Y][N_TIME_LEVELS + 1][FLU_NIN][HYB_NXT];
   __shared__ int   s_HasWaveCounterpart [CGPU_FLU_BLOCK_SIZE_Y]                            [HYB_NXT];        // booleans indicating where to switch to first-order

#  ifdef CONSERVE_MASS
   __shared__ real  s_Flux               [CGPU_FLU_BLOCK_SIZE_Y]                            [HYB_NXT];         // the average density fluxes
#  else  // #  ifdef CONSERVE_MASS
              real (*s_Flux)                                                                [HYB_NXT] = NULL;  // useless if CONSERVE_MASS is off
#  endif // #  ifdef CONSERVE_MASS ... # else

   const int NPatchGroup = 0;

#  else // #  ifdef __CUDACC__
// allocate memory on stack within loop for CPU run
   real (*s_In)     [N_TIME_LEVELS + 1][FLU_NIN][HYB_NXT] = NULL;
   int  (*s_HasWaveCounterpart)                 [HYB_NXT] = NULL;
   real (*s_Flux)                               [HYB_NXT] = NULL;
   const real _dh                                         = real(1.0)/dh;
#  endif

   if ( XYZ )
   {
      CUFLU_Advance( g_Fluid_In, g_Fluid_Out, g_Flux, g_IsCompletelyRefined, g_HasWaveCounterpart, NPatchGroup, dt, _dh, Eta, StoreFlux,
                                  0,              0, s_In, s_Flux, s_HasWaveCounterpart, false, 0, MinDens );
      CUFLU_Advance( g_Fluid_In, g_Fluid_Out, g_Flux, g_IsCompletelyRefined, g_HasWaveCounterpart, NPatchGroup, dt, _dh, Eta, StoreFlux,
                     HYB_GHOST_SIZE,              0, s_In, s_Flux, s_HasWaveCounterpart, false, 3, MinDens );
      CUFLU_Advance( g_Fluid_In, g_Fluid_Out, g_Flux, g_IsCompletelyRefined, g_HasWaveCounterpart, NPatchGroup, dt, _dh, Eta, StoreFlux,
                     HYB_GHOST_SIZE, HYB_GHOST_SIZE, s_In, s_Flux, s_HasWaveCounterpart,  true, 6, MinDens );
   }
   else
   {
      CUFLU_Advance( g_Fluid_In, g_Fluid_Out, g_Flux, g_IsCompletelyRefined, g_HasWaveCounterpart, NPatchGroup, dt, _dh, Eta, StoreFlux,
                                  0,              0, s_In, s_Flux, s_HasWaveCounterpart, false, 6, MinDens );
      CUFLU_Advance( g_Fluid_In, g_Fluid_Out, g_Flux, g_IsCompletelyRefined, g_HasWaveCounterpart, NPatchGroup, dt, _dh, Eta, StoreFlux,
                                  0, HYB_GHOST_SIZE, s_In, s_Flux, s_HasWaveCounterpart, false, 3, MinDens );
      CUFLU_Advance( g_Fluid_In, g_Fluid_Out, g_Flux, g_IsCompletelyRefined, g_HasWaveCounterpart, NPatchGroup, dt, _dh, Eta, StoreFlux,
                     HYB_GHOST_SIZE, HYB_GHOST_SIZE, s_In, s_Flux, s_HasWaveCounterpart,  true, 0, MinDens );
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
// Parameter   :  g_Fluid_In            : Global memory array storing the input variables
//                g_Fluid_Out           : Global memory array to store the output variables
//                g_Flux                : Global memory array to store the output fluxes (useful only if StoreFlux == true)
//                g_IsCompletelyRefined : Global memory array storing whether PG is completely refined
//                g_HasWaveCounterpart  : Global memory array storing which cells have wave counterpart
//                dt                    : Time interval to advance solution
//                _dh                   : 1 / grid size
//                Eta                   : Particle mass / Planck constant
//                StoreFlux             : true --> store the coarse-fine fluxes
//                                          --> useful only if CONSERVE_MASS is defined
//                j_gap                 : Number of useless grids on each side in the j direction (j may not be equal to y)
//                k_gap                 : Number of useless grids on each side in the k direction (k mya not be equal to z)
//                s_In                  : Shared memory array to store the input data and the solutions at different times
//                s_Flux                : Shared memory array to store the boundary fluxes
//                s_HasWaveCounterpart  : Shared memory array to store where cells have wave counterpart
//                FinalOut              : true --> store the updated data to g_Fluid_Out
//                XYZ                   : 0 : Update the solution in the x direction
//                                        3 : Update the solution in the y direction
//                                        6 : Update the solution in the z direction
//                                        --> This parameter is also used to determine the place to store the output fluxes
//                MinDens               : Minimum allowed density
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE
void CUFLU_Advance(  real g_Fluid_In [][FLU_NIN  ][ CUBE(HYB_NXT) ],
                     #ifdef GAMER_DEBUG
                     real g_Fluid_Out[][FLU_NOUT ][ CUBE(PS2) ],
                     #else
                     real g_Fluid_Out[][FLU_NIN  ][ CUBE(PS2) ],
                     #endif
                     real g_Flux     [][9][NFLUX_TOTAL][ SQR(PS2) ],
                     const bool g_IsCompletelyRefined[],
                     const bool g_HasWaveCounterpart [][ CUBE(HYB_NXT) ],
                     int NPatchGroup,
                     const real dt, const real _dh, const real Eta, const bool StoreFlux,
                     const uint j_gap, const uint k_gap,
                     real s_In                 [][N_TIME_LEVELS + 1][FLU_NIN][HYB_NXT],
                     real s_Flux               [][HYB_NXT],
                     int  s_HasWaveCounterpart [][HYB_NXT],
                     const bool FinalOut,
                     const int XYZ, const real MinDens )
{

   const real dh           = real(1.0)/_dh;                   // grid spacing
   const real Coeff1       = real(1.0) * dt /(dh * Eta);      // coefficient for continuity equation
   const real Coeff2       = real(0.5) * dt /(dh * dh * Eta); // coefficient for HJ-equation

   const uint size_j       = HYB_NXT - 2 * j_gap;             // number of y-columns to be updated
   const uint size_k       = HYB_NXT - 2 * k_gap;             // number of z-columns to be updated
   const uint NColumnTotal = size_j * size_k;                 // total number of data columns to be updated

// openmp pragma for the CPU solver
#  ifndef __CUDACC__
#  pragma omp parallel
#  endif
   {
#     ifdef __CUDACC__
      const int bx = blockIdx.x;
#     else
//    in CPU mode, every thread works on one patch group at a time and corresponds to one block in the grid of the GPU solver
#     pragma omp for schedule( runtime ) private ( s_In, s_Flux, s_HasWaveCounterpart )
      for (int bx=0; bx<NPatchGroup; bx++)
#     endif
      {

#        ifdef GAMER_DEBUG
/*
if (XYZ == 0)
{
         int counter = 0;
         for (int k=0; k<CUBE(HYB_NXT); k++)
         {
            if (g_HasWaveCounterpart[bx][k]) counter++;
         }

         printf("PG with bx = %d counter %d\n", bx, counter);
         int i = 10;
         for (int j=0; j<HYB_NXT; j++)
         {
            for (int k=0; k<HYB_NXT; k++)
            {
               const int t1 = to1D1( i, j, k );

               if (g_HasWaveCounterpart[bx][t1])
                  printf("X ");
               else
                  printf("O ");
            }
            printf("\n");
         }
         printf("\n\n");
}*/
#        endif

         uint Column0 = 0;                // the total number of columns that have been updated
         uint Idx, Idx1, Idx2;            // temporary indices used for indexing column updates, writing data to g_Fluid_In, g_Fluid_Out
#        ifdef CONSERVE_MASS
         uint Idx3;                       // temporary index used for writing data to g_Flux
#        endif // # ifdef CONSERVE_MASS
         uint si, sj;                     // array indices used in the shared memory array
         uint NStep;                      // number of iterations for updating each column

#        ifdef __CUDACC__
//       use two-dimensional thread blocks in GPU mode
         const uint tx            = threadIdx.x;
         const uint ty            = threadIdx.y;
#        else  // # ifdef __CUDACC__
//       every block just has a single thread with temporary memory on the stack in CPU mode
         const uint tx            = 0;
         const uint ty            = 0;

//       create arrays for columns of various intermediate fields on the stack
         real s_In_1PG                      [CGPU_FLU_BLOCK_SIZE_Y][N_TIME_LEVELS + 1][FLU_NIN][HYB_NXT]; // density and phase fields at all RK stages
         int  s_HasWaveCounterpart_1PG      [CGPU_FLU_BLOCK_SIZE_Y]                            [HYB_NXT]; // use RK1

#        ifdef CONSERVE_MASS
         real s_Flux_1PG                    [CGPU_FLU_BLOCK_SIZE_Y]                            [HYB_NXT]; // boundary fluxes for fixup flux
#        endif // #  ifdef CONSERVE_MASS

         s_In                      = s_In_1PG;
         s_HasWaveCounterpart      = s_HasWaveCounterpart_1PG;

#        ifdef CONSERVE_MASS
         s_Flux     = s_Flux_1PG;
#        endif // #  ifdef CONSERVE_MASS

#        endif // # ifdef __CUDACC__ ... # else

         const uint tid                 = ty * CGPU_FLU_BLOCK_SIZE_X + tx;                // thread ID within block
               uint j,k;                                                             // (j,k): array indices used in g_Fluid_In
         uint NColumnOnce               = MIN( NColumnTotal, CGPU_FLU_BLOCK_SIZE_Y );     // number of columns updated per iteration
         const uint NThread             = CGPU_FLU_BLOCK_SIZE_X * CGPU_FLU_BLOCK_SIZE_Y;  // total number of threads within block
         const bool IsCompletelyRefined = g_IsCompletelyRefined[bx];


//       loop over all data columns
         while ( Column0 < NColumnTotal )
         {
//          1. load data into shared memory
            CELL_LOOP(HYB_NXT, 0, 0) {
               j = j_gap + ( sj + Column0 ) % size_j;
               k = k_gap + ( sj + Column0 ) / size_j;

//             1.1 determine the array indices for loading global memory data along different directions
               Idx1 = get1D1( k, j, si, XYZ );

//             1.2 load the interior data into shared memory at time_level 0
               s_In[sj][0][DENS][si] = g_Fluid_In[bx][DENS][Idx1];
               s_In[sj][0][PHAS][si] = g_Fluid_In[bx][PHAS][Idx1];

#              ifdef CONSERVE_MASS
               s_Flux[sj][si] = 0;
#              endif

               s_HasWaveCounterpart[sj][si]  = g_HasWaveCounterpart[bx][Idx1];
            }

//          1.4 sync data read into S_In
#           ifdef __CUDACC__
            __syncthreads();
#           endif // # ifdef __CUDACC_

//          2. Runge-Kutta iterations
            for (uint time_level = 0; time_level < N_TIME_LEVELS; ++time_level)
            {
               const uint ghost = GHOST_ZONE_PER_STAGE * ( time_level + 1 ) ;

//             3. update density and phase
               CELL_LOOP(HYB_NXT, ghost, ghost)
               {

                  real De_New, Ph_New, vp, vm, fm, fp;

//                compute density logarithms
                  const real LogRho_c   = LOG(s_In[sj][time_level][DENS][si    ]);
                  const real LogRho_p1  = LOG(s_In[sj][time_level][DENS][si + 1]);
                  const real LogRho_m1  = LOG(s_In[sj][time_level][DENS][si - 1]);

//                compute quantum pressure
                  const real LogRho_Vel = LogRho_p1 - LogRho_m1;
                  const real LogRho_Lap = LogRho_p1 - 2 * LogRho_c + LogRho_m1;
                  const real QP         = real(1.0/2.0) * LogRho_Lap  + real(1.0/16.0) * SQR(LogRho_Vel);


//                compute kinetic velocity
//                make sure fluxes computed at wave-fluid boundaries agree on left and right side
                  if ( s_HasWaveCounterpart[sj][si - 1] && s_HasWaveCounterpart[sj][si] && s_HasWaveCounterpart[sj][si + 1]) {
                     vp = _dh * UNWRAPGRADF1(s_In[sj][time_level][PHAS], si);
                     vm = _dh * UNWRAPGRADB1(s_In[sj][time_level][PHAS], si);
                  } else {
                     vp = _dh * GRADF1(s_In[sj][time_level][PHAS], si);
                     vm = _dh * GRADB1(s_In[sj][time_level][PHAS], si);
                  }

#                 if ( HYBRID_SCHEME == HYBRID_UPWIND )
//                access Rc[time_level][i, i-1], Pc[time_level][i, i-1]
                  fm = UPWIND_FM(s_In[sj][time_level][DENS], vm, si    );
                  fp = UPWIND_FM(s_In[sj][time_level][DENS], vp, si + 1);
#                 elif ( HYBRID_SCHEME == HYBRID_FROMM )
//                access Rc[time_level][i, i-1, i-2], Pc[time_level][i, i-1]
                  fm = FROMM_FM (s_In[sj][time_level][DENS], vm, si    , Coeff1);
                  fp = FROMM_FM (s_In[sj][time_level][DENS], vp, si + 1, Coeff1);
#                 elif ( HYBRID_SCHEME == HYBRID_MUSCL )
//                access Rc[time_level][i, i-1, i-2], Pc[time_level][i, i-1]
                  fm = MUSCL_FM (s_In[sj][time_level][DENS], vm, si    , Coeff1);
                  fp = MUSCL_FM (s_In[sj][time_level][DENS], vp, si + 1, Coeff1);
#                 endif


#                 ifdef CONSERVE_MASS
//                2.2.1 update density fluxes
                  s_Flux[sj][si] += FLUX_COEFFS[time_level] * fm;
                  if ( si == HYB_NXT - ghost - 1) {
                     s_Flux[sj][si + 1] += FLUX_COEFFS[time_level] * fp;
                  }
#                 endif


//                solve wave equation to get well-defined phase update even in regions where density vanishes
                  if ( s_HasWaveCounterpart[sj][si] )
                  {
                     vp = UNWRAPGRADF1(s_In[sj][time_level][PHAS], si);
                     vm = UNWRAPGRADB1(s_In[sj][time_level][PHAS], si);
                  } else if ( FABS(QP) > 0.15 )
                  {
                     vp = GRADF1(s_In[sj][time_level][PHAS], si);
                     vm = GRADB1(s_In[sj][time_level][PHAS], si);
                  }  else
                  {
                     vp = GRADF3(s_In[sj][time_level][PHAS], si);
                     vm = GRADB3(s_In[sj][time_level][PHAS], si);
                  }

//                evolve continuity equation with density fluxes
//                evolve Hamilton-Jacobi equation with Osher-Sethian flux and quantum pressure discretisation
                  De_New = TIME_COEFFS[time_level] * Coeff1 * ( fm - fp );
                  Ph_New = TIME_COEFFS[time_level] * Coeff2 * ( - SQR(MIN(vp, 0)) - SQR(MAX(vm, 0)) + QP );

//                3.3 use N_TIME_LEVELS-stages RK-algorithm
                  for (uint tl = 0; tl < time_level + 1; ++tl) {
                     De_New += RK_COEFFS[time_level][tl] * s_In[sj][tl][DENS][si];
                     Ph_New += RK_COEFFS[time_level][tl] * s_In[sj][tl][PHAS][si];
                  }

                  if ( De_New < 0 || De_New != De_New || Ph_New != Ph_New ) {
                     De_New = s_In[sj][0][DENS][si];
                     Ph_New = s_In[sj][0][PHAS][si];
                  }

//                3.5 while computing the temporary results in RK algorithm, just write them to s_In
                  s_In[sj][time_level + 1][DENS][si] = De_New;
                  s_In[sj][time_level + 1][PHAS][si] = Ph_New;
               }

#              ifdef  __CUDACC__
               __syncthreads();
#              endif

            } // if ( time_level < N_TIME_LEVELS - 1 ) {

//          4. write back final results to g_Fluid_In or g_Fluid_Out to save memory

//          4.1 write shared memory array back to global memory
            CELL_LOOP(HYB_NXT, HYB_GHOST_SIZE, HYB_GHOST_SIZE)
            {

               j = j_gap + ( sj + Column0 ) % size_j ;
               k = k_gap + ( sj + Column0 ) / size_j;

               real De_New = s_In[sj][N_TIME_LEVELS][DENS][si];
               real Ph_New = s_In[sj][N_TIME_LEVELS][PHAS][si];

               if ( FinalOut )
               {
//                apply the the minimum density check
                  De_New = (De_New < MinDens) ? MinDens : De_New;

                  Idx2 = get1D2( k, j, si, XYZ );

                  g_Fluid_Out[bx][DENS][Idx2] = De_New;
                  g_Fluid_Out[bx][PHAS][Idx2] = Ph_New;
#                 ifdef GAMER_DEBUG
                  g_Fluid_Out[bx][STUB][Idx2] = s_HasWaveCounterpart[sj][si];
#                 endif
               } else {
                  Idx1 = get1D1( k, j, si, XYZ );

                  g_Fluid_In[bx][DENS][Idx1] = De_New;
                  g_Fluid_In[bx][PHAS][Idx1] = Ph_New;
               }


//             4.2 fluxes (for the flux-correction operation)
#              ifdef CONSERVE_MASS
               if ( StoreFlux  &&  tx == 0 )
               if ( k >= HYB_GHOST_SIZE  &&  k < HYB_NXT-HYB_GHOST_SIZE )
               if ( j >= HYB_GHOST_SIZE  &&  j < HYB_NXT-HYB_GHOST_SIZE )
               {
                  Idx3 = ( k - HYB_GHOST_SIZE ) * PS2 + ( j - HYB_GHOST_SIZE );

                  g_Flux[bx][XYZ+0][0][Idx3] = s_Flux[ty][  0 + HYB_GHOST_SIZE] / Eta;
                  g_Flux[bx][XYZ+1][0][Idx3] = s_Flux[ty][PS1 + HYB_GHOST_SIZE] / Eta;
                  g_Flux[bx][XYZ+2][0][Idx3] = s_Flux[ty][PS2 + HYB_GHOST_SIZE] / Eta;
               }
#              endif // # ifdef CONSERVE_MASS
            }

#           ifdef __CUDACC__
            __syncthreads();
#           endif

//          4.3 update remaining number of columns
            Column0     += NColumnOnce;
            NColumnOnce  = MIN( NColumnTotal - Column0, CGPU_FLU_BLOCK_SIZE_Y );

         } // while ( Column0 < NColumnTotal )
      }
   }
} // FUNCTION : CUFLU_Advance

#endif // #if ( defined GPU  &&  MODEL == ELBDM && ELBDM_SCHEME == ELBDM_HYBRID)
