#include "CUFLU.h"
#include "Macro.h"

#if ( ELBDM_SCHEME == ELBDM_HYBRID )


// convert to 1D index with ghost boundary
# define to1D1(z,y,x) (  (z)                 * HYB_NXT * HYB_NXT +  (y)                 * HYB_NXT +  (x)                  )
// convert to 1D index without ghost boundary
# define to1D2(z,y,x) ( ((z)-HYB_GHOST_SIZE) * PS2     * PS2     + ((y)-HYB_GHOST_SIZE) * PS2     + ((x)-HYB_GHOST_SIZE)  )

# define to1D3(z,y,x) (  (z) * PS2     * PS2                     +  (y) * PS2                     +  (x)                  )

#ifdef __CUDACC__
# define CGPU_FLU_BLOCK_SIZE_X   FLU_BLOCK_SIZE_X
# define CGPU_FLU_BLOCK_SIZE_Y   FLU_HJ_BLOCK_SIZE_Y
#else
# define CGPU_FLU_BLOCK_SIZE_X   1
# define CGPU_FLU_BLOCK_SIZE_Y   1
#endif

#ifdef __CUDACC__
#define CGPU_SHARED __shared__
#else
#define CGPU_SHARED
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

// Osher-Sethian flux for Hamilton-Jacobi equation
# define OSHER_SETHIAN_FLUX(Vel_p, Vel_m) ((real) 0.5 * (pow(MIN(Vel_p, 0), 2) + pow(MAX(Vel_m, 0), 2)))


// use third-order Runge-Kutta scheme since it provides a good compromise between the required ghost zone
// (6 for the second-order spatial discretisation) and the achievable time steps (close to finite-difference scheme)
// --> lower-order Runge-Kutta methods suffer from small time steps (CFL condition < 0.01 for first-order
//     and <= 0.05 for second-order methods according to empirical tests)
// --> using double precision for the following three constant coefficients improves the mass conservation errors
//     from ~1e-5 to ~1e-8 but deteriorates the performance by ~10%
#define ELBDM_HJ_RK_ORDER 3
GPU_DEVICE_VARIABLE
const static double TIME_COEFFS[ELBDM_HJ_RK_ORDER]                    = {1.0, 1.0/4.0, 2.0/3.0};

GPU_DEVICE_VARIABLE
const static double RK_COEFFS  [ELBDM_HJ_RK_ORDER][ELBDM_HJ_RK_ORDER] = { {1.0, 0.0, 0.0}, {3.0/4.0, 1.0/4.0, 0.0}, {1.0/3.0, 0.0, 2.0/3.0} };

#ifdef CONSERVE_MASS
GPU_DEVICE_VARIABLE
const static double FLUX_COEFFS[ELBDM_HJ_RK_ORDER]                    = {1.0/6.0, 1.0/6.0, 2.0/3.0};
#endif

// density floor for computation of quantum pressure from input density
// should be small so as to not affect the accuracy of the quantum pressure
// --> must be positive and larger than machine precision to ensure that density close to machine precision
//     does not lead to nan in logarithm
#define QP_DENSITY_FLOOR TINY_NUMBER

GPU_DEVICE
static void CUFLU_Advance( real g_Fluid_In [][FLU_NIN ][ CUBE(HYB_NXT) ],
#                          ifdef GAMER_DEBUG
                           real g_Fluid_Out[][FLU_NOUT][ CUBE(PS2) ],
#                          else
                           real g_Fluid_Out[][FLU_NIN ][ CUBE(PS2) ],
#                          endif
                           real g_Flux     [][9][NFLUX_TOTAL][ SQR(PS2) ],
                           const bool g_IsCompletelyRefined [],
                           const bool g_HasWaveCounterpart[][ CUBE(HYB_NXT) ],
                           int NPatchGroup,
                           const real dt, const real _dh, const real Eta, const bool StoreFlux,
                           const uint j_gap, const uint k_gap,
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
// Parameter   :  g_Fluid_In  : Global memory array storing the input variables
//                g_Fluid_Out : Global memory array to store the output variables
//                g_Flux      : Global memory array to store the output fluxes (useful only if StoreFlux == true)
//                dt          : Time interval to advance solution
//                _dh         : 1 / grid size
//                Eta         : Particle mass / Planck constant
//                StoreFlux   : true --> store the coarse-fine fluxes
//                                   --> useful only if CONSERVE_MASS is defined
//                XYZ         : true  : x->y->z ( forward sweep)
//                              false : z->y->x (backward sweep)
//                              --> Meaningless if CONSERVE_MASS is off since the operators along different directions
//                                  commute
//                              --> Meaningful if CONSERVE_MASS is on, in which the symmetry along different directions
//                                  are broken ...
//                MinDens     : Minimum allowed density
//-------------------------------------------------------------------------------------------------------
#ifdef __CUDACC__
__global__
void CUFLU_ELBDMSolver_HamiltonJacobi( real g_Fluid_In [][FLU_NIN ][ CUBE(HYB_NXT) ],
#                                      ifdef GAMER_DEBUG
                                       real g_Fluid_Out[][FLU_NOUT][ CUBE(PS2) ],
#                                      else
                                       real g_Fluid_Out[][FLU_NIN ][ CUBE(PS2) ],
#                                      endif
                                       real g_Flux     [][9][NFLUX_TOTAL][ SQR(PS2) ],
                                       const bool g_IsCompletelyRefined[],
                                       const bool g_HasWaveCounterpart[][ CUBE(HYB_NXT) ],
                                       const real dt, const real _dh, const real Eta, const bool StoreFlux,
                                       const bool XYZ, const real MinDens )
#else
void CPU_ELBDMSolver_HamiltonJacobi(   real g_Fluid_In [][FLU_NIN ][ CUBE(HYB_NXT)],
#                                      ifdef GAMER_DEBUG
                                       real g_Fluid_Out[][FLU_NOUT][ CUBE(PS2) ],
#                                      else
                                       real g_Fluid_Out[][FLU_NIN ][ CUBE(PS2) ],
#                                      endif
                                       real g_Flux     [][9][NFLUX_TOTAL][ SQR(PS2) ],
                                       const bool g_IsCompletelyRefined[],
                                       const bool g_HasWaveCounterpart[][ CUBE(HYB_NXT) ],
                                       const int NPatchGroup,
                                       const real dt, const real dh, const real Eta, const bool StoreFlux,
                                       const bool XYZ, const real MinDens )
#endif
{

#  ifdef __CUDACC__
// parameter useless when GPU is used
   const int NPatchGroup = NULL_INT;
#  else
   const real _dh = real(1.0)/dh;
#  endif

   if ( XYZ )
   {
      CUFLU_Advance( g_Fluid_In, g_Fluid_Out, g_Flux, g_IsCompletelyRefined, g_HasWaveCounterpart, NPatchGroup, dt, _dh, Eta, StoreFlux,
                                  0,              0, false, 0, MinDens );
      CUFLU_Advance( g_Fluid_In, g_Fluid_Out, g_Flux, g_IsCompletelyRefined, g_HasWaveCounterpart, NPatchGroup, dt, _dh, Eta, StoreFlux,
                     HYB_GHOST_SIZE,              0, false, 3, MinDens );
      CUFLU_Advance( g_Fluid_In, g_Fluid_Out, g_Flux, g_IsCompletelyRefined, g_HasWaveCounterpart, NPatchGroup, dt, _dh, Eta, StoreFlux,
                     HYB_GHOST_SIZE, HYB_GHOST_SIZE,  true, 6, MinDens );
   }
   else
   {
      CUFLU_Advance( g_Fluid_In, g_Fluid_Out, g_Flux, g_IsCompletelyRefined, g_HasWaveCounterpart, NPatchGroup, dt, _dh, Eta, StoreFlux,
                                  0,              0, false, 6, MinDens );
      CUFLU_Advance( g_Fluid_In, g_Fluid_Out, g_Flux, g_IsCompletelyRefined, g_HasWaveCounterpart, NPatchGroup, dt, _dh, Eta, StoreFlux,
                                  0, HYB_GHOST_SIZE, false, 3, MinDens );
      CUFLU_Advance( g_Fluid_In, g_Fluid_Out, g_Flux, g_IsCompletelyRefined, g_HasWaveCounterpart, NPatchGroup, dt, _dh, Eta, StoreFlux,
                     HYB_GHOST_SIZE, HYB_GHOST_SIZE,  true, 0, MinDens );
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
#                    ifdef GAMER_DEBUG
                     real g_Fluid_Out[][FLU_NOUT ][ CUBE(PS2) ],
#                    else
                     real g_Fluid_Out[][FLU_NIN  ][ CUBE(PS2) ],
#                    endif
                     real g_Flux     [][9][NFLUX_TOTAL][ SQR(PS2) ],
                     const bool g_IsCompletelyRefined[],
                     const bool g_HasWaveCounterpart [][ CUBE(HYB_NXT) ],
                     int NPatchGroup,
                     const real dt, const real _dh, const real Eta, const bool StoreFlux,
                     const uint j_gap, const uint k_gap,
                     const bool FinalOut,
                     const int XYZ, const real MinDens )
{

   const real dh           = real(1.0)/_dh;                   // grid spacing
   const real Coeff1       = real(1.0) * dt /(dh * Eta);      // coefficient for continuity equation
   const real Coeff2       = real(0.5) * dt /(dh * dh * Eta); // coefficient for HJ-equation

   const uint size_j       = HYB_NXT - 2 * j_gap;             // number of y-columns to be updated
   const uint size_k       = HYB_NXT - 2 * k_gap;             // number of z-columns to be updated
   const uint NColumnTotal = size_j * size_k;                 // total number of data columns to be updated
   const uint NThread      = CGPU_FLU_BLOCK_SIZE_X * CGPU_FLU_BLOCK_SIZE_Y;  // total number of threads in thread block


// openmp pragma for the CPU solver
#  ifndef __CUDACC__
#  pragma omp parallel
#  endif
   {
//###OPTIMIZATION: change the order of different dimensions to [N_TIME_LEVELS + 1][FLU_NIN][CGPU_FLU_BLOCK_SIZE_Y][HYB_NXT]
//    create memories for columns of various intermediate fields on stack or shared GPU memory
      CGPU_SHARED real s_In                [CGPU_FLU_BLOCK_SIZE_Y][ELBDM_HJ_RK_ORDER+1][FLU_NIN][HYB_NXT];
      CGPU_SHARED int  s_HasWaveCounterpart[CGPU_FLU_BLOCK_SIZE_Y]                              [HYB_NXT];  // booleans indicating where to switch to first-order

#     ifdef CONSERVE_MASS
      CGPU_SHARED real s_Flux              [CGPU_FLU_BLOCK_SIZE_Y]                              [HYB_NXT];  // the average density fluxes
#     endif

#     ifdef __CUDACC__
//    use two-dimensional thread blocks in GPU mode
      const uint tx            = threadIdx.x;
      const uint ty            = threadIdx.y;
      const uint tid           = ty * CGPU_FLU_BLOCK_SIZE_X + tx;    // thread ID within block
      const int  bx            = blockIdx.x;
#     else
//    every block just has a single thread with temporary memory on the stack in CPU mode
      const uint tx            = 0;
      const uint ty            = 0;
      const uint tid           = 0;
//    in CPU mode, every thread works on one patch group at a time and corresponds to one block in the grid of the GPU solver
#     pragma omp for schedule( runtime )
      for (int bx=0; bx<NPatchGroup; bx++)
#     endif // #ifdef __CUDACC__ ... else ...
      {

         uint Idx, Idx1, Idx2;            // temporary indices used for indexing column updates, writing data to g_Fluid_In, g_Fluid_Out
#        ifdef CONSERVE_MASS
         uint Idx3;                       // temporary index used for writing data to g_Flux
#        endif
         uint si, sj;                     // array indices used in the shared memory array
         uint NStep;                      // number of iterations for updating each column
         uint j,k;                        // (j,k): array indices used in g_Fluid_In
         uint Column0                   = 0;                                              // the total number of columns that have been updated
         uint NColumnOnce               = MIN( NColumnTotal, CGPU_FLU_BLOCK_SIZE_Y );     // number of columns updated per iteration
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
#           endif

//          first-order (in time and space) upwind update in completely refine patches
//          ensures that information from wave patches only propagates one cell per time step
//          part of several efforts to ensure stability of the solver when vortices form in wave patches
            if ( IsCompletelyRefined )
            {
               const uint time_level = 0;

               CELL_LOOP(HYB_NXT, HYB_GHOST_SIZE, HYB_GHOST_SIZE - 1)
               {

//                compute quantum pressure
                  const real LogRho_c   = LOG(MAX(s_In[sj][time_level][DENS][si    ], QP_DENSITY_FLOOR));
                  const real LogRho_p1  = LOG(MAX(s_In[sj][time_level][DENS][si + 1], QP_DENSITY_FLOOR));
                  const real LogRho_m1  = LOG(MAX(s_In[sj][time_level][DENS][si - 1], QP_DENSITY_FLOOR));
                  const real LogRho_Vel = LogRho_p1 - LogRho_m1;
                  const real LogRho_Lap = LogRho_p1 - 2 * LogRho_c + LogRho_m1;
                  const real QP         = real(1.0/2.0) * LogRho_Lap  + real(1.0/16.0) * SQR(LogRho_Vel);

//                compute density
                  const real vp         = GRADF1(s_In[sj][time_level][PHAS], si);
                  const real vm         = GRADB1(s_In[sj][time_level][PHAS], si);

#                 if   ( HYBRID_SCHEME == HYBRID_UPWIND )
//                access Rc[time_level][i, i-1], Pc[time_level][i, i-1]
                  const real fm         = UPWIND_FM(s_In[sj][time_level][DENS], _dh * vm, si    );
                  const real fp         = UPWIND_FM(s_In[sj][time_level][DENS], _dh * vp, si + 1);
#                 elif ( HYBRID_SCHEME == HYBRID_FROMM )
//                access Rc[time_level][i, i-1, i-2], Pc[time_level][i, i-1]
                  const real fm         = FROMM_FM (s_In[sj][time_level][DENS], _dh * vm, si    , Coeff1);
                  const real fp         = FROMM_FM (s_In[sj][time_level][DENS], _dh * vp, si + 1, Coeff1);
#                 elif ( HYBRID_SCHEME == HYBRID_MUSCL )
//                access Rc[time_level][i, i-1, i-2], Pc[time_level][i, i-1]
                  const real fm         = MUSCL_FM (s_In[sj][time_level][DENS], _dh * vm, si    , Coeff1);
                  const real fp         = MUSCL_FM (s_In[sj][time_level][DENS], _dh * vp, si + 1, Coeff1);
#                 endif

#                 ifdef CONSERVE_MASS
                  s_Flux[sj][si]        = fm;
#                 endif
                  s_In[sj][ELBDM_HJ_RK_ORDER][DENS][si] = s_In[sj][0][DENS][si] + Coeff1 * ( fm - fp );
                  s_In[sj][ELBDM_HJ_RK_ORDER][PHAS][si] = s_In[sj][0][PHAS][si] + Coeff2 * ( - SQR(MIN(vp, 0)) - SQR(MAX(vm, 0)) + QP );
               }
#              ifdef  __CUDACC__
               __syncthreads();
#              endif

            } else { // if ( IsCompletelyRefined )
//             2. Runge-Kutta iterations
               for (uint time_level=0; time_level<ELBDM_HJ_RK_ORDER; ++time_level)
               {
                  const uint ghost = GHOST_ZONE_PER_STAGE * ( time_level + 1 ) ;

//                3. update density and phase
                  CELL_LOOP(HYB_NXT, ghost, ghost)
                  {
                     real De_New, Ph_New, vp, vm, fm, fp;

//                   compute quantum pressure
                     const real LogRho_c   = LOG(MAX(s_In[sj][time_level][DENS][si    ], QP_DENSITY_FLOOR));
                     const real LogRho_p1  = LOG(MAX(s_In[sj][time_level][DENS][si + 1], QP_DENSITY_FLOOR));
                     const real LogRho_m1  = LOG(MAX(s_In[sj][time_level][DENS][si - 1], QP_DENSITY_FLOOR));
                     const real LogRho_Vel = LogRho_p1 - LogRho_m1;
                     const real LogRho_Lap = LogRho_p1 - 2 * LogRho_c + LogRho_m1;
                     const real QP         = real(1.0/2.0) * LogRho_Lap  + real(1.0/16.0) * SQR(LogRho_Vel);

//                   compute kinetic velocity
//                   make sure fluxes computed at wave-fluid boundaries agree on left and right side
                     vp = _dh * GRADF1(s_In[sj][time_level][PHAS], si);
                     vm = _dh * GRADB1(s_In[sj][time_level][PHAS], si);

#                    if   ( HYBRID_SCHEME == HYBRID_UPWIND )
//                   access Rc[time_level][i, i-1], Pc[time_level][i, i-1]
                     fm = UPWIND_FM(s_In[sj][time_level][DENS], vm, si    );
                     fp = UPWIND_FM(s_In[sj][time_level][DENS], vp, si + 1);
#                    elif ( HYBRID_SCHEME == HYBRID_FROMM )
//                   access Rc[time_level][i, i-1, i-2], Pc[time_level][i, i-1]
                     fm = FROMM_FM (s_In[sj][time_level][DENS], vm, si    , Coeff1);
                     fp = FROMM_FM (s_In[sj][time_level][DENS], vp, si + 1, Coeff1);
#                    elif ( HYBRID_SCHEME == HYBRID_MUSCL )
//                   access Rc[time_level][i, i-1, i-2], Pc[time_level][i, i-1]
                     fm = MUSCL_FM (s_In[sj][time_level][DENS], vm, si    , Coeff1);
                     fp = MUSCL_FM (s_In[sj][time_level][DENS], vp, si + 1, Coeff1);
#                    endif

#                    ifdef CONSERVE_MASS
//                   2.2.1 update density fluxes
                     s_Flux[sj][si] += FLUX_COEFFS[time_level] * fm;
                     if ( si == HYB_NXT - ghost - 1) {
                        s_Flux[sj][si + 1] += FLUX_COEFFS[time_level] * fp;
                     }
#                    endif

//                   solve wave equation to get well-defined phase update even in regions where density vanishes
                     vp = GRADF3(s_In[sj][time_level][PHAS], si);
                     vm = GRADB3(s_In[sj][time_level][PHAS], si);

//                   evolve continuity equation with density fluxes
//                   evolve Hamilton-Jacobi equation with Osher-Sethian flux and quantum pressure discretisation
                     De_New = Coeff1 * ( fm - fp );
                     Ph_New = Coeff2 * ( - SQR(MIN(vp, 0)) - SQR(MAX(vm, 0)) + QP );

                     De_New *= TIME_COEFFS[time_level];
                     Ph_New *= TIME_COEFFS[time_level];

//                   3.3 use ELBDM_HJ_RK_ORDER-stages RK-algorithm
                     for (uint tl=0; tl<time_level+1; ++tl) {
                        De_New += RK_COEFFS[time_level][tl] * s_In[sj][tl][DENS][si];
                        Ph_New += RK_COEFFS[time_level][tl] * s_In[sj][tl][PHAS][si];
                     }

//                   3.5 while computing the temporary results in RK algorithm, just write them to s_In
                     s_In[sj][time_level + 1][DENS][si] = De_New;
                     s_In[sj][time_level + 1][PHAS][si] = Ph_New;
                  } // CELL_LOOP(HYB_NXT, ghost, ghost)

#                 ifdef  __CUDACC__
                  __syncthreads();
#                 endif

               } // for (uint time_level=0; time_level<ELBDM_HJ_RK_ORDER; ++time_level)
            } // if ( isCompletelyRefined ) ... else ...

//          4. write back final results to g_Fluid_In or g_Fluid_Out to save memory

//          4.1 write shared memory fluid arrays back to global memory
            CELL_LOOP(HYB_NXT, HYB_GHOST_SIZE, HYB_GHOST_SIZE)
            {
               j = j_gap + ( sj + Column0 ) % size_j ;
               k = k_gap + ( sj + Column0 ) / size_j;

               const real Ph_Old = s_In[sj][0][PHAS][si];
               real De_New       = s_In[sj][ELBDM_HJ_RK_ORDER][DENS][si];
               real Ph_New       = s_In[sj][ELBDM_HJ_RK_ORDER][PHAS][si];

//             4.1.1 detect failure of fluid scheme
//             if cell has wave counterpart and density is negative or there is nan, use second-order finite-difference forward-in-time discretisation of wave equation
//             if cell does not have wave counterpart, unphysical field values are unexpected and should lead to the termination of the run
               if ( s_HasWaveCounterpart[sj][si] && (De_New < 0 || De_New != De_New || Ph_New != Ph_New || FABS(Ph_New - Ph_Old) > M_PI )) {
//                compute real and imaginary parts
                  const real Re_c   = SQRT(s_In[sj][0][DENS][si    ]) * COS(s_In[sj][0][PHAS][si    ]);
                  const real Re_p1  = SQRT(s_In[sj][0][DENS][si + 1]) * COS(s_In[sj][0][PHAS][si + 1]);
                  const real Re_m1  = SQRT(s_In[sj][0][DENS][si - 1]) * COS(s_In[sj][0][PHAS][si - 1]);
                  const real Im_c   = SQRT(s_In[sj][0][DENS][si    ]) * SIN(s_In[sj][0][PHAS][si    ]);
                  const real Im_p1  = SQRT(s_In[sj][0][DENS][si + 1]) * SIN(s_In[sj][0][PHAS][si + 1]);
                  const real Im_m1  = SQRT(s_In[sj][0][DENS][si - 1]) * SIN(s_In[sj][0][PHAS][si - 1]);
//                compute laplacians
                  const real Re_Lap = Re_p1 - 2 * Re_c + Re_m1;
                  const real Im_Lap = Im_p1 - 2 * Im_c + Im_m1;
//                second-order centered in space forward in time wave equation update
                  const real Re_New = Re_c - Coeff2 * Im_Lap;
                  const real Im_New = Im_c + Coeff2 * Re_Lap;
//                convert back to phase and match to old phase
                  Ph_New            = EQUALISE( Ph_Old, SATAN2(Im_New, Re_New) );
                  De_New            = SQR(Re_New) + SQR(Im_New);


#                 ifdef CONSERVE_MASS
//                4.1.2 set fluxes of cells where fluid scheme fails to zero
                  s_Flux[sj][si] = 0;
                  if ( si == HYB_NXT - HYB_GHOST_SIZE - 1) {
                     s_Flux[sj][si + 1] = 0;
                  }
#                 endif
               } // if ( s_HasWaveCounterpart[sj][si] && ... )

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

            } // CELL_LOOP(HYB_NXT, HYB_GHOST_SIZE, HYB_GHOST_SIZE)

#           ifdef  __CUDACC__
            __syncthreads();
#           endif

//          4.2 write shared memory flux arrays back to global memory
#           ifdef CONSERVE_MASS
            CELL_LOOP(HYB_NXT, HYB_GHOST_SIZE, HYB_GHOST_SIZE)
            {
               if ( StoreFlux  &&  tx == 0 )
               if ( k >= HYB_GHOST_SIZE  &&  k < HYB_NXT-HYB_GHOST_SIZE )
               if ( j >= HYB_GHOST_SIZE  &&  j < HYB_NXT-HYB_GHOST_SIZE )
               {
                  Idx3 = ( k - HYB_GHOST_SIZE ) * PS2 + ( j - HYB_GHOST_SIZE );

                  g_Flux[bx][XYZ+0][0][Idx3] = s_Flux[ty][  0 + HYB_GHOST_SIZE] / Eta;
                  g_Flux[bx][XYZ+1][0][Idx3] = s_Flux[ty][PS1 + HYB_GHOST_SIZE] / Eta;
                  g_Flux[bx][XYZ+2][0][Idx3] = s_Flux[ty][PS2 + HYB_GHOST_SIZE] / Eta;
               }
            }

#           ifdef  __CUDACC__
            __syncthreads();
#           endif
#           endif // # ifdef CONSERVE_MASS

//          4.3 update remaining number of columns
            Column0     += NColumnOnce;
            NColumnOnce  = MIN( NColumnTotal - Column0, CGPU_FLU_BLOCK_SIZE_Y );

         } // while ( Column0 < NColumnTotal )
      } // for (int bx=0; bx<NPatchGroup; bx++)
   } // #pragma omp parallel

} // FUNCTION : CUFLU_Advance



#endif // #if ( ELBDM_SCHEME == ELBDM_HYBRID )
