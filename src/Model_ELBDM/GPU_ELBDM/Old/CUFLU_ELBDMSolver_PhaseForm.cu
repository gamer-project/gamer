#include "Macro.h"
#include "CUFLU.h"
#include "stdio.h"

#if ( defined GPU  &&  MODEL == ELBDM && ELBDM_SCHEME == HYBRID && ( HYBRID_SCHEME == HYBRID_UPWIND || HYBRID_SCHEME == HYBRID_MUSCL ) )



// useful macros
# define to1D1(z,y,x) ( __umul24(z, FLU_NXT*FLU_NXT) + __umul24(y, FLU_NXT) + x )
# define to1D2(z,y,x) ( __umul24(z-FLU_GHOST_SIZE, PS2*PS2) + __umul24(y-FLU_GHOST_SIZE, PS2) + x-FLU_GHOST_SIZE )


# define CELL_LOOP( NCell, leftGhost, rightGhost )    for ( (Idx   = tid, (NStep = (NCell) - (leftGhost) - (rightGhost), (si = Idx % NStep + (leftGhost), sj = Idx / NStep))); \
                                                             Idx   < NColumnOnce * NStep; \
                                                            (Idx  += NThread, (si = Idx % NStep + (leftGhost), sj = Idx / NStep)) )

# define GTR( a, b )     (  ( (a) > (b) ) ? (1) : (0)  )
# define LSS( a, b )     (  ( (a) < (b) ) ? (1) : (0)  )
# define SGN( a )        (  ( (a) > (0) ) ? (1) : ( (a) < (0) ) ? (-1) : (0) )
            
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
# if ( HYBRID_SCHEME == HYBRID_UPWIND )

# define UPWIND_FM(Rc, Vb, t)        (   FMAX(Vb, 0) * Rc[t-1] \
                                       + FMIN(Vb, 0) * Rc[t  ] )

# endif // #if ( HYBRID_SCHEME == HYBRID_UPWIND )

//Second-order MUSCL flux reconstruction
#if ( HYBRID_SCHEME == HYBRID_MUSCL )

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

# define MUSCL_FM(Rc, Vb, t, dx, dt) (   FMAX(Vb, 0) * Rc[t-1] \
                                       + FMIN(Vb, 0) * Rc[t  ] \
                                       +  real(0.5) * FABS(Vb) * (1. - FABS(Vb * dt/dx)) * LIMITER(UPWIND_GRADIENT_RATIO(Rc, Vb, t)) * (Rc[t] - Rc[t - 1]) )
# endif // #if ( HYBRID_SCHEME == HYBRID_MUSCL )

//Third-order PPM flux reconstruction
# if ( HYBRID_SCHEME == HYBRID_PPM )
void PPM_INTERPOLATION(real* a_array, real* a_L_array, real* a_R_array, int i, int densityLimiter);
void PPM_LIMITER      (real* a_array, real* a_L_array, real* a_R_array, int i, int densityLimiter);
real PPM_FM           (real* a_array, real* a_L_array, real* a_R_array, real* v_L_array, int i, real dh, real dt);
# endif // # if ( HYBRID_SCHEME == HYBRID_PPM )

//Second- and fourth-order quantum pressure depending on the scheme used
# if (HYBRID_SCHEME == HYBRID_MUSCL || HYBRID_SCHEME == HYBRID_UPWIND)
# define QUANTUM_PRESSURE(Sr, t)  ((real) 0.5 * (pow(GRADC2(Sr, t), 2) + LAP2(Sr, t)))
# else 
# define QUANTUM_PRESSURE(Sr, t)  ((real) 0.5 * (pow(GRADC4(Sr, t), 2) + LAP4(Sr, t)))
# endif // # if (HYBRID_SCHEME == HYBRID_MUSCL || HYBRID_SCHEME == HYBRID_UPWIND)

//Osher-Sethian flux for Hamilton-Jacobi equation
# define OSHER_SETHIAN_FLUX(vp, vm) ((real) 0.5 * (pow(MIN(vp, 0), 2) + pow(MAX(vm, 0), 2)))

//#define N_TIME_LEVELS 1
//const real TIME_COEFFS[N_TIME_LEVELS]                = {1.0};
//const real RK_COEFFS  [N_TIME_LEVELS][N_TIME_LEVELS] = {{1.0}};


//#define N_TIME_LEVELS 2
//const real TIME_COEFFS[N_TIME_LEVELS]                = {1.0/2.0, 1.0};
//const real RK_COEFFS  [N_TIME_LEVELS][N_TIME_LEVELS] = {{1.0, 0.0}, {1.0, 0.0}};

#define N_TIME_LEVELS 3
const static __device__ real TIME_COEFFS[N_TIME_LEVELS]                = {1.0, 1.0/4.0, 2.0/3.0};
const static __device__ real RK_COEFFS  [N_TIME_LEVELS][N_TIME_LEVELS] = {{1.0, 0.0, 0.0}, {3.0/4.0, 1.0/4.0, 0.0}, {1.0/3.0, 0.0, 2.0/3.0}};

#ifdef CONSERVE_MASS
const static __device__ real FLUX_COEFFS[N_TIME_LEVELS]                = {1.0/6.0, 1.0/6.0, 2.0/3.0};
#endif 

#define NO_LIMITER    0
#define SOME_LIMITER  1
#define FULL_LIMITER  2 

static __device__ void CUFLU_Advance( real g_Fluid_In [][FLU_NIN ][ CUBE(FLU_NXT) ],
                                      real g_Fluid_Out[][FLU_NOUT][ CUBE(PS2) ],
                                      real g_Flux     [][9][NFLUX_TOTAL][ SQR(PS2) ],
                                      const real dt, const real _dh, const real Eta, const bool StoreFlux,
                                      const uint j_gap, const uint k_gap,         
                                      real s_In[][N_TIME_LEVELS + 1][FLU_NIN][FLU_NXT],
                                      real s_Sr[][FLU_NXT],
                                      real s_Fm[][FLU_NXT],
                                      real s_Flux[][FLU_NXT], 
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



__global__ void CUFLU_ELBDMSolver_PhaseForm( 
                                    real g_Fluid_In [][FLU_NIN ][ CUBE(FLU_NXT) ],
                                    real g_Fluid_Out[][FLU_NOUT][ CUBE(PS2) ],
                                    real g_Flux     [][9][NFLUX_TOTAL][ SQR(PS2) ],
                                    const real dt, const real _dh, const real Eta, const bool StoreFlux,
                                    const bool XYZ, const real MinDens )
{


   __shared__ real s_In  [FLU_BLOCK_SIZE_Y][N_TIME_LEVELS + 1][FLU_NIN][FLU_NXT];
   __shared__ real s_Sr  [FLU_BLOCK_SIZE_Y][FLU_NXT]; // one column of the 0.5 * log(rho) for every thread block
   __shared__ real s_Fm  [FLU_BLOCK_SIZE_Y][FLU_NXT]; // one column of the fluxes for every thread block

#  ifdef CONSERVE_MASS
   __shared__ real s_Flux [FLU_BLOCK_SIZE_Y][FLU_NXT];
#  else
              real (*s_Flux)[FLU_NXT] = NULL;  // useless if CONSERVE_MASS is off
#  endif

   if ( XYZ )
   {
      CUFLU_Advance( g_Fluid_In, g_Fluid_Out, g_Flux, dt, _dh, Eta, StoreFlux,
                                  0,              0, s_In, s_Sr, s_Fm, s_Flux, false, 0, MinDens );
      CUFLU_Advance( g_Fluid_In, g_Fluid_Out, g_Flux, dt, _dh, Eta, StoreFlux,
                     FLU_GHOST_SIZE,              0, s_In, s_Sr, s_Fm, s_Flux, false, 3, MinDens );
      CUFLU_Advance( g_Fluid_In, g_Fluid_Out, g_Flux, dt, _dh, Eta, StoreFlux,
                     FLU_GHOST_SIZE, FLU_GHOST_SIZE, s_In, s_Sr, s_Fm, s_Flux,  true, 6, MinDens );
   }

   else
   {
      CUFLU_Advance( g_Fluid_In, g_Fluid_Out, g_Flux, dt, _dh, Eta, StoreFlux,
                                  0,              0, s_In, s_Sr, s_Fm, s_Flux, false, 6, MinDens );
      CUFLU_Advance( g_Fluid_In, g_Fluid_Out, g_Flux, dt, _dh, Eta, StoreFlux,
                                  0, FLU_GHOST_SIZE, s_In, s_Sr, s_Fm, s_Flux, false, 3, MinDens );
      CUFLU_Advance( g_Fluid_In, g_Fluid_Out, g_Flux, dt, _dh, Eta, StoreFlux,
                     FLU_GHOST_SIZE, FLU_GHOST_SIZE, s_In, s_Sr, s_Fm, s_Flux,  true, 0, MinDens );
   }

} // FUNCTION : CUFLU_ELBDMSolver_PhaseForm_MUSCL



//-------------------------------------------------------------------------------------------------------
// Function    :  CUFLU_Advance
// Description :  Use GPU to advance solutions by one time-step
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
//                s_Sr           : Shared memory array to store the density logarithms
//                s_Fm           : Shared memory array to store the boundary density fluxes during the computation
//                s_Flux         : Shared memory array to store the boundary fluxes
//                FinalOut       : true --> store the updated data to g_Fluid_Out
//                XYZ            : 0 : Update the solution in the x direction
//                                 3 : Update the solution in the y direction
//                                 6 : Update the solution in the z direction
//                                 --> This parameter is also used to determine the place to store the output fluxes
//                MinDens        : Minimum allowed density
//-------------------------------------------------------------------------------------------------------

__device__ void CUFLU_Advance( real g_Fluid_In [][FLU_NIN ][ CUBE(FLU_NXT) ],
                               real g_Fluid_Out[][FLU_NOUT][ CUBE(PS2) ],
                               real g_Flux     [][9][NFLUX_TOTAL][ SQR(PS2) ],
                               const real dt, const real _dh, const real Eta, const bool StoreFlux,
                               const uint j_gap, const uint k_gap, 
                               real s_In    [][N_TIME_LEVELS + 1][FLU_NIN][FLU_NXT],
                               real s_Sr    [][FLU_NXT],
                               real s_Fm    [][FLU_NXT],
                               real s_Flux  [][FLU_NXT], 
                               const bool FinalOut,
                               const int XYZ, const real MinDens )
{

   const real dh      = 1./_dh; 
   const real Coeff1 = dt/(dh * Eta);
   const real Coeff2 = dt/(dh * dh * Eta);

   const uint bx           = blockIdx.x;
   const uint tx           = threadIdx.x;
   const uint ty           = threadIdx.y;
   const uint tid          = __umul24(ty,FLU_BLOCK_SIZE_X) + tx;
   const uint size_j       = FLU_NXT - (j_gap<<1);
   const uint size_k       = FLU_NXT - (k_gap<<1);
   const uint NColumnTotal = __umul24( size_j, size_k );    // total number of data columns to be updated
   const uint i            = tx + FLU_GHOST_SIZE;           // (i,j,k): array indices used in g_Fluid_In
   const uint j_end        = FLU_NXT - j_gap;
         uint j            = j_gap + ty%size_j;
         uint k            = k_gap + ty/size_j;
         uint Column0      = 0;                              // the total number of columns that have been updated
         uint NColumnOnce  = MIN( NColumnTotal, FLU_BLOCK_SIZE_Y );

   uint   Idx, Idx1, Idx2, delta_k;

   real FluidMinDens = FMAX(1e-10, MinDens); 

   const uint NThread     = FLU_BLOCK_SIZE_X*FLU_BLOCK_SIZE_Y;

   real vp, vm, fp, fm, qp, osf;
   //change of density and phase in time step
   real De_New, Ph_New;

   int time_level; 

#  if ( HYBRID_SCHEME == HYBRID_MUSCL || HYBRID_SCHEME == HYBRID_UPWIND)
   const int ghostZonePerStage = 2; 
#  else 
   return; 
#  endif 

   uint   si, sj; // array indices used in the shared memory array
   uint   g1, g2; 

   uint NStep;

#  ifdef CONSERVE_MASS
   uint Idx3, f;
#  endif

// determine the array indices for loading the ghost-zone data
   bool LoadGhost = false;                                  // true --> load the ghost-zone data
   uint LoadGhost_i;
   int  LoadGhost_di, LoadGhost_dIdx1;


// use the first 2*FLU_GHOST_SIZE threads to load the ghost zones in addition to their regular cells 
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


// loop over all data columns
   while ( Column0 < NColumnTotal )
   {
//    1. load data into shared memory
      if ( tid < NColumnOnce*PS2 )
      {
//       1.1 determine the array indices for loading global memory data along different directions
         switch ( XYZ )
         {
            case 0:  Idx1 = to1D1( k, j, i );    break;
            case 3:  Idx1 = to1D1( k, i, j );    break;
            case 6:  Idx1 = to1D1( i, k, j );    break;
         }

         time_level = 0; 
      
//       1.2 load the interior data into shared memory at time_level 0
         s_In[ty][time_level][DENS][i] = g_Fluid_In[bx][DENS][Idx1];
         s_In[ty][time_level][PHAS][i] = g_Fluid_In[bx][PHAS][Idx1];

//       1.3 load the ghost-zone data into shared memory
         if ( LoadGhost )
         {
            s_In[ty][time_level][DENS][LoadGhost_i] = g_Fluid_In[bx][DENS][ (int)Idx1 + LoadGhost_dIdx1 ];
            s_In[ty][time_level][PHAS][LoadGhost_i] = g_Fluid_In[bx][PHAS][ (int)Idx1 + LoadGhost_dIdx1 ];
         }
      } // if ( tid < NColumnOnce*PS2 )

#     ifdef __CUDACC__ 
      __syncthreads();
#     endif 

      // RK iterations
      for (time_level = 0; time_level < N_TIME_LEVELS; ++time_level) 
      {
         g1 = ghostZonePerStage *   time_level       ;
         g2 = ghostZonePerStage * ( time_level + 1 ) ;         

//       compute density logarithms
         CELL_LOOP(FLU_NXT, g1, g1)
         { 
            s_Sr[sj][si] = real(0.5) * log(FMAX(s_In[sj][time_level][DENS][si], FluidMinDens));
         } 

//       compute backward density fluxes at all cell faces of real cells
         CELL_LOOP(FLU_NXT, g2, g2 - 1)
         {
#           if ( HYBRID_SCHEME == HYBRID_UPWIND )
//             access Rc[time_level][i, i-1], Pc[time_level][i, i-1]
               s_Fm[sj][si] = UPWIND_FM(s_In[sj][time_level][DENS], _dh * GRADB1 (s_In[sj][time_level][PHAS], si), si); 
#           elif ( HYBRID_SCHEME == HYBRID_MUSCL )
//             access Rc[time_level][i, i-1, i-2], Pc[time_level][i, i-1]
               s_Fm[sj][si] = MUSCL_FM (s_In[sj][time_level][DENS], _dh * GRADB1 (s_In[sj][time_level][PHAS], si), si, dh, dt); 
#           elif ( HYBRID_SCHEME == HYBRID_PPM ) 
//             access rho_L[i, i-1], rho_R[i, i-1], v_L[i]
               s_Fm[sj][si] = PPM_FM   (s_In[sj][time_level][DENS], rho_L, rho_R, v_L, si, dh, dt);
#           endif

#           ifdef CONSERVE_MASS
//             update density fluxes
               if ( time_level == 0 ) s_Flux[sj][si] = 0;
               s_Flux[sj][si] += FLUX_COEFFS[time_level] * s_Fm[sj][si];
#           endif            
         }

#        ifdef __CUDACC__ 
         __syncthreads();
#        endif 


         // Compute density and phase changes
         CELL_LOOP(FLU_NXT, g2, g2)
         {
            fp  = s_Fm[sj][si+1];
            fm  = s_Fm[sj][si  ];
            qp  = QUANTUM_PRESSURE  (s_Sr[sj], si);
            vp  = GRADF3(s_In[sj][time_level][PHAS], si);
            vm  = GRADB3(s_In[sj][time_level][PHAS], si);
            osf = OSHER_SETHIAN_FLUX(vp, vm); 

            De_New = TIME_COEFFS[time_level] * Coeff1 * (fm - fp);
            Ph_New = TIME_COEFFS[time_level] * Coeff2 * (qp - osf);
            for (int tl = 0; tl < time_level + 1; ++tl) {
               De_New += RK_COEFFS[time_level][tl] * s_In[sj][tl][DENS][si];
               Ph_New += RK_COEFFS[time_level][tl] * s_In[sj][tl][PHAS][si];
            }

//          Write density and phase change as well as density fluxes after RK1 update to buffer
            if ( time_level == 0 ) {
               s_In[sj][N_TIME_LEVELS][DENS][si] = FMAX(De_New, FluidMinDens);
               s_In[sj][N_TIME_LEVELS][PHAS][si] =      Ph_New;   
            }

//          While computing the temporary results in RK algorithm, just write them to s_In
            if ( time_level < N_TIME_LEVELS - 1 ) {
               s_In[sj][time_level+1][DENS][si] = De_New;
               s_In[sj][time_level+1][PHAS][si] = Ph_New;    
            } 
//          Write back final results to g_Fluid_In[0] or g_Fluid_Out to save memory
            else if ( time_level == N_TIME_LEVELS - 1 ) {

//             handle the case that we have negative densities or if the velocity timestep criterion is not met -> switch to RK1
               if ( De_New < 0 || De_New != De_New || Ph_New != Ph_New ) {             
                  De_New = s_In[sj][N_TIME_LEVELS][DENS][si];
                  Ph_New = s_In[sj][N_TIME_LEVELS][PHAS][si];
               }
               
//             5.1 data
               if ( FinalOut )
               {
//                apply the the minimum density check
//                --> to be consistent with the CPU solver, we apply it just before storing the output results to g_Fluid_Out
                  De_New = (De_New < MinDens) ? MinDens : De_New;

                  switch ( XYZ )
                  {
                     case 0:  Idx2 = to1D2( k, j, i );    break;
                     case 3:  Idx2 = to1D2( k, i, j );    break;
                     case 6:  Idx2 = to1D2( i, k, j );    break;
                  }

                  g_Fluid_Out[bx][DENS][Idx2] = De_New;
                  g_Fluid_Out[bx][PHAS][Idx2] = Ph_New;
                  g_Fluid_Out[bx][STUB][Idx2] = 0;
               }
               else
               {
                  g_Fluid_In[bx][DENS][Idx1] = De_New;
                  g_Fluid_In[bx][PHAS][Idx1] = Ph_New;
               }

#              ifdef CONSERVE_MASS
//             5.2 fluxes (for the flux-correction operation)
               if ( StoreFlux  &&  tx == 0 )
               if ( k >= FLU_GHOST_SIZE  &&  k < FLU_NXT-FLU_GHOST_SIZE )
               if ( j >= FLU_GHOST_SIZE  &&  j < FLU_NXT-FLU_GHOST_SIZE )
               {
                  Idx3 = __umul24( k-FLU_GHOST_SIZE, PS2 ) + (j-FLU_GHOST_SIZE);

                  g_Flux[bx][XYZ+0][0][Idx3] = s_Flux[ty][  0 + FLU_GHOST_SIZE] / Eta;
                  g_Flux[bx][XYZ+1][0][Idx3] = s_Flux[ty][PS1 + FLU_GHOST_SIZE] / Eta;
                  g_Flux[bx][XYZ+2][0][Idx3] = s_Flux[ty][PS2 + FLU_GHOST_SIZE] / Eta;
               }
#              endif 

//             5.3 reset the target array indices
               j += NColumnOnce;

               if ( j >= j_end )
               {
                  delta_k  = ( j - j_end )/size_j + 1;
                  k       += delta_k;
                  j       -= __umul24( size_j, delta_k );
               }
            }              
         }
#        ifdef __CUDACC__ 
         __syncthreads();
#        endif 
         
      }
      

      Column0     += NColumnOnce;
      NColumnOnce  = MIN( NColumnTotal - Column0, FLU_BLOCK_SIZE_Y );

   } // while ( Column0 < NColumnTotal )

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