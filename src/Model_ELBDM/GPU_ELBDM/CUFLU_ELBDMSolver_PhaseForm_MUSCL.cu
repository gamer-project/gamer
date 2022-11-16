#include "Macro.h"
#include "CUFLU.h"

#if ( defined GPU  &&  MODEL == ELBDM && ELBDM_SCHEME == HYBRID )



// useful macros
#define to1D1(z,y,x) ( __umul24(z, FLU_NXT*FLU_NXT) + __umul24(y, FLU_NXT) + x )
#define to1D2(z,y,x) ( __umul24(z-FLU_GHOST_SIZE, PS2*PS2) + __umul24(y-FLU_GHOST_SIZE, PS2) + x-FLU_GHOST_SIZE )

# define GTR( a, b )     (  ( (a) > (b) ) ? (1) : (0)  )
# define LSS( a, b )     (  ( (a) < (b) ) ? (1) : (0)  )
            
# define CENTERED_GRADIENT(In, t) ( real(1.0/2.0 ) * (   In[t + 1] - In[t - 1] ) )
# define BACKWARD_GRADIENT(In, t) (                      In[t    ] - In[t - 1]   )
# define FORWARD_GRADIENT(In, t)  (                      In[t + 1] - In[t    ]   )
# define LAPLACIAN(In, t)         ( In[t - 1] - (real)2.0*In[t] + In[t + 1] )

// VAN ALBADA LIMITER
# define LIMITER(In) ( (SQR(In) + In)/((real)1. + SQR(In)) )
// MC
//# define LIMITER(r) MAX(0, MIN(MIN((1. + r) / 2., 2.), 2. * r))
// VAN LEER
//# define LIMITER(r) ((r + FABS(r))/(1.0 + FABS(r)))

# define UPWIND_GRADIENT_RATIO(Rc, Vb, t) ( ((Rc[t-1] - Rc[t-2]) * GTR(Vb, 0)  \
                                           + (Rc[t+1] - Rc[t  ]) * LSS(Vb, 0)) \
                                           / (Rc[t  ] - Rc[t-1] + (((Rc[t] - Rc[t-1]) == 0) ? 1e-8 : 0)))

//Second-order MUSCL flux reconstruction
# define MUSCL_FLUX(Rc, Vb, t, dx, dt) ( FMAX(Vb, 0) * Rc[t-1] \
                                       + FMIN(Vb, 0) * Rc[t  ] \
                                       + real(0.5) * FABS(Vb) * (1. - FABS(Vb * dt/dx)) * LIMITER(UPWIND_GRADIENT_RATIO(Rc, Vb, t)) * (Rc[t] - Rc[t - 1]) )

//Ratio of subsequent backward gradients
# define BACKWARD_GRADIENT_RATIO(Pc, t) ((Pc[t] - Pc[t-1]) / (Pc[t-1] - Pc[t-2] + (((Pc[t-1] - Pc[t-2]) == 0) ?  1e-8 : 0)))

#define N_TIME_LEVELS 3
static const __device__ real TIME_COEFFS[N_TIME_LEVELS] = {1., 1./4, 2./3};
static const __device__ real RK_COEFFS [N_TIME_LEVELS][N_TIME_LEVELS] = {{1., 0., 0.}, {3./4, 1./4, 0.}, {1./3, 0, 2./3}};

#ifdef CONSERVE_MASS
static const __device__ real FLUX_COEFFS[N_TIME_LEVELS] = {1./6, 1./6, 2./3};
#endif 

static __device__ void CUFLU_Advance( real g_Fluid_In [][FLU_NIN ][ CUBE(FLU_NXT) ],
                                      real g_Fluid_Out[][FLU_NOUT][ CUBE(PS2) ],
                                      real g_Flux     [][9][NFLUX_TOTAL][ SQR(PS2) ],
                                      const real dt, const real _dh, const real Eta, const bool StoreFlux,
                                      const uint j_gap, const uint k_gap,         
                                      real s_In[][FLU_NIN][FLU_BLOCK_SIZE_Y][FLU_NXT],
                                      real s_Ql[][FLU_NXT],
                                      real s_LogRho[][FLU_NXT], 
                                      real s_Fm[][FLU_NXT],
                                      bool s_RK1[][FLU_NXT],
                                      real s_Flux[][PS2+1], const bool FinalOut, const int XYZ, const real MinDens );




//-------------------------------------------------------------------------------------------------------
// Function    :  CUFLU_ELBDMSolver_PhaseForm_MUSCL
// Description :  GPU solver for kinetic term in Hamilton Jacobi-Madelung equations
//
// Note        :  1. The three-dimensional evolution is achieved by applying x, y, and z operators successively.
//                   Since these operators commute, the order of applying them are irrelevant.
//                   --> Input parameter "XYZ" is actually useless
//                2. The implementation is very similar to the function " CUFLU_FluidSolver_RTVD"
//                4. Prefix "g" for pointers pointing to the "Global" memory space
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



__global__ void CUFLU_ELBDMSolver_PhaseForm_MUSCL( real g_Fluid_In [][FLU_NIN ][ CUBE(FLU_NXT) ],
                                   real g_Fluid_Out[][FLU_NOUT][ CUBE(PS2) ],
                                   real g_Flux     [][9][NFLUX_TOTAL][ SQR(PS2) ],
                                   const real dt, const real _dh, const real Eta, const bool StoreFlux,
                                   const bool XYZ, const real MinDens )
{

   __shared__ real s_In  [N_TIME_LEVELS+1][FLU_NIN][FLU_BLOCK_SIZE_Y][FLU_NXT];

   __shared__ real s_Ql                          [FLU_BLOCK_SIZE_Y][FLU_NXT]; // one column of the gradient ratios for phase for every thread block
   __shared__ real s_LogRho                      [FLU_BLOCK_SIZE_Y][FLU_NXT]; // one column of the log(rho) for every thread block
   __shared__ real s_Fm                          [FLU_BLOCK_SIZE_Y][FLU_NXT]; // one column of the fluxes for every thread block
   __shared__ bool s_RK1                         [FLU_BLOCK_SIZE_Y][FLU_NXT]; // one column of the flags with control RK1 vs RK3 for every thread block

#  ifdef CONSERVE_MASS
   __shared__ real s_Flux[FLU_BLOCK_SIZE_Y][PS2+1];
#  else
   real (*s_Flux)[PS2+1]                     = NULL;  // useless if CONSERVE_MASS is off
#  endif

   if ( XYZ )
   {
      CUFLU_Advance( g_Fluid_In, g_Fluid_Out, g_Flux, dt, _dh, Eta, StoreFlux,
                                  0,              0, s_In, s_Ql, s_LogRho, s_Fm, s_RK1, s_Flux, false, 0, MinDens );
      CUFLU_Advance( g_Fluid_In, g_Fluid_Out, g_Flux, dt, _dh, Eta, StoreFlux,
                     FLU_GHOST_SIZE,              0, s_In, s_Ql, s_LogRho, s_Fm, s_RK1, s_Flux, false, 3, MinDens );
      CUFLU_Advance( g_Fluid_In, g_Fluid_Out, g_Flux, dt, _dh, Eta, StoreFlux,
                     FLU_GHOST_SIZE, FLU_GHOST_SIZE, s_In, s_Ql, s_LogRho, s_Fm, s_RK1, s_Flux,  true, 6, MinDens );
   }

   else
   {
      CUFLU_Advance( g_Fluid_In, g_Fluid_Out, g_Flux, dt, _dh, Eta, StoreFlux,
                                  0,              0, s_In, s_Ql, s_LogRho, s_Fm, s_RK1, s_Flux, false, 6, MinDens );
      CUFLU_Advance( g_Fluid_In, g_Fluid_Out, g_Flux, dt, _dh, Eta, StoreFlux,
                                  0, FLU_GHOST_SIZE, s_In, s_Ql, s_LogRho, s_Fm, s_RK1, s_Flux, false, 3, MinDens );
      CUFLU_Advance( g_Fluid_In, g_Fluid_Out, g_Flux, dt, _dh, Eta, StoreFlux,
                     FLU_GHOST_SIZE, FLU_GHOST_SIZE, s_In, s_Ql, s_LogRho, s_Fm,s_RK1,  s_Flux,  true, 0, MinDens );
   }

} // FUNCTION : CUFLU_ELBDMSolver_PhaseForm_MUSCL



//-------------------------------------------------------------------------------------------------------
// Function    :  CUFLU_Advance
// Description :  Use GPU to advance solutions by one time-step
//
// Note        :  1. Based on expanding the kinematic propagator to 3rd order
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
//                s_Ql           : Shared memory array to store the slope limiters
//                s_Vm           : Shared memory array to store the backward velocities
//                s_LogRho       : Shared memory array to store the density logarithms
//                s_Fm           : Shared memory array to store the boundary fluxes during the computation
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
                               real s_In    [][FLU_NIN][FLU_BLOCK_SIZE_Y][FLU_NXT],
                               real s_Ql    [][FLU_NXT],
                               real s_LogRho[][FLU_NXT], 
                               real s_Fm    [][FLU_NXT],
                               bool s_RK1   [][FLU_NXT],
                               real s_Flux  [][PS2+1], const bool FinalOut,
                               const int XYZ, const real MinDens )
{

   const real  dh      = 1./_dh; 
   const real _dh2     = SQR(_dh);

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

   uint   Idx1, Idx2, Idx3, delta_k;

   const uint NThread     = FLU_BLOCK_SIZE_X*FLU_BLOCK_SIZE_Y;

   //slope-limiter
   real ql, qc, Qc, Qr;
   //velocities dS/dx = v and density fluxes f at i - 1/2, i + 1/2 
   real vm, vp;

   //minimum timestep allowed by velocity-dependent CFL criterion
   real dt_min; 

   //change of density and phase in time step
   real De_New, Ph_New;
   real ddensity, dphase;

   uint   Idx;
   uint   si, sj;                                           // array indices used in the shared memory array

#  ifdef CONSERVE_MASS
   uint   f;                                           // array indices used in the s_Flux array
#  endif 

   uint NStep;

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

//       1.2 load the interior data into shared memory at time_level 0
         s_In[0][0][ty][i] = g_Fluid_In[bx][0][Idx1];
         s_In[0][1][ty][i] = g_Fluid_In[bx][1][Idx1];
         s_RK1[ty][i]      = false; 

//       1.3 load the ghost-zone data into shared memory
         if ( LoadGhost )
         {
            s_In[0][0][ty][LoadGhost_i] = g_Fluid_In[bx][0][ (int)Idx1 + LoadGhost_dIdx1 ];
            s_In[0][1][ty][LoadGhost_i] = g_Fluid_In[bx][1][ (int)Idx1 + LoadGhost_dIdx1 ];
            s_RK1 [ty][LoadGhost_i]     = false; 
         }
      } // if ( tid < NColumnOnce*PS2 )

      __syncthreads();


      // RK iterations

      for (int time_level = 0; time_level < N_TIME_LEVELS; ++time_level) 
      {

         // Compute second-order backward velocities and backward density fluxes
         Idx = tid;
         // Backward fluxes need a larger ghost zone at the beginning of the array
         NStep = FLU_NXT - 2*(time_level + 1) - 2*time_level;

         while ( Idx < NColumnOnce * NStep )
         {
            si = Idx % NStep + 2*(time_level + 1);
            sj = Idx / NStep;

            s_Ql  [sj][si]  =       BACKWARD_GRADIENT_RATIO  (s_In[time_level][1][sj], si);
            vm              = _dh * BACKWARD_GRADIENT        (s_In[time_level][1][sj], si);
            s_Fm  [sj][si]  =       MUSCL_FLUX               (s_In[time_level][0][sj], vm, si, dh, dt);    

//          check the velocity-dependent CFL-condition and switch to forward-Euler for updating the density wherever the CFL-condition is not met
//          dt = 1 / MaxdS_dx * 0.5 * ELBDM_ETA * DT__VELOCITY;
//          compute CFL condition timestep
            dt_min = dh / FABS(vm) * 0.5 * Eta * 2.0 ;

//          if the time step adopted in solver is larger than what velocity-dependent CFL condition allows, we switch to first-order RK
            if ( dt > dt_min ) {
//             compute how far wrong information can propagate to determine where we need to switch to forward Euler
               int l_min = i - 2;
               int l_max = i + 2 + 1;
               if (l_min < 0)        l_min = 0;
               if (l_max > FLU_NXT ) l_max = FLU_NXT; 
               for (int l = l_min; l < l_max; ++l) s_RK1[sj][l] = true;  
            }
            Idx += NThread;
         } // while ( Idx < NColumnOnce*NStep )
         

         // Compute density logarithms
         Idx = tid;
         NStep = FLU_NXT - 4*time_level;

         while ( Idx < NColumnOnce * NStep )
         {
            si = Idx % NStep + 2*time_level;
            sj = Idx / NStep;

            //Compute density logarithms
            s_LogRho[sj][si] = log(FMAX(s_In[time_level][0][sj][si], MinDens));

            Idx += NThread;
         } // while ( Idx < NColumnOnce*NStep )

         __syncthreads();

#        ifdef CONSERVE_MASS
         // Update density fluxes
         Idx = tid;
         while ( Idx < NColumnOnce*(PS2+1) )
         {
            si  = Idx % (PS2+1);
            sj  = Idx / (PS2+1);
            f   = si + FLU_GHOST_SIZE;
            if ( time_level == 0 ) s_Flux[sj][si] = 0;
            s_Flux[sj][si] += FLUX_COEFFS[time_level] * s_Fm[sj][f];

            Idx += NThread;
         } // while ( Idx < NColumnOnce*(PS2+1) )

         __syncthreads();
#        endif            
         
         // Compute density and phase changes
         Idx = tid;
         NStep = FLU_NXT - 4*(time_level + 1);

         while ( Idx < NColumnOnce * NStep )
         {
            si = Idx % NStep + 2*(time_level + 1);
            sj = Idx / NStep;


            //Slope-limited second order gradients
            ql = s_Ql[sj][si  ];
            qc = s_Ql[sj][si+1];
            Qc = 1./(s_Ql[sj][si+1] + ((s_Ql[sj][si+1] == 0) ? 1e-8 : 0));
            Qr = 1./(s_Ql[sj][si+2] + ((s_Ql[sj][si+2] == 0) ? 1e-8 : 0));

            vm = _dh * BACKWARD_GRADIENT(s_In[time_level][1][sj], si) * ( (real) 1. + (real) 0.5 * LIMITER(qc) - (real) 0.5 * LIMITER(ql)/(ql + ((ql==0) ? 1e-8 : 0)));
            vp = _dh * FORWARD_GRADIENT (s_In[time_level][1][sj], si) * ( (real) 1. + (real) 0.5 * LIMITER(Qc) - (real) 0.5 * LIMITER(Qr)/(Qr + ((Qr==0) ? 1e-8 : 0)));
   
            ddensity = _dh * (s_Fm[sj][si+1] - s_Fm[sj][si]);

            dphase   = (SQR(FMIN(vp, 0)) + SQR(FMAX(vm, 0)))/2;
            dphase  += - _dh2/4 * LAPLACIAN(s_LogRho[sj], si) - _dh2/8 * SQR(CENTERED_GRADIENT(s_LogRho[sj], si));


            De_New = - TIME_COEFFS[time_level] * dt / Eta * ddensity;
            Ph_New = - TIME_COEFFS[time_level] * dt / Eta * dphase;
            for (int tl = 0; tl < time_level + 1; ++tl) {
               De_New += RK_COEFFS[time_level][tl] * s_In[tl][0][sj][si];
               Ph_New += RK_COEFFS[time_level][tl] * s_In[tl][1][sj][si];
            }

//          Write density and phase change as well as density fluxes after RK1 update to buffer
            if ( time_level == 0 ) {
               //s_Fm[1][sj][si] = s_Fm[0][sj][si];
               s_In[N_TIME_LEVELS][0][sj][si] = FMAX(s_In[0][0][sj][si] - dt / Eta * ddensity, MinDens);
               s_In[N_TIME_LEVELS][1][sj][si] =      s_In[0][1][sj][si] - dt / Eta * dphase;   
            }

//          While computing the temporary results in RK algorithm, just write them to s_In
            if ( time_level < N_TIME_LEVELS - 1 ) {
               s_In[time_level + 1][0][sj][si] = De_New;
               s_In[time_level + 1][1][sj][si] = Ph_New;    
               Idx += NThread;
            } 
//          Write back final results to g_Fluid_In[0] or g_Fluid_Out to save memory
            else if ( time_level == N_TIME_LEVELS - 1 ) {


//             handle the case that we have negative densities or if the velocity timestep criterion is not met -> switch to RK1
               if ( De_New < 0 || De_New != De_New ) {             
                  De_New = s_In[N_TIME_LEVELS][0][sj][si];
                  Ph_New = s_In[N_TIME_LEVELS][1][sj][si];
                  //s_Flux[sj][si - FLU_GHOST_SIZE] = s_Fm[1][sj][si];
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

                  g_Fluid_Out[bx][0][Idx2] = De_New;
                  g_Fluid_Out[bx][1][Idx2] = Ph_New;
                  g_Fluid_Out[bx][2][Idx2] = 0;
               }
               else
               {
                  g_Fluid_In[bx][0][Idx1] = De_New;
                  g_Fluid_In[bx][1][Idx1] = Ph_New;
               }

#              ifdef CONSERVE_MASS
//             5.2 fluxes (for the flux-correction operation)
               if ( StoreFlux  &&  tx == 0 )
               if ( k >= FLU_GHOST_SIZE  &&  k < FLU_NXT-FLU_GHOST_SIZE )
               if ( j >= FLU_GHOST_SIZE  &&  j < FLU_NXT-FLU_GHOST_SIZE )
               {
                  Idx3 = __umul24( k-FLU_GHOST_SIZE, PS2 ) + (j-FLU_GHOST_SIZE);

                  g_Flux[bx][XYZ+0][0][Idx3] = s_Flux[ty][  0] / Eta;
                  g_Flux[bx][XYZ+1][0][Idx3] = s_Flux[ty][PS1] / Eta;
                  g_Flux[bx][XYZ+2][0][Idx3] = s_Flux[ty][PS2] / Eta;
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
               Idx += NThread;
            }                  
         }
         __syncthreads();
         
      }
      

      Column0     += NColumnOnce;
      NColumnOnce  = MIN( NColumnTotal - Column0, FLU_BLOCK_SIZE_Y );

   } // while ( Column0 < NColumnTotal )

} // FUNCTION : CUFLU_Advance



#endif // #if ( defined GPU  &&  MODEL == ELBDM && ELBDM_SCHEME == HYBRID)