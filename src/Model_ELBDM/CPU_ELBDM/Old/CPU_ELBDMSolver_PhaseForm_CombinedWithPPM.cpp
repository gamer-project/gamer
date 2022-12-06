#include "GAMER.h"
#include "CUFLU.h"

//###################################################################################################
#if (  !defined GPU  &&  MODEL == ELBDM  &&  ELBDM_SCHEME == HYBRID )
//###################################################################################################



// useful macros
#define to1D(z,y,x)     ( z*FLU_NXT*FLU_NXT + y*FLU_NXT + x )

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
# if ( HYBRID_SCHEME == UPWIND )

# define UPWIND_FM(Rc, Vb, t)        (   FMAX(Vb, 0) * Rc[t-1] \
                                       + FMIN(Vb, 0) * Rc[t  ] )

# endif // #if ( HYBRID_SCHEME == UPWIND )

//Second-order MUSCL flux reconstruction
#if ( HYBRID_SCHEME == MUSCL )

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
# endif // #if ( HYBRID_SCHEME == MUSCL )

//Third-order PPM flux reconstruction
# if ( HYBRID_SCHEME == TOS )
void PPM_INTERPOLATION(real* a_array, real* a_L_array, real* a_R_array, int i, int densityLimiter);
void PPM_LIMITER      (real* a_array, real* a_L_array, real* a_R_array, int i, int densityLimiter);
real PPM_FM           (real* a_array, real* a_L_array, real* a_R_array, real* v_L_array, int i, real dh, real dt);
# endif // # if ( HYBRID_SCHEME == TOS )

//Second- and fourth-order quantum pressure depending on the scheme used
# if (HYBRID_SCHEME == MUSCL || HYBRID_SCHEME == UPWIND)
# define QUANTUM_PRESSURE(Sr, t)  ((real) 0.5 * (pow(GRADC2(Sr, t), 2) + LAP2(Sr, t)))
# else 
# define QUANTUM_PRESSURE(Sr, t)  ((real) 0.5 * (pow(GRADC4(Sr, t), 2) + LAP4(Sr, t)))
# endif // # if (HYBRID_SCHEME == MUSCL || HYBRID_SCHEME == UPWIND)

//Osher-Sethian flux for Hamilton-Jacobi equation
# define OSHER_SETHIAN_FLUX(vp, vm) ((real) 0.5 * (pow(MIN(vp, 0), 2) + pow(MAX(vm, 0), 2)))


static void CPU_AdvanceX( real u[][ CUBE(FLU_NXT)], real Flux_Array[][NFLUX_TOTAL][ SQR(PS2) ], 
                          const real dt, const real dh, const real Eta, const bool StoreFlux, const real MinDens, 
                          const int j_gap, const int k_gap, const int Flux_XYZ );
static void TransposeXY( real u[][ CUBE(FLU_NXT)] );
static void TransposeXZ( real u[][ CUBE(FLU_NXT)] );


//-------------------------------------------------------------------------------------------------------
// Function    :  CPU_ELBDMSolver
// Description :  CPU ELBDM kinematic solver based on expanding the propagator to the 3rd order 
//
// Note        :  1. The three-dimensional evolution is achieved by applying x, y, and z operators successively.
//                   Since these operators commute, the order of applying them are irrelevant.
//                   --> Input pamameter "XYZ" is actually meaningless (if CONSERVE_MASS is off)
//                   --> Nevertheless, the symmetry in different directions will be broken if CONSERVE_MASS is on
//                2. The implementation is very similar to the function "CPU_FluidSolver_RTVD"
//
// Parameter   :  Flu_Array_In   : Array storing the input variables (only REAL/IMAG)
//                Flu_Array_Out  : Array to store the output variables (DENS/REAL/IMAG)
//                Flux_Array     : Array to store the output flux
//                NPatchGroup    : Number of patch groups to be evaluated
//                dt             : Time interval to advance solution
//                dh             : Grid size
//                Eta            : Particle mass / Planck constant
//                StoreFlux      : true --> store the coarse-fine fluxes
//                                      --> useful only if CONSERVE_MASS is defined
//                Taylor3_Coeff  : Coefficient in front of the third term in the Taylor expansion
//                XYZ            : true  : x->y->z ( forward sweep)
//                                 false : z->y->x (backward sweep)
//                                 --> Meaningless if CONSERVE_MASS is off since the operators along different directions 
//                                     commute
//                                 --> Meaningful if CONSERVE_MASS is on, in which the symmetry along different directions
//                                     are broken ...
//-------------------------------------------------------------------------------------------------------
void CPU_ELBDMSolver_PhaseForm( real Flu_Array_In [][FLU_NIN ][ CUBE(FLU_NXT)], 
                      real Flu_Array_Out[][FLU_NOUT][ SQR(PS2)*PS2 ], 
                      real Flux_Array[][9][NFLUX_TOTAL][ SQR(PS2) ], 
                      const int NPatchGroup, const real dt, const real dh, const real Eta, const bool StoreFlux,
                      const bool XYZ, const real MinDens )
{

   const real FluidMinDens = FMAX(1e-10, MinDens); 

   if ( XYZ )
   {
#     pragma omp parallel for schedule( runtime )
      for (int P=0; P<NPatchGroup; P++)
      {
         CPU_AdvanceX( Flu_Array_In[P], Flux_Array[P], dt, dh, Eta, StoreFlux, FluidMinDens,
                                    0,              0, 0 );                    
                                                                               
         TransposeXY ( Flu_Array_In[P] );                                      
                                                                               
         CPU_AdvanceX( Flu_Array_In[P], Flux_Array[P], dt, dh, Eta, StoreFlux, FluidMinDens,
                       FLU_GHOST_SIZE,              0, 3 );                    
                                                                               
         TransposeXZ ( Flu_Array_In[P] );                                      
                                                                               
         CPU_AdvanceX( Flu_Array_In[P], Flux_Array[P], dt, dh, Eta, StoreFlux, FluidMinDens,
                       FLU_GHOST_SIZE, FLU_GHOST_SIZE, 6 );
   
         TransposeXZ ( Flu_Array_In[P] );
         TransposeXY ( Flu_Array_In[P] );
      }
   }

   else
   {
#     pragma omp parallel for schedule( runtime )
      for (int P=0; P<NPatchGroup; P++)
      {
         TransposeXY ( Flu_Array_In[P] );
         TransposeXZ ( Flu_Array_In[P] );

         CPU_AdvanceX( Flu_Array_In[P], Flux_Array[P], dt, dh, Eta, StoreFlux, FluidMinDens,
                                    0,              0, 6 );         
                                                                    
         TransposeXZ ( Flu_Array_In[P] );                           
                                                                    
         CPU_AdvanceX( Flu_Array_In[P], Flux_Array[P], dt, dh, Eta, StoreFlux, FluidMinDens,
                                    0, FLU_GHOST_SIZE, 3 );         
                                                                    
         TransposeXY ( Flu_Array_In[P] );                           
                                                                    
         CPU_AdvanceX( Flu_Array_In[P], Flux_Array[P], dt, dh, Eta, StoreFlux, FluidMinDens,
                       FLU_GHOST_SIZE, FLU_GHOST_SIZE, 0 );
      }
   }


// copy the updated data to Flu_Array_Out
   int Idx1, Idx2;

#  pragma omp parallel for private( Idx1, Idx2) schedule( runtime )
   for (int P=0; P<NPatchGroup; P++)
   {
//    copy data
      for (int v=0; v<FLU_NOUT; v++) 
      {
         Idx1 = 0;

         for (int k=FLU_GHOST_SIZE; k<FLU_GHOST_SIZE+PS2; k++)
         for (int j=FLU_GHOST_SIZE; j<FLU_GHOST_SIZE+PS2; j++)
         for (int i=FLU_GHOST_SIZE; i<FLU_GHOST_SIZE+PS2; i++)
         {
            Idx2 = to1D(k,j,i);

            Flu_Array_Out[P][v][ Idx1++ ] = Flu_Array_In[P][v][Idx2];
         }
      }
   } // for (int P=0; P<NPatchGroup; P++)

} // FUNCTION : CPU_ELBDMSolver


//#define N_TIME_LEVELS 1
//const real TIME_COEFFS[N_TIME_LEVELS]                = {1.0};
//const real RK_COEFFS  [N_TIME_LEVELS][N_TIME_LEVELS] = {{1.0}};


//#define N_TIME_LEVELS 2
//const real TIME_COEFFS[N_TIME_LEVELS]                = {1.0/2.0, 1.0};
//const real RK_COEFFS  [N_TIME_LEVELS][N_TIME_LEVELS] = {{1.0, 0.0}, {1.0, 0.0}};

#define N_TIME_LEVELS 3
const double TIME_COEFFS[N_TIME_LEVELS]                = {1.0, 1.0/4.0, 2.0/3.0};
const double FLUX_COEFFS[N_TIME_LEVELS]                = {1.0/6.0, 1.0/6.0, 2.0/3.0};
const double RK_COEFFS  [N_TIME_LEVELS][N_TIME_LEVELS] = {{1.0, 0.0, 0.0}, {3.0/4.0, 1.0/4.0, 0.0}, {1.0/3.0, 0.0, 2.0/3.0}};

#define NO_LIMITER    0
#define SOME_LIMITER  1
#define FULL_LIMITER  2 


//-------------------------------------------------------------------------------------------------------
// Function    :  CPU_AdvanceX
// Description :  Use CPU to advance a single patch group by one time-step in the x direction
//
// Note        :  Based on expanding the kinematic propagator to the 3rd order
//
// Parameter   :  u              : Array storing the input variables (only REAL/IMAG)
//                Flux_Array     : Array to store the output flux (only density)
//                dt             : Time interval to advance solution
//                dh             : Grid size
//                Eta            : Particle mass / Planck constant
//                StoreFlux      : true --> store the coarse-fine fluxes
//                                      --> useful only if CONSERVE_MASS is defined
//                Taylor3_Coeff  : Coefficient in front of the third term in the Taylor expansion
//                j_gap          : Number of cells to be skipped on each side in the y direction
//                k_gap          : Number of cells to be skipped on each side in the z direction
//                Flux_XYZ       : Parameter used to determine the place to store the output fluxes
//                                 --> (0,3,6) <-> (x/y/z) fluxes
//                                 --> useful only if CONSERVE_MASS is defined
//-------------------------------------------------------------------------------------------------------
void CPU_AdvanceX( real u[][ CUBE(FLU_NXT)], real Flux_Array[][NFLUX_TOTAL][ SQR(PS2) ], 
                   const real dt, const real dh, const real Eta, const bool StoreFlux, const real MinDens, 
                   const int j_gap, const int k_gap, const int Flux_XYZ )
{
   const real Coeff1 = dt/(dh * Eta);
   const real Coeff2 = dt/(dh*dh * Eta);
   const real _dh    = 1./dh;
   const int j_start = j_gap;
   const int k_start = k_gap;
   const int j_end   = FLU_NXT - j_gap;
   const int k_end   = FLU_NXT - k_gap;
   int Idx;

   real vp, vm, fp, fm, qp, osf;

   real Rc[N_TIME_LEVELS][FLU_NXT];  // one column of the density in the input array "u" at all time levels
   real Pc[N_TIME_LEVELS][FLU_NXT];  // one column of the phase   in the input array "u" at all time levels
   real *Rc_target = NULL;           // pointer to set density array in update loop
   real *Pc_target = NULL;           // pointer to set phase   array in update loop
   real *Rc_N = NULL;                // pointer to store the full-step density
   real *Pc_N = NULL;                // pointer to store the full-step phase  

   real Fm[FLU_NXT], Sr[FLU_NXT], Fm_avg[FLU_NXT], Rc_RK1[FLU_NXT], Pc_RK1[FLU_NXT];
  
#  if ( HYBRID_SCHEME == MUSCL || HYBRID_SCHEME == UPWIND)
   const int ghostZonePerStage = 2; 
#  else
   const int ghostZonePerStage = 4; 
   real v_C[FLU_NXT], v_L[FLU_NXT], v_R[FLU_NXT], rho_L[FLU_NXT], rho_R[FLU_NXT];
#  endif 


#  ifdef CONSERVE_MASS
   uint Idx3;
#  endif

// loop over all targeted columns
   for (int k=k_start; k<k_end; k++)
   for (int j=j_start; j<j_end; j++)
   {

//    1. backup one column of data
//    ------------------------------------------------------------------------------------------------------------
      Idx    = to1D(k,j,0);

//    set pointers to output arrays
      Rc_N = &u[DENS][Idx];
      Pc_N = &u[PHAS][Idx];
      
      memcpy( Rc[0], Rc_N, FLU_NXT*sizeof(real) );
      memcpy( Pc[0], Pc_N, FLU_NXT*sizeof(real) );
      memset( Fm_avg,   0, FLU_NXT*sizeof(real) );

#     ifdef GAMER_DEBUG
      if ( Rc[0][i] != Rc[0][i] || Pc[0][i] != Pc[0][i] || Rc[0][i] < 0)
          Aux_Error( ERROR_INFO, "nan in input array of fluid solver k %d j %d i %d Rho %f Pha %f !\n\n", k, j, i, Rc[0][i], Pc[0][i]); 
#     endif


      for (int time_level = 0; time_level < N_TIME_LEVELS; ++time_level) 
      {
         int g1 = ghostZonePerStage *   time_level       ;
         int g2 = ghostZonePerStage * ( time_level + 1 ) ;

//       exclude ghost zones from previous stages from computation of quantum pressure
         for (int i = g1; i < FLU_NXT - g1; i++) {
            Sr[i] = real(0.5) * log(FMAX(Rc[time_level][i], MinDens));
         }

//       prepare arrays of velocity and density interpolated to cell-face for PPM scheme
#        if ( HYBRID_SCHEME == TOS )

//       GRADC4 requires ghost boundary of size 2 in each direction
//       Respect ghost boundary of size 2 from GRADC4 since PPM_INTERPOLATION accesses v_C[i-1, i, i+1, i+2]
//       PPM interpolation fills v_R from [i_begin to i_end] v_L from [i_begin + 1 to i_end + 1]
//       Therefore, it fills v_L[g1 + 4 to FLU_NXT - g1 - 3 ] and v_R[g1 + 3 to FLU_NXT - g1 - 4 ]
         for (int i = g1 + 2; i < FLU_NXT - g1 - 2; i++) {
            v_C[i] = _dh * GRADC4(Pc[time_level], i);
         }
         for (int i = g1 - 3; i < FLU_NXT - g1 - 4; i++) {
            PPM_INTERPOLATION(v_C, v_L, v_R, i, NO_LIMITER);
         }

         //Fill rho_L[g2 - 1] as PPM_INTERPOATION only fills a_L[i+1]
         PPM_INTERPOLATION(Rc[time_level], rho_L, rho_R, g2 - 2, FULL_LIMITER);

         for (int i = g2 - 1; i < FLU_NXT - g2 + 2; i++) {
            PPM_INTERPOLATION(Rc[time_level], rho_L, rho_R, i, FULL_LIMITER);
            PPM_LIMITER      (Rc[time_level], rho_L, rho_R, i, FULL_LIMITER);
         } 
#        endif // #        if ( HYBRID_SCHEME == TOS )

//       compute backward density fluxes at all cell faces of real cells
         for (int i = g2; i < FLU_NXT - g2 + 1; i++)
         {
#        if ( HYBRID_SCHEME == UPWIND )
//          Access Rc[time_level][i, i-1], Pc[time_level][i, i-1]
            Fm[i] = UPWIND_FM(Rc[time_level], _dh * GRADB1 (Pc[time_level], i), i); 
#        elif ( HYBRID_SCHEME == MUSCL )
//          Access Rc[time_level][i, i-1, i-2], Pc[time_level][i, i-1]
            Fm[i] = MUSCL_FM (Rc[time_level], _dh * GRADB1 (Pc[time_level], i), i, dh, dt); 
#        elif ( HYBRID_SCHEME == TOS ) 
//          Access rho_L[i, i-1], rho_R[i, i-1], v_L[i]
            Fm[i] = PPM_FM   (Rc[time_level], rho_L, rho_R, v_L, i, dh, dt);
#        endif

#        ifdef CONSERVE_MASS
            Fm_avg[i] += FLUX_COEFFS[time_level] * Fm[i];
#        endif 
         }

//       update density and phase of real cells
         for (int i = g2; i < FLU_NXT - g2; i++)
         {
            fp  = Fm[i+1];
            fm  = Fm[i  ];
            qp  = QUANTUM_PRESSURE  (Sr, i);
            vp  = GRADF3(Pc[time_level], i);
            vm  = GRADB3(Pc[time_level], i);
            osf = OSHER_SETHIAN_FLUX(vp, vm); 

            if (time_level + 1 < N_TIME_LEVELS)
            {
               Rc_target = &Rc[time_level + 1][i];
               Pc_target = &Pc[time_level + 1][i];
            }
            else 
            {
               Rc_target = &Rc_N[i];
               Pc_target = &Pc_N[i];
            }

            *Rc_target = TIME_COEFFS[time_level] * Coeff1 * (fm - fp);
            *Pc_target = TIME_COEFFS[time_level] * Coeff2 * (qp - osf);
            for (int tx = 0; tx < time_level + 1; ++tx) {
               *Rc_target += RK_COEFFS[time_level][tx] * Rc[tx][i];
               *Pc_target += RK_COEFFS[time_level][tx] * Pc[tx][i];
            }

//          overwrite nans and negative densities with RK1 solution
            if ( time_level + 1 ==  N_TIME_LEVELS ) {
               if ( *Rc_target < 0 || *Rc_target  != *Rc_target || *Pc_target != *Pc_target ) {             
                  *Rc_target = FMIN(Rc_RK1[i], MinDens);
                  *Pc_target =      Pc_RK1[i];               
               }
#              ifdef GAMER_DEBUG
               if (*Rc_target != *Rc_target || *Pc_target != *Pc_target) {
                  Aux_Error( ERROR_INFO, "nan for dt = %f input dens %f input phas %f rk1 dens = %f rk1 phas = %f!\n\n", dt, Rc[0][i], Pc[0][i], Rc_N[i], Pc_N[i]); 
               }
#              endif
            }
         }


//       store RK1 solutions until end of RKN iterations
         if (time_level == 0) {
            memcpy( Rc_RK1, Rc_target, FLU_NXT*sizeof(real));
            memcpy( Pc_RK1, Pc_target, FLU_NXT*sizeof(real));
         }
      }

//    4 save the fluxes across all patch boundaries
//    ------------------------------------------------------------------------------------------------------------
#     ifdef CONSERVE_MASS
      if ( StoreFlux )
      if (  ( j>=FLU_GHOST_SIZE && j<FLU_NXT-FLU_GHOST_SIZE )  &&  ( k>=FLU_GHOST_SIZE && k<FLU_NXT-FLU_GHOST_SIZE )  )
      {
         Idx3 = (k-FLU_GHOST_SIZE)*PS2 + (j-FLU_GHOST_SIZE);

         Flux_Array[Flux_XYZ+0][0][Idx3] = Fm_avg[      FLU_GHOST_SIZE] / Eta;
         Flux_Array[Flux_XYZ+1][0][Idx3] = Fm_avg[PS1 + FLU_GHOST_SIZE] / Eta;
         Flux_Array[Flux_XYZ+2][0][Idx3] = Fm_avg[PS2 + FLU_GHOST_SIZE] / Eta;
      }
#     endif // #ifdef CONSERVE_MASS

   } // for j,k





} // FUNCTION : CPU_AdvanceX




//-------------------------------------------------------------------------------------------------------
// Function    :  TrasposeXY
// Description :  Transpose the x and y directions
//
// Parameter   :  u : Input wave function (density, real, imaginary)
//-------------------------------------------------------------------------------------------------------
void TransposeXY( real u[][ CUBE(FLU_NXT)] )
{

   real (*u_xy)[FLU_NXT*FLU_NXT] = new real [FLU_NIN][FLU_NXT*FLU_NXT];
   int Idx1, Idx2;

   for (int k=0; k<FLU_NXT; k++)
   {
      for (int j=0; j<FLU_NXT; j++)
      for (int i=0; i<FLU_NXT; i++)
      {
         Idx1 = to1D(k,j,i);
         Idx2 = j + i*FLU_NXT;

         u_xy[0][Idx2] = u[0][Idx1];
         u_xy[1][Idx2] = u[1][Idx1];
      }

      for (int v=0; v<FLU_NIN; v++)    memcpy( &u[v][to1D(k,0,0)], u_xy[v], FLU_NXT*FLU_NXT*sizeof(real) );
   }

   delete [] u_xy;

} // FUNCTION : TrasposeXY



//-------------------------------------------------------------------------------------------------------
// Function    :  TrasposeXZ
// Description :  Transpose the x and z directions
//
// Parameter   :  u : Input wave function (density, real, imaginary)
//-------------------------------------------------------------------------------------------------------
void TransposeXZ( real u[][ CUBE(FLU_NXT)] )
{

   real u_temp[FLU_NIN];
   int Idx1, Idx2;

   for (int j=0; j<FLU_NXT; j++)
   for (int k=0; k<FLU_NXT; k++)
   {
      for (int i=0; i<k; i++)
      {
         Idx1 = to1D(k,j,i);
         Idx2 = to1D(i,j,k);
         
         u_temp[0] = u[0][Idx1];
         u_temp[1] = u[1][Idx1];
         
         u[0][Idx1] = u[0][Idx2];
         u[1][Idx1] = u[1][Idx2];
         
         u[0][Idx2] = u_temp[0];
         u[1][Idx2] = u_temp[1];
      }

      Idx1 = to1D(k,j,k);
   } // j,k

} // FUNCTION : TrasposeXZ


# if ( HYBRID_SCHEME == TOS )


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
#endif



#endif // #if ( MODEL == ELBDM  &&  ELBDM_SCHEME == HYBRID )
