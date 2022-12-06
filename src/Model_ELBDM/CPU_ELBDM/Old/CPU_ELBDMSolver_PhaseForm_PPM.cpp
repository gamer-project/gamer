#include "GAMER.h"
#include "CUFLU.h"


#if ( !defined GPU  &&  MODEL == ELBDM && ELBDM_SCHEME == HYBRID && HYBRID_SCHEME == TOS )


// useful macros
#define to1D(z,y,x)     ( z*FLU_NXT*FLU_NXT + y*FLU_NXT + x )

# define GTR( a, b )     (  ( (a) > (b) ) ? (1) : (0)  )
# define LSS( a, b )     (  ( (a) < (b) ) ? (1) : (0)  )
            
# define BACKWARD_GRADIENT(In, t) ( In[t    ] - In[t - 1] )
# define FORWARD_GRADIENT( In, t) ( In[t + 1] - In[t    ] )

# define F2_GRADIENT(In, t) ( real(1.0/2.0 ) * ( -3*In[t] + 4*In[t+1] - In[t+2]))
# define B2_GRADIENT(In, t) ( real(1.0/2.0 ) * (  3*In[t] - 4*In[t-1] + In[t-2]))

# define F3_GRADIENT(In, t) ( real(1.0/6.0 ) * ( -2*In[t-1] - 3*In[t] + 6*In[t+1] - In[t+2]))
# define B3_GRADIENT(In, t) ( real(1.0/6.0 ) * (  2*In[t+1] + 3*In[t] - 6*In[t-1] + In[t-2]))

# define LAP1(In, t)   (                            +                 In[t-1] - real(2.0/1.0) * In[t  ] +                 In[t+1]                            )
# define LAP2(In, t)   ( - real(1.0/12.0) * In[t-2] + real(4.0/3.0) * In[t-1] - real(5.0/2.0) * In[t  ] + real(4.0/3.0) * In[t+1] - real(1.0/12.0) * In[t+2] )

# define GRAD1(In, t)  (                            - real(1.0/2.0) * In[t-1]                           + real(1.0/2.0) * In[t+1]                            )
# define GRAD2(In, t)  ( + real(1.0/12.0) * In[t-2] - real(2.0/3.0) * In[t-1]                           + real(2.0/3.0) * In[t+1] - real(1.0/12.0) * In[t+2] )


static void CPU_AdvanceX( real u[][ CUBE(FLU_NXT) ], real Flux_Array[][NFLUX_TOTAL][ SQR(PS2) ],
                          const real dt, const real dh, const real Eta, const bool StoreFlux, const real Taylor3_Coeff, const real MinDens, 
                          const int j_gap, const int k_gap, const int Flux_XYZ );
static void TransposeXY( real u[][ CUBE(FLU_NXT) ] );
static void TransposeXZ( real u[][ CUBE(FLU_NXT) ] );
static void computePPMInterpolation(real* a_array, real* a_L_array, real* a_R_array, bool fix2, bool fix3, int ghostBoundary);



//-------------------------------------------------------------------------------------------------------
// Function    :  CPU_ELBDMSolver_PhaseForm
// Description :  Solve kinetic term in Hamilton Jacobi-Madelung equations
//
// Note        :  1. The three-dimensional evolution is achieved by applying x, y, and z operators successively.
//                   Since these operators commute, the order of applying them are irrelevant.
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
//-------------------------------------------------------------------------------------------------------
void CPU_ELBDMSolver_PhaseForm_PPM( real Flu_Array_In [][FLU_NIN ][ CUBE(FLU_NXT) ],
                      real Flu_Array_Out[][FLU_NOUT][ CUBE(PS2) ],
                      real Flux_Array[][9][NFLUX_TOTAL][ SQR(PS2) ],
                      const int NPatchGroup, const real dt, const real dh, const real Eta, const bool StoreFlux,
                      const real Taylor3_Coeff, const bool XYZ, const real MinDens )
{
   const real FluidMinDens = FMAX(1e-10, MinDens); 

   if ( true )
   {
#     pragma omp parallel for schedule( runtime )
      for (int P=0; P<NPatchGroup; P++)
      {
         CPU_AdvanceX( Flu_Array_In[P], Flux_Array[P], dt, dh, Eta, StoreFlux, Taylor3_Coeff, FluidMinDens,
                                    0,              0, 0 );

         TransposeXY ( Flu_Array_In[P] );

         CPU_AdvanceX( Flu_Array_In[P], Flux_Array[P], dt, dh, Eta, StoreFlux, Taylor3_Coeff, FluidMinDens,
                       FLU_GHOST_SIZE,              0, 3 );

         TransposeXZ ( Flu_Array_In[P] );

         CPU_AdvanceX( Flu_Array_In[P], Flux_Array[P], dt, dh, Eta, StoreFlux, Taylor3_Coeff, FluidMinDens,
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

         CPU_AdvanceX( Flu_Array_In[P], Flux_Array[P], dt, dh, Eta, StoreFlux, Taylor3_Coeff, FluidMinDens,
                                    0,              0, 6 );

         TransposeXZ ( Flu_Array_In[P] );

         CPU_AdvanceX( Flu_Array_In[P], Flux_Array[P], dt, dh, Eta, StoreFlux, Taylor3_Coeff, FluidMinDens,
                                    0, FLU_GHOST_SIZE, 3 );

         TransposeXY ( Flu_Array_In[P] );

         CPU_AdvanceX( Flu_Array_In[P], Flux_Array[P], dt, dh, Eta, StoreFlux, Taylor3_Coeff, FluidMinDens,
                       FLU_GHOST_SIZE, FLU_GHOST_SIZE, 0 );
      }
   }


// copy the updated data to Flu_Array_Out
   int  Idx1, Idx2;
   
#  pragma omp parallel for private( Idx1, Idx2 ) schedule( runtime )
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

            if (v == DENS || v == PHAS)
               Flu_Array_Out[P][v][ Idx1++ ] = Flu_Array_In[P][v][Idx2];
            else 
               Flu_Array_Out[P][v][ Idx1++ ] = 0;
         }
      }

//    evaluate the new density (and apply the minimum density check)
      for (int t=0; t<CUBE(PS2); t++)
      {
         if ( Flu_Array_Out[P][0][t] < MinDens )
            Flu_Array_Out[P][0][t] = MinDens;
      }
   } // for (int P=0; P<NPatchGroup; P++)

} // FUNCTION : CPU_ELBDMSolver_PhaseForm_MUSCL



//-------------------------------------------------------------------------------------------------------
// Function    :  CPU_AdvanceX
// Description :  Use CPU to advance a single patch group by one time-step in the x direction
//
// Note        :  Based on solving the Hamilton-Jacobi-Madelung equations
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
#define N_TIME_LEVELS 3
const real TIME_COEFFS[N_TIME_LEVELS] = {1., 1./4, 2./3};
const real FLUX_COEFFS[N_TIME_LEVELS] = {1./6, 1./6, 2./3};
const real RK_COEFFS  [N_TIME_LEVELS][N_TIME_LEVELS] = {{1., 0., 0.}, {3./4, 1./4, 0.}, {1./3, 0, 2./3}};

//#define N_TIME_LEVELS 2
//const real TIME_COEFFS[N_TIME_LEVELS] = {1/2., 1};
//const real RK_COEFFS [N_TIME_LEVELS][N_TIME_LEVELS] = {{1., 0.}, {0, 1.}};

void CPU_AdvanceX( real u[][ FLU_NXT*FLU_NXT*FLU_NXT ], real Flux_Array[][NFLUX_TOTAL][ PS2*PS2 ], 
                   const real dt, const real dh, const real Eta, const bool StoreFlux, const real Taylor3_Coeff,  const real MinDens,
                   const int j_gap, const int k_gap, const int Flux_XYZ )
{

   const real _dh  = 1./dh; 
   const real _dh2 = SQR(_dh);
   const real _Eta2_dh = (real)0.5*_dh/Eta;

   const int j_start   = j_gap;
   const int k_start   = k_gap;
   const int j_end     = FLU_NXT - j_gap;
   const int k_end     = FLU_NXT - k_gap;

   real Rc[N_TIME_LEVELS][FLU_NXT];  // one column of the density in the input array "u" at all time levels
   real Pc[N_TIME_LEVELS][FLU_NXT];  // one column of the phase   in the input array "u" at all time levels
   real *Rc_target = NULL; // pointer to set density array in update loop
   real *Pc_target = NULL; // pointer to set phase   array in update loop
   real *Rc_N = NULL;    // pointer to store the full-step density
   real *Pc_N = NULL;    // pointer to store the full-step phase  
   int Idx;
   real Fm[N_TIME_LEVELS + 1][FLU_NXT];  // one column of the density flux at all time levels
   
   real *Fm_target = NULL; // pointer to set phase   array in update loop
   //change of density and phase in time step
   real ddensity, dphase;
   real v_C[FLU_NXT], v_L[FLU_NXT], v_R[FLU_NXT], rho_L[FLU_NXT], rho_R[FLU_NXT], log_density[FLU_NXT]; 
   real fp_L[FLU_NXT], fp_R[FLU_NXT], fm_L[FLU_NXT], fm_R[FLU_NXT]; 
   bool use_rk1[FLU_NXT]; 
   real rk1_density[FLU_NXT], rk1_phase[FLU_NXT]; 


   //indices to determine how far wrong information from evolving velocity field with wrong timestep can propagate
   int l_min, l_max; 

#  ifdef CONSERVE_MASS
   int Idx3;
   int Idx2;
   //real   R, I, dR, dI, Flux[PS2+1], R1, R2, I1, I2, dR1, dR2, dI1, dI2;
#  endif

// loop over all targeted columns
   for (int k=k_start; k<k_end; k++)
   for (int j=j_start; j<j_end; j++)
   {

//    1. backup one column of data
//    ------------------------------------------------------------------------------------------------------------
      Idx    = to1D(k,j,0);

      //Set pointers to output arrays
      Rc_N = &u[DENS][Idx];
      Pc_N = &u[PHAS][Idx];
      
      memcpy( Rc[0], Rc_N, FLU_NXT*sizeof(real) );
      memcpy( Pc[0], Pc_N, FLU_NXT*sizeof(real) );

      for (int i=0; i<FLU_NXT; i++) {
         use_rk1[i] = false;

#        ifdef GAMER_DEBUG
         if ( Rc[0][i] != Rc[0][i] || Pc[0][i] != Pc[0][i] || Rc[0][i] < 0)
            printf("Nan in input array of fluid solver k %d j %d i %d Rho %f Pha %f\n", k, j, i, Rc[0][i], Pc[0][i]); 
#        endif
      }

      for (int time_level = 0; time_level < N_TIME_LEVELS; ++time_level) 
      {
         int ghostBoundary = 4*time_level;

         for (int i = ghostBoundary    ; i < FLU_NXT - ghostBoundary    ; i++) {
            log_density[i] = 0.5 * log(FMAX(Rc[time_level][i], MinDens));
         }

         
         for (int i = ghostBoundary + 2; i < FLU_NXT - ghostBoundary - 2; i++) { 
            v_C[i]         = _dh * GRAD2(Pc[time_level], i);
         }
         
         
         computePPMInterpolation(v_C,             v_L,    v_R, false, true, ghostBoundary + 2);
         computePPMInterpolation(Rc[time_level], rho_L, rho_R, true, false, ghostBoundary);


//       computePPMFlux

         for (int i = ghostBoundary + 4; i < FLU_NXT - ghostBoundary - 4; i++) {
            real a    = Rc[time_level][i];
            real ap   = Rc[time_level][i+1];
            real am   = Rc[time_level][i-1];
            real a_L  = rho_L[i];
            real a_R  = rho_R[i];
            real a_Lp = rho_L[i+1];
            real a_Rp = rho_R[i+1];
            real a_Lm = rho_L[i-1];
            real a_Rm = rho_R[i-1];
            real vp2  = v_R[i]; 
            real vm2  = v_L[i]; 

            real d_a  =                      a_R  - a_L;
            real a_6  = real(6.0) * (a  - real(1.0/2.0) * ( a_R  + a_L ));

            real d_ap =                      a_Rp - a_Lp;
            real a_6p = real(6.0) * (ap - real(1.0/2.0) * ( a_Rp + a_Lp));

            real d_am =                      a_Rm - a_Lm;
            real a_6m = real(6.0) * (am - real(1.0/2.0) * ( a_Rm + a_Lm));

   //      Compute density fluxes at i+1/2 as seen by cells centered at i (fp_R) and i + 1 (fp_L)
            real y  =  vp2 * dt;
            real x  =  y / dh;
            fp_L[i] = a_R  - x/2.0 * (d_a   - ( 1.0 - 2.0/3.0 * x) * a_6  );

            y       = -vp2 * dt;
            x       =  y / dh;
            fp_R[i] = a_Lp + x/2.0 * (d_ap  + ( 1.0 - 2.0/3.0 * x) * a_6p );


            y       =  vm2 * dt;
            x       =  y / dh;
            fm_L[i] = a_Rm - x/2.0 * (d_am  - ( 1.0 - 2.0/3.0 * x) * a_6m);

            y       = -vm2 * dt;
            x       =  y / dh;
            fm_R[i] = a_L  + x/2.0 * (d_a   + ( 1.0 - 2.0/3.0 * x) * a_6 );



//          Enforce upwinding for density fluxes abar
            real a_bar_p2 = fp_L[i] * FMAX( vp2, 0.0 ) + fp_R[i] * FMIN( vp2, 0.0 );
            real a_bar_m2 = fm_L[i] * FMAX( vm2, 0.0 ) + fm_R[i] * FMIN( vm2, 0.0 );


            ddensity = _dh * (a_bar_p2 - a_bar_m2);



            real vm  = _dh * B3_GRADIENT(Pc[time_level], i);
            real vp  = _dh * F3_GRADIENT(Pc[time_level], i);
            dphase   = (SQR(MIN(vp, 0)) + SQR(MAX(vm, 0)))/2.0;

            dphase  -= _dh2/real(2.0) * (LAP2(log_density, i) +  SQR(GRAD2(log_density, i)));

            if (time_level == 0) {
               rk1_density[i] = FMAX(Rc[0][i] - dt / Eta * ddensity, MinDens); 
               rk1_phase[i]   =      Pc[0][i] - dt / Eta * dphase; 
            }

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

            *Rc_target = - TIME_COEFFS[time_level] * dt / Eta * ddensity;
            *Pc_target = - TIME_COEFFS[time_level] * dt / Eta * dphase;
            for (int tx = 0; tx < time_level + 1; ++tx) {
               *Rc_target += RK_COEFFS[time_level][tx] * Rc[tx][i];
               *Pc_target += RK_COEFFS[time_level][tx] * Pc[tx][i];
            }


//          check the velocity-dependent CFL-condition and switch to forward-Euler for updating the density wherever the CFL-condition is not met
//          dt = 1 / MaxdS_dx * 0.5 * ELBDM_ETA * DT__VELOCITY;
            if ( time_level + 1 ==  N_TIME_LEVELS ) {
               if ( (Rc_N[i] < 0 || Rc_N[i] != Rc_N[i] || Pc_N[i] != Pc_N[i]) && use_rk1[i] ) {             
                  Rc_N[i] = rk1_density[i];
                  Pc_N[i] = rk1_phase[i];               
               }
#              ifdef GAMER_DEBUG
               if (Rc_N[i] != Rc_N[i] || Pc_N[i] != Pc_N[i]) {
                  Rc_N[i] = MinDens;
                  Pc_N[i] = 0;
                  printf("FE after nan for dt = %f and dt_min = %f org. dens %f org. phas %f new dens = %f new phas = %f!\n", dt, dt_min, Rc[0][i], Pc[0][i], Rc_N[i], Pc_N[i]); 
               }
#              endif
            }

         }
      }


//    4 save the fluxes across all patch boundaries
//    ------------------------------------------------------------------------------------------------------------
#     ifdef CONSERVE_MASS
      for (int i=0; i<FLU_NXT; i++)
      {
         Fm[N_TIME_LEVELS][i] = 0;
         for (int time_level = 0; time_level < N_TIME_LEVELS; ++time_level) 
         {
         Fm[N_TIME_LEVELS][i] += FLUX_COEFFS[time_level] * Fm[time_level][i];
         }

         if ( use_rk1[i] )
            Fm[N_TIME_LEVELS][i] = Fm[0][i];
      }
      
      if ( StoreFlux )
      if (  ( j>=FLU_GHOST_SIZE && j<FLU_NXT-FLU_GHOST_SIZE )  &&  ( k>=FLU_GHOST_SIZE && k<FLU_NXT-FLU_GHOST_SIZE )  )
      {
         Idx3 = (k-FLU_GHOST_SIZE)*PS2 + (j-FLU_GHOST_SIZE);

         Flux_Array[Flux_XYZ+0][0][Idx3] = Fm[N_TIME_LEVELS][      FLU_GHOST_SIZE] / Eta;//Fm[0][      FLU_GHOST_SIZE];
         Flux_Array[Flux_XYZ+1][0][Idx3] = Fm[N_TIME_LEVELS][PS1 + FLU_GHOST_SIZE] / Eta;//Fm[0][PS1 + FLU_GHOST_SIZE];
         Flux_Array[Flux_XYZ+2][0][Idx3] = Fm[N_TIME_LEVELS][PS2 + FLU_GHOST_SIZE] / Eta;//Fm[0][PS2 + FLU_GHOST_SIZE];
      }
#     endif // #ifdef CONSERVE_MASS


   } // for j,k

} // FUNCTION : CPU_AdvanceX


//-------------------------------------------------------------------------------------------------------
// Function    :  TrasposeXY
// Description :  Transpose the x and y directions
//
// Parameter   :  u : Input wave function (density, phase)
//-------------------------------------------------------------------------------------------------------
void TransposeXY( real u[][ CUBE(FLU_NXT) ] )
{

   real (*u_xy)[ SQR(FLU_NXT) ] = new real [FLU_NIN][ SQR(FLU_NXT) ];
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

      for (int v=0; v<FLU_NIN; v++)    memcpy( &u[v][to1D(k,0,0)], u_xy[v], SQR(FLU_NXT)*sizeof(real) );
   }

   delete [] u_xy;

} // FUNCTION : TrasposeXY


template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

//-------------------------------------------------------------------------------------------------------
// Function    :  TrasposeXZ
// Description :  Transpose the x and z directions
//
// Parameter   :  u : Input wave function (density, phase)
//-------------------------------------------------------------------------------------------------------
void TransposeXZ( real u[][ CUBE(FLU_NXT) ] )
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

//Requires ghost zone of size 2 on either side
//Fills a_R[1 : -2], a_L[2 : -1]
void computePPMInterpolation(real* a_array, real* a_L_array, real* a_R_array, bool fix2, bool fix3, int ghostBoundary) {
   real a, ap, app, am, amm, delta_a, delta_ap, delta_m, delta_mp, ap2, a_L, a_R; 

// -------------------------------------------------------------------------
//  interpolate the cell-centered data to the edges
// -------------------------------------------------------------------------
   for (int i=1+ghostBoundary; i<FLU_NXT-2-ghostBoundary; i++) {
      a   = a_array[i];
      ap  = a_array[i+1];
      app = a_array[i+2];
      am  = a_array[i-1];

//    Average slope of parabola
      delta_a  = 0.5 * (ap - am);
      delta_ap = 0.5 * (app - a); 

      if ((ap - a) * ( a - am )   > 0.0)
         delta_m  = FMIN(FABS(delta_a ), 2.0 * FMIN(FABS(a - am), FABS(ap  - a ))) * sgn(delta_a);
      else
         delta_m = 0.0; 


      if ((app - ap) * ( ap - a ) > 0.0)
         delta_mp = FMIN(FABS(delta_ap), 2.0 * FMIN(FABS(ap - a), FABS(app - ap))) * sgn(delta_ap);
      else
         delta_mp = 0.0; 


      if (fix3)
         ap2 = real(7.0/12.0) * ( a + ap ) - real(1.0/12.0) * ( app + am );
      else
         ap2 = a + real(1.0/2.0) * ( ap - a ) - real(1.0/6.0) * (delta_mp - delta_m);


      a_R_array[i  ] = ap2; 
      a_L_array[i+1] = ap2; 


      if ( i == 1 + ghostBoundary) continue;

      a_R = a_R_array[i];
      a_L = a_L_array[i]; 


// -------------------------------------------------------------------------
//    Face-centered density approximations that take into account monotonicity
//    Potentially discontinuous, shock-resolving approximation to density
// -------------------------------------------------------------------------

/*
      if (fix1) {
//       Switch to different interpolation if we detect discontinuities in a 
//       Second derivative as measure for discontinuities
         d2_a  = 1/(6*dh**2) * (ap - 2*a + am)
         d2_ap = np.roll(d2_a, ROLL_R, axis = axis)
         d2_am = np.roll(d2_a, ROLL_L, axis = axis)

         if (-d2_ap * d2_am <= 0)
            eta_bar = 0
         else if (FABS(ap - am) - epsilon * FMIN(FABS(ap), FABS(am)) <= 0) 
            eta_bar = 0
         else 
            eta_bar = - ( (d2_ap - d2_am ) / ( 2 * dh) ) * ( 2*dh**3 / (ap - am) )


         eta = FMAX(0, FMIN(eta1 * (eta_bar - eta2), 1))


         a_Ld = am + 1/2 * delta_mm
         a_Rd = ap - 1/2 * delta_mp

         a_L = a_L * ( 1 - eta ) + a_Ld * eta
         a_R = a_R * ( 1 - eta ) + a_Rd * eta
*/
      if (fix2) {
//       Set coefficients of the interpolating parabola such that it does not overshoot
         
//        1. If a is local extremum, set the interpolation function to be constant
         if ((a_R - a)*(a - a_L) <= (0.0)) {
            a_L = a;
            a_R = a;
         }
         
//       2. If a between a_R and a_L, but very close to them, the parabola might still overshoot
         if (+ (a_R - a_L)*(a_R - a_L) / real(6.0) < (a_R - a_L) * (a - real(1.0/2.0) * (a_L + a_R)))
            a_L = real(3.0) * a - real(2.0) * a_R;
         if (- (a_R - a_L)*(a_R - a_L) / real(6.0) > (a_R - a_L) * (a - real(1.0/2.0) * (a_L + a_R)))
            a_R = real(3.0) * a - real(2.0) * a_L;

//       make sure that we didn't over or undersoot -- this may not
//       be needed, but is discussed in Colella & Sekora (2008)
         a_R = FMAX(a_R, FMIN(a, ap));
         a_R = FMIN(a_R, FMAX(a, ap));
         a_L = FMAX(a_L, FMIN(am, a));
         a_L = FMIN(a_L, FMAX(am, a));
      }


      a_R_array[i] = a_R;
      a_L_array[i] = a_L; 
   }

}





#endif // #if ( !defined GPU  &&  MODEL == ELBDM && ELBDM_SCHEME == HYBRID && HYBRID_SCHEME == PPM )
