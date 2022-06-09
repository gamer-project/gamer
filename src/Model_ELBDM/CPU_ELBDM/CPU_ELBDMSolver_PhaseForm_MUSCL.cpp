#include "GAMER.h"
#include "CUFLU.h"

#include <iostream>


#if ( !defined GPU  &&  MODEL == ELBDM  && ELBDM_SCHEME == PHASE)



// useful macros
#define to1D(z,y,x)     ( z*FLU_NXT*FLU_NXT + y*FLU_NXT + x )

#ifdef LAPLACIAN_4TH
#  define LAP1(In,t)    (  real(1.0/ 12.0)*( - In[t-2] + (real)16.0*In[t-1] - (real)30.0*In[t] + \
                                             - In[t+2] + (real)16.0*In[t+1] )  )
#  define LAP2(In,t)    (  real(1.0/144.0)*(  In[t-4] - (real)32.0*In[t-3] + (real)316.0*In[t-2] - (real)992.0*In[t-1] + \
                                              In[t+4] - (real)32.0*In[t+3] + (real)316.0*In[t+2] - (real)992.0*In[t+1] + \
                                              (real)1414.0*In[t] )  )
#else
#  define LAP1(In,t)    ( In[t-1] - (real)2.0*In[t] + In[t+1] )
#  define LAP2(In,t)    ( In[t-2] - (real)4.0*In[t-1] + (real)6.0*In[t] - (real)4.0*In[t+1] + In[t+2] )
#endif

# define GTR( a, b )     (  ( (a) > (b) ) ? (1) : (0)  )
# define LSS( a, b )     (  ( (a) < (b) ) ? (1) : (0)  )
             

# define CENTERED_GRADIENT(In, t) ( real(1.0/2.0 ) * (   In[t + 1] - In[t - 1] ) )
# define BACKWARD_GRADIENT(In, t) ( In[t    ] - In[t - 1] )
# define FORWARD_GRADIENT(In, t)  ( In[t + 1] - In[t    ] )
# define LAPLACIAN(In, t)         ( In[t - 1] - (real)2.0*In[t] + In[t + 1] )

// VAN ALBADA LIMITER
# define LIMITER(In) ( (SQR(In) + In)/((real)1. + SQR(In)) )

# define UPWIND_GRADIENT_RATIO(Rc, Vb, t) ( ((Rc[t-1] - Rc[t-2]) * GTR(Vb[t], 0)  \
                                           + (Rc[t+1] - Rc[t  ]) * LSS(Vb[t], 0)) \
                                           / (Rc[t] - Rc[t-1] + (((Rc[t] - Rc[t-1]) == 0) ? 1e-8 : 0)))

//Second-order MUSCL flux reconstruction
# define MUSCL_FLUX(Rc, Vb, t, dx, dt) (  MAX(Vb[t], 0) * Rc[t-1] \
                                       +  MIN(Vb[t], 0) * Rc[t  ] \
                                       +  real(0.5) * FABS(Vb[t]) * (1. - FABS(Vb[t] * dt/dx)) * LIMITER(UPWIND_GRADIENT_RATIO(Rc, Vb, t)) * (Rc[t] -Rc[t - 1]) )

//Ratio of subsequent backward gradients
# define BACKWARD_GRADIENT_RATIO(Pc, t) ((Pc[t] - Pc[t-1]) / (Pc[t-1] - Pc[t-2] + (((Pc[t-1] - Pc[t-2]) == 0) ?  1e-8 : 0)))



static void CPU_AdvanceX( real u[][ CUBE(FLU_NXT) ], real Flux_Array[][NFLUX_TOTAL][ SQR(PS2) ],
                          const real dt, const real dh, const real Eta, const bool StoreFlux, const real Taylor3_Coeff,
                          const int j_gap, const int k_gap, const int Flux_XYZ );
static void TransposeXY( real u[][ CUBE(FLU_NXT) ] );
static void TransposeXZ( real u[][ CUBE(FLU_NXT) ] );




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
//                MinDens        : Minimum allowed density
//-------------------------------------------------------------------------------------------------------
void CPU_ELBDMSolver_PhaseForm_MUSCL( real Flu_Array_In [][FLU_NIN ][ CUBE(FLU_NXT) ],
                      real Flu_Array_Out[][FLU_NOUT][ CUBE(PS2) ],
                      real Flux_Array[][9][NFLUX_TOTAL][ SQR(PS2) ],
                      const int NPatchGroup, const real dt, const real dh, const real Eta, const bool StoreFlux,
                      const real Taylor3_Coeff, const bool XYZ, const real MinDens )
{

   if ( XYZ )
   {
#     pragma omp parallel for schedule( runtime )
      for (int P=0; P<NPatchGroup; P++)
      {
         CPU_AdvanceX( Flu_Array_In[P], Flux_Array[P], dt, dh, Eta, StoreFlux, Taylor3_Coeff,
                                    0,              0, 0 );

         TransposeXY ( Flu_Array_In[P] );

         CPU_AdvanceX( Flu_Array_In[P], Flux_Array[P], dt, dh, Eta, StoreFlux, Taylor3_Coeff,
                       FLU_GHOST_SIZE,              0, 3 );

         TransposeXZ ( Flu_Array_In[P] );

         CPU_AdvanceX( Flu_Array_In[P], Flux_Array[P], dt, dh, Eta, StoreFlux, Taylor3_Coeff,
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

         CPU_AdvanceX( Flu_Array_In[P], Flux_Array[P], dt, dh, Eta, StoreFlux, Taylor3_Coeff,
                                    0,              0, 6 );

         TransposeXZ ( Flu_Array_In[P] );

         CPU_AdvanceX( Flu_Array_In[P], Flux_Array[P], dt, dh, Eta, StoreFlux, Taylor3_Coeff,
                                    0, FLU_GHOST_SIZE, 3 );

         TransposeXY ( Flu_Array_In[P] );

         CPU_AdvanceX( Flu_Array_In[P], Flux_Array[P], dt, dh, Eta, StoreFlux, Taylor3_Coeff,
                       FLU_GHOST_SIZE, FLU_GHOST_SIZE, 0 );
      }
   }


// copy the updated data to Flu_Array_Out
   int  Idx1, Idx2, v_m1;
   real Amp, Rescale;   // not using double precision since MinDens will break the mass conservation anyway

#  pragma omp parallel for private( Idx1, Idx2, v_m1, Amp, Rescale ) schedule( runtime )
   for (int P=0; P<NPatchGroup; P++)
   {
//    copy data
      for (int v=1; v<FLU_NOUT; v++)
      {
         v_m1 = v-1;
         Idx1 = 0;

         for (int k=FLU_GHOST_SIZE; k<FLU_GHOST_SIZE+PS2; k++)
         for (int j=FLU_GHOST_SIZE; j<FLU_GHOST_SIZE+PS2; j++)
         for (int i=FLU_GHOST_SIZE; i<FLU_GHOST_SIZE+PS2; i++)
         {
            Idx2 = to1D(k,j,i);

            Flu_Array_Out[P][v][ Idx1++ ] = Flu_Array_In[P][v_m1][Idx2];
         }
      }

//    evaluate the new density (and apply the minimum density check)
      for (int t=0; t<CUBE(PS2); t++)
      {
         Amp = SQR( Flu_Array_Out[P][1][t] ) + SQR( Flu_Array_Out[P][2][t] );

         if ( Amp < MinDens )
         {
            Rescale                 = SQRT( MinDens / Amp );
            Flu_Array_Out[P][1][t] *= Rescale;
            Flu_Array_Out[P][2][t] *= Rescale;
            Amp                     = MinDens;
         }

         Flu_Array_Out[P][0][t] = Amp;
      }
   } // for (int P=0; P<NPatchGroup; P++)

} // FUNCTION : CPU_ELBDMSolver



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
const real RK_COEFFS [N_TIME_LEVELS][N_TIME_LEVELS] = {{1., 0, 0}, {3./4, 1./4, 0}, {1./3, 0, 2./3}};

void CPU_AdvanceX( real u[][ FLU_NXT*FLU_NXT*FLU_NXT ], real Flux_Array[][NFLUX_TOTAL][ PS2*PS2 ], 
                   const real dt, const real dh, const real Eta, const bool StoreFlux, const real Taylor3_Coeff, 
                   const int j_gap, const int k_gap, const int Flux_XYZ )
{

   const real _dh  = 1./dh; 
   const real _dh2 = SQR(_dh);

   const int j_start   = j_gap;
   const int k_start   = k_gap;
   const int j_end     = FLU_NXT - j_gap;
   const int k_end     = FLU_NXT - k_gap;

   const real dT       = (real)0.5*dt/Eta;
   const real _Eta2_dh = (real)0.5*_dh/Eta;
   const real Coeff1   = dT*_dh*_dh;
   const real Coeff2   = Taylor3_Coeff*Coeff1*Coeff1;

   real Rc[N_TIME_LEVELS][FLU_NXT];  // one column of the density in the input array "u" at all time levels
   real Pc[N_TIME_LEVELS][FLU_NXT];  // one column of the phase   in the input array "u" at all time levels
   real *Rc_target = NULL; // pointer to set density array in update loop
   real *Pc_target = NULL; // pointer to set phase   array in update loop
   real *Re_N = NULL;    // pointer to store the full-step density
   real *Im_N = NULL;    // pointer to store the full-step phase  
   int Idx;

   //slope-limite
   real ql, qc, Qc, Qr;
   //velocities dS/dx = v and density fluxes f at i - 1/2, i + 1/2 
   real vm, vp, fm, fp;
   //change of density and phase in time step
   real ddensity, dphase;
   //temporary variables to convert density and phase to real and imaginary parts
   real re, im;
   real ql_ratios [FLU_NXT], backward_velocities[FLU_NXT], log_density[FLU_NXT];  // one column of the gradient ratios for phase, velocities dS/dx and log(rho)


   real Re_Old [FLU_NXT];  // one column of the real      part in the input array "u"
   real Im_Old [FLU_NXT];  // one column of the imaginary part in the input array "u"
   real Re_Half[FLU_NXT];  // one column of the real      part at the half time-step
   real Im_Half[FLU_NXT];  // one column of the imaginary part at the half time-step
   real *Re_New = NULL;    // pointer to store the full-step real      part
   real *Im_New = NULL;    // pointer to store the full-step imaginary part

   real vmax = 0;
   
// loop over all targeted columns
   for (int k=k_start; k<k_end; k++)
   for (int j=j_start; j<j_end; j++)
   {

//    1. backup one column of data
//    ------------------------------------------------------------------------------------------------------------
      Idx    = to1D(k,j,0);

      //Set pointers to output arrays
      Re_N = &u[0][Idx];
      Im_N = &u[1][Idx];
      
      //for (int i=0; i<FLU_NXT; i++) {
      //   std::cout<<"initial i: " << i << " Re : " << Re_N[i] << " Im: "<< Im_N[i] << "\n";
      //}
      
      for (int i=0; i<FLU_NXT; i++) {
         //Set input data by converting real and imaginary parts to density Rc and phase Pc
         Rc[0][i] = SQR(Re_N[i]) + SQR(Im_N[i]);
         Pc[0][i] = atan2(Im_N[i], Re_N[i]);
         //Make sure that phase is continuous in input array by requiring that it may no jump by more than 2 pi between neighbouring points
         if (i > 0)
            Pc[0][i] = ELBDM_UnwrapPhase(Pc[0][i-1], Pc[0][i]);
      }

      /*
      for (int i=0; i<FLU_NXT; i++) {
         re = sqrt(Rc[0][i]) * cos(Pc[0][i]);
         im = sqrt(Rc[0][i]) * sin(Pc[0][i]);
         Re_Old[i] = re;
         Im_Old[i] = im;
      }

      //memcpy( Rc[0], Rc_N, FLU_NXT*sizeof(real) );
      //memcpy( Pc[0], Pc_N, FLU_NXT*sizeof(real) );


//    2. half-step solution
//    ------------------------------------------------------------------------------------------------------------
      for (int i=2; i<FLU_NXT-2; i++)
      {
         Re_Half[i] = Re_Old[i] - (real)0.5*Coeff1*LAP1( Im_Old, i ) - Coeff2*LAP2( Re_Old, i );
         Im_Half[i] = Im_Old[i] + (real)0.5*Coeff1*LAP1( Re_Old, i ) - Coeff2*LAP2( Im_Old, i );
      }


//    3. full-step solution (equivalent to the 3rd-order Taylor expansion)
//    ------------------------------------------------------------------------------------------------------------
      for (int i=FLU_GHOST_SIZE; i<FLU_NXT-FLU_GHOST_SIZE; i++)
      {
         Re_N[i] = Re_Old[i] - Coeff1*LAP1( Im_Half, i );
         Im_N[i] = Im_Old[i] + Coeff1*LAP1( Re_Half, i );
      }*/


      
      for (int time_level = 0; time_level < N_TIME_LEVELS; ++time_level) 
      {
         
         //for (int i=2*(time_level + 1); i<FLU_NXT-2*(time_level + 1); i++)
         //   std::cout<<"time level :" << time_level << "i: " << i << " Rho : " << Rc[time_level][i] << " Phase: "<< Pc[time_level][i] << "\n";
         
         
         //Compute second-order backward velocities
         for (int i=2*(time_level + 1); i<FLU_NXT-time_level*2; i++)
         {
            ql_ratios[i]           =       BACKWARD_GRADIENT_RATIO(Pc[time_level], i);
            backward_velocities[i] = _dh * BACKWARD_GRADIENT      (Pc[time_level], i);
            vmax = MAX(vmax, FABS(backward_velocities[i]));
         }

         //Compute density logarithms
         for (int i=2*time_level; i<FLU_NXT-2*time_level; i++)
         {
            log_density[i] = log(Rc[time_level][i]);
         }
               
         //Update density and phase fields
         for (int i=2*(time_level + 1); i<FLU_NXT-2*(time_level + 1); i++)
         {
            fm = MUSCL_FLUX(Rc[time_level], backward_velocities, i    , dh, dt);
            fp = MUSCL_FLUX(Rc[time_level], backward_velocities, i + 1, dh, dt);

            //Slope-limited second order gradients (useful if phase develops discontinuities, could happen in hybrid scheme with 2 pi jump)
            ql = ql_ratios[i  ];
            qc = ql_ratios[i+1];
            Qc = 1./(ql_ratios[i+1] + ((ql_ratios[i+1] == 0) ? 1e-8 : 0));
            Qr = 1./(ql_ratios[i+2] + ((ql_ratios[i+2] == 0) ? 1e-8 : 0));

            vm = _dh * BACKWARD_GRADIENT(Pc[time_level], i) * ( (real) 1. + (real) 0.5 * LIMITER(qc) - (real) 0.5 * LIMITER(ql)/(ql + ((ql==0) ? 1e-8 : 0)));
            vp = _dh * FORWARD_GRADIENT (Pc[time_level], i) * ( (real) 1. + (real) 0.5 * LIMITER(Qc) - (real) 0.5 * LIMITER(Qr)/(Qr + ((Qr==0) ? 1e-8 : 0)));
   
            ddensity = _dh * (fp - fm);
            dphase   = (SQR(MIN(vp, 0)) + SQR(MAX(vm, 0)))/2;
            dphase  += -_dh2/4 * LAPLACIAN(log_density, i) - _dh2/8 * SQR(CENTERED_GRADIENT(log_density, i));


            if (time_level + 1 < N_TIME_LEVELS)
            {
               Rc_target = &Rc[time_level + 1][i];
               Pc_target = &Pc[time_level + 1][i];
            }
            else 
            {
               Rc_target = &Re_N[i];
               Pc_target = &Im_N[i];
            }

            *Rc_target = - TIME_COEFFS[time_level] * dt * Eta * ddensity;
            *Pc_target = - TIME_COEFFS[time_level] * dt * Eta * dphase;
            for (int tx = 0; tx < time_level + 1; ++tx) {
               *Rc_target += RK_COEFFS[time_level][tx] * Rc[tx][i];
               *Pc_target += RK_COEFFS[time_level][tx] * Pc[tx][i];
            }

            if (time_level + 1 == N_TIME_LEVELS) {
               re = sqrt(Re_N[i]) * cos(Im_N[i]);
               im = sqrt(Re_N[i]) * sin(Im_N[i]);
               Re_N[i] = re;
               Im_N[i] = im;
            }
         }
      }
      
      //for (int i=0; i<FLU_NXT; i++) {
      //   std::cout<<"final i: " << i << " Re : " << Re_N[i] << " Im: "<< Im_N[i] << "\n";
      //}


   } // for j,k

} // FUNCTION : CPU_AdvanceX


//-------------------------------------------------------------------------------------------------------
// Function    :  TrasposeXY
// Description :  Transpose the x and y directions
//
// Parameter   :  u : Input wave function (density, real, imaginary)
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



//-------------------------------------------------------------------------------------------------------
// Function    :  TrasposeXZ
// Description :  Transpose the x and z directions
//
// Parameter   :  u : Input wave function (density, real, imaginary)
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



#endif // #if ( !defined GPU  &&  MODEL == ELBDM )
