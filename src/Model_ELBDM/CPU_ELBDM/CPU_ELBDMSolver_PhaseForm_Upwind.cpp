#include "GAMER.h"
#include "CUFLU.h"

//###################################################################################################
#if (  !defined GPU  &&  MODEL == ELBDM  &&  ELBDM_SCHEME == HYBRID && HYBRID_SCHEME == UPWIND )
//###################################################################################################



// useful macros
#define to1D(z,y,x)     ( z*FLU_NXT*FLU_NXT + y*FLU_NXT + x )

#  define LAP1(In,t)    ( In[t-1] - (real)2.0*In[t] + In[t+1] )
#  define LAP2(In,t)    ( - (real)1./12 * In[t-2] + (real)4./3*In[t-1] - (real)5./2*In[t] - (real)4./3*In[t+1] - (real)1./12*In[t+2] )
#  define GRADF(In, t)  ( In[t+1] - In[t]   )
#  define GRADB(In, t)  ( In[t]   - In[t-1] )
#  define GRAD1(In, t)  ( In[t+1] - In[t-1] )
#  define GRAD2(In, t)  ( -1./6 * In[t+2] + 4./3 * In[t+1] - 4./3 * In[t-1] + 1./6 * In[t-2] )

static void CPU_AdvanceX( real u[][ CUBE(FLU_NXT)], real Flux_Array[][NFLUX_TOTAL][ SQR(PS2) ], 
                          const real dt, const real dh, const real Eta, const bool StoreFlux, const real Taylor3_Coeff, 
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
void CPU_ELBDMSolver_PhaseForm_Upwind( real Flu_Array_In [][FLU_NIN ][ CUBE(FLU_NXT)], 
                      real Flu_Array_Out[][FLU_NOUT][ SQR(PS2)*PS2 ], 
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
                   const real dt, const real dh, const real Eta, const bool StoreFlux, const real Taylor3_Coeff, 
                   const int j_gap, const int k_gap, const int Flux_XYZ )
{
   const real Coeff1 = dt/(dh*dh);

   const int j_start   = j_gap;
   const int k_start   = k_gap;
   const int j_end     = FLU_NXT - j_gap;
   const int k_end     = FLU_NXT - k_gap;

   real Rh_Old [FLU_NXT];  // one column of the real      part in the input array "u" 
   real Ph_Old [FLU_NXT];  // one column of the imaginary part in the input array "u"
   real Sr_Old [FLU_NXT];  // one column of the imaginary part in the input array "u"
   real Rh_Half[FLU_NXT];  // one column of the real      part at the half time-step 
   real Ph_Half[FLU_NXT];  // one column of the imaginary part at the half time-step 
   real Sr_Half[FLU_NXT];  // one column of the imaginary part at the half time-step 
   real *Rh_New = NULL;    // pointer to store the full-step real      part
   real *Ph_New = NULL;    // pointer to store the full-step imaginary part
   int Idx;

   real vp, vm, fp, fm, qp, osf;
   

// loop over all targeted columns
   for (int k=k_start; k<k_end; k++)
   for (int j=j_start; j<j_end; j++)
   {

//    1. backup one column of data
//    ------------------------------------------------------------------------------------------------------------
      Idx    = to1D(k,j,0);
      Rh_New = &u[0][Idx];
      Ph_New = &u[1][Idx];

      memcpy( Rh_Old, Rh_New, FLU_NXT*sizeof(real) );
      memcpy( Ph_Old, Ph_New, FLU_NXT*sizeof(real) );

      for (int i=0; i<FLU_NXT; i++) {
         Sr_Old[i] = 0.5 * log(Rh_Old[i]);
      }


//    2. half-step solution
//    ------------------------------------------------------------------------------------------------------------
      for (int i=2; i<FLU_NXT-2; i++)
      {
         vp  = GRADF(Ph_Old, i);
         vm  = GRADB(Ph_Old, i);
         fp  = MAX(vp, 0) * Rh_Old[i    ] + MIN(vp, 0) * Rh_Old[i + 1];
         fm  = MAX(vm, 0) * Rh_Old[i - 1] + MIN(vm, 0) * Rh_Old[i    ];
         qp  = (real) 0.5 * ((real) 0.25 * pow(GRAD1(Sr_Old, i), 2) + LAP1(Sr_Old, i));
         osf = (real) 0.5 * (pow(MIN(vp, 0), 2) + pow(MAX(vm, 0), 2));
         
         Rh_Half[i] = Rh_Old[i] - (real) 0.5 * Coeff1 * (fp - fm);
         Ph_Half[i] = Ph_Old[i] - (real) 0.5 * Coeff1 * (osf - qp);
         Sr_Half[i] = 0.5 * log(Rh_Half[i]);
      }


//    3. full-step solution (equivalent to the 3rd-order Taylor expansion)
//    ------------------------------------------------------------------------------------------------------------
      for (int i=FLU_GHOST_SIZE; i<FLU_NXT-FLU_GHOST_SIZE; i++)
      {
         vp  = GRADF(Ph_Half, i);
         vm  = GRADB(Ph_Half, i);
         fp  = MAX(vp, 0) * Rh_Half[i    ] + MIN(vp, 0) * Rh_Half[i + 1];
         fm  = MAX(vm, 0) * Rh_Half[i - 1] + MIN(vm, 0) * Rh_Half[i    ];
         qp  = (real) 0.5 * ((real) 0.25 * pow(GRAD1(Sr_Half, i), 2) + LAP1(Sr_Half, i));
         osf = (real) 0.5 * (pow(MIN(vp, 0), 2) + pow(MAX(vm, 0), 2));
         
         Rh_New[i] = Rh_Old[i] - (real) 1.0 * Coeff1 * (fp - fm);
         Ph_New[i] = Ph_Old[i] - (real) 1.0 * Coeff1 * (osf - qp);
      }
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



#endif // #if ( MODEL == ELBDM  &&  ELBDM_SCHEME == HYBRID && HYBRID_SCHEME == UPWIND )
