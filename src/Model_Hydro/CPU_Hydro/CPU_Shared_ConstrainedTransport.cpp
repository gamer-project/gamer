#ifndef __CUFLU_CONSTRAINEDTRANSPORT__
#define __CUFLU_CONSTRAINEDTRANSPORT__



#include "CUFLU.h"

#if ( MODEL == HYDRO  &&  defined MHD )


// internal functions
GPU_DEVICE
static real dE_Upwind( const real FC_Ele_R, const real FC_Ele_L, const real FC_Mom,
                       const real PriVar_R[NCOMP_TOTAL_PLUS_MAG], const real PriVar_L[NCOMP_TOTAL_PLUS_MAG],
                       const int XYZ, const real dt_dh );




//-------------------------------------------------------------------------------------------------------
// Function    :  MHD_ComputeElectric
// Description :  Compute the edge-centered line-averaged electric field E=B x V (electromotive force; EMF)
//                for the constrained-transport algorithm
//
// Note        :  1. Ref : (a) Gardiner & Stone, J. Comput. Phys., 227, 4123 (2008)
//                         (b) Stone et al., ApJS, 178, 137 (2008)
//                2. This function is shared by MHM_RP and CTU schemes
//                3. g_EC_Ele [] has the size of N_EC_ELE^3  but is accessed with a stride "NEle"
//                   g_FC_Flux[] has the size of N_FC_FLUX^3 but is accessed with a stride "NFlux"
//                   g_PriVar [] has the size of FLU_NXT^3   but is accessed with a stride "NPri"
//                4. EMF-x/y/z( i, j, k ) are defined at the lower-left edge center of
//                   g_PriVar( i+OffsetPri+1, j+OffsetPri+1, k+OffsetPri+1 )
//
// Parameter   :  g_EC_Ele  : Array to store the output electric field
//                g_FC_Flux : Array storing the input face-centered fluxes
//                g_PriVar  : Array storing the input cell-centered primitive variables
//                NEle      : Stride for accessing g_EC_Ele[]
//                NFlux     : Stride for accessing g_FC_Flux[]
//                NPri      : Stride for accessing g_PriVar[]
//                OffsetPri : Offset for accessing g_PriVar[]
//                dt        : Time interval to advance solution
//                dh        : Cell size
//
// Return      :  g_EC_Ele[]
//------------------------------------------------------------------------------------------------------
GPU_DEVICE
void MHD_ComputeElectric(       real g_EC_Ele[][ CUBE(N_EC_ELE) ],
                          const real g_FC_Flux[][NCOMP_TOTAL_PLUS_MAG][ CUBE(N_FC_FLUX) ],
                          const real g_PriVar[][ CUBE(FLU_NXT) ],
                          const int NEle, const int NFlux, const int NPri, const int OffsetPri,
                          const real dt, const real dh )
{

   const int  NEleM1       = NEle - 1;
   const int  didx_flux[3] = { 1, NFlux, SQR(NFlux) };
   const int  didx_pri [3] = { 1, NPri,  SQR(NPri)  };
   const real dt_dh        = dt / dh;

   for (int d=0; d<3; d++)
   {
      const int TDir1 = (d+1)%3;             // transverse direction 1
      const int TDir2 = (d+2)%3;             // transverse direction 2
      const int TB1   = TDir1 + MAG_OFFSET;  // B flux component along the transverse direction 1
      const int TB2   = TDir2 + MAG_OFFSET;  // B flux component along the transverse direction 2

      int idx_ele_e[3], idx_flux_s[3];

      switch ( d )
      {
         case 0 : idx_ele_e [0] = NEleM1;  idx_ele_e [1] = NEle;    idx_ele_e [2] = NEle;
                  idx_flux_s[0] = 1;       idx_flux_s[1] = 0;       idx_flux_s[2] = 0;
                  break;

         case 1 : idx_ele_e [0] = NEle;    idx_ele_e [1] = NEleM1;  idx_ele_e [2] = NEle;
                  idx_flux_s[0] = 0;       idx_flux_s[1] = 1;       idx_flux_s[2] = 0;
                  break;

         case 2 : idx_ele_e [0] = NEle;    idx_ele_e [1] = NEle;    idx_ele_e [2] = NEleM1;
                  idx_flux_s[0] = 0;       idx_flux_s[1] = 0;       idx_flux_s[2] = 1;
                  break;
      }

      const int size_ij = idx_ele_e[0]*idx_ele_e[1];
      CGPU_LOOP( idx_ele, idx_ele_e[0]*idx_ele_e[1]*idx_ele_e[2] )
      {
         const int i_ele    = idx_ele % idx_ele_e[0];
         const int j_ele    = idx_ele % size_ij / idx_ele_e[0];
         const int k_ele    = idx_ele / size_ij;

         const int i_flux   = i_ele + idx_flux_s[0];
         const int j_flux   = j_ele + idx_flux_s[1];
         const int k_flux   = k_ele + idx_flux_s[2];
         const int idx_flux = IDX321( i_flux, j_flux, k_flux, NFlux, NFlux );

         const int i_pri    = i_flux + OffsetPri;
         const int j_pri    = j_flux + OffsetPri;
         const int k_pri    = k_flux + OffsetPri;
         const int idx_pri  = IDX321( i_pri, j_pri, k_pri, NPri, NPri );

         g_EC_Ele[d][idx_ele] = ( - g_FC_Flux[TDir1][TB2][ idx_flux + didx_flux[TDir2] ]
                                  - g_FC_Flux[TDir1][TB2][ idx_flux                    ]
                                  + g_FC_Flux[TDir2][TB1][ idx_flux + didx_flux[TDir1] ]
                                  + g_FC_Flux[TDir2][TB1][ idx_flux                    ] );

         g_EC_Ele[d][idx_ele] += dE_Upwind( -g_FC_Flux[ idx_flux + didx_flux[TDir2] ][TDir1][TB2],
                                            -g_FC_Flux[ idx_flux                    ][TDir1][TB2],
                                             g_FC_Flux[ idx_flux                    ][TDir2][  0],
                                             g_PriVar[ idx_pri + didx_pri[TDir2] ],
                                             g_PriVar[ idx_pri                   ],
                                             d, dt_dh );

         g_EC_Ele[d][idx_ele] += dE_Upwind( -g_FC_Flux[ idx_flux + didx_flux[TDir2] ][TDir1][TB2],
                                            -g_FC_Flux[ idx_flux                    ][TDir1][TB2],
                                             g_FC_Flux[ idx_flux + didx_flux[TDir1] ][TDir2][  0],
                                             g_PriVar[ idx_pri + didx_pri[TDir1] + didx_pri[TDir2] ],
                                             g_PriVar[ idx_pri + didx_pri[TDir1]                   ],
                                             d, dt_dh );


         g_EC_Ele[d][idx_ele] += dE_Upwind( +g_FC_Flux[ idx_flux + didx_flux[TDir1] ][TDir2][TB1],
                                            +g_FC_Flux[ idx_flux                    ][TDir2][TB1],
                                             g_FC_Flux[ idx_flux                    ][TDir1][  0],
                                             g_PriVar[ idx_pri + didx_pri[TDir1] ],
                                             g_PriVar[ idx_pri                   ],
                                             d, dt_dh );

         g_EC_Ele[d][idx_ele] += dE_Upwind( +g_FC_Flux[ idx_flux + didx_flux[TDir1] ][TDir2][TB1],
                                            +g_FC_Flux[ idx_flux                    ][TDir2][TB1],
                                             g_FC_Flux[ idx_flux + didx_flux[TDir2] ][TDir1][  0],
                                             g_PriVar[ idx_pri + didx_pri[TDir2] + didx_pri[TDir1] ],
                                             g_PriVar[ idx_pri + didx_pri[TDir2]                   ],
                                             d, dt_dh );

         g_EC_Ele[d][idx_ele] *= (real)0.25;

      } // CGPU_LOOP( idx_ele, idx_ele_e[0]*idx_ele_e[1]*idx_ele_e[2] )
   } // for ( int d=0; d<3; d++)


#  ifdef __CUDACC__
   __syncthreads();
#  endif

} // FUNCTION : MHD_ComputeElectric



//-------------------------------------------------------------------------------------------------------
// Function    :  dE_Upwind
// Description :  Calculate the first partial derivative of electric field with the upwind scheme
//
// Note        :  1. Ref : Gardiner & Stone, J. Comput. Phys., 227, 4123 (2008)
//                2. Invoked by MHD_ComputeElectric()
//
// Parameter   :  FC_Ele_R/L : Right/left face-centered electric field
//                FC_Mom     : Face-centered momentum for determining the upwind direction
//                PriVar_R/L : Cell-centered primitive variable for computing the reference electric field
//                XYZ        : Target spatial direction : (0/1/2) --> (x/y/z)
//                dt_dh      : dt/dh --> for normalizing velocity only
//
// Return      :  dE/dx
//------------------------------------------------------------------------------------------------------
GPU_DEVICE
real dE_Upwind( const real FC_Ele_R, const real FC_Ele_L, const real FC_Mom,
                const real PriVar_R[NCOMP_TOTAL_PLUS_MAG], const real PriVar_L[NCOMP_TOTAL_PLUS_MAG],
                const int XYZ, const real dt_dh )
{

   const int TDir1 = (XYZ+1)%3;           // transverse direction 1
   const int TDir2 = (XYZ+2)%3;           // transverse direction 2
   const int TV1   = TDir1 + 1;           // velocity component along the transverse direction 1
   const int TV2   = TDir2 + 1;           // velocity component along the transverse direction 2
   const int TB1   = TDir1 + MAG_OFFSET;  // B field  component along the transverse direction 1
   const int TB2   = TDir2 + MAG_OFFSET;  // B field  component along the transverse direction 2

// convert dimensional momentum to dimensionless velocity to reduce the effect of round-off errors
   const real FC_Vel = (real)2.0*dt_dh*FC_Mom/( PriVar_R[0] + PriVar_L[0] );

   real dE, CC_Ele_R, CC_Ele_L;  // CC_Ele_R/L: right/left cell-centered electric field

// MAX_ERROR is defined in CUFLU.h
   if ( FABS(FC_Vel) <= MAX_ERROR )
   {
      CC_Ele_R = PriVar_R[TB1]*PriVar_R[TV2] - PriVar_R[TB2]*PriVar_R[TV1];
      CC_Ele_L = PriVar_L[TB1]*PriVar_L[TV2] - PriVar_L[TB2]*PriVar_L[TV1];
      dE       = (real)0.5*( FC_Ele_R - CC_Ele_R + FC_Ele_L - CC_Ele_L );
   }

   else if ( FC_Vel > (real)0.0 )
   {
      CC_Ele_L = PriVar_L[TB1]*PriVar_L[TV2] - PriVar_L[TB2]*PriVar_L[TV1];
      dE       = FC_Ele_L - CC_Ele_L;
   }

   else
   {
      CC_Ele_R = PriVar_R[TB1]*PriVar_R[TV2] - PriVar_R[TB2]*PriVar_R[TV1];
      dE       = FC_Ele_R - CC_Ele_R;
   }

   return dE;

} // FUNCTION : dE_Upwind



//-------------------------------------------------------------------------------------------------------
// Function    :  MHD_UpdataMagnetic
// Description :  Update magnetic field with the constrained transport algorithm
//
// Note        :  1. This function is shared by MHM_RP and CTU schemes
//                2. g_FC_Bx/y/z_Out[] are accessed with a stride "NOut"
//                   g_FC_B_In[] has the size of FLU_NXT_P1^3 and is also accessed with the same stride
//                   g_EC_Ele[] has the size of N_EC_ELE^3 but is accessed with a stride "NEle"
//
// Parameter   :  g_FC_B_Out  : Array to store the output face-centered B field
//                              --> Separate into three arrays since the array dimension is different
//                                  during the half- and full-step updates
//                g_FC_B_In   : Array storing the input face-centered B field
//                g_EC_Ele    : Array storing the input edge-centered electric field
//                dt          : Time interval to advance solution
//                dh          : Cell size
//                NOut        : Stride for accessing g_FC_Bx/y/z_Out[]
//                NEle        : Stride for accessing g_EC_Ele[]
//                Offset_B_In : Offset for accessing g_FC_B_In[]
//
// Return      :  g_FC_B_Out[]
//------------------------------------------------------------------------------------------------------
GPU_DEVICE
void MHD_UpdateMagnetic( real *g_FC_Bx_Out, real *g_FC_By_Out, real *g_FC_Bz_Out,
                         const real g_FC_B_In[][ FLU_NXT_P1*SQR(FLU_NXT) ],
                         const real g_EC_Ele[][ CUBE(N_EC_ELE) ],
                         const real dt, const real dh, const int NOut, const int NEle, const int Offset_B_In )
{

   const int  NOutP1      = NOut + 1;
   const int  didx_ele[3] = { 1, NEle, SQR(NEle) };
   const real dt_dh       = dt / dh;

   real *g_FC_B_Out[3] = { g_FC_Bx_Out, g_FC_By_Out, g_FC_Bz_Out };
   real dE1, dE2;

   for (int d=0; d<3; d++)
   {
      const int TDir1 = (d+1)%3;    // transverse direction 1
      const int TDir2 = (d+2)%3;    // transverse direction 2

      int idx_out_e_i, idx_out_e_j, idx_out_e_k, stride_in_i, stride_in_j;

      switch ( d )
      {
         case 0 : idx_out_e_i = NOutP1;      idx_out_e_j = NOut;        idx_out_e_k = NOut;
                  stride_in_i = FLU_NXT_P1;  stride_in_j = FLU_NXT;
                  break;

         case 1 : idx_out_e_i = NOut;        idx_out_e_j = NOutP1;      idx_out_e_k = NOut;
                  stride_in_i = FLU_NXT;     stride_in_j = FLU_NXT_P1;
                  break;

         case 2 : idx_out_e_i = NOut;        idx_out_e_j = NOut;        idx_out_e_k = NOutP1;
                  stride_in_i = FLU_NXT;     stride_in_j = FLU_NXT;
                  break;
      }

      const int size_ij = idx_out_e_i*idx_out_e_j;
      CGPU_LOOP( idx_out, idx_out_e_i*idx_out_e_j*idx_out_e_k )
      {
         const int i_out   = idx_out % idx_out_e_i;
         const int j_out   = idx_out % size_ij / idx_out_e_i;
         const int k_out   = idx_out / size_ij;

         const int idx_ele = IDX321( i_out, j_out, k_out, NEle, NEle );

         const int i_in    = i_out + Offset_B_In;
         const int j_in    = j_out + Offset_B_In;
         const int k_in    = k_out + Offset_B_In;
         const int idx_in  = IDX321( i_in, j_in, k_in, stride_in_i, stride_in_j );

          dE1 = g_EC_Ele[TDir1][ idx_ele + didx_ele[TDir2] ] - g_EC_Ele[TDir1][idx_ele];
          dE2 = g_EC_Ele[TDir2][ idx_ele + didx_ele[TDir1] ] - g_EC_Ele[TDir2][idx_ele];

          g_FC_B_Out[d][idx_out] = g_FC_B_In[d][idx_in] + dt_dh*( dE1 - dE2 );
      } // CGPU_LOOP( idx_out, idx_out_e_i*idx_out_e_j*idx_out_e_k )
   } // for (int d=0; d<3; d++)

} // FUNCTION : MHD_UpdateMagnetic



#endif // #if ( MODEL == HYDRO  &&  defined MHD )



#endif // #ifndef __CUFLU_CONSTRAINEDTRANSPORT__
