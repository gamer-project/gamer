#ifndef __CUFLU_CONSTRAINEDTRANSPORT__
#define __CUFLU_CONSTRAINEDTRANSPORT__



#include "CUFLU.h"

#if ( MODEL == HYDRO  &&  defined MHD )



// external functions
#ifdef __CUDACC__
#include "CUFLU_Shared_FluUtility.cu"
#endif


// internal functions
GPU_DEVICE
static real dE_Upwind( const real FC_Ele_L, const real FC_Ele_R, const real FC_Mom, const real D_L, const real D_R,
                       const real V_L1, const real V_L2, const real V_R1, const real V_R2,
                       const real B_L1, const real B_L2, const real B_R1, const real B_R2,
                       const real dt_dh );
#ifdef UNSPLIT_GRAVITY
GPU_DEVICE
void UpdateVelocityByGravity( real &v1, real &v2, const int TDir1, const int TDir2,
                              const int i_usg, const int j_usg, const int k_usg,
                              const real dt_half, const double dh_f8, const real GraConst,
                              const real g_Pot_USG[], const double Corner_USG[], const double Time,
                              const bool UsePot, const OptExtAcc_t ExtAcc, const ExtAcc_t ExtAcc_Func,
                              const double ExtAcc_AuxArray[] );
#endif




//-------------------------------------------------------------------------------------------------------
// Function    :  MHD_ComputeElectric
// Description :  Compute the edge-centered line-averaged electric field E=B x V (electromotive force; EMF)
//                for the constrained-transport algorithm
//
// Note        :  1. Ref : (a) Gardiner & Stone, J. Comput. Phys., 227, 4123 (2008)
//                         (b) Stone et al., ApJS, 178, 137 (2008)
//                2. This function is shared by MHM, MHM_RP and CTU schemes
//                3. g_EC_Ele [] has the size of N_EC_ELE^3  but is accessed with a stride "NEle"
//                   --> But there are only NEle-1 useful elements along x/y/z for Ex/Ey/Ez, respectively
//                   g_FC_Flux[] has the size of N_FC_FLUX^3 but is accessed with a stride "NFlux"
//                   g_PriVar [] has the size of FLU_NXT^3   but is accessed with a stride "NPri"
//                4. EMF-x/y/z( i, j, k ) are defined at the lower-left edge center of
//                   g_PriVar( i+OffsetPri+1, j+OffsetPri+1, k+OffsetPri+1 )
//                5. Store the electric field at the patch boundaries for correcting the coarse-grid B field
//                   --> Option "DumpIntEle"
//                   --> Structure of g_IntEle[] = [face index][E field index][cell index]
//                          Face index: [0/1/2] = left/middle/right x faces
//                                      [3/4/5] = left/middle/right y faces
//                                      [6/7/8] = left/middle/right z faces
//                          E field index on x faces: [0/1] = Ey/Ez
//                                           y faces: [0/1] = Ez/Ex
//                                           z faces: [0/1] = Ex/Ey
//                          Cell index on x faces: [Nz][Ny] (= [PS2+1][PS2] for Ey and [PS2][PS2+1] for Ez)
//                                     on y faces: [Nx][Nz] (= [PS2+1][PS2] for Ez and [PS2][PS2+1] for Ex)
//                                     on z faces: [Ny][Nx] (= [PS2+1][PS2] for Ex and [PS2][PS2+1] for Ey)
//                6. For the unsplitting scheme in gravity (i.e., UNSPLIT_GRAVITY), this function also corrects the half-step
//                   velocity by gravity when CorrHalfVel==true
//
// Parameter   :  g_EC_Ele        : Array to store the output electric field
//                g_FC_Flux       : Array storing the input face-centered fluxes
//                g_PriVar        : Array storing the input cell-centered primitive variables
//                NEle            : Stride for accessing g_EC_Ele[]
//                NFlux           : Stride for accessing g_FC_Flux[]
//                NPri            : Stride for accessing g_PriVar[]
//                OffsetPri       : Offset for accessing g_PriVar[]
//                dt              : Time interval to advance solution
//                dh              : Cell size
//                DumpIntEle      : Store the inter-patch electric field (i.e., E field at the patch boundaries)
//                                  in g_IntEle[]
//                g_IntEle        : Array for DumpIntEle
//                CorrHalfVel     : true --> correct the half-step velocity by gravity       (for UNSPLIT_GRAVITY only)
//                g_Pot_USG       : Array storing the input potential for CorrHalfVel        (for UNSPLIT_GRAVITY only)
//                g_Corner        : Array storing the corner coordinates of each patch group (for UNSPLIT_GRAVITY only)
//                Time            : Current physical time                                    (for UNSPLIT_GRAVITY only)
//                UsePot          : Add self-gravity and/or external potential               (for UNSPLIT_GRAVITY only)
//                ExtAcc          : Add external acceleration                                (for UNSPLIT_GRAVITY only)
//                ExtAcc_Func     : Function pointer to the external acceleration routine    (for UNSPLIT_GRAVITY only)
//                ExtAcc_AuxArray : Auxiliary array for external acceleration                (for UNSPLIT_GRAVITY only)
//
// Return      :  g_EC_Ele[], g_IntEle[]
//------------------------------------------------------------------------------------------------------
GPU_DEVICE
void MHD_ComputeElectric(       real g_EC_Ele[][ CUBE(N_EC_ELE) ],
                          const real g_FC_Flux[][NCOMP_TOTAL_PLUS_MAG][ CUBE(N_FC_FLUX) ],
                          const real g_PriVar[][ CUBE(FLU_NXT) ],
                          const int NEle, const int NFlux, const int NPri, const int OffsetPri,
                          const real dt, const real dh,
                          const bool DumpIntEle, real g_IntEle[][NCOMP_ELE][ PS2P1*PS2 ],
                          const bool CorrHalfVel, const real g_Pot_USG[], const double g_Corner[], const double Time,
                          const bool UsePot, const OptExtAcc_t ExtAcc, const ExtAcc_t ExtAcc_Func,
                          const double ExtAcc_AuxArray[] )
{

// check
#  ifdef GAMER_DEBUG
#  ifdef UNSPLIT_GRAVITY
   if ( CorrHalfVel )
   {
      if ( UsePot  &&  g_Pot_USG == NULL )
         printf( "ERROR : g_Pot_USG == NULL !!\n" );

      if ( ExtAcc  &&  g_Corner == NULL )
         printf( "ERROR : g_Corner == NULL !!\n" );
   }
#  else
   if ( CorrHalfVel )
      printf( "ERROR : CorrHalfVel is NOT supported when UNSPLIT_GRAVITY is off !!\n" );
#  endif
#  endif // #ifdef GAMER_DEBUG


   const int  NEleM1       = NEle - 1;
   const int  didx_flux[3] = { 1, NFlux, SQR(NFlux) };
   const int  didx_pri [3] = { 1, NPri,  SQR(NPri)  };
   const real dt_dh        = dt / dh;

#  ifdef UNSPLIT_GRAVITY
   const double dh_f8       = (double)dh;
   const real   GraConst    = -(real)0.25*dt_dh;
   const real   dt_half     = (real)0.5*dt;
   const int    pri_ghost   = ( NPri - PS2 )/2;                // number of ghost zones on each side for g_PriVar[]
   const int    idx_pri2usg = USG_GHOST_SIZE_F - pri_ghost;    // index difference between g_PriVar[] and g_Pot_USG[]

   int    ijk_usg[3];
   double Corner_USG[3];   // central coordinates of the 0th cell in g_Pot_USG[]
   if ( CorrHalfVel  &&  ExtAcc )
      for (int d=0; d<3; d++)    Corner_USG[d] = g_Corner[d] - dh_f8*USG_GHOST_SIZE_F;

// check
#  ifdef GAMER_DEBUG
   if ( CorrHalfVel  &&  idx_pri2usg + OffsetPri < 1 )
         printf( "ERROR : idx_pri2usg (%d) + OffsetPri (%d) < 1 (USG_GHOST_SIZE_F %d, NPri %d) !!\n",
                 idx_pri2usg, OffsetPri, USG_GHOST_SIZE_F, NPri );
#  endif
#  endif // #ifdef UNSPLIT_GRAVITY


   for (int d=0; d<3; d++)
   {
      const int TDir1 = (d+1)%3;             // transverse direction 1
      const int TDir2 = (d+2)%3;             // transverse direction 2
      const int TV1   = TDir1 + 1;           // velocity component along the transverse direction 1
      const int TV2   = TDir2 + 1;           // velocity component along the transverse direction 2
      const int TB1   = TDir1 + MAG_OFFSET;  // B flux   component along the transverse direction 1
      const int TB2   = TDir2 + MAG_OFFSET;  // B flux   component along the transverse direction 2

      int idx_ele_e[2], idx_flux_s[3];

      switch ( d )
      {
         case 0 : idx_ele_e [0] = NEleM1;  idx_ele_e [1] = NEle;
                  idx_flux_s[0] = 1;       idx_flux_s[1] = 0;       idx_flux_s[2] = 0;
                  break;

         case 1 : idx_ele_e [0] = NEle;    idx_ele_e [1] = NEleM1;
                  idx_flux_s[0] = 0;       idx_flux_s[1] = 1;       idx_flux_s[2] = 0;
                  break;

         case 2 : idx_ele_e [0] = NEle;    idx_ele_e [1] = NEle;
                  idx_flux_s[0] = 0;       idx_flux_s[1] = 0;       idx_flux_s[2] = 1;
                  break;
      }

      const int size_ij = idx_ele_e[0]*idx_ele_e[1];
      CGPU_LOOP( idx0, NEleM1*SQR(NEle)  )
      {
         const int i_ele    = idx0 % idx_ele_e[0];
         const int j_ele    = idx0 % size_ij / idx_ele_e[0];
         const int k_ele    = idx0 / size_ij;
         const int idx_ele  = IDX321( i_ele, j_ele, k_ele, NEle, NEle );

         const int i_flux   = i_ele + idx_flux_s[0];
         const int j_flux   = j_ele + idx_flux_s[1];
         const int k_flux   = k_ele + idx_flux_s[2];
         const int idx_flux = IDX321( i_flux, j_flux, k_flux, NFlux, NFlux );

         const int i_pri    = i_flux + OffsetPri;
         const int j_pri    = j_flux + OffsetPri;
         const int k_pri    = k_flux + OffsetPri;
         const int idx_pri  = IDX321( i_pri, j_pri, k_pri, NPri, NPri );

         real D_L, D_R, V_L1, V_L2, V_R1, V_R2, B_L1, B_L2, B_R1, B_R2;
         int  idx_L, idx_R;
         real Ele_Out;

//       compute the edge-centered electric field
         Ele_Out = ( - g_FC_Flux[TDir1][TB2][ idx_flux + didx_flux[TDir2] ]
                     - g_FC_Flux[TDir1][TB2][ idx_flux                    ]
                     + g_FC_Flux[TDir2][TB1][ idx_flux + didx_flux[TDir1] ]
                     + g_FC_Flux[TDir2][TB1][ idx_flux                    ] );


         idx_L = idx_pri;
         idx_R = idx_L + didx_pri[TDir2];
         D_L   = g_PriVar[  0][ idx_L ];
         V_L1  = g_PriVar[TV1][ idx_L ];
         V_L2  = g_PriVar[TV2][ idx_L ];
         B_L1  = g_PriVar[TB1][ idx_L ];
         B_L2  = g_PriVar[TB2][ idx_L ];
         D_R   = g_PriVar[  0][ idx_R ];
         V_R1  = g_PriVar[TV1][ idx_R ];
         V_R2  = g_PriVar[TV2][ idx_R ];
         B_R1  = g_PriVar[TB1][ idx_R ];
         B_R2  = g_PriVar[TB2][ idx_R ];

//       correct the half-step velocity by gravity for the unsplitting scheme
#        ifdef UNSPLIT_GRAVITY
         if ( CorrHalfVel )
         {
            ijk_usg[0] = i_pri + idx_pri2usg;
            ijk_usg[1] = j_pri + idx_pri2usg;
            ijk_usg[2] = k_pri + idx_pri2usg;
            UpdateVelocityByGravity( V_L1, V_L2, TDir1, TDir2, ijk_usg[0], ijk_usg[1], ijk_usg[2], dt_half, dh_f8,
                                     GraConst, g_Pot_USG, Corner_USG, Time, UsePot, ExtAcc, ExtAcc_Func, ExtAcc_AuxArray );

            ijk_usg[TDir2] ++;
            UpdateVelocityByGravity( V_R1, V_R2, TDir1, TDir2, ijk_usg[0], ijk_usg[1], ijk_usg[2], dt_half, dh_f8,
                                     GraConst, g_Pot_USG, Corner_USG, Time, UsePot, ExtAcc, ExtAcc_Func, ExtAcc_AuxArray );
         } // if ( CorrHalfVel )
#        endif // #ifdef UNSPLIT_GRAVITY

         Ele_Out += dE_Upwind( -g_FC_Flux[TDir1][TB2][ idx_flux                    ],
                               -g_FC_Flux[TDir1][TB2][ idx_flux + didx_flux[TDir2] ],
                                g_FC_Flux[TDir2][  0][ idx_flux                    ],
                               D_L, D_R, V_L1, V_L2, V_R1, V_R2, B_L1, B_L2, B_R1, B_R2, dt_dh );


         idx_L = idx_pri + didx_pri[TDir1];
         idx_R = idx_L   + didx_pri[TDir2];
         D_L   = g_PriVar[  0][ idx_L ];
         V_L1  = g_PriVar[TV1][ idx_L ];
         V_L2  = g_PriVar[TV2][ idx_L ];
         B_L1  = g_PriVar[TB1][ idx_L ];
         B_L2  = g_PriVar[TB2][ idx_L ];
         D_R   = g_PriVar[  0][ idx_R ];
         V_R1  = g_PriVar[TV1][ idx_R ];
         V_R2  = g_PriVar[TV2][ idx_R ];
         B_R1  = g_PriVar[TB1][ idx_R ];
         B_R2  = g_PriVar[TB2][ idx_R ];

//       correct the half-step velocity by gravity for the unsplitting scheme
#        ifdef UNSPLIT_GRAVITY
         if ( CorrHalfVel )
         {
            ijk_usg[0] = i_pri + idx_pri2usg;
            ijk_usg[1] = j_pri + idx_pri2usg;
            ijk_usg[2] = k_pri + idx_pri2usg;
            ijk_usg[TDir1] ++;
            UpdateVelocityByGravity( V_L1, V_L2, TDir1, TDir2, ijk_usg[0], ijk_usg[1], ijk_usg[2], dt_half, dh_f8,
                                     GraConst, g_Pot_USG, Corner_USG, Time, UsePot, ExtAcc, ExtAcc_Func, ExtAcc_AuxArray );

            ijk_usg[TDir2] ++;
            UpdateVelocityByGravity( V_R1, V_R2, TDir1, TDir2, ijk_usg[0], ijk_usg[1], ijk_usg[2], dt_half, dh_f8,
                                     GraConst, g_Pot_USG, Corner_USG, Time, UsePot, ExtAcc, ExtAcc_Func, ExtAcc_AuxArray );
         } // if ( CorrHalfVel )
#        endif // #ifdef UNSPLIT_GRAVITY

         Ele_Out += dE_Upwind( -g_FC_Flux[TDir1][TB2][ idx_flux                    ],
                               -g_FC_Flux[TDir1][TB2][ idx_flux + didx_flux[TDir2] ],
                                g_FC_Flux[TDir2][  0][ idx_flux + didx_flux[TDir1] ],
                               D_L, D_R, V_L1, V_L2, V_R1, V_R2, B_L1, B_L2, B_R1, B_R2, dt_dh );


         idx_L = idx_pri;
         idx_R = idx_L + didx_pri[TDir1];
         D_L   = g_PriVar[  0][ idx_L ];
         V_L1  = g_PriVar[TV1][ idx_L ];
         V_L2  = g_PriVar[TV2][ idx_L ];
         B_L1  = g_PriVar[TB1][ idx_L ];
         B_L2  = g_PriVar[TB2][ idx_L ];
         D_R   = g_PriVar[  0][ idx_R ];
         V_R1  = g_PriVar[TV1][ idx_R ];
         V_R2  = g_PriVar[TV2][ idx_R ];
         B_R1  = g_PriVar[TB1][ idx_R ];
         B_R2  = g_PriVar[TB2][ idx_R ];

//       correct the half-step velocity by gravity for the unsplitting scheme
#        ifdef UNSPLIT_GRAVITY
         if ( CorrHalfVel )
         {
            ijk_usg[0] = i_pri + idx_pri2usg;
            ijk_usg[1] = j_pri + idx_pri2usg;
            ijk_usg[2] = k_pri + idx_pri2usg;
            UpdateVelocityByGravity( V_L1, V_L2, TDir1, TDir2, ijk_usg[0], ijk_usg[1], ijk_usg[2], dt_half, dh_f8,
                                     GraConst, g_Pot_USG, Corner_USG, Time, UsePot, ExtAcc, ExtAcc_Func, ExtAcc_AuxArray );

            ijk_usg[TDir1] ++;
            UpdateVelocityByGravity( V_R1, V_R2, TDir1, TDir2, ijk_usg[0], ijk_usg[1], ijk_usg[2], dt_half, dh_f8,
                                     GraConst, g_Pot_USG, Corner_USG, Time, UsePot, ExtAcc, ExtAcc_Func, ExtAcc_AuxArray );
         } // if ( CorrHalfVel )
#        endif // #ifdef UNSPLIT_GRAVITY

         Ele_Out += dE_Upwind( +g_FC_Flux[TDir2][TB1][ idx_flux                    ],
                               +g_FC_Flux[TDir2][TB1][ idx_flux + didx_flux[TDir1] ],
                                g_FC_Flux[TDir1][  0][ idx_flux                    ],
                               D_L, D_R, V_L1, V_L2, V_R1, V_R2, B_L1, B_L2, B_R1, B_R2, dt_dh );


         idx_L = idx_pri + didx_pri[TDir2];
         idx_R = idx_L   + didx_pri[TDir1];
         D_L   = g_PriVar[  0][ idx_L ];
         V_L1  = g_PriVar[TV1][ idx_L ];
         V_L2  = g_PriVar[TV2][ idx_L ];
         B_L1  = g_PriVar[TB1][ idx_L ];
         B_L2  = g_PriVar[TB2][ idx_L ];
         D_R   = g_PriVar[  0][ idx_R ];
         V_R1  = g_PriVar[TV1][ idx_R ];
         V_R2  = g_PriVar[TV2][ idx_R ];
         B_R1  = g_PriVar[TB1][ idx_R ];
         B_R2  = g_PriVar[TB2][ idx_R ];

//       correct the half-step velocity by gravity for the unsplitting scheme
#        ifdef UNSPLIT_GRAVITY
         if ( CorrHalfVel )
         {
            ijk_usg[0] = i_pri + idx_pri2usg;
            ijk_usg[1] = j_pri + idx_pri2usg;
            ijk_usg[2] = k_pri + idx_pri2usg;
            ijk_usg[TDir2] ++;
            UpdateVelocityByGravity( V_L1, V_L2, TDir1, TDir2, ijk_usg[0], ijk_usg[1], ijk_usg[2], dt_half, dh_f8,
                                     GraConst, g_Pot_USG, Corner_USG, Time, UsePot, ExtAcc, ExtAcc_Func, ExtAcc_AuxArray );

            ijk_usg[TDir1] ++;
            UpdateVelocityByGravity( V_R1, V_R2, TDir1, TDir2, ijk_usg[0], ijk_usg[1], ijk_usg[2], dt_half, dh_f8,
                                     GraConst, g_Pot_USG, Corner_USG, Time, UsePot, ExtAcc, ExtAcc_Func, ExtAcc_AuxArray );
         } // if ( CorrHalfVel )
#        endif // #ifdef UNSPLIT_GRAVITY

         Ele_Out += dE_Upwind( +g_FC_Flux[TDir2][TB1][ idx_flux                    ],
                               +g_FC_Flux[TDir2][TB1][ idx_flux + didx_flux[TDir1] ],
                                g_FC_Flux[TDir1][  0][ idx_flux + didx_flux[TDir2] ],
                               D_L, D_R, V_L1, V_L2, V_R1, V_R2, B_L1, B_L2, B_R1, B_R2, dt_dh );


         Ele_Out *= (real)0.25;

//       store the electric field of all cells in g_EC_Ele[]
         g_EC_Ele[d][idx_ele] = Ele_Out;

//       store the inter-patch electric field in g_IntEle[]
         if ( DumpIntEle )
         {
//          sanity check: this function assumes N_FL_ELE == PS2+1
#           if ( N_FL_ELE != PS2+1 )
#              error : ERROR : N_FL_ELE != PS2+1 !!
#           endif

            if      ( d == 0 ) {
               if ( j_ele == 0 || j_ele == PS1 || j_ele == PS2 )  g_IntEle[ 3 + j_ele/PS1 ][1][ i_ele*PS2P1 + k_ele ] = Ele_Out;
               if ( k_ele == 0 || k_ele == PS1 || k_ele == PS2 )  g_IntEle[ 6 + k_ele/PS1 ][0][ j_ele*PS2   + i_ele ] = Ele_Out;
            } // d == 0

            else if ( d == 1 ) {
               if ( k_ele == 0 || k_ele == PS1 || k_ele == PS2 )  g_IntEle[ 6 + k_ele/PS1 ][1][ j_ele*PS2P1 + i_ele ] = Ele_Out;
               if ( i_ele == 0 || i_ele == PS1 || i_ele == PS2 )  g_IntEle[ 0 + i_ele/PS1 ][0][ k_ele*PS2   + j_ele ] = Ele_Out;
            } // d == 1

            else {
               if ( i_ele == 0 || i_ele == PS1 || i_ele == PS2 )  g_IntEle[ 0 + i_ele/PS1 ][1][ k_ele*PS2P1 + j_ele ] = Ele_Out;
               if ( j_ele == 0 || j_ele == PS1 || j_ele == PS2 )  g_IntEle[ 3 + j_ele/PS1 ][0][ i_ele*PS2   + k_ele ] = Ele_Out;
            } // d == 2
         } // if ( DumpIntEle )
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
// Parameter   :  FC_Ele_L/R : Left/right face-centered electric field
//                FC_Mom     : Face-centered momentum for determining the upwind direction
//                D_L/R      : Left/right cell-centered density
//                V_L/R_1/2  : Left/right cell-centered velocity along the transverse direction 1/2
//                B_L/R_1/2  : Left/right cell-centered B field  along the transverse direction 1/2
//                dt_dh      : dt/dh --> for normalizing velocity only
//
// Return      :  dE
//------------------------------------------------------------------------------------------------------
GPU_DEVICE
real dE_Upwind( const real FC_Ele_L, const real FC_Ele_R, const real FC_Mom, const real D_L, const real D_R,
                const real V_L1, const real V_L2, const real V_R1, const real V_R2,
                const real B_L1, const real B_L2, const real B_R1, const real B_R2,
                const real dt_dh )
{

// convert dimensional momentum to dimensionless velocity to reduce the effect of round-off errors
   const real FC_Vel = (real)2.0*dt_dh*FC_Mom/( D_L + D_R );

   real dE, CC_Ele_L, CC_Ele_R;  // CC_Ele_L/R: left/right cell-centered electric field

// MAX_ERROR is defined in Macro.h
   if ( FABS(FC_Vel) <= MAX_ERROR )
   {
      CC_Ele_R = B_R1*V_R2 - B_R2*V_R1;
      CC_Ele_L = B_L1*V_L2 - B_L2*V_L1;
      dE       = (real)0.5*( FC_Ele_R - CC_Ele_R + FC_Ele_L - CC_Ele_L );
   }

   else if ( FC_Vel > (real)0.0 )
   {
      CC_Ele_L = B_L1*V_L2 - B_L2*V_L1;
      dE       = FC_Ele_L - CC_Ele_L;
   }

   else
   {
      CC_Ele_R = B_R1*V_R2 - B_R2*V_R1;
      dE       = FC_Ele_R - CC_Ele_R;
   }

   return dE;

} // FUNCTION : dE_Upwind



#ifdef UNSPLIT_GRAVITY
//-------------------------------------------------------------------------------------------------------
// Function    :  UpdateVelocityByGravity
// Description :  Update the half-step velocity by gravity for the unsplitting scheme
//
// Note        :  1. Invoked by MHD_ComputeElectric()
//                2. Only used when enabling UNSPLIT_GRAVITY
//
// Parameter   :  v1/2            : Velocity along the 1st/2nd transverse directions to be updated
//                                  --> call-by-reference
//                TDir1/2         : 1st/2nd transverse directions
//                i/j/k_usg       : Array indices of g_Pot_USG[]
//                dt_half         : 0.5*dt
//                dh_f8           : Cell size in double precision
//                GraConst        : -0.25*dt/dh
//                g_Pot_USG       : Input potential array
//                Corner_USG      : Cell-centered coordinates of the 0th cell in g_Pot_USG[]
//                Time            : Current physical time
//                UsePot          : Add self-gravity and/or external potential
//                ExtAcc          : Add external acceleration
//                ExtAcc_Func     : Function pointer to the external acceleration routine
//                ExtAcc_AuxArray : Auxiliary array for external acceleration
//
// Return      :  v1, v2
//------------------------------------------------------------------------------------------------------
GPU_DEVICE
void UpdateVelocityByGravity( real &v1, real &v2, const int TDir1, const int TDir2,
                              const int i_usg, const int j_usg, const int k_usg,
                              const real dt_half, const double dh_f8, const real GraConst,
                              const real g_Pot_USG[], const double Corner_USG[], const double Time,
                              const bool UsePot, const OptExtAcc_t ExtAcc, const ExtAcc_t ExtAcc_Func,
                              const double ExtAcc_AuxArray[] )
{

   const int didx_usg[3] = { 1, USG_NXT_F, SQR(USG_NXT_F) };
   real Acc[3] = { (real)0.0, (real)0.0, (real)0.0 };

// external acceleration
   if ( ExtAcc )
   {
      double xyz[3]; // cell-centered coordinates

      xyz[0] = Corner_USG[0] + i_usg*dh_f8;  // always use double precision to calculate coordinates
      xyz[1] = Corner_USG[1] + j_usg*dh_f8;
      xyz[2] = Corner_USG[2] + k_usg*dh_f8;

      ExtAcc_Func( Acc, xyz[0], xyz[1], xyz[2], Time, ExtAcc_AuxArray );

      Acc[TDir1] *= dt_half;
      Acc[TDir2] *= dt_half;
   }

// self-gravity and external potential
   if ( UsePot )
   {
      const int idx_usg = IDX321( i_usg, j_usg, k_usg, USG_NXT_F, USG_NXT_F );

      Acc[TDir1] += GraConst*( g_Pot_USG[ idx_usg + didx_usg[TDir1] ] - g_Pot_USG[ idx_usg - didx_usg[TDir1] ] );
      Acc[TDir2] += GraConst*( g_Pot_USG[ idx_usg + didx_usg[TDir2] ] - g_Pot_USG[ idx_usg - didx_usg[TDir2] ] );
   }

// advance velocity by gravity
   v1 += Acc[TDir1];
   v2 += Acc[TDir2];

} // FUNCTION : UpdateVelocityByGravity
#endif // #ifdef UNSPLIT_GRAVITY



//-------------------------------------------------------------------------------------------------------
// Function    :  MHD_UpdateMagnetic
// Description :  Update magnetic field with the constrained transport algorithm
//
// Note        :  1. This function is shared by MHM, MHM_RP and CTU schemes
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

      int idx_out_e_i, idx_out_e_j, idx_out_e_k, size_in_i, size_in_j;

      switch ( d )
      {
         case 0 : idx_out_e_i = NOutP1;      idx_out_e_j = NOut;        idx_out_e_k = NOut;
                  size_in_i   = FLU_NXT_P1;  size_in_j   = FLU_NXT;
                  break;

         case 1 : idx_out_e_i = NOut;        idx_out_e_j = NOutP1;      idx_out_e_k = NOut;
                  size_in_i   = FLU_NXT;     size_in_j   = FLU_NXT_P1;
                  break;

         case 2 : idx_out_e_i = NOut;        idx_out_e_j = NOut;        idx_out_e_k = NOutP1;
                  size_in_i   = FLU_NXT;     size_in_j   = FLU_NXT;
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
         const int idx_in  = IDX321( i_in, j_in, k_in, size_in_i, size_in_j );

          dE1 = g_EC_Ele[TDir1][ idx_ele + didx_ele[TDir2] ] - g_EC_Ele[TDir1][idx_ele];
          dE2 = g_EC_Ele[TDir2][ idx_ele + didx_ele[TDir1] ] - g_EC_Ele[TDir2][idx_ele];

          g_FC_B_Out[d][idx_out] = g_FC_B_In[d][idx_in] + dt_dh*( dE1 - dE2 );
      } // CGPU_LOOP( idx_out, idx_out_e_i*idx_out_e_j*idx_out_e_k )
   } // for (int d=0; d<3; d++)


#  ifdef __CUDACC__
   __syncthreads();
#  endif

} // FUNCTION : MHD_UpdateMagnetic



//-------------------------------------------------------------------------------------------------------
// Function    :  MHD_HalfStepPrimitive
// Description :  Evaluate the half-step cell-centered primitive variables
//
// Note        :  1. Used by the CTU scheme
//                2. Use face-centered fluxes for the conservative update and then convert momentum to velocity
//                3. Cell-centered B field is simply obtained by averaging the half-step face-centered B field
//                4. Cell-centered primitive variables are only used for computing the edge-centered
//                   electric field, which is then used for the full-step CT update
//                   --> Only need to calculate velocity and B field
//                   --> Skip energy and passive scalars
//                       --> No need to apply the dual-energy formalism
//                5. g_Flu_In[]     has the size of FLU_NXT^3 and is accessed with the same stride
//                   g_FC_B_Half[]  has the size of FLU_NXT_P1*SQR(FLU_NXT) but is accessed with the dimension
//                                  (N_HF_VAR+1)*SQR(N_HF_VAR)
//                   g_PriVar_Out[] has the size of FLU_NXT^3 but is accessed with a stride N_HF_VAR
//                                  --> g_FC_B_Half[] contains all the face-centered B field for g_PriVar_Out[],
//                                      nothing more, nothing less
//                                  --> Just like the relation between g_Flu_Array_In[] and g_Mag_Array_In[]
//                                  --> One can invoke MHD_GetCellCenteredBField() to compute the cell-centered
//                                      B field directly
//                   g_Flux[]       has the size of N_FC_FLUX^3 but is accessed with a stride N_HF_FLUX
//                                  --> Although currently we have N_FC_FLUX == N_HF_FLUX
//
// Parameter   :  g_Flu_In     : Array storing the input initial cell-centered fluid data
//                g_FC_B_Half  : Array storing the input half-step face-centered B field
//                g_PriVar_Out : Array to store the output half-step primitive variables
//                g_Flux       : Array storing the input face-centered fluxes
//                dt           : Full-step time interval
//                dh           : Cell size
//                MinDens      : Minimum allowed density
//
// Return      :  g_PriVar_Out[]
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE
void MHD_HalfStepPrimitive( const real g_Flu_In[][ CUBE(FLU_NXT) ],
                            const real g_FC_B_Half[][ FLU_NXT_P1*SQR(FLU_NXT) ],
                                  real g_PriVar_Out[][ CUBE(FLU_NXT) ],
                            const real g_Flux[][NCOMP_TOTAL_PLUS_MAG][ CUBE(N_FC_FLUX) ],
                            const real dt, const real dh, const real MinDens )
{

   const int  didx_flux[3] = { 1, N_HF_FLUX, SQR(N_HF_FLUX) };
   const real dt_dh2       = (real)0.5*dt/dh;
   const int  NFluVar      = NCOMP_FLUID - 1;   // density + momentum*3

   real dFlux[3][NFluVar], Output_1Cell[ NFluVar + NCOMP_MAG ];

   const int size_ij = SQR(N_HF_VAR);
   CGPU_LOOP( idx_out, CUBE(N_HF_VAR) )
   {
      const int i_out      = idx_out % N_HF_VAR;
      const int j_out      = idx_out % size_ij / N_HF_VAR;
      const int k_out      = idx_out / size_ij;

      const int i_flux     = i_out + 1;
      const int j_flux     = j_out + 1;
      const int k_flux     = k_out + 1;
      const int idx_flux   = IDX321( i_flux, j_flux, k_flux, N_HF_FLUX, N_HF_FLUX );

      const int i_flu_in   = i_out + FLU_GHOST_SIZE - 1;    // assuming N_HF_VAR = PS2+2
      const int j_flu_in   = j_out + FLU_GHOST_SIZE - 1;
      const int k_flu_in   = k_out + FLU_GHOST_SIZE - 1;
      const int idx_flu_in = IDX321( i_flu_in, j_flu_in, k_flu_in, FLU_NXT, FLU_NXT );


//    1. calculate flux difference to update the fluid data by 0.5*dt
      for (int d=0; d<3; d++)
      for (int v=0; v<NFluVar; v++)
         dFlux[d][v] = g_Flux[d][v][idx_flux] - g_Flux[d][v][ idx_flux - didx_flux[d] ];

      for (int v=0; v<NFluVar; v++)
         Output_1Cell[v] = g_Flu_In[v][idx_flu_in] - dt_dh2*( dFlux[0][v] + dFlux[1][v] + dFlux[2][v] );

//    apply density floor
      Output_1Cell[DENS] = FMAX( Output_1Cell[DENS], MinDens );

//    check negative density
#     ifdef CHECK_UNPHYSICAL_IN_FLUID
      Hydro_IsUnphysical_Single( Output_1Cell[DENS], "density", TINY_NUMBER, HUGE_NUMBER, ERROR_INFO, UNPHY_VERBOSE );
#     endif


//    2. momentum --> velocity
      const real _Dens = (real)1.0 / Output_1Cell[DENS];
      for (int v=1; v<NFluVar; v++)    Output_1Cell[v] *= _Dens;


//    3. compute the cell-centered half-step B field
      MHD_GetCellCenteredBField( Output_1Cell+NFluVar, g_FC_B_Half[0], g_FC_B_Half[1], g_FC_B_Half[2],
                                 N_HF_VAR, N_HF_VAR, N_HF_VAR, i_out, j_out, k_out );


//    4. store results to the output array
//       --> variable indices in g_PriVar_Out[] remain consistent with other arrays even though
//           energy and passive scalars have not been skipped here
      for (int v=0; v<NFluVar; v++)    g_PriVar_Out[ v              ][idx_out] = Output_1Cell[ v           ];
      for (int v=0; v<NCOMP_MAG; v++)  g_PriVar_Out[ v + MAG_OFFSET ][idx_out] = Output_1Cell[ v + NFluVar ];

   } // CGPU_LOOP( idx_out, CUBE(N_HF_VAR) )


#  ifdef __CUDACC__
   __syncthreads();
#  endif

} // FUNCTION : MHD_HalfStepPrimitive



# if ( FLU_SCHEME == MHM )
//-------------------------------------------------------------------------------------------------------
// Function    :  MHD_ComputeElectric_Half
// Description :  Compute the edge-centered line-averaged electric field E=B x V (electromotive force; EMF)
//                for the constrained-transport algorithm
//
// Note        :  1. This function can only be used by MHM scheme
//                2. Ref : (a) Teyssier R., 2002, A&A, 385, 337. doi:10.1051/0004-6361:20011817
//                         (b) Fromang, S., Hennebelle, P., & Teyssier, R. 2006, aap, 457, 371. doi:10.1051/0004-6361:20065371
//                3. g_EC_Ele [] has the size of N_EC_ELE^3 but is accessed with a stride "NEle"
//                   --> But there are only NEle-1 useful elements along x/y/z for Ex/Ey/Ez, respectively
//                   g_ConVar [] has the size of FLU_NXT^3  but is accessed with a stride "NCon"
//                4. EMF-x/y/z( i, j, k ) are defined at the lower-left edge center of
//                   g_ConVar( i+OffsetCon, j+OffsetCon, k+OffsetCon )
//
// Parameter   :  g_EC_Ele  : Array to store the output electric field
//                g_ConVar  : Array storing the input cell-centered conserved variables
//                g_FC_B    : Array storing the input face-centered magnetic field
//                NEle      : Stride for accessing g_EC_Ele[]
//                NCon      : Stride for accessing g_ConVar[]
//                OffsetCon : Offset for accessing g_ConVar[]
//
// Return      :  g_EC_Ele[]
//------------------------------------------------------------------------------------------------------
GPU_DEVICE
void MHD_ComputeElectric_Half(       real g_EC_Ele[][ CUBE(N_EC_ELE) ],
                               const real g_ConVar[][ CUBE(FLU_NXT) ],
                               const real g_FC_B  [][ SQR(FLU_NXT)*FLU_NXT_P1 ],
                               const int NEle, const int NCon, const int OffsetCon )
{

#  ifdef GAMER_DEBUG
   if ( NCon != FLU_NXT )
      printf( "ERROR : NCon != FLU_NXT (NCon %d, FLU_NXT %d)\n", NCon, FLU_NXT );

   if ( NEle != N_HF_ELE )
      printf( "ERROR : NEle != N_HF_ELE (NEle %d, N_HF_ELE %d)\n", NEle, N_HF_ELE );

   if ( NEle != NCon - 2*OffsetCon + 1 )
      printf( "ERROR : NEle != NCon - 2*OffsetCon + 1 (NEle %d, NCon %d, OffestCon %d)", NEle, NCon, OffsetCon );
#  endif // # ifdef GAMER_DEBUG


   const int NEleM1       = NEle - 1;
   const int NConP1       = NCon + 1;
   const int didx_con [3] = { 1, NCon, SQR(NCon) };

   for (int d=0; d<3; d++)
   {
      const int TDir1 = (d+1)%3;       // transverse direction 1
      const int TDir2 = (d+2)%3;       // transverse direction 2
      const int TV1   = TDir1 + 1;     // velocity component along the transverse direction 1
      const int TV2   = TDir2 + 1;     // velocity component along the transverse direction 2

      int idx_ele_e[2];
      int sizeB_i, sizeB_j, sizeB_k;   // B field stride along each direction

      switch ( d )
      {
         case 0 : idx_ele_e[0] = NEleM1;   idx_ele_e[1] = NEle;
                  sizeB_i      = NConP1;   sizeB_j      = NCon;     sizeB_k = NCon;
                  break;

         case 1 : idx_ele_e[0] = NEle;     idx_ele_e[1] = NEleM1;
                  sizeB_i      = NCon;     sizeB_j      = NConP1;   sizeB_k = NCon;
                  break;

         case 2 : idx_ele_e[0] = NEle;     idx_ele_e[1] = NEle;
                  sizeB_i      = NCon;     sizeB_j      = NCon;     sizeB_k = NConP1;
                  break;
      }

      const int size_ij = idx_ele_e[0]*idx_ele_e[1];

      CGPU_LOOP( idx, NEleM1*SQR(NEle) )
      {
         const int i_ele        = idx % idx_ele_e[0];
         const int j_ele        = idx % size_ij / idx_ele_e[0];
         const int k_ele        = idx / size_ij;
         const int idx_ele      = IDX321( i_ele, j_ele, k_ele, NEle, NEle );

         const int i_con        = i_ele + OffsetCon;
         const int j_con        = j_ele + OffsetCon;
         const int k_con        = k_ele + OffsetCon;
         const int idx_con      = IDX321( i_con, j_con, k_con, NCon, NCon );
         const int idx_con_m1   = idx_con - didx_con[TDir1];
         const int idx_con_m2   = idx_con - didx_con[TDir2];
         const int idx_con_m1m2 = idx_con - didx_con[TDir1] - didx_con[TDir2];

//       magnetic field index of the transverse direction 1
         const int i_B1   = i_con;
         const int j_B1   = j_con;
         const int k_B1   = k_con;
         const int idx_B1 = IDX321( i_B1, j_B1, k_B1, sizeB_k, sizeB_i );

//       magnetic field index of the transverse direction 2
         const int i_B2   = i_B1;
         const int j_B2   = j_B1;
         const int k_B2   = k_B1;
         const int idx_B2 = IDX321( i_B2, j_B2, k_B2, sizeB_j, sizeB_k );

         const int didx_B1[3] = { 1, sizeB_k, sizeB_k*sizeB_i };
         const int didx_B2[3] = { 1, sizeB_j, sizeB_j*sizeB_k };

         real V1, V2, B1, B2;

//       get the velocity average
         V1 = (real)0.25*( g_ConVar[TV1][idx_con     ] / g_ConVar[DENS][idx_con     ] +
                           g_ConVar[TV1][idx_con_m1  ] / g_ConVar[DENS][idx_con_m1  ] +
                           g_ConVar[TV1][idx_con_m2  ] / g_ConVar[DENS][idx_con_m2  ] +
                           g_ConVar[TV1][idx_con_m1m2] / g_ConVar[DENS][idx_con_m1m2] );

         V2 = (real)0.25*( g_ConVar[TV2][idx_con     ] / g_ConVar[DENS][idx_con     ] +
                           g_ConVar[TV2][idx_con_m1  ] / g_ConVar[DENS][idx_con_m1  ] +
                           g_ConVar[TV2][idx_con_m2  ] / g_ConVar[DENS][idx_con_m2  ] +
                           g_ConVar[TV2][idx_con_m1m2] / g_ConVar[DENS][idx_con_m1m2] );

//       get the magnetic field average
         B1 = (real)0.5*( g_FC_B[TDir1][ idx_B1 - didx_B1[TDir2] ] + g_FC_B[TDir1][ idx_B1 ] );
         B2 = (real)0.5*( g_FC_B[TDir2][ idx_B2 - didx_B2[TDir1] ] + g_FC_B[TDir2][ idx_B2 ] );

//       store the electric field of all cells in g_EC_Ele[]
         g_EC_Ele[d][idx_ele] = B1*V2 - B2*V1;

      } // CGPU_LOOP( idx, NEleM1*SQR(NEle) )
   } // for ( int d=0; d<3; d++)

#  ifdef __CUDACC__
   __syncthreads();
#  endif

} // FUNCTION : MHD_ComputeElectric_Half



//-------------------------------------------------------------------------------------------------------
// Function    :  MHD_UpdateMagnetic_Half
// Description :  Update magnetic field with the constrained transport algorithm
//
// Note        :  1. This function can only be used by MHM scheme
//                2. g_EC_Ele[] has the size of N_EC_ELE^3 but is accessed with a stride "NEle"
//                3. Ref : (a) Teyssier R., 2002, A&A, 385, 337. doi:10.1051/0004-6361:20011817
//                         (b) Fromang, S., Hennebelle, P., & Teyssier, R. 2006, aap, 457, 371. doi:10.1051/0004-6361:20065371
//
// Parameter   :  fc          : Array to store the updated face-centered B field and the other fluid variables
//                g_EC_Ele    : Array storing the input edge-centered electric field
//                dt          : Time interval to advance solution
//                dh          : Cell size
//                idx_{i,j,k} : Indices for accessing g_EC_Ele[]
//                NEle        : Stride for accessing g_EC_Ele[]
//
// Return      :  fc[][NCOMP_LR]
//------------------------------------------------------------------------------------------------------
GPU_DEVICE
void MHD_UpdateMagnetic_Half(       real fc[][NCOMP_LR],
                              const real g_EC_Ele[][ CUBE(N_EC_ELE) ],
                              const real dt, const real dh,
                              const int idx_i, const int idx_j, const int idx_k,
                              const int NEle )
{

   const real dt_dh2    = (real)0.5*dt/dh;
   const int  fL        = 0;
   const int  fR        = 1;
   const int  idx_E     = IDX321( idx_i, idx_j, idx_k, NEle, NEle );
   const int  didx_E[3] = { 1, NEle, SQR(NEle) };


// 1. calculate the source terms of magnetic field at all faces
   real B_source[2][3];

   for (int d=0; d<3; d++)
   {
      const int TDir1 = (d+1)%3;    // transverse direction 1
      const int TDir2 = (d+2)%3;    // transverse direction 2

      B_source[fL][d] = ( g_EC_Ele[TDir1][ idx_E + didx_E[TDir2] ] - g_EC_Ele[TDir1][ idx_E ] ) * dt_dh2 -
                        ( g_EC_Ele[TDir2][ idx_E + didx_E[TDir1] ] - g_EC_Ele[TDir2][ idx_E ] ) * dt_dh2 ;

      B_source[fR][d] = ( g_EC_Ele[TDir1][ idx_E + didx_E[TDir2] + didx_E[d] ] - g_EC_Ele[TDir1][ idx_E + didx_E[d] ] ) * dt_dh2 -
                        ( g_EC_Ele[TDir2][ idx_E + didx_E[TDir1] + didx_E[d] ] - g_EC_Ele[TDir2][ idx_E + didx_E[d] ] ) * dt_dh2 ;
   } // for (int d=0; d<3; d++)


// 2. update the magnetic field
//    --> face-centered normal B field is updated by the edge-centered EMF defined at the boundaries of the target face
   for (int d=0; d<3; d++)
   {
      const int faceL = d*2;
      const int faceR = faceL+1;
      const int TDir1 = (d+1)%3;    // transverse direction 1
      const int TDir2 = (d+2)%3;    // transverse direction 2

//    2-a. update the normal magnetic field
      fc[faceL][MAG_OFFSET+d] += B_source[fL][d];
      fc[faceR][MAG_OFFSET+d] += B_source[fR][d];

//    2-b. update the transverse magnetic field
//         --> source term of the transverse component should take the average along the transverse direction
      const real B_source_TDir1 = (real)0.5*( B_source[fL][TDir1] + B_source[fR][TDir1] );
      const real B_source_TDir2 = (real)0.5*( B_source[fL][TDir2] + B_source[fR][TDir2] );

      fc[faceL][MAG_OFFSET+TDir1] += B_source_TDir1;
      fc[faceR][MAG_OFFSET+TDir1] += B_source_TDir1;
      fc[faceL][MAG_OFFSET+TDir2] += B_source_TDir2;
      fc[faceR][MAG_OFFSET+TDir2] += B_source_TDir2;
   } // for (int d=0; d<3; d++)

} // FUNCTION : MHD_UpdateMagnetic_Half
#endif // # if ( FLU_SCHEME == MHM )



#endif // #if ( MODEL == HYDRO  &&  defined MHD )



#endif // #ifndef __CUFLU_CONSTRAINEDTRANSPORT__
