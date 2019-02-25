#ifndef __CUFLU_CONSTRAINEDTRANSPORT__
#define __CUFLU_CONSTRAINEDTRANSPORT__



#include "CUFLU.h"

#if ( MODEL == HYDRO  &&  defined MHD )



// external functions
#ifdef __CUDACC__

#include "CUFLU_Shared_FluUtility.cu"

#else

void MHD_GetCellCenteredB( real B_CC[], const real Bx_FC[], const real By_FC[], const real Bz_FC[],
                           const int Width_FC, const int i, const int j, const int k );

#endif // #ifdef __CUDACC__ ... else ...


// internal functions
GPU_DEVICE
static real dE_Upwind( const real FC_Ele_L, const real FC_Ele_R, const real FC_Mom, const real D_L, const real D_R,
                       const real V_L1, const real V_L2, const real V_R1, const real V_R2,
                       const real B_L1, const real B_L2, const real B_R1, const real B_R2,
                       const real dt_dh );




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
      const int TV1   = TDir1 + 1;           // velocity component along the transverse direction 1
      const int TV2   = TDir2 + 1;           // velocity component along the transverse direction 2
      const int TB1   = TDir1 + MAG_OFFSET;  // B flux   component along the transverse direction 1
      const int TB2   = TDir2 + MAG_OFFSET;  // B flux   component along the transverse direction 2

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

         real D_L, D_R, V_L1, V_L2, V_R1, V_R2, B_L1, B_L2, B_R1, B_R2;
         int  idx_L, idx_R;

         g_EC_Ele[d][idx_ele] = ( - g_FC_Flux[TDir1][TB2][ idx_flux + didx_flux[TDir2] ]
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

         g_EC_Ele[d][idx_ele] += dE_Upwind( -g_FC_Flux[TDir1][TB2][ idx_flux                    ],
                                            -g_FC_Flux[TDir1][TB2][ idx_flux + didx_flux[TDir2] ],
                                             g_FC_Flux[TDir2][  0][ idx_flux                    ],
                                             D_L, D_R, V_L1, V_L2, V_R1, V_R2, B_L1, B_L2, B_R1, B_R2,
                                             dt_dh );

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

         g_EC_Ele[d][idx_ele] += dE_Upwind( -g_FC_Flux[TDir1][TB2][ idx_flux                    ],
                                            -g_FC_Flux[TDir1][TB2][ idx_flux + didx_flux[TDir2] ],
                                             g_FC_Flux[TDir2][  0][ idx_flux + didx_flux[TDir1] ],
                                             D_L, D_R, V_L1, V_L2, V_R1, V_R2, B_L1, B_L2, B_R1, B_R2,
                                             dt_dh );

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

         g_EC_Ele[d][idx_ele] += dE_Upwind( +g_FC_Flux[TDir2][TB1][ idx_flux                    ],
                                            +g_FC_Flux[TDir2][TB1][ idx_flux + didx_flux[TDir1] ],
                                             g_FC_Flux[TDir1][  0][ idx_flux                    ],
                                             D_L, D_R, V_L1, V_L2, V_R1, V_R2, B_L1, B_L2, B_R1, B_R2,
                                             dt_dh );

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

         g_EC_Ele[d][idx_ele] += dE_Upwind( +g_FC_Flux[TDir2][TB1][ idx_flux                    ],
                                            +g_FC_Flux[TDir2][TB1][ idx_flux + didx_flux[TDir1] ],
                                             g_FC_Flux[TDir1][  0][ idx_flux + didx_flux[TDir2] ],
                                             D_L, D_R, V_L1, V_L2, V_R1, V_R2, B_L1, B_L2, B_R1, B_R2,
                                             dt_dh );

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

// MAX_ERROR is defined in CUFLU.h
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
//                   g_Flux[]       has the size of N_FC_FLUX^3 and is accessed with the same stride
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

   const int  didx_flux[3] = { 1, N_FC_FLUX, SQR(N_FC_FLUX) };
   const real dt_dh2       = (real)0.5*dt/dh;
   const int  NFluVar      = NCOMP_FLUID - 1;   // density + momentum*3

   real dFlux[3][NFluVar], Output_1Cell[ NFluVar + NCOMP_MAG ];

   const int size_ij = SQR(N_HF_VAR);
   CGPU_LOOP( idx_out, CUBE(N_HF_VAR) )
   {
      const int i_out      = idx_out % N_HF_VAR;
      const int j_out      = idx_out % size_ij / N_HF_VAR;
      const int k_out      = idx_out / size_ij;
      const int idx_flux   = IDX321( i_out, j_out, k_out, N_FC_FLUX, N_FC_FLUX );

      const int i_flu_in   = i_out + FLU_GHOST_SIZE - 1;    // assuming N_HF_VAR = PS2+2
      const int j_flu_in   = j_out + FLU_GHOST_SIZE - 1;
      const int k_flu_in   = k_out + FLU_GHOST_SIZE - 1;
      const int idx_flu_in = IDX321( i_flu_in, j_flu_in, k_flu_in, FLU_NXT, FLU_NXT );


//    1. calculate flux difference to update the fluid data by 0.5*dt
      for (int d=0; d<3; d++)
      for (int v=0; v<NFluVar; v++)
         dFlux[d][v] = g_Flux[d][v][ idx_flux + didx_flux[d] ] - g_Flux[d][v][idx_flux];

      for (int v=0; v<NFluVar; v++)
         Output_1Cell[v] = g_Flu_In[v][idx_flu_in] - dt_dh2*( dFlux[0][v] + dFlux[1][v] + dFlux[2][v] );

//    apply density floor
      Output_1Cell[DENS] = FMAX( Output_1Cell[DENS], MinDens );

//    check negative density
#     ifdef CHECK_NEGATIVE_IN_FLUID
      if ( Hydro_CheckNegative(Output_1Cell[DENS]) )
         printf( "WARNING : invalid density (%14.7e) at file <%s>, line <%d>, function <%s>\n",
                 Output_1Cell[DENS], __FILE__, __LINE__, __FUNCTION__ );
#     endif


//    2. momentum --> velocity
      const real _Dens = (real)1.0 / Output_1Cell[DENS];
      for (int v=1; v<NFluVar; v++)    Output_1Cell[v] *= _Dens;


//    3. compute the cell-centered half-step B field
      MHD_GetCellCenteredB( Output_1Cell+NFluVar, g_FC_B_Half[0], g_FC_B_Half[1], g_FC_B_Half[2],
                            N_HF_VAR, i_out, j_out, k_out );


//    4. store results to the output array
//       --> variable indices in g_PriVar_Out[] remain consistent with other arrays even though
//           energy and passive scalars have not been skipped here
      for (int v=0; v<NFluVar; v++)    g_PriVar_Out[ v              ][idx_out] = Output_1Cell[ v           ];
      for (int v=0; v<NCOMP_MAG; v++)  g_PriVar_Out[ v + MAG_OFFSET ][idx_out] = Output_1Cell[ v + NFluVar ];

   } // CGPU_LOOP( idx_out, CUBE(N_HF_VAR) )

} // FUNCTION : MHD_HalfStepPrimitive



#endif // #if ( MODEL == HYDRO  &&  defined MHD )



#endif // #ifndef __CUFLU_CONSTRAINEDTRANSPORT__
