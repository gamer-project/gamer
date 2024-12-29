#ifndef __CUFLU_ADDCONDUCTIVEFLUX__
#define __CUFLU_ADDCONDUCTIVEFLUX__



#include "CUFLU.h"

#ifdef CONDUCTION



// external functions
#ifdef __CUDACC__

# include "CUFLU_ComputeConduction.cu"
# include "../CUFLU_Microphysics_SharedUtility.cu"

#else // #ifdef __CUDACC__

void Hydro_ComputeConduction( real &cond_kappa, real &cond_chi, const MicroPhy_t *MicroPhy,
                              const real Dens, const real Temp );
real MC_limiter( const real a, const real b );

#endif // #ifdef __CUDACC__ ... else ...


//-----------------------------------------------------------------------------------------
// Function    : Hydro_AddConductiveFlux_HalfStep
//
// Description : Compute the half-step conductive fluxes
//
// Note        : 1. Must enable CONDUCTION
//               2. Must enable MHD for anisotropic thermal conduction
//               3. Invoked by CPU/CUFLU_FluidSolver_MHM()
//
// Reference   :
//
// Parameter   : g_Con_Var   : Array storing the input cell-centered temperature
//               g_Flux_Half : Array with hydrodynamic fluxes for adding the conductive fluxes
//               g_FC_B      : Array storing the input face-centered B field
//               g_CC_B      : Array storing the input cell-centered B field
//               dh          : Cell size
//               MicroPhy    : Microphysics object
//
// Return      : g_Flux_Half[]
//-----------------------------------------------------------------------------------------
GPU_DEVICE
void Hydro_AddConductiveFlux_HalfStep( const real g_ConVar[][ CUBE(FLU_NXT) ],
                                       const real Temp[ CUBE(FLU_NXT) ],
                                             real g_Flux_Half[][NCOMP_TOTAL_PLUS_MAG][ CUBE(N_FC_FLUX) ],
                                       const real g_FC_B[][ SQR(FLU_NXT)*FLU_NXT_P1 ],
                                       const real g_CC_B[][ CUBE(FLU_NXT) ],
                                       const real dh, const MicroPhy_t *MicroPhy )
{

   const int  didx_cvar[3] = { 1, FLU_NXT, SQR(FLU_NXT) };
   const int  flux_offset  = 1;  // skip the additional fluxes along the transverse directions for computing the CT electric field
   const real _dh          = (real)1.0 / dh;

   for (int d=0; d<3; d++)
   {
      const int TDir1 = (d+1)%3;    // transverse direction 1
      const int TDir2 = (d+2)%3;    // transverse direction 2

      int sizeB_i, sizeB_j, stride_fc_B;
      int size_i, size_j, size_k;
      int i_offset, j_offset, k_offset;

      switch ( d )
      {
         case 0 : size_i   = N_HF_FLUX-1;              size_j   = N_HF_FLUX-2*flux_offset;  size_k      = N_HF_FLUX-2*flux_offset;
                  i_offset = 0;                        j_offset = flux_offset;              k_offset    = flux_offset;
#                 ifdef MHD
                  sizeB_i  = FLU_NXT_P1;               sizeB_j  = FLU_NXT;                  stride_fc_B = 1;
#                 endif
                  break;

         case 1 : size_i   = N_HF_FLUX-2*flux_offset;  size_j   = N_HF_FLUX-1;              size_k      = N_HF_FLUX-2*flux_offset;
                  i_offset = flux_offset;              j_offset = 0;                        k_offset    = flux_offset;
#                 ifdef MHD
                  sizeB_i  = FLU_NXT;                  sizeB_j  = FLU_NXT_P1;               stride_fc_B = FLU_NXT;
#                 endif
                  break;

         case 2 : size_i   = N_HF_FLUX-2*flux_offset;  size_j   = N_HF_FLUX-2*flux_offset;  size_k      = N_HF_FLUX-1;
                  i_offset = flux_offset;              j_offset = flux_offset;              k_offset    = 0;
#                 ifdef MHD
                  sizeB_i  = FLU_NXT;                  sizeB_j  = FLU_NXT;                  stride_fc_B = SQR(FLU_NXT);
#                 endif
                  break;
      } // switch ( d )

      const int size_ij = size_i*size_j;

      CGPU_LOOP( idx, size_i*size_j*size_k )
      {
//       flux index
         const int i_flux   = idx % size_i           + i_offset;
         const int j_flux   = idx % size_ij / size_i + j_offset;
         const int k_flux   = idx / size_ij          + k_offset;
         const int idx_flux = IDX321( i_flux, j_flux, k_flux, N_HF_FLUX, N_HF_FLUX );

//       conserved variable and cell-centered magnetic field index
         const int i_cvar   = i_flux;
         const int j_cvar   = j_flux;
         const int k_cvar   = k_flux;
         const int idx_cvar = IDX321( i_cvar, j_cvar, k_cvar, FLU_NXT, FLU_NXT );

#        ifdef MHD
//       face-centered magnetic field index
         const int idx_fc_B = IDX321( i_cvar, j_cvar, k_cvar, sizeB_i, sizeB_j ) + stride_fc_B;
#        endif

#        ifdef MHD
//       1. compute the mean magnetic field
//       ---------------------
//       |         |         |
//       |    ^    |    ^    |
//       -----1---------2-----
//       |    |    |    |    |
//       |         |         |
//       |  i j k -->        |
//       |         |         |
//       |    ^    |    ^    |
//       -----3---------4-----
//       |    |    |    |    |
//       |         |         |
//       ---------------------
         real B_N_mean, B_T1_mean, B_T2_mean, B_amp;
         if ( MicroPhy->CondFluxType == ANISOTROPIC_CONDUCTION )
         {
            B_N_mean  =             g_FC_B[    d][ idx_fc_B                ];
            B_T1_mean = (real)0.5*( g_CC_B[TDir1][ idx_cvar                ] +
                                    g_CC_B[TDir1][ idx_cvar + didx_cvar[d] ]   );
            B_T2_mean = (real)0.5*( g_CC_B[TDir2][ idx_cvar                ] +
                                    g_CC_B[TDir2][ idx_cvar + didx_cvar[d] ]   );
            B_amp     = SQRT( SQR(B_N_mean) + SQR(B_T1_mean) + SQR(B_T2_mean) );
//          normalize magnetic field
            B_N_mean  /= B_amp;
            B_T1_mean /= B_amp;
            B_T2_mean /= B_amp;
         }
#        endif // #ifdef MHD

//       2. compute temperature slope
//       ---------------------
//       |         |         |
//       ----bl--------br-----
//       |         |         |
//       |      N_slope      |
//       |         |         |
//       ----al--------ar-----
//       |         |         |
//       ---------------------
         real N_slope, T1_slope, T2_slope;
         real al, bl, ar, br;

//       normal direction
         N_slope = ( Temp[ idx_cvar + didx_cvar[d] ] - Temp[ idx_cvar ] ) * _dh;

         if ( MicroPhy->CondFluxType == ANISOTROPIC_CONDUCTION )
         {
//          transverse direction 1
            al = Temp[ idx_cvar                                   ] -
                 Temp[ idx_cvar                - didx_cvar[TDir1] ];
            bl = Temp[ idx_cvar                + didx_cvar[TDir1] ] -
                 Temp[ idx_cvar                                   ];
            ar = Temp[ idx_cvar + didx_cvar[d]                    ] -
                 Temp[ idx_cvar + didx_cvar[d] - didx_cvar[TDir1] ];
            br = Temp[ idx_cvar + didx_cvar[d] + didx_cvar[TDir1] ] -
                 Temp[ idx_cvar + didx_cvar[d]                    ];
            T1_slope = (  MC_limiter( MC_limiter(al,bl), MC_limiter(ar,br) )  ) * _dh;

//          transverse direction 2
            al = Temp[ idx_cvar                                   ] -
                 Temp[ idx_cvar                - didx_cvar[TDir2] ];
            bl = Temp[ idx_cvar                + didx_cvar[TDir2] ] -
                 Temp[ idx_cvar                                   ];
            ar = Temp[ idx_cvar + didx_cvar[d]                    ] -
                 Temp[ idx_cvar + didx_cvar[d] - didx_cvar[TDir2] ];
            br = Temp[ idx_cvar + didx_cvar[d] + didx_cvar[TDir2] ] -
                 Temp[ idx_cvar + didx_cvar[d]                    ];
            T2_slope = (  MC_limiter( MC_limiter(al,bl), MC_limiter(ar,br) )  ) * _dh;
         } else if ( CONDUCTION_SATURATION ) {
            T1_slope = ( Temp[ idx_cvar + didx_cvar[TDir1] ] -
                         Temp[ idx_cvar                    ] ) * _dh;
            T2_slope = ( Temp[ idx_cvar + didx_cvar[TDir2] ] -
                         Temp[ idx_cvar                    ] ) * _dh;
         }

//       3. compute conductive flux
         real gradient, gradT, Total_Flux;
#        ifdef MHD
         if ( MicroPhy->CondFluxType == ISOTROPIC_CONDUCTION )
            gradient = -N_slope;
         else if ( MicroPhy->CondFluxType == ANISOTROPIC_CONDUCTION )
            gradient = -B_N_mean*( B_N_mean*N_slope + B_T1_mean*T1_slope + B_T2_mean*T2_slope );
#        else
         gradient = -N_slope;
#        endif
//       get the conductivity
//       --> for non-constant conduction coefficients, take the spatial average along the normal direction
//           to get the face-centered coefficients
         real kappa_l, kappa_r, kappa, chi;
         Hydro_ComputeConduction( kappa_l, chi, MicroPhy, g_ConVar[DENS][ idx_cvar ],
                                  Temp[ idx_cvar ] );
         Hydro_ComputeConduction( kappa_r, chi, MicroPhy, g_ConVar[DENS][ idx_cvar + didx_cvar[d] ],
                                  Temp[ idx_cvar + didx_cvar[d] ] );
         kappa = 0.5*( kappa_l + kappa_r );
         if ( MicroPhy->CondSaturation ) {
            real gradT = SQRT( SQR(N_slope) + SQR(T1_slope) + SQR(T2_slope) );
            real Thalf = 0.5*( Temp[ idx_cvar ] + Temp[ idx_cvar + didx_cvar[d] ] );
            real Dhalf = 0.5*( g_ConVar[DENS][ idx_cvar ] + g_ConVar[DENS][ idx_cvar + didx_cvar[d] ] );
            real l_e = MicroPhy->CondMFPConst * Thalf * Thalf / Dhalf;
            real l_T = Thalf/gradT;
            real sat_factor = 4.0;
#           ifdef MHD
            if ( MicroPhy->CondSatWhistler )
            {
               real Phalf = MicroPhy->CondPresConv*Dhalf*Thalf;
               real beta = (real)2.0*Phalf/( B_amp * B_amp );
               sat_factor += beta/(real)3.0;
            }
#           endif
            kappa /= 1.0+sat_factor*l_e/l_T;
         }
         Total_Flux = kappa*gradient;

//       5. flux add-up
         g_Flux_Half[d][ENGY][idx_flux] += Total_Flux;

      } // CGPU_LOOP( idx, size_i*size_j*size_k )
   } // for (int d=0; d<3; d++)

#  ifdef __CUDACC__
   __syncthreads();
#  endif

} // FUNCTION : Hydro_AddConductiveFlux_HalfStep



//-----------------------------------------------------------------------------------------
// Function    : Hydro_AddConductiveFlux_FullStep
//
// Description : Compute the full-step conductive fluxes
//
// Note        : 1. Must enable CONDUCTION
//               2. Must enable MHD for anisotropic thermal conduction
//               3. Invoked by CPU/CUFLU_FluidSolver_MHM()
//
// Reference   : Yang et al., ApJ 761, 185 (2012); doi:10.1088/0004-637X/761/2/185
//
// Parameter   : g_PriVar_Half : Array storing the input cell-centered, half-step primitive fluid variables
//               g_FC_Flux     : Array with hydrodynamic fluxes for adding the conductive fluxes
//               g_FC_B_Half   : Array storing the input face-centered, half-step magnetic field
//               NFlux         : Stride for accessing g_FC_Flux[]
//               dh            : Cell size
//               MicroPhy      : Microphysics object
//
// Return      : g_FC_Flux[]
//-----------------------------------------------------------------------------------------
GPU_DEVICE
void Hydro_AddConductiveFlux_FullStep( const real g_PriVar_Half[][ CUBE(FLU_NXT) ],
                                       const real Temp[ CUBE(FLU_NXT) ],
                                             real g_FC_Flux[][NCOMP_TOTAL_PLUS_MAG][ CUBE(N_FC_FLUX) ],
                                       const real g_FC_B_Half[][ FLU_NXT_P1*SQR(FLU_NXT) ],
                                       const int NFlux, const real dh, const MicroPhy_t *MicroPhy )
{

#  ifdef MHD
   const int  mag_offset   = ( N_HF_VAR - PS2 ) / 2;
#  endif
   const int  didx_half[3] = { 1, N_HF_VAR, SQR(N_HF_VAR) };
   const int  half_offset  = ( N_HF_VAR - NFlux ) / 2;
   const int  flux_offset  = 1;  // skip the additional fluxes along the transverse directions for computing the CT electric field
   const real _dh          = (real)1.0 / dh;

   for (int d=0; d<3; d++)
   {
      const int TDir1 = (d+1)%3;    // transverse direction 1
      const int TDir2 = (d+2)%3;    // transverse direction 2

      int sizeB_i, sizeB_j, sizeB_k;
      int mag_offset_i, mag_offset_j, mag_offset_k;
      int size_i, size_j, size_k;
      int i_offset, j_offset, k_offset;

      switch ( d )
      {
         case 0 : size_i       = NFlux-1;               size_j        = NFlux-2*flux_offset;  size_k        = NFlux-2*flux_offset;
#                 ifdef MHD
                  sizeB_i      = N_HF_VAR+1;            sizeB_j       = N_HF_VAR;             sizeB_k       = N_HF_VAR;
                  mag_offset_i = mag_offset;            mag_offset_j  = mag_offset-1;         mag_offset_k  = mag_offset-1;
#                 endif
                  i_offset     = 0;                     j_offset      = flux_offset;          k_offset      = flux_offset;
                  break;

         case 1 : size_i        = NFlux-2*flux_offset;  size_j        = NFlux-1;              size_k        = NFlux-2*flux_offset;
#                 ifdef MHD
                  sizeB_i       = N_HF_VAR;             sizeB_j       = N_HF_VAR+1;           sizeB_k       = N_HF_VAR;
                  mag_offset_i  = mag_offset-1;         mag_offset_j  = mag_offset;           mag_offset_k  = mag_offset-1;
#                 endif
                  i_offset      = flux_offset;          j_offset      = 0;                    k_offset      = flux_offset;
                  break;

         case 2 : size_i        = NFlux-2*flux_offset;  size_j        = NFlux-2*flux_offset;  size_k        = NFlux-1;
#                 ifdef MHD
                  sizeB_i       = N_HF_VAR;             sizeB_j       = N_HF_VAR;             sizeB_k       = N_HF_VAR+1;
                  mag_offset_i  = mag_offset-1;         mag_offset_j  = mag_offset-1;         mag_offset_k  = mag_offset;
#                 endif
                  i_offset      = flux_offset;          j_offset      = flux_offset;          k_offset      = 0;
                  break;
      } // switch ( d )

#     ifdef MHD
      const int stride_fc_BT1[3] = { 1, sizeB_k, sizeB_k*sizeB_i };
      const int stride_fc_BT2[3] = { 1, sizeB_j, sizeB_j*sizeB_k };
#     endif
      const int size_ij          = size_i*size_j;

      CGPU_LOOP( idx, size_i*size_j*size_k )
      {
//       flux index
         const int i_flux     = idx % size_i           + i_offset;
         const int j_flux     = idx % size_ij / size_i + j_offset;
         const int k_flux     = idx / size_ij          + k_offset;
         const int idx_flux   = IDX321( i_flux, j_flux, k_flux, NFlux, NFlux );

//       half-step primitive variable index
         const int i_half     = i_flux + half_offset;
         const int j_half     = j_flux + half_offset;
         const int k_half     = k_flux + half_offset;
         const int idx_half   = IDX321( i_half, j_half, k_half, N_HF_VAR, N_HF_VAR );

#        ifdef MHD
//       magnetic field indices
         int i_fc, j_fc, k_fc, idx_fc_BN, idx_fc_BT1, idx_fc_BT2;
         if ( MicroPhy->CondFluxType == ANISOTROPIC_CONDUCTION )
         {
            i_fc       = i_flux + mag_offset_i;
            j_fc       = j_flux + mag_offset_j;
            k_fc       = k_flux + mag_offset_k;
            idx_fc_BN  = IDX321( i_fc, j_fc, k_fc, sizeB_i, sizeB_j );
            idx_fc_BT1 = IDX321( i_fc, j_fc, k_fc, sizeB_k, sizeB_i );
            idx_fc_BT2 = IDX321( i_fc, j_fc, k_fc, sizeB_j, sizeB_k );
         }
#        endif

#        ifdef MHD
//       1. compute the mean magnetic field
//       ---------------------
//       |         |         |
//       |    ^    |    ^    |
//       -----1---------2-----
//       |    |    |    |    |
//       |         |         |
//       |  i j k -->        |
//       |         |         |
//       |    ^    |    ^    |
//       -----3---------4-----
//       |    |    |    |    |
//       |         |         |
//       ---------------------
         real B_N_mean, B_T1_mean, B_T2_mean, B_amp;
         if ( MicroPhy->CondFluxType == ANISOTROPIC_CONDUCTION )
         {
            B_N_mean  =              g_FC_B_Half[    d][ idx_fc_BN                                            ];
            B_T1_mean = (real)0.25*( g_FC_B_Half[TDir1][ idx_fc_BT1                                           ] +
                                     g_FC_B_Half[TDir1][ idx_fc_BT1                    + stride_fc_BT1[TDir1] ] +
                                     g_FC_B_Half[TDir1][ idx_fc_BT1 - stride_fc_BT1[d]                        ] +
                                     g_FC_B_Half[TDir1][ idx_fc_BT1 - stride_fc_BT1[d] + stride_fc_BT1[TDir1] ]   );
            B_T2_mean = (real)0.25*( g_FC_B_Half[TDir2][ idx_fc_BT2                                           ] +
                                     g_FC_B_Half[TDir2][ idx_fc_BT2                    + stride_fc_BT2[TDir2] ] +
                                     g_FC_B_Half[TDir2][ idx_fc_BT2 - stride_fc_BT2[d]                        ] +
                                     g_FC_B_Half[TDir2][ idx_fc_BT2 - stride_fc_BT2[d] + stride_fc_BT2[TDir2] ]   );
            B_amp     = SQRT( SQR(B_N_mean) + SQR(B_T1_mean) + SQR(B_T2_mean) );

//          normalize magnetic field
            B_N_mean  /= B_amp;
            B_T1_mean /= B_amp;
            B_T2_mean /= B_amp;
         }
#        endif

//       2. compute temperature slope
//       ---------------------
//       |         |         |
//       ----bl--------br-----
//       |         |         |
//       |      N_slope      |
//       |         |         |
//       ----al--------ar-----
//       |         |         |
//       ---------------------
         real N_slope, T1_slope, T2_slope;
         real al, bl, ar, br;

//       normal direction
         N_slope = ( Temp[ idx_half + didx_half[d] ] - Temp[ idx_half ] ) * _dh;

         if ( MicroPhy->CondFluxType == ANISOTROPIC_CONDUCTION ) 
         {
//          transverse direction 1
            al = Temp[ idx_half                                   ] -
                 Temp[ idx_half                - didx_half[TDir1] ];
            bl = Temp[ idx_half                + didx_half[TDir1] ] -
                 Temp[ idx_half                                   ];
            ar = Temp[ idx_half + didx_half[d]                    ] -
                 Temp[ idx_half + didx_half[d] - didx_half[TDir1] ];
            br = Temp[ idx_half + didx_half[d] + didx_half[TDir1] ] -
                 Temp[ idx_half + didx_half[d]                    ];
            T1_slope = (  MC_limiter( MC_limiter(al,bl), MC_limiter(ar,br) )  ) * _dh;

//          transverse direction 2
            al = Temp[ idx_half                                   ] -
                 Temp[ idx_half                - didx_half[TDir2] ];
            bl = Temp[ idx_half                + didx_half[TDir2] ] -
                 Temp[ idx_half                                   ];
            ar = Temp[ idx_half + didx_half[d]                    ] -
                 Temp[ idx_half + didx_half[d] - didx_half[TDir2] ];
            br = Temp[ idx_half + didx_half[d] + didx_half[TDir2] ] -
                 Temp[ idx_half + didx_half[d]                    ];
            T2_slope = (  MC_limiter( MC_limiter(al,bl), MC_limiter(ar,br) )  ) * _dh;
         } else if ( CONDUCTION_SATURATION ) {
            T1_slope = ( Temp[ idx_half + didx_half[TDir1] ] -
                         Temp[ idx_half                    ] ) * _dh;
            T2_slope = ( Temp[ idx_half + didx_half[TDir2] ] -
                         Temp[ idx_half                    ] ) * _dh;
         }

//       3. compute conductive flux
         real gradient, Total_Flux;
#        ifdef MHD
         if ( MicroPhy->CondFluxType == ISOTROPIC_CONDUCTION )
            gradient = -N_slope;
         else if ( MicroPhy->CondFluxType == ANISOTROPIC_CONDUCTION )
            gradient = -B_N_mean*( B_N_mean*N_slope + B_T1_mean*T1_slope + B_T2_mean*T2_slope );
#        else
	 gradient = -N_slope;
#        endif
//       get the conductivity
//       --> for non-constant conduction coefficients, take the spatial average along the normal direction
//           to get the face-centered coefficients
         real kappa_l, kappa_r, kappa, chi;
         Hydro_ComputeConduction( kappa_l, chi, MicroPhy, g_PriVar_Half[DENS][ idx_half ],
                                  Temp[ idx_half ] );
         Hydro_ComputeConduction( kappa_r, chi, MicroPhy, g_PriVar_Half[DENS][ idx_half + didx_half[d] ],
                                  Temp[ idx_half + didx_half[d] ] );
         kappa = 0.5*( kappa_l + kappa_r );
         if ( MicroPhy->CondSaturation ) {
            real gradT = SQRT( SQR(N_slope) + SQR(T1_slope) + SQR(T2_slope) );
            real Thalf = 0.5*( Temp[ idx_half ] + Temp[ idx_half + didx_half[d] ] );
            real Dhalf = 0.5*( g_PriVar_Half[DENS][ idx_half ] + g_PriVar_Half[DENS][ idx_half + didx_half[d] ] );
            real l_e = MicroPhy->CondMFPConst * Thalf * Thalf / Dhalf;
            real l_T = Thalf/gradT;
            real sat_factor = 4.0;
#           ifdef MHD
            if ( MicroPhy->CondSatWhistler )
            {
               real Phalf = MicroPhy->CondPresConv*Dhalf*Thalf;
               real beta = (real)2.0*Phalf/( B_amp * B_amp );
               sat_factor += beta/(real)3.0;
            }
#           endif
            kappa /= 1.0+sat_factor*l_e/l_T;
         }

         Total_Flux = kappa*gradient;

//       4. flux add-up
         g_FC_Flux[d][ENGY][idx_flux] += Total_Flux;

      } // CGPU_LOOP( idx, size_i*size_j*size_k )
   } // for (int d=0; d<3; d++)

#  ifdef __CUDACC__
   __syncthreads();
#  endif

} // FUNCTION : Hydro_AddConductiveFlux_FullStep



//-----------------------------------------------------------------------------------------
// Function    : Hydro_AddConductiveFlux
//
// Description : Compute the conductive fluxes
//
// Note        : 1. Must enable CONDUCTION
//               2. Must enable MHD for anisotropic thermal conduction
//               3. Invoked by CPU/CUFLU_FluidSolver_MHM() and Hydro_DataReconstruction()
//
// Reference   :
//
// Parameter   : g_ConVar       : Array storing the input cell-centered temperature
//               g_PriVar       : Array storing the input cell-centered primitive fluid variables (NCOMP_TOTAL_PLUS_MAG)
//               g_Flux         : Array with hydrodynamic fluxes for adding the extra fluxes
//               g_FC_B         : Array storing the input face-centered B field
//               N_Var          : Size of g_ConVar/g_PriVar
//               N_Ghost        : Ghost zone size of data-reconstruction
//               N_Flux         : Size of flux
//               NSkip_N        : Empty size of flux on normal     direction
//               NSkip_T        : Empty size of flux on transverse direction
//               NSkip_MHM_Half : Skip size of g_ConVar for MHM scheme half-step only (TODO: need a better name)
//               dh             : Cell size
//               initialize     : initialize g_Flux to zero or not
//               MicroPhy       : Microphysics object
//
// Return      : g_Flux[], initialize
//-----------------------------------------------------------------------------------------
GPU_DEVICE
void Hydro_AddConductiveFlux( const real g_ConVar[][ CUBE(FLU_NXT) ],
                              const real g_PriVar[][ CUBE(FLU_NXT) ],
                              const real Temp[ CUBE(FLU_NXT) ],
                                    real g_Flux[][NCOMP_TOTAL_PLUS_MAG][ CUBE(N_FC_FLUX) ],
                              const real g_FC_B[][ SQR(FLU_NXT)*FLU_NXT_P1 ],
                              const int N_Var, const int N_Ghost, const int N_Flux, const int NSkip_N,
                              const int NSkip_T, const int NSkip_MHM_Half, const real dh,
                              bool &initialize, const MicroPhy_t *MicroPhy )
{
#  ifdef GAMER_DEBUG
   if ( g_ConVar == NULL  &&  g_PriVar == NULL )   Aux_Error( ERROR_INFO, "Both g_ConVar and g_PriVar are NULL!\n");

#  ifdef MHD
   if ( g_PriVar == NULL  &&  g_FC_B == NULL )   Aux_Error( ERROR_INFO, "Both g_PriVar and g_FC_B are NULL!\n");
#  endif
#  endif // #ifdef GAMER_DEBUG

   const int  didx_cvar[3] = { 1, N_Var, SQR(N_Var) };
   const real _dh          = (real)1.0 / dh;
// index skip update of flux arrray
#  ifdef MHD
   const int skip_update_T = 1;
#  else
   const int skip_update_T = 0;
#  endif

   for (int d=0; d<3; d++)
   {
      const int TDir1 = (d+1)%3;    // transverse direction 1
      const int TDir2 = (d+2)%3;    // transverse direction 2

      int flux_size_i, flux_size_j, flux_size_k; // flux need to be updated only, instead of N_Flux
      int flux_offset_i, flux_offset_j, flux_offset_k;
      int var_offset_i, var_offset_j, var_offset_k;

      switch ( d )
      {
         case 0 : flux_size_i   = N_Var-2*N_Ghost+2*NSkip_MHM_Half-NSkip_N-1;                 flux_size_j   = N_Var-2*N_Ghost+2*NSkip_MHM_Half-2*NSkip_T-2*skip_update_T; flux_size_k   = N_Var-2*N_Ghost+2*NSkip_MHM_Half-2*NSkip_T-2*skip_update_T;
                  flux_offset_i = 0;                                                          flux_offset_j = skip_update_T;                                              flux_offset_k = skip_update_T;
                  var_offset_i  = N_Ghost-NSkip_MHM_Half+NSkip_N;                             var_offset_j  = N_Ghost-NSkip_MHM_Half+NSkip_T;                             var_offset_k  = N_Ghost+NSkip_T-NSkip_MHM_Half+NSkip_T;
                  // not sure if var_offset_i works for CTU
                  break;

         case 1 : flux_size_i   = N_Var-2*N_Ghost+2*NSkip_MHM_Half-2*NSkip_T-2*skip_update_T; flux_size_j   = N_Var-2*N_Ghost+2*NSkip_MHM_Half-NSkip_N-1;                 flux_size_k   = N_Var-2*N_Ghost+2*NSkip_MHM_Half-2*NSkip_T-2*skip_update_T;
                  flux_offset_i = skip_update_T;                                              flux_offset_j = 0;                                                          flux_offset_k = skip_update_T;
                  var_offset_i  = N_Ghost-NSkip_MHM_Half+NSkip_T;                             var_offset_j  = N_Ghost-NSkip_MHM_Half+NSkip_N;                             var_offset_k  = N_Ghost-NSkip_MHM_Half+NSkip_T;
                  break;

         case 2 : flux_size_i   = N_Var-2*N_Ghost+2*NSkip_MHM_Half-2*NSkip_T-2*skip_update_T; flux_size_j   = N_Var-2*N_Ghost+2*NSkip_MHM_Half-2*NSkip_T-2*skip_update_T; flux_size_k   = N_Var-2*N_Ghost+2*NSkip_MHM_Half-NSkip_N-1;
                  flux_offset_i = skip_update_T;                                              flux_offset_j = skip_update_T;                                              flux_offset_k = 0;
                  var_offset_i  = N_Ghost-NSkip_MHM_Half+NSkip_T;                             var_offset_j  = N_Ghost-NSkip_MHM_Half+NSkip_T;                             var_offset_k  = N_Ghost-NSkip_MHM_Half+NSkip_N;
                  break;
      } // switch ( d )

#     ifdef MHD
      int sizeB_i, sizeB_j, sizeB_k;
      int mag_offset_i, mag_offset_j, mag_offset_k;
      switch ( d )
      {
         case 0 : sizeB_i      = N_Var+1;                  sizeB_j      = N_Var;                    sizeB_k      = N_Var;
                  mag_offset_i = N_Ghost-NSkip_MHM_Half+1; mag_offset_j = N_Ghost-NSkip_MHM_Half;   mag_offset_k = N_Ghost-NSkip_MHM_Half;
                  break;
         case 1 : sizeB_i      = N_Var;                    sizeB_j      = N_Var+1;                  sizeB_k      = N_Var;
                  mag_offset_i = N_Ghost-NSkip_MHM_Half;   mag_offset_j = N_Ghost-NSkip_MHM_Half+1; mag_offset_k = N_Ghost-NSkip_MHM_Half;
                  break;
         case 2 : sizeB_i      = N_Var;                    sizeB_j      = N_Var;                    sizeB_k      = N_Var+1;
                  mag_offset_i = N_Ghost-NSkip_MHM_Half;   mag_offset_j = N_Ghost-NSkip_MHM_Half;   mag_offset_k = N_Ghost-NSkip_MHM_Half+1;
                  break;
      } // switch ( d )

      const int stride_fc_BT1[3] = { 1, sizeB_k, sizeB_k*sizeB_i };
      const int stride_fc_BT2[3] = { 1, sizeB_j, sizeB_j*sizeB_k };
#     endif

      const int flux_size_ij = flux_size_i*flux_size_j;

      CGPU_LOOP( idx, flux_size_i*flux_size_j*flux_size_k )
      {
//       1. calcualte indices
//       flux index
//       --------------------
//       |        |         |
//       --------------------
//       |        |         |
//       |    -(i j k)->    |
//       |        |         |
//       --------------------
//       |        |         |
//       --------------------
         const int i_flux   = idx % flux_size_i                + flux_offset_i;
         const int j_flux   = idx % flux_size_ij / flux_size_i + flux_offset_j;
         const int k_flux   = idx / flux_size_ij               + flux_offset_k;
         const int idx_flux = IDX321( i_flux, j_flux, k_flux, N_Flux, N_Flux );

//       cell-centered conserved/primitive variable index
//       -----------------
//       |       |       |
//       -----------------
//       |       |       |
//       | i j k |       |
//       |       |       |
//       -----------------
//       |       |       |
//       -----------------
         const int i_cvar   = i_flux + var_offset_i;
         const int j_cvar   = j_flux + var_offset_j;
         const int k_cvar   = k_flux + var_offset_k;
         const int idx_cvar = IDX321( i_cvar, j_cvar, k_cvar, N_Var, N_Var );

#        ifdef MHD
//       magnetic field indices
//       --------------------
//       |        |         |
//       --------------------
//       |        |         |
//       |    -(i j k)->    |
//       |        |         |
//       --------------------
//       |        |         |
//       --------------------
         const int i_fc       = i_flux + mag_offset_i;
         const int j_fc       = j_flux + mag_offset_j;
         const int k_fc       = k_flux + mag_offset_k;
         const int idx_fc_BN  = IDX321( i_fc, j_fc, k_fc, sizeB_i, sizeB_j );
         const int idx_fc_BT1 = IDX321( i_fc, j_fc, k_fc, sizeB_k, sizeB_i );
         const int idx_fc_BT2 = IDX321( i_fc, j_fc, k_fc, sizeB_j, sizeB_k );
#        endif

//       2. compute the mean magnetic field at the face-centered flux location
#        ifdef MHD
         real B_N_mean, B_T1_mean, B_T2_mean, B_amp;
         if ( MicroPhy->CondFluxType == ANISOTROPIC_CONDUCTION )
         {
            if ( g_FC_B == NULL )
            {
//             MHM full-step only
               B_N_mean  = (real)0.5 * ( g_PriVar[NCOMP_TOTAL+d    ][idx_cvar] + g_PriVar[NCOMP_TOTAL+d    ][idx_cvar + didx_cvar[d]] );
               B_T1_mean = (real)0.5 * ( g_PriVar[NCOMP_TOTAL+TDir1][idx_cvar] + g_PriVar[NCOMP_TOTAL+TDir1][idx_cvar + didx_cvar[d]] );
               B_T2_mean = (real)0.5 * ( g_PriVar[NCOMP_TOTAL+TDir2][idx_cvar] + g_PriVar[NCOMP_TOTAL+TDir2][idx_cvar + didx_cvar[d]] );
            }
            else
            {
               B_N_mean  =              g_FC_B[    d][ idx_fc_BN                                            ];
               B_T1_mean = (real)0.25*( g_FC_B[TDir1][ idx_fc_BT1                                           ] +
                                        g_FC_B[TDir1][ idx_fc_BT1                    + stride_fc_BT1[TDir1] ] +
                                        g_FC_B[TDir1][ idx_fc_BT1 - stride_fc_BT1[d]                        ] +
                                        g_FC_B[TDir1][ idx_fc_BT1 - stride_fc_BT1[d] + stride_fc_BT1[TDir1] ]   );
               B_T2_mean = (real)0.25*( g_FC_B[TDir2][ idx_fc_BT2                                           ] +
                                        g_FC_B[TDir2][ idx_fc_BT2                    + stride_fc_BT2[TDir2] ] +
                                        g_FC_B[TDir2][ idx_fc_BT2 - stride_fc_BT2[d]                        ] +
                                        g_FC_B[TDir2][ idx_fc_BT2 - stride_fc_BT2[d] + stride_fc_BT2[TDir2] ]   );
            } // if ( g_FC_B == NULL )

            B_amp      = SQRT( SQR(B_N_mean) + SQR(B_T1_mean) + SQR(B_T2_mean) );
//          normalize magnetic field
            B_N_mean  /= B_amp;
            B_T1_mean /= B_amp;
            B_T2_mean /= B_amp;
         } // if ( MicroPhy->CondFluxType == ANISOTROPIC_CONDUCTION )
#        endif // #ifdef MHD

//       3. compute the mean temperature slope at the face-centered flux location
//       ---------------------
//       |         |         |
//       ----bl--------br-----
//       |         |         |
//       |      N_slope      |
//       |         |         |
//       ----al--------ar-----
//       |         |         |
//       ---------------------
         real N_slope, T1_slope, T2_slope; // normal, transverse 1, and transverse 2 direction
         real al, bl, ar, br; // temporary slope variables, see the graph above

         N_slope = ( Temp[ idx_cvar + didx_cvar[d] ] - Temp[ idx_cvar ] ) * _dh;

         if ( MicroPhy->CondFluxType == ANISOTROPIC_CONDUCTION )
         {
            al = Temp[ idx_cvar                                   ] -
                 Temp[ idx_cvar                - didx_cvar[TDir1] ];
            bl = Temp[ idx_cvar                + didx_cvar[TDir1] ] -
                 Temp[ idx_cvar                                   ];
            ar = Temp[ idx_cvar + didx_cvar[d]                    ] -
                 Temp[ idx_cvar + didx_cvar[d] - didx_cvar[TDir1] ];
            br = Temp[ idx_cvar + didx_cvar[d] + didx_cvar[TDir1] ] -
                 Temp[ idx_cvar + didx_cvar[d]                    ];
            T1_slope = (  MC_limiter( MC_limiter(al,bl), MC_limiter(ar,br) )  ) * _dh;

            al = Temp[ idx_cvar                                   ] -
                 Temp[ idx_cvar                - didx_cvar[TDir2] ];
            bl = Temp[ idx_cvar                + didx_cvar[TDir2] ] -
                 Temp[ idx_cvar                                   ];
            ar = Temp[ idx_cvar + didx_cvar[d]                    ] -
                 Temp[ idx_cvar + didx_cvar[d] - didx_cvar[TDir2] ];
            br = Temp[ idx_cvar + didx_cvar[d] + didx_cvar[TDir2] ] -
                 Temp[ idx_cvar + didx_cvar[d]                    ];
            T2_slope = (  MC_limiter( MC_limiter(al,bl), MC_limiter(ar,br) )  ) * _dh;
         }
         else if ( CONDUCTION_SATURATION )
         {
            T1_slope = ( Temp[ idx_cvar + didx_cvar[TDir1] ] -
                         Temp[ idx_cvar                    ] ) * _dh;
            T2_slope = ( Temp[ idx_cvar + didx_cvar[TDir2] ] -
                         Temp[ idx_cvar                    ] ) * _dh;
         }
         else
         {
            // TODO: Add error message here?
         } // if ( MicroPhy->CondFluxType == ANISOTROPIC_CONDUCTION ) ... else if ... else ...

//       3. compute conductive flux
         real gradient, gradT, Total_Flux;
#        ifdef MHD
         if      ( MicroPhy->CondFluxType ==   ISOTROPIC_CONDUCTION )
            gradient = -N_slope;
         else if ( MicroPhy->CondFluxType == ANISOTROPIC_CONDUCTION )
            gradient = -B_N_mean*( B_N_mean*N_slope + B_T1_mean*T1_slope + B_T2_mean*T2_slope );
         // TODO: Add error message for else condition here?
#        else
         gradient = -N_slope;
#        endif

//       get the conductivity
//       --> for non-constant conduction coefficients, take the spatial average along the normal direction
//           to get the face-centered coefficients
         real kappa_l, kappa_r, kappa, chi;
         real dens_L, dens_R, temp_L, temp_R;
         if ( g_PriVar == NULL )
         {
            dens_L = g_ConVar[DENS][ idx_cvar                ];
            dens_R = g_ConVar[DENS][ idx_cvar + didx_cvar[d] ];
         }
         else
         {
            dens_L = g_PriVar[DENS][ idx_cvar                ];
            dens_R = g_PriVar[DENS][ idx_cvar + didx_cvar[d] ];
         } // if ( g_PriVar == NULL ) ... else ...
         temp_L = Temp[ idx_cvar                ];
         temp_R = Temp[ idx_cvar + didx_cvar[d] ];

         Hydro_ComputeConduction( kappa_l, chi, MicroPhy, dens_L, temp_L );
         Hydro_ComputeConduction( kappa_r, chi, MicroPhy, dens_R, temp_R );
         kappa = 0.5*( kappa_l + kappa_r );

         if ( MicroPhy->CondSaturation )
         {
            real gradT      = SQRT( SQR(N_slope) + SQR(T1_slope) + SQR(T2_slope) );
            real Thalf      = 0.5 * ( temp_L + temp_R );
            real Dhalf      = 0.5 * ( dens_L + dens_R );
            real l_e        = MicroPhy->CondMFPConst * Thalf * Thalf / Dhalf;
            real l_T        = Thalf / gradT;
            real sat_factor = 4.0;
#           ifdef MHD
            if ( MicroPhy->CondSatWhistler )
            {
               real Phalf   = MicroPhy->CondPresConv * Dhalf * Thalf;
               real beta    = (real)2.0 * Phalf / ( B_amp * B_amp );
               sat_factor  += beta / (real)3.0;
            }
#           endif
            kappa /= 1.0 + sat_factor * l_e / l_T;
         } // if ( MicroPhy->CondSaturation )

         Total_Flux = kappa * gradient;

//       5. initialize flux if need
         if ( initialize )
         {
            for (int v=0; v<NCOMP_TOTAL_PLUS_MAG; v++)   g_Flux[d][v][idx_flux] = (real)0.0;
         }

//       6. flux add-up
         g_Flux[d][ENGY][idx_flux] += Total_Flux;

      } // CGPU_LOOP( idx, flux_size_i*flux_size_j*flux_size_k )
   } // for (int d=0; d<3; d++)

   if ( initialize ) initialize = false;

#  ifdef __CUDACC__
   __syncthreads();
#  endif

} // FUNCTION : Hydro_AddConductiveFlux



#endif // #ifdef CONDUCTION



#endif // #ifndef __CUFLU_ADDCONDUCTIVEFLUX__
