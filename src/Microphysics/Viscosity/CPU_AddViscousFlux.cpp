#ifndef __CUFLU_ADDVISCOUSFLUX__
#define __CUFLU_ADDVISCOUSFLUX__


#include "CUFLU.h"

#ifdef VISCOSITY


// external functions
#ifdef __CUDACC__

# include "CUFLU_ComputeViscosity.cu"
# include "../CUFLU_Microphysics_SharedUtility.cu"

#else // #ifdef __CUDACC__

void Hydro_ComputeViscosity( real &visc_mu, real &visc_nu, const MicroPhy_t *MicroPhy,
                             const real Dens, const real Temp );
real MC_limiter( const real a, const real b );

#endif // #ifdef __CUDACC__ ... else ...


//-----------------------------------------------------------------------------------------
// Function    : Hydro_AddViscousFlux_HalfStep
//
// Description : Compute the half-step conductive fluxes
//
// Note        : 1. Must enable VISCOSITY
//               2. Must enable MHD for anisotropic (Braginskii) viscosity
//               3. Invoked by CPU/CUFLU_FluidSolver_MHM()
//
// Reference   :
//
// Parameter   : g_PriVar    : Array storing the input cell-centered primitive variables
//               g_Flux_Half : Array with hydrodynamic fluxes for adding the viscous fluxes
//               g_FC_B      : Array storing the input face-centered B field
//               g_CC_B      : Array storing the input cell-centered B field
//               dh          : Cell size
//               MicroPhy    : Microphysics object
//
// Return      : g_Flux_Half[]
//-----------------------------------------------------------------------------------------
GPU_DEVICE
void Hydro_AddViscousFlux_HalfStep( const real g_PriVar[][ CUBE(FLU_NXT) ],
                                    const real Temp[],
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
         real B_N_mean, B_T1_mean, B_T2_mean, B_amp, B2;
         if ( MicroPhy.ViscFluxType == ANISOTROPIC_VISCOSITY )
         {
            B_N_mean  =             g_FC_B[    d][ idx_fc_B                ];
            B_T1_mean = (real)0.5*( g_CC_B[TDir1][ idx_cvar                ] +
                                    g_CC_B[TDir1][ idx_cvar + didx_cvar[d] ]   );
            B_T2_mean = (real)0.5*( g_CC_B[TDir2][ idx_cvar                ] +
                                    g_CC_B[TDir2][ idx_cvar + didx_cvar[d] ]   );
            B2        = SQR(B_N_mean) + SQR(B_T1_mean) + SQR(B_T2_mean);
            B_amp     = SQRT( B2 );
//          normalize magnetic field
            B_N_mean  /= B_amp;
            B_T1_mean /= B_amp;
            B_T2_mean /= B_amp;
         }
#        endif // #ifdef MHD

//       2. compute jacobian matrix for velocity gradients
//       ---------------------
//       |         |         |
//       ----bl--------br-----
//       |         |         |
//       |      N_slope      |
//       |         |         |
//       ----al--------ar-----
//       |         |         |
//       ---------------------
         real v_N, v_T1, v_T2;
         real N_slope_N, T1_slope_N, T2_slope_N;
         real N_slope_T1, T1_slope_T1, T2_slope_T1;
         real N_slope_T2, T1_slope_T2, T2_slope_T2;
         real divV_l, divV_r;

//       face-centered velocities
         v_N  = 0.5*( g_PriVar[d+1]    [ idx_cvar ] + g_PriVar[d+1]    [ idx_cvar + didx_cvar[d] ] );
         v_T1 = 0.5*( g_PriVar[TDir1+1][ idx_cvar ] + g_PriVar[TDir1+1][ idx_cvar + didx_cvar[d] ] );
         v_T2 = 0.5*( g_PriVar[TDIr2+1][ idx_cvar ] + g_PriVar[TDir2+1][ idx_cvar + didx_cvar[d] ] );

//       normal direction
         N_slope_N  = ( g_PriVar[d+1]    [ idx_cvar + didx_cvar[d] ] - g_PriVar[d+1]    [ idx_cvar ] ) * _dh;
         T1_slope_N = ( g_PriVar[TDir1+1][ idx_cvar + didx_cvar[d] ] - g_PriVar[TDir1+1][ idx_cvar ] ) * _dh;
         T2_slope_N = ( g_PriVar[TDIr2+1][ idx_cvar + didx_cvar[d] ] - g_PriVar[TDir2+1][ idx_cvar ] ) * _dh;

         if ( MicroPhy.ViscFluxType == ANISOTROPIC_VISCOSITY )
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
         }
         else // ISOTROPIC_VISCOSITY
         {
            divV_l = N_slope_N + T1_slope + T2_slope;
            divV_r = N_slope_N + T1_slope + T2_slope;
         }

//       3. compute viscous flux
//       get the viscosity
//       --> for non-constant viscosity coefficients, take the spatial average along the normal direction
//           to get the face-centered coefficients
         real mu_l, mu_r, mu, visc_nu, delta_p, Total_Flux;

         Hydro_ComputeViscosity( mu_l, visc_nu, MicroPhy, g_PriVar[DENS][ idx_cvar ],
                                 Temp[ idx_cvar ] );
         Hydro_ComputeViscosity( mu_r, visc_nu, MicroPhy, g_PriVar[DENS][ idx_cvar + didx_cvar[d] ],
                                 Temp[ idx_cvar + didx_cvar[d] ] );
         mu = 0.5*( mu_l + mu_r );

         real stress_N, stress_T1, stress_T2;
         if ( MicroPhy.ViscFluxType == ISOTROPIC_VISCOSITY )
            stress_N  = ;
            stress_T1 = ;
            stress_T2 = ;
         else if ( MicroPhy.ViscFluxType == ANISOTROPIC_VISCOSITY )
            real BBdV, delta_p;
            delta_p = mu*(3.0*BBdV - divV);
            if ( MicroPhy.ViscFluxType == ANISOTROPIC_VISCOSITY && MicroPhy.ViscBounds )
               delta_p = FMIN( FMAX( delta_p, -B2 ) 0.5*B2 );
            stress_N  = -delta_p*(B_N_mean*B_N_mean - 1./3.);
            stress_T1 = -delta_p*B_T1_mean*B_N_mean;
            stress_T2 = -delta_p*B_T2_mean*B_N_mean;

//       5. flux add-up
         g_Flux_Half[d][d+1][idx_flux] += stress_N;
         g_Flux_Half[d][TDir1+1][idx_flux] += stress_T1;
         g_Flux_Half[d][TDir2+1][idx_flux] += stress_T2;
#        ifndef BAROTROPIC_EOS
         real Total_Flux = stress_N*v_N + stress_T1*v_T1 + stress_T2*v_T2;
         g_Flux_Half[d][ENGY][idx_flux] += Total_Flux;
#        endif


      } // CGPU_LOOP( idx, size_i*size_j*size_k )
   } // for (int d=0; d<3; d++)

#  ifdef __CUDACC__
   __syncthreads();
#  endif

} // FUNCTION : Hydro_AddViscousFlux_HalfStep



//-----------------------------------------------------------------------------------------
// Function    : Hydro_AddViscousFlux_FullStep
//
// Description : Compute the full-step viscosity fluxes
//
// Note        : 1. Must enable VISCOSITY
//               2. Must enable MHD for anisotropic (Braginskii) viscosity
//               3. Invoked by CPU/CUFLU_FluidSolver_MHM()
//
// Reference   :
//
// Parameter   : g_PriVar_Half : Array storing the input cell-centered, half-step primitive fluid variables
//               g_FC_Flux     : Array with hydrodynamic fluxes for adding the viscous fluxes
//               g_FC_B_Half   : Array storing the input face-centered, half-step magnetic field
//               NFlux         : Stride for accessing g_FC_Flux[]
//               dh            : Cell size
//               MicroPhy      : Microphysics object
//
// Return      : g_FC_Flux[]
//-----------------------------------------------------------------------------------------
GPU_DEVICE
void Hydro_AddViscousFlux_FullStep( const real g_PriVar_Half[][ CUBE(FLU_NXT) ],`
                                    const real Temp[],
                                          real g_FC_Flux[][NCOMP_TOTAL_PLUS_MAG][ CUBE(N_FC_FLUX) ],
                                    const real g_FC_B_Half[][ FLU_NXT_P1*SQR(FLU_NXT) ],
                                    const int NFlux, const real dh, const MicroPhy_t *MicroPhy )
{

   const int  didx_half[3] = { 1, N_HF_VAR, SQR(N_HF_VAR) };
   const int  mag_offset   = ( N_HF_VAR - PS2 ) / 2;
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
         if ( MicroPhy.CondFluxType == ANISOTROPIC_VISCOSITY )
         {
            const int i_fc       = i_flux + mag_offset_i;
            const int j_fc       = j_flux + mag_offset_j;
            const int k_fc       = k_flux + mag_offset_k;
            const int idx_fc_BN  = IDX321( i_fc, j_fc, k_fc, sizeB_i, sizeB_j );
            const int idx_fc_BT1 = IDX321( i_fc, j_fc, k_fc, sizeB_k, sizeB_i );
            const int idx_fc_BT2 = IDX321( i_fc, j_fc, k_fc, sizeB_j, sizeB_k );
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
         real B_N_mean, B_T1_mean, B_T2_mean, B_amp, B2;
         if ( MicroPhy.ViscFluxType == ANISOTROPIC_VISCOSITY )
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
            B2        = SQR(B_N_mean) + SQR(B_T1_mean) + SQR(B_T2_mean);
            B_amp     = SQRT( B2 );

//          normalize magnetic field
            B_N_mean  /= B_amp;
            B_T1_mean /= B_amp;
            B_T2_mean /= B_amp;
         }
#        endif

//       2. compute jacobian matrix for velocity gradients
//       ---------------------
//       |         |         |
//       ----bl--------br-----
//       |         |         |
//       |      N_slope      |
//       |         |         |
//       ----al--------ar-----
//       |         |         |
//       ---------------------
         real v_N, v_T1, v_T2;
         real N_slope_N, T1_slope_N, T2_slope_N;
         real N_slope_T1, T1_slope_T1, T2_slope_T1;
         real N_slope_T2, T1_slope_T2, T2_slope_T2;

//       face-centered velocities
         v_N  = 0.5*( g_PriVar_Half[d+1]    [ idx_half ] + g_PriVar_Half[d+1]    [ idx_half + didx_half[d] ] );
         v_T1 = 0.5*( g_PriVar_Half[TDir1+1][ idx_half ] + g_PriVar_Half[TDir1+1][ idx_half + didx_half[d] ] );
         v_T2 = 0.5*( g_PriVar_Half[TDir2+1][ idx_half ] + g_PriVar_Half[TDir2+1][ idx_half + didx_half[d] ] );

//       normal direction
         N_slope_N  = ( g_PriVar_Half[d+1]    [ idx_half + didx_cvar[d] ] - g_PriVar_Half[d+1]    [ idx_half ] ) * _dh;
         T1_slope_N = ( g_PriVar_Half[TDir1+1][ idx_half + didx_half[d] ] - g_PriVar_Half[TDir1+1][ idx_half ] ) * _dh;
         T2_slope_N = ( g_PriVar_Half[TDIr2+1][ idx_half + didx_half[d] ] - g_PriVar_Half[TDir2+1][ idx_half ] ) * _dh;

         if ( MicroPhy.ViscFluxType == ANISOTROPIC_VISCOSITY )
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
         }
         else // ISOTROPIC_VISCOSITY
         {

         }

//       3. compute viscous flux
//       get the viscosity
//       --> for non-constant viscosity coefficients, take the spatial average along the normal direction
//           to get the face-centered coefficients
         real mu_l, mu_r, mu, visc_nu;
         Hydro_ComputeViscosity( mu_l, visc_nu, MicroPhy, g_PriVar_Half[DENS][ idx_half ],
                                 Temp[ idx_half ] );
         Hydro_ComputeViscosity( mu_r, visc_nu, MicroPhy, g_PriVar_Half[DENS][ idx_half + didx_half[d] ],
                                 Temp[ idx_half + didx_half[d] ] );
         mu = 0.5*( mu_l + mu_r );
         real stress_N, stress_T1, stress_T2;
         if ( MicroPhy.ViscFluxType == ISOTROPIC_VISCOSITY )
            stress_N  = ;
            stress_T1 = ;
            stress_T2 = ;
         else if ( MicroPhy.ViscFluxType == ANISOTROPIC_VISCOSITY )
            real BBdV, delta_p;
            delta_p = mu*(3.0*BBdV - divV);
            if ( MicroPhy.ViscFluxType == ANISOTROPIC_VISCOSITY && MicroPhy.ViscBounds )
               delta_p = FMIN( FMAX( delta_p, -B2 ) 0.5*B2 );
            stress_N  = -delta_p*(B_N_mean*B_N_mean - 1./3.);
            stress_T1 = -delta_p*B_T1_mean*B_N_mean;
            stress_T2 = -delta_p*B_T2_mean*B_N_mean;

//       4. flux add-up
         g_FC_Flux[d][d+1][idx_flux] += stress_N;
         g_FC_Flux[d][TDir1+1][idx_flux] += stress_T1;
         g_FC_Flux[d][TDir2+1][idx_flux] += stress_T2;
#        ifndef BAROTROPIC_EOS
         real Total_Flux = stress_N*v_N + stress_T1*v_T1 + stress_T2*v_T2;
         g_FC_Flux[d][ENGY][idx_flux] += Total_Flux;
#        endif

      } // CGPU_LOOP( idx, size_i*size_j*size_k )
   } // for (int d=0; d<3; d++)

#  ifdef __CUDACC__
   __syncthreads();
#  endif

} // FUNCTION : Hydro_AddViscousFlux_FullStep



//-----------------------------------------------------------------------------------------
// Function    : Hydro_AddViscousFlux
//
// Description : Compute the viscousity fluxes
//
// Note        : 1. Must enable VISCOSITY
//               2. Must enable MHD for anisotropic (Braginskii) viscosity
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
// Return      : g_Flux, initialize
//-----------------------------------------------------------------------------------------
GPU_DEVICE
void Hydro_AddViscousFlux( const real g_ConVar[][ CUBE(FLU_NXT) ],
                           const real g_PriVar[][ CUBE(FLU_NXT) ],
                           const real Temp[],
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
         real B_N_mean, B_T1_mean, B_T2_mean, B_amp, B2;
         if ( MicroPhy.ViscFluxType == ANISOTROPIC_VISCOSITY )
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
            } // if ( g_FC_B == NULL ) ...  else ...

            B2        = SQR(B_N_mean) + SQR(B_T1_mean) + SQR(B_T2_mean);
            B_amp     = SQRT( B2 );
//          normalize magnetic field
            B_N_mean  /= B_amp;
            B_T1_mean /= B_amp;
            B_T2_mean /= B_amp;
         } // if ( MicroPhy.ViscFluxType == ANISOTROPIC_VISCOSITY )
#        endif // #ifdef MHD

//       2. compute jacobian matrix for velocity gradients
//       ---------------------
//       |         |         |
//       ----bl--------br-----
//       |         |         |
//       |      N_slope      |
//       |         |         |
//       ----al--------ar-----
//       |         |         |
//       ---------------------
         real v_N, v_T1, v_T2;
         real N_slope_N, T1_slope_N, T2_slope_N;
         real N_slope_T1, T1_slope_T1, T2_slope_T1;
         real N_slope_T2, T1_slope_T2, T2_slope_T2;
         real divV_l, divV_r;

//       face-centered velocities
         if ( g_PriVar == NULL )
         {
            v_N  = 0.5*( g_ConVar[    d+1][ idx_cvar                ] / g_ConVar[DENS][ idx_cvar                ] +
                         g_ConVar[    d+1][ idx_cvar + didx_cvar[d] ] / g_ConVar[DENS][ idx_cvar + didx_cvar[d] ] );
            v_T1 = 0.5*( g_ConVar[TDir1+1][ idx_cvar                ] / g_ConVar[DENS][ idx_cvar                ] +
                         g_ConVar[TDir1+1][ idx_cvar + didx_cvar[d] ] / g_ConVar[DENS][ idx_cvar + didx_cvar[d] ] );
            v_T2 = 0.5*( g_ConVar[TDIr2+1][ idx_cvar                ] / g_ConVar[DENS][ idx_cvar                ] +
                         g_ConVar[TDir2+1][ idx_cvar + didx_cvar[d] ] / g_ConVar[DENS][ idx_cvar + didx_cvar[d] ] );

//          normal direction
            N_slope_N  = ( g_ConVar[    d+1][ idx_cvar + didx_cvar[d] ] / g_ConVar[DENS][ idx_cvar + didx_cvar[d] ] -
                           g_ConVar[    d+1][ idx_cvar                ] / g_ConVar[DENS][ idx_cvar                ] ) * _dh;
            T1_slope_N = ( g_ConVar[TDir1+1][ idx_cvar + didx_cvar[d] ] / g_ConVar[DENS][ idx_cvar + didx_cvar[d] ] -
                           g_ConVar[TDir1+1][ idx_cvar                ] / g_ConVar[DENS][ idx_cvar                ] ) * _dh;
            T2_slope_N = ( g_ConVar[TDIr2+1][ idx_cvar + didx_cvar[d] ] / g_ConVar[DENS][ idx_cvar + didx_cvar[d] ] -
                           g_ConVar[TDir2+1][ idx_cvar                ] / g_ConVar[DENS][ idx_cvar                ] ) * _dh;
         }
         else
         {
            v_N  = 0.5*( g_PriVar[    d+1][ idx_cvar ] + g_PriVar[    d+1][ idx_cvar + didx_cvar[d] ] );
            v_T1 = 0.5*( g_PriVar[TDir1+1][ idx_cvar ] + g_PriVar[TDir1+1][ idx_cvar + didx_cvar[d] ] );
            v_T2 = 0.5*( g_PriVar[TDIr2+1][ idx_cvar ] + g_PriVar[TDir2+1][ idx_cvar + didx_cvar[d] ] );

//          normal direction
            N_slope_N  = ( g_PriVar[d+1]    [ idx_cvar + didx_cvar[d] ] - g_PriVar[d+1]    [ idx_cvar ] ) * _dh;
            T1_slope_N = ( g_PriVar[TDir1+1][ idx_cvar + didx_cvar[d] ] - g_PriVar[TDir1+1][ idx_cvar ] ) * _dh;
            T2_slope_N = ( g_PriVar[TDIr2+1][ idx_cvar + didx_cvar[d] ] - g_PriVar[TDir2+1][ idx_cvar ] ) * _dh;
         } // if ( g_PriVar == NULL ) ... else ...

         if ( MicroPhy.ViscFluxType == ANISOTROPIC_VISCOSITY )
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
         }
         else // ISOTROPIC_VISCOSITY
         {
            divV_l = N_slope_N + T1_slope + T2_slope;
            divV_r = N_slope_N + T1_slope + T2_slope;
         } // if ( MicroPhy.ViscFluxType == ANISOTROPIC_VISCOSITY ) ... else ...

//       3. compute viscous flux
//       get the viscosity
//       --> for non-constant viscosity coefficients, take the spatial average along the normal direction
//           to get the face-centered coefficients
         real mu_l, mu_r, mu, visc_nu, delta_p, Total_Flux;
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

         Hydro_ComputeViscosity( mu_l, visc_nu, MicroPhy, dens_L, temp_L );
         Hydro_ComputeViscosity( mu_r, visc_nu, MicroPhy, dens_R, temp_L );
         mu = 0.5*( mu_l + mu_r );

         real stress_N, stress_T1, stress_T2;
         if ( MicroPhy.ViscFluxType == ISOTROPIC_VISCOSITY )
         {
            stress_N  = ;
            stress_T1 = ;
            stress_T2 = ;
         }
         else if ( MicroPhy.ViscFluxType == ANISOTROPIC_VISCOSITY )
         {
            real BBdV, delta_p;
            delta_p = mu*(3.0*BBdV - divV);
            if ( MicroPhy.ViscFluxType == ANISOTROPIC_VISCOSITY && MicroPhy.ViscBounds )
               delta_p = FMIN( FMAX( delta_p, -B2 ) 0.5*B2 );
            stress_N  = -delta_p*(B_N_mean*B_N_mean - 1./3.);
            stress_T1 = -delta_p*B_T1_mean*B_N_mean;
            stress_T2 = -delta_p*B_T2_mean*B_N_mean;
         }
         else
         {
            // TODO: some error message?
         }

//       5. initialize flux if need
         if ( initialize )
         {
            for (int v=0; v<NCOMP_TOTAL_PLUS_MAG; v++)   g_Flux[d][v][idx_flux] = (real)0.0;
         }

//       6. flux add-up
         g_Flux_Half[d][    d+1][idx_flux] += stress_N;
         g_Flux_Half[d][TDir1+1][idx_flux] += stress_T1;
         g_Flux_Half[d][TDir2+1][idx_flux] += stress_T2;
#        ifndef BAROTROPIC_EOS
         real Total_Flux = stress_N*v_N + stress_T1*v_T1 + stress_T2*v_T2;
         g_Flux_Half[d][ENGY][idx_flux] += Total_Flux;
#        endif

      } // CGPU_LOOP( idx, flux_size_i*flux_size_j*flux_size_k )
   } // for (int d=0; d<3; d++)
   
   if ( initialize ) initialize = false;

#  ifdef __CUDACC__
   __syncthreads();
#  endif

} // FUNCTION : Hydro_AddViscousFlux



#endif // #ifdef VISCOSITY

#endif // #ifndef __CUFLU_ADDVISCOUSFLUX__
