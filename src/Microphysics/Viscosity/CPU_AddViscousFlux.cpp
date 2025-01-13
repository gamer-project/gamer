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
real van_leer2( const real a, const real b, const real c, const real d );

#endif // #ifdef __CUDACC__ ... else ...

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
   const real two_thirds   = 2./3.;
#  ifdef MHD
   const int skip_update_T = 1;
#  else
   const int skip_update_T = 0;
#  endif

// pre-compute velocity for simplicity, depending on if we have conservative or primitive variables
   real Vel[3][ CUBE(FLU_NXT) ];
   
   for (int d=0; d<3; d++)
   {
      CGPU_LOOP( idx, CUBE(FLU_NXT) )
      {

         if ( g_PriVar == NULL ) 
            Vel[d][idx] = g_ConVar[d+1][idx] / g_ConVar[DENS][idx];
         else
            Vel[d][idx] = g_PriVar[d+1][idx];

      }

   }

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
//       1. calculate indices
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
//       These are only needed if viscosity is anisotropic
         if ( MicroPhy->ViscFluxType == ANISOTROPIC_VISCOSITY )
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
         } // if ( MicroPhy->ViscFluxType == ANISOTROPIC_VISCOSITY )
#        endif // #ifdef MHD

//       2. compute jacobian matrix for velocity gradients
         real v_N, v_T1, v_T2;
         real N_slope_N, T1_slope_N, T2_slope_N;
         real N_slope_T1, T1_slope_T1, T2_slope_T1;
         real N_slope_T2, T1_slope_T2, T2_slope_T2;
         real divV;

//       face-centered velocities
         v_N  = 0.5 * ( Vel[    d][ idx_cvar                ] + Vel[    d][ idx_cvar + didx_cvar[d] ] );
         v_T1 = 0.5 * ( Vel[TDir1][ idx_cvar                ] + Vel[TDir1][ idx_cvar + didx_cvar[d] ] );
         v_T2 = 0.5 * ( Vel[TDir2][ idx_cvar                ] + Vel[TDir2][ idx_cvar + didx_cvar[d] ] );
   
//       normal direction derivative
         N_slope_N  = ( Vel[    d][ idx_cvar + didx_cvar[d] ] - Vel[    d][ idx_cvar ] ) * _dh;
         T1_slope_N = ( Vel[TDir1][ idx_cvar + didx_cvar[d] ] - Vel[TDir1][ idx_cvar ] ) * _dh;
         T2_slope_N = ( Vel[TDir2][ idx_cvar + didx_cvar[d] ] - Vel[TDir2][ idx_cvar ] ) * _dh;

         if ( MicroPhy->ViscFluxType == ANISOTROPIC_VISCOSITY )
         {

//          transverse direction 1 derivative
            N_slope_T1  = van_leer2( 
               Vel[    d][ idx_cvar + didx_cvar[d] + didx_cvar[TDir1] ] - 
               Vel[    d][ idx_cvar + didx_cvar[d]                    ],
               Vel[    d][ idx_cvar + didx_cvar[d]                    ] - 
               Vel[    d][ idx_cvar + didx_cvar[d] - didx_cvar[TDir1] ],
               Vel[    d][ idx_cvar                + didx_cvar[TDir1] ] - 
               Vel[    d][ idx_cvar                                   ],
               Vel[    d][ idx_cvar                                   ] - 
               Vel[    d][ idx_cvar                - didx_cvar[TDir1] ]
            ) * _dh;
            T1_slope_T1 = van_leer2( 
               Vel[TDir1][ idx_cvar + didx_cvar[d] + didx_cvar[TDir1] ] - 
               Vel[TDir1][ idx_cvar + didx_cvar[d]                    ],
               Vel[TDir1][ idx_cvar + didx_cvar[d]                    ] - 
               Vel[TDir1][ idx_cvar + didx_cvar[d] - didx_cvar[TDir1] ],
               Vel[TDir1][ idx_cvar                + didx_cvar[TDir1] ] - 
               Vel[TDir1][ idx_cvar                                   ],
               Vel[TDir1][ idx_cvar                                   ] - 
               Vel[TDir1][ idx_cvar                - didx_cvar[TDir1] ]
            ) * _dh;
            T2_slope_T1  = van_leer2( 
               Vel[TDir2][ idx_cvar + didx_cvar[d] + didx_cvar[TDir1] ] - 
               Vel[TDir2][ idx_cvar + didx_cvar[d]                    ],
               Vel[TDir2][ idx_cvar + didx_cvar[d]                    ] - 
               Vel[TDir2][ idx_cvar + didx_cvar[d] - didx_cvar[TDir1] ],
               Vel[TDir2][ idx_cvar                + didx_cvar[TDir1] ] - 
               Vel[TDir2][ idx_cvar                                   ],
               Vel[TDir2][ idx_cvar                                   ] - 
               Vel[TDir2][ idx_cvar                - didx_cvar[TDir1] ]
            ) * _dh;

//          transverse direction 2 derivative
            N_slope_T2  = van_leer2( 
               Vel[    d][ idx_cvar + didx_cvar[d] + didx_cvar[TDir2] ] - 
               Vel[    d][ idx_cvar + didx_cvar[d]                    ],
               Vel[    d][ idx_cvar + didx_cvar[d]                    ] - 
               Vel[    d][ idx_cvar + didx_cvar[d] - didx_cvar[TDir2] ],
               Vel[    d][ idx_cvar                + didx_cvar[TDir2] ] - 
               Vel[    d][ idx_cvar                                   ],
               Vel[    d][ idx_cvar                                   ] - 
               Vel[    d][ idx_cvar                - didx_cvar[TDir2] ]
            ) * _dh;
            T1_slope_T2 = van_leer2( 
               Vel[TDir1][ idx_cvar + didx_cvar[d] + didx_cvar[TDir2] ] - 
               Vel[TDir1][ idx_cvar + didx_cvar[d]                    ],
               Vel[TDir1][ idx_cvar + didx_cvar[d]                    ] - 
               Vel[TDir1][ idx_cvar + didx_cvar[d] - didx_cvar[TDir2] ],
               Vel[TDir1][ idx_cvar                + didx_cvar[TDir2] ] - 
               Vel[TDir1][ idx_cvar                                   ],
               Vel[TDir1][ idx_cvar                                   ] - 
               Vel[TDir1][ idx_cvar                - didx_cvar[TDir2] ]
            ) * _dh;
            T2_slope_T2 = van_leer2( 
               Vel[TDir2][ idx_cvar + didx_cvar[d] + didx_cvar[TDir2] ] - 
               Vel[TDir2][ idx_cvar + didx_cvar[d]                    ],
               Vel[TDir2][ idx_cvar + didx_cvar[d]                    ] - 
               Vel[TDir2][ idx_cvar + didx_cvar[d] - didx_cvar[TDir2] ],
               Vel[TDir2][ idx_cvar                + didx_cvar[TDir2] ] - 
               Vel[TDir2][ idx_cvar                                   ],
               Vel[TDir2][ idx_cvar                                   ] - 
               Vel[TDir2][ idx_cvar                - didx_cvar[TDir2] ]
            ) * _dh;

//          divergence
            divV = N_slope_N + T1_slope_T1 + T2_slope_T2;     

         }
         else // ISOTROPIC_VISCOSITY
         {

//          transverse direction 1 and 2 derivatives

            N_slope_T1 = 0.25 * ( Vel[    d][ idx_cvar + didx_cvar[d] + didx_cvar[TDir1] ] +
                                  Vel[    d][ idx_cvar                + didx_cvar[TDir1] ] -
                                  Vel[    d][ idx_cvar + didx_cvar[d] - didx_cvar[TDir1] ] - 
                                  Vel[    d][ idx_cvar                - didx_cvar[TDir1] ] ) * _dh;
            N_slope_T2 = 0.25 * ( Vel[    d][ idx_cvar + didx_cvar[d] + didx_cvar[TDir2] ] +
                                  Vel[    d][ idx_cvar                + didx_cvar[TDir2] ] -
                                  Vel[    d][ idx_cvar + didx_cvar[d] - didx_cvar[TDir2] ] - 
                                  Vel[    d][ idx_cvar                - didx_cvar[TDir2] ] ) * _dh;

//          divergence
            divV = 0.25 * ( Vel[    d][ idx_cvar                    + didx_cvar[    d] ] -
                            Vel[    d][ idx_cvar                    - didx_cvar[    d] ] +
                            Vel[    d][ idx_cvar + didx_cvar[    d] + didx_cvar[    d] ] -
                            Vel[    d][ idx_cvar + didx_cvar[    d] - didx_cvar[    d] ] +
                            Vel[TDir1][ idx_cvar                    + didx_cvar[TDir1] ] -
                            Vel[TDir1][ idx_cvar                    - didx_cvar[TDir1] ] +
                            Vel[TDir1][ idx_cvar + didx_cvar[    d] + didx_cvar[TDir1] ] -
                            Vel[TDir1][ idx_cvar + didx_cvar[    d] - didx_cvar[TDir1] ] +
                            Vel[TDir2][ idx_cvar                    + didx_cvar[TDir2] ] - 
                            Vel[TDir2][ idx_cvar                    - didx_cvar[TDir2] ] +
                            Vel[TDir2][ idx_cvar + didx_cvar[    d] + didx_cvar[TDir2] ] - 
                            Vel[TDir2][ idx_cvar + didx_cvar[    d] - didx_cvar[TDir2] ] ) * _dh;

         } // if ( MicroPhy->ViscFluxType == ANISOTROPIC_VISCOSITY ) ... else ...

//       3. compute viscous flux
//       get the viscosity
//       --> for non-constant viscosity coefficients, take the spatial average along the normal direction
//           to get the face-centered coefficients
         real mu_l, mu_r, mu, visc_nu, delta_p;
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
         Hydro_ComputeViscosity( mu_r, visc_nu, MicroPhy, dens_R, temp_R );
         mu = 0.5*( mu_l + mu_r );

         real stress_N, stress_T1, stress_T2;
         if ( MicroPhy->ViscFluxType == ISOTROPIC_VISCOSITY )
         {

            stress_N  = -mu * ( 2.0*N_slope_N - two_thirds*divV );
            stress_T1 = -mu * ( T1_slope_N );
            stress_T2 = -mu * ( T2_slope_N );

         }
         else if ( MicroPhy->ViscFluxType == ANISOTROPIC_VISCOSITY )
         {

            real BBdV, delta_p;
            BBdV  = B_N_mean  * ( B_N_mean *  N_slope_N + B_T1_mean *  N_slope_T1 + B_T2_mean *  N_slope_T2 );
            BBdV += B_T1_mean * ( B_N_mean * T1_slope_N + B_T1_mean * T1_slope_T1 + B_T2_mean * T1_slope_T2 );
            BBdV += B_T2_mean * ( B_N_mean * T2_slope_N + B_T1_mean * T2_slope_T1 + B_T2_mean * T2_slope_T2 );
            delta_p = mu*(3.0*BBdV - divV);
            if ( MicroPhy->ViscFluxType == ANISOTROPIC_VISCOSITY && MicroPhy->ViscBounds )
               delta_p = FMIN( FMAX( delta_p, -B2 ), 0.5*B2 );
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
         g_Flux[d][    d+1][idx_flux] += stress_N;
         g_Flux[d][TDir1+1][idx_flux] += stress_T1;
         g_Flux[d][TDir2+1][idx_flux] += stress_T2;
#        ifndef BAROTROPIC_EOS
         g_Flux[d][ENGY][idx_flux] += stress_N*v_N + stress_T1*v_T1 + stress_T2*v_T2;
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
