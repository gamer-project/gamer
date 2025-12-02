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
real van_leer2( const real a, const real b, const real c, const real d );
real compute_temperature( const real ConVar[][ CUBE(FLU_NXT) ],
                          const real PriVar[][ CUBE(FLU_NXT) ],
                          const real   FC_B[][ SQR(FLU_NXT)*FLU_NXT_P1 ],
                          const int idx, const real MinTemp,
                          const long PassiveFloor, const EoS_t *EoS );

#endif // #ifdef __CUDACC__ ... else ...

GPU_DEVICE
real get_vel_op( const real ConVar[][ CUBE(FLU_NXT) ],
                 const real PriVar[][ CUBE(FLU_NXT) ],
                 const int comp,  // 0,1,2 for vx,vy,vz
                 const int idxL, const int idxR, const int op )
{
   real v_L, v_R, result;

   if ( PriVar == NULL )
   {
      // conservative -> primitive
      v_L = ConVar[comp+1][idxL] / ConVar[DENS][idxL];
      v_R = ConVar[comp+1][idxR] / ConVar[DENS][idxR];
   }
   else
   {
      // primitive already
      v_L = PriVar[comp+1][idxL];
      v_R = PriVar[comp+1][idxR];
   }

   if ( op == 0 )
      result = 0.5*( v_L + v_R );
   else
      result = v_R - v_L;

   return result;

} // FUNCTION: get_vel_op


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
// Parameter   : g_ConVar       : Array storing the input cell-centered conserved variables (NCOMP_TOTAL_PLUS_MAG)
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
//               MinTemp        : Temperature floor
//               PassiveFloor   : Bitwise flag to specify the passive scalars to be floored
//               EoS            : Equation of state object
//               MicroPhy       : Microphysics object
//
// Return      : g_Flux, initialize
//-----------------------------------------------------------------------------------------
GPU_DEVICE
void Hydro_AddViscousFlux( const real g_ConVar[][ CUBE(FLU_NXT) ],
                           const real g_PriVar[][ CUBE(FLU_NXT) ],
                                 real g_Flux[][NCOMP_TOTAL_PLUS_MAG][ CUBE(N_FC_FLUX) ],
                           const real g_FC_B[][ SQR(FLU_NXT)*FLU_NXT_P1 ],
                           const int N_Var, const int N_Ghost, const int N_Flux,
                           const int NSkip_N, const int NSkip_T, const int NSkip_MHM_Half,
                           const real dh, bool &initialize, const real MinTemp,
                           const long PassiveFloor, const EoS_t *EoS, const MicroPhy_t *MicroPhy )
{
#  ifdef GAMER_DEBUG
   if ( g_ConVar == NULL  &&  g_PriVar == NULL )   Aux_Error( ERROR_INFO, "Both g_ConVar and g_PriVar are NULL!\n");

#  ifdef MHD
   if ( g_PriVar == NULL  &&  g_FC_B == NULL )   Aux_Error( ERROR_INFO, "Both g_PriVar and g_FC_B are NULL!\n");
#  endif
#  endif // #ifdef GAMER_DEBUG

   const int  didx_cvar[3] = { 1, N_Var, SQR(N_Var) };
   const real _dh          = (real)1.0 / dh;
   const real one_third    = (real)(1./3.);
   const real two_thirds   = (real)(2./3.);
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
         const int i_cvar      = i_flux + var_offset_i;
         const int j_cvar      = j_flux + var_offset_j;
         const int k_cvar      = k_flux + var_offset_k;
         const int idx_cvar    = IDX321( i_cvar, j_cvar, k_cvar, N_Var, N_Var );
         const int idx_cvar_dd = idx_cvar + didx_cvar[d];

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

//       2. compute the mean magnetic field at the face-centered flux location
         real B_N_mean, B_T1_mean, B_T2_mean, B_amp, B2;
         if ( g_FC_B == NULL )
         {
//          MHM full-step only
            B_N_mean  = (real)0.5 * ( g_PriVar[NCOMP_TOTAL+d    ][idx_cvar] + g_PriVar[NCOMP_TOTAL+d    ][idx_cvar_dd] );
            B_T1_mean = (real)0.5 * ( g_PriVar[NCOMP_TOTAL+TDir1][idx_cvar] + g_PriVar[NCOMP_TOTAL+TDir1][idx_cvar_dd] );
            B_T2_mean = (real)0.5 * ( g_PriVar[NCOMP_TOTAL+TDir2][idx_cvar] + g_PriVar[NCOMP_TOTAL+TDir2][idx_cvar_dd] );
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
//       normalize magnetic field
         B_N_mean  /= B_amp;
         B_T1_mean /= B_amp;
         B_T2_mean /= B_amp;
#        endif // #ifdef MHD

//       2. compute jacobian matrix for velocity gradients

         real N_slope_N, T1_slope_N, T2_slope_N;
         real N_slope_T1, T1_slope_T1, T2_slope_T1;
         real N_slope_T2, T1_slope_T2, T2_slope_T2;
         real divV, BBdV;

//       face-centered velocities
         real v_N  = get_vel_op( g_ConVar, g_PriVar,     d, idx_cvar, idx_cvar_dd, 0 );
         real v_T1 = get_vel_op( g_ConVar, g_PriVar, TDir1, idx_cvar, idx_cvar_dd, 0 );
         real v_T2 = get_vel_op( g_ConVar, g_PriVar, TDir2, idx_cvar, idx_cvar_dd, 0 );

//       normal direction derivative
         N_slope_N  = get_vel_op( g_ConVar, g_PriVar,     d, idx_cvar, idx_cvar_dd, 1 ) * _dh;
         T1_slope_N = get_vel_op( g_ConVar, g_PriVar, TDir1, idx_cvar, idx_cvar_dd, 1 ) * _dh;
         T2_slope_N = get_vel_op( g_ConVar, g_PriVar, TDir2, idx_cvar, idx_cvar_dd, 1 ) * _dh;

         if ( MicroPhy->ViscFluxType == ANISOTROPIC_VISCOSITY )
         {

//          transverse direction 1 derivative
            N_slope_T1   = van_leer2(
               get_vel_op( g_ConVar, g_PriVar,     d, idx_cvar_dd,                    idx_cvar_dd + didx_cvar[TDir1], 1 ),
               get_vel_op( g_ConVar, g_PriVar,     d, idx_cvar_dd - didx_cvar[TDir1], idx_cvar_dd,                    1 ),
               get_vel_op( g_ConVar, g_PriVar,     d, idx_cvar,                       idx_cvar    + didx_cvar[TDir1], 1 ),
               get_vel_op( g_ConVar, g_PriVar,     d, idx_cvar    - didx_cvar[TDir1], idx_cvar,                       1 )
            ) * _dh;
            T1_slope_T1  = van_leer2(
               get_vel_op( g_ConVar, g_PriVar, TDir1, idx_cvar_dd,                    idx_cvar_dd + didx_cvar[TDir1], 1 ),
               get_vel_op( g_ConVar, g_PriVar, TDir1, idx_cvar_dd - didx_cvar[TDir1], idx_cvar_dd,                    1 ),
               get_vel_op( g_ConVar, g_PriVar, TDir1, idx_cvar,                       idx_cvar    + didx_cvar[TDir1], 1 ),
               get_vel_op( g_ConVar, g_PriVar, TDir1, idx_cvar    - didx_cvar[TDir1], idx_cvar,                       1 )
            ) * _dh;
            T2_slope_T1  = van_leer2(
               get_vel_op( g_ConVar, g_PriVar, TDir2, idx_cvar_dd,                    idx_cvar_dd + didx_cvar[TDir1], 1 ),
               get_vel_op( g_ConVar, g_PriVar, TDir2, idx_cvar_dd - didx_cvar[TDir1], idx_cvar_dd,                    1 ),
               get_vel_op( g_ConVar, g_PriVar, TDir2, idx_cvar,                       idx_cvar    + didx_cvar[TDir1], 1 ),
               get_vel_op( g_ConVar, g_PriVar, TDir2, idx_cvar    - didx_cvar[TDir1], idx_cvar,                       1 )
            ) * _dh;

//          transverse direction 2 derivative
            N_slope_T2   = van_leer2(
               get_vel_op( g_ConVar, g_PriVar,     d, idx_cvar_dd,                    idx_cvar_dd + didx_cvar[TDir2], 1 ),
               get_vel_op( g_ConVar, g_PriVar,     d, idx_cvar_dd - didx_cvar[TDir1], idx_cvar_dd,                    1 ),
               get_vel_op( g_ConVar, g_PriVar,     d, idx_cvar,                       idx_cvar    + didx_cvar[TDir2], 1 ),
               get_vel_op( g_ConVar, g_PriVar,     d, idx_cvar    - didx_cvar[TDir1], idx_cvar,                       1 )
            ) * _dh;
            T1_slope_T2  = van_leer2(
               get_vel_op( g_ConVar, g_PriVar, TDir1, idx_cvar_dd,                    idx_cvar_dd + didx_cvar[TDir2], 1 ),
               get_vel_op( g_ConVar, g_PriVar, TDir1, idx_cvar_dd - didx_cvar[TDir1], idx_cvar_dd,                    1 ),
               get_vel_op( g_ConVar, g_PriVar, TDir1, idx_cvar,                       idx_cvar    + didx_cvar[TDir2], 1 ),
               get_vel_op( g_ConVar, g_PriVar, TDir1, idx_cvar    - didx_cvar[TDir1], idx_cvar,                       1 )
            ) * _dh;
            T2_slope_T2  = van_leer2(
               get_vel_op( g_ConVar, g_PriVar, TDir2, idx_cvar_dd,                    idx_cvar_dd + didx_cvar[TDir2], 1 ),
               get_vel_op( g_ConVar, g_PriVar, TDir2, idx_cvar_dd - didx_cvar[TDir1], idx_cvar_dd,                    1 ),
               get_vel_op( g_ConVar, g_PriVar, TDir2, idx_cvar,                       idx_cvar    + didx_cvar[TDir2], 1 ),
               get_vel_op( g_ConVar, g_PriVar, TDir2, idx_cvar    - didx_cvar[TDir1], idx_cvar,                       1 )
            ) * _dh;

            BBdV  = B_N_mean  * ( B_N_mean *  N_slope_N + B_T1_mean *  N_slope_T1 + B_T2_mean *  N_slope_T2 );
            BBdV += B_T1_mean * ( B_N_mean * T1_slope_N + B_T1_mean * T1_slope_T1 + B_T2_mean * T1_slope_T2 );
            BBdV += B_T2_mean * ( B_N_mean * T2_slope_N + B_T1_mean * T2_slope_T1 + B_T2_mean * T2_slope_T2 );

         }
         else // ISOTROPIC_VISCOSITY
         {

//          transverse direction 1 and 2 derivatives
            N_slope_T1  = 0.5 * (
               get_vel_op( g_ConVar, g_PriVar,     d, idx_cvar + didx_cvar[TDir1], idx_cvar_dd + didx_cvar[TDir1], 0 ) -
               get_vel_op( g_ConVar, g_PriVar,     d, idx_cvar + didx_cvar[TDir1], idx_cvar_dd - didx_cvar[TDir1], 0 )
            ) * _dh;
            N_slope_T2  = 0.5 * (
               get_vel_op( g_ConVar, g_PriVar,     d, idx_cvar + didx_cvar[TDir2], idx_cvar_dd + didx_cvar[TDir2], 0 ) -
               get_vel_op( g_ConVar, g_PriVar,     d, idx_cvar + didx_cvar[TDir2], idx_cvar_dd - didx_cvar[TDir2], 0 )
            ) * _dh;
            T1_slope_T1 = 0.5 * (
               get_vel_op( g_ConVar, g_PriVar, TDir1, idx_cvar + didx_cvar[TDir1], idx_cvar_dd + didx_cvar[TDir1], 0 ) -
               get_vel_op( g_ConVar, g_PriVar, TDir1, idx_cvar + didx_cvar[TDir1], idx_cvar_dd - didx_cvar[TDir1], 0 )
            ) * _dh;
            T2_slope_T2 = 0.5 * (
               get_vel_op( g_ConVar, g_PriVar, TDir2, idx_cvar + didx_cvar[TDir2], idx_cvar_dd + didx_cvar[TDir2], 0 ) -
               get_vel_op( g_ConVar, g_PriVar, TDir2, idx_cvar + didx_cvar[TDir2], idx_cvar_dd - didx_cvar[TDir2], 0 )
            ) * _dh;

            if ( MicroPhy->ViscSaturation )
            {

               T1_slope_T2 = 0.5 * (
                  get_vel_op( g_ConVar, g_PriVar, TDir1, idx_cvar + didx_cvar[TDir2], idx_cvar_dd + didx_cvar[TDir2], 0 ) -
                  get_vel_op( g_ConVar, g_PriVar, TDir1, idx_cvar + didx_cvar[TDir2], idx_cvar_dd - didx_cvar[TDir2], 0 )
               ) * _dh;
               T2_slope_T1 = 0.5 * (
                  get_vel_op( g_ConVar, g_PriVar, TDir2, idx_cvar + didx_cvar[TDir1], idx_cvar_dd + didx_cvar[TDir1], 0 ) -
                  get_vel_op( g_ConVar, g_PriVar, TDir2, idx_cvar + didx_cvar[TDir1], idx_cvar_dd - didx_cvar[TDir1], 0 )
               ) * _dh;

            }


         } // if ( MicroPhy->ViscFluxType == ANISOTROPIC_VISCOSITY ) ... else ...

//       divergence
         divV = N_slope_N + T1_slope_T1 + T2_slope_T2;

//       3. compute viscous flux
//       get the viscosity
//       --> for non-constant viscosity coefficients, take the spatial average along the normal direction
//           to get the face-centered coefficients
         real mu_l, mu_r, mu, visc_nu;
         real dens_L, dens_R, temp_L, temp_R;
         if ( g_PriVar == NULL )
         {
            dens_L  = g_ConVar[DENS][ idx_cvar    ];
            dens_R  = g_ConVar[DENS][ idx_cvar_dd ];
         }
         else
         {
            dens_L = g_PriVar[DENS][ idx_cvar    ];
            dens_R = g_PriVar[DENS][ idx_cvar_dd ];
         }
         temp_L = compute_temperature( g_ConVar, g_PriVar, g_FC_B, idx_cvar,
                                       MinTemp, PassiveFloor, EoS );
         temp_R = compute_temperature( g_ConVar, g_PriVar, g_FC_B, idx_cvar_dd,
                                       MinTemp, PassiveFloor, EoS );

         Hydro_ComputeViscosity( mu_l, visc_nu, MicroPhy, dens_L, temp_L );
         Hydro_ComputeViscosity( mu_r, visc_nu, MicroPhy, dens_R, temp_R );
         mu = 0.5*( mu_l + mu_r );

         if ( MicroPhy->ViscSaturation )
         {

            real norm_pi;

            if ( MicroPhy->ViscFluxType == ANISOTROPIC_VISCOSITY )
            {

               norm_pi = two_thirds * SQR( (3.0*BBdV - divV ) );

            }
            else // ISOTROPIC_VISCOSITY
            {

               norm_pi   = SQR( 2.0*N_slope_N   - two_thirds*divV );
               norm_pi  += SQR( 2.0*T1_slope_T1 - two_thirds*divV );
               norm_pi  += SQR( 2.0*T2_slope_T2 - two_thirds*divV );
               norm_pi  += 2.0 * SQR( N_slope_T1 + T1_slope_N );
               norm_pi  += 2.0 * SQR( N_slope_T2 + T2_slope_N );
               norm_pi  += 2.0 * SQR( T1_slope_T2 + T2_slope_T1 );

            }

            norm_pi = SQRT( norm_pi );
			if ( norm_pi > 0.0 )
			{
               real Thalf = 0.5 * ( temp_L + temp_R );
               real Dhalf = 0.5 * ( dens_L + dens_R );
               real v_th  = SQRT( MicroPhy->ViscThermalSpeedConv*Thalf );
               real l_i   = MicroPhy->ViscMFPConst * Thalf * Thalf / Dhalf;
               real l_v   = v_th / norm_pi;
               mu /= 1.0 + 2.0 * l_i / l_v;
            }

         }

         real stress_N, stress_T1, stress_T2;

         if ( MicroPhy->ViscFluxType == ISOTROPIC_VISCOSITY )
         {

            stress_N  = -mu * ( 2.0*N_slope_N - two_thirds*divV );
            stress_T1 = -mu * ( T1_slope_N + N_slope_T1 );
            stress_T2 = -mu * ( T2_slope_N + N_slope_T2 );

         }
         else // ANISOTROPIC_VISCOSITY
         {

#           ifdef MHD
            real delta_p;

            delta_p = mu*(3.0*BBdV - divV);

            if ( MicroPhy->ViscFluxType == ANISOTROPIC_VISCOSITY && MicroPhy->ViscBounds )
               delta_p = FMIN( FMAX( delta_p, -B2 ), 0.5*B2 );
            stress_N  = -delta_p*(B_N_mean*B_N_mean - 1./3.);
            stress_T1 = -delta_p*B_T1_mean*B_N_mean;
            stress_T2 = -delta_p*B_T2_mean*B_N_mean;
#           endif
	    
         } // if ( MicroPhy->ViscFluxType == ISOTROPIC_VISCOSITY ) ... else ...

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
