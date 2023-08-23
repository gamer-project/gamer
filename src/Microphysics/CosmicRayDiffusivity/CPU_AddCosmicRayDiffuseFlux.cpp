#ifndef __CPU_COSMICRAYDIFFUSE_FLUXES__
#define __CPU_COSMICRAYDIFFUSE_FLUXES__

#include "CUFLU.h"

#if ( ( MODEL == HYDRO ) && defined COSMIC_RAY && defined MICROPHYSICS && defined CR_DIFFUSION )

//external functions
#ifdef __CUDACC__

# include "CUFLU_ComputeCosmicRayDiffusivity.cu"

#else // #ifdef __CUDACC__

void CR_ComputeDiffusivity( real &diff_cr_para, real &diff_cr_perp, const MicroPhy_t *Mic );

#endif // #ifdef __CUDACC__ ... else ...


// internal funciton
real MC_limiter(real a, real b);
real minmod(real a, real b);



//-----------------------------------------------------------------------------------------
// Function    : CR_DiffuseFlux_HalfStep
//
// Description : Compute the half-step cosmic ray diffusive flux.
//
// Note        : 1. Must enable MHD and MICROPHYSICS.
//
// Reference   : Yang, H.-Y.~K., Ruszkowski, M., Ricker, P.~M., et al. 2012, apj, 761, 185. doi:10.1088/0004-637X/761/2/185
//
// Parameter   : g_Con_Var      : Array storing the input conserved fluid variables
//               g_Flux_Half    : Array storing the fluxes, and to store the cosmic rays diffusive flux.
//               g_FC_B         : Array storing the face-centered input B field
//               g_CC_B         : Array storing the cell-centered input B field
//               dh             : Cell size
//               Mic            : Microphysics object
//-----------------------------------------------------------------------------------------
GPU_DEVICE
void CR_DiffuseFlux_HalfStep( const real g_ConVar[][ CUBE(FLU_NXT) ],
                                    real g_Flux_Half[][NCOMP_TOTAL_PLUS_MAG][ CUBE(N_FC_FLUX) ],
                              const real g_FC_B[][ SQR(FLU_NXT)*FLU_NXT_P1 ],
                              const real g_CC_B[][ CUBE(FLU_NXT) ],
                              const real dh, const MicroPhy_t *Mic )
{
   const int didx_cvar[3] = { 1, FLU_NXT, SQR(FLU_NXT) };
   const int flux_offset = 1;
   const real small_B    = 1.e-30;

   for ( int d=0; d<3; d++ )
   {
      const int TDir1          = (d+1)%3;    // transverse direction 1
      const int TDir2          = (d+2)%3;    // transverse direction 2
      const int stride_fc_B[3] = { 1, FLU_NXT, SQR(FLU_NXT) };
      int sizeB_i, sizeB_j;

      int i_cvar_s=0, j_cvar_s=0, k_cvar_s=0, size_i, size_j, size_k;
      int i_offset, j_offset, k_offset;

      switch ( d )
      {
         case 0 : size_i  = N_HF_FLUX-1;  size_j  = N_HF_FLUX-0-2*flux_offset;  size_k = N_HF_FLUX-0-2*flux_offset;
                  i_offset = 0; j_offset = 1; k_offset = 1;
                  sizeB_i = FLU_NXT_P1;   sizeB_j = FLU_NXT;
                  break;

         case 1 : size_i  = N_HF_FLUX-0-2*flux_offset;  size_j  = N_HF_FLUX-1;  size_k = N_HF_FLUX-0-2*flux_offset;
                  i_offset = 1; j_offset = 0; k_offset = 1;
                  sizeB_i = FLU_NXT;      sizeB_j = FLU_NXT_P1;
                  break;

         case 2 : size_i  = N_HF_FLUX-0-2*flux_offset;  size_j  = N_HF_FLUX-0-2*flux_offset;  size_k = N_HF_FLUX-1;
                  i_offset = 1; j_offset = 1; k_offset = 0;
                  sizeB_i = FLU_NXT;      sizeB_j = FLU_NXT;
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
         const int i_cvar   = i_flux + i_cvar_s;
         const int j_cvar   = j_flux + j_cvar_s;
         const int k_cvar   = k_flux + k_cvar_s;
         const int idx_cvar = IDX321( i_cvar, j_cvar, k_cvar, FLU_NXT, FLU_NXT );

//       face-centered magnetic field index
         const int idx_fc_B = IDX321( i_cvar, j_cvar, k_cvar, sizeB_i, sizeB_j ) + stride_fc_B[d];

//       get the cell size
         real dh_N, dh_T1, dh_T2;
         dh_N  = dh;
         dh_T1 = dh;
         dh_T2 = dh;


//       1. get the diffusivity
         real diff_cr_eff_para, diff_cr_eff_perp;
         CR_ComputeDiffusivity( diff_cr_eff_para, diff_cr_eff_perp, Mic );

//       2. compute the mean magnetic field
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
         real B_N_mean, B_T1_mean, B_T2_mean, B_tot;
         B_N_mean  =       g_FC_B[    d][ idx_fc_B                ];
         B_T1_mean = 0.5*( g_CC_B[TDir1][ idx_cvar                ] +
                           g_CC_B[TDir1][ idx_cvar + didx_cvar[d] ] );
         B_T2_mean = 0.5*( g_CC_B[TDir2][ idx_cvar                ] +
                           g_CC_B[TDir2][ idx_cvar + didx_cvar[d] ] );
         B_tot    =  SQRT( B_N_mean*B_N_mean + B_T1_mean*B_T1_mean + B_T2_mean*B_T2_mean );

         //if ( B_tot < small_B && diff_cr_eff_perp != diff_cr_eff_para ) {
            // Error
            // xFlux = 0;
         //}

//       normalize magnetic field
         B_N_mean  = B_N_mean  / B_tot;
         B_T1_mean = B_T1_mean / B_tot;
         B_T2_mean = B_T2_mean / B_tot;


//       3. compute slope
//       ---------------------
//       |         |         |
//       ----al--------ar-----
//       |         |         |
//       |      N_slope      |
//       |         |         |
//       ----bl--------br-----
//       |         |         |
//       ---------------------
         real N_slope, T1_slope, T2_slope;
         real al, bl, ar, br;

//       normal direction
         N_slope = ( g_ConVar[CRAY][ idx_cvar + didx_cvar[d] ] -  g_ConVar[CRAY][ idx_cvar ] ) / dh_N;

//       transverse direction 1
         al = g_ConVar[CRAY][ idx_cvar                                   ] -
              g_ConVar[CRAY][ idx_cvar                - didx_cvar[TDir1] ]   ;
         bl = g_ConVar[CRAY][ idx_cvar                + didx_cvar[TDir1] ] -
              g_ConVar[CRAY][ idx_cvar                                   ]   ;
         ar = g_ConVar[CRAY][ idx_cvar + didx_cvar[d]                    ] -
              g_ConVar[CRAY][ idx_cvar + didx_cvar[d] - didx_cvar[TDir1] ]   ;
         br = g_ConVar[CRAY][ idx_cvar + didx_cvar[d] + didx_cvar[TDir1] ] -
              g_ConVar[CRAY][ idx_cvar + didx_cvar[d]                    ]   ;
         T1_slope = ( MC_limiter( MC_limiter(al, bl), MC_limiter(ar, br) ) ) / dh_T1;

//       transverse direction 2
         al = g_ConVar[CRAY][ idx_cvar                                   ] -
              g_ConVar[CRAY][ idx_cvar                - didx_cvar[TDir2] ]   ;
         bl = g_ConVar[CRAY][ idx_cvar                + didx_cvar[TDir2] ] -
              g_ConVar[CRAY][ idx_cvar                                   ]   ;
         ar = g_ConVar[CRAY][ idx_cvar + didx_cvar[d]                    ] -
              g_ConVar[CRAY][ idx_cvar + didx_cvar[d] - didx_cvar[TDir2] ]   ;
         br = g_ConVar[CRAY][ idx_cvar + didx_cvar[d] + didx_cvar[TDir2] ] -
              g_ConVar[CRAY][ idx_cvar + didx_cvar[d]                    ]   ;
         T2_slope = ( MC_limiter( MC_limiter(al, bl), MC_limiter(ar, br) ) ) / dh_T2;



//       4. compute Flux
         real Flux_1Face, Flux_Para, Flux_Perp;
         real common = -B_N_mean * ( B_N_mean*N_slope + B_T1_mean*T1_slope + B_T2_mean*T2_slope );
         Flux_Para = diff_cr_eff_para * common;
         Flux_Perp = diff_cr_eff_perp * (-N_slope - common);
         Flux_1Face = Flux_Para + Flux_Perp;

//       5. flux add-up
         g_Flux_Half[d][CRAY][idx_flux] += Flux_1Face;
         g_Flux_Half[d][ENGY][idx_flux] += Flux_1Face;
         // TODO: dual energy (internal energy) update?

      } // CGPU_LOOP( idx, size_i*size_j*size_k )
   } // for (int d=0; d<3; d++)

#  ifdef __CUDACC__
     __syncthreads();
#  endif

} // FUNCTION : CR_DiffuseFlux_HalfStep




//-----------------------------------------------------------------------------------------
// Function    : CR_DiffuseFlux_FullStep
//
// Description : Compute the full-step cosmic ray diffusive flux.
//
// Note        : 1. Must enable MHD and MICROPHYSICS.
//
// Reference   : Yang, H.-Y.~K., Ruszkowski, M., Ricker, P.~M., et al. 2012, apj, 761, 185. doi:10.1088/0004-637X/761/2/185
//
// Parameter   : g_PriVar_Half  : Array storing the half-step primitive fluid variables.
//               g_FC_Flux      : Array storing the fluxes, and to store the cosmic rays diffusive flux.
//               g_FC_B_Half    : Array storing the half-step face-centered magnetic field
//               NFlux          : Stride for accessing g_FC_Flux
//               dh             : Cell size
//               Mic            : Microphysics object
//-----------------------------------------------------------------------------------------
GPU_DEVICE
void CR_DiffuseFlux_FullStep( const real g_PriVar_Half[][ CUBE(FLU_NXT) ],
                                    real g_FC_Flux[][NCOMP_TOTAL_PLUS_MAG][ CUBE(N_FC_FLUX) ],
                              const real g_FC_B_Half[][ FLU_NXT_P1*SQR(FLU_NXT) ],
                              const int NFlux, const real dh, const MicroPhy_t *Mic )
{

   const int didx_half[3] = { 1, N_HF_VAR, SQR(N_HF_VAR) };
   const int mag_offset  = (N_HF_VAR - PS2)/2;
   const int cell_offset = (N_HF_VAR - N_FC_VAR)/2;
   const real small_B    = 1.e-30;

   for ( int d=0; d<3; d++ )
   {
      const int TDir1          = (d+1)%3;    // transverse direction 1
      const int TDir2          = (d+2)%3;    // transverse direction 2
      int sizeB_i, sizeB_j, sizeB_k;
      int mag_offset_i, mag_offset_j, mag_offset_k;

      int idx_flux_e[3];

      switch ( d )
      {
         case 0 : idx_flux_e[0] = N_FC_VAR-1;   idx_flux_e[1] = N_FC_VAR;     idx_flux_e[2] = N_FC_VAR;
                  sizeB_i       = N_HF_VAR+1;   sizeB_j       = N_HF_VAR;     sizeB_k       = N_HF_VAR;
                  mag_offset_i  = mag_offset;   mag_offset_j  = mag_offset-1; mag_offset_k  = mag_offset-1;
                  break;

         case 1 : idx_flux_e[0] = N_FC_VAR;     idx_flux_e[1] = N_FC_VAR-1;   idx_flux_e[2] = N_FC_VAR;
                  sizeB_i       = N_HF_VAR;     sizeB_j       = N_HF_VAR+1;   sizeB_k       = N_HF_VAR;
                  mag_offset_i  = mag_offset-1; mag_offset_j  = mag_offset;   mag_offset_k  = mag_offset-1;
                  break;

         case 2 : idx_flux_e[0] = N_FC_VAR;     idx_flux_e[1] = N_FC_VAR;     idx_flux_e[2] = N_FC_VAR-1;
                  sizeB_i       = N_HF_VAR;     sizeB_j       = N_HF_VAR;     sizeB_k       = N_HF_VAR+1;
                  mag_offset_i  = mag_offset-1; mag_offset_j  = mag_offset-1; mag_offset_k  = mag_offset;
                  break;
      } // switch ( d )

      const int size_ij = idx_flux_e[0]*idx_flux_e[1];
      CGPU_LOOP( idx, idx_flux_e[0]*idx_flux_e[1]*idx_flux_e[2] )
      {
//       flux index
         const int i_flux   = idx % idx_flux_e[0];
         const int j_flux   = idx % size_ij / idx_flux_e[0];
         const int k_flux   = idx / size_ij;
         const int idx_flux = IDX321( i_flux, j_flux, k_flux, NFlux, NFlux );

//       half variable index
         const int i_half   = i_flux + cell_offset;
         const int j_half   = j_flux + cell_offset;
         const int k_half   = k_flux + cell_offset;
         const int idx_half = IDX321( i_half, j_half, k_half, N_HF_VAR, N_HF_VAR );

//       magnetic field indexes
         const int i_fc     = i_flux + mag_offset_i;
         const int j_fc     = j_flux + mag_offset_j;
         const int k_fc     = k_flux + mag_offset_k;
         const int idx_fc_BN  = IDX321( i_fc, j_fc, k_fc, sizeB_i, sizeB_j );
         const int idx_fc_BT1 = IDX321( i_fc, j_fc, k_fc, sizeB_k, sizeB_i );
         const int idx_fc_BT2 = IDX321( i_fc, j_fc, k_fc, sizeB_j, sizeB_k );
         const int stride_fc_BN[3]  = { 1, sizeB_i, sizeB_i*sizeB_j };
         const int stride_fc_BT1[3] = { 1, sizeB_k, sizeB_k*sizeB_i };
         const int stride_fc_BT2[3] = { 1, sizeB_j, sizeB_j*sizeB_k };

//       get the cell size
         real dh_N, dh_T1, dh_T2;
         dh_N  = dh;
         dh_T1 = dh;
         dh_T2 = dh;


//       1. get the diffusivity
         real diff_cr_eff_para, diff_cr_eff_perp;
         CR_ComputeDiffusivity( diff_cr_eff_para, diff_cr_eff_perp, Mic );

//       2. compute the mean magnetic field
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
         real B_N_mean, B_T1_mean, B_T2_mean, B_tot;
         B_N_mean =         g_FC_B_Half[    d][ idx_fc_BN                                            ];
         B_T1_mean = 0.25*( g_FC_B_Half[TDir1][ idx_fc_BT1                                           ] +
                            g_FC_B_Half[TDir1][ idx_fc_BT1                    + stride_fc_BT1[TDir1] ] +
                            g_FC_B_Half[TDir1][ idx_fc_BT1 - stride_fc_BT1[d]                        ] +
                            g_FC_B_Half[TDir1][ idx_fc_BT1 - stride_fc_BT1[d] + stride_fc_BT1[TDir1] ] );
         B_T2_mean = 0.25*( g_FC_B_Half[TDir2][ idx_fc_BT2                                           ] +
                            g_FC_B_Half[TDir2][ idx_fc_BT2                    + stride_fc_BT2[TDir2] ] +
                            g_FC_B_Half[TDir2][ idx_fc_BT2 - stride_fc_BT2[d]                        ] +
                            g_FC_B_Half[TDir2][ idx_fc_BT2 - stride_fc_BT2[d] + stride_fc_BT2[TDir2] ] );
         B_tot     = SQRT( B_N_mean*B_N_mean + B_T1_mean*B_T1_mean + B_T2_mean*B_T2_mean );

         //if ( B_tot < small_B && diff_cr_eff_perp != diff_cr_eff_para ) {
            // Error
            // xFlux = 0;
         //}

//       normalize magnetic field
         B_N_mean  = B_N_mean  / B_tot;
         B_T1_mean = B_T1_mean / B_tot;
         B_T2_mean = B_T2_mean / B_tot;


//       3. compute cosmic ray slope
//       ---------------------
//       |         |         |
//       ----al--------ar-----
//       |         |         |
//       |      N_slope      |
//       |         |         |
//       ----bl--------br-----
//       |         |         |
//       ---------------------
         real N_slope, T1_slope, T2_slope;
         real al, bl, ar, br;

//       normal direction
         N_slope = ( g_PriVar_Half[CRAY][ idx_half + didx_half[d] ] -  g_PriVar_Half[CRAY][ idx_half ] ) / dh_N;

//       transverse direction 1
         al = g_PriVar_Half[CRAY][ idx_half                                   ] -
              g_PriVar_Half[CRAY][ idx_half                - didx_half[TDir1] ]   ;
         bl = g_PriVar_Half[CRAY][ idx_half                + didx_half[TDir1] ] -
              g_PriVar_Half[CRAY][ idx_half                                   ]   ;
         ar = g_PriVar_Half[CRAY][ idx_half + didx_half[d]                    ] -
              g_PriVar_Half[CRAY][ idx_half + didx_half[d] - didx_half[TDir1] ]   ;
         br = g_PriVar_Half[CRAY][ idx_half + didx_half[d] + didx_half[TDir1] ] -
              g_PriVar_Half[CRAY][ idx_half + didx_half[d]                    ]   ;
         T1_slope = ( MC_limiter( MC_limiter(al, bl), MC_limiter(ar, br) ) ) / dh_T1;

//       transverse direction 2
         al = g_PriVar_Half[CRAY][ idx_half                                   ] -
              g_PriVar_Half[CRAY][ idx_half                - didx_half[TDir2] ]   ;
         bl = g_PriVar_Half[CRAY][ idx_half                + didx_half[TDir2] ] -
              g_PriVar_Half[CRAY][ idx_half                                   ]   ;
         ar = g_PriVar_Half[CRAY][ idx_half + didx_half[d]                    ] -
              g_PriVar_Half[CRAY][ idx_half + didx_half[d] - didx_half[TDir2] ]   ;
         br = g_PriVar_Half[CRAY][ idx_half + didx_half[d] + didx_half[TDir2] ] -
              g_PriVar_Half[CRAY][ idx_half + didx_half[d]                    ]   ;
         T2_slope = ( MC_limiter( MC_limiter(al, bl), MC_limiter(ar, br) ) ) / dh_T2;



//       4. compute Flux
         real Flux_1Face, Flux_Para, Flux_Perp;
         real common = -B_N_mean * ( B_N_mean*N_slope + B_T1_mean*T1_slope + B_T2_mean*T2_slope );
         Flux_Para = diff_cr_eff_para * common;
         Flux_Perp = diff_cr_eff_perp * (-N_slope - common);
         Flux_1Face = Flux_Para + Flux_Perp;

//       5. Flux add-up
         g_FC_Flux[d][CRAY][idx_flux] += Flux_1Face;
         g_FC_Flux[d][ENGY][idx_flux] += Flux_1Face;
         // TODO: dual energy (internal energy) update?

      } // CGPU_LOOP( idx, idx_flux_e[0]*idx_flux_e[1]*idx_flux_e[2] )
   } // for (int d=0; d<3; d++)

#  ifdef __CUDACC__
     __syncthreads();
#  endif

} // FUNCTION : CR_DiffuseFlux_FullStep



//-----------------------------------------------------------------------------------------
// Function    : MC_limiter
// Description :
// Note        :
// Parameter   : a :
//               b :
// Return      :
// Reference   :
//-----------------------------------------------------------------------------------------
real MC_limiter( real a, real b )
{
    return minmod( 2.0*minmod(a, b), 0.5*(a+b) );
} // FUNCTION : MC_limiter



//-----------------------------------------------------------------------------------------
// Function    : minmod
// Description :
// Note        :
// Parameter   : a :
//               b :
// Return      :
// Reference   :
//-----------------------------------------------------------------------------------------
real minmod( real a, real b )
{
    if ( a > 0.0 && b > 0.0 ) return MIN(a, b);
    if ( a < 0.0 && b < 0.0 ) return MAX(a, b);
    if ( a * b <= 0.0 )       return 0.0;
} // FUNCTION : minmod

#endif // #if ( ( MODEL == HYDRO ) && defined COSMIC_RAY && defined MICROPHYSICS && defined CR_DIFFUSION )

#endif // #ifndef __CPU_COSMICRAYDIFFUSE_FLUXES__
