#ifndef __CUFLU_HYDRO_ADDEXTRAFLUX__
#define __CUFLU_HYDRO_ADDEXTRAFLUX__



#include "GAMER.h"
#include "CUFLU.h"



#if (  MODEL == HYDRO  &&  ( FLU_SCHEME == MHM || FLU_SCHEME == MHM_RP )  )



//-----------------------------------------------------------------------------------------
// Function    : AddExtraFlux_Template
//
// Description : Add extra flux
//
// Note        : 1. Only for MHM_RP half-step, MHM_RP full-step, and MHM full-step
//               2. Does not upadate the additional flux along transverse direction for computing the CT electric field
//
// Usage       :     MHM_RP half-step: AddExtraFlux_Template( g_Flu_Array_In[P],              NULL, g_Flux_Half_1PG, g_Mag_Array_In[P],  FLU_NXT,             0, N_HF_FLUX, NSkip_N, NSkip_T, dh );
//               MHM/MHM_RP full-step: AddExtraFlux_Template(              NULL, g_PriVar_Half_1PG,   g_FC_Flux_1PG, g_FC_Mag_Half_1PG, N_HF_VAR, LR_GHOST_SIZE, N_FL_FLUX, NSkip_N, NSkip_T, dh );
//
// Parameter   : g_ConVar : Array storing the input cell-centered conserved fluid variables (NCOMP_TOTAL)
//               g_PriVar : Array storing the input cell-centered primitive fluid variables (NCOMP_TOTAL_PLUS_MAG)
//               g_Flux   : Array with hydrodynamic fluxes for adding the extra fluxes
//               g_FC_B   : Array storing the input face-centered B field
//               N_Var    : Size of g_ConVar/g_PriVar
//               N_Ghost  : Ghost zone size of data-reconstruction
//               N_Flux   : Size of flux
//               NSkip_N  : Empty size of flux on normal     direction
//               NSkip_T  : Empty size of flux on transverse direction
//               dh       : Cell size
//
// Return      : g_Flux[]
//-----------------------------------------------------------------------------------------
void AddExtraFlux_Template( const real g_ConVar[][ CUBE(FLU_NXT) ],
                            const real g_PriVar[][ CUBE(FLU_NXT) ],
                                  real g_Flux[][NCOMP_TOTAL_PLUS_MAG][ CUBE(N_FC_FLUX) ],
                            const real g_FC_B[][ SQR(FLU_NXT)*FLU_NXT_P1 ],
                            const int N_Var,
                            const int N_Ghost,
                            const int N_Flux,
                            const int NSkip_N,
                            const int NSkip_T,
                            const real dh )
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
         case 0 : flux_size_i   = N_Var-2*N_Ghost-NSkip_N-1;                 flux_size_j   = N_Var-2*N_Ghost-2*NSkip_T-2*skip_update_T; flux_size_k   = N_Var-2*N_Ghost-2*NSkip_T-2*skip_update_T;
                  flux_offset_i = 0;                                         flux_offset_j = skip_update_T;                             flux_offset_k = skip_update_T;
                  var_offset_i  = N_Ghost+NSkip_N;                           var_offset_j  = N_Ghost+NSkip_T;                           var_offset_k  = N_Ghost+NSkip_T;
                  // not sure if var_offset_i works for CTU
                  break;

         case 1 : flux_size_i   = N_Var-2*N_Ghost-2*NSkip_T-2*skip_update_T; flux_size_j   = N_Var-2*N_Ghost-NSkip_N-1;                 flux_size_k   = N_Var-2*N_Ghost-2*NSkip_T-2*skip_update_T;
                  flux_offset_i = skip_update_T;                             flux_offset_j = 0;                                         flux_offset_k = skip_update_T;
                  var_offset_i  = N_Ghost+NSkip_T;                           var_offset_j  = N_Ghost+NSkip_N;                           var_offset_k  = N_Ghost+NSkip_T;
                  break;

         case 2 : flux_size_i   = N_Var-2*N_Ghost-2*NSkip_T-2*skip_update_T; flux_size_j   = N_Var-2*N_Ghost-2*NSkip_T-2*skip_update_T; flux_size_k   = N_Var-2*N_Ghost-NSkip_N-1;
                  flux_offset_i = skip_update_T;                             flux_offset_j = skip_update_T;                             flux_offset_k = 0;
                  var_offset_i  = N_Ghost+NSkip_T;                           var_offset_j  = N_Ghost+NSkip_T;                           var_offset_k  = N_Ghost+NSkip_N;
                  break;
      } // switch ( d )

#     ifdef MHD
      int sizeB_i, sizeB_j, sizeB_k;
      int mag_offset_i, mag_offset_j, mag_offset_k;
      switch ( d )
      {
         case 0 : sizeB_i      = N_Var+1;   sizeB_j      = N_Var;     sizeB_k      = N_Var;
                  mag_offset_i = N_Ghost+1; mag_offset_j = N_Ghost;   mag_offset_k = N_Ghost;
                  break;
         case 1 : sizeB_i      = N_Var;     sizeB_j      = N_Var+1;   sizeB_k      = N_Var;
                  mag_offset_i = N_Ghost;   mag_offset_j = N_Ghost+1; mag_offset_k = N_Ghost;
                  break;
         case 2 : sizeB_i      = N_Var;     sizeB_j      = N_Var;     sizeB_k      = N_Var+1;
                  mag_offset_i = N_Ghost;   mag_offset_j = N_Ghost;   mag_offset_k = N_Ghost+1;
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

/*
#ifdef MHD
         printf("d=%d, flux: (%02d, %02d, %02d, %05d), cvar: (%02d, %02d, %02d, %05d), fc: (%02d, %02d, %02d) => (%04d, %04d, %04d), N_Flux=%02d, N_Var=%02d, sizeB=(%02d %02d %02d)\n",
                 d,
                 i_flux, j_flux, k_flux, idx_flux,
                 i_cvar, j_cvar, k_cvar, idx_cvar,
                 i_fc, j_fc, k_fc,
                 idx_fc_BN, idx_fc_BT1, idx_fc_BT2,
                 N_Flux, N_Var, sizeB_i, sizeB_j, sizeB_k);
#else
         printf("d=%d, flux: (%02d, %02d, %02d, %05d), cvar: (%02d, %02d, %02d, %05d), N_Flux=%02d, N_Var=%02d\n",
                 d,
                 i_flux, j_flux, k_flux, idx_flux,
                 i_cvar, j_cvar, k_cvar, idx_cvar,
                 N_Flux, N_Var);
#endif
*/

//       2. compute the mean magnetic field at the face-centered flux location
#        ifdef MHD
         real B_N_mean, B_T1_mean, B_T2_mean;
         if ( g_FC_B == NULL )
         {
            B_N_mean  = g_PriVar[NCOMP_TOTAL+d    ][idx_cvar];
            B_T1_mean = g_PriVar[NCOMP_TOTAL+TDir1][idx_cvar];
            B_T2_mean = g_PriVar[NCOMP_TOTAL+TDir2][idx_cvar];
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
#        endif

//       3. compute the mean density slope at the face-centered flux location
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
         N_slope = ( g_ConVar[DENS][ idx_cvar + didx_cvar[d] ] - g_ConVar[DENS][ idx_cvar ] ) * _dh;

//       transverse direction 1
         al = g_ConVar[DENS][ idx_cvar                                   ] -
              g_ConVar[DENS][ idx_cvar                - didx_cvar[TDir1] ];
         bl = g_ConVar[DENS][ idx_cvar                + didx_cvar[TDir1] ] -
              g_ConVar[DENS][ idx_cvar                                   ];
         ar = g_ConVar[DENS][ idx_cvar + didx_cvar[d]                    ] -
              g_ConVar[DENS][ idx_cvar + didx_cvar[d] - didx_cvar[TDir1] ];
         br = g_ConVar[DENS][ idx_cvar + didx_cvar[d] + didx_cvar[TDir1] ] -
              g_ConVar[DENS][ idx_cvar + didx_cvar[d]                    ];
         T1_slope = 0.25 * (al + bl + ar + br) * _dh;

//       transverse direction 2
         al = g_ConVar[DENS][ idx_cvar                                   ] -
              g_ConVar[DENS][ idx_cvar                - didx_cvar[TDir2] ];
         bl = g_ConVar[DENS][ idx_cvar                + didx_cvar[TDir2] ] -
              g_ConVar[DENS][ idx_cvar                                   ];
         ar = g_ConVar[DENS][ idx_cvar + didx_cvar[d]                    ] -
              g_ConVar[DENS][ idx_cvar + didx_cvar[d] - didx_cvar[TDir2] ];
         br = g_ConVar[DENS][ idx_cvar + didx_cvar[d] + didx_cvar[TDir2] ] -
              g_ConVar[DENS][ idx_cvar + didx_cvar[d]                    ];
         T2_slope = 0.25 * (al + bl + ar + br) * _dh;

//       4. do your calculation of the extra flux here
         real Extra_Flux = (real)0.0;

//       5. flux add-up
         g_Flux[d][DENS][idx_flux] += Extra_Flux;
      } // CGPU_LOOP( idx, size_i*size_j*size_k )
   } // for (int d=0; d<3; d++)

#  ifdef __CUDACC__
   __syncthreads();
#  endif

} // FUNCTION : AddExtraFlux_Template



//-----------------------------------------------------------------------------------------
// Function    : AddExtraFlux_AddExtraFlux_HancockPredict_TemplateTemplate
//
// Description : Add extra flux only in Hancock predict
//
// Note        : 1. Only for MHM half-step
//
// Usage       : MHM half-step : AddExtraFlux_HancockPredict_Template( Flux, g_cc_array, g_FC_B, cc_idx, cc_i, cc_j, cc_k, NGhost, dh );
//
// Parameter   : Flux         : Array with hydrodynamic fluxes for adding the extra fluxes
//               g_g_cc_array : Array storing the input cell-centered conserved fluid variables
//               g_FC_B       : Array storing the input face-centered B field
//               cc_idx/i/j/k : Index of g_cc_array
//               NGhost       : Ghost zone size of data-reconstruction
//               dh           : Cell size
//
// Return      : Flux[]
//-----------------------------------------------------------------------------------------
GPU_DEVICE
void AddExtraFlux_HancockPredict_Template(       real Flux[][NCOMP_TOTAL_PLUS_MAG],
                                           const real g_cc_array[][ CUBE(FLU_NXT) ],
                                           const real g_FC_B[][ FLU_NXT_P1*SQR(FLU_NXT) ],
                                           const int cc_idx, const int cc_i, const int cc_j, const int cc_k,
                                           const int NGhost, const real dh )
{
   const int didx_cvar[3] = { 1, FLU_NXT, SQR(FLU_NXT) };
   const real _dh          = (real)1.0 / dh;

   for (int f=0; f<6; f++)
   {
//    -------------
//    |     ^     |
//    ------3------
//    |     |     |
//  --0-> i j k --1->
//    |     ^     |
//    ------2------
//    |     |     |
//    -------------
      const int d     = f / 2;
      const int TDir1 = (d+1)%3;    // transverse direction 1
      const int TDir2 = (d+2)%3;    // transverse direction 2
      const int LorR  = f % 2; // 0: Left, 1: Right

#     ifdef MHD
      int sizeB_i, sizeB_j, sizeB_k;
      switch ( d )
      {
         case 0 : sizeB_i = FLU_NXT_P1; sizeB_j = FLU_NXT;    sizeB_k = FLU_NXT;    break;
         case 1 : sizeB_i = FLU_NXT;    sizeB_j = FLU_NXT_P1; sizeB_k = FLU_NXT;    break;
         case 2 : sizeB_i = FLU_NXT;    sizeB_j = FLU_NXT;    sizeB_k = FLU_NXT_P1; break;
      } // switch ( d )
      const int stride_fc_BT1[3] = { 1, sizeB_k, sizeB_k*sizeB_i };
      const int stride_fc_BT2[3] = { 1, sizeB_j, sizeB_j*sizeB_k };

      const int idx_fc_BN  = IDX321( cc_i+LorR, cc_j, cc_j, sizeB_i, sizeB_j );
      const int idx_fc_BT1 = IDX321( cc_i+LorR, cc_j, cc_j, sizeB_k, sizeB_i );
      const int idx_fc_BT2 = IDX321( cc_i+LorR, cc_j, cc_j, sizeB_j, sizeB_k );
#     endif

/*
#ifdef MHD
      printf("f: %d, d: %d, LR: %d, cc: (%02d, %02d, %02d, %05d), fc: (%02d, %02d, %02d) => (%05d, %05d, %05d)\n",
              f, d, LorR,
              cc_i, cc_j, cc_k, cc_idx,
              cc_i+LorR, cc_j, cc_k,
              idx_fc_BN, idx_fc_BT1, idx_fc_BT2);
#else
      printf("f: %d, d: %d, LR: %d, cc: (%02d, %02d, %02d, %05d)\n",
              f, d, LorR,
              cc_i, cc_j, cc_k, cc_idx
              );
#endif
*/

//    1. compute the mean magnetic field at the face-centered flux location
#     ifdef MHD
      real B_N_mean, B_T1_mean, B_T2_mean;
      B_N_mean  =              g_FC_B[    d][ idx_fc_BN                                            ];
      B_T1_mean = (real)0.25*( g_FC_B[TDir1][ idx_fc_BT1                                           ] +
                               g_FC_B[TDir1][ idx_fc_BT1                    + stride_fc_BT1[TDir1] ] +
                               g_FC_B[TDir1][ idx_fc_BT1 - stride_fc_BT1[d]                        ] +
                               g_FC_B[TDir1][ idx_fc_BT1 - stride_fc_BT1[d] + stride_fc_BT1[TDir1] ]   );
      B_T2_mean = (real)0.25*( g_FC_B[TDir2][ idx_fc_BT2                                           ] +
                               g_FC_B[TDir2][ idx_fc_BT2                    + stride_fc_BT2[TDir2] ] +
                               g_FC_B[TDir2][ idx_fc_BT2 - stride_fc_BT2[d]                        ] +
                               g_FC_B[TDir2][ idx_fc_BT2 - stride_fc_BT2[d] + stride_fc_BT2[TDir2] ]   );
#     endif

//    2. calculate the mean density slope at flux location
//    ---------------------
//    |         |         |
//    ----bl--------br-----
//    |         |         |
//    |      N_slope      |
//    |         |         |
//    ----al--------ar-----
//    |         |         |
//    ---------------------
      real N_slope, T1_slope, T2_slope;
      real al, bl, ar, br;

//    normal direction
      N_slope = ( g_cc_array[DENS][ cc_idx + LorR*didx_cvar[d] ] - g_cc_array[DENS][ cc_idx + LorR*didx_cvar[d] - didx_cvar[d] ] ) * _dh;

//    transverse direction 1
      al = g_cc_array[DENS][ cc_idx + LorR*didx_cvar[d] - didx_cvar[d]                    ] -
           g_cc_array[DENS][ cc_idx + LorR*didx_cvar[d] - didx_cvar[d] - didx_cvar[TDir1] ];
      bl = g_cc_array[DENS][ cc_idx + LorR*didx_cvar[d] - didx_cvar[d] + didx_cvar[TDir1] ] -
           g_cc_array[DENS][ cc_idx + LorR*didx_cvar[d] - didx_cvar[d]                    ];
      ar = g_cc_array[DENS][ cc_idx + LorR*didx_cvar[d]                                   ] -
           g_cc_array[DENS][ cc_idx + LorR*didx_cvar[d]                - didx_cvar[TDir1] ];
      br = g_cc_array[DENS][ cc_idx + LorR*didx_cvar[d]                + didx_cvar[TDir1] ] -
           g_cc_array[DENS][ cc_idx + LorR*didx_cvar[d]                                   ];
      T1_slope = 0.25 * (al + bl + ar + br) * _dh;

//    transverse direction 2
      al = g_cc_array[DENS][ cc_idx + LorR*didx_cvar[d] - didx_cvar[d]                    ] -
           g_cc_array[DENS][ cc_idx + LorR*didx_cvar[d] - didx_cvar[d] - didx_cvar[TDir2] ];
      bl = g_cc_array[DENS][ cc_idx + LorR*didx_cvar[d] - didx_cvar[d] + didx_cvar[TDir2] ] -
           g_cc_array[DENS][ cc_idx + LorR*didx_cvar[d] - didx_cvar[d]                    ];
      ar = g_cc_array[DENS][ cc_idx + LorR*didx_cvar[d]                                   ] -
           g_cc_array[DENS][ cc_idx + LorR*didx_cvar[d]                - didx_cvar[TDir2] ];
      br = g_cc_array[DENS][ cc_idx + LorR*didx_cvar[d]                + didx_cvar[TDir2] ] -
           g_cc_array[DENS][ cc_idx + LorR*didx_cvar[d]                                   ];
      T2_slope = 0.25 * (al + bl + ar + br) * _dh;

//    3. do your calculation of the extra flux here
      real Extra_Flux = (real)0.0;

//    4. flux add-up
      Flux[f][DENS] += Extra_Flux;
   } // for (int f=0; f<6; f++)
} // FUNCTION : AddExtraFlux_HancockPredict_Template



#endif // #if (  MODEL == HYDRO  &&  ( FLU_SCHEME == MHM || FLU_SCHEME == MHM_RP )  )



#endif // #ifdef __CUFLU_HYDRO_ADDEXTRAFLUX__
