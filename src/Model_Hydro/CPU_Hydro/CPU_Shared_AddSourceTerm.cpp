#ifndef __CUFLU_HYDRO_ADDSOURCETERM__
#define __CUFLU_HYDRO_ADDSOURCETERM__




#include "CUFLU.h"



#if (  MODEL == HYDRO  &&  ( FLU_SCHEME == MHM  ||  FLU_SCHEME == MHM_RP )  )



#if ( FLU_SCHEME == MHM )
//-------------------------------------------------------------------------------------------------------
// Function    :  Hydro_AddSourceTerm_HalfStep_MHM
//
// Description :  Add source term on the face-centered conserved varibles at half-step update
//
// Note        :  1. This function can only be used by MHM
//                2. Invoked by Hydro_HancockPredict()
//                3. Must be used after MHD_UpdateMagnetic_Half()
//
// Parameter   :  g_ConVar_In : Array storing the input conserved variables
//                g_FC_B_In   : Array storing the input face-centered magnetic field (for MHD only)
//                fcCon       : Face-centered conserved variables to be updated
//                idx_in      : Index of accessing g_ConVar_In[]
//                dt_dh2      : 0.5 * dt / dh
//                EoS         : EoS object
//
// Return      :  fcCon
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE
void Hydro_AddSourceTerm_HalfStep_MHM( const real g_ConVar_In[][ CUBE(FLU_NXT) ],
                                       const real g_FC_B_In[][ FLU_NXT_P1*SQR(FLU_NXT) ],
                                             real fcCon[][NCOMP_LR],
                                       const int idx_in, const real dt_dh2, const EoS_t *EoS )
{

   const int didx_in[3] = { 1, FLU_NXT, SQR(FLU_NXT) };
   const int i_in       = idx_in % FLU_NXT;
   const int j_in       = idx_in % SQR(FLU_NXT) / FLU_NXT;
   const int k_in       = idx_in / SQR(FLU_NXT);

// |    |            |    |
// ----------f=3-----------
// |    |            |    |
// |   f=0  idx_in  f=1   |
// |    |            |    |
// ----------f=2-----------
// |    |            |    |
   for (int f=0; f<6; f++)
   {
      const int d     = f/2;
      const int TDir1 = (d+1)%3;    // transverse direction 1
      const int TDir2 = (d+2)%3;    // transverse direction 2
      const int LR    = 2*(f%2) - 1;

#     ifdef MHD
      int size_B_i, size_B_j, size_B_k;
      switch ( d )
      {
         case 0 : size_B_i = FLU_NXT_P1; size_B_j = FLU_NXT;    size_B_k = FLU_NXT;    break;
         case 1 : size_B_i = FLU_NXT;    size_B_j = FLU_NXT_P1; size_B_k = FLU_NXT;    break;
         case 2 : size_B_i = FLU_NXT;    size_B_j = FLU_NXT;    size_B_k = FLU_NXT_P1; break;
      } // switch ( d )

//    index of magnetic field
      const int idx_B_N  = IDX321( i_in, j_in, k_in, size_B_i, size_B_j );
      const int idx_B_T1 = IDX321( i_in, j_in, k_in, size_B_k, size_B_i );
      const int idx_B_T2 = IDX321( i_in, j_in, k_in, size_B_j, size_B_k );

//    index increment of magnetic field
      const int didx_B_N [3] = { 1, size_B_i, size_B_i*size_B_j };
      const int didx_B_T1[3] = { 1, size_B_k, size_B_k*size_B_i };
      const int didx_B_T2[3] = { 1, size_B_j, size_B_j*size_B_k };
#     endif

//    1. calculate extra term
      real ExtraTerm = 0.0;
#     ifdef MHD
//    magnetic field at face center
//    const real B_N  =                g_FC_B_In[d    ][ idx_B_N  + (f%2)*didx_B_N[d]                  ];
//    const real B_T1 = (real)0.25 * ( g_FC_B_In[TDir1][ idx_B_T1                                      ] +
//                                     g_FC_B_In[TDir1][ idx_B_T1 + LR*didx_B_T1[d]                    ] +
//                                     g_FC_B_In[TDir1][ idx_B_T1                   + didx_B_T1[TDir1] ] +
//                                     g_FC_B_In[TDir1][ idx_B_T1 + LR*didx_B_T1[d] + didx_B_T1[TDir1] ] );
//    const real B_T2 = (real)0.25 * ( g_FC_B_In[TDir2][ idx_B_T2                                      ] +
//                                     g_FC_B_In[TDir2][ idx_B_T2 + LR*didx_B_T2[d]                    ] +
//                                     g_FC_B_In[TDir2][ idx_B_T2                   + didx_B_T2[TDir2] ] +
//                                     g_FC_B_In[TDir2][ idx_B_T2 + LR*didx_B_T2[d] + didx_B_T2[TDir2] ] );
#     endif // #ifdef MHD


//    2. update
//    fcCon[f][DENS] += ExtraTerm;

   } // for (int f=0; f<6; f++)

} // FUNCTION : Hydro_AddSourceTerm_HalfStep_MHM
#endif // #if ( FLU_SCHEME == MHM )



#if ( FLU_SCHEME == MHM_RP )
//-------------------------------------------------------------------------------------------------------
// Function    :  Hydro_AddSourceTerm_CCVar_HalfStep_MHM_RP
//
// Description :  Add source term on the cell-centered primitive varibles at half-step update
//
// Note        :  1. This function can only be used by MHM_RP
//                2. Invoked by Hydro_RiemannPredict()
//                3. Do not update magnetic field here, please update at Hydro_AddSourceTerm_FCVar_HalfStep_MHM_RP()
//
// Parameter   :  g_ConVar_In  : Array storing the input conserved variables
//                g_FC_B_In    : Array storing the input face-centered magnetic field (for MHD only)
//                OneCell      : Single-cell fluid array to store the updated cell-centered conserved variables
//                idx_in       : Index of accessing g_ConVar_In[]
//                didx_cvar_in : Index increment of g_ConVar_In[]
//                dt_dh2       : 0.5 * dt / dh
//                EoS          : EoS object
//
// Return      :  OneCell
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE
void Hydro_AddSourceTerm_CCVar_HalfStep_MHM_RP( const real g_ConVar_In[][ CUBE(FLU_NXT) ],
                                                const real g_FC_B_In[][ FLU_NXT_P1*SQR(FLU_NXT) ],
                                                      real OneCell[NCOMP_TOTAL_PLUS_MAG],
                                                const int idx_in, const int didx_cvar_in[3],
                                                const real dt_dh2, const EoS_t *EoS )
{

#  ifdef GAMER_DEBUG
   if ( didx_cvar_in[0] != 1  ||  didx_cvar_in[1] != FLU_NXT  ||  didx_cvar_in[2] != SQR(FLU_NXT) )
      printf( "ERROR : didx_cvar_in {%d, %d, %d} != {%d, %d, %d} !!\n", didx_cvar_in[0], didx_cvar_in[1], didx_cvar_in[2],
              1, FLU_NXT, SQR(FLU_NXT) );
#  endif

#  ifdef MHD
   const int i_in       = idx_in % FLU_NXT;
   const int j_in       = idx_in % SQR(FLU_NXT) / FLU_NXT;
   const int k_in       = idx_in / SQR(FLU_NXT);

// index of magnetic field
   const int idx_B_xL   = IDX321( i_in,   j_in,   k_in,   FLU_NXT_P1, FLU_NXT    );
   const int idx_B_xR   = IDX321( i_in+1, j_in,   k_in,   FLU_NXT_P1, FLU_NXT    );
   const int idx_B_yL   = IDX321( i_in,   j_in,   k_in,   FLU_NXT,    FLU_NXT_P1 );
   const int idx_B_yR   = IDX321( i_in,   j_in+1, k_in,   FLU_NXT,    FLU_NXT_P1 );
   const int idx_B_zL   = IDX321( i_in,   j_in,   k_in,   FLU_NXT,    FLU_NXT    );
   const int idx_B_zR   = IDX321( i_in,   j_in,   k_in+1, FLU_NXT,    FLU_NXT    );

// index increment of magnetic field
   const int didx_Bx[3] = { 1, FLU_NXT_P1, FLU_NXT_P1*FLU_NXT    };
   const int didx_By[3] = { 1, FLU_NXT,    FLU_NXT   *FLU_NXT_P1 };
   const int didx_Bz[3] = { 1, FLU_NXT,    FLU_NXT   *FLU_NXT    };
#  endif

// -----------------------------
// |    |                |     |
// |    |        ^       |     |
// |    |        |       |     |
// -------------B_yR------------
// |    |        |       |     |
// |    |                |     |
// | -B_xL->  OneCell  -B_xR-> |
// |    |                |     |
// |    |        ^       |     |
// |    |        |       |     |
// -------------B_yL------------
// |    |        |       |     |
// |    |                |     |
// -----------------------------

#  ifdef MHD
// magnetic field (see the figure)
// const real B_xL = g_FC_B_In[0][idx_B_xL];
// const real B_xR = g_FC_B_In[0][idx_B_xR];
// const real B_yL = g_FC_B_In[1][idx_B_yL];
// const real B_yR = g_FC_B_In[1][idx_B_yR];
// const real B_zL = g_FC_B_In[2][idx_B_zL];
// const real B_zR = g_FC_B_In[2][idx_B_zR];
#  endif

// example: cosmic-ray adiabatic work
// #  ifdef COSMIC_RAY
// // 1. calculate the cosmic-ray pressure
//    const real pCR_old = EoS->CREint2CRPres_FuncPtr( g_ConVar_In[CRAY][idx_in], EoS->AuxArrayDevPtr_Flt,
//                                                     EoS->AuxArrayDevPtr_Int, EoS->Table );
//
//
// // 2. compute \div V
//    real div_V[3];
//
//    for (int d=0; d<3; d++)
//    {
//       div_V[d] = (real)0.5 * ( g_ConVar_In[MOMX+d][idx_in + didx_cvar_in[d]] / g_ConVar_In[DENS][idx_in + didx_cvar_in[d]] -
//                                g_ConVar_In[MOMX+d][idx_in - didx_cvar_in[d]] / g_ConVar_In[DENS][idx_in - didx_cvar_in[d]] );
//    } // for (int d=0; d<3; d++)
//
//
// // 3. update the cosmic-ray energy
//    OneCell[CRAY] -= pCR_old*dt_dh2*( div_V[0] + div_V[1] + div_V[2] );
// #  endif // #ifdef COSMIC_RAY

} // FUNCTION : Hydro_AddSourceTerm_CCVar_HalfStep_MHM_RP



#ifdef MHD
//-------------------------------------------------------------------------------------------------------
// Function    :  Hydro_AddSourceTerm_FCVar_HalfStep_MHM_RP
//
// Description :  Add source term on the face-centered magnetic field at half-step update
//
// Note        :  1. This function can only be used by MHM_RP
//                2. Invoked by CPU/CUFLU_FluidSolver_MHM()
//                3. Must be used after MHD_UpdateMagnetic()
//
// Parameter   :  g_ConVar_In : Array storing the input conserved variables
//                g_FC_B_In   : Array storing the input     face-centered magnetic field (for MHD only)
//                g_FC_B_Half : Array storing the half-step face-centered magnetic field (for MHD only)
//                dt          : Time interval to advance solution
//                dh          : Cell size
//                EoS         : EoS object
//
// Return      :  g_FC_B_Half
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE
void Hydro_AddSourceTerm_FCVar_HalfStep_MHM_RP( const real g_ConVar_In[][ CUBE(FLU_NXT) ],
                                                const real g_FC_B_In[][ FLU_NXT_P1*SQR(FLU_NXT) ],
                                                      real g_FC_B_Half[][ FLU_NXT_P1*SQR(FLU_NXT) ],
                                                const real dt, const real dh, const EoS_t *EoS )
{

#  ifdef GAMER_DEBUG
   if ( FLU_NXT-N_HF_VAR != 2 )   printf( "FLU_NXT(%d) - N_HF_VAR(%d) != 2 !\n", FLU_NXT, N_HF_VAR );
#  endif

   const real dt_dh2          = (real)0.5 * dt / dh;
   const int  N_HF_VAR_P1     = N_HF_VAR + 1;
   const int  didx_cvar_in[3] = { 1, FLU_NXT, SQR(FLU_NXT) };

   for (int d=0; d<3; d++)
   {
      const int TDir1 = (d+1)%3;    // transverse direction 1
      const int TDir2 = (d+2)%3;    // transverse direction 2

      int size_BIn_i, size_BIn_j, size_BIn_k;
      int size_i, size_j;
      switch ( d )
      {
         case 0 : size_BIn_i = FLU_NXT_P1;  size_BIn_j = FLU_NXT;     size_BIn_k = FLU_NXT;
                  size_i     = N_HF_VAR_P1; size_j     = N_HF_VAR;    break;
         case 1 : size_BIn_i = FLU_NXT;     size_BIn_j = FLU_NXT_P1;  size_BIn_k = FLU_NXT;
                  size_i     = N_HF_VAR;    size_j     = N_HF_VAR_P1; break;
         case 2 : size_BIn_i = FLU_NXT;     size_BIn_j = FLU_NXT;     size_BIn_k = FLU_NXT_P1;
                  size_i     = N_HF_VAR;    size_j     = N_HF_VAR;    break;
      }
      const int didx_in_BN [3] = { 1, size_BIn_i, size_BIn_i*size_BIn_j };
      const int didx_in_BT1[3] = { 1, size_BIn_k, size_BIn_k*size_BIn_i };
      const int didx_in_BT2[3] = { 1, size_BIn_j, size_BIn_j*size_BIn_k };

      const int size_ij = size_i*size_j;
      CGPU_LOOP( idx_half, N_HF_VAR_P1*SQR(N_HF_VAR) )
      {
//       index of the half-step array
         const int i_half     = idx_half % size_i;
         const int j_half     = idx_half % size_ij / size_i;
         const int k_half     = idx_half / size_ij;

//       index of the cell-centered input variables
         const int i_in       = i_half + 1;
         const int j_in       = j_half + 1;
         const int k_in       = k_half + 1;
         const int idx_in     = IDX321( i_in, j_in, k_in, FLU_NXT, FLU_NXT );

//       index of the face-centered input variables
         const int i_BIn      = i_half + 1;
         const int j_BIn      = j_half + 1;
         const int k_BIn      = k_half + 1;
         const int idx_BIn_N  = IDX321( i_BIn, j_BIn, k_BIn, size_BIn_i, size_BIn_j );
         const int idx_BIn_T1 = IDX321( i_BIn, j_BIn, k_BIn, size_BIn_k, size_BIn_i );
         const int idx_BIn_T2 = IDX321( i_BIn, j_BIn, k_BIn, size_BIn_j, size_BIn_k );

//       |                      |              |
//       |      ^               |        ^     |
//       -------4------------------------3------
//       |      |               |        |     |
//       |                      |              |
//       |   idx_in             |              |
//       |      -            -BIn_N->  idx_in  |
//       |  didx_cvar_in[d]     |              |
//       |                      |              |
//       |      ^               |        ^     |
//       -------2------------------------1------
//       |      |               |        |     |
//       |                      |              |

//       1. calculate extra term
         real ExtraTerm = 0.0;
//       magnetic field at face-centered BIn_N
//       const real BIn_N  =                g_FC_B_In[d    ][ idx_BIn_N                                        ];
//       const real BIn_T1 = (real)0.25 * ( g_FC_B_In[TDir1][ idx_BIn_T1                                       ] +
//                                          g_FC_B_In[TDir1][ idx_BIn_T1 - didx_in_BT1[d]                      ] +
//                                          g_FC_B_In[TDir1][ idx_BIn_T1                  + didx_in_BT1[TDir1] ] +
//                                          g_FC_B_In[TDir1][ idx_BIn_T1 - didx_in_BT1[d] + didx_in_BT1[TDir1] ] );
//       const real BIn_T2 = (real)0.25 * ( g_FC_B_In[TDir2][ idx_BIn_T2                                       ] +
//                                          g_FC_B_In[TDir2][ idx_BIn_T2 - didx_in_BT2[d]                      ] +
//                                          g_FC_B_In[TDir2][ idx_BIn_T2                  + didx_in_BT2[TDir2] ] +
//                                          g_FC_B_In[TDir2][ idx_BIn_T2 - didx_in_BT2[d] + didx_in_BT2[TDir2] ] );

//       2. update
//       g_FC_B_Half[d][idx_half] += ExtraTerm;
      } // CGPU_LOOP( idx_half, N_HF_VAR_P1*SQR(N_HF_VAR) )
   } // for (int d=0; d<3; d++)

#  ifdef __CUDACC__
   __syncthreads();
#  endif

} // FUNCTION : Hydro_AddSourceTerm_FCVar_HalfStep_MHM_RP
#endif // #ifdef MHD
#endif // #if ( FLU_SCHEME == MHM_RP )



//-------------------------------------------------------------------------------------------------------
// Function    :  Hydro_AddSourceTerm_CCVar_FullStep
//
// Description :  Add source term on the cell-centered conserved varibles at full-step update
//
// Note        :  1. Shared by both MHM and MHM_RP
//                2. Invoked by Hydro_FullStepUpdate()
//
// Parameter   :  g_PriVar_Half : Array storing the input cell-centered primitive variables
//                                --> Accessed with the stride N_HF_VAR
//                                --> Although its actually allocated size is FLU_NXT^3 since it points to g_PriVar_1PG[]
//                OutCell       : Array to store the updated conserved variables
//                dt_dh         : Time interval to advance solution / cell size
//                EoS           : EoS object
//
// Return      :  OutputCell
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE
void Hydro_AddSourceTerm_CCVar_FullStep( const real g_PriVar_Half[][ CUBE(FLU_NXT) ],
                                         real OutCell[], const int idx_hf, const int didx_hf[3],
                                         const real dt_dh, const EoS_t *EoS )
{

// Example: cosmic ray adiabatic work done term
// #  ifdef COSMIC_RAY
// // 1. calculate the cosmic-ray pressure
//    const real pCR_half = EoS->CREint2CRPres_FuncPtr( g_PriVar_Half[CRAY][idx_hf], EoS->AuxArrayDevPtr_Flt,
//                                                      EoS->AuxArrayDevPtr_Int, EoS->Table );
//
//
// // 2. compute \div V
//    real div_V[3];
//    for (int d=0; d<3; d++)
//    {
//       div_V[d] = (real)0.5 * ( g_PriVar_Half[MOMX+d][idx_hf + didx_hf[d]] -
//                                g_PriVar_Half[MOMX+d][idx_hf - didx_hf[d]] );
//    } // for (int d=0; d<3; d++)
//
//
// // 3. update the cosmic-ray energy
//    OutCell[CRAY] -= pCR_half*dt_dh*( div_V[0] + div_V[1] + div_V[2] );
// #  endif

} // FUNCTION : Hydro_AddSourceTerm_CCVar_FullStep



#ifdef MHD
//-------------------------------------------------------------------------------------------------------
// Function    :  Hydro_AddSourceTerm_FCVar_FullStep
//
// Description :  Add source term on the face-centered magnetic field at full-step update
//
// Note        :  1. Shared by both MHM and MHM_RP
//                2. Invoked by CPU/CUFLU_FluidSolver_MHM()
//                3. Must be used after MHD_UpdateMagnetic()
//
// Parameter   :  g_PriVar_Half : Array storing the input cell-centered primitive variables
//                                --> Accessed with the stride N_HF_VAR
//                                --> Although its actually allocated size is FLU_NXT^3 since it points to g_PriVar_1PG[]
//                g_Output      : Array to store the updated fluid data
//                dt            : Time interval to advance solution
//                dh            : Cell size
//                EoS           : EoS object
//
// Return      :  g_FC_B_Out
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE
void Hydro_AddSourceTerm_FCVar_FullStep( const real g_PriVar_Half[][ CUBE(FLU_NXT) ],
                                               real g_FC_B_Out[][ PS2P1*SQR(PS2) ],
                                         const real dt, const real dh, const EoS_t *EoS )
{

   const real dt_dh      = dt/dh;
   const int  didx_hf[3] = { 1, N_HF_VAR, SQR(N_HF_VAR) };

   for (int d=0; d<3; d++)
   {
      int size_i, size_j;
      switch ( d )
      {
         case 0 : size_i = PS2P1; size_j = PS2;   break;
         case 1 : size_i = PS2;   size_j = PS2P1; break;
         case 2 : size_i = PS2;   size_j = PS2;   break;
      }

      const int size_ij = size_i*size_j;
      CGPU_LOOP( idx_out, PS2P1*SQR(PS2) )
      {
//       index of the output array
         const int i_out  = idx_out % size_i;
         const int j_out  = idx_out % size_ij / size_i;
         const int k_out  = idx_out / size_ij;

//       index of the half-step variables
         const int i_hf   = i_out + (N_HF_VAR-PS2)/2;
         const int j_hf   = j_out + (N_HF_VAR-PS2)/2;
         const int k_hf   = k_out + (N_HF_VAR-PS2)/2;
         const int idx_hf = IDX321( i_hf, j_hf, k_hf, N_HF_VAR, N_HF_VAR );

//       ------------------------------------
//       |                  |               |
//       |   idx_hf         |               |
//       |      -        idx_out    idx_hf  |
//       |  didx_hf[d]      |               |
//       |                  |               |
//       ------------------------------------

//       1. calculate extra term
         real ExtraTerm = (real)0.0;

//       2. update
//          g_FC_B_Out[d][idx_out] += ExtraTerm;
      } // CGPU_LOOP( idx_out, PS2P1*SQR(PS2) )
   } // for (int d=0; d<3; d++)

#  ifdef __CUDACC__
   __syncthreads();
#  endif

} // FUNCTION : Hydro_AddSourceTerm_FCVar_FullStep
#endif // #ifdef MHD



#endif // #if (  MODEL == HYDRO  &&  ( FLU_SCHEME == MHM  ||  FLU_SCHEME == MHM_RP )  )




#endif // #ifdef __CUFLU_HYDRO_ADDSOURCETERM__
