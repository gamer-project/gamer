#ifndef __CUFLU_HYDRO_SCANALLCELLS__
#define __CUFLU_HYDRO_SCANALLCELLS__




#include "CUFLU.h"



#if (  MODEL == HYDRO  &&  ( FLU_SCHEME == MHM  ||  FLU_SCHEME == MHM_RP )  )



#if ( FLU_SCHEME == MHM )
//-------------------------------------------------------------------------------------------------------
// Function    :  Hydro_Scan_HalfStep_MHM
//
// Description :  Scan through all face-centered conserved varibles after half-step update
//
// Note        :  1. This function can only be used by MHM
//                2. Invoked by Hydro_HancockPredict()
//
// Parameter   :  fcCon       : Face-centered conserved variables to be updated
//                g_ConVar_In : Array storing the input conserved variables
//                g_FC_B_In   : Array storing the input face-centered magnetic field (for MHD only)
//                idx_in      : Index of accessing g_ConVar_In[]
//                didx_in     : Index increment of g_ConVar_In[]
//                dt_dh2      : 0.5 * dt / dh
//                EoS         : EoS object
//
// Return      :  fcCon
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE
void Hydro_Scan_HalfStep_MHM( const real g_ConVar_In[][ CUBE(FLU_NXT) ],
                              const real g_FC_B_In[][ FLU_NXT_P1*SQR(FLU_NXT) ],
                                    real fcCon[][NCOMP_LR],
                              const int idx_in, const int didx_in[3],
                              const real dt_dh2, const EoS_t *EoS )
{

#  ifdef GAMER_DEBUG
   if ( didx_in[0] != 1  ||  didx_in[1] != FLU_NXT  ||  didx_in[2] != SQR(FLU_NXT) )
      printf( "didx_in {%d, %d, %d} != {%d, %d, %d} !\n", didx_in[0], didx_in[1], didx_in[2],
              1, FLU_NXT, SQR(FLU_NXT) );
#  endif

   const int idx_i = idx_in % FLU_NXT;
   const int idx_j = idx_in % SQR(FLU_NXT) / FLU_NXT;
   const int idx_k = idx_in / SQR(FLU_NXT);

#  ifdef MHD
// index of magnetic field
   const int idx_B_xL   = IDX321( idx_i,   idx_j,   idx_k,   FLU_NXT_P1, FLU_NXT    );
   const int idx_B_xR   = IDX321( idx_i+1, idx_j,   idx_k,   FLU_NXT_P1, FLU_NXT    );
   const int idx_B_yL   = IDX321( idx_i,   idx_j,   idx_k,   FLU_NXT,    FLU_NXT_P1 );
   const int idx_B_yR   = IDX321( idx_i,   idx_j+1, idx_k,   FLU_NXT,    FLU_NXT_P1 );
   const int idx_B_zL   = IDX321( idx_i,   idx_j,   idx_k,   FLU_NXT,    FLU_NXT    );
   const int idx_B_zR   = IDX321( idx_i,   idx_j,   idx_k+1, FLU_NXT,    FLU_NXT    );
   const int didx_Bx[3] = { 1, FLU_NXT_P1, FLU_NXT_P1*FLU_NXT    };
   const int didx_By[3] = { 1, FLU_NXT,    FLU_NXT   *FLU_NXT_P1 };
   const int didx_Bz[3] = { 1, FLU_NXT,    FLU_NXT   *FLU_NXT    };
#  endif

// 1. calculate extra term
   real B_xL = g_FC_B_In[0][idx_B_xL];
   real B_xR = g_FC_B_In[0][idx_B_xR];
   real B_yL = g_FC_B_In[1][idx_B_yL];
   real B_yR = g_FC_B_In[1][idx_B_yR];
   real B_zL = g_FC_B_In[2][idx_B_zL];
   real B_zR = g_FC_B_In[2][idx_B_zR];

   real ExtraTerm = 0.0;
   ExtraTerm = g_ConVar_In[DENS][idx_in];

// 2. update
   fcCon[0][DENS] += ExtraTerm;

} // FUNCTION : Hydro_Scan_HalfStep_MHM
#endif // #if ( FLU_SCHEME == MHM )



#if ( FLU_SCHEME == MHM_RP )
//-------------------------------------------------------------------------------------------------------
// Function    :  Hydro_Scan_HalfStep_MHM_RP
//
// Description :  Scan through all cell-centered primitive varibles after half-step update
//
// Note        :  1. This function can only be used by MHM_RP
//                2. Invoked by Hydro_RiemannPredict()
//                3. Must Invoked after the cell-centered magnetic field calculation
//
// Parameter   :  g_ConVar_In : Array storing the input conserved variables
//                g_FC_B_In   : Array storing the input face-centered magnetic field (for MHD only)
//                OneCell     : Single-cell fluid array to store the updated cell-centered primitive variables
//                idx_in      : Index of accessing g_ConVar_In[]
//                didx_in     : Index increment of g_ConVar_In[]
//                dt_dh2      : 0.5 * dt / dh
//                EoS         : EoS object
//
// Return      :  OneCell
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE
void Hydro_Scan_HalfStep_MHM_RP( const real g_ConVar_In[][ CUBE(FLU_NXT) ],
                                 const real g_FC_B_In[][ FLU_NXT_P1*SQR(FLU_NXT) ],
                                       real OneCell[NCOMP_TOTAL_PLUS_MAG],
                                 const int idx_in, const int didx_in[3],
                                 const real dt_dh2, const EoS_t *EoS )
{

#  ifdef GAMER_DEBUG
   if ( didx_in[0] != 1  ||  didx_in[1] != FLU_NXT  ||  didx_in[2] != SQR(FLU_NXT) )
      printf( "didx_in {%d, %d, %d} != {%d, %d, %d} !\n", didx_in[0], didx_in[1], didx_in[2],
              1, FLU_NXT, SQR(FLU_NXT) );
#  endif

   const int idx_i = idx_in % FLU_NXT;
   const int idx_j = idx_in % SQR(FLU_NXT) / FLU_NXT;
   const int idx_k = idx_in / SQR(FLU_NXT);

#  ifdef MHD
// index of magnetic field
   const int idx_B_xL   = IDX321( idx_i,   idx_j,   idx_k,   FLU_NXT_P1, FLU_NXT    );
   const int idx_B_xR   = IDX321( idx_i+1, idx_j,   idx_k,   FLU_NXT_P1, FLU_NXT    );
   const int idx_B_yL   = IDX321( idx_i,   idx_j,   idx_k,   FLU_NXT,    FLU_NXT_P1 );
   const int idx_B_yR   = IDX321( idx_i,   idx_j+1, idx_k,   FLU_NXT,    FLU_NXT_P1 );
   const int idx_B_zL   = IDX321( idx_i,   idx_j,   idx_k,   FLU_NXT,    FLU_NXT    );
   const int idx_B_zR   = IDX321( idx_i,   idx_j,   idx_k+1, FLU_NXT,    FLU_NXT    );
   const int didx_Bx[3] = { 1, FLU_NXT_P1, FLU_NXT_P1*FLU_NXT    };
   const int didx_By[3] = { 1, FLU_NXT,    FLU_NXT   *FLU_NXT_P1 };
   const int didx_Bz[3] = { 1, FLU_NXT,    FLU_NXT   *FLU_NXT    };
#  endif

#  ifdef COSMIC_RAY
// 1. calculate the cosmic-ray pressure
   const real pCR_old = EoS->CREint2CRPres_FuncPtr( g_ConVar_In[CRAY][idx_in], EoS->AuxArrayDevPtr_Flt,
                                                    EoS->AuxArrayDevPtr_Int, EoS->Table );


// 2. compute \div V using the upwind data; reference: [2]
   real div_V[3];

   for (int d=0; d<3; d++)
   {
      div_V[d] = (real)0.5 * ( g_ConVar_In[DENS+d][idx_in + didx_in[d]] / g_ConVar_In[DENS][idx_in + didx_in[d]] -
                               g_ConVar_In[DENS+d][idx_in - didx_in[d]] / g_ConVar_In[DENS][idx_in - didx_in[d]] );
   } // for (int d=0; d<3; d++)


// 3. update the cosmic-ray energy
   OneCell[CRAY] -= pCR_old*dt_dh2*( div_V[0] + div_V[1] + div_V[2] );
#  endif // #ifdef COSMIC_RAY

} // FUMCTION : Hydro_Scan_HalfStep_MHM_RP
#endif // #if ( FLU_SCHEME == MHM_RP )




//-------------------------------------------------------------------------------------------------------
// Function    :  Hydro_Scan_CCVar_FullStep
//
// Description :  Scan through all cell-centered conserved varibles after full-step update
//
// Note        :  1. Shared by both MHM and MHM_RP
//                2. Invoked by CPU/CUFLU_FluidSolver_MHM()
//
// Parameter   :  g_PriVar_Half : Array storing the input cell-centered primitive variables
//                                --> Accessed with the stride N_HF_VAR
//                                --> Although its actually allocated size is FLU_NXT^3 since it points to g_PriVar_1PG[]
//                g_Output      : Array to store the updated fluid data
//                dt            : Time interval to advance solution
//                dh            : Cell size
//                EoS           : EoS object
//
// Return      :  g_Output
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE
void Hydro_Scan_CCVar_FullStep( const real g_PriVar_Half[][ CUBE(FLU_NXT) ],
                                      real g_Output[][ CUBE(PS2) ],
                                const real dt, const real dh, const EoS_t *EoS )
{

#  ifdef GAMER_DEBUG
#  endif

   const real dt_dh      = dt/dh;
   const int  didx_hf[3] = { 1, N_HF_VAR, SQR(N_HF_VAR) };
   const int  size_ij    = SQR(PS2);

   CGPU_LOOP( idx_out, CUBE(PS2) )
   {
//    index of the output array
      const int i_out    = idx_out % PS2;
      const int j_out    = idx_out % size_ij / PS2;
      const int k_out    = idx_out / size_ij;

//    index of the half-step variables
      const int i_hf     = i_out + (N_HF_VAR-PS2)/2;
      const int j_hf     = j_out + (N_HF_VAR-PS2)/2;
      const int k_hf     = k_out + (N_HF_VAR-PS2)/2;
      const int idx_hf   = IDX321( i_hf, j_hf, k_hf, N_HF_VAR, N_HF_VAR );

#     ifdef COSMIC_RAY
//    1. calculate the cosmic-ray pressure
      const real pCR_half = EoS->CREint2CRPres_FuncPtr( g_PriVar_Half[CRAY][idx_hf], EoS->AuxArrayDevPtr_Flt,
                                                        EoS->AuxArrayDevPtr_Int, EoS->Table );


//    2. compute \div V using the upwind data; reference: [2]
      real div_V[3];
      for (int d=0; d<3; d++)
      {
         div_V[d] = (real)0.5 * ( g_PriVar_Half[DENS+d][idx_hf + didx_hf[d]] -
                                  g_PriVar_Half[DENS+d][idx_hf - didx_hf[d]] );
      } // for (int d=0; d<3; d++)


//    3. update the cosmic-ray energy
      g_Output[CRAY][idx_out] -= pCR_half*dt_dh*( div_V[0] + div_V[1] + div_V[2] );
#     endif

   } // CGPU_LOOP( idx_out, CUBE(PS2) )

#  ifdef __CUDACC__
   __syncthreads();
#  endif

} // FUNCTION : Hydro_Scan_CCVar_FullStep



//-------------------------------------------------------------------------------------------------------
// Function    :  Hydro_Scan_FCVar_FullStep
//
// Description :  Scan through all face-centered magnetic field after full-step update
//
// Note        :  1. Shared by both MHM and MHM_RP
//                2. Invoked by CPU/CUFLU_FluidSolver_MHM()
//
// Parameter   :  g_PriVar_Half : Array storing the input cell-centered primitive variables
//                                --> Accessed with the stride N_HF_VAR
//                                --> Although its actually allocated size is FLU_NXT^3 since it points to g_PriVar_1PG[]
//                g_Output      : Array to store the updated fluid data
//                dt            : Time interval to advance solution
//                dh            : Cell size
//                EoS         : EoS object
//
// Return      :  g_FC_B_Out
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE
void Hydro_Scan_FCVar_FullStep( const real g_PriVar_Half[][ CUBE(FLU_NXT) ],
                                      real g_FC_B_Out[][ PS2P1*SQR(PS2) ],
                                const real dt, const real dh, const EoS_t *EoS )
{

#  ifdef GAMER_DEBUG
#  endif

   const real dt_dh      = dt/dh;
   const int  didx_hf[3] = { 1, N_HF_VAR, SQR(N_HF_VAR) };

   for (int d=0; d<3; d++)
   {
      const int TDir1 = (d+1)%3;    // transverse direction 1
      const int TDir2 = (d+2)%3;    // transverse direction 2

      int size_i, size_j, size_k;

      switch ( d )
      {
         case 0 : size_i = PS2P1; size_j = PS2;   size_k = PS2;
                  break;

         case 1 : size_i = PS2;   size_j = PS2P1; size_k = PS2;
                  break;

         case 2 : size_i = PS2;   size_j = PS2;   size_k = PS2P1;
                  break;
      }

      const int size_ij = size_i*size_j;
      CGPU_LOOP( idx_out, PS2P1*SQR(PS2) )
      {
//       index of the output array
         const int i_out    = idx_out % size_i;
         const int j_out    = idx_out % size_ij / size_i;
         const int k_out    = idx_out / size_ij;

//       index of the half-step variables
         const int i_hf     = i_out + (N_HF_VAR-PS2)/2;
         const int j_hf     = j_out + (N_HF_VAR-PS2)/2;
         const int k_hf     = k_out + (N_HF_VAR-PS2)/2;
         const int idx_hf   = IDX321( i_hf, j_hf, k_hf, N_HF_VAR, N_HF_VAR );

//       1. calculate extra term
         real ExtraTerm = 0.0;
         ExtraTerm = g_PriVar_Half[DENS][idx_hf];


//       2. update
         // g_FC_B_Out[d][idx_out] += ExtraTerm;
      } // CGPU_LOOP( idx_out, PS2P1*SQR(PS2) )
   } // for (int d=0; d<3; d++)

#  ifdef __CUDACC__
   __syncthreads();
#  endif

} // FUNCTION : Hydro_Scan_FCVar_FullStep



#endif // #if (  MODEL == HYDRO  &&  ( FLU_SCHEME == MHM  ||  FLU_SCHEME == MHM_RP )  )




#endif // #ifdef __CUFLU_HYDRO_SCANALLCELLS__
