#include "CUFLU.h"

#ifdef COSMIC_RAY

// external functions
#ifdef __CUDACC__

#include "CUFLU_Shared_FluUtility.cu"

#else // #ifdef __CUDACC__

void Hydro_Con2Pri( const real In[], real Out[], const real MinPres,
                    const bool FracPassive, const int NFrac, const int FracIdx[],
                    const bool JeansMinPres, const real JeansMinPres_Coeff,
                    const EoS_DE2P_t EoS_DensEint2Pres, const EoS_DP2E_t EoS_DensPres2Eint,
                    const EoS_GUESS_t EoS_GuessHTilde, const EoS_H2TEM_t EoS_HTilde2Temp,
                    const double EoS_AuxArray_Flt[], const int EoS_AuxArray_Int[],
                    const real *const EoS_Table[EOS_NTABLE_MAX], real* const EintOut, real* LorentzFactorPtr );

#endif // #ifdef __CUDACC__ ... else ...



//-------------------------------------------------------------------------------------------------------
// Function    : CR_AdiabaticWork_HalfStep_MHM_RP
//
// Description : Add the adiabatic work term to update the cosmic-ray energy for the half-step solution of MHM_RP
//
// Note        : 1. MHM should not use this function
//               2. Work w/ and w/o MHD
//               3. Invoked by Hydro_RiemannPredict()
//
// Reference   : [1] Yang et al., ApJ 761, 185 (2012); doi:10.1088/0004-637X/761/2/185
//               [2] A simple dual implementation to track pressure accurately, S. Li, Astronum Proceeding, 385, 273 (2007)
//
// Parameter   : OneCell     : Single-cell fluid array to store the updated cell-centered cosmic-ray energy
//               g_ConVar_In : Array storing the input conserved variables
//               g_Flux_Half : Array storing the input face-centered fluxes
//                             --> Accessed with the stride didx_flux
//               idx_in      : Index of accessing g_ConVar_In[]
//               didx_in     : Index increment of g_ConVar_In[]
//               idx_flux    : Index of accessing g_flux_Half[]
//               didx_flux   : Index increment of g_Flux_Half[]
//               dt_dh2      : 0.5 * dt / dh
//               EoS         : EoS object
//
// Return      : OneCell[CRAY]
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE
void CR_AdiabaticWork_HalfStep_MHM_RP( real OneCell[NCOMP_TOTAL_PLUS_MAG],
                                       const real g_ConVar_In[][ CUBE(FLU_NXT) ],
                                       const real g_Flux_Half[][NCOMP_TOTAL_PLUS_MAG][ CUBE(N_FC_FLUX) ],
                                       const int idx_in, const int didx_in[3],
                                       const int idx_flux, const int didx_flux[3],
                                       const real dt_dh2, const EoS_t *EoS )
{
#  ifdef SRHD
// Although SRHD does not support the magnetic field yet, we still declare the size as 
// NCOMP_TOTAL_PLUS_MAG in case the magnetic field is supported someday.
   real Con_L[NCOMP_TOTAL_PLUS_MAG], Con_C[NCOMP_TOTAL_PLUS_MAG], Con_R[NCOMP_TOTAL_PLUS_MAG];
   real Pri_L[NCOMP_TOTAL_PLUS_MAG], Pri_C[NCOMP_TOTAL_PLUS_MAG], Pri_R[NCOMP_TOTAL_PLUS_MAG];
#  endif

// 1. compute \div V using the upwind data; reference: [2]
   real div_V[3], LorentzFactor = 1.0;

   for (int d=0; d<3; d++)
   {
#     ifdef MHD
      const real DensFlux_L = g_Flux_Half[d][DENS][ idx_flux - didx_flux[d] ];
      const real DensFlux_R = g_Flux_Half[d][DENS][ idx_flux                ];
#     else
      const real DensFlux_L = g_Flux_Half[d][DENS][ idx_flux                ];
      const real DensFlux_R = g_Flux_Half[d][DENS][ idx_flux + didx_flux[d] ];
#     endif

#     ifdef SRHD
      for (int v=0; v<NCOMP_TOTAL_PLUS_MAG; v++)
      {
         Con_L[v] = g_ConVar_In[v][ idx_in - didx_in[d] ];
         Con_C[v] = g_ConVar_In[v][ idx_in              ];
         Con_R[v] = g_ConVar_In[v][ idx_in + didx_in[d] ];
      }

      const real minPres = TINY_NUMBER;
      const real minJeansPres =  TINY_NUMBER;
      Hydro_Con2Pri( Con_L, Pri_L, minPres, NULL_BOOL, NULL_INT, NULL, NULL_BOOL, NULL_REAL,
                     EoS->DensEint2Pres_FuncPtr, EoS->DensPres2Eint_FuncPtr, EoS->GuessHTilde_FuncPtr,
                     EoS->HTilde2Temp_FuncPtr, EoS->AuxArrayDevPtr_Flt, EoS->AuxArrayDevPtr_Int,
                     EoS->Table, NULL, NULL );
      Hydro_Con2Pri( Con_C, Pri_C, minPres, NULL_BOOL, NULL_INT, NULL, NULL_BOOL, NULL_REAL,
                     EoS->DensEint2Pres_FuncPtr, EoS->DensPres2Eint_FuncPtr, EoS->GuessHTilde_FuncPtr,
                     EoS->HTilde2Temp_FuncPtr, EoS->AuxArrayDevPtr_Flt, EoS->AuxArrayDevPtr_Int,
                     EoS->Table, NULL, &LorentzFactor );
      Hydro_Con2Pri( Con_R, Pri_R, minPres, NULL_BOOL, NULL_INT, NULL, NULL_BOOL, NULL_REAL,
                     EoS->DensEint2Pres_FuncPtr, EoS->DensPres2Eint_FuncPtr, EoS->GuessHTilde_FuncPtr,
                     EoS->HTilde2Temp_FuncPtr, EoS->AuxArrayDevPtr_Flt, EoS->AuxArrayDevPtr_Int,
                     EoS->Table, NULL, NULL );

      const real Rho_L = Pri_L[DENS];
      const real Rho_C = Pri_C[DENS];
      const real Rho_R = Pri_R[DENS];
#     else
      const real Rho_L = g_ConVar_In[DENS][ idx_in - didx_in[d] ];
      const real Rho_C = g_ConVar_In[DENS][ idx_in              ];
      const real Rho_R = g_ConVar_In[DENS][ idx_in + didx_in[d] ];
#     endif

      div_V[d]  = ( DensFlux_R > (real)0.0 ) ? DensFlux_R / Rho_C : DensFlux_R / Rho_R;
      div_V[d] -= ( DensFlux_L > (real)0.0 ) ? DensFlux_L / Rho_L : DensFlux_L / Rho_C;
   } // for (int d=0; d<3; d++)

// 2. calculate the cosmic-ray pressure
   const real pCR_old = EoS->CREint2CRPres_FuncPtr( g_ConVar_In[CRAY][idx_in]/LorentzFactor, EoS->AuxArrayDevPtr_Flt,
                                                    EoS->AuxArrayDevPtr_Int, EoS->Table );

// 3. update the cosmic-ray energy
   OneCell[CRAY] -= pCR_old * dt_dh2 * ( div_V[0] + div_V[1] + div_V[2] );

} // FUMCTION : CR_AdiabaticWork_HalfStep_MHM_RP




//-------------------------------------------------------------------------------------------------------
// Function    : CR_AdiabaticWork_FullStep
//
// Description : Add the adiabatic work term to update the cosmic-ray energy for the full-step solution of MHM_RP/MHM
//
// Note        : 1. Shared by both MHM and MHM_RP (but it hasn't been tested for MHM yet)
//               2. Work w/ and w/o MHD
//               3. Invoked by CPU/CUFLU_FluidSolver_MHM()
//
// Reference   : [1] Yang et al., ApJ 761, 185 (2012); doi:10.1088/0004-637X/761/2/185
//               [2] A simple dual implementation to track pressure accurately, S. Li, Astronum Proceeding, 385, 273 (2007)
//
// Parameter   : g_PriVar_Half : Array storing the input cell-centered primitive variables
//                               --> Accessed with the stride N_HF_VAR
//                               --> Although its actually allocated size is FLU_NXT^3 since it points to g_PriVar_1PG[]
//               g_Output      : Array to store the updated fluid data
//               g_Flux        : Array storing the input face-centered fluxes
//                               --> Accessed with the array stride N_FL_FLUX even thought its actually
//                                   allocated size is N_FC_FLUX^3
//               g_FC_Var      : Array storing the input face-centered conserved variables
//                               --> Accessed with the array stride N_FC_VAR^3
//               dt            : Time interval to advance solution
//               dh            : Cell size
//               EoS           : EoS object
//
// Return      : g_Output[CRAY][]
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE
void CR_AdiabaticWork_FullStep( const real g_PriVar_Half[][ CUBE(FLU_NXT) ],
                                      real g_Output[][ CUBE(PS2) ],
                                const real g_Flux[][NCOMP_TOTAL_PLUS_MAG][ CUBE(N_FC_FLUX) ],
                                const real g_FC_Var[][NCOMP_TOTAL_PLUS_MAG][ CUBE(N_FC_VAR) ],
                                const real dt, const real dh, const EoS_t *EoS )
{

   const int  didx_flux[3] = { 1, N_FL_FLUX, SQR(N_FL_FLUX) };
   const int  didx_fc  [3] = { 1, N_FC_VAR,  SQR(N_FC_VAR)  };
   const real dt_dh        = dt/dh;

   real div_V[3];
#  ifdef SRHD
   real Con_LR[NCOMP_TOTAL_PLUS_MAG], Con_CL[NCOMP_TOTAL_PLUS_MAG], Con_CR[NCOMP_TOTAL_PLUS_MAG], Con_RL[NCOMP_TOTAL_PLUS_MAG];
   real Pri_LR[NCOMP_TOTAL_PLUS_MAG], Pri_CL[NCOMP_TOTAL_PLUS_MAG], Pri_CR[NCOMP_TOTAL_PLUS_MAG], Pri_RL[NCOMP_TOTAL_PLUS_MAG];
#  endif

   const int size_ij = SQR(PS2);
   CGPU_LOOP( idx_out, CUBE(PS2) )
   {
//    index of the output array
      const int i_out    = idx_out % PS2;
      const int j_out    = idx_out % size_ij / PS2;
      const int k_out    = idx_out / size_ij;

//    index of the flux array
//    --> for MHD, one additional flux is evaluated along each transverse direction for computing the CT electric field
#     ifdef MHD
      const int i_flux   = i_out + 1;
      const int j_flux   = j_out + 1;
      const int k_flux   = k_out + 1;
#     else
      const int i_flux   = i_out;
      const int j_flux   = j_out;
      const int k_flux   = k_out;
#     endif
      const int idx_flux = IDX321( i_flux, j_flux, k_flux, N_FL_FLUX, N_FL_FLUX );

//    index of the half-step variables
      const int i_hf     = i_out + (N_HF_VAR-PS2)/2;
      const int j_hf     = j_out + (N_HF_VAR-PS2)/2;
      const int k_hf     = k_out + (N_HF_VAR-PS2)/2;
      const int idx_hf   = IDX321( i_hf, j_hf, k_hf, N_HF_VAR, N_HF_VAR );

//    index of the face-centered variables
      const int i_fc     = i_out + 1;
      const int j_fc     = j_out + 1;
      const int k_fc     = k_out + 1;
      const int idx_fc   = IDX321( i_fc, j_fc, k_fc, N_FC_VAR, N_FC_VAR );

//    1. compute \div V using the upwind data; reference: [2]
      for (int d=0; d<3; d++)
      {
         const int faceL = 2*d;
         const int faceR = faceL+1;

#        ifdef MHD
         const real DensFlux_L = g_Flux[d][DENS][ idx_flux - didx_flux[d] ];
         const real DensFlux_R = g_Flux[d][DENS][ idx_flux                ];
#        else
         const real DensFlux_L = g_Flux[d][DENS][ idx_flux                ];
         const real DensFlux_R = g_Flux[d][DENS][ idx_flux + didx_flux[d] ];
#        endif

#        ifdef SRHD
         for (int v=0; v<NCOMP_TOTAL_PLUS_MAG; v++)
         {
            Con_LR[v] = g_FC_Var[faceR][v][ idx_fc - didx_fc[d] ];
            Con_CL[v] = g_FC_Var[faceL][v][ idx_fc              ];
            Con_CR[v] = g_FC_Var[faceR][v][ idx_fc              ];
            Con_RL[v] = g_FC_Var[faceL][v][ idx_fc + didx_fc[d] ];
         }

         const real minPres = TINY_NUMBER;
         const real minJeansPres =  TINY_NUMBER;
         Hydro_Con2Pri( Con_LR, Pri_LR, minPres, NULL_BOOL, NULL_INT, NULL, NULL_BOOL, NULL_REAL,
                        EoS->DensEint2Pres_FuncPtr, EoS->DensPres2Eint_FuncPtr, EoS->GuessHTilde_FuncPtr,
                        EoS->HTilde2Temp_FuncPtr, EoS->AuxArrayDevPtr_Flt, EoS->AuxArrayDevPtr_Int,
                        EoS->Table, NULL, NULL );
         Hydro_Con2Pri( Con_CL, Pri_CL, minPres, NULL_BOOL, NULL_INT, NULL, NULL_BOOL, NULL_REAL,
                        EoS->DensEint2Pres_FuncPtr, EoS->DensPres2Eint_FuncPtr, EoS->GuessHTilde_FuncPtr,
                        EoS->HTilde2Temp_FuncPtr, EoS->AuxArrayDevPtr_Flt, EoS->AuxArrayDevPtr_Int,
                        EoS->Table, NULL, NULL );
         Hydro_Con2Pri( Con_CR, Pri_CR, minPres, NULL_BOOL, NULL_INT, NULL, NULL_BOOL, NULL_REAL,
                        EoS->DensEint2Pres_FuncPtr, EoS->DensPres2Eint_FuncPtr, EoS->GuessHTilde_FuncPtr,
                        EoS->HTilde2Temp_FuncPtr, EoS->AuxArrayDevPtr_Flt, EoS->AuxArrayDevPtr_Int,
                        EoS->Table, NULL, NULL );
         Hydro_Con2Pri( Con_RL, Pri_RL, minPres, NULL_BOOL, NULL_INT, NULL, NULL_BOOL, NULL_REAL,
                        EoS->DensEint2Pres_FuncPtr, EoS->DensPres2Eint_FuncPtr, EoS->GuessHTilde_FuncPtr,
                        EoS->HTilde2Temp_FuncPtr, EoS->AuxArrayDevPtr_Flt, EoS->AuxArrayDevPtr_Int,
                        EoS->Table, NULL, NULL );

         const real Rho_LR = Pri_LR[DENS];
         const real Rho_CL = Pri_CL[DENS];
         const real Rho_CR = Pri_CR[DENS];
         const real Rho_RL = Pri_RL[DENS];
#        else
         const real Rho_LR = g_FC_Var[faceR][DENS][ idx_fc - didx_fc[d] ];
         const real Rho_CL = g_FC_Var[faceL][DENS][ idx_fc              ];
         const real Rho_CR = g_FC_Var[faceR][DENS][ idx_fc              ];
         const real Rho_RL = g_FC_Var[faceL][DENS][ idx_fc + didx_fc[d] ];
#        endif

         div_V[d]  = ( DensFlux_R > (real)0.0 ) ? DensFlux_R / Rho_CR : DensFlux_R / Rho_RL;
         div_V[d] -= ( DensFlux_L > (real)0.0 ) ? DensFlux_L / Rho_LR : DensFlux_L / Rho_CL;
      } // for (int d=0; d<3; d++)

//    2. calculate the cosmic-ray pressure
      const real pCR_half = EoS->CREint2CRPres_FuncPtr( g_PriVar_Half[CRAY][idx_hf], EoS->AuxArrayDevPtr_Flt,
                                                        EoS->AuxArrayDevPtr_Int, EoS->Table );

//    3. update the cosmic-ray energy
      g_Output[CRAY][idx_out] -= pCR_half * dt_dh * ( div_V[0] + div_V[1] + div_V[2] );

//    4. apply floor
      g_Output[CRAY][idx_out] = FMAX( g_Output[CRAY][idx_out], TINY_NUMBER );

   } // CGPU_LOOP( idx_out, CUBE(PS2) )

#  ifdef __CUDACC__
   __syncthreads();
#  endif

} // FUNCTION : CR_AdiabaticWork_FullStep



#endif // #ifdef COSMIC_RAY
