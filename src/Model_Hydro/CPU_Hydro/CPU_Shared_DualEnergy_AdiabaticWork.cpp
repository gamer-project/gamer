#ifndef __CUFLU_SHARED_DUALENERGY_ADIABATICWORK__
#define __CUFLU_SHARED_DUALENERGY_ADIABATICWORK__



#include "CUFLU.h"

#if ( MODEL == HYDRO  &&  DUAL_ENERGY == DE_EINT  &&  !defined SRHD )



// external functions
#ifdef __CUDACC__
#include "CUFLU_Shared_FluUtility.cu"
#include "CUFLU_Shared_DualEnergy.cu"
#endif




#if ( FLU_SCHEME == MHM_RP )
//-------------------------------------------------------------------------------------------------------
// Function    :  Hydro_DualEnergy_AdiabaticWork_HalfStep_MHM_RP
//
// Description :  Add the adiabatic work term to update the dual energy for the half-step solution of MHM_RP
//
// Note        :  1. MHM should not use this function
//                2. Work w/ and w/o MHD
//                3. Invoked by Hydro_RiemannPredict()
//
// Reference   :  [1] Bryan et al., ApJS 211, 19 (2012); doi:10.1088/0067-0049/211/2/19
//                [2] A simple dual implementation to track pressure accurately, S. Li, Astronum Proceeding, 385, 273 (2007)
//
// Parameter   :  OneCell     : Single-cell fluid array to store the updated cell-centered dual energy
//                g_ConVar_In : Array storing the input conserved variables
//                g_Flux_Half : Array storing the input face-centered fluxes
//                              --> Accessed with the stride didx_flux
//                idx_in      : Index of accessing g_ConVar_In[]
//                didx_in     : Index increment of g_ConVar_In[]
//                idx_flux    : Index of accessing g_flux_Half[]
//                didx_flux   : Index increment of g_Flux_Half[]
//                dt_dh2      : 0.5 * dt / dh
//                EoS         : EoS object
//
// Return      :  OneCell[DUAL]
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE
void Hydro_DualEnergy_AdiabaticWork_HalfStep_MHM_RP( real OneCell[NCOMP_TOTAL_PLUS_MAG],
                                                     const real g_ConVar_In[][ CUBE(FLU_NXT) ],
                                                     const real g_Flux_Half[][NCOMP_TOTAL_PLUS_MAG][ CUBE(N_FC_FLUX) ],
                                                     const int idx_in, const int didx_in[3],
                                                     const int idx_flux, const int didx_flux[3],
                                                     const real dt_dh2, const EoS_t *EoS )
{

// 1. calculate the pressure from the dual-energy variable
   real Passive[NCOMP_PASSIVE];
   for (int v=0; v<NCOMP_PASSIVE; v++)   Passive[v] = g_ConVar_In[ NCOMP_FLUID + v ][idx_in];

   const bool CheckMinPres_No = false;
   const real pDual_old = Hydro_DensDual2Pres( g_ConVar_In[DENS][idx_in], g_ConVar_In[DUAL][idx_in], Passive,
                                               CheckMinPres_No, NULL_REAL, EoS->DensEint2Pres_FuncPtr,
                                               EoS->AuxArrayDevPtr_Flt, EoS->AuxArrayDevPtr_Int, EoS->Table );


// 2. compute \div V using the upwind data; reference: [2]
   real div_V[3];

   for (int d=0; d<3; d++)
   {
#     ifdef MHD
      const real DensFlux_L = g_Flux_Half[d][DENS][ idx_flux - didx_flux[d] ];
      const real DensFlux_R = g_Flux_Half[d][DENS][ idx_flux                ];
#     else
      const real DensFlux_L = g_Flux_Half[d][DENS][ idx_flux                ];
      const real DensFlux_R = g_Flux_Half[d][DENS][ idx_flux + didx_flux[d] ];
#     endif

      div_V[d]  = ( DensFlux_R > (real)0.0 ) ?
                  ( DensFlux_R / g_ConVar_In[DENS][ idx_in              ] ) :
                  ( DensFlux_R / g_ConVar_In[DENS][ idx_in + didx_in[d] ] );

      div_V[d] -= ( DensFlux_L > (real)0.0 ) ?
                  ( DensFlux_L / g_ConVar_In[DENS][ idx_in - didx_in[d] ] ) :
                  ( DensFlux_L / g_ConVar_In[DENS][ idx_in              ] );
   } // for (int d=0; d<3; d++)


// 3. update the dual energy
   OneCell[DUAL] -= pDual_old*dt_dh2*( div_V[0] + div_V[1] + div_V[2] );

} // FUNCTION : Hydro_DualEnergy_AdiabaticWork_HalfStep_MHM_RP
#endif // #if ( FLU_SCHEME == MHM_RP )



#if ( FLU_SCHEME == MHM  ||  FLU_SCHEME == MHM_RP )
//-------------------------------------------------------------------------------------------------------
// Function    :  Hydro_DualEnergy_AdiabaticWork_FullStep
//
// Description :  Add the adiabatic work term to update the dual energy for the full-step solution of MHM/MHM_RP
//
// Note        :  1. Shared by both MHM and MHM_RP
//                2. Work w/ and w/o MHD
//                3. Invoked by Hydro_FullStepUpdate()
//
// Reference   :  [1] Bryan et al., ApJS 211, 19 (2012); doi:10.1088/0067-0049/211/2/19
//                [2] A simple dual implementation to track pressure accurately, S. Li, Astronum Proceeding, 385, 273 (2007)
//
// Parameter   :  Edual       : Dual energy to be updated
//                g_PriVar    : Array storing the input cell-centered primitive variables
//                              --> MHM without MHD: original-time data, stride FLU_NXT, offset FLU_GHOST_SIZE
//                              --> MHM with MHD   : half-step data, stride N_HF_VAR=N_FC_VAR, offset 1
//                              --> MHM_RP         : half-step data, stride N_HF_VAR=FLU_NXT-2, offset FLU_GHOST_SIZE-1
//                              --> Although its actual allocated size is FLU_NXT^3 since it points to g_PriVar_1PG[]
//                g_Flux      : Array storing the input face-centered fluxes
//                              --> Accessed with the array stride N_FL_FLUX even though its actually
//                                  allocated size is N_FC_FLUX^3
//                g_FC_Var    : Array storing the input face-centered conserved variables
//                              --> Accessed with the array stride N_FC_VAR^3
//                FracPassive : true --> input passive scalars are mass fraction instead of density
//                NFrac       : Number of passive scalars for the option "FracPassive"
//                FracIdx     : Target variable indices for the option "FracPassive"
//                dt          : Time interval to advance solution
//                dh          : Cell size
//                EoS         : EoS object
//                idx_out     : Array index associated with Edual
//
// Return      :  Edual
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE
void Hydro_DualEnergy_AdiabaticWork_FullStep( real &Edual,
                                              const real g_PriVar[][ CUBE(FLU_NXT) ],
                                              const real g_Flux[][NCOMP_TOTAL_PLUS_MAG][ CUBE(N_FC_FLUX) ],
                                              const real g_FC_Var[][NCOMP_TOTAL_PLUS_MAG][ CUBE(N_FC_VAR) ],
                                              const bool FracPassive, const int NFrac, const int FracIdx[],
                                              const real dt, const real dh, const EoS_t *EoS, const int idx_out )
{

   const int  size_ij      = SQR(PS2);
   const int  didx_flux[3] = { 1, N_FL_FLUX, SQR(N_FL_FLUX) };
   const int  didx_fc[3]   = { 1, N_FC_VAR,  SQR(N_FC_VAR)  };
   const real dt_dh        = dt/dh;

// index of the output array
   const int i_out    = idx_out % PS2;
   const int j_out    = idx_out % size_ij / PS2;
   const int k_out    = idx_out / size_ij;

// index of the flux array
// --> for MHD, one additional flux is evaluated along each transverse direction for computing the CT electric field
#  ifdef MHD
   const int i_flux   = i_out + 1;
   const int j_flux   = j_out + 1;
   const int k_flux   = k_out + 1;
#  else
   const int i_flux   = i_out;
   const int j_flux   = j_out;
   const int k_flux   = k_out;
#  endif
   const int idx_flux = IDX321( i_flux, j_flux, k_flux, N_FL_FLUX, N_FL_FLUX );

// index of the g_PriVar array
// --> both PLM and PPM retain the original layout for MHM without MHD
// --> MHM+MHD repacks all fields in Hydro_ConFC2PriCC_MHM(); MHM_RP uses Hydro_RiemannPredict()
#  if ( FLU_SCHEME == MHM  &&  !defined MHD )
   const int i_hf     = i_out + FLU_GHOST_SIZE;
   const int j_hf     = j_out + FLU_GHOST_SIZE;
   const int k_hf     = k_out + FLU_GHOST_SIZE;
   const int idx_hf   = IDX321( i_hf, j_hf, k_hf, FLU_NXT, FLU_NXT );
#  else // MHM_RP or MHM+MHD
   const int i_hf     = i_out + (N_HF_VAR-PS2)/2;
   const int j_hf     = j_out + (N_HF_VAR-PS2)/2;
   const int k_hf     = k_out + (N_HF_VAR-PS2)/2;
   const int idx_hf   = IDX321( i_hf, j_hf, k_hf, N_HF_VAR, N_HF_VAR );
#  endif

// index of the face-centered variables (N_FC_VAR=PS2+2 for all supported combinations)
   const int i_fc     = i_out + 1;
   const int j_fc     = j_out + 1;
   const int k_fc     = k_out + 1;
   const int idx_fc   = IDX321( i_fc, j_fc, k_fc, N_FC_VAR, N_FC_VAR );


// 1. calculate the pressure
   const bool CheckMinPres_No = false;
   real Passive[NCOMP_PASSIVE];

   for (int v=0; v<NCOMP_PASSIVE; v++)   Passive[v] = g_PriVar[ NCOMP_FLUID + v ][idx_hf];

// convert the mass fraction of target passive scalars to mass density
   if ( FracPassive )
      for (int v=0; v<NFrac; v++)   Passive[ FracIdx[v] ] *= g_PriVar[DENS][idx_hf];

   const real pDual_half = Hydro_DensDual2Pres( g_PriVar[DENS][idx_hf], g_PriVar[DUAL][idx_hf], Passive,
                                                CheckMinPres_No, NULL_REAL, EoS->DensEint2Pres_FuncPtr,
                                                EoS->AuxArrayDevPtr_Flt, EoS->AuxArrayDevPtr_Int, EoS->Table );


// 2. compute \div V using the upwind data; reference: [2]
   real div_V[3];

   for (int d=0; d<3; d++)
   {
      const int faceL = 2*d;
      const int faceR = faceL+1;

#     ifdef MHD
      const real DensFlux_L = g_Flux[d][DENS][ idx_flux - didx_flux[d] ];
      const real DensFlux_R = g_Flux[d][DENS][ idx_flux                ];
#     else
      const real DensFlux_L = g_Flux[d][DENS][ idx_flux                ];
      const real DensFlux_R = g_Flux[d][DENS][ idx_flux + didx_flux[d] ];
#     endif

      div_V[d]  = ( DensFlux_R > (real)0.0 ) ?
                  ( DensFlux_R / g_FC_Var[faceR][DENS][ idx_fc              ] ) :
                  ( DensFlux_R / g_FC_Var[faceL][DENS][ idx_fc + didx_fc[d] ] );

      div_V[d] -= ( DensFlux_L > (real)0.0 ) ?
                  ( DensFlux_L / g_FC_Var[faceR][DENS][ idx_fc - didx_fc[d] ] ) :
                  ( DensFlux_L / g_FC_Var[faceL][DENS][ idx_fc              ] );
   } // for (int d=0; d<3; d++)


// 3. update the dual energy
   Edual -= pDual_half*dt_dh*( div_V[0] + div_V[1] + div_V[2] );

} // FUNCTION : Hydro_DualEnergy_AdiabaticWork_FullStep
#endif // #if ( FLU_SCHEME == MHM  ||  FLU_SCHEME == MHM_RP )



#endif // #if ( MODEL == HYDRO  &&  DUAL_ENERGY == DE_EINT  &&  !defined SRHD )



#endif // #ifndef __CUFLU_SHARED_DUALENERGY_ADIABATICWORK__
