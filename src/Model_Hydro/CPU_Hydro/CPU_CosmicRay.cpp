#include "CUFLU.h"

#ifdef COSMIC_RAY




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

// 1. calculate the cosmic-ray pressure
   const real pCR_old = EoS->CREint2CRPres_FuncPtr( g_ConVar_In[CRAY][idx_in], EoS->AuxArrayDevPtr_Flt,
                                                    EoS->AuxArrayDevPtr_Int, EoS->Table );


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


// 3. update the cosmic-ray energy
   OneCell[CRAY] -= pCR_old*dt_dh2*( div_V[0] + div_V[1] + div_V[2] );

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

   const int didx_flux[3] = { 1, N_FL_FLUX, SQR(N_FL_FLUX) };
   const int didx_fc[3]   = { 1, N_FC_VAR,  SQR(N_FC_VAR)  };
   const real dt_dh       = dt/dh;

   real div_V[3];

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


//    1. calculate the cosmic-ray pressure
      const real pCR_half = EoS->CREint2CRPres_FuncPtr( g_PriVar_Half[CRAY][idx_hf], EoS->AuxArrayDevPtr_Flt,
                                                        EoS->AuxArrayDevPtr_Int, EoS->Table );


//    2. compute \div V using the upwind data; reference: [2]
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

         div_V[d]  = ( DensFlux_R > (real)0.0 ) ?
                     ( DensFlux_R / g_FC_Var[faceR][DENS][ idx_fc              ] ) :
                     ( DensFlux_R / g_FC_Var[faceL][DENS][ idx_fc + didx_fc[d] ] );

         div_V[d] -= ( DensFlux_L > (real)0.0 ) ?
                     ( DensFlux_L / g_FC_Var[faceR][DENS][ idx_fc - didx_fc[d] ] ) :
                     ( DensFlux_L / g_FC_Var[faceL][DENS][ idx_fc              ] );
      } // for (int d=0; d<3; d++)


//    3. update the cosmic-ray energy
      g_Output[CRAY][idx_out] -= pCR_half*dt_dh*( div_V[0] + div_V[1] + div_V[2] );

   } // CGPU_LOOP( idx_out, CUBE(PS2) )

#  ifdef __CUDACC__
   __syncthreads();
#  endif

} // FUNCTION : CR_AdiabaticWork_FullStep



#endif // #ifdef COSMIC_RAY
