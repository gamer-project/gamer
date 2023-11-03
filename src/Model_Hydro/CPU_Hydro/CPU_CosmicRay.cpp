#include "CUFLU.h"

#ifdef COSMIC_RAY



//-------------------------------------------------------------------------------------------------------
// Function    : CR_AdiabaticWork_HalfStepUpdate
//
// Description : Update the adiabatic work done term of cosmic ray
//
// Note        :
//
// Reference   : Yang, H.-Y.~K., Ruszkowski, M., Ricker, P.~M., et al. 2012, apj, 761, 185. doi:10.1088/0004-637X/761/2/185
//
// Parameter   :  OneCell     : Array storing the cell conserved variables
//                g_ConVar_In : Array storing the input conserved variables
//                g_Flux_Half : Array storing the input face-centered fluxes
//                              --> Accessed with the stride N_HF_FLUX
//                idx_flux    : Index of accessing g_flux_Half
//                didx_flux   : Index increment of g_Flux_Half
//                idx_fc      : Index of accessing g_ConVar_In
//                didx_fc     : Index increment of g_ConVar_In
//                dt_dh2      : 0.5 * dt / dh
//                EoS         : EoS object
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE
void CR_AdiabaticWork_HalfStepUpdate(       real OneCell[NCOMP_TOTAL_PLUS_MAG],
                                      const real g_ConVar_In[][ CUBE(FLU_NXT) ],
                                      const real g_Flux_Half[][NCOMP_TOTAL_PLUS_MAG][ CUBE(N_FC_FLUX) ],
                                      const int idx_fc, const int didx_fc[3],
                                      const int idx_flux, const int didx_flux[3],
                                      const real dt_dh2, const EoS_t *EoS )
{
   real div_V[3];

// 1. store the cosmic ray and calculate the pressure
   const real pCR_old = EoS->CREint2CRPres_FuncPtr( g_ConVar_In[CRAY][idx_fc], EoS->AuxArrayDevPtr_Flt,
                                                    EoS->AuxArrayDevPtr_Int, EoS->Table );

// 2. \div V term, reference: "Simple Method to Track Pressure Accurately", S. Li, Astronum Proceeding, 2007
   for (int d=0; d<3; d++)
   {
#     ifdef MHD
      div_V[d]  = ( g_Flux_Half[d][DENS][ idx_flux                ] > 0 ) ?
                  ( g_Flux_Half[d][DENS][ idx_flux                ] / g_ConVar_In[DENS][ idx_fc              ] ) :
                  ( g_Flux_Half[d][DENS][ idx_flux                ] / g_ConVar_In[DENS][ idx_fc + didx_fc[d] ] );

      div_V[d] -= ( g_Flux_Half[d][DENS][ idx_flux - didx_flux[d] ] > 0 ) ?
                  ( g_Flux_Half[d][DENS][ idx_flux - didx_flux[d] ] / g_ConVar_In[DENS][ idx_fc - didx_fc[d] ] ) :
                  ( g_Flux_Half[d][DENS][ idx_flux - didx_flux[d] ] / g_ConVar_In[DENS][ idx_fc              ] );
#     else
      div_V[d]  = ( g_Flux_Half[d][DENS][ idx_flux + didx_flux[d] ] > 0 ) ?
                  ( g_Flux_Half[d][DENS][ idx_flux + didx_flux[d] ] / g_ConVar_In[DENS][ idx_fc              ] ) :
                  ( g_Flux_Half[d][DENS][ idx_flux + didx_flux[d] ] / g_ConVar_In[DENS][ idx_fc + didx_fc[d] ] );

      div_V[d] -= ( g_Flux_Half[d][DENS][ idx_flux                ] > 0 ) ?
                  ( g_Flux_Half[d][DENS][ idx_flux                ] / g_ConVar_In[DENS][ idx_fc - didx_fc[d] ] ) :
                  ( g_Flux_Half[d][DENS][ idx_flux                ] / g_ConVar_In[DENS][ idx_fc              ] );
#     endif

   } // for (int d=0; d<3; d++)

// 3. update the cosmic ray
   OneCell[CRAY] = OneCell[CRAY] - pCR_old * dt_dh2 * ( div_V[0] + div_V[1] + div_V[2] );

} // FUMCTION : CR_AdiabaticWork_HalfStepUpdate




//-------------------------------------------------------------------------------------------------------
// Function    : CR_AdiabaticWork_FullStepUpdate
//
// Description : Update the adiabatic work done term of cosmic ray
//
// Note        :
//
// Reference   : Yang, H.-Y.~K., Ruszkowski, M., Ricker, P.~M., et al. 2012, apj, 761, 185. doi:10.1088/0004-637X/761/2/185
//
// Parameter   : g_PriVar_Half : Array to store the output primitive variables
//                               --> Accessed with the stride N_HF_VAR
//                               --> Although its actually allocated size is FLU_NXT^3 since it points to g_PriVar_1PG[]
//               g_Output      : Array to store the updated fluid data
//               g_Flux        : Array to store the output fluxes
//               g_FC_Var      : Array to store the half-step variables
//               dt            : Time interval to advance solution
//               dh            : Cell size
//               EoS           : EoS object
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE
void CR_AdiabaticWork_FullStepUpdate( const real g_PriVar_Half[][ CUBE(FLU_NXT) ],
                                            real g_Output[][ CUBE(PS2) ],
                                      const real g_Flux[][NCOMP_TOTAL_PLUS_MAG][ CUBE(N_FC_FLUX) ],
                                      const real g_FC_Var[][NCOMP_TOTAL_PLUS_MAG][ CUBE(N_FC_VAR) ],
                                      const real dt, const real dh, const EoS_t *EoS )
{
   const int didx_flux[3] = { 1, N_FL_FLUX, SQR(N_FL_FLUX) };
   const int didx_fc[3]   = { 1, N_FC_VAR,  SQR(N_FC_VAR)  };
   const real dt_dh       = dt/dh;

   real div_V[3], Output_1Cell[NCOMP_TOTAL_PLUS_MAG];

   const int size_ij = SQR(PS2);
   CGPU_LOOP( idx_out, CUBE(PS2) )
   {
      const int i_out    = idx_out % PS2;
      const int j_out    = idx_out % size_ij / PS2;
      const int k_out    = idx_out / size_ij;

//    for MHD, one additional flux is evaluated along each transverse direction for computing the CT electric field
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

//    index of the half step variable
      const int i_hf     = i_out + (N_HF_VAR-PS2)/2;
      const int j_hf     = j_out + (N_HF_VAR-PS2)/2;
      const int k_hf     = k_out + (N_HF_VAR-PS2)/2;
      const int idx_hf   = IDX321( i_hf, j_hf, k_hf, N_HF_VAR, N_HF_VAR );

//    index for the face center velocity
      const int i_fc     = i_out + 1;
      const int j_fc     = j_out + 1;
      const int k_fc     = k_out + 1;
      const int idx_fc   = IDX321( i_fc, j_fc, k_fc, N_FC_VAR, N_FC_VAR);

//    1. store the cosmic ray and calculate the pressure
      Output_1Cell[CRAY]  = g_Output[CRAY][idx_out];
      const real pCR_half = EoS->CREint2CRPres_FuncPtr( g_PriVar_Half[CRAY][idx_hf], EoS->AuxArrayDevPtr_Flt,
                                                        EoS->AuxArrayDevPtr_Int, EoS->Table );

//    2. \div V term, reference: "Simple Method to Track Pressure Accurately", S. Li, Astronum Proceeding, 2007
      for (int d=0; d<3; d++)
      {
         const int faceL = 2*d;
         const int faceR = faceL+1;

#        ifdef MHD
         div_V[d]  = ( g_Flux[d][DENS][ idx_flux                ] > 0 ) ?
                     ( g_Flux[d][DENS][ idx_flux                ] / g_FC_Var[faceR][DENS][ idx_fc              ] ) :
                     ( g_Flux[d][DENS][ idx_flux                ] / g_FC_Var[faceL][DENS][ idx_fc + didx_fc[d] ] );

         div_V[d] -= ( g_Flux[d][DENS][ idx_flux - didx_flux[d] ] > 0 ) ?
                     ( g_Flux[d][DENS][ idx_flux - didx_flux[d] ] / g_FC_Var[faceR][DENS][ idx_fc - didx_fc[d] ] ) :
                     ( g_Flux[d][DENS][ idx_flux - didx_flux[d] ] / g_FC_Var[faceL][DENS][ idx_fc              ] );
#        else
         div_V[d]  = ( g_Flux[d][DENS][ idx_flux + didx_flux[d] ] > 0 ) ?
                     ( g_Flux[d][DENS][ idx_flux + didx_flux[d] ] / g_FC_Var[faceR][DENS][ idx_fc              ] ) :
                     ( g_Flux[d][DENS][ idx_flux + didx_flux[d] ] / g_FC_Var[faceL][DENS][ idx_fc + didx_fc[d] ] );

         div_V[d] -= ( g_Flux[d][DENS][ idx_flux                ] > 0 ) ?
                     ( g_Flux[d][DENS][ idx_flux                ] / g_FC_Var[faceR][DENS][ idx_fc - didx_fc[d] ] ) :
                     ( g_Flux[d][DENS][ idx_flux                ] / g_FC_Var[faceL][DENS][ idx_fc              ] );
#        endif

      } // for (int d=0; d<3; d++)

//    3. update the cosmic ray
      Output_1Cell[CRAY] = Output_1Cell[CRAY] - pCR_half * dt_dh * ( div_V[0] + div_V[1] + div_V[2] );
      g_Output[CRAY][idx_out] = Output_1Cell[CRAY];

      } // CGPU_LOOP( idx_out, CUBE(PS2) )

#  ifdef __CUDACC__
   __syncthreads();
#  endif

} // FUNCTION : CR_AdiabaticWork_FullStepUpdate



#endif // #ifdef COSMIC_RAY
