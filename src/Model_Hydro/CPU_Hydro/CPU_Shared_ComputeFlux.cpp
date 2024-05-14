#ifndef __CUFLU_COMPUTEFLUX__
#define __CUFLU_COMPUTEFLUX__



#include "CUFLU.h"

#if ( MODEL == HYDRO  &&  (FLU_SCHEME == MHM || FLU_SCHEME == MHM_RP || FLU_SCHEME == CTU) )



// external functions
#ifdef __CUDACC__

#if ( RSOLVER == EXACT  ||  RSOLVER_RESCUE == EXACT )
# include "CUFLU_Shared_RiemannSolver_Exact.cu"
#endif
#if ( RSOLVER == ROE    ||  RSOLVER_RESCUE == ROE   )
# include "CUFLU_Shared_RiemannSolver_Roe.cu"
#endif
#if ( RSOLVER == HLLE   ||  RSOLVER_RESCUE == HLLE  )
# include "CUFLU_Shared_RiemannSolver_HLLE.cu"
#endif
#if ( RSOLVER == HLLC   ||  RSOLVER_RESCUE == HLLC  )
# include "CUFLU_Shared_RiemannSolver_HLLC.cu"
#endif
#if ( RSOLVER == HLLD   ||  RSOLVER_RESCUE == HLLD  )
# include "CUFLU_Shared_RiemannSolver_HLLD.cu"
#endif

#else // #ifdef __CUDACC__

#if ( RSOLVER == EXACT  ||  RSOLVER_RESCUE == EXACT )
void Hydro_RiemannSolver_Exact( const int XYZ, real Flux_Out[], const real L_In[], const real R_In[],
                                const real MinDens, const real MinPres, const EoS_DE2P_t EoS_DensEint2Pres,
                                const EoS_DP2C_t EoS_DensPres2CSqr, const double EoS_AuxArray_Flt[],
                                const int EoS_AuxArray_Int[], const real* const EoS_Table[EOS_NTABLE_MAX] );
#endif
#if ( RSOLVER == ROE    ||  RSOLVER_RESCUE == ROE   )
void Hydro_RiemannSolver_Roe( const int XYZ, real Flux_Out[], const real L_In[], const real R_In[],
                              const real MinDens, const real MinPres, const EoS_DE2P_t EoS_DensEint2Pres,
                              const EoS_DP2C_t EoS_DensPres2CSqr, const double EoS_AuxArray_Flt[],
                              const int EoS_AuxArray_Int[], const real* const EoS_Table[EOS_NTABLE_MAX] );
#endif
#if ( RSOLVER == HLLE   ||  RSOLVER_RESCUE == HLLE  )
void Hydro_RiemannSolver_HLLE( const int XYZ, real Flux_Out[], const real L_In[], const real R_In[],
                               const real MinDens, const real MinPres, const EoS_DE2P_t EoS_DensEint2Pres,
                               const EoS_DP2C_t EoS_DensPres2CSqr, const EoS_GUESS_t EoS_GuessHTilde,
                               const EoS_H2TEM_t EoS_HTilde2Temp,
                               const double EoS_AuxArray_Flt[], const int EoS_AuxArray_Int[],
                               const real* const EoS_Table[EOS_NTABLE_MAX] );
#endif
#if ( RSOLVER == HLLC   ||  RSOLVER_RESCUE == HLLC  )
void Hydro_RiemannSolver_HLLC( const int XYZ, real Flux_Out[], const real L_In[], const real R_In[],
                               const real MinDens, const real MinPres, const EoS_DE2P_t EoS_DensEint2Pres,
                               const EoS_DP2C_t EoS_DensPres2CSqr, const EoS_GUESS_t EoS_GuessHTilde,
                               const EoS_H2TEM_t EoS_HTilde2Temp,
                               const double EoS_AuxArray_Flt[], const int EoS_AuxArray_Int[],
                               const real* const EoS_Table[EOS_NTABLE_MAX] );
#endif
#if ( RSOLVER == HLLD   ||  RSOLVER_RESCUE == HLLD  )
void Hydro_RiemannSolver_HLLD( const int XYZ, real Flux_Out[], const real L_In[], const real R_In[],
                               const real MinDens, const real MinPres, const EoS_DE2P_t EoS_DensEint2Pres,
                               const EoS_DP2C_t EoS_DensPres2CSqr, const double EoS_AuxArray_Flt[],
                               const int EoS_AuxArray_Int[], const real* const EoS_Table[EOS_NTABLE_MAX] );
#endif

#endif // #ifdef __CUDACC__ ... else ...




//-------------------------------------------------------------------------------------------------------
// Function    :  Hydro_ComputeFlux
// Description :  Compute the face-centered fluxes by Riemann solver
//
// Note        :  1. Currently support the exact, HLLC, HLLE, HLLD, and Roe solvers
//                2. g_FC_Var[] has the size of N_FC_VAR^3
//                   --> (N_FC_VAR-1-2*NSkip_N)*(N_FC_VAR-2*NSkip_T)^2 fluxes will be computed
//                   --> See below for the definitions of NSkip_N and NSkip_T
//                3. g_FC_Flux[] has the size of N_FC_FLUX^3
//                   --> But (i,j,k) flux will be stored in the "(k*NFlux+j)*NFlux+i" element in g_FC_Flux[]
//                       --> We have assumed that NFlux <= N_FC_FLUX
//                   --> (i,j,k) in g_FC_Flux_x[] is defined on the +x surface of the cell (i+NSkip_N, j+NSkip_T, k+NSkip_T) in g_FC_Var[]
//                       (i,j,k) in g_FC_Flux_y[] is defined on the +y surface of the cell (i+NSkip_T, j+NSkip_N, k+NSkip_T) in g_FC_Var[]
//                       (i,j,k) in g_FC_Flux_z[] is defined on the +z surface of the cell (i+NSkip_T, j+NSkip_T, k+NSkip_N) in g_FC_Var[]
//                4. This function is shared by MHM, MHM_RP, and CTU schemes
//                5. For the unsplitting scheme in gravity (i.e., UNSPLIT_GRAVITY), this function also corrects the half-step
//                   velocity by gravity when CorrHalfVel==true
//
// Parameter   :  g_FC_Var        : Array storing the input face-centered conserved variables
//                g_FC_Flux       : Array to store the output face-centered fluxes
//                NFlux           : Stride for accessing g_FC_Flux[]
//                NSkip_N         : Number of cells to be skipped in the normal directions
//                                  --> "(N_FC_VAR-1-2*NSkip_N)" fluxes will be computed along the normal direction
//                NSkip_T         : Number of cells to be skipped in the transverse directions
//                                  --> "(N_FC_VAR-2*NSkip_T)^2" fluxes will be computed along the transverse direction
//                CorrHalfVel     : true --> correct the half-step velocity by gravity       (for UNSPLIT_GRAVITY only)
//                g_Pot_USG       : Array storing the input potential for CorrHalfVel        (for UNSPLIT_GRAVITY only)
//                g_Corner        : Array storing the corner coordinates of each patch group (for UNSPLIT_GRAVITY only)
//                dt              : Time interval to advance the full-step solution          (for UNSPLIT_GRAVITY only)
//                dh              : Cell size                                                (for UNSPLIT_GRAVITY only)
//                Time            : Current physical time                                    (for UNSPLIT_GRAVITY only)
//                UsePot          : Add self-gravity and/or external potential               (for UNSPLIT_GRAVITY only)
//                ExtAcc          : Add external acceleration                                (for UNSPLIT_GRAVITY only)
//                ExtAcc_Func     : Function pointer to the external acceleration routine    (for UNSPLIT_GRAVITY only)
//                ExtAcc_AuxArray : Auxiliary array for external acceleration                (for UNSPLIT_GRAVITY only)
//                MinDens/Pres    : Density and pressure floors
//                EoS             : EoS object
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE
void Hydro_ComputeFlux( const real g_FC_Var [][NCOMP_TOTAL_PLUS_MAG][ CUBE(N_FC_VAR) ],
                              real g_FC_Flux[][NCOMP_TOTAL_PLUS_MAG][ CUBE(N_FC_FLUX) ],
                        const int NFlux, const int NSkip_N, const int NSkip_T,
                        const bool CorrHalfVel, const real g_Pot_USG[], const double g_Corner[],
                        const real dt, const real dh, const double Time, const bool UsePot,
                        const OptExtAcc_t ExtAcc, const ExtAcc_t ExtAcc_Func, const double ExtAcc_AuxArray[],
                        const real MinDens, const real MinPres, const EoS_t *EoS )
{

// check
#  ifdef GAMER_DEBUG
#  ifdef UNSPLIT_GRAVITY
   if ( CorrHalfVel )
   {
      if ( UsePot  &&  g_Pot_USG == NULL )
         printf( "ERROR : g_Pot_USG == NULL !!\n" );

      if ( ExtAcc  &&  g_Corner == NULL )
         printf( "ERROR : g_Corner == NULL !!\n" );
   }
#  else
   if ( CorrHalfVel )
      printf( "ERROR : CorrHalfVel is NOT supported when UNSPLIT_GRAVITY is off !!\n" );
#  endif

   if ( NFlux > N_FC_FLUX )
      printf( "ERROR : NFlux (%d) > N_FC_FLUX (%d) !!\n", NFlux, N_FC_FLUX );
#  endif // #ifdef GAMER_DEBUG


   const int didx_fc[3] = { 1, N_FC_VAR, N_FC_VAR*N_FC_VAR };

   real ConVar_L[NCOMP_TOTAL_PLUS_MAG], ConVar_R[NCOMP_TOTAL_PLUS_MAG], Flux_1Face[NCOMP_TOTAL_PLUS_MAG];

#  ifdef UNSPLIT_GRAVITY
   const real   GraConst    = -(real)0.5*dt/dh;
   const int    didx_usg[3] = { 1, USG_NXT_F, SQR(USG_NXT_F) };
   const int    fc_ghost    = ( N_FC_VAR - PS2 )/2;         // number of ghost zones on each side for g_FC_Var[]
   const int    idx_fc2usg  = USG_GHOST_SIZE_F - fc_ghost;  // index difference between g_FC_Var[] and g_Pot_USG[]
   const double dh_half     = 0.5*(double)dh;               // always use double precision to calculate the cell position
   const real   dt_half     = (real)0.5*dt;

   double CrShift[3] = { 0.0 };

// CrShift[]: central coordinates of the 0th cell in g_FC_Var[]
   if ( CorrHalfVel && ExtAcc )
      for (int d=0; d<3; d++)    CrShift[d] = g_Corner[d] - double(dh*fc_ghost);

// check
#  ifdef GAMER_DEBUG
   if ( CorrHalfVel )
   {
      if ( idx_fc2usg + NSkip_N < 0 )
         printf( "ERROR : idx_fc2usg (%d) + NSkip_N (%d) < 0 (USG_GHOST_SIZE_F %d, N_FC_VAR %d) !!\n",
                 idx_fc2usg, NSkip_N, USG_GHOST_SIZE_F, N_FC_VAR );

//    one additional cell is required to calculate the derivative along the transverse direction
      if ( idx_fc2usg + NSkip_T < 1 )
         printf( "ERROR : idx_fc2usg (%d) + NSkip_T (%d) < 1 (USG_GHOST_SIZE_F %d, N_FC_VAR %d) !!\n",
                 idx_fc2usg, NSkip_T, USG_GHOST_SIZE_F, N_FC_VAR );
   }
#  endif
#  endif // #ifdef UNSPLIT_GRAVITY


// loop over different spatial directions
   for (int d=0; d<3; d++)
   {
      const int faceL = 2*d;
      const int faceR = faceL+1;

#     ifdef UNSPLIT_GRAVITY
      int d1, d2, d3;

      if ( CorrHalfVel )
      {
         d1 =  d;
         d2 = (d+1)%3;
         d3 = (d+2)%3;
      }
#     endif

      int idx_fc_s[3], idx_flux_e[3];

      switch ( d )
      {
         case 0 : idx_fc_s  [0] = NSkip_N;              idx_fc_s  [1] = NSkip_T;              idx_fc_s  [2] = NSkip_T;
                  idx_flux_e[0] = N_FC_VAR-1-2*NSkip_N; idx_flux_e[1] = N_FC_VAR-2*NSkip_T;   idx_flux_e[2] = N_FC_VAR-2*NSkip_T;
                  break;

         case 1 : idx_fc_s  [0] = NSkip_T;              idx_fc_s  [1] = NSkip_N;              idx_fc_s  [2] = NSkip_T;
                  idx_flux_e[0] = N_FC_VAR-2*NSkip_T;   idx_flux_e[1] = N_FC_VAR-1-2*NSkip_N; idx_flux_e[2] = N_FC_VAR-2*NSkip_T;
                  break;

         case 2 : idx_fc_s  [0] = NSkip_T;              idx_fc_s  [1] = NSkip_T;              idx_fc_s  [2] = NSkip_N;
                  idx_flux_e[0] = N_FC_VAR-2*NSkip_T;   idx_flux_e[1] = N_FC_VAR-2*NSkip_T;   idx_flux_e[2] = N_FC_VAR-1-2*NSkip_N;
                  break;
      }

      const int size_ij = idx_flux_e[0]*idx_flux_e[1];
      CGPU_LOOP( idx, idx_flux_e[0]*idx_flux_e[1]*idx_flux_e[2] )
      {
         const int i_flux   = idx % idx_flux_e[0];
         const int j_flux   = idx % size_ij / idx_flux_e[0];
         const int k_flux   = idx / size_ij;
         const int idx_flux = IDX321( i_flux, j_flux, k_flux, NFlux, NFlux );

         const int i_fc     = i_flux + idx_fc_s[0];
         const int j_fc     = j_flux + idx_fc_s[1];
         const int k_fc     = k_flux + idx_fc_s[2];
         const int idx_fc   = IDX321( i_fc, j_fc, k_fc, N_FC_VAR, N_FC_VAR );

         for (int v=0; v<NCOMP_TOTAL_PLUS_MAG; v++)
         {
            ConVar_L[v] = g_FC_Var[faceR][v][ idx_fc            ];
            ConVar_R[v] = g_FC_Var[faceL][v][ idx_fc+didx_fc[d] ];
         }


//       1. correct the half-step velocity by gravity
#        ifdef UNSPLIT_GRAVITY
         if ( CorrHalfVel )
         {
            real Acc[3] = { (real)0.0, (real)0.0, (real)0.0 };
            real Enki_L, Enki_R;

//          external acceleration
            if ( ExtAcc )
            {
               double xyz[3]; // face-centered coordinates

               xyz[0]  = CrShift[0] + (double)(i_fc*dh);
               xyz[1]  = CrShift[1] + (double)(j_fc*dh);
               xyz[2]  = CrShift[2] + (double)(k_fc*dh);
               xyz[d] += dh_half;

               ExtAcc_Func( Acc, xyz[0], xyz[1], xyz[2], Time, ExtAcc_AuxArray );

               for (int t=0; t<3; t++)    Acc[t] *= dt_half;
            }

//          self-gravity and external potential
            if ( UsePot )
            {
               const int idx_usg = IDX321( i_fc+idx_fc2usg, j_fc+idx_fc2usg, k_fc+idx_fc2usg, USG_NXT_F, USG_NXT_F );

               Acc[d1] +=            GraConst*( g_Pot_USG[ idx_usg+didx_usg[d1] ] - g_Pot_USG[ idx_usg                           ] );
               Acc[d2] += (real)0.25*GraConst*( g_Pot_USG[ idx_usg+didx_usg[d2] ] + g_Pot_USG[ idx_usg+didx_usg[d2]+didx_usg[d1] ]
                                               -g_Pot_USG[ idx_usg-didx_usg[d2] ] - g_Pot_USG[ idx_usg-didx_usg[d2]+didx_usg[d1] ] );
               Acc[d3] += (real)0.25*GraConst*( g_Pot_USG[ idx_usg+didx_usg[d3] ] + g_Pot_USG[ idx_usg+didx_usg[d3]+didx_usg[d1] ]
                                               -g_Pot_USG[ idx_usg-didx_usg[d3] ] - g_Pot_USG[ idx_usg-didx_usg[d3]+didx_usg[d1] ] );
            }

//          store the "non"-kinetic energy (i.e. total energy - kinetic energy)
            Enki_L = ConVar_L[4] - (real)0.5*( SQR(ConVar_L[1]) + SQR(ConVar_L[2]) + SQR(ConVar_L[3]) )/ConVar_L[0];
            Enki_R = ConVar_R[4] - (real)0.5*( SQR(ConVar_R[1]) + SQR(ConVar_R[2]) + SQR(ConVar_R[3]) )/ConVar_R[0];

//          advance velocity by gravity
            for (int t=0; t<3; t++)
            {
               ConVar_L[t+1] += ConVar_L[0]*Acc[t];
               ConVar_R[t+1] += ConVar_R[0]*Acc[t];
            }

//          update total energy density with the non-kinetic energy fixed
            ConVar_L[4] = Enki_L + (real)0.5*( SQR(ConVar_L[1]) + SQR(ConVar_L[2]) + SQR(ConVar_L[3]) )/ConVar_L[0];
            ConVar_R[4] = Enki_R + (real)0.5*( SQR(ConVar_R[1]) + SQR(ConVar_R[2]) + SQR(ConVar_R[3]) )/ConVar_R[0];
         } // if ( CorrHalfVel )
#        endif // #ifdef UNSPLIT_GRAVITY


//       2. invoke Riemann solver
#        if   ( RSOLVER == EXACT  &&  !defined MHD )
         Hydro_RiemannSolver_Exact( d, Flux_1Face, ConVar_L, ConVar_R, MinDens, MinPres,
                                    EoS->DensEint2Pres_FuncPtr, EoS->DensPres2CSqr_FuncPtr,
                                    EoS->AuxArrayDevPtr_Flt, EoS->AuxArrayDevPtr_Int, EoS->Table );
#        elif ( RSOLVER == ROE )
         Hydro_RiemannSolver_Roe  ( d, Flux_1Face, ConVar_L, ConVar_R, MinDens, MinPres,
                                    EoS->DensEint2Pres_FuncPtr, EoS->DensPres2CSqr_FuncPtr,
                                    EoS->AuxArrayDevPtr_Flt, EoS->AuxArrayDevPtr_Int, EoS->Table );
#        elif ( RSOLVER == HLLE )
         Hydro_RiemannSolver_HLLE ( d, Flux_1Face, ConVar_L, ConVar_R, MinDens, MinPres,
                                    EoS->DensEint2Pres_FuncPtr, EoS->DensPres2CSqr_FuncPtr,
                                    EoS->GuessHTilde_FuncPtr, EoS->HTilde2Temp_FuncPtr,
                                    EoS->AuxArrayDevPtr_Flt, EoS->AuxArrayDevPtr_Int, EoS->Table );
#        elif ( RSOLVER == HLLC  &&  !defined MHD )
         Hydro_RiemannSolver_HLLC ( d, Flux_1Face, ConVar_L, ConVar_R, MinDens, MinPres,
                                    EoS->DensEint2Pres_FuncPtr, EoS->DensPres2CSqr_FuncPtr,
                                    EoS->GuessHTilde_FuncPtr, EoS->HTilde2Temp_FuncPtr,
                                    EoS->AuxArrayDevPtr_Flt, EoS->AuxArrayDevPtr_Int, EoS->Table );
#        elif ( RSOLVER == HLLD  &&  defined MHD )
         Hydro_RiemannSolver_HLLD ( d, Flux_1Face, ConVar_L, ConVar_R, MinDens, MinPres,
                                    EoS->DensEint2Pres_FuncPtr, EoS->DensPres2CSqr_FuncPtr,
                                    EoS->AuxArrayDevPtr_Flt, EoS->AuxArrayDevPtr_Int, EoS->Table );
#        else
#        error : ERROR : unsupported Riemann solver (EXACT/ROE/HLLE/HLLC/HLLD) !!
#        endif


//       3. switch to a different Riemann solver if the default one fails
#        if ( RSOLVER_RESCUE != NONE )
         for (int v=0; v<NCOMP_TOTAL_PLUS_MAG; v++)
         {
//          only check NaN for now
            if ( Flux_1Face[v] != Flux_1Face[v] )
            {
#              ifdef CHECK_UNPHYSICAL_IN_FLUID
               printf( "WARNING : default Riemann solver failed in Hydro_ComputeFlux() --> switch to RSOLVER_RESCUE (%d) !!\n", RSOLVER_RESCUE );
#              endif

#              if   ( RSOLVER_RESCUE == EXACT  &&  !defined MHD )
               Hydro_RiemannSolver_Exact( d, Flux_1Face, ConVar_L, ConVar_R, MinDens, MinPres,
                                          EoS->DensEint2Pres_FuncPtr, EoS->DensPres2CSqr_FuncPtr,
                                          EoS->AuxArrayDevPtr_Flt, EoS->AuxArrayDevPtr_Int, EoS->Table );
#              elif ( RSOLVER_RESCUE == ROE )
               Hydro_RiemannSolver_Roe  ( d, Flux_1Face, ConVar_L, ConVar_R, MinDens, MinPres,
                                          EoS->DensEint2Pres_FuncPtr, EoS->DensPres2CSqr_FuncPtr,
                                          EoS->AuxArrayDevPtr_Flt, EoS->AuxArrayDevPtr_Int, EoS->Table );
#              elif ( RSOLVER_RESCUE == HLLE )
               Hydro_RiemannSolver_HLLE ( d, Flux_1Face, ConVar_L, ConVar_R, MinDens, MinPres,
                                          EoS->DensEint2Pres_FuncPtr, EoS->DensPres2CSqr_FuncPtr,
                                          EoS->GuessHTilde_FuncPtr, EoS->HTilde2Temp_FuncPtr,
                                          EoS->AuxArrayDevPtr_Flt, EoS->AuxArrayDevPtr_Int, EoS->Table );
#              elif ( RSOLVER_RESCUE == HLLC  &&  !defined MHD )
               Hydro_RiemannSolver_HLLC ( d, Flux_1Face, ConVar_L, ConVar_R, MinDens, MinPres,
                                          EoS->DensEint2Pres_FuncPtr, EoS->DensPres2CSqr_FuncPtr,
                                          EoS->GuessHTilde_FuncPtr, EoS->HTilde2Temp_FuncPtr,
                                          EoS->AuxArrayDevPtr_Flt, EoS->AuxArrayDevPtr_Int, EoS->Table );
#              elif ( RSOLVER_RESCUE == HLLD  &&  defined MHD )
               Hydro_RiemannSolver_HLLD ( d, Flux_1Face, ConVar_L, ConVar_R, MinDens, MinPres,
                                          EoS->DensEint2Pres_FuncPtr, EoS->DensPres2CSqr_FuncPtr,
                                          EoS->AuxArrayDevPtr_Flt, EoS->AuxArrayDevPtr_Int, EoS->Table );
#              else
#              error : ERROR : unsupported Riemann solver (EXACT/ROE/HLLE/HLLC/HLLD) !!
#              endif

//             check again
#              ifdef CHECK_UNPHYSICAL_IN_FLUID
               for (int w=0; w<NCOMP_TOTAL_PLUS_MAG; w++) {
                  if ( Flux_1Face[w] != Flux_1Face[w] ) {
                     printf( "ERROR : RSOLVER_RESCUE still failed !!\n" );
                     break;
                  }
               }
#              endif

               break;
            } // if ( Flux_1Face[v] != Flux_1Face[v] )
         } // for (int v=0; v<NCOMP_TOTAL_PLUS_MAG; v++)
#        endif // #if ( RSOLVER_RESCUE != NONE )


//       4. store the fluxes of all cells in g_FC_Flux[]
//       --> including the magnetic components since they are required for CT
         for (int v=0; v<NCOMP_TOTAL_PLUS_MAG; v++)   g_FC_Flux[d][v][idx_flux] = Flux_1Face[v];
      } // i,j,k
   } // for (int d=0; d<3; d++)


#  ifdef __CUDACC__
   __syncthreads();
#  endif

} // FUNCTION : Hydro_ComputeFlux



//-------------------------------------------------------------------------------------------------------
// Function    :  Hydro_StoreIntFlux
// Description :  Store the inter-patch fluxes in g_IntFlux[]
//
// Note        :  1. No need to store the magnetic components since this array is only for the fluid flux fix-up operation
//
// Parameter   :  g_FC_Flux : Array storing the face-centered fluxes
//                g_IntFlux : Array to store the face-center fluxes at the patch interfaces
//                NFlux     : Stride for accessing g_FC_Flux[]
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE
void Hydro_StoreIntFlux( const real g_FC_Flux[][NCOMP_TOTAL_PLUS_MAG][ CUBE(N_FC_FLUX) ],
                               real g_IntFlux[][NCOMP_TOTAL][ SQR(PS2) ],
                         const int NFlux )
{

// check
#  ifdef GAMER_DEBUG
   if ( NFlux > N_FC_FLUX )
      printf( "ERROR : NFlux (%d) > N_FC_FLUX (%d) !!\n", NFlux, N_FC_FLUX );
#  endif


   CGPU_LOOP( idx_out, SQR(PS2) )
   {
//    indices of the 2 transverse directions
//    --> for MHD, one additional flux is evaluated along each transverse direction for computing the CT electric field
#     ifdef MHD
      const int idx_flux_0  = idx_out % PS2 + 1;
      const int idx_flux_1  = idx_out / PS2 + 1;
#     else
      const int idx_flux_0  = idx_out % PS2;
      const int idx_flux_1  = idx_out / PS2;
#     endif

//    flux index
//    x-face
      const int idx_in_0  = IDX321(          0, idx_flux_0, idx_flux_1, NFlux, NFlux );
      const int idx_in_1  = IDX321(        PS1, idx_flux_0, idx_flux_1, NFlux, NFlux );
      const int idx_in_2  = IDX321(        PS2, idx_flux_0, idx_flux_1, NFlux, NFlux );
//    y-face
      const int idx_in_3  = IDX321( idx_flux_0,          0, idx_flux_1, NFlux, NFlux );
      const int idx_in_4  = IDX321( idx_flux_0,        PS1, idx_flux_1, NFlux, NFlux );
      const int idx_in_5  = IDX321( idx_flux_0,        PS2, idx_flux_1, NFlux, NFlux );
//    z-face
      const int idx_in_6  = IDX321( idx_flux_0, idx_flux_1,          0, NFlux, NFlux );
      const int idx_in_7  = IDX321( idx_flux_0, idx_flux_1,        PS1, NFlux, NFlux );
      const int idx_in_8  = IDX321( idx_flux_0, idx_flux_1,        PS2, NFlux, NFlux );

//    store the flux
      for (int v=0; v<NCOMP_TOTAL; v++) {
//       x-face
         g_IntFlux[0][v][idx_out] = g_FC_Flux[0][v][idx_in_0];
         g_IntFlux[1][v][idx_out] = g_FC_Flux[0][v][idx_in_1];
         g_IntFlux[2][v][idx_out] = g_FC_Flux[0][v][idx_in_2];
//       y-face
         g_IntFlux[3][v][idx_out] = g_FC_Flux[1][v][idx_in_3];
         g_IntFlux[4][v][idx_out] = g_FC_Flux[1][v][idx_in_4];
         g_IntFlux[5][v][idx_out] = g_FC_Flux[1][v][idx_in_5];
//       z-face
         g_IntFlux[6][v][idx_out] = g_FC_Flux[2][v][idx_in_6];
         g_IntFlux[7][v][idx_out] = g_FC_Flux[2][v][idx_in_7];
         g_IntFlux[8][v][idx_out] = g_FC_Flux[2][v][idx_in_8];
      }
   } // CGPU_LOOP( idx_out, SQR(PS2) )

#  ifdef __CUDACC__
   __syncthreads();
#  endif

} // FUNCTION : Hydro_StoreIntFlux



#endif // #if ( MODEL == HYDRO  &&  (FLU_SCHEME == MHM || FLU_SCHEME == MHM_RP || FLU_SCHEME == CTU) )



#endif // #ifndef __CUFLU_COMPUTEFLUX__
