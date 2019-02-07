#ifndef __CUFLU_COMPUTEFLUX__
#define __CUFLU_COMPUTEFLUX__



#include "CUFLU.h"

#if ( MODEL == SR_HYDRO  &&  (FLU_SCHEME == MHM || FLU_SCHEME == MHM_RP) )



// external functions
#ifdef __CUDACC__

#if ( RSOLVER == HLLE )
# include "CUFLU_Shared_RiemannSolver_HLLE.cu"
#elif ( RSOLVER == HLLC )
# include "CUFLU_Shared_RiemannSolver_HLLC.cu"
#endif

#ifdef UNSPLIT_GRAVITY
# include "../../SelfGravity/GPU_Gravity/CUPOT_ExternalAcc.cu"
#endif

#else // #ifdef __CUDACC__

#if ( RSOLVER == HLLE )
void SRHydro_RiemannSolver_HLLE( const int XYZ, real Flux_Out[], const real L_In[], const real R_In[],
                               const real Gamma, const real MinTemp );
#elif ( RSOLVER == HLLC )
void SRHydro_RiemannSolver_HLLC( const int XYZ, real Flux_Out[], const real L_In[], const real R_In[],
                               const real Gamma, const real MinTemp );
#endif

#ifdef UNSPLIT_GRAVITY
void ExternalAcc( real Acc[], const double x, const double y, const double z, const double Time, const double UserArray[] );
#endif

#endif // #ifdef __CUDACC__ ... else ...




//-------------------------------------------------------------------------------------------------------
// Function    :  SRHydro_ComputeFlux
// Description :  Compute the face-centered fluxes by Riemann solver
//
// Note        :  1. Currently support the exact, HLLC, HLLE, and Roe solvers
//                2. g_FC_Var[] has the size of N_FC_VAR^3
//                   --> "N_FC_VAR-1" fluxes will be computed along each normal direction
//                   --> But "2*Gap" cells will skipped along the transverse direction (Gap cells on each side)
//                3. g_FC_Flux[] has the size of N_FC_FLUX^3
//                   --> But (i,j,k) flux will be stored in the "(k*N_FL_FLUX+j)*N_FL_FLUX+i" element in g_FC_Flux[]
//                       --> We have assumed that N_FL_FLUX <= N_FC_FLUX
//                   --> (i,j,k) in g_FC_Flux_x is defined on the +x surface of the cell (i,     j+Gap, k+Gap) in g_FC_Var[]
//                       (i,j,k) in g_FC_Flux_y is defined on the +y surface of the cell (i+Gap, j,     k+Gap) in g_FC_Var[]
//                       (i,j,k) in g_FC_Flux_z is defined on the +z surface of the cell (i+Gap, j+Gap, k    ) in g_FC_Var[]
//                4. This function is shared by MHM, MHM_RP, and CTU schemes
//                5. For the performance consideration, this function will also be responsible for storing the
//                   inter-patch fluxes
//                   --> Option "DumpIntFlux"
//
// Parameter   :  g_FC_Var        : Array storing the input face-centered conserved variables
//                g_FC_Flux       : Array to store the output face-centered fluxes
//                Gap             : Number of cells to be skipped in the transverse directions
//                                  --> "(N_FC_VAR-2*Gap)^2" fluxes will be computed on each surface
//                Gamma           : Ratio of specific heats
//                MinTemp         : Minimum allowed temerature
//                DumpIntFlux     : true --> store the inter-patch fluxes in g_IntFlux[]
//                g_IntFlux       : Array for DumpIntFlux
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE
void SRHydro_ComputeFlux( const real g_FC_Var [][NCOMP_TOTAL][ CUBE(N_FC_VAR) ],
                                real g_FC_Flux[][NCOMP_TOTAL][ CUBE(N_FC_FLUX) ],
                          const int Gap, const real Gamma, const real MinTemp, 
                          const bool DumpIntFlux, real g_IntFlux[][NCOMP_TOTAL][ SQR(PS2) ] )
{

// check
#  ifdef GAMER_DEBUG
   if ( CorrHalfVel )
      printf( "ERROR : CorrHalfVel is NOT supported when UNSPLIT_GRAVITY is off !!\n" );
#  endif // #ifdef GAMER_DEBUG

   const int didx_fc[3] = { 1, N_FC_VAR, N_FC_VAR*N_FC_VAR };

   real ConVar_L[NCOMP_TOTAL], ConVar_R[NCOMP_TOTAL], Flux_1Face[NCOMP_TOTAL];

// loop over different spatial directions
   for (int d=0; d<3; d++)
   {
      const int faceL = 2*d;
      const int faceR = faceL+1;

      int idx_fc_s[3], idx_flux_e[3];

      switch ( d )
      {
         case 0 : idx_fc_s  [0] = 0;              idx_fc_s  [1] = Gap;            idx_fc_s  [2] = Gap;
                  idx_flux_e[0] = N_FC_VAR-1;     idx_flux_e[1] = N_FC_VAR-2*Gap; idx_flux_e[2] = N_FC_VAR-2*Gap; break;

         case 1 : idx_fc_s  [0] = Gap;            idx_fc_s  [1] = 0;              idx_fc_s  [2] = Gap;
                  idx_flux_e[0] = N_FC_VAR-2*Gap; idx_flux_e[1] = N_FC_VAR-1;     idx_flux_e[2] = N_FC_VAR-2*Gap; break;

         case 2 : idx_fc_s  [0] = Gap;            idx_fc_s  [1] = Gap;            idx_fc_s  [2] = 0;
                  idx_flux_e[0] = N_FC_VAR-2*Gap; idx_flux_e[1] = N_FC_VAR-2*Gap; idx_flux_e[2] = N_FC_VAR-1;     break;
      }

      const int size_ij = idx_flux_e[0]*idx_flux_e[1];
      CGPU_LOOP( idx, idx_flux_e[0]*idx_flux_e[1]*idx_flux_e[2] )
      {
         const int i_flux   = idx % idx_flux_e[0];
         const int j_flux   = idx % size_ij / idx_flux_e[0];
         const int k_flux   = idx / size_ij;
         const int idx_flux = IDX321( i_flux, j_flux, k_flux, N_FL_FLUX, N_FL_FLUX );

         const int i_fc     = i_flux + idx_fc_s[0];
         const int j_fc     = j_flux + idx_fc_s[1];
         const int k_fc     = k_flux + idx_fc_s[2];
         const int idx_fc   = IDX321( i_fc, j_fc, k_fc, N_FC_VAR, N_FC_VAR );

         for (int v=0; v<NCOMP_TOTAL; v++)
         {
            ConVar_L[v] = g_FC_Var[faceR][v][ idx_fc            ];
            ConVar_R[v] = g_FC_Var[faceL][v][ idx_fc+didx_fc[d] ];
         }

#        ifdef CHECK_NEGATIVE_IN_FLUID
         if( SRHydro_CheckUnphysical( ConVar_L, NULL, Gamma, MinTemp, __FUNCTION__, __LINE__, true )
          || SRHydro_CheckUnphysical( ConVar_R, NULL, Gamma, MinTemp, __FUNCTION__, __LINE__, true ) )
         exist(EXIT_FAILURE);
#        endif

//       2. invoke Riemann solver
#        if ( RSOLVER == HLLE )
         SRHydro_RiemannSolver_HLLE ( d, Flux_1Face, ConVar_L, ConVar_R, Gamma, MinTemp );
#        elif ( RSOLVER == HLLC )
         SRHydro_RiemannSolver_HLLC ( d, Flux_1Face, ConVar_L, ConVar_R, Gamma, MinTemp );
#        else
#        error : ERROR : unsupported Riemann solver !!
#        endif


//       3. store results
//       store the fluxes of all cells in g_FC_Flux[]
         for (int v=0; v<NCOMP_TOTAL; v++)   g_FC_Flux[d][v][idx_flux] = Flux_1Face[v];

//       store the inter-patch fluxes in g_IntFlux[]
         if ( DumpIntFlux )
         {
            int int_face, int_idx;

//          we have assumed N_FC_VAR=PS2+2
            if      (  d == 0  &&  ( i_flux == 0 || i_flux == PS1 || i_flux == PS2 )  )
            {
               int_face = i_flux/PS1;
               int_idx  = k_flux*PS2 + j_flux;
               for (int v=0; v<NCOMP_TOTAL; v++)   g_IntFlux[int_face][v][int_idx] = Flux_1Face[v];
            }

            else if (  d == 1  &&  ( j_flux == 0 || j_flux == PS1 || j_flux == PS2 )  )
            {
               int_face = j_flux/PS1 + 3;
               int_idx  = k_flux*PS2 + i_flux;
               for (int v=0; v<NCOMP_TOTAL; v++)   g_IntFlux[int_face][v][int_idx] = Flux_1Face[v];
            }

            else if (  d == 2  &&  ( k_flux == 0 || k_flux == PS1 || k_flux == PS2 )  )
            {
               int_face = k_flux/PS1 + 6;
               int_idx  = j_flux*PS2 + i_flux;
               for (int v=0; v<NCOMP_TOTAL; v++)   g_IntFlux[int_face][v][int_idx] = Flux_1Face[v];
            }
         } // if ( DumpIntFlux )
      } // i,j,k
   } // for (int d=0; d<3; d++)


#  ifdef __CUDACC__
   __syncthreads();
#  endif

} // FUNCTION : SRHydro_ComputeFlux



#endif // #if ( MODEL == SR_HYDRO  &&  (FLU_SCHEME == MHM || FLU_SCHEME == MHM_RP || FLU_SCHEME == CTU) )



#endif // #ifndef __CUFLU_COMPUTEFLUX__
