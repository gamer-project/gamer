#ifndef __CUFLU_COMPUTEFLUX__
#define __CUFLU_COMPUTEFLUX__



#include "CUFLU.h"

#if ( MODEL == HYDRO  &&  (FLU_SCHEME == MHM || FLU_SCHEME == MHM_RP || FLU_SCHEME == CTU) )



// external functions
#ifdef __CUDACC__

#if   ( RSOLVER == EXACT )
# include "CUFLU_Shared_RiemannSolver_Exact.cu"
#elif ( RSOLVER == ROE )
# include "CUFLU_Shared_RiemannSolver_Roe.cu"
#elif ( RSOLVER == HLLE )
# include "CUFLU_Shared_RiemannSolver_HLLE.cu"
#elif ( RSOLVER == HLLC )
# include "CUFLU_Shared_RiemannSolver_HLLC.cu"
#elif ( RSOLVER == HLLD )
# include "CUFLU_Shared_RiemannSolver_HLLD.cu"
#endif

#ifdef UNSPLIT_GRAVITY
# include "../../SelfGravity/GPU_Gravity/CUPOT_ExternalAcc.cu"
#endif

#else // #ifdef __CUDACC__

#if   ( RSOLVER == EXACT )
void Hydro_Con2Pri( const real In[], real Out[], const real Gamma_m1, const real MinPres,
                    const bool NormPassive, const int NNorm, const int NormIdx[],
                    const bool JeansMinPres, const real JeansMinPres_Coeff );
void Hydro_RiemannSolver_Exact( const int XYZ, real Flux_Out[], const real L_In[], const real R_In[], const real Gamma );
#elif ( RSOLVER == ROE )
void Hydro_RiemannSolver_Roe( const int XYZ, real Flux_Out[], const real L_In[], const real R_In[],
                              const real Gamma, const real MinPres );
#elif ( RSOLVER == HLLE )
void Hydro_RiemannSolver_HLLE( const int XYZ, real Flux_Out[], const real L_In[], const real R_In[],
                               const real Gamma, const real MinPres );
#elif ( RSOLVER == HLLC )
void Hydro_RiemannSolver_HLLC( const int XYZ, real Flux_Out[], const real L_In[], const real R_In[],
                               const real Gamma, const real MinPres );
#elif ( RSOLVER == HLLD )
void Hydro_RiemannSolver_HLLD( const int XYZ, real Flux_Out[], const real L_In[], const real R_In[],
                               const real Gamma, const real MinPres );
#endif
#ifdef UNSPLIT_GRAVITY
void ExternalAcc( real Acc[], const double x, const double y, const double z, const double Time, const double UserArray[] );
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
//                5. For the performance consideration, this function will also be responsible for storing the
//                   inter-patch fluxes
//                   --> Option "DumpIntFlux"
//                6. For the unsplitting scheme in gravity (i.e., UNSPLIT_GRAVITY), this function also corrects the half-step
//                   velocity by gravity when CorrHalfVel==true
//
// Parameter   :  g_FC_Var        : Array storing the input face-centered conserved variables
//                g_FC_Flux       : Array to store the output face-centered fluxes
//                NFlux           : Stride for accessing g_FC_Flux[]
//                NSkip_N         : Number of cells to be skipped in the normal directions
//                                  --> "(N_FC_VAR-1-2*NSkip_N)" fluxes will be computed along the normal direction
//                NSkip_T         : Number of cells to be skipped in the transverse directions
//                                  --> "(N_FC_VAR-2*NSkip_T)^2" fluxes will be computed along the transverse direction
//                Gamma           : Ratio of specific heats
//                CorrHalfVel     : true --> correct the half-step velocity by gravity        (for UNSPLIT_GRAVITY only)
//                g_Pot_USG       : Array storing the input potential for CorrHalfVel         (for UNSPLIT_GRAVITY only)
//                g_Corner        : Array storing the corner coordinates of each patch group  (for UNSPLIT_GRAVITY only)
//                dt              : Time interval to advance the full-step solution           (for UNSPLIT_GRAVITY only)
//                dh              : Cell size                                                 (for UNSPLIT_GRAVITY only)
//                Time            : Current physical time                                     (for UNSPLIT_GRAVITY only)
//                GravityType     : Types of gravity --> self-gravity, external gravity, both (for UNSPLIT_GRAVITY only)
//                ExtAcc_AuxArray : Auxiliary array for external acceleration                 (for UNSPLIT_GRAVITY only)
//                MinPres         : Minimum allowed pressure
//                DumpIntFlux     : true --> store the inter-patch fluxes in g_IntFlux[]
//                g_IntFlux       : Array for DumpIntFlux
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE
void Hydro_ComputeFlux( const real g_FC_Var [][NCOMP_TOTAL_PLUS_MAG][ CUBE(N_FC_VAR) ],
                              real g_FC_Flux[][NCOMP_TOTAL_PLUS_MAG][ CUBE(N_FC_FLUX) ],
                        const int NFlux, const int NSkip_N, const int NSkip_T, const real Gamma,
                        const bool CorrHalfVel, const real g_Pot_USG[], const double g_Corner[],
                        const real dt, const real dh, const double Time,
                        const OptGravityType_t GravityType, const double ExtAcc_AuxArray[],
                        const real MinPres, const bool DumpIntFlux, real g_IntFlux[][NCOMP_TOTAL][ SQR(PS2) ] )
{

// check
#  ifdef GAMER_DEBUG
#  ifdef UNSPLIT_GRAVITY
   if ( CorrHalfVel )
   {
      if (  ( GravityType == GRAVITY_SELF || GravityType == GRAVITY_BOTH )  &&  g_Pot_USG == NULL  )
         printf( "ERROR : g_Pot_USG == NULL !!\n" );

      if (  ( GravityType == GRAVITY_EXTERNAL || GravityType == GRAVITY_BOTH )  &&  g_Corner == NULL  )
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

#  if ( RSOLVER == EXACT )
   const real Gamma_m1 = Gamma - (real)1.0;
   real PriVar_L[NCOMP_TOTAL], PriVar_R[NCOMP_TOTAL];
#  endif

#  ifdef UNSPLIT_GRAVITY
   const real   GraConst    = -(real)0.5*dt/dh;
   const int    didx_usg[3] = { 1, USG_NXT_F, SQR(USG_NXT_F) };
   const int    fc_ghost    = ( N_FC_VAR - PS2 )/2;         // number of ghost zones on each side for g_FC_Var[]
   const int    idx_fc2usg  = USG_GHOST_SIZE_F - fc_ghost;  // index difference between g_FC_Var[] and g_Pot_USG[]
   const double dh_half     = 0.5*(double)dh;               // always use double precision to calculate the cell position
   const real   dt_half     = (real)0.5*dt;

   double CrShift[3];

// CrShift[]: central coordinates of the 0th cell in g_FC_Var[]
   if (  CorrHalfVel  &&  ( GravityType == GRAVITY_EXTERNAL || GravityType == GRAVITY_BOTH )  )
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
            real   Acc[3], eL, eR;
            double xyz[3];

            Acc[0] = (real)0.0;
            Acc[1] = (real)0.0;
            Acc[2] = (real)0.0;

//          external gravity
            if ( GravityType == GRAVITY_EXTERNAL  ||  GravityType == GRAVITY_BOTH )
            {
//             xyz[]: face-centered coordinates
               xyz[0]  = CrShift[0] + (double)(i_fc*dh);
               xyz[1]  = CrShift[1] + (double)(j_fc*dh);
               xyz[2]  = CrShift[2] + (double)(k_fc*dh);
               xyz[d] += dh_half;

               ExternalAcc( Acc, xyz[0], xyz[1], xyz[2], Time, ExtAcc_AuxArray );

               for (int t=0; t<3; t++)    Acc[t] *= dt_half;
            }

//          self-gravity
            if ( GravityType == GRAVITY_SELF  ||  GravityType == GRAVITY_BOTH )
            {
               const int idx_usg = IDX321( i_fc+idx_fc2usg, j_fc+idx_fc2usg, k_fc+idx_fc2usg, USG_NXT_F, USG_NXT_F );

               Acc[d1] +=            GraConst*( g_Pot_USG[ idx_usg+didx_usg[d1] ] - g_Pot_USG[ idx_usg                           ] );
               Acc[d2] += (real)0.25*GraConst*( g_Pot_USG[ idx_usg+didx_usg[d2] ] + g_Pot_USG[ idx_usg+didx_usg[d2]+didx_usg[d1] ]
                                               -g_Pot_USG[ idx_usg-didx_usg[d2] ] - g_Pot_USG[ idx_usg-didx_usg[d2]+didx_usg[d1] ] );
               Acc[d3] += (real)0.25*GraConst*( g_Pot_USG[ idx_usg+didx_usg[d3] ] + g_Pot_USG[ idx_usg+didx_usg[d3]+didx_usg[d1] ]
                                               -g_Pot_USG[ idx_usg-didx_usg[d3] ] - g_Pot_USG[ idx_usg-didx_usg[d3]+didx_usg[d1] ] );
            }

//          store the internal energy density (plus the magnetic energy in MHD)
            eL = ConVar_L[4] - (real)0.5*( SQR(ConVar_L[1]) + SQR(ConVar_L[2]) + SQR(ConVar_L[3]) )/ConVar_L[0];
            eR = ConVar_R[4] - (real)0.5*( SQR(ConVar_R[1]) + SQR(ConVar_R[2]) + SQR(ConVar_R[3]) )/ConVar_R[0];

//          advance velocity by gravity
            for (int t=0; t<3; t++)
            {
               ConVar_L[t+1] += ConVar_L[0]*Acc[t];
               ConVar_R[t+1] += ConVar_R[0]*Acc[t];
            }

//          update total energy density with the internal energy density (plus the magnetic energy in MHD) fixed
            ConVar_L[4] = eL + (real)0.5*( SQR(ConVar_L[1]) + SQR(ConVar_L[2]) + SQR(ConVar_L[3]) )/ConVar_L[0];
            ConVar_R[4] = eR + (real)0.5*( SQR(ConVar_R[1]) + SQR(ConVar_R[2]) + SQR(ConVar_R[3]) )/ConVar_R[0];
         } // if ( CorrHalfVel )
#        endif // #ifdef UNSPLIT_GRAVITY


//       2. invoke Riemann solver
#        if   ( RSOLVER == EXACT  &&  !defined MHD )
         const bool NormPassive_No  = false; // do NOT convert any passive variable to mass fraction for the Riemann solvers
         const bool JeansMinPres_No = false;

         Hydro_Con2Pri( ConVar_L, PriVar_L, Gamma_m1, MinPres, NormPassive_No, NULL_INT, NULL, JeansMinPres_No, NULL_REAL );
         Hydro_Con2Pri( ConVar_R, PriVar_R, Gamma_m1, MinPres, NormPassive_No, NULL_INT, NULL, JeansMinPres_No, NULL_REAL );

         Hydro_RiemannSolver_Exact( d, Flux_1Face, PriVar_L, PriVar_R, Gamma );
#        elif ( RSOLVER == ROE )
         Hydro_RiemannSolver_Roe  ( d, Flux_1Face, ConVar_L, ConVar_R, Gamma, MinPres );
#        elif ( RSOLVER == HLLE )
         Hydro_RiemannSolver_HLLE ( d, Flux_1Face, ConVar_L, ConVar_R, Gamma, MinPres );
#        elif ( RSOLVER == HLLC  &&  !defined MHD )
         Hydro_RiemannSolver_HLLC ( d, Flux_1Face, ConVar_L, ConVar_R, Gamma, MinPres );
#        elif ( RSOLVER == HLLD  &&  defined MHD )
         Hydro_RiemannSolver_HLLD ( d, Flux_1Face, ConVar_L, ConVar_R, Gamma, MinPres );
#        else
#        error : ERROR : unsupported Riemann solver (EXACT/ROE/HLLE/HLLC/HLLD) !!
#        endif


//       3. store the fluxes of all cells in g_FC_Flux[]
//       --> including the magnetic components since they are required for CT
         for (int v=0; v<NCOMP_TOTAL_PLUS_MAG; v++)   g_FC_Flux[d][v][idx_flux] = Flux_1Face[v];


//       4. store the inter-patch fluxes in g_IntFlux[]
//       --> no need to store the magnetic components since this array is only for the flux fix-up operation
         if ( DumpIntFlux )
         {
            int int_face, int_idx;

//          we have assumed N_FC_VAR=PS2+2 for pure hydro
//          --> for MHD, one additional flux is evaluated along each transverse direction for computing the CT electric field
//          --> must exclude it when storing the inter-patch fluxes
            if (  d == 0  &&  ( i_flux == 0 || i_flux == PS1 || i_flux == PS2 )  )
            {
#              ifdef MHD
               if ( j_flux > 0  &&  j_flux < PS2+1  &&  k_flux > 0  &&  k_flux < PS2+1 )
#              endif
               {
                  int_face = i_flux/PS1;
#                 ifdef MHD
                  int_idx  = (k_flux-1)*PS2 + j_flux-1;
#                 else
                  int_idx  = (k_flux  )*PS2 + j_flux;
#                 endif
                  for (int v=0; v<NCOMP_TOTAL; v++)   g_IntFlux[int_face][v][int_idx] = Flux_1Face[v];
               }
            }

            else if (  d == 1  &&  ( j_flux == 0 || j_flux == PS1 || j_flux == PS2 )  )
            {
#              ifdef MHD
               if ( i_flux > 0  &&  i_flux < PS2+1  &&  k_flux > 0  &&  k_flux < PS2+1 )
#              endif
               {
                  int_face = j_flux/PS1 + 3;
#                 ifdef MHD
                  int_idx  = (k_flux-1)*PS2 + i_flux-1;
#                 else
                  int_idx  = (k_flux  )*PS2 + i_flux;
#                 endif
                  for (int v=0; v<NCOMP_TOTAL; v++)   g_IntFlux[int_face][v][int_idx] = Flux_1Face[v];
               }
            }

            else if (  d == 2  &&  ( k_flux == 0 || k_flux == PS1 || k_flux == PS2 )  )
            {
#              ifdef MHD
               if ( i_flux > 0  &&  i_flux < PS2+1  &&  j_flux > 0  &&  j_flux < PS2+1 )
#              endif
               {
                  int_face = k_flux/PS1 + 6;
#                 ifdef MHD
                  int_idx  = (j_flux-1)*PS2 + i_flux-1;
#                 else
                  int_idx  = (j_flux  )*PS2 + i_flux;
#                 endif
                  for (int v=0; v<NCOMP_TOTAL; v++)   g_IntFlux[int_face][v][int_idx] = Flux_1Face[v];
               }
            }
         } // if ( DumpIntFlux )
      } // i,j,k
   } // for (int d=0; d<3; d++)


#  ifdef __CUDACC__
   __syncthreads();
#  endif

} // FUNCTION : Hydro_ComputeFlux



#endif // #if ( MODEL == HYDRO  &&  (FLU_SCHEME == MHM || FLU_SCHEME == MHM_RP || FLU_SCHEME == CTU) )



#endif // #ifndef __CUFLU_COMPUTEFLUX__
