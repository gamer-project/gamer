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
// Parameter   :  [ 1] g_FC_Var        : Array storing the input face-centered conserved variables
//                [ 2] g_FC_Flux       : Array to store the output face-centered fluxes
//                [ 3] Gap             : Number of cells to be skipped in the transverse directions
//                                       --> "(N_FC_VAR-2*Gap)^2" fluxes will be computed on each surface
//                [ 4] Gamma           : Ratio of specific heats
//                [ 5] CorrHalfVel     : true --> correct the half-step velocity by gravity        (for UNSPLIT_GRAVITY only)
//                [ 6] g_Pot_USG       : Array storing the input potential for CorrHalfVel         (for UNSPLIT_GRAVITY only)
//                                       --> must have the same size as g_FC_Var[] --> (PS2+2)^3
//                [ 7] g_Corner        : Array storing the corner coordinates of each patch group  (for UNSPLIT_GRAVITY only)
//                [ 8] dt              : Time interval to advance the full-step solution           (for UNSPLIT_GRAVITY only)
//                [ 9] dh              : Cell size                                                 (for UNSPLIT_GRAVITY only)
//                [10] Time            : Current physical time                                     (for UNSPLIT_GRAVITY only)
//                [11] GravityType     : Types of gravity --> self-gravity, external gravity, both (for UNSPLIT_GRAVITY only)
//                [12] ExtAcc_AuxArray : Auxiliary array for external acceleration                 (for UNSPLIT_GRAVITY only)
//                [13] MinTemp         : Minimum allowed temerature
//                [14] DumpIntFlux     : true --> store the inter-patch fluxes in g_IntFlux[]
//                [15] g_IntFlux       : Array for DumpIntFlux
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE
void SRHydro_ComputeFlux( const real g_FC_Var [][NCOMP_TOTAL][ CUBE(N_FC_VAR) ],
                                real g_FC_Flux[][NCOMP_TOTAL][ CUBE(N_FC_FLUX) ],
                          const int Gap, const real Gamma, const bool CorrHalfVel, const real g_Pot_USG[],
                          const double g_Corner[], const real dt, const real dh, const double Time,
                          const OptGravityType_t GravityType, const double ExtAcc_AuxArray[], const real MinTemp,
                          const bool DumpIntFlux, real g_IntFlux[][NCOMP_TOTAL][ SQR(PS2) ] )
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

      if ( N_FC_VAR != PS2+2 )
         printf( "ERROR : N_FC_VAR (%d) != PS2+2 (%d) !!\n", N_FC_VAR, PS2+2 );
   }
#  else
   if ( CorrHalfVel )
      printf( "ERROR : CorrHalfVel is NOT supported when UNSPLIT_GRAVITY is off !!\n" );
#  endif
#  endif // #ifdef GAMER_DEBUG

   const int didx_fc[3] = { 1, N_FC_VAR, N_FC_VAR*N_FC_VAR };

   real ConVar_L[NCOMP_TOTAL], ConVar_R[NCOMP_TOTAL], Flux_1Face[NCOMP_TOTAL];

#  ifdef UNSPLIT_GRAVITY
   real PriVar_L[NCOMP_TOTAL], PriVar_R[NCOMP_TOTAL];

   const real   GraConst    = -(real)0.5*dt/dh;
   const int    didx_usg[3] = { 1, USG_NXT_F, SQR(USG_NXT_F) };
   const int    didx        = USG_GHOST_SIZE - 1;  // assuming g_FC_Var[] has one ghost zone on each side
   const double dh_half     = 0.5*(double)dh;      // always use double precision to calculate the cell position
   const real   dt_half     = (real)0.5*dt;

   double CrShift[3];

// assuming N_FC_VAR = PS2+2 (ghost-zone=1 on each side)
// --> CrShift is the central coordinates of the 1st cell in g_FC_Var[]
   if (  CorrHalfVel  &&  ( GravityType == GRAVITY_EXTERNAL || GravityType == GRAVITY_BOTH )  )
      for (int d=0; d<3; d++)    CrShift[d] = g_Corner[d] - (double)dh;
#  endif


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
#        ifdef UNSPLIT_GRAVITY
         real n_R, Ux_R, Uy_R, Uz_R, P_L, Uxx_R, Uyy_R, Uzz_R, Uxy_R, Uxz_R, Uyz_R, LorentzFactor_R;
         real n_L, Ux_L, Uy_L, Uz_L, P_R, Uxx_L, Uyy_L, Uzz_L, Uxy_L, Uxz_L, Uyz_L, LorentzFactor_L;

         real Const1, Const2, Const3;
#        endif

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


//       1. correct the half-step velocity by gravity
#        ifdef UNSPLIT_GRAVITY
         if ( CorrHalfVel )
         {
            real   Acc[3];
            double xyz[3];

            Acc[0] = (real)0.0;
            Acc[1] = (real)0.0;
            Acc[2] = (real)0.0;

//          external gravity
            if ( GravityType == GRAVITY_EXTERNAL  ||  GravityType == GRAVITY_BOTH )
            {
               xyz[0]  = CrShift[0] + (double)(i_fc*dh);
               xyz[1]  = CrShift[1] + (double)(j_fc*dh);
               xyz[2]  = CrShift[2] + (double)(k_fc*dh);
               xyz[d] += dh_half;

               ExternalAcc( Acc, xyz[0], xyz[1], xyz[2], Time, ExtAcc_AuxArray );

               for (int d=0; d<3; d++)    Acc[d] *= dt_half;
            }

//          self-gravity
            if ( GravityType == GRAVITY_SELF  ||  GravityType == GRAVITY_BOTH )
            {
               const int idx_usg = IDX321( i_fc+didx, j_fc+didx, k_fc+didx, USG_NXT_F, USG_NXT_F );

               Acc[d1] +=            GraConst*( g_Pot_USG[ idx_usg+didx_usg[d1] ] - g_Pot_USG[ idx_usg                           ] );
               Acc[d2] += (real)0.25*GraConst*( g_Pot_USG[ idx_usg+didx_usg[d2] ] + g_Pot_USG[ idx_usg+didx_usg[d2]+didx_usg[d1] ]
                                               -g_Pot_USG[ idx_usg-didx_usg[d2] ] - g_Pot_USG[ idx_usg-didx_usg[d2]+didx_usg[d1] ] );
               Acc[d3] += (real)0.25*GraConst*( g_Pot_USG[ idx_usg+didx_usg[d3] ] + g_Pot_USG[ idx_usg+didx_usg[d3]+didx_usg[d1] ]
                                               -g_Pot_USG[ idx_usg-didx_usg[d3] ] - g_Pot_USG[ idx_usg-didx_usg[d3]+didx_usg[d1] ] );
            }


//          conserved --> primitive variables
            LorentzFactor_L = SRHydro_Con2Pri( ConVar_L, PriVar_L, Gamma, MinTemp);
            LorentzFactor_R = SRHydro_Con2Pri( ConVar_R, PriVar_R, Gamma, MinTemp);

#           ifdef USE_3_VELOCITY
#           error: UNSPLIT_GRAVITY do not support 3-velocities
#           endif

//          update the momentum density (ConVar_L)
            n_L    = PriVar_L[0];
            Ux_L   = PriVar_L[1];
            Uy_L   = PriVar_L[2];
            Uz_L   = PriVar_L[3];
            P_L    = PriVar_L[4];

            Const1 = n_L*( (real)2.0 * SQR(LorentzFactor_L)-(real)1.0 );
            Const2 = ( SQR(Ux_L) + SQR(Uy_L) + SQR(Uz_L) ) * P_L;
            
            Uxx_L = Ux_L*Ux_L;
            Uyy_L = Uy_L*Uy_L;
            Uzz_L = Uz_L*Uz_L;
            Uxy_L = Ux_L*Uy_L;
            Uxz_L = Ux_L*Uz_L;
            Uyz_L = Uy_L*Uz_L;

            ConVar_L[MOMX] += - Uxx_L * Acc[0] - Uxy_L * Acc[1] - Uxz_L * Acc[2] + Const1 * Acc[0] + Const2 * Acc[0];
            ConVar_L[MOMY] += - Uxy_L * Acc[0] - Uyy_L * Acc[1] - Uyz_L * Acc[2] + Const1 * Acc[1] + Const2 * Acc[1];
            ConVar_L[MOMZ] += - Uxz_L * Acc[0] - Uyz_L * Acc[1] - Uzz_L * Acc[2] + Const1 * Acc[2] + Const2 * Acc[2];


//          update the momentum density (ConVar_R)
            n_R    = PriVar_R[0];
            Ux_R   = PriVar_R[1];
            Uy_R   = PriVar_R[2];
            Uz_R   = PriVar_R[3];
            P_R    = PriVar_R[4];

            Const1 = n_R*( (real)2.0 * SQR(LorentzFactor_R)-(real)1.0 );
            Const2 = ( SQR(Ux_R) + SQR(Uy_R) + SQR(Uz_R) ) * P_R;
            
            Uxx_R = Ux_R*Ux_R;
            Uyy_R = Uy_R*Uy_R;
            Uzz_R = Uz_R*Uz_R;
            Uxy_R = Ux_R*Uy_R;
            Uxz_R = Ux_R*Uz_R;
            Uyz_R = Uy_R*Uz_R;

            ConVar_R[MOMX] += - Uxx_R * Acc[0] - Uxy_R * Acc[1] - Uxz_R * Acc[2] + Const1 * Acc[0] + Const2 * Acc[0];
            ConVar_R[MOMY] += - Uxy_R * Acc[0] - Uyy_R * Acc[1] - Uyz_R * Acc[2] + Const1 * Acc[1] + Const2 * Acc[1];
            ConVar_R[MOMZ] += - Uxz_R * Acc[0] - Uyz_R * Acc[1] - Uzz_R * Acc[2] + Const1 * Acc[2] + Const2 * Acc[2];


//          update total energy density (ConVar_R)
            Const3 = LorentzFactor_L * ( n_L + P_L );
			ConVar_L[ENGY] += Const3 * ( Ux_L * Acc[0] + Uy_L * Acc[1] + Uz_L * Acc[2] );


//          update total energy density (ConVar_R)
            Const3 = LorentzFactor_R * ( n_R + P_R );
			ConVar_R[ENGY] += Const3 * ( Ux_R * Acc[0] + Uy_R * Acc[1] + Uz_R * Acc[2] );

         } // if ( CorrHalfVel )
#        endif // #ifdef UNSPLIT_GRAVITY

#        ifdef CHECK_FAILED_CELL_IN_FLUID
         SRHydro_CheckUnphysical( ConVar_L, NULL, Gamma, MinTemp, __FUNCTION__, __LINE__, true );
         SRHydro_CheckUnphysical( ConVar_R, NULL, Gamma, MinTemp, __FUNCTION__, __LINE__, true );
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
