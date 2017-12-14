#ifndef __CUFLU_COMPUTEFLUX_CU__
#define __CUFLU_COMPUTEFLUX_CU__



#include "Macro.h"
#include "CUFLU.h"
#if   ( RSOLVER == EXACT )
#include "CUFLU_Shared_RiemannSolver_Exact.cu"
#elif ( RSOLVER == ROE )
#include "CUFLU_Shared_RiemannSolver_Roe.cu"
#elif ( RSOLVER == HLLE )
#include "CUFLU_Shared_RiemannSolver_HLLE.cu"
#elif ( RSOLVER == HLLC )
#include "CUFLU_Shared_RiemannSolver_HLLC.cu"
#endif

#ifdef UNSPLIT_GRAVITY
#include "../../SelfGravity/GPU_Gravity/CUPOT_ExternalAcc.cu"
#endif


static __device__ void CUFLU_ComputeFlux( const real g_FC_Var_xL[][NCOMP_TOTAL][ N_FC_VAR*N_FC_VAR*N_FC_VAR ],
                                          const real g_FC_Var_xR[][NCOMP_TOTAL][ N_FC_VAR*N_FC_VAR*N_FC_VAR ],
                                          const real g_FC_Var_yL[][NCOMP_TOTAL][ N_FC_VAR*N_FC_VAR*N_FC_VAR ],
                                          const real g_FC_Var_yR[][NCOMP_TOTAL][ N_FC_VAR*N_FC_VAR*N_FC_VAR ],
                                          const real g_FC_Var_zL[][NCOMP_TOTAL][ N_FC_VAR*N_FC_VAR*N_FC_VAR ],
                                          const real g_FC_Var_zR[][NCOMP_TOTAL][ N_FC_VAR*N_FC_VAR*N_FC_VAR ],
                                          real g_FC_Flux_x[][NCOMP_TOTAL][ N_FC_FLUX*N_FC_FLUX*N_FC_FLUX ],
                                          real g_FC_Flux_y[][NCOMP_TOTAL][ N_FC_FLUX*N_FC_FLUX*N_FC_FLUX ],
                                          real g_FC_Flux_z[][NCOMP_TOTAL][ N_FC_FLUX*N_FC_FLUX*N_FC_FLUX ],
                                          real g_Flux[][9][NCOMP_TOTAL][ PS2*PS2 ], const bool DumpFlux,
                                          const uint Gap, const real Gamma, const real MinPres );




//-------------------------------------------------------------------------------------------------------
// Function    :  CUFLU_ComputeFlux
// Description :  Compute the face-centered fluxes by Riemann solver
//
// Note        :  1. Currently support the exact, Roe, HLLE, and HLLC solvers
//                2. Prefix "g" for pointers pointing to the "Global" memory space
//                   Prefix "s" for pointers pointing to the "Shared" memory space
//                3. The function is asynchronous
//                   --> "__syncthreads" must be called before using the output data
//                4. The sizes of the array g_FC_Var_x/y/z are assumed to be "N_FC_VAR"
//                   --> "N_FC_VAR-1" fluxes will be computed along each normal direction
//                5. The sizes of the arrays g_FC_Flux_XX in each direction are assumed to be "N_FL_FLUX"
//                   --> The (i,j,k) flux will be stored in the array g_FC_Flux_XX with
//                       the index "(k*N_FL_FLUX+j)*N_FL_FLUX+i"
//                   --> The (i,j,k) FC_Flux_x is defined at the +x surfaces of the cell (i,     j-Gap, k-Gap)
//                       The (i,j,k) FC_Flux_y is defined at the +y surfaces of the cell (i-Gap, j,     k-Gap)
//                       The (i,j,k) FC_Flux_z is defined at the +z surfaces of the cell (i-Gap, j-Gap, k    )
//                6. This function is shared by MHM, MHM_RP, and CTU schemes
//                7. For the performance consideration, this function will also be responsible for store the
//                   inter-patch fluxes
//                   --> Option "DumpFlux"
//                8. The "__forceinline__" qualifier is added for higher performance
//
// Parameter   :  g_FC_Var_xL     : Global memory array storing the face-centered variables on the -x surface
//                g_FC_Var_xR     : Global memory array storing the face-centered variables on the +x surface
//                g_FC_Var_yL     : Global memory array storing the face-centered variables on the -y surface
//                g_FC_Var_yR     : Global memory array storing the face-centered variables on the +y surface
//                g_FC_Var_zL     : Global memory array storing the face-centered variables on the -z surface
//                g_FC_Var_zR     : Global memory array storing the face-centered variables on the +z surface
//                g_FC_Flux_x     : Global memory array to store the face-centered fluxes in the x direction
//                g_FC_Flux_y     : Global memory array to store the face-centered fluxes in the y direction
//                g_FC_Flux_z     : Global memory array to store the face-centered fluxes in the z direction
//                g_Flux          : Global memory array to store the output fluxes
//                DumpFlux        : true --> store the inter-patch fluxes to "g_Flux" for the AMR correction
//                Gap             : Number of grids to be skipped in the transverse direction
//                                  --> "(N_FC_VAR-2*Gap)^2" fluxes will be computed in each surface
//                Gamma           : Ratio of specific heats
//                CorrHalfVel     : true --> correcting the half-step velocity by gravity (for UNSPLIT_GRAVITY only)
//                g_Pot_USG       : Global memory array storing the input potential for CorrHalfVel (for UNSPLIT_GRAVITY only)
//                                  --> must have the same size as FC_Var ( (PS2+2)^3 )
//                g_Corner        : Global memory array storing the physical corner coordinates of each patch group (for UNSPLIT_GRAVITY)
//                dt              : Time interval to advance the full-step solution (for UNSPLIT_GRAVITY only)
//                _dh             : 1 / Grid size                                   (for UNSPLIT_GRAVITY only)
//                Time            : Current physical time                           (for UNSPLIT_GRAVITY only)
//                GravityType     : Types of gravity --> self-gravity, external gravity, both (for UNSPLIT_GRAVITY only)
//                ExtAcc_AuxArray : Auxiliary array for adding external acceleration (for UNSPLIT_GRAVITY only)
//                MinPres         : Minimum allowed pressure
//-------------------------------------------------------------------------------------------------------
__forceinline__
__device__ void CUFLU_ComputeFlux( const real g_FC_Var_xL[][NCOMP_TOTAL][ N_FC_VAR*N_FC_VAR*N_FC_VAR ],
                                   const real g_FC_Var_xR[][NCOMP_TOTAL][ N_FC_VAR*N_FC_VAR*N_FC_VAR ],
                                   const real g_FC_Var_yL[][NCOMP_TOTAL][ N_FC_VAR*N_FC_VAR*N_FC_VAR ],
                                   const real g_FC_Var_yR[][NCOMP_TOTAL][ N_FC_VAR*N_FC_VAR*N_FC_VAR ],
                                   const real g_FC_Var_zL[][NCOMP_TOTAL][ N_FC_VAR*N_FC_VAR*N_FC_VAR ],
                                   const real g_FC_Var_zR[][NCOMP_TOTAL][ N_FC_VAR*N_FC_VAR*N_FC_VAR ],
                                   real g_FC_Flux_x[][NCOMP_TOTAL][ N_FC_FLUX*N_FC_FLUX*N_FC_FLUX ],
                                   real g_FC_Flux_y[][NCOMP_TOTAL][ N_FC_FLUX*N_FC_FLUX*N_FC_FLUX ],
                                   real g_FC_Flux_z[][NCOMP_TOTAL][ N_FC_FLUX*N_FC_FLUX*N_FC_FLUX ],
                                   real g_Flux[][9][NCOMP_TOTAL][ PS2*PS2 ], const bool DumpFlux,
                                   const uint Gap, const real Gamma, const bool CorrHalfVel,
                                   const real g_Pot_USG[][ USG_NXT_F*USG_NXT_F*USG_NXT_F ],
                                   const double g_Corner[][3], const real dt, const real _dh, const double Time,
                                   const OptGravityType_t GravityType, const double ExtAcc_AuxArray[], const real MinPres )
{

// check
#  ifdef UNSPLIT_GRAVITY
#  if ( N_FC_VAR != PS2+2 )
#     error : ERROR : N_FC_VAR != PS2+2 !!
#  endif

#  if ( USG_GHOST_SIZE < 1 )
#     error : ERROR : USG_GHOST_SIZE < 1 !!
#  endif
#  endif // #ifdef UNSPLIT_GRAVITY


   const uint  bx       = blockIdx.x;
   const uint  tx       = threadIdx.x;
   const uint  dID      = blockDim.x;
   const uint3 dID_In   = make_uint3( 1, N_FC_VAR, N_FC_VAR*N_FC_VAR );
   const uint  N_N      = N_FC_VAR - 1;      // number of fluxes to be calculated in the normal direction
   const uint  N_T      = N_FC_VAR - 2*Gap;  // number of fluxes to be calculated in the transverse direction

#  if ( RSOLVER == EXACT )
   const real  Gamma_m1 = Gamma - (real)1.0;
#  endif

   uint   ID, ID_In, ID_Out, Nxy, SurfID;
   uint3  ID3d;
   FluVar VarL, VarR, FC_Flux;

#  ifdef UNSPLIT_GRAVITY
   const uint   dID_USG[3] = { 1, USG_NXT_F, USG_NXT_F*USG_NXT_F };
   const real   GraConst   = -(real)0.5*dt*_dh;
   const uint   didx       = USG_GHOST_SIZE - 1;   // assuming FC_Var has one ghost zone on each side
   const double dh         = 1.0/_dh;
   const double dh_half    = 0.5*dh;
   const real   dt_half    = (real)0.5*dt;

   uint   d1, d2, d3, ID_USG;
   real   eL, eR, Acc[3];
   double xyz[3];
#  endif


// macro to correct the half-step velocity for unsplitting update
#  ifdef UNSPLIT_GRAVITY
#  define CorrHalfVel( Gap_x, Gap_y, Gap_z )                                                       \
   {                                                                                               \
      if ( CorrHalfVel )                                                                           \
      {                                                                                            \
         Acc[0] = (real)0.0;                                                                       \
         Acc[1] = (real)0.0;                                                                       \
         Acc[2] = (real)0.0;                                                                       \
                                                                                                   \
         /* external gravity */                                                                    \
         if ( GravityType == GRAVITY_EXTERNAL  ||  GravityType == GRAVITY_BOTH )                   \
         {                                                                                         \
            /* do not write "(ID3d.?+Gap_?-1)*dh" since ID3d and Gap are uint and 0U+0U-1 != -1 */ \
            xyz[0 ]  = g_Corner[bx][0] + (ID3d.x+Gap_x)*dh - dh;                                   \
            xyz[1 ]  = g_Corner[bx][1] + (ID3d.y+Gap_y)*dh - dh;                                   \
            xyz[2 ]  = g_Corner[bx][2] + (ID3d.z+Gap_z)*dh - dh;                                   \
            xyz[d1] += dh_half;                                                                    \
                                                                                                   \
            CUPOT_ExternalAcc( Acc, xyz[0], xyz[1], xyz[2], Time, ExtAcc_AuxArray );               \
                                                                                                   \
            Acc[0] *= dt_half;                                                                     \
            Acc[1] *= dt_half;                                                                     \
            Acc[2] *= dt_half;                                                                     \
         }                                                                                         \
                                                                                                   \
         /* self-gravity */                                                                        \
         if ( GravityType == GRAVITY_SELF  ||  GravityType == GRAVITY_BOTH )                       \
         {                                                                                         \
            ID_USG   = ( (ID3d.z+Gap_z+didx)*USG_NXT_F + (ID3d.y+Gap_y+didx) )*USG_NXT_F           \
                       + (ID3d.x+Gap_x+didx);                                                      \
                                                                                                   \
            Acc[d1] +=            GraConst*(  g_Pot_USG[bx][ ID_USG+dID_USG[d1]             ]      \
                                            - g_Pot_USG[bx][ ID_USG                         ] );   \
                                                                                                   \
            Acc[d2] += (real)0.25*GraConst*(  g_Pot_USG[bx][ ID_USG+dID_USG[d2]             ]      \
                                            + g_Pot_USG[bx][ ID_USG+dID_USG[d2]+dID_USG[d1] ]      \
                                            - g_Pot_USG[bx][ ID_USG-dID_USG[d2]             ]      \
                                            - g_Pot_USG[bx][ ID_USG-dID_USG[d2]+dID_USG[d1] ] );   \
                                                                                                   \
            Acc[d3] += (real)0.25*GraConst*(  g_Pot_USG[bx][ ID_USG+dID_USG[d3]             ]      \
                                            + g_Pot_USG[bx][ ID_USG+dID_USG[d3]+dID_USG[d1] ]      \
                                            - g_Pot_USG[bx][ ID_USG-dID_USG[d3]             ]      \
                                            - g_Pot_USG[bx][ ID_USG-dID_USG[d3]+dID_USG[d1] ] );   \
         }                                                                                         \
                                                                                                   \
         /* store the internal energy density */                                                   \
         eL = VarL.Egy - (real)0.5*( SQR(VarL.Px) + SQR(VarL.Py) + SQR(VarL.Pz) )/VarL.Rho;        \
         eR = VarR.Egy - (real)0.5*( SQR(VarR.Px) + SQR(VarR.Py) + SQR(VarR.Pz) )/VarR.Rho;        \
                                                                                                   \
         /* advance velocity by gravity */                                                         \
         VarL.Px += VarL.Rho*Acc[0];                                                               \
         VarR.Px += VarR.Rho*Acc[0];                                                               \
         VarL.Py += VarL.Rho*Acc[1];                                                               \
         VarR.Py += VarR.Rho*Acc[1];                                                               \
         VarL.Pz += VarL.Rho*Acc[2];                                                               \
         VarR.Pz += VarR.Rho*Acc[2];                                                               \
                                                                                                   \
         /* update total energy density with the internal energy density fixed */                  \
         VarL.Egy = eL + (real)0.5*( SQR(VarL.Px) + SQR(VarL.Py) + SQR(VarL.Pz) )/VarL.Rho;        \
         VarR.Egy = eR + (real)0.5*( SQR(VarR.Px) + SQR(VarR.Py) + SQR(VarR.Pz) )/VarR.Rho;        \
      } /* if ( CorrHalfVel ) */                                                                   \
   } // CorrHalfVel

#  else

#  define CorrHalfVel( Gap_x, Gap_y, Gap_z )    /* nothing to do */

#  endif // #ifdef UNSPLIT_GRAVITY ... else ...


// macro to load data from the global memory
#  define Load( Input, Output, ID )                            \
   {                                                           \
      Output.Rho = Input[bx][0][ID];                           \
      Output.Px  = Input[bx][1][ID];                           \
      Output.Py  = Input[bx][2][ID];                           \
      Output.Pz  = Input[bx][3][ID];                           \
      Output.Egy = Input[bx][4][ID];                           \
                                                               \
      for (int v=0; v<NCOMP_PASSIVE; v++)                      \
      Output.Passive[v] = Input[bx][ NCOMP_FLUID + v ][ID];    \
   } // Load


// macro to dump face-centered flux to the global memory
#  define Dump_FC_Flux( Input, Output, ID )                    \
   {                                                           \
      Output[bx][0][ID] = Input.Rho;                           \
      Output[bx][1][ID] = Input.Px;                            \
      Output[bx][2][ID] = Input.Py;                            \
      Output[bx][3][ID] = Input.Pz;                            \
      Output[bx][4][ID] = Input.Egy;                           \
                                                               \
      for (int v=0; v<NCOMP_PASSIVE; v++)                      \
      Output[bx][ NCOMP_FLUID + v ][ID] = Input.Passive[v];    \
   } // Dump_FC_Flux


// macro to dump fluxes across patch boundaries to the global memory
#  define Dump_InterPatch_Flux( FC_Flux, SurfID, ID_Flux )                    \
   {                                                                          \
      g_Flux[bx][SurfID][0][ID_Flux] = FC_Flux.Rho;                           \
      g_Flux[bx][SurfID][1][ID_Flux] = FC_Flux.Px;                            \
      g_Flux[bx][SurfID][2][ID_Flux] = FC_Flux.Py;                            \
      g_Flux[bx][SurfID][3][ID_Flux] = FC_Flux.Pz;                            \
      g_Flux[bx][SurfID][4][ID_Flux] = FC_Flux.Egy;                           \
                                                                              \
      for (int v=0; v<NCOMP_PASSIVE; v++)                                     \
      g_Flux[bx][SurfID][ NCOMP_FLUID + v ][ID_Flux] = FC_Flux.Passive[v];    \
   } // Dump_InterPatch_Flux


// macro for different Riemann solvers
#  if   ( RSOLVER == EXACT )

      /* exact solver */
      #define RiemannSolver( Dir, VarL, VarR )                                               \
      {                                                                                      \
         /* do NOT convert any passive variable to mass fraction for the Riemann solvers */  \
         const bool NormPassive_No  = false;                                                 \
         const bool JeansMinPres_No = false;                                                 \
                                                                                             \
         VarL = CUFLU_Con2Pri( VarL, Gamma_m1, MinPres, NormPassive_No, NULL_INT, NULL, JeansMinPres_No, NULL_REAL ); \
         VarR = CUFLU_Con2Pri( VarR, Gamma_m1, MinPres, NormPassive_No, NULL_INT, NULL, JeansMinPres_No, NULL_REAL ); \
                                                                                             \
         FC_Flux = CUFLU_RiemannSolver_Exact( Dir, NULL, NULL, NULL, VarL, VarR, Gamma );    \
      } // RiemannSolver

#  elif ( RSOLVER == ROE )

      /* Roe solver */
      #define RiemannSolver( Dir, VarL, VarR )                                               \
      {                                                                                      \
         FC_Flux = CUFLU_RiemannSolver_Roe( Dir, VarL, VarR, Gamma, MinPres );               \
      } // RiemannSolver

#  elif ( RSOLVER == HLLE )

      /* HLLE solver */
      #define RiemannSolver( Dir, VarL, VarR )                                               \
      {                                                                                      \
         FC_Flux = CUFLU_RiemannSolver_HLLE( Dir, VarL, VarR, Gamma, MinPres );              \
      } // RiemannSolver

#  elif ( RSOLVER == HLLC )

      /* HLLC solver */
      #define RiemannSolver( Dir, VarL, VarR )                                               \
      {                                                                                      \
         FC_Flux = CUFLU_RiemannSolver_HLLC( Dir, VarL, VarR, Gamma, MinPres );              \
      } // RiemannSolver

#  else

#     error : ERROR : unsupported Riemann solver (EXACT/ROE/HLLE/HLLC) !!

#  endif


// key macro to invoke other macros to calculate and store the fluxes
#  define GetFlux( Dir, Nx, Ny, Gap_x, Gap_y, Gap_z, dID_In, g_FC_Var_L, g_FC_Var_R, g_FC_Flux )               \
   {                                                                                                           \
      ID  = tx;                                                                                                \
      Nxy = (Nx)*(Ny);                                                                                         \
                                                                                                               \
      while ( ID < N_N*N_T*N_T )                                                                               \
      {                                                                                                        \
         ID3d.x = ID%(Nx);                                                                                     \
         ID3d.y = ID%Nxy/(Nx);                                                                                 \
         ID3d.z = ID/Nxy;                                                                                      \
         ID_In  = __umul24( __umul24( ID3d.z+Gap_z, N_FC_VAR  ) + ID3d.y+Gap_y, N_FC_VAR  ) + ID3d.x+Gap_x;    \
         ID_Out = __umul24( __umul24( ID3d.z,       N_FL_FLUX ) + ID3d.y,       N_FL_FLUX ) + ID3d.x;          \
                                                                                                               \
         Load( g_FC_Var_R, VarL, ID_In        );                                                               \
         Load( g_FC_Var_L, VarR, ID_In+dID_In );                                                               \
                                                                                                               \
         CorrHalfVel( Gap_x, Gap_y, Gap_z );                                                                   \
                                                                                                               \
         RiemannSolver( Dir, VarL, VarR );                                                                     \
                                                                                                               \
         Dump_FC_Flux( FC_Flux, g_FC_Flux, ID_Out );                                                           \
                                                                                                               \
         /* store the inter-patch fluxes */                                                                    \
         if ( DumpFlux )                                                                                       \
         {                                                                                                     \
            if      (  Dir == 0  &&  ( ID3d.x == 0 || ID3d.x == PS1 || ID3d.x == PS2 )  )                      \
            {                                                                                                  \
               SurfID = 3*Dir + ID3d.x/PS1;                                                                    \
               Dump_InterPatch_Flux( FC_Flux, SurfID, __umul24( ID3d.z, PS2 ) + ID3d.y );                      \
            }                                                                                                  \
                                                                                                               \
            else if (  Dir == 1  &&  ( ID3d.y == 0 || ID3d.y == PS1 || ID3d.y == PS2 )  )                      \
            {                                                                                                  \
               SurfID = 3*Dir + ID3d.y/PS1;                                                                    \
               Dump_InterPatch_Flux( FC_Flux, SurfID, __umul24( ID3d.z, PS2 ) + ID3d.x );                      \
            }                                                                                                  \
                                                                                                               \
            else if (  Dir == 2  &&  ( ID3d.z == 0 || ID3d.z == PS1 || ID3d.z == PS2 )  )                      \
            {                                                                                                  \
               SurfID = 3*Dir + ID3d.z/PS1;                                                                    \
               Dump_InterPatch_Flux( FC_Flux, SurfID, __umul24( ID3d.y, PS2 ) + ID3d.x );                      \
            }                                                                                                  \
         }                                                                                                     \
                                                                                                               \
         ID += dID;                                                                                            \
                                                                                                               \
      } /* while ( ID < N_N*N_T*N_T ) */                                                                       \
   } // GetFlux


// actual code lines to invoke various macros
#  ifdef UNSPLIT_GRAVITY
   if ( CorrHalfVel )   { d1 = 0;   d2 = 1;   d3 = 2; };
#  endif
   GetFlux( 0, N_N, N_T,   0, Gap, Gap, dID_In.x, g_FC_Var_xL, g_FC_Var_xR, g_FC_Flux_x );

#  ifdef UNSPLIT_GRAVITY
   if ( CorrHalfVel )   { d1 = 1;   d2 = 2;   d3 = 0; };
#  endif
   GetFlux( 1, N_T, N_N, Gap,   0, Gap, dID_In.y, g_FC_Var_yL, g_FC_Var_yR, g_FC_Flux_y );

#  ifdef UNSPLIT_GRAVITY
   if ( CorrHalfVel )   { d1 = 2;   d2 = 0;   d3 = 1; };
#  endif
   GetFlux( 2, N_T, N_T, Gap, Gap,   0, dID_In.z, g_FC_Var_zL, g_FC_Var_zR, g_FC_Flux_z );


// remove all macros used only in this function
#  undef CorrHalfVel
#  undef Load
#  undef Dump_FC_Flux
#  undef Dump_InterPatch_Flux
#  undef RiemannSolver
#  undef GetFlux

} // FUNCTION : CUFLU_ComputeFlux



#endif // #ifndef __CUFLU_COMPUTEFLUX_CU__
