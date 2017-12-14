#ifndef __CUFLU_DATARECONSTRUCTION_CU__
#define __CUFLU_DATARECONSTRUCTION_CU__



#include "Macro.h"
#include "CUFLU.h"
#include "CUFLU_Shared_FluUtility.cu"

static __device__ FluVar LimitSlope( const FluVar L2, const FluVar L1, const FluVar C0, const FluVar R1,
                                     const FluVar R2, const LR_Limiter_t LR_Limiter, const real MinMod_Coeff,
                                     const real EP_Coeff, const real Gamma, const int XYZ );
static __device__ void CUFLU_DataReconstruction( const real g_PriVar[][NCOMP_TOTAL][ FLU_NXT*FLU_NXT*FLU_NXT ],
                                                 real g_Slope_PPM_x[][NCOMP_TOTAL][ N_SLOPE_PPM*N_SLOPE_PPM*N_SLOPE_PPM ],
                                                 real g_Slope_PPM_y[][NCOMP_TOTAL][ N_SLOPE_PPM*N_SLOPE_PPM*N_SLOPE_PPM ],
                                                 real g_Slope_PPM_z[][NCOMP_TOTAL][ N_SLOPE_PPM*N_SLOPE_PPM*N_SLOPE_PPM ],
                                                 real g_FC_Var_xL[][NCOMP_TOTAL][ N_FC_VAR*N_FC_VAR*N_FC_VAR ],
                                                 real g_FC_Var_xR[][NCOMP_TOTAL][ N_FC_VAR*N_FC_VAR*N_FC_VAR ],
                                                 real g_FC_Var_yL[][NCOMP_TOTAL][ N_FC_VAR*N_FC_VAR*N_FC_VAR ],
                                                 real g_FC_Var_yR[][NCOMP_TOTAL][ N_FC_VAR*N_FC_VAR*N_FC_VAR ],
                                                 real g_FC_Var_zL[][NCOMP_TOTAL][ N_FC_VAR*N_FC_VAR*N_FC_VAR ],
                                                 real g_FC_Var_zR[][NCOMP_TOTAL][ N_FC_VAR*N_FC_VAR*N_FC_VAR ],
                                                 const int NIn, const int NGhost, const real Gamma,
                                                 const LR_Limiter_t LR_Limiter, const real MinMod_Coeff,
                                                 const real EP_Coeff, const real dt, const real _dh,
                                                 const real MinDens, const real MinPres,
                                                 const bool NormPassive, const int NNorm, const int NormIdx[] );
#ifdef CHAR_RECONSTRUCTION
static __device__ FluVar Pri2Char( FluVar InOut, const real Gamma, const real Rho, const real Pres,
                                   const int XYZ );
static __device__ FluVar Char2Pri( FluVar InOut, const real Gamma, const real Rho, const real Pres,
                                   const int XYZ );
#endif
#if ( FLU_SCHEME == MHM )
static __device__ void HancockPredict( real g_FC_Var_xL[][NCOMP_TOTAL][ N_FC_VAR*N_FC_VAR*N_FC_VAR ],
                                       real g_FC_Var_xR[][NCOMP_TOTAL][ N_FC_VAR*N_FC_VAR*N_FC_VAR ],
                                       real g_FC_Var_yL[][NCOMP_TOTAL][ N_FC_VAR*N_FC_VAR*N_FC_VAR ],
                                       real g_FC_Var_yR[][NCOMP_TOTAL][ N_FC_VAR*N_FC_VAR*N_FC_VAR ],
                                       real g_FC_Var_zL[][NCOMP_TOTAL][ N_FC_VAR*N_FC_VAR*N_FC_VAR ],
                                       real g_FC_Var_zR[][NCOMP_TOTAL][ N_FC_VAR*N_FC_VAR*N_FC_VAR ],
                                       FluVar FC_Var_xL, FluVar FC_Var_xR,
                                       FluVar FC_Var_yL, FluVar FC_Var_yR,
                                       FluVar FC_Var_zL, FluVar FC_Var_zR,
                                       const real dt, const real _dh, const real Gamma,
                                       const uint ID_Out, const FluVar P_In,
                                       const real MinDens, const real MinPres,
                                       const bool NormPassive, const int NNorm, const int NormIdx[] );
#endif
#if ( FLU_SCHEME == CTU )
static __device__ void Get_EigenSystem( const FluVar CC_Var, FluVar5 &EVal_x, FluVar5 &EVal_y, FluVar5 &EVal_z,
                                        FluVar5 &LEVec1, FluVar5 &LEVec2, FluVar5 &LEVec3, FluVar5 &LEVec4,
                                        FluVar5 &LEVec5, FluVar5 &REVec1, FluVar5 &REVec2, FluVar5 &REVec3,
                                        FluVar5 &REVec4, FluVar5 &REVec5, const real Gamma );
#endif




#if ( LR_SCHEME == PLM )
//-------------------------------------------------------------------------------------------------------
// Function    :  CUFLU_DataReconstruction
// Description :  Reconstruct the face-centered variables by the piecewise-linear method (PLM)
//
// Note        :  1. Use the parameter "LR_Limiter" to choose different slope limiters
//                2. The input data should be primitive variables
//                3. For the performance consideration, the output data will be conserved variables
//                4. The face-centered variables will be advanced by half time-step for the MHM and
//                   CTU schemes
//                5. The function is asynchronous
//                   --> "__syncthreads" must be called before using the output data
//                6. The PLM and PPM data reconstruction functions share the same function name
//                7. The data reconstruction can be applied to characteristic variables by
//                   defining "CHAR_RECONSTRUCTION"
//                8. This function is shared by MHM, MHM_RP, and CTU schemes
//
// Parameter   :  g_PriVar          : Global memory array storing the input primitive variables
//                g_Slope_PPM_x     : Global memory array to store the x-slope for the PPM reconstruction
//                g_Slope_PPM_y     : Global memory array to store the y-slope for the PPM reconstruction
//                g_Slope_PPM_z     : Global memory array to store the z-slope for the PPM reconstruction
//                g_FC_Var_xL       : Global memory array to store the face-centered variables on the -x surface
//                g_FC_Var_xR       : Global memory array to store the face-centered variables on the +x surface
//                g_FC_Var_yL       : Global memory array to store the face-centered variables on the -y surface
//                g_FC_Var_yR       : Global memory array to store the face-centered variables on the +y surface
//                g_FC_Var_zL       : Global memory array to store the face-centered variables on the -z surface
//                g_FC_Var_zR       : Global memory array to store the face-centered variables on the +z surface
//                NIn               : Size of the input array "PriVar" in one direction
//                                    (can be smaller than FLU_NXT)
//                NGhost            : Size of the ghost zone
//                                    --> "NIn-2*NGhost" cells will be computed along each direction
//                                    --> The size of the output array "g_FC_Var" is assumed to be
//                                        "(NIn-2*NGhost)^3"
//                                    --> The reconstructed data at cell (i,j,k) will be stored in the
//                                        array "g_FC_Var" with the index "(i-NGhost,j-NGhost,k-NGhost)
//                Gamma             : Ratio of specific heats
//                LR_Limiter        : Slope limiter for the data reconstruction in the MHM/MHM_RP/CTU schemes
//                                    (0/1/2/3/4) = (vanLeer/generalized MinMod/vanAlbada/
//                                                   vanLeer + generalized MinMod/extrema-preserving) limiter
//                MinMod_Coeff      : Coefficient of the generalized MinMod limiter
//                EP_Coeff          : Coefficient of the extrema-preserving limiter
//                dt                : Time interval to advance solution (useful only for CTU and MHM schemes)
//                _dh               : 1 / grid size (useful only for CTU and MHM schemes)
//                MinDens/Pres      : Minimum allowed density and pressure
//                NormPassive       : true --> convert passive scalars to mass fraction
//                NNorm             : Number of passive scalars for the option "NormPassive"
//                                    --> Should be set to the global variable "PassiveNorm_NVar"
//                NormIdx           : Target variable indices for the option "NormPassive"
//                                    --> Should be set to the global variable "PassiveNorm_VarIdx"
//-------------------------------------------------------------------------------------------------------
__device__ void CUFLU_DataReconstruction( const real g_PriVar[][NCOMP_TOTAL][ FLU_NXT*FLU_NXT*FLU_NXT ],
                                          real g_Slope_PPM_x[][NCOMP_TOTAL][ N_SLOPE_PPM*N_SLOPE_PPM*N_SLOPE_PPM ],
                                          real g_Slope_PPM_y[][NCOMP_TOTAL][ N_SLOPE_PPM*N_SLOPE_PPM*N_SLOPE_PPM ],
                                          real g_Slope_PPM_z[][NCOMP_TOTAL][ N_SLOPE_PPM*N_SLOPE_PPM*N_SLOPE_PPM ],
                                          real g_FC_Var_xL[][NCOMP_TOTAL][ N_FC_VAR*N_FC_VAR*N_FC_VAR ],
                                          real g_FC_Var_xR[][NCOMP_TOTAL][ N_FC_VAR*N_FC_VAR*N_FC_VAR ],
                                          real g_FC_Var_yL[][NCOMP_TOTAL][ N_FC_VAR*N_FC_VAR*N_FC_VAR ],
                                          real g_FC_Var_yR[][NCOMP_TOTAL][ N_FC_VAR*N_FC_VAR*N_FC_VAR ],
                                          real g_FC_Var_zL[][NCOMP_TOTAL][ N_FC_VAR*N_FC_VAR*N_FC_VAR ],
                                          real g_FC_Var_zR[][NCOMP_TOTAL][ N_FC_VAR*N_FC_VAR*N_FC_VAR ],
                                          const int NIn, const int NGhost, const real Gamma,
                                          const LR_Limiter_t LR_Limiter, const real MinMod_Coeff,
                                          const real EP_Coeff, const real dt, const real _dh,
                                          const real MinDens, const real MinPres,
                                          const bool NormPassive, const int NNorm, const int NormIdx[] )
{

   const uint  bx       = blockIdx.x;
   const uint  tx       = threadIdx.x;
   const uint  dID2     = blockDim.x;
   const uint  NOut     = NIn - 2*NGhost;
   const uint  NOut2    = __umul24( NOut, NOut );
   const uint3 dr1      = make_uint3( 1, NIn, NIn*NIn );
   const real _Gamma_m1 = (real)1.0 / ( Gamma - (real)1.0 );

   uint   ID2 = tx;
   uint   ID1, ID1_L, ID1_R;
   uint3  ID1_3d;
   real   Min, Max;
   FluVar Slope_Limiter;
   FluVar CC_L, CC_C, CC_R;                           // cell-centered variables
   FluVar FC_xL, FC_xR, FC_yL, FC_yR, FC_zL, FC_zR;   // face-centered variables

// variables for the extrema-preserving limiter
   uint   ID1_L2, ID1_R2;
   FluVar CC_L2, CC_R2;

// variables for the CTU scheme
#  if ( FLU_SCHEME == CTU )
   const real dt_dh2 = (real)0.5*dt*_dh;

   FluVar5 EVal_x, EVal_y, EVal_z;
   FluVar  Correct_L, Correct_R, dFC;
   real    Coeff_L, Coeff_R;

// initialize the constant components of the matrices of the left and right eigenvectors
   FluVar5 LEVec1 = { 0.0, NULL_REAL, 0.0, 0.0, NULL_REAL };
   FluVar5 LEVec2 = { 1.0,       0.0, 0.0, 0.0, NULL_REAL };
   FluVar5 LEVec3 = { 0.0,       0.0, 1.0, 0.0,       0.0 };
   FluVar5 LEVec4 = { 0.0,       0.0, 0.0, 1.0,       0.0 };
   FluVar5 LEVec5 = { 0.0, NULL_REAL, 0.0, 0.0, NULL_REAL };

   FluVar5 REVec1 = { 1.0, NULL_REAL, 0.0, 0.0, NULL_REAL };
   FluVar5 REVec2 = { 1.0,       0.0, 0.0, 0.0,       0.0 };
   FluVar5 REVec3 = { 0.0,       0.0, 1.0, 0.0,       0.0 };
   FluVar5 REVec4 = { 0.0,       0.0, 0.0, 1.0,       0.0 };
   FluVar5 REVec5 = { 1.0, NULL_REAL, 0.0, 0.0, NULL_REAL };
#  endif // CTU


// loop over all cells
   while ( ID2 < NOut*NOut*NOut )
   {
      ID1_3d.x = ID2%NOut       + NGhost;
      ID1_3d.y = ID2%NOut2/NOut + NGhost;
      ID1_3d.z = ID2/NOut2      + NGhost;
      ID1      = __umul24(  __umul24( ID1_3d.z, NIn ) + ID1_3d.y, NIn  ) + ID1_3d.x;

      CC_C.Rho        = g_PriVar[bx][ 0][ID1];
      CC_C.Px         = g_PriVar[bx][ 1][ID1];
      CC_C.Py         = g_PriVar[bx][ 2][ID1];
      CC_C.Pz         = g_PriVar[bx][ 3][ID1];
      CC_C.Egy        = g_PriVar[bx][ 4][ID1];
#     if ( NCOMP_PASSIVE > 0 )
      for (int v=0, vv=NCOMP_FLUID; v<NCOMP_PASSIVE; v++, vv++)
      CC_C.Passive[v] = g_PriVar[bx][vv][ID1];
#     endif


//    1. evaluate the eigenvalues and eigenvectors for the CTU integrator
#     if ( FLU_SCHEME == CTU )
      Get_EigenSystem( CC_C, EVal_x, EVal_y, EVal_z, LEVec1, LEVec2, LEVec3, LEVec4, LEVec5,
                       REVec1, REVec2, REVec3, REVec4, REVec5, Gamma );
#     endif


#     define Interpolate_PLM( comp, FC_L, FC_R )                                                      \
      {                                                                                               \
         FC_L.comp = CC_C.comp - (real)0.5*Slope_Limiter.comp;                                        \
         FC_R.comp = CC_C.comp + (real)0.5*Slope_Limiter.comp;                                        \
      } // Interpolate_PLM

#     define Monotonic_PLM( comp, FC_L, FC_R, CC_L, CC_C, CC_R )                                      \
      {                                                                                               \
         Min       = ( CC_C.comp < CC_L.comp ) ? CC_C.comp : CC_L.comp;                               \
         Max       = ( CC_C.comp > CC_L.comp ) ? CC_C.comp : CC_L.comp;                               \
         FC_L.comp = ( FC_L.comp > Min       ) ? FC_L.comp : Min;                                     \
         FC_L.comp = ( FC_L.comp < Max       ) ? FC_L.comp : Max;                                     \
         FC_R.comp = (real)2.0*CC_C.comp - FC_L.comp;                                                 \
                                                                                                      \
         Min       = ( CC_C.comp < CC_R.comp ) ? CC_C.comp : CC_R.comp;                               \
         Max       = ( CC_C.comp > CC_R.comp ) ? CC_C.comp : CC_R.comp;                               \
         FC_R.comp = ( FC_R.comp > Min       ) ? FC_R.comp : Min;                                     \
         FC_R.comp = ( FC_R.comp < Max       ) ? FC_R.comp : Max;                                     \
         FC_L.comp = (real)2.0*CC_C.comp - FC_R.comp;                                                 \
      } // Monotonic_PLM

#     define GetFCVar_PLM( XYZ, Dir, FC_L, FC_R )                                                     \
      {                                                                                               \
         ID1_L = ID1 - dr1.XYZ;                                                                       \
         ID1_R = ID1 + dr1.XYZ;                                                                       \
                                                                                                      \
         /* 2. load the cell-centered values */                                                       \
         CC_L.Rho        = g_PriVar[bx][ 0][ID1_L];                                                   \
         CC_L.Px         = g_PriVar[bx][ 1][ID1_L];                                                   \
         CC_L.Py         = g_PriVar[bx][ 2][ID1_L];                                                   \
         CC_L.Pz         = g_PriVar[bx][ 3][ID1_L];                                                   \
         CC_L.Egy        = g_PriVar[bx][ 4][ID1_L];                                                   \
         for (int v=0, vv=NCOMP_FLUID; v<NCOMP_PASSIVE; v++, vv++)                                    \
         CC_L.Passive[v] = g_PriVar[bx][vv][ID1_L];                                                   \
                                                                                                      \
         CC_R.Rho        = g_PriVar[bx][ 0][ID1_R];                                                   \
         CC_R.Px         = g_PriVar[bx][ 1][ID1_R];                                                   \
         CC_R.Py         = g_PriVar[bx][ 2][ID1_R];                                                   \
         CC_R.Pz         = g_PriVar[bx][ 3][ID1_R];                                                   \
         CC_R.Egy        = g_PriVar[bx][ 4][ID1_R];                                                   \
         for (int v=0, vv=NCOMP_FLUID; v<NCOMP_PASSIVE; v++, vv++)                                    \
         CC_R.Passive[v] = g_PriVar[bx][vv][ID1_R];                                                   \
                                                                                                      \
         if ( LR_Limiter == EXTPRE )                                                                  \
         {                                                                                            \
            ID1_L2 = ID1 - 2*dr1.XYZ;                                                                 \
            ID1_R2 = ID1 + 2*dr1.XYZ;                                                                 \
                                                                                                      \
            CC_L2.Rho        = g_PriVar[bx][ 0][ID1_L2];                                              \
            CC_L2.Px         = g_PriVar[bx][ 1][ID1_L2];                                              \
            CC_L2.Py         = g_PriVar[bx][ 2][ID1_L2];                                              \
            CC_L2.Pz         = g_PriVar[bx][ 3][ID1_L2];                                              \
            CC_L2.Egy        = g_PriVar[bx][ 4][ID1_L2];                                              \
            for (int v=0, vv=NCOMP_FLUID; v<NCOMP_PASSIVE; v++, vv++)                                 \
            CC_L2.Passive[v] = g_PriVar[bx][vv][ID1_L2];                                              \
                                                                                                      \
            CC_R2.Rho        = g_PriVar[bx][ 0][ID1_R2];                                              \
            CC_R2.Px         = g_PriVar[bx][ 1][ID1_R2];                                              \
            CC_R2.Py         = g_PriVar[bx][ 2][ID1_R2];                                              \
            CC_R2.Pz         = g_PriVar[bx][ 3][ID1_R2];                                              \
            CC_R2.Egy        = g_PriVar[bx][ 4][ID1_R2];                                              \
            for (int v=0, vv=NCOMP_FLUID; v<NCOMP_PASSIVE; v++, vv++)                                 \
            CC_R2.Passive[v] = g_PriVar[bx][vv][ID1_R2];                                              \
         }                                                                                            \
                                                                                                      \
                                                                                                      \
         /* 3. evaluate the monotonic slope */                                                        \
         Slope_Limiter = LimitSlope( CC_L2, CC_L, CC_C, CC_R, CC_R2, LR_Limiter, MinMod_Coeff,        \
                                     EP_Coeff, Gamma, Dir );                                          \
                                                                                                      \
                                                                                                      \
         /* 4. evaluate the face-centered values */                                                   \
         Interpolate_PLM( Rho,        FC_L, FC_R );                                                   \
         Interpolate_PLM( Px,         FC_L, FC_R );                                                   \
         Interpolate_PLM( Py,         FC_L, FC_R );                                                   \
         Interpolate_PLM( Pz,         FC_L, FC_R );                                                   \
         Interpolate_PLM( Egy,        FC_L, FC_R );                                                   \
         for (int v=0; v<NCOMP_PASSIVE; v++)                                                          \
         Interpolate_PLM( Passive[v], FC_L, FC_R );                                                   \
                                                                                                      \
                                                                                                      \
         /* ensure the face-centered variables lie between neighboring cell-centered values */        \
         if ( LR_Limiter != EXTPRE )                                                                  \
         {                                                                                            \
            Monotonic_PLM( Rho,        FC_L, FC_R, CC_L, CC_C, CC_R );                                \
            Monotonic_PLM( Px,         FC_L, FC_R, CC_L, CC_C, CC_R );                                \
            Monotonic_PLM( Py,         FC_L, FC_R, CC_L, CC_C, CC_R );                                \
            Monotonic_PLM( Pz,         FC_L, FC_R, CC_L, CC_C, CC_R );                                \
            Monotonic_PLM( Egy,        FC_L, FC_R, CC_L, CC_C, CC_R );                                \
            for (int v=0; v<NCOMP_PASSIVE; v++)                                                       \
            Monotonic_PLM( Passive[v], FC_L, FC_R, CC_L, CC_C, CC_R );                                \
         }                                                                                            \
                                                                                                      \
         else /* for the extrema-preserving limiter --> ensure positive density and pressure */       \
         {                                                                                            \
            FC_L.Rho = FMAX( FC_L.Rho, MinDens );                                                     \
            FC_R.Rho = FMAX( FC_R.Rho, MinDens );                                                     \
                                                                                                      \
            FC_L.Egy = CUFLU_CheckMinPres( FC_L.Egy, MinPres );                                       \
            FC_R.Egy = CUFLU_CheckMinPres( FC_R.Egy, MinPres );                                       \
                                                                                                      \
            /* ensure positive mass fractions for passive scalars */                                  \
            for (int v=0; v<NCOMP_PASSIVE; v++)                                                       \
            {                                                                                         \
               FC_L.Passive[v] = FMAX( FC_L.Passive[v], TINY_NUMBER );                                \
               FC_R.Passive[v] = FMAX( FC_R.Passive[v], TINY_NUMBER );                                \
            }                                                                                         \
         }                                                                                            \
                                                                                                      \
      } // GetFCVar_PLM

      GetFCVar_PLM( x, 0, FC_xL, FC_xR );
      GetFCVar_PLM( y, 1, FC_yL, FC_yR );
      GetFCVar_PLM( z, 2, FC_zL, FC_zR );

#     undef Interpolate_PLM
#     undef Monotonic_PLM
#     undef GetFCVar_PLM


//    6. advance the face-centered variables by half time-step for the CTU integrator
#     if ( FLU_SCHEME == CTU )

//    =====================================================================================
//    a. for the HLL solvers (HLLE/HLLC)
//    =====================================================================================
#     if (  ( RSOLVER == HLLE  ||  RSOLVER == HLLC )  &&  defined HLL_NO_REF_STATE  )

#     define GetRefState_PLM( EVal )                                                                  \
      {                                                                                               \
         Correct_L.Rho = (real)0.0;                                                                   \
         Correct_L.Px  = (real)0.0;                                                                   \
         Correct_L.Py  = (real)0.0;                                                                   \
         Correct_L.Pz  = (real)0.0;                                                                   \
         Correct_L.Egy = (real)0.0;                                                                   \
                                                                                                      \
         Correct_R.Rho = (real)0.0;                                                                   \
         Correct_R.Px  = (real)0.0;                                                                   \
         Correct_R.Py  = (real)0.0;                                                                   \
         Correct_R.Pz  = (real)0.0;                                                                   \
         Correct_R.Egy = (real)0.0;                                                                   \
      } // GetRefState_PLM


//    include waves both from left and right directions during the data reconstruction, as suggested in ATHENA
#     ifdef HLL_INCLUDE_ALL_WAVES

#     define Corr_L_PLM( comp, EVal, LEVec, REVec )                                                   \
      {                                                                                               \
         Coeff_L = (real)0.0;                                                                         \
                                                                                                      \
         Coeff_L += LEVec.Rho*dFC.Rho;                                                                \
         Coeff_L += LEVec.Px *dFC.Px;                                                                 \
         Coeff_L += LEVec.Py *dFC.Py;                                                                 \
         Coeff_L += LEVec.Pz *dFC.Pz;                                                                 \
         Coeff_L += LEVec.Egy*dFC.Egy;                                                                \
                                                                                                      \
         Coeff_L *= -dt_dh2*EVal.comp;                                                                \
                                                                                                      \
         Correct_L.Rho += Coeff_L*REVec.Rho;                                                          \
         Correct_L.Px  += Coeff_L*REVec.Px;                                                           \
         Correct_L.Py  += Coeff_L*REVec.Py;                                                           \
         Correct_L.Pz  += Coeff_L*REVec.Pz;                                                           \
         Correct_L.Egy += Coeff_L*REVec.Egy;                                                          \
      } // Corr_L_PLM

#     define Corr_R_PLM( comp, EVal, LEVec, REVec )                                                   \
      {                                                                                               \
         Correct_R.comp = Correct_L.comp;                                                             \
      } // Corr_R_PLM


//    only include waves propagating toward the target cell interface during the data reconstruction
#     else // ifndef HLL_INCLUDE_ALL_WAVES

#     define Corr_L_PLM( comp, EVal, LEVec, REVec )                                                   \
      {                                                                                               \
         Coeff_L = (real)0.0;                                                                         \
                                                                                                      \
         if ( EVal.comp <= (real)0.0 )                                                                \
         {                                                                                            \
            Coeff_L += LEVec.Rho*dFC.Rho;                                                             \
            Coeff_L += LEVec.Px *dFC.Px;                                                              \
            Coeff_L += LEVec.Py *dFC.Py;                                                              \
            Coeff_L += LEVec.Pz *dFC.Pz;                                                              \
            Coeff_L += LEVec.Egy*dFC.Egy;                                                             \
                                                                                                      \
            Coeff_L *= -dt_dh2*EVal.comp;                                                             \
                                                                                                      \
            Correct_L.Rho += Coeff_L*REVec.Rho;                                                       \
            Correct_L.Px  += Coeff_L*REVec.Px;                                                        \
            Correct_L.Py  += Coeff_L*REVec.Py;                                                        \
            Correct_L.Pz  += Coeff_L*REVec.Pz;                                                        \
            Correct_L.Egy += Coeff_L*REVec.Egy;                                                       \
         }                                                                                            \
      } // Corr_L_PLM

#     define Corr_R_PLM( comp, EVal, LEVec, REVec )                                                   \
      {                                                                                               \
         Coeff_R = (real)0.0;                                                                         \
                                                                                                      \
         if ( EVal.comp >= (real)0.0 )                                                                \
         {                                                                                            \
            Coeff_R += LEVec.Rho*dFC.Rho;                                                             \
            Coeff_R += LEVec.Px *dFC.Px;                                                              \
            Coeff_R += LEVec.Py *dFC.Py;                                                              \
            Coeff_R += LEVec.Pz *dFC.Pz;                                                              \
            Coeff_R += LEVec.Egy*dFC.Egy;                                                             \
                                                                                                      \
            Coeff_R *= -dt_dh2*EVal.comp;                                                             \
                                                                                                      \
            Correct_R.Rho += Coeff_R*REVec.Rho;                                                       \
            Correct_R.Px  += Coeff_R*REVec.Px;                                                        \
            Correct_R.Py  += Coeff_R*REVec.Py;                                                        \
            Correct_R.Pz  += Coeff_R*REVec.Pz;                                                        \
            Correct_R.Egy += Coeff_R*REVec.Egy;                                                       \
         }                                                                                            \
      } // Corr_R_PLM

#     endif // ifdef HLL_INCLUDE_ALL_WAVES ... else ...


//    =====================================================================================
//    b. for the Roe's and exact solvers
//    =====================================================================================
#     else // ( RSOLVER == ROE/EXACT && ifndef HLL_NO_REF_STATE )

#     define GetRefState_PLM( EVal )                                                                  \
      {                                                                                               \
         Coeff_L = -dt_dh2*FMIN( EVal.Rho, (real)0.0 );                                               \
         Coeff_R = -dt_dh2*FMAX( EVal.Egy, (real)0.0 );                                               \
                                                                                                      \
         Correct_L.Rho = Coeff_L*dFC.Rho;                                                             \
         Correct_L.Px  = Coeff_L*dFC.Px;                                                              \
         Correct_L.Py  = Coeff_L*dFC.Py;                                                              \
         Correct_L.Pz  = Coeff_L*dFC.Pz;                                                              \
         Correct_L.Egy = Coeff_L*dFC.Egy;                                                             \
                                                                                                      \
         Correct_R.Rho = Coeff_R*dFC.Rho;                                                             \
         Correct_R.Px  = Coeff_R*dFC.Px;                                                              \
         Correct_R.Py  = Coeff_R*dFC.Py;                                                              \
         Correct_R.Pz  = Coeff_R*dFC.Pz;                                                              \
         Correct_R.Egy = Coeff_R*dFC.Egy;                                                             \
      } // GetRefState_PLM

#     define Corr_L_PLM( comp, EVal, LEVec, REVec )                                                   \
      {                                                                                               \
         Coeff_L = (real)0.0;                                                                         \
                                                                                                      \
         if ( EVal.comp <= (real)0.0 )                                                                \
         {                                                                                            \
            Coeff_L += LEVec.Rho*dFC.Rho;                                                             \
            Coeff_L += LEVec.Px *dFC.Px;                                                              \
            Coeff_L += LEVec.Py *dFC.Py;                                                              \
            Coeff_L += LEVec.Pz *dFC.Pz;                                                              \
            Coeff_L += LEVec.Egy*dFC.Egy;                                                             \
                                                                                                      \
            Coeff_L *= dt_dh2*( EVal.Rho - EVal.comp );                                               \
                                                                                                      \
            Correct_L.Rho += Coeff_L*REVec.Rho;                                                       \
            Correct_L.Px  += Coeff_L*REVec.Px;                                                        \
            Correct_L.Py  += Coeff_L*REVec.Py;                                                        \
            Correct_L.Pz  += Coeff_L*REVec.Pz;                                                        \
            Correct_L.Egy += Coeff_L*REVec.Egy;                                                       \
         }                                                                                            \
      } // Corr_L_PLM

#     define Corr_R_PLM( comp, EVal, LEVec, REVec )                                                   \
      {                                                                                               \
         Coeff_R = (real)0.0;                                                                         \
                                                                                                      \
         if ( EVal.comp >= (real)0.0 )                                                                \
         {                                                                                            \
            Coeff_R += LEVec.Rho*dFC.Rho;                                                             \
            Coeff_R += LEVec.Px *dFC.Px;                                                              \
            Coeff_R += LEVec.Py *dFC.Py;                                                              \
            Coeff_R += LEVec.Pz *dFC.Pz;                                                              \
            Coeff_R += LEVec.Egy*dFC.Egy;                                                             \
                                                                                                      \
            Coeff_R *= dt_dh2*( EVal.Egy - EVal.comp );                                               \
                                                                                                      \
            Correct_R.Rho += Coeff_R*REVec.Rho;                                                       \
            Correct_R.Px  += Coeff_R*REVec.Px;                                                        \
            Correct_R.Py  += Coeff_R*REVec.Py;                                                        \
            Correct_R.Pz  += Coeff_R*REVec.Pz;                                                        \
            Correct_R.Egy += Coeff_R*REVec.Egy;                                                       \
         }                                                                                            \
      } // Corr_R_PLM


#     endif // if (  ( RSOLVER == HLLE  ||  RSOLVER == HLLC )  &&  defined HLL_NO_REF_STATE  ) ... else ...


//    macro for correcting the passive scalars
#     if ( NCOMP_PASSIVE > 0 )

#     define Corr_Passive( EVal )                                                                     \
      {                                                                                               \
         Coeff_L = -dt_dh2*FMIN( EVal.Px, (real)0.0 );                                                \
         Coeff_R = -dt_dh2*FMAX( EVal.Px, (real)0.0 );                                                \
                                                                                                      \
         for (int v=0; v<NCOMP_PASSIVE; v++)                                                          \
         {                                                                                            \
            Correct_L.Passive[v] = Coeff_L*dFC.Passive[v];                                            \
            Correct_R.Passive[v] = Coeff_R*dFC.Passive[v];                                            \
         }                                                                                            \
      }                                                                                               \

#     else

#     define Corr_Passive( EVal )   /* nothing to do */

#     endif // #if ( NCOMP_PASSIVE > 0 ) ... else ...


#     define CTU_Predict_PLM( dir, FC_L, FC_R, EVal )                                                 \
      {                                                                                               \
         /* (6-1) evaluate the slope */                                                               \
         dFC.Rho        = FC_R.Rho        - FC_L.Rho;                                                 \
         dFC.Px         = FC_R.Px         - FC_L.Px;                                                  \
         dFC.Py         = FC_R.Py         - FC_L.Py;                                                  \
         dFC.Pz         = FC_R.Pz         - FC_L.Pz;                                                  \
         dFC.Egy        = FC_R.Egy        - FC_L.Egy;                                                 \
         for (int v=0; v<NCOMP_PASSIVE; v++)                                                          \
         dFC.Passive[v] = FC_R.Passive[v] - FC_L.Passive[v];                                          \
                                                                                                      \
                                                                                                      \
         /* (6-2) re-order variables for the y/z directions */                                        \
         dFC = CUFLU_Rotate3D( dFC, dir, true );                                                      \
                                                                                                      \
                                                                                                      \
         /* (6-3) evaluate the reference states */                                                    \
         GetRefState_PLM( EVal );                                                                     \
                                                                                                      \
                                                                                                      \
         /* (6-4) evaluate the corrections to the left and right face-centered active variables */    \
         Corr_L_PLM( Rho, EVal, LEVec1, REVec1 );                                                     \
         Corr_L_PLM( Px,  EVal, LEVec2, REVec2 );                                                     \
         Corr_L_PLM( Py,  EVal, LEVec3, REVec3 );                                                     \
         Corr_L_PLM( Pz,  EVal, LEVec4, REVec4 );                                                     \
         Corr_L_PLM( Egy, EVal, LEVec5, REVec5 );                                                     \
                                                                                                      \
         Corr_R_PLM( Rho, EVal, LEVec1, REVec1 );                                                     \
         Corr_R_PLM( Px,  EVal, LEVec2, REVec2 );                                                     \
         Corr_R_PLM( Py,  EVal, LEVec3, REVec3 );                                                     \
         Corr_R_PLM( Pz,  EVal, LEVec4, REVec4 );                                                     \
         Corr_R_PLM( Egy, EVal, LEVec5, REVec5 );                                                     \
                                                                                                      \
                                                                                                      \
         /* (6-5) evaluate the corrections to the left and right face-centered passive scalars */     \
         /*       --> passive scalars travel with fluid velocity (i.e., entropy mode) */              \
         Corr_Passive( EVal );                                                                        \
                                                                                                      \
                                                                                                      \
         /* (6-6) evolve the face-centered variables by half time-step */                             \
         Correct_L = CUFLU_Rotate3D( Correct_L, dir, false );                                         \
         Correct_R = CUFLU_Rotate3D( Correct_R, dir, false );                                         \
                                                                                                      \
         FC_L.Rho        += Correct_L.Rho;                                                            \
         FC_L.Px         += Correct_L.Px;                                                             \
         FC_L.Py         += Correct_L.Py;                                                             \
         FC_L.Pz         += Correct_L.Pz;                                                             \
         FC_L.Egy        += Correct_L.Egy;                                                            \
         for (int v=0; v<NCOMP_PASSIVE; v++)                                                          \
         FC_L.Passive[v] += Correct_L.Passive[v];                                                     \
                                                                                                      \
         FC_R.Rho        += Correct_R.Rho;                                                            \
         FC_R.Px         += Correct_R.Px;                                                             \
         FC_R.Py         += Correct_R.Py;                                                             \
         FC_R.Pz         += Correct_R.Pz;                                                             \
         FC_R.Egy        += Correct_R.Egy;                                                            \
         for (int v=0; v<NCOMP_PASSIVE; v++)                                                          \
         FC_R.Passive[v] += Correct_R.Passive[v];                                                     \
                                                                                                      \
      } // CTU_Predict_PLM

      CTU_Predict_PLM( 0, FC_xL, FC_xR, EVal_x );
      CTU_Predict_PLM( 1, FC_yL, FC_yR, EVal_y );
      CTU_Predict_PLM( 2, FC_zL, FC_zR, EVal_z );


//    (6-7) ensure positive density and pressure
      FC_xL.Rho = FMAX( FC_xL.Rho, MinDens );
      FC_xR.Rho = FMAX( FC_xR.Rho, MinDens );
      FC_yL.Rho = FMAX( FC_yL.Rho, MinDens );
      FC_yR.Rho = FMAX( FC_yR.Rho, MinDens );
      FC_zL.Rho = FMAX( FC_zL.Rho, MinDens );
      FC_zR.Rho = FMAX( FC_zR.Rho, MinDens );

      FC_xL.Egy = CUFLU_CheckMinPres( FC_xL.Egy, MinPres );
      FC_xR.Egy = CUFLU_CheckMinPres( FC_xR.Egy, MinPres );
      FC_yL.Egy = CUFLU_CheckMinPres( FC_yL.Egy, MinPres );
      FC_yR.Egy = CUFLU_CheckMinPres( FC_yR.Egy, MinPres );
      FC_zL.Egy = CUFLU_CheckMinPres( FC_zL.Egy, MinPres );
      FC_zR.Egy = CUFLU_CheckMinPres( FC_zR.Egy, MinPres );

#     if ( NCOMP_PASSIVE > 0 )
      for (int v=0; v<NCOMP_PASSIVE; v++)
      {
         FC_xL.Passive[v] = FMAX( FC_xL.Passive[v], TINY_NUMBER );
         FC_xR.Passive[v] = FMAX( FC_xR.Passive[v], TINY_NUMBER );
         FC_yL.Passive[v] = FMAX( FC_yL.Passive[v], TINY_NUMBER );
         FC_yR.Passive[v] = FMAX( FC_yR.Passive[v], TINY_NUMBER );
         FC_zL.Passive[v] = FMAX( FC_zL.Passive[v], TINY_NUMBER );
         FC_zR.Passive[v] = FMAX( FC_zR.Passive[v], TINY_NUMBER );
      }
#     endif

#     undef GetRefState_PLM
#     undef Corr_L_PLM
#     undef Corr_R_PLM
#     undef Corr_Passive
#     undef CTU_Predict_PLM

#     endif // #if ( FLU_SCHEME == CTU )


//    7. primitive variables --> conserved variables
      FC_xL = CUFLU_Pri2Con( FC_xL, _Gamma_m1, NormPassive, NNorm, NormIdx );
      FC_xR = CUFLU_Pri2Con( FC_xR, _Gamma_m1, NormPassive, NNorm, NormIdx );
      FC_yL = CUFLU_Pri2Con( FC_yL, _Gamma_m1, NormPassive, NNorm, NormIdx );
      FC_yR = CUFLU_Pri2Con( FC_yR, _Gamma_m1, NormPassive, NNorm, NormIdx );
      FC_zL = CUFLU_Pri2Con( FC_zL, _Gamma_m1, NormPassive, NNorm, NormIdx );
      FC_zR = CUFLU_Pri2Con( FC_zR, _Gamma_m1, NormPassive, NNorm, NormIdx );


#     if ( FLU_SCHEME == MHM )
//    8. advance the face-centered variables by half time-step for the MHM integrator
      HancockPredict( g_FC_Var_xL, g_FC_Var_xR, g_FC_Var_yL, g_FC_Var_yR, g_FC_Var_zL, g_FC_Var_zR,
                      FC_xL, FC_xR, FC_yL, FC_yR, FC_zL, FC_zR, dt, _dh, Gamma, ID2, CC_C, MinDens, MinPres,
                      NormPassive, NNorm, NormIdx );

#     else  // for MHM_RP and CTU

//    8. store the face-centered values to the output global arrays
#     define Store( g_array, FC )                                    \
      {                                                              \
         g_array[bx][ 0][ID2] = FC.Rho;                              \
         g_array[bx][ 1][ID2] = FC.Px;                               \
         g_array[bx][ 2][ID2] = FC.Py;                               \
         g_array[bx][ 3][ID2] = FC.Pz;                               \
         g_array[bx][ 4][ID2] = FC.Egy;                              \
                                                                     \
         for (int v=0, vv=NCOMP_FLUID; v<NCOMP_PASSIVE; v++, vv++)   \
         g_array[bx][vv][ID2] = FC.Passive[v];                       \
      } // Store

      Store( g_FC_Var_xL, FC_xL );
      Store( g_FC_Var_xR, FC_xR );
      Store( g_FC_Var_yL, FC_yL );
      Store( g_FC_Var_yR, FC_yR );
      Store( g_FC_Var_zL, FC_zL );
      Store( g_FC_Var_zR, FC_zR );

#     undef Store

#     endif // if ( FLU_SCHEME == MHM ) ... else ...

      ID2 += dID2;

   } // while ( ID2 < NOut*NOut*NOut )

} // FUNCTION : CUFLU_DataReconstruction (PLM)
#endif // #if ( LR_SCHEME == PLM )



#if ( LR_SCHEME == PPM )
//-------------------------------------------------------------------------------------------------------
// Function    :  CUFLU_DataReconstruction
// Description :  Reconstruct the face-centered variables by the piecewise-parabolic method (PPM)
//
// Note        :  1. Use the parameter "LR_Limiter" to choose different slope limiters
//                2. The input data should be primitive variables
//                3. For the performance consideration, the output data will be conserved variables
//                4. The face-centered variables will be advanced by half time-step for MHM and CTU schemes
//                5. The PPM reconstruction does NOT support the extrema-preserving limiter
//                6. The function is asynchronous
//                   --> "__syncthreads" must be called before using the output data
//                7. The PLM and PPM data reconstruction functions share the same function name
//                8. The data reconstruction can be applied to characteristic variables by
//                   defining "CHAR_RECONSTRUCTION"
//                9. This function is shared by MHM, MHM_RP, and CTU schemes
//
// Parameter   :  g_PriVar          : Global memory array storing the input primitive variables
//                g_Slope_PPM_x     : Global memory array to store the x-slope for the PPM reconstruction
//                g_Slope_PPM_y     : Global memory array to store the y-slope for the PPM reconstruction
//                g_Slope_PPM_z     : Global memory array to store the z-slope for the PPM reconstruction
//                g_FC_Var_xL       : Global memory array to store the face-centered variables on the -x surface
//                g_FC_Var_xR       : Global memory array to store the face-centered variables on the +x surface
//                g_FC_Var_yL       : Global memory array to store the face-centered variables on the -y surface
//                g_FC_Var_yR       : Global memory array to store the face-centered variables on the +y surface
//                g_FC_Var_zL       : Global memory array to store the face-centered variables on the -z surface
//                g_FC_Var_zR       : Global memory array to store the face-centered variables on the +z surface
//                NIn               : Size of the input array "PriVar" in one direction
//                                    (can be smaller than FLU_NXT)
//                NGhost            : Size of the ghost zone
//                                    --> "NIn-2*NGhost" cells will be computed along each direction
//                                    --> The size of the output array "g_FC_Var" is assumed to be
//                                        "(NIn-2*NGhost)^3"
//                                    --> The reconstructed data at cell (i,j,k) will be stored in the
//                                        array "g_FC_Var" with the index "(i-NGhost,j-NGhost,k-NGhost)
//                Gamma             : Ratio of specific heats
//                LR_Limiter        : Slope limiter for the data reconstruction in the MHM/MHM_RP/CTU schemes
//                                    (0/1/2/3/4) = (vanLeer/generalized MinMod/vanAlbada/
//                                                   vanLeer + generalized MinMod/extrema-preserving) limiter
//                MinMod_Coeff      : Coefficient of the generalized MinMod limiter
//                EP_Coeff          : Coefficient of the extrema-preserving limiter
//                dt                : Time interval to advance solution (useful only for CTU and MHM schemes)
//                _dh               : 1 / grid size (useful only for CTU and MHM schemes)
//                MinDens/Pres      : Minimum allowed density and pressure
//                NormPassive       : true --> convert passive scalars to mass fraction
//                NNorm             : Number of passive scalars for the option "NormPassive"
//                                    --> Should be set to the global variable "PassiveNorm_NVar"
//                NormIdx           : Target variable indices for the option "NormPassive"
//                                    --> Should be set to the global variable "PassiveNorm_VarIdx"
//-------------------------------------------------------------------------------------------------------
__device__ void CUFLU_DataReconstruction( const real g_PriVar[][NCOMP_TOTAL][ FLU_NXT*FLU_NXT*FLU_NXT ],
                                          real g_Slope_PPM_x[][NCOMP_TOTAL][ N_SLOPE_PPM*N_SLOPE_PPM*N_SLOPE_PPM ],
                                          real g_Slope_PPM_y[][NCOMP_TOTAL][ N_SLOPE_PPM*N_SLOPE_PPM*N_SLOPE_PPM ],
                                          real g_Slope_PPM_z[][NCOMP_TOTAL][ N_SLOPE_PPM*N_SLOPE_PPM*N_SLOPE_PPM ],
                                          real g_FC_Var_xL[][NCOMP_TOTAL][ N_FC_VAR*N_FC_VAR*N_FC_VAR ],
                                          real g_FC_Var_xR[][NCOMP_TOTAL][ N_FC_VAR*N_FC_VAR*N_FC_VAR ],
                                          real g_FC_Var_yL[][NCOMP_TOTAL][ N_FC_VAR*N_FC_VAR*N_FC_VAR ],
                                          real g_FC_Var_yR[][NCOMP_TOTAL][ N_FC_VAR*N_FC_VAR*N_FC_VAR ],
                                          real g_FC_Var_zL[][NCOMP_TOTAL][ N_FC_VAR*N_FC_VAR*N_FC_VAR ],
                                          real g_FC_Var_zR[][NCOMP_TOTAL][ N_FC_VAR*N_FC_VAR*N_FC_VAR ],
                                          const int NIn, const int NGhost, const real Gamma,
                                          const LR_Limiter_t LR_Limiter, const real MinMod_Coeff,
                                          const real EP_Coeff, const real dt, const real _dh,
                                          const real MinDens, const real MinPres,
                                          const bool NormPassive, const int NNorm, const int NormIdx[] )
{

   const uint  bx       = blockIdx.x;
   const uint  tx       = threadIdx.x;
   const uint  dID2     = blockDim.x;
   const uint  NOut     = NIn - 2*NGhost;
   const uint  NOut2    = __umul24( NOut, NOut );
   const uint  NSlope   = NOut + 2;
   const uint  NSlope2  = __umul24( NSlope, NSlope );
   const uint3 dr1      = make_uint3( 1, NIn, NIn*NIn );
   const uint3 dr3      = make_uint3( 1, NSlope, NSlope2 );
   const real _Gamma_m1 = (real)1.0 / ( Gamma - (real)1.0 );

   uint   ID2 = tx;
   uint   ID1, ID3, ID1_L, ID1_R, ID3_L, ID3_R;
   uint3  ID1_3d, ID3_3d;
   FluVar Slope_Limiter, dFC, dFC6;
   FluVar CC_L, CC_C, CC_R;                              // cell-centered variables
   FluVar FC_xL, FC_xR, FC_yL, FC_yR, FC_zL, FC_zR;      // face-centered variables
   real   CC1_L, CC1_R;                                  // cell-centered 1-element variables
   real   dCC1_L, dCC1_R, dCC1_C;                        // cell-centered 1-element slopes
   real   Min, Max;

// variables for the CTU scheme
#  if ( FLU_SCHEME == CTU )
   const real dt_dh2 = (real)0.5*dt*_dh;

// include waves both from left and right directions during the data reconstruction, as suggested in ATHENA
#  if (  ( RSOLVER == HLLE  ||  RSOLVER == HLLC )  &&  defined HLL_NO_REF_STATE  )
#  ifdef HLL_INCLUDE_ALL_WAVES
   const bool HLL_Include_All_Waves = true;
#  else
   const bool HLL_Include_All_Waves = false;
#  endif
#  endif // if (  ( RSOLVER == HLLE  ||  RSOLVER == HLLC )  &&  defined HLL_NO_REF_STATE  )

   FluVar5 EVal_x, EVal_y, EVal_z;
   FluVar  Correct_L, Correct_R;
   real    Coeff_L, Coeff_R, Coeff_A, Coeff_B, Coeff_C, Coeff_D;

// initialize the constant components of the matrices of the left and right eigenvectors
   FluVar5 LEVec1 = { 0.0, NULL_REAL, 0.0, 0.0, NULL_REAL };
   FluVar5 LEVec2 = { 1.0,       0.0, 0.0, 0.0, NULL_REAL };
   FluVar5 LEVec3 = { 0.0,       0.0, 1.0, 0.0,       0.0 };
   FluVar5 LEVec4 = { 0.0,       0.0, 0.0, 1.0,       0.0 };
   FluVar5 LEVec5 = { 0.0, NULL_REAL, 0.0, 0.0, NULL_REAL };

   FluVar5 REVec1 = { 1.0, NULL_REAL, 0.0, 0.0, NULL_REAL };
   FluVar5 REVec2 = { 1.0,       0.0, 0.0, 0.0,       0.0 };
   FluVar5 REVec3 = { 0.0,       0.0, 1.0, 0.0,       0.0 };
   FluVar5 REVec4 = { 0.0,       0.0, 0.0, 1.0,       0.0 };
   FluVar5 REVec5 = { 1.0, NULL_REAL, 0.0, 0.0, NULL_REAL };
#  endif // CTU


// 1. evaluate the monotonic slope and store the results to the global arrays
// load over all cells
   while ( ID2 < NSlope*NSlope*NSlope )
   {
      ID1_3d.x = ID2%NSlope         + NGhost-1;
      ID1_3d.y = ID2%NSlope2/NSlope + NGhost-1;
      ID1_3d.z = ID2/NSlope2        + NGhost-1;
      ID1      = __umul24(  __umul24( ID1_3d.z, NIn ) + ID1_3d.y, NIn  ) + ID1_3d.x;

      CC_C.Rho        = g_PriVar[bx][ 0][ID1];
      CC_C.Px         = g_PriVar[bx][ 1][ID1];
      CC_C.Py         = g_PriVar[bx][ 2][ID1];
      CC_C.Pz         = g_PriVar[bx][ 3][ID1];
      CC_C.Egy        = g_PriVar[bx][ 4][ID1];
#     if ( NCOMP_PASSIVE > 0 )
      for (int v=0, vv=NCOMP_FLUID; v<NCOMP_PASSIVE; v++, vv++)
      CC_C.Passive[v] = g_PriVar[bx][vv][ID1];
#     endif


#     define GetSlope_PPM( XYZ, Dir, g_Slope_PPM )                                                    \
      {                                                                                               \
         ID1_L  = ID1 - dr1.XYZ;                                                                      \
         ID1_R  = ID1 + dr1.XYZ;                                                                      \
                                                                                                      \
         /* (1-1) load the cell-centered values */                                                    \
         CC_L.Rho        = g_PriVar[bx][ 0][ID1_L];                                                   \
         CC_L.Px         = g_PriVar[bx][ 1][ID1_L];                                                   \
         CC_L.Py         = g_PriVar[bx][ 2][ID1_L];                                                   \
         CC_L.Pz         = g_PriVar[bx][ 3][ID1_L];                                                   \
         CC_L.Egy        = g_PriVar[bx][ 4][ID1_L];                                                   \
         for (int v=0, vv=NCOMP_FLUID; v<NCOMP_PASSIVE; v++, vv++)                                    \
         CC_L.Passive[v] = g_PriVar[bx][vv][ID1_L];                                                   \
                                                                                                      \
         CC_R.Rho        = g_PriVar[bx][ 0][ID1_R];                                                   \
         CC_R.Px         = g_PriVar[bx][ 1][ID1_R];                                                   \
         CC_R.Py         = g_PriVar[bx][ 2][ID1_R];                                                   \
         CC_R.Pz         = g_PriVar[bx][ 3][ID1_R];                                                   \
         CC_R.Egy        = g_PriVar[bx][ 4][ID1_R];                                                   \
         for (int v=0, vv=NCOMP_FLUID; v<NCOMP_PASSIVE; v++, vv++)                                    \
         CC_R.Passive[v] = g_PriVar[bx][vv][ID1_R];                                                   \
                                                                                                      \
         /* (1-2) evaluate the monotonic slope */                                                     \
         Slope_Limiter = LimitSlope( CC_L, CC_L, CC_C, CC_R, CC_R, LR_Limiter, MinMod_Coeff,          \
                                     NULL_REAL, Gamma, Dir );                                         \
                                                                                                      \
         /* (1-3) store the results to the global arrays */                                           \
         g_Slope_PPM[bx][ 0][ID2] = Slope_Limiter.Rho;                                                \
         g_Slope_PPM[bx][ 1][ID2] = Slope_Limiter.Px;                                                 \
         g_Slope_PPM[bx][ 2][ID2] = Slope_Limiter.Py;                                                 \
         g_Slope_PPM[bx][ 3][ID2] = Slope_Limiter.Pz;                                                 \
         g_Slope_PPM[bx][ 4][ID2] = Slope_Limiter.Egy;                                                \
         for (int v=0, vv=NCOMP_FLUID; v<NCOMP_PASSIVE; v++, vv++)                                    \
         g_Slope_PPM[bx][vv][ID2] = Slope_Limiter.Passive[v];                                         \
      }

      GetSlope_PPM( x, 0, g_Slope_PPM_x );
      GetSlope_PPM( y, 1, g_Slope_PPM_y );
      GetSlope_PPM( z, 2, g_Slope_PPM_z );

#     undef GetSlope_PPM

      ID2 += dID2;

   } // while ( ID2 < NSlope*NSlope*NSlope )

   __syncthreads();



// 2. evaluate the face-centered values
   ID2 = tx;

// load over all cells
   while ( ID2 < NOut*NOut*NOut )
   {
      ID1_3d.x = ID2%NOut       + NGhost;
      ID1_3d.y = ID2%NOut2/NOut + NGhost;
      ID1_3d.z = ID2/NOut2      + NGhost;
      ID1      = __umul24(  __umul24( ID1_3d.z, NIn    ) + ID1_3d.y, NIn     ) + ID1_3d.x;

      ID3_3d.x = ID2%NOut       + 1;
      ID3_3d.y = ID2%NOut2/NOut + 1;
      ID3_3d.z = ID2/NOut2      + 1;
      ID3      = __umul24(  __umul24( ID3_3d.z, NSlope ) + ID3_3d.y, NSlope  ) + ID3_3d.x;

      CC_C.Rho        = g_PriVar[bx][ 0][ID1];
      CC_C.Px         = g_PriVar[bx][ 1][ID1];
      CC_C.Py         = g_PriVar[bx][ 2][ID1];
      CC_C.Pz         = g_PriVar[bx][ 3][ID1];
      CC_C.Egy        = g_PriVar[bx][ 4][ID1];
#     if ( NCOMP_PASSIVE > 0 )
      for (int v=0, vv=NCOMP_FLUID; v<NCOMP_PASSIVE; v++, vv++)
      CC_C.Passive[v] = g_PriVar[bx][vv][ID1];
#     endif


//    (2-1) evaluate the eigenvalues and eigenvectors for the CTU integrator
#     if ( FLU_SCHEME == CTU )
      Get_EigenSystem( CC_C, EVal_x, EVal_y, EVal_z, LEVec1, LEVec2, LEVec3, LEVec4, LEVec5,
                       REVec1, REVec2, REVec3, REVec4, REVec5, Gamma );
#     endif


//    (2-2) get the face-centered primitive variables
#     define GetFCVar_PPM( XYZ, comp, v, CC1_C, FC1_L, FC1_R, g_Slope_PPM )                           \
      {                                                                                               \
         /* parabolic interpolation */                                                                \
         ID1_L  = ID1 - dr1.XYZ;                                                                      \
         ID1_R  = ID1 + dr1.XYZ;                                                                      \
         ID3_L  = ID3 - dr3.XYZ;                                                                      \
         ID3_R  = ID3 + dr3.XYZ;                                                                      \
                                                                                                      \
         CC1_L  = g_PriVar[bx][v][ID1_L];                                                             \
         CC1_R  = g_PriVar[bx][v][ID1_R];                                                             \
                                                                                                      \
         dCC1_L = g_Slope_PPM[bx][v][ID3_L];                                                          \
         dCC1_C = g_Slope_PPM[bx][v][ID3  ];                                                          \
         dCC1_R = g_Slope_PPM[bx][v][ID3_R];                                                          \
                                                                                                      \
         FC1_L = (real)0.5*( CC1_C + CC1_L ) - (real)1.0/(real)6.0*( dCC1_C - dCC1_L );               \
         FC1_R = (real)0.5*( CC1_C + CC1_R ) - (real)1.0/(real)6.0*( dCC1_R - dCC1_C );               \
                                                                                                      \
                                                                                                      \
         /* monotonicity constraint */                                                                \
         dFC.comp  = FC1_R - FC1_L;                                                                   \
         dFC6.comp = (real)6.0*(  CC1_C - (real)0.5*( FC1_L + FC1_R )  );                             \
                                                                                                      \
         if (  ( FC1_R - CC1_C )*( CC1_C - FC1_L ) <= (real)0.0  )                                    \
         {                                                                                            \
            FC1_L = CC1_C;                                                                            \
            FC1_R = CC1_C;                                                                            \
         }                                                                                            \
         else if ( dFC.comp*dFC6.comp > +dFC.comp*dFC.comp )                                          \
            FC1_L = (real)3.0*CC1_C - (real)2.0*FC1_R;                                                \
         else if ( dFC.comp*dFC6.comp < -dFC.comp*dFC.comp )                                          \
            FC1_R = (real)3.0*CC1_C - (real)2.0*FC1_L;                                                \
                                                                                                      \
                                                                                                      \
         /* ensure the face-centered variables lie between neighboring cell-centered values */        \
         Min   = ( CC1_C < CC1_L ) ? CC1_C : CC1_L;                                                   \
         Max   = ( CC1_C > CC1_L ) ? CC1_C : CC1_L;                                                   \
         FC1_L = ( FC1_L > Min   ) ? FC1_L : Min;                                                     \
         FC1_L = ( FC1_L < Max   ) ? FC1_L : Max;                                                     \
                                                                                                      \
         Min   = ( CC1_C < CC1_R ) ? CC1_C : CC1_R;                                                   \
         Max   = ( CC1_C > CC1_R ) ? CC1_C : CC1_R;                                                   \
         FC1_R = ( FC1_R > Min   ) ? FC1_R : Min;                                                     \
         FC1_R = ( FC1_R < Max   ) ? FC1_R : Max;                                                     \
                                                                                                      \
      } // GetFCVar_PPM

      GetFCVar_PPM( x, Rho,         0, CC_C.Rho,        FC_xL.Rho,        FC_xR.Rho,        g_Slope_PPM_x );
      GetFCVar_PPM( x, Px,          1, CC_C.Px,         FC_xL.Px,         FC_xR.Px,         g_Slope_PPM_x );
      GetFCVar_PPM( x, Py,          2, CC_C.Py,         FC_xL.Py,         FC_xR.Py,         g_Slope_PPM_x );
      GetFCVar_PPM( x, Pz,          3, CC_C.Pz,         FC_xL.Pz,         FC_xR.Pz,         g_Slope_PPM_x );
      GetFCVar_PPM( x, Egy,         4, CC_C.Egy,        FC_xL.Egy,        FC_xR.Egy,        g_Slope_PPM_x );
#     if ( NCOMP_PASSIVE > 0 )
      for (int v=0, vv=NCOMP_FLUID; v<NCOMP_PASSIVE; v++, vv++)
      GetFCVar_PPM( x, Passive[v], vv, CC_C.Passive[v], FC_xL.Passive[v], FC_xR.Passive[v], g_Slope_PPM_x );
#     endif

      GetFCVar_PPM( y, Rho,         0, CC_C.Rho,        FC_yL.Rho,        FC_yR.Rho,        g_Slope_PPM_y );
      GetFCVar_PPM( y, Px,          1, CC_C.Px,         FC_yL.Px,         FC_yR.Px,         g_Slope_PPM_y );
      GetFCVar_PPM( y, Py,          2, CC_C.Py,         FC_yL.Py,         FC_yR.Py,         g_Slope_PPM_y );
      GetFCVar_PPM( y, Pz,          3, CC_C.Pz,         FC_yL.Pz,         FC_yR.Pz,         g_Slope_PPM_y );
      GetFCVar_PPM( y, Egy,         4, CC_C.Egy,        FC_yL.Egy,        FC_yR.Egy,        g_Slope_PPM_y );
#     if ( NCOMP_PASSIVE > 0 )
      for (int v=0, vv=NCOMP_FLUID; v<NCOMP_PASSIVE; v++, vv++)
      GetFCVar_PPM( y, Passive[v], vv, CC_C.Passive[v], FC_yL.Passive[v], FC_yR.Passive[v], g_Slope_PPM_y );
#     endif

      GetFCVar_PPM( z, Rho,         0, CC_C.Rho,        FC_zL.Rho,        FC_zR.Rho,        g_Slope_PPM_z );
      GetFCVar_PPM( z, Px,          1, CC_C.Px,         FC_zL.Px,         FC_zR.Px,         g_Slope_PPM_z );
      GetFCVar_PPM( z, Py,          2, CC_C.Py,         FC_zL.Py,         FC_zR.Py,         g_Slope_PPM_z );
      GetFCVar_PPM( z, Pz,          3, CC_C.Pz,         FC_zL.Pz,         FC_zR.Pz,         g_Slope_PPM_z );
      GetFCVar_PPM( z, Egy,         4, CC_C.Egy,        FC_zL.Egy,        FC_zR.Egy,        g_Slope_PPM_z );
#     if ( NCOMP_PASSIVE > 0 )
      for (int v=0, vv=NCOMP_FLUID; v<NCOMP_PASSIVE; v++, vv++)
      GetFCVar_PPM( z, Passive[v], vv, CC_C.Passive[v], FC_zL.Passive[v], FC_zR.Passive[v], g_Slope_PPM_z );
#     endif

#     undef GetFCVar_PPM


//    (2-3) advance the face-centered variables by half time-step for the CTU integrator
#     if ( FLU_SCHEME == CTU )

//    =====================================================================================
//    a. for the HLL solvers (HLLE/HLLC)
//    =====================================================================================
#     if (  ( RSOLVER == HLLE  ||  RSOLVER == HLLC )  &&  defined HLL_NO_REF_STATE  )

#     define GetRefState_PPM( EVal )                                                                  \
      {                                                                                               \
         Correct_L.Rho = (real)0.0;                                                                   \
         Correct_L.Px  = (real)0.0;                                                                   \
         Correct_L.Py  = (real)0.0;                                                                   \
         Correct_L.Pz  = (real)0.0;                                                                   \
         Correct_L.Egy = (real)0.0;                                                                   \
                                                                                                      \
         Correct_R.Rho = (real)0.0;                                                                   \
         Correct_R.Px  = (real)0.0;                                                                   \
         Correct_R.Py  = (real)0.0;                                                                   \
         Correct_R.Pz  = (real)0.0;                                                                   \
         Correct_R.Egy = (real)0.0;                                                                   \
      } // GetRefState_PPM

#     define Corr_L_PPM( comp, EVal, LEVec, REVec )                                                   \
      {                                                                                               \
         Coeff_L = (real)0.0;                                                                         \
                                                                                                      \
         if ( HLL_Include_All_Waves  ||  EVal.comp <= (real)0.0 )                                     \
         {                                                                                            \
            Coeff_C = -EVal.comp;                                                                     \
            Coeff_D = real(-4.0/3.0)*dt_dh2*SQR(Coeff_C);                                             \
                                                                                                      \
            Coeff_L += LEVec.Rho*(  Coeff_C*( dFC.Rho + dFC6.Rho ) + Coeff_D*( dFC6.Rho )  );         \
            Coeff_L += LEVec.Px *(  Coeff_C*( dFC.Px  + dFC6.Px  ) + Coeff_D*( dFC6.Px  )  );         \
            Coeff_L += LEVec.Py *(  Coeff_C*( dFC.Py  + dFC6.Py  ) + Coeff_D*( dFC6.Py  )  );         \
            Coeff_L += LEVec.Pz *(  Coeff_C*( dFC.Pz  + dFC6.Pz  ) + Coeff_D*( dFC6.Pz  )  );         \
            Coeff_L += LEVec.Egy*(  Coeff_C*( dFC.Egy + dFC6.Egy ) + Coeff_D*( dFC6.Egy )  );         \
                                                                                                      \
            Coeff_L *= dt_dh2;                                                                        \
                                                                                                      \
            Correct_L.Rho += Coeff_L*REVec.Rho;                                                       \
            Correct_L.Px  += Coeff_L*REVec.Px;                                                        \
            Correct_L.Py  += Coeff_L*REVec.Py;                                                        \
            Correct_L.Pz  += Coeff_L*REVec.Pz;                                                        \
            Correct_L.Egy += Coeff_L*REVec.Egy;                                                       \
         }                                                                                            \
      } // Corr_L_PPM

#     define Corr_R_PPM( comp, EVal, LEVec, REVec )                                                   \
      {                                                                                               \
         Coeff_R = (real)0.0;                                                                         \
                                                                                                      \
         if ( HLL_Include_All_Waves  ||  EVal.comp >= (real)0.0 )                                     \
         {                                                                                            \
            Coeff_A = -EVal.comp;                                                                     \
            Coeff_B = real(-4.0/3.0)*dt_dh2*SQR(Coeff_A);                                             \
                                                                                                      \
            Coeff_R += LEVec.Rho*(  Coeff_A*( dFC.Rho - dFC6.Rho ) + Coeff_B*( dFC6.Rho )  );         \
            Coeff_R += LEVec.Px *(  Coeff_A*( dFC.Px  - dFC6.Px  ) + Coeff_B*( dFC6.Px  )  );         \
            Coeff_R += LEVec.Py *(  Coeff_A*( dFC.Py  - dFC6.Py  ) + Coeff_B*( dFC6.Py  )  );         \
            Coeff_R += LEVec.Pz *(  Coeff_A*( dFC.Pz  - dFC6.Pz  ) + Coeff_B*( dFC6.Pz  )  );         \
            Coeff_R += LEVec.Egy*(  Coeff_A*( dFC.Egy - dFC6.Egy ) + Coeff_B*( dFC6.Egy )  );         \
                                                                                                      \
            Coeff_R *= dt_dh2;                                                                        \
                                                                                                      \
            Correct_R.Rho += Coeff_R*REVec.Rho;                                                       \
            Correct_R.Px  += Coeff_R*REVec.Px;                                                        \
            Correct_R.Py  += Coeff_R*REVec.Py;                                                        \
            Correct_R.Pz  += Coeff_R*REVec.Pz;                                                        \
            Correct_R.Egy += Coeff_R*REVec.Egy;                                                       \
         }                                                                                            \
      } // Corr_R_PPM


//    =====================================================================================
//    b. for the Roe's and exact solvers
//    =====================================================================================
#     else // ( RSOLVER == ROE/EXACT && ifndef HLL_NO_REF_STATE )

#     define GetRefState_PPM( EVal )                                                                  \
      {                                                                                               \
         Coeff_L = -dt_dh2*FMIN( EVal.Rho, (real)0.0 );                                               \
         Coeff_R = -dt_dh2*FMAX( EVal.Egy, (real)0.0 );                                               \
                                                                                                      \
         Correct_L.Rho = Coeff_L*(  dFC.Rho + ( (real)1.0 - real(4.0/3.0)*Coeff_L )*dFC6.Rho  );      \
         Correct_L.Px  = Coeff_L*(  dFC.Px  + ( (real)1.0 - real(4.0/3.0)*Coeff_L )*dFC6.Px   );      \
         Correct_L.Py  = Coeff_L*(  dFC.Py  + ( (real)1.0 - real(4.0/3.0)*Coeff_L )*dFC6.Py   );      \
         Correct_L.Pz  = Coeff_L*(  dFC.Pz  + ( (real)1.0 - real(4.0/3.0)*Coeff_L )*dFC6.Pz   );      \
         Correct_L.Egy = Coeff_L*(  dFC.Egy + ( (real)1.0 - real(4.0/3.0)*Coeff_L )*dFC6.Egy  );      \
                                                                                                      \
         Correct_R.Rho = Coeff_R*(  dFC.Rho - ( (real)1.0 + real(4.0/3.0)*Coeff_R )*dFC6.Rho  );      \
         Correct_R.Px  = Coeff_R*(  dFC.Px  - ( (real)1.0 + real(4.0/3.0)*Coeff_R )*dFC6.Px   );      \
         Correct_R.Py  = Coeff_R*(  dFC.Py  - ( (real)1.0 + real(4.0/3.0)*Coeff_R )*dFC6.Py   );      \
         Correct_R.Pz  = Coeff_R*(  dFC.Pz  - ( (real)1.0 + real(4.0/3.0)*Coeff_R )*dFC6.Pz   );      \
         Correct_R.Egy = Coeff_R*(  dFC.Egy - ( (real)1.0 + real(4.0/3.0)*Coeff_R )*dFC6.Egy  );      \
      } // GetRefState_PPM

#     define Corr_L_PPM( comp, EVal, LEVec, REVec )                                                   \
      {                                                                                               \
         Coeff_L = (real)0.0;                                                                         \
                                                                                                      \
         if ( EVal.comp <= (real)0.0 )                                                                \
         {                                                                                            \
            Coeff_C = EVal.Rho - EVal.comp;                                                           \
            /* write as (a-b)*(a+b) instead of a^2-b^2 to ensure that Coeff_D=0 when Coeff_C=0 */     \
            /* Coeff_D = real(4.0/3.0)*dt_dh2*( EVal.Rho*EVal.Rho - EVal.comp*EVal.comp ); */         \
            Coeff_D = real(4.0/3.0)*dt_dh2*Coeff_C*( EVal.Rho + EVal.comp );                          \
                                                                                                      \
            Coeff_L += LEVec.Rho*(  Coeff_C*( dFC.Rho + dFC6.Rho ) + Coeff_D*( dFC6.Rho )  );         \
            Coeff_L += LEVec.Px *(  Coeff_C*( dFC.Px  + dFC6.Px  ) + Coeff_D*( dFC6.Px  )  );         \
            Coeff_L += LEVec.Py *(  Coeff_C*( dFC.Py  + dFC6.Py  ) + Coeff_D*( dFC6.Py  )  );         \
            Coeff_L += LEVec.Pz *(  Coeff_C*( dFC.Pz  + dFC6.Pz  ) + Coeff_D*( dFC6.Pz  )  );         \
            Coeff_L += LEVec.Egy*(  Coeff_C*( dFC.Egy + dFC6.Egy ) + Coeff_D*( dFC6.Egy )  );         \
                                                                                                      \
            Coeff_L *= dt_dh2;                                                                        \
                                                                                                      \
            Correct_L.Rho += Coeff_L*REVec.Rho;                                                       \
            Correct_L.Px  += Coeff_L*REVec.Px;                                                        \
            Correct_L.Py  += Coeff_L*REVec.Py;                                                        \
            Correct_L.Pz  += Coeff_L*REVec.Pz;                                                        \
            Correct_L.Egy += Coeff_L*REVec.Egy;                                                       \
         }                                                                                            \
      } // Corr_L_PPM

#     define Corr_R_PPM( comp, EVal, LEVec, REVec )                                                   \
      {                                                                                               \
         Coeff_R = (real)0.0;                                                                         \
                                                                                                      \
         if ( EVal.comp >= (real)0.0 )                                                                \
         {                                                                                            \
            Coeff_A = EVal.Egy - EVal.comp;                                                           \
            /* write as (a-b)*(a+b) instead of a^2-b^2 to ensure that Coeff_B=0 when Coeff_A=0 */     \
            /* Coeff_B = real(4.0/3.0)*dt_dh2*( EVal.Egy*EVal.Egy - EVal.comp*EVal.comp ); */         \
            Coeff_B = real(4.0/3.0)*dt_dh2*Coeff_A*( EVal.Egy + EVal.comp );                          \
                                                                                                      \
            Coeff_R += LEVec.Rho*(  Coeff_A*( dFC.Rho - dFC6.Rho ) + Coeff_B*( dFC6.Rho )  );         \
            Coeff_R += LEVec.Px *(  Coeff_A*( dFC.Px  - dFC6.Px  ) + Coeff_B*( dFC6.Px  )  );         \
            Coeff_R += LEVec.Py *(  Coeff_A*( dFC.Py  - dFC6.Py  ) + Coeff_B*( dFC6.Py  )  );         \
            Coeff_R += LEVec.Pz *(  Coeff_A*( dFC.Pz  - dFC6.Pz  ) + Coeff_B*( dFC6.Pz  )  );         \
            Coeff_R += LEVec.Egy*(  Coeff_A*( dFC.Egy - dFC6.Egy ) + Coeff_B*( dFC6.Egy )  );         \
                                                                                                      \
            Coeff_R *= dt_dh2;                                                                        \
                                                                                                      \
            Correct_R.Rho += Coeff_R*REVec.Rho;                                                       \
            Correct_R.Px  += Coeff_R*REVec.Px;                                                        \
            Correct_R.Py  += Coeff_R*REVec.Py;                                                        \
            Correct_R.Pz  += Coeff_R*REVec.Pz;                                                        \
            Correct_R.Egy += Coeff_R*REVec.Egy;                                                       \
         }                                                                                            \
      } // Corr_R_PPM


#     endif // if (  ( RSOLVER == HLLE  ||  RSOLVER == HLLC )  &&  defined HLL_NO_REF_STATE  ) ... else ...


//    macro for correcting the passive scalars
#     if ( NCOMP_PASSIVE > 0 )

#     define Corr_Passive( EVal )                                                                     \
      {                                                                                               \
         Coeff_L = -dt_dh2*FMIN( EVal.Px, (real)0.0 );                                                \
         Coeff_R = -dt_dh2*FMAX( EVal.Px, (real)0.0 );                                                \
                                                                                                      \
         for (int v=0; v<NCOMP_PASSIVE; v++)                                                          \
         {                                                                                            \
            Correct_L.Passive[v]                                                                      \
               = Coeff_L*(  dFC.Passive[v] + ( (real)1.0 - real(4.0/3.0)*Coeff_L )*dFC6.Passive[v]  );\
            Correct_R.Passive[v]                                                                      \
               = Coeff_R*(  dFC.Passive[v] - ( (real)1.0 + real(4.0/3.0)*Coeff_R )*dFC6.Passive[v]  );\
         }                                                                                            \
      }                                                                                               \

#     else

#     define Corr_Passive( EVal )   /* nothing to do */

#     endif // #if ( NCOMP_PASSIVE > 0 ) ... else ...


#     define CTU_Predict_PPM( dir, FC_L, FC_R, EVal )                                                          \
      {                                                                                                        \
         /* (2-3-1) compute the PPM coefficient */                                                             \
         dFC.Rho         = FC_R.Rho        - FC_L.Rho;                                                         \
         dFC.Px          = FC_R.Px         - FC_L.Px;                                                          \
         dFC.Py          = FC_R.Py         - FC_L.Py;                                                          \
         dFC.Pz          = FC_R.Pz         - FC_L.Pz;                                                          \
         dFC.Egy         = FC_R.Egy        - FC_L.Egy;                                                         \
         for (int v=0; v<NCOMP_PASSIVE; v++)                                                                   \
         dFC.Passive[v]  = FC_R.Passive[v] - FC_L.Passive[v];                                                  \
                                                                                                               \
         dFC6.Rho        = (real)6.0*(  CC_C.Rho        - (real)0.5*( FC_L.Rho        + FC_R.Rho        )  );  \
         dFC6.Px         = (real)6.0*(  CC_C.Px         - (real)0.5*( FC_L.Px         + FC_R.Px         )  );  \
         dFC6.Py         = (real)6.0*(  CC_C.Py         - (real)0.5*( FC_L.Py         + FC_R.Py         )  );  \
         dFC6.Pz         = (real)6.0*(  CC_C.Pz         - (real)0.5*( FC_L.Pz         + FC_R.Pz         )  );  \
         dFC6.Egy        = (real)6.0*(  CC_C.Egy        - (real)0.5*( FC_L.Egy        + FC_R.Egy        )  );  \
         for (int v=0; v<NCOMP_PASSIVE; v++)                                                                   \
         dFC6.Passive[v] = (real)6.0*(  CC_C.Passive[v] - (real)0.5*( FC_L.Passive[v] + FC_R.Passive[v] )  );  \
                                                                                                               \
                                                                                                               \
         /* (2-3-2) re-order variables for the y/z directions */                                               \
         dFC  = CUFLU_Rotate3D( dFC,  dir, true );                                                             \
         dFC6 = CUFLU_Rotate3D( dFC6, dir, true );                                                             \
                                                                                                               \
                                                                                                               \
         /* (2-3-3) evaluate the reference states */                                                           \
         GetRefState_PPM( EVal );                                                                              \
                                                                                                               \
                                                                                                               \
         /* (2-3-4) evaluate the corrections to the left and right face-centered active variables */           \
         Corr_L_PPM( Rho, EVal, LEVec1, REVec1 );                                                              \
         Corr_L_PPM( Px,  EVal, LEVec2, REVec2 );                                                              \
         Corr_L_PPM( Py,  EVal, LEVec3, REVec3 );                                                              \
         Corr_L_PPM( Pz,  EVal, LEVec4, REVec4 );                                                              \
         Corr_L_PPM( Egy, EVal, LEVec5, REVec5 );                                                              \
                                                                                                               \
         Corr_R_PPM( Rho, EVal, LEVec1, REVec1 );                                                              \
         Corr_R_PPM( Px,  EVal, LEVec2, REVec2 );                                                              \
         Corr_R_PPM( Py,  EVal, LEVec3, REVec3 );                                                              \
         Corr_R_PPM( Pz,  EVal, LEVec4, REVec4 );                                                              \
         Corr_R_PPM( Egy, EVal, LEVec5, REVec5 );                                                              \
                                                                                                               \
                                                                                                               \
         /* (2-3-5) evaluate the corrections to the left and right face-centered passive scalars */            \
         /*         --> passive scalars travel with fluid velocity (i.e., entropy mode) */                     \
         Corr_Passive( EVal );                                                                                 \
                                                                                                               \
                                                                                                               \
         /* (2-3-6) evolve the face-centered variables by half time-step */                                    \
         Correct_L = CUFLU_Rotate3D( Correct_L, dir, false );                                                  \
         Correct_R = CUFLU_Rotate3D( Correct_R, dir, false );                                                  \
                                                                                                               \
         FC_L.Rho        += Correct_L.Rho;                                                                     \
         FC_L.Px         += Correct_L.Px;                                                                      \
         FC_L.Py         += Correct_L.Py;                                                                      \
         FC_L.Pz         += Correct_L.Pz;                                                                      \
         FC_L.Egy        += Correct_L.Egy;                                                                     \
         for (int v=0; v<NCOMP_PASSIVE; v++)                                                                   \
         FC_L.Passive[v] += Correct_L.Passive[v];                                                              \
                                                                                                               \
         FC_R.Rho        += Correct_R.Rho;                                                                     \
         FC_R.Px         += Correct_R.Px;                                                                      \
         FC_R.Py         += Correct_R.Py;                                                                      \
         FC_R.Pz         += Correct_R.Pz;                                                                      \
         FC_R.Egy        += Correct_R.Egy;                                                                     \
         for (int v=0; v<NCOMP_PASSIVE; v++)                                                                   \
         FC_R.Passive[v] += Correct_R.Passive[v];                                                              \
      } // CTU_Predict_PPM

      CTU_Predict_PPM( 0, FC_xL, FC_xR, EVal_x );
      CTU_Predict_PPM( 1, FC_yL, FC_yR, EVal_y );
      CTU_Predict_PPM( 2, FC_zL, FC_zR, EVal_z );

//    (2-3-5) ensure positive density and pressure
      FC_xL.Rho = FMAX( FC_xL.Rho, MinDens );
      FC_xR.Rho = FMAX( FC_xR.Rho, MinDens );
      FC_yL.Rho = FMAX( FC_yL.Rho, MinDens );
      FC_yR.Rho = FMAX( FC_yR.Rho, MinDens );
      FC_zL.Rho = FMAX( FC_zL.Rho, MinDens );
      FC_zR.Rho = FMAX( FC_zR.Rho, MinDens );

      FC_xL.Egy = CUFLU_CheckMinPres( FC_xL.Egy, MinPres );
      FC_xR.Egy = CUFLU_CheckMinPres( FC_xR.Egy, MinPres );
      FC_yL.Egy = CUFLU_CheckMinPres( FC_yL.Egy, MinPres );
      FC_yR.Egy = CUFLU_CheckMinPres( FC_yR.Egy, MinPres );
      FC_zL.Egy = CUFLU_CheckMinPres( FC_zL.Egy, MinPres );
      FC_zR.Egy = CUFLU_CheckMinPres( FC_zR.Egy, MinPres );

#     if ( NCOMP_PASSIVE > 0 )
      for (int v=0; v<NCOMP_PASSIVE; v++)
      {
         FC_xL.Passive[v] = FMAX( FC_xL.Passive[v], TINY_NUMBER );
         FC_xR.Passive[v] = FMAX( FC_xR.Passive[v], TINY_NUMBER );
         FC_yL.Passive[v] = FMAX( FC_yL.Passive[v], TINY_NUMBER );
         FC_yR.Passive[v] = FMAX( FC_yR.Passive[v], TINY_NUMBER );
         FC_zL.Passive[v] = FMAX( FC_zL.Passive[v], TINY_NUMBER );
         FC_zR.Passive[v] = FMAX( FC_zR.Passive[v], TINY_NUMBER );
      }
#     endif

#     undef GetRefState_PPM
#     undef Corr_L_PPM
#     undef Corr_R_PPM
#     undef Corr_Passive
#     undef CTU_Predict_PPM

#     endif // #if ( FLU_SCHEME == CTU )


//    (2-4) primitive variables --> conserved variables
      FC_xL = CUFLU_Pri2Con( FC_xL, _Gamma_m1, NormPassive, NNorm, NormIdx );
      FC_xR = CUFLU_Pri2Con( FC_xR, _Gamma_m1, NormPassive, NNorm, NormIdx );
      FC_yL = CUFLU_Pri2Con( FC_yL, _Gamma_m1, NormPassive, NNorm, NormIdx );
      FC_yR = CUFLU_Pri2Con( FC_yR, _Gamma_m1, NormPassive, NNorm, NormIdx );
      FC_zL = CUFLU_Pri2Con( FC_zL, _Gamma_m1, NormPassive, NNorm, NormIdx );
      FC_zR = CUFLU_Pri2Con( FC_zR, _Gamma_m1, NormPassive, NNorm, NormIdx );


#     if ( FLU_SCHEME == MHM )
//    (2-5) advance the face-centered variables by half time-step for the MHM scheme
      HancockPredict( g_FC_Var_xL, g_FC_Var_xR, g_FC_Var_yL, g_FC_Var_yR, g_FC_Var_zL, g_FC_Var_zR,
                      FC_xL, FC_xR, FC_yL, FC_yR, FC_zL, FC_zR, dt, _dh, Gamma, ID2, CC_C, MinDens, MinPres,
                      NormPassive, NNorm, NormIdx );

#     else  // for MHM_RP and CTU

//    (2-5) store the face-centered values to the output global arrays
#     define Store( g_array, FC )                                    \
      {                                                              \
         g_array[bx][ 0][ID2] = FC.Rho;                              \
         g_array[bx][ 1][ID2] = FC.Px;                               \
         g_array[bx][ 2][ID2] = FC.Py;                               \
         g_array[bx][ 3][ID2] = FC.Pz;                               \
         g_array[bx][ 4][ID2] = FC.Egy;                              \
                                                                     \
         for (int v=0, vv=NCOMP_FLUID; v<NCOMP_PASSIVE; v++, vv++)   \
         g_array[bx][vv][ID2] = FC.Passive[v];                       \
      } // Store

      Store( g_FC_Var_xL, FC_xL );
      Store( g_FC_Var_xR, FC_xR );
      Store( g_FC_Var_yL, FC_yL );
      Store( g_FC_Var_yR, FC_yR );
      Store( g_FC_Var_zL, FC_zL );
      Store( g_FC_Var_zR, FC_zR );

#     undef Store

#     endif // if ( FLU_SCHEME == MHM ) ... else ...

      ID2 += dID2;

   } // while ( ID2 < NOut*NOut*NOut )

} // FUNCTION : CUFLU_DataReconstruction (PPM)
#endif // #if ( LR_SCHEME == PPM )



#ifdef CHAR_RECONSTRUCTION
//-------------------------------------------------------------------------------------------------------
// Function    :  Pri2Char
// Description :  Convert the primitive variables to the characteristic variables
//
// Note           1. Passive scalars require no conversion
//                   --> Their eigenmatrices are just identity matrix
//
// Parameter   :  InOut : Array storing both the input primitive variables and output characteristic variables
//                Gamma : Ratio of specific heats
//                Rho   : Density
//                Pres  : Pressure
//                XYZ   : Target spatial direction : (0/1/2) --> (x/y/z)
//-------------------------------------------------------------------------------------------------------
__device__ FluVar Pri2Char( FluVar InOut, const real Gamma, const real Rho, const real Pres, const int XYZ )
{

#  ifdef CHECK_NEGATIVE_IN_FLUID
   if ( CUFLU_CheckNegative(Pres) )
      printf( "ERROR : negative pressure (%14.7e) at file <%s>, line <%d>, function <%s>\n",
              Pres, __FILE__, __LINE__, __FUNCTION__ );

   if ( CUFLU_CheckNegative(Rho) )
      printf( "ERROR : negative density (%14.7e) at file <%s>, line <%d>, function <%s>\n",
              Rho,  __FILE__, __LINE__, __FUNCTION__ );
#  endif

   const real _Cs2 = (real)1.0 / ( Gamma*Pres/Rho );
   const real _Cs  = SQRT( _Cs2 );

   FluVar Temp = CUFLU_Rotate3D( InOut, XYZ, true );

   InOut.Rho = -(real)0.5*Rho*_Cs*Temp.Px + (real)0.5*_Cs2*Temp.Egy;
   InOut.Px  = Temp.Rho - _Cs2*Temp.Egy;
   InOut.Py  = Temp.Py;
   InOut.Pz  = Temp.Pz;
   InOut.Egy = +(real)0.5*Rho*_Cs*Temp.Px + (real)0.5*_Cs2*Temp.Egy;

   return InOut;

} // FUNCTION : Pri2Char



//-------------------------------------------------------------------------------------------------------
// Function    :  Char2Pri
// Description :  Convert the characteristic variables to the primitive variables
//
// Note           1. Passive scalars require no conversion
//                   --> Their eigenmatrices are just identity matrix
//
// Parameter   :  InOut : Array storing both the input characteristic variables and output primitive variables
//                Gamma : Ratio of specific heats
//                Rho   : Density
//                Pres  : Pressure
//                XYZ   : Target spatial direction : (0/1/2) --> (x/y/z)
//-------------------------------------------------------------------------------------------------------
__device__ FluVar Char2Pri( FluVar InOut, const real Gamma, const real Rho, const real Pres, const int XYZ )
{

#  ifdef CHECK_NEGATIVE_IN_FLUID
   if ( CUFLU_CheckNegative(Pres) )
      printf( "ERROR : negative pressure (%14.7e) at file <%s>, line <%d>, function <%s>\n",
              Pres, __FILE__, __LINE__, __FUNCTION__ );

   if ( CUFLU_CheckNegative(Rho) )
      printf( "ERROR : negative density (%14.7e) at file <%s>, line <%d>, function <%s>\n",
              Rho,  __FILE__, __LINE__, __FUNCTION__ );
#  endif

   const real _Rho = (real)1.0 / Rho;
   const real Cs2  = Gamma*Pres*_Rho;
   const real Cs   = SQRT( Cs2 );
   FluVar Temp;

   Temp.Rho = InOut.Rho + InOut.Px + InOut.Egy;
   Temp.Px  = Cs*_Rho*( -InOut.Rho + InOut.Egy );
   Temp.Py  = InOut.Py;
   Temp.Pz  = InOut.Pz;
   Temp.Egy = Cs2*( InOut.Rho + InOut.Egy );

#  if ( NCOMP_PASSIVE > 0 )
   for (int v=0; v<NCOMP_PASSIVE; v++)    Temp.Passive[v] = InOut.Passive[v];
#  endif

   InOut = CUFLU_Rotate3D( Temp, XYZ, false );

   return InOut;

} // FUNCTION : Char2Pri
#endif // #ifdef CHAR_RECONSTRUCTION



//-------------------------------------------------------------------------------------------------------
// Function    :  LimitSlope
// Description :  Evaluate the monotonic slope by applying slope limiters
//
// Note        :  1. The input data should be primitive variables
//                2. The parameters "L2, R2, EP_Coeff" are useful only for the extrema-preserving limiter
//
// Parameter   :  L2             : Element x-2
//                L1             : Element x-1
//                C0             : Element x
//                R1             : Element x+1
//                R2             : Element x+2
//                LR_Limiter     : Slope limiter for the data reconstruction in the MHM/MHM_RP/CTU schemes
//                                 (0/1/2/3/4) = (vanLeer/generalized MinMod/vanAlbada/
//                                                vanLeer + generalized MinMod/extrema-preserving) limiter
//                MinMod_Coeff   : Coefficient of the generalized MinMod limiter
//                EP_Coeff       : Coefficient of the extrema-preserving limiter
//                Gamma          : Ratio of specific heats
//                                 --> Useful only if the option "CHAR_RECONSTRUCTION" is turned on
//                XYZ            : Target spatial direction : (0/1/2) --> (x/y/z)
//                                 --> Useful only if the option "CHAR_RECONSTRUCTION" is turned on
//-------------------------------------------------------------------------------------------------------
__device__ FluVar LimitSlope( const FluVar L2, const FluVar L1, const FluVar C0, const FluVar R1, const FluVar R2,
                              const LR_Limiter_t LR_Limiter, const real MinMod_Coeff, const real EP_Coeff,
                              const real Gamma, const int XYZ )
{

   FluVar Slope_L, Slope_C, Slope_R, Slope_A, Slope_Limiter;
   real   Slope_LR, Sign;

// evaluate different slopes
#  define GetSlope( comp )                                                                               \
   {                                                                                                     \
      Slope_L.comp = C0.comp - L1.comp;                                                                  \
      Slope_R.comp = R1.comp - C0.comp;                                                                  \
      Slope_C.comp = (real)0.5*( Slope_L.comp + Slope_R.comp );                                          \
   } // GetSlope

#  define GetSlope_vL_GMinMod( comp )                                                                    \
   {                                                                                                     \
      Slope_LR = Slope_L.comp*Slope_R.comp;                                                              \
                                                                                                         \
      if ( Slope_LR > (real)0.0 )                                                                        \
         Slope_A.comp = Slope_LR / Slope_C.comp;                                                         \
      else                                                                                               \
         Slope_A.comp = (real)0.0;                                                                       \
   } // GetSlope_vL_GMinMod

#  define GetSlope_ExtPre( comp )                                                                        \
   {                                                                                                     \
      Slope_L2.comp = L1.comp - L2.comp;                                                                 \
      Slope_R2.comp = R2.comp - R1.comp;                                                                 \
   } // GetSlope_ExtPre

   GetSlope( Rho );
   GetSlope( Px  );
   GetSlope( Py  );
   GetSlope( Pz  );
   GetSlope( Egy );
#  if ( NCOMP_PASSIVE > 0 )
   for (int v=0; v<NCOMP_PASSIVE; v++)
   GetSlope( Passive[v] );
#  endif

   if ( LR_Limiter == VL_GMINMOD )
   {
      GetSlope_vL_GMinMod( Rho );
      GetSlope_vL_GMinMod( Px  );
      GetSlope_vL_GMinMod( Py  );
      GetSlope_vL_GMinMod( Pz  );
      GetSlope_vL_GMinMod( Egy );
#     if ( NCOMP_PASSIVE > 0 )
      for (int v=0; v<NCOMP_PASSIVE; v++)
      GetSlope_vL_GMinMod( Passive[v] );
#     endif
   }

// for the extrema-preserving limiter
#  if ( LR_SCHEME == PLM )
   FluVar Slope_L2, Slope_R2;
   real   D2_L, D2_R, D2_C, D2_Sign, D2_Limiter, Slope_Sign;

   if ( LR_Limiter == EXTPRE )
   {
      GetSlope_ExtPre( Rho );
      GetSlope_ExtPre( Px  );
      GetSlope_ExtPre( Py  );
      GetSlope_ExtPre( Pz  );
      GetSlope_ExtPre( Egy );
#     if ( NCOMP_PASSIVE > 0 )
      for (int v=0; v<NCOMP_PASSIVE; v++)
      GetSlope_ExtPre( Passive[v] );
#     endif
   }
#  endif // # if ( LR_SCHEME == PLM )

// primitive variables --> characteristic variables
#  ifdef CHAR_RECONSTRUCTION
   Slope_L = Pri2Char( Slope_L, Gamma, C0.Rho, C0.Egy, XYZ );
   Slope_R = Pri2Char( Slope_R, Gamma, C0.Rho, C0.Egy, XYZ );
   Slope_C = Pri2Char( Slope_C, Gamma, C0.Rho, C0.Egy, XYZ );

   if ( LR_Limiter == VL_GMINMOD )
      Slope_A = Pri2Char( Slope_A, Gamma, C0.Rho, C0.Egy, XYZ );

#  if ( LR_SCHEME == PLM )
   if ( LR_Limiter == EXTPRE )
   {
      Slope_L2 = Pri2Char( Slope_L2, Gamma, C0.Rho, C0.Egy, XYZ );
      Slope_R2 = Pri2Char( Slope_R2, Gamma, C0.Rho, C0.Egy, XYZ );
   }
#  endif
#  endif // #ifdef CHAR_RECONSTRUCTION


#  undef GetSlope
#  undef GetSlope_vL_GMinMod
#  undef GetSlope_ExtPre


// for the extrema-preserving limiter
#  if ( LR_SCHEME == PLM )

#     define EXTPRE_CONDITION( comp )                                                                    \
         ( LR_Limiter != EXTPRE || Slope_L2.comp*Slope_R2.comp > (real)0.0 )

#     define EXTPRE_LIMITER( comp )                                                                      \
      {                                                                                                  \
         if ( LR_Limiter == EXTPRE )                                                                     \
         {                                                                                               \
            D2_L = Slope_L .comp - Slope_L2.comp;                                                        \
            D2_R = Slope_R2.comp - Slope_R .comp;                                                        \
            D2_C = Slope_R .comp - Slope_L .comp;                                                        \
                                                                                                         \
            D2_Sign    = SIGN( D2_C );                                                                   \
            Slope_Sign = SIGN( Slope_C.comp );                                                           \
                                                                                                         \
            D2_Limiter = FMIN(  FABS(D2_C),                                                              \
                                FMIN( FMAX(D2_Sign*D2_L, (real)0.0), FMAX(D2_Sign*D2_R, (real)0.0) )  ); \
                                                                                                         \
            if ( D2_Sign*Slope_Sign < (real)0.0 )                                                        \
               Slope_Limiter.comp = FMIN( (real)1.5*EP_Coeff*D2_Limiter,                                 \
                                          MinMod_Coeff*FABS(Slope_L.comp) );                             \
            else                                                                                         \
               Slope_Limiter.comp = FMIN( (real)1.5*EP_Coeff*D2_Limiter,                                 \
                                          MinMod_Coeff*FABS(Slope_R.comp) );                             \
                                                                                                         \
            Slope_Limiter.comp = Slope_Sign * FMIN( FABS(Slope_C.comp), Slope_Limiter.comp );            \
         }                                                                                               \
      } // EXTPRE_LIMITER

#  else // ( LR_SCHEME == PPM )

#     define EXTPRE_CONDITION( comp ) true
#     define EXTPRE_LIMITER( comp )   /* nothing to do */

#  endif // #if ( LR_SCHEME == PLM ) ... else ...


// apply the slope limiter
#  define Limiter( comp )                                                                                \
   {                                                                                                     \
      Slope_Limiter.comp = (real)0.0;                                                                    \
      Slope_LR           = Slope_R.comp*Slope_L.comp;                                                    \
                                                                                                         \
      if (  Slope_LR > (real)0.0  &&  EXTPRE_CONDITION( comp )  )                                        \
      {                                                                                                  \
         switch ( LR_Limiter )                                                                           \
         {                                                                                               \
            case VANLEER :               /* van-Leer */                                                  \
               Slope_Limiter.comp = (real)2.0*Slope_LR/( Slope_L.comp + Slope_R.comp );                  \
               break;                                                                                    \
                                                                                                         \
            case GMINMOD : case EXTPRE : /* generalized MinMod & extrema-preserving */                   \
               Sign                = SIGN( Slope_C.comp );                                               \
                                                                                                         \
               Slope_L.comp       *= MinMod_Coeff;                                                       \
               Slope_R.comp       *= MinMod_Coeff;                                                       \
                                                                                                         \
               Slope_L.comp        = FABS( Slope_L.comp );                                               \
               Slope_C.comp        = FABS( Slope_C.comp );                                               \
               Slope_R.comp        = FABS( Slope_R.comp );                                               \
                                                                                                         \
               Slope_Limiter.comp  = FMIN( Slope_L.comp, Slope_R.comp );                                 \
               Slope_Limiter.comp  = FMIN( Slope_C.comp, Slope_Limiter.comp );                           \
               Slope_Limiter.comp *= Sign;                                                               \
               break;                                                                                    \
                                                                                                         \
            case ALBADA :                /* van-Albada */                                                \
               Slope_Limiter.comp = Slope_LR*( Slope_L.comp + Slope_R.comp ) /                           \
                                    ( Slope_L.comp*Slope_L.comp + Slope_R.comp*Slope_R.comp );           \
               break;                                                                                    \
                                                                                                         \
            case VL_GMINMOD :            /* van-Leer + generalized MinMod */                             \
               Sign                = SIGN( Slope_C.comp );                                               \
                                                                                                         \
               Slope_L.comp       *= MinMod_Coeff;                                                       \
               Slope_R.comp       *= MinMod_Coeff;                                                       \
                                                                                                         \
               Slope_L.comp        = FABS( Slope_L.comp );                                               \
               Slope_C.comp        = FABS( Slope_C.comp );                                               \
               Slope_R.comp        = FABS( Slope_R.comp );                                               \
               Slope_A.comp        = FABS( Slope_A.comp );                                               \
                                                                                                         \
               Slope_Limiter.comp  = FMIN( Slope_L.comp, Slope_R.comp );                                 \
               Slope_Limiter.comp  = FMIN( Slope_C.comp, Slope_Limiter.comp );                           \
               Slope_Limiter.comp  = FMIN( Slope_A.comp, Slope_Limiter.comp );                           \
               Slope_Limiter.comp *= Sign;                                                               \
               break;                                                                                    \
         }                                                                                               \
      } /* if (  Slope_LR > (real)0.0  &&  EXTPRE_CONDITION( comp )  ) */                                \
                                                                                                         \
      else                                                                                               \
         EXTPRE_LIMITER( comp );                                                                         \
   } // LimitSlope

   Limiter( Rho );
   Limiter( Px  );
   Limiter( Py  );
   Limiter( Pz  );
   Limiter( Egy );
#  if ( NCOMP_PASSIVE > 0 )
   for (int v=0; v<NCOMP_PASSIVE; v++)
   Limiter( Passive[v] );
#  endif

#  undef EXTPRE_CONDITION
#  undef EXTPRE_LIMITER
#  undef Limiter


// characteristic variables --> primitive variables
#  ifdef CHAR_RECONSTRUCTION
   Slope_Limiter = Char2Pri( Slope_Limiter, Gamma, C0.Rho, C0.Egy, XYZ );
#  endif


   return Slope_Limiter;

} // FUNCTION : LimitSlope



#if ( FLU_SCHEME == MHM )
//-------------------------------------------------------------------------------------------------------
// Function    :  HancockPredict
// Description :  Evolve the face-centered variables by half time-step by calculating the face-centered fluxes
//                (no Riemann solver is required)
//
// Note        :  1. Work for the MHM scheme
//                2. Do NOT require data in the neighboring cells
//
// Parameter   :  g_FC_Var_xL  : Face-centered conserved global variables on the -x surface
//                g_FC_Var_xR  : Face-centered conserved global variables on the +x surface
//                g_FC_Var_yL  : Face-centered conserved global variables on the -y surface
//                g_FC_Var_yR  : Face-centered conserved global variables on the +y surface
//                g_FC_Var_zL  : Face-centered conserved global variables on the -z surface
//                g_FC_Var_zR  : Face-centered conserved global variables on the +z surface
//                FC_Var_xL    : Face-centered conserved local  variables on the -x surface
//                FC_Var_xR    : Face-centered conserved local  variables on the +x surface
//                FC_Var_yL    : Face-centered conserved local  variables on the -y surface
//                FC_Var_yR    : Face-centered conserved local  variables on the +y surface
//                FC_Var_zL    : Face-centered conserved local  variables on the -z surface
//                FC_Var_zR    : Face-centered conserved local  variables on the +z surface
//                dt           : Time interval to advance solution
//                _dh          : 1 / grid size
//                Gamma        : Ratio of specific heats
//                ID_Out       : Array index to store the output data
//                P_In         : Input primitive variables (data before update)
//                               --> For checking negative density and pressure
//                MinDens/Pres : Minimum allowed density and pressure
//                NormPassive  : true --> convert passive scalars to mass fraction
//                NNorm        : Number of passive scalars for the option "NormPassive"
//                               --> Should be set to the global variable "PassiveNorm_NVar"
//                NormIdx      : Target variable indices for the option "NormPassive"
//                               --> Should be set to the global variable "PassiveNorm_VarIdx"
//-------------------------------------------------------------------------------------------------------
__device__ void HancockPredict( real g_FC_Var_xL[][NCOMP_TOTAL][ N_FC_VAR*N_FC_VAR*N_FC_VAR ],
                                real g_FC_Var_xR[][NCOMP_TOTAL][ N_FC_VAR*N_FC_VAR*N_FC_VAR ],
                                real g_FC_Var_yL[][NCOMP_TOTAL][ N_FC_VAR*N_FC_VAR*N_FC_VAR ],
                                real g_FC_Var_yR[][NCOMP_TOTAL][ N_FC_VAR*N_FC_VAR*N_FC_VAR ],
                                real g_FC_Var_zL[][NCOMP_TOTAL][ N_FC_VAR*N_FC_VAR*N_FC_VAR ],
                                real g_FC_Var_zR[][NCOMP_TOTAL][ N_FC_VAR*N_FC_VAR*N_FC_VAR ],
                                FluVar FC_Var_xL, FluVar FC_Var_xR,
                                FluVar FC_Var_yL, FluVar FC_Var_yR,
                                FluVar FC_Var_zL, FluVar FC_Var_zR,
                                const real dt, const real _dh, const real Gamma, const uint ID_Out,
                                const FluVar P_In, const real MinDens, const real MinPres,
                                const bool NormPassive, const int NNorm, const int NormIdx[] )
{

   const uint  bx       = blockIdx.x;
   const real  Gamma_m1 = Gamma - (real)1.0;
   const real _Gamma_m1 = (real)1.0 / Gamma_m1;
   const real  dt_dh2   = (real)0.5*dt*_dh;

   FluVar FC_Flux_xL, FC_Flux_xR, FC_Flux_yL, FC_Flux_yR, FC_Flux_zL, FC_Flux_zR;
   real   FluxDiff;


// 1. evaluate the fluxes
   FC_Flux_xL = CUFLU_Con2Flux( FC_Var_xL, Gamma_m1, 0, MinPres );
   FC_Flux_xR = CUFLU_Con2Flux( FC_Var_xR, Gamma_m1, 0, MinPres );
   FC_Flux_yL = CUFLU_Con2Flux( FC_Var_yL, Gamma_m1, 1, MinPres );
   FC_Flux_yR = CUFLU_Con2Flux( FC_Var_yR, Gamma_m1, 1, MinPres );
   FC_Flux_zL = CUFLU_Con2Flux( FC_Var_zL, Gamma_m1, 2, MinPres );
   FC_Flux_zR = CUFLU_Con2Flux( FC_Var_zR, Gamma_m1, 2, MinPres );


// 2. get the half-step solution
#  define HalfStepUpdate( comp )                                           \
   {                                                                       \
      FluxDiff = dt_dh2 * (  FC_Flux_xR.comp - FC_Flux_xL.comp +           \
                             FC_Flux_yR.comp - FC_Flux_yL.comp +           \
                             FC_Flux_zR.comp - FC_Flux_zL.comp  );         \
                                                                           \
      FC_Var_xL.comp -= FluxDiff;                                          \
      FC_Var_xR.comp -= FluxDiff;                                          \
      FC_Var_yL.comp -= FluxDiff;                                          \
      FC_Var_yR.comp -= FluxDiff;                                          \
      FC_Var_zL.comp -= FluxDiff;                                          \
      FC_Var_zR.comp -= FluxDiff;                                          \
   } // HalfStepUpdate

   HalfStepUpdate( Rho );
   HalfStepUpdate( Px  );
   HalfStepUpdate( Py  );
   HalfStepUpdate( Pz  );
   HalfStepUpdate( Egy );
#  if ( NCOMP_PASSIVE > 0 )
   for (int v=0; v<NCOMP_PASSIVE; v++)
   HalfStepUpdate( Passive[v] );
#  endif

#  undef HalfStepUpdate


// 3. enforce the positive density and pressure
// check negative density and energy
   if ( FC_Var_xL.Rho <= (real)0.0  ||  FC_Var_xR.Rho <= (real)0.0  ||  FC_Var_yL.Rho <= (real)0.0  ||
        FC_Var_yR.Rho <= (real)0.0  ||  FC_Var_zL.Rho <= (real)0.0  ||  FC_Var_zR.Rho <= (real)0.0  ||
        FC_Var_xL.Egy <= (real)0.0  ||  FC_Var_xR.Egy <= (real)0.0  ||  FC_Var_yL.Egy <= (real)0.0  ||
        FC_Var_yR.Egy <= (real)0.0  ||  FC_Var_zL.Egy <= (real)0.0  ||  FC_Var_zR.Egy <= (real)0.0   )
   {
      const FluVar C_In = CUFLU_Pri2Con( P_In, _Gamma_m1, NormPassive, NNorm, NormIdx );

      FC_Var_xL.Rho = FC_Var_xR.Rho = FC_Var_yL.Rho = FC_Var_yR.Rho = FC_Var_zL.Rho = FC_Var_zR.Rho = C_In.Rho;
      FC_Var_xL.Px  = FC_Var_xR.Px  = FC_Var_yL.Px  = FC_Var_yR.Px  = FC_Var_zL.Px  = FC_Var_zR.Px  = C_In.Px;
      FC_Var_xL.Py  = FC_Var_xR.Py  = FC_Var_yL.Py  = FC_Var_yR.Py  = FC_Var_zL.Py  = FC_Var_zR.Py  = C_In.Py;
      FC_Var_xL.Pz  = FC_Var_xR.Pz  = FC_Var_yL.Pz  = FC_Var_yR.Pz  = FC_Var_zL.Pz  = FC_Var_zR.Pz  = C_In.Pz;
      FC_Var_xL.Egy = FC_Var_xR.Egy = FC_Var_yL.Egy = FC_Var_yR.Egy = FC_Var_zL.Egy = FC_Var_zR.Egy = C_In.Egy;

#     if ( NCOMP_PASSIVE > 0 )
      for (int v=0; v<NCOMP_PASSIVE; v++)
      {
         FC_Var_xL.Passive[v] = FC_Var_xR.Passive[v] = FC_Var_yL.Passive[v] = FC_Var_yR.Passive[v]
                              = FC_Var_zL.Passive[v] = FC_Var_zR.Passive[v] = C_In.Passive[v];
      }
#     endif
   }

// ensure positive density and pressure
   FC_Var_xL.Rho = FMAX( FC_Var_xL.Rho, MinDens );
   FC_Var_xR.Rho = FMAX( FC_Var_xR.Rho, MinDens );
   FC_Var_yL.Rho = FMAX( FC_Var_yL.Rho, MinDens );
   FC_Var_yR.Rho = FMAX( FC_Var_yR.Rho, MinDens );
   FC_Var_zL.Rho = FMAX( FC_Var_zL.Rho, MinDens );
   FC_Var_zR.Rho = FMAX( FC_Var_zR.Rho, MinDens );

   FC_Var_xL.Egy = CUFLU_CheckMinPresInEngy( FC_Var_xL, Gamma_m1, _Gamma_m1, MinPres );
   FC_Var_xR.Egy = CUFLU_CheckMinPresInEngy( FC_Var_xR, Gamma_m1, _Gamma_m1, MinPres );
   FC_Var_yL.Egy = CUFLU_CheckMinPresInEngy( FC_Var_yL, Gamma_m1, _Gamma_m1, MinPres );
   FC_Var_yR.Egy = CUFLU_CheckMinPresInEngy( FC_Var_yR, Gamma_m1, _Gamma_m1, MinPres );
   FC_Var_zL.Egy = CUFLU_CheckMinPresInEngy( FC_Var_zL, Gamma_m1, _Gamma_m1, MinPres );
   FC_Var_zR.Egy = CUFLU_CheckMinPresInEngy( FC_Var_zR, Gamma_m1, _Gamma_m1, MinPres );

#  if ( NCOMP_PASSIVE > 0 )
   for (int v=0; v<NCOMP_PASSIVE; v++)
   {
      FC_Var_xL.Passive[v] = FMAX( FC_Var_xL.Passive[v], TINY_NUMBER );
      FC_Var_xR.Passive[v] = FMAX( FC_Var_xR.Passive[v], TINY_NUMBER );
      FC_Var_yL.Passive[v] = FMAX( FC_Var_yL.Passive[v], TINY_NUMBER );
      FC_Var_yR.Passive[v] = FMAX( FC_Var_yR.Passive[v], TINY_NUMBER );
      FC_Var_zL.Passive[v] = FMAX( FC_Var_zL.Passive[v], TINY_NUMBER );
      FC_Var_zR.Passive[v] = FMAX( FC_Var_zR.Passive[v], TINY_NUMBER );
   }
#  endif


// 4. store the half-step face-centered variables to the global memory
#  define Store( g_array, FC_Var )                                \
   {                                                              \
      g_array[bx][ 0][ID_Out] = FC_Var.Rho;                       \
      g_array[bx][ 1][ID_Out] = FC_Var.Px;                        \
      g_array[bx][ 2][ID_Out] = FC_Var.Py;                        \
      g_array[bx][ 3][ID_Out] = FC_Var.Pz;                        \
      g_array[bx][ 4][ID_Out] = FC_Var.Egy;                       \
                                                                  \
      for (int v=0, vv=NCOMP_FLUID; v<NCOMP_PASSIVE; v++, vv++)   \
      g_array[bx][vv][ID_Out] = FC_Var.Passive[v];                \
   } // Store

   Store( g_FC_Var_xL, FC_Var_xL );
   Store( g_FC_Var_xR, FC_Var_xR );
   Store( g_FC_Var_yL, FC_Var_yL );
   Store( g_FC_Var_yR, FC_Var_yR );
   Store( g_FC_Var_zL, FC_Var_zL );
   Store( g_FC_Var_zR, FC_Var_zR );

#  undef Store

} // FUNCTION : HancockPredict
#endif // #if ( FLU_SCHEME == MHM )



#if ( FLU_SCHEME == CTU )
//-------------------------------------------------------------------------------------------------------
// Function    :  Get_EigenSystem
// Description :  Evaluate the eigenvalues and left/right eigenvectors
//
// Note        :  1. The input data should be primitive variables
//                2. The constant components of eigenvectors should be initialized in advance
//                3. Work for the CTU scheme
//
// Parameter   :  CC_Var   : Array storing the input cell-centered primitive variables
//                EVal     : Array to store the output eigenvalues (in three spatial directions)
//                LEVec    : Array to store the output left eigenvectors
//                REVec    : Array to store the output right eigenvectors
//                Gamma    : Ratio of specific heats
//-------------------------------------------------------------------------------------------------------
__device__ void Get_EigenSystem( const FluVar CC_Var, FluVar5 &EVal_x, FluVar5 &EVal_y, FluVar5 &EVal_z,
                                 FluVar5 &LEVec1, FluVar5 &LEVec2, FluVar5 &LEVec3, FluVar5 &LEVec4,
                                 FluVar5 &LEVec5, FluVar5 &REVec1, FluVar5 &REVec2, FluVar5 &REVec3,
                                 FluVar5 &REVec4, FluVar5 &REVec5, const real Gamma )
{

#  ifdef CHECK_NEGATIVE_IN_FLUID
   if ( CUFLU_CheckNegative(CC_Var.Egy) )
      printf( "ERROR : negative pressure (%14.7e) at file <%s>, line <%d>, function <%s>\n",
              CC_Var.Egy, __FILE__, __LINE__, __FUNCTION__ );

   if ( CUFLU_CheckNegative(CC_Var.Rho) )
      printf( "ERROR : negative density (%14.7e) at file <%s>, line <%d>, function <%s>\n",
              CC_Var.Rho, __FILE__, __LINE__, __FUNCTION__ );
#  endif

   const real  Rho = CC_Var.Rho;
   const real _Rho = (real)1.0/Rho;
   const real  Vx  = CC_Var.Px;
   const real  Vy  = CC_Var.Py;
   const real  Vz  = CC_Var.Pz;
   const real  Cs2 = Gamma*CC_Var.Egy*_Rho;
   const real  Cs  = SQRT( Cs2 );
   const real _Cs  = (real)1.0/Cs;
   const real _Cs2 = _Cs*_Cs;


// a. eigenvalues in three spatial directions
   EVal_x.Rho = Vx - Cs;
   EVal_x.Px  = Vx;
   EVal_x.Py  = Vx;
   EVal_x.Pz  = Vx;
   EVal_x.Egy = Vx + Cs;

   EVal_y.Rho = Vy - Cs;
   EVal_y.Px  = Vy;
   EVal_y.Py  = Vy;
   EVal_y.Pz  = Vy;
   EVal_y.Egy = Vy + Cs;

   EVal_z.Rho = Vz - Cs;
   EVal_z.Px  = Vz;
   EVal_z.Py  = Vz;
   EVal_z.Pz  = Vz;
   EVal_z.Egy = Vz + Cs;


// NOTE : the left and right eigenvectors have the same form in different directions
// b. left eigenvectors (rows of the matrix LEVec)
   LEVec1.Px  = -(real)0.5*Rho*_Cs;
   LEVec1.Egy = +(real)0.5*_Cs2;
   LEVec2.Egy = -_Cs2;
   LEVec5.Px  = -LEVec1.Px;
   LEVec5.Egy = +LEVec1.Egy;


// c. right eigenvectors (rows of the matrix REVec)
   REVec1.Px  = -Cs*_Rho;
   REVec1.Egy = +Cs2;
   REVec5.Px  = -REVec1.Px;
   REVec5.Egy = +Cs2;

} // FUNCTION : Get_EigenSystem
#endif // #if ( FLU_SCHEME == CTU )



#endif // #ifndef __CUFLU_DATARECONSTRUCTION_CU__
