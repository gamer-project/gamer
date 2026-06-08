#include "CUFLU.h"

#if ( MODEL == HYDRO  &&  FLU_SCHEME == CTU  &&  !defined SRHD )



// external functions
#ifdef __CUDACC__

#include "CUFLU_Shared_FluUtility.cu"
#include "CUFLU_Shared_DataReconstruction.cu"
#include "CUFLU_Shared_ComputeFlux.cu"
#include "CUFLU_Shared_FullStepUpdate.cu"
#ifdef MHD
#include "CUFLU_Shared_ConstrainedTransport.cu"
#endif

#include "CUDA_ConstMemory.h"

#else // #ifdef __CUDACC__

void Hydro_Rotate3D( real InOut[], const int XYZ, const bool Forward, const int Mag_Offset );
void Hydro_DataReconstruction( const real g_ConVar   [][ CUBE(FLU_NXT) ],
                               const real g_FC_B     [][ SQR(FLU_NXT)*FLU_NXT_P1 ],
                                     real g_PriVar   [][ CUBE(FLU_NXT) ],
                                     real g_FC_Var   [][NCOMP_TOTAL_PLUS_MAG][ CUBE(N_FC_VAR) ],
                                     real g_Slope_PPM[][NCOMP_LR            ][ CUBE(N_SLOPE_PPM) ],
                                     real g_EC_Ele   [][ CUBE(N_EC_ELE) ],
                               const bool Con2Pri, const LR_Limiter_t LR_Limiter, const real MinMod_Coeff,
                               const real dt, const real dh,
                               const real MinDens, const real MinPres, const real MinEint,
                               const long PassiveFloor, const bool FracPassive, const int NFrac, const int FracIdx[],
                               const bool JeansMinPres, const real JeansMinPres_Coeff,
                               const EoS_t *EoS );
void Hydro_ComputeFlux( const real g_FC_Var [][NCOMP_TOTAL_PLUS_MAG][ CUBE(N_FC_VAR) ],
                              real g_FC_Flux[][NCOMP_TOTAL_PLUS_MAG][ CUBE(N_FC_FLUX) ],
                        const int NFlux, const int NSkip_N, const int NSkip_T,
                        const bool CorrHalfVel, const real g_Pot_USG[], const double g_Corner[],
                        const real dt, const real dh, const double Time, const bool UsePot,
                        const OptExtAcc_t ExtAcc, const ExtAcc_t ExtAcc_Func, const double ExtAcc_AuxArray[],
                        const real MinDens, const real MinPres, const long PassiveFloor, const EoS_t *EoS );
void Hydro_StoreIntFlux( const real g_FC_Flux[][NCOMP_TOTAL_PLUS_MAG][ CUBE(N_FC_FLUX) ],
                               real g_IntFlux[][NCOMP_TOTAL][ SQR(PS2) ],
                         const int NFlux );
void Hydro_FullStepUpdate( const real g_Input[][ CUBE(FLU_NXT) ], real g_Output[][ CUBE(PS2) ], char g_DE_Status[],
                           const real g_FC_B[][ PS2P1*SQR(PS2) ], const real g_Flux[][NCOMP_TOTAL_PLUS_MAG][ CUBE(N_FC_FLUX) ],
                           const real dt, const real dh, const real MinDens, const real MinEint, const real DualEnergySwitch,
                           const long PassiveFloor, const bool NormPassive, const int NNorm, const int NormIdx[],
                           const EoS_t *EoS, int *s_FullStepFailure, const int Iteration, const int MinMod_MaxIter );
#ifdef MHD
void MHD_ComputeElectric(       real g_EC_Ele[][ CUBE(N_EC_ELE) ],
                          const real g_FC_Flux[][NCOMP_TOTAL_PLUS_MAG][ CUBE(N_FC_FLUX) ],
                          const real g_PriVar[][ CUBE(FLU_NXT) ],
                          const int NEle, const int NFlux, const int NPri, const int OffsetPri,
                          const real dt, const real dh,
                          const bool DumpIntEle, real g_IntEle[][NCOMP_ELE][ PS2P1*PS2 ],
                          const bool CorrHalfVel, const real g_Pot_USG[], const double g_Corner[], const double Time,
                          const bool UsePot, const OptExtAcc_t ExtAcc, const ExtAcc_t ExtAcc_Func,
                          const double ExtAcc_AuxArray[] );
void MHD_UpdateMagnetic( real *g_FC_Bx_Out, real *g_FC_By_Out, real *g_FC_Bz_Out,
                         const real g_FC_B_In[][ FLU_NXT_P1*SQR(FLU_NXT) ],
                         const real g_EC_Ele[][ CUBE(N_EC_ELE) ],
                         const real dt, const real dh, const int NOut, const int NEle, const int Offset_B_In );
void MHD_HalfStepPrimitive( const real g_Flu_In[][ CUBE(FLU_NXT) ],
                            const real g_FC_B_Half[][ FLU_NXT_P1*SQR(FLU_NXT) ],
                                  real g_PriVar_Out[][ CUBE(FLU_NXT) ],
                            const real g_Flux[][NCOMP_TOTAL_PLUS_MAG][ CUBE(N_FC_FLUX) ],
                            const real dt, const real dh, const real MinDens );
#endif // #ifdef MHD

#endif // #ifdef __CUDACC__ ... else ...


// internal functions
GPU_DEVICE
void Hydro_TGradientCorrection(       real g_FC_Var   [][NCOMP_TOTAL_PLUS_MAG][ CUBE(N_FC_VAR)  ],
                                const real g_FC_Flux  [][NCOMP_TOTAL_PLUS_MAG][ CUBE(N_FC_FLUX) ],
                                const real g_FC_B_In  [][ FLU_NXT_P1*SQR(FLU_NXT) ],
                                const real g_FC_B_Half[][ FLU_NXT_P1*SQR(FLU_NXT) ],
                                const real g_EC_Ele   [][ CUBE(N_EC_ELE) ],
                                const real g_PriVar   [][ CUBE(FLU_NXT) ],
                                const real dt, const real dh, const real MinDens, const real MinEint,
                                const long PassiveFloor );




//-------------------------------------------------------------------------------------------------------
// Function    :  CPU/CUFLU_FluidSolver_CTU
// Description :  CPU/GPU fluid solver based on the Corner-Transport-Upwind (CTU) scheme
//
// Note        :  1. Ref: (a) Stone et al., ApJS, 178, 137 (2008)
//                        (b) Gardiner & Stone, J. Comput. Phys., 227, 4123 (2008)
//                2. See include/CUFLU.h for the values and description of different symbolic constants
//                   such as N_FC_VAR, N_FC_FLUX, N_SLOPE_PPM, N_FL_FLUX, N_HF_VAR
//                3. Arrays with a prefix "g_" are stored in the global memory of GPU
//
// Parameter   :  g_Flu_Array_In     : Array storing the input fluid variables
//                g_Flu_Array_Out    : Array to store the output fluid variables
//                g_Mag_Array_In     : Array storing the input B field (for MHD only)
//                g_Mag_Array_Out    : Array to store the output B field (for MHD only)
//                g_DE_Array_Out     : Array to store the dual-energy status
//                g_Flux_Array       : Array to store the output fluxes
//                g_Ele_Array        : Array to store the output electric field (for MHD only)
//                g_Corner_Array     : Array storing the physical corner coordinates of each patch group (for UNSPLIT_GRAVITY)
//                g_Pot_Array_USG    : Array storing the input potential for UNSPLIT_GRAVITY
//                g_PriVar           : Array to store the primitive variables
//                g_Slope_PPM        : Array to store the slope for the PPM reconstruction
//                g_FC_Var           : Array to store the half-step variables
//                g_FC_Flux          : Array to store the face-centered fluxes
//                g_FC_Mag_Half      : Array to store the half-step B field (for MHD only)
//                g_EC_Ele           : Array to store the edge-centered electric field (for MHD only)
//                NPatchGroup        : Number of patch groups to be evaluated
//                dt                 : Time interval to advance solution
//                dh                 : Cell size
//                StoreFlux          : true --> store the coarse-fine fluxes
//                StoreElectric      : true --> store the coarse-fine electric field
//                LR_Limiter         : Slope limiter for the data reconstruction in the MHM/MHM_RP/CTU schemes
//                                     (0/1/2/3/4) = (vanLeer/generalized MinMod/vanAlbada/
//                                                    vanLeer + generalized MinMod/extrema-preserving) limiter
//                MinMod_Coeff       : Coefficient of the generalized MinMod limiter
//                Time               : Current physical time                                 (for UNSPLIT_GRAVITY only)
//                UsePot             : Add self-gravity and/or external potential            (for UNSPLIT_GRAVITY only)
//                ExtAcc             : Add external acceleration                             (for UNSPLIT_GRAVITY only)
//                ExtAcc_Func        : Function pointer to the external acceleration routine (for UNSPLIT_GRAVITY only)
//                c_ExtAcc_AuxArray  : Auxiliary array for adding external acceleration      (for UNSPLIT_GRAVITY and CPU only)
//                                     --> When using GPU, this array is stored in the constant memory header
//                                         CUDA_ConstMemory.h and does not need to be passed as a function argument
//                MinDens/Pres/Eint  : Density, pressure, and internal energy floors
//                DualEnergySwitch   : Use the dual-energy formalism if E_int/E_kin < DualEnergySwitch
//                PassiveFloor       : Bitwise flag to specify the passive scalars to be floored
//                NormPassive        : true --> normalize passive scalars so that the sum of their mass density
//                                              is equal to the gas mass density
//                NNorm              : Number of passive scalars to be normalized
//                                     --> Should be set to the global variable "PassiveNorm_NVar"
//                c_NormIdx          : Target variable indices to be normalized
//                                     --> Should be set to the global variable "PassiveNorm_VarIdx"
//                                     --> When using GPU, this array is stored in the constant memory and does
//                                         not need to be passed as a function argument
//                                         --> Declared in CUDA_ConstMemory.h with the prefix "c_" to
//                                             highlight that this is a constant variable on GPU
//                FracPassive        : true --> convert passive scalars to mass fraction during data reconstruction
//                NFrac              : Number of passive scalars for the option "FracPassive"
//                                     --> Should be set to the global variable "PassiveIntFrac_NVar"
//                c_FracIdx          : Target variable indices for the option "FracPassive"
//                                     --> Should be set to the global variable "PassiveIntFrac_VarIdx"
//                                     --> When using GPU, this array is stored in the constant memory and does
//                                         not need to be passed as a function argument
//                                         --> Declared in CUDA_ConstMemory.h with the prefix "c_" to
//                                             highlight that this is a constant variable on GPU
//                JeansMinPres       : Apply minimum pressure estimated from the Jeans length
//                JeansMinPres_Coeff : Coefficient used by JeansMinPres = G*(Jeans_NCell*Jeans_dh)^2/(Gamma*pi);
//                EoS                : EoS object
//-------------------------------------------------------------------------------------------------------
#ifdef __CUDACC__
#ifndef MHD
// scalar slope limiter (component-wise; no CHAR_RECONSTRUCTION needed)
GPU_DEVICE
real Hydro_LimitSlope_Scalar( const real vL, const real vC, const real vR,
                               const LR_Limiter_t LR_Limiter, const real MinMod_Coeff )
{
   const real SL = vC - vL;
   const real SR = vR - vC;
   if ( SL * SR <= (real)0.0 ) return (real)0.0;
   const real SC = (real)0.5 * ( SL + SR );
   switch ( LR_Limiter )
   {
      case LR_LIMITER_CENTRAL:
         return SC;
      case LR_LIMITER_VANLEER:
         return (real)2.0 * SL * SR / ( SL + SR );
      case LR_LIMITER_GMINMOD:
         return SIGN(SC) * FMIN( MinMod_Coeff*FABS(SL), FMIN( MinMod_Coeff*FABS(SR), FABS(SC) ) );
      case LR_LIMITER_ALBADA:
         return SL * SR * ( SL + SR ) / ( SL*SL + SR*SR );
      case LR_LIMITER_VL_GMINMOD: {
         const real SA = (real)2.0 * SL * SR / ( SL + SR );
         return SIGN(SC) * FMIN( MinMod_Coeff*FABS(SL),
                           FMIN( MinMod_Coeff*FABS(SR),
                           FMIN( FABS(SC), FABS(SA) ) ) );
      }
      default: return (real)0.0;
   }
}

// fused kernel: Con2Pri + PPM slope + face-centered data reconstruction (replaces CUFLU_PPM_Slope + CUFLU_DR_FCVar)
// eliminates g_PriVar and g_Slope_PPM global roundtrip by computing Con2Pri inline during slab loading
// and computing PPM limited slopes on-the-fly from the 5-slab shared-memory window
__global__
__launch_bounds__( FLU_BLOCK_SIZE_X )
void CUFLU_DR_FCVar(
   const real   g_Flu_Array_In [][NCOMP_TOTAL][ CUBE(FLU_NXT) ],
   const real   g_Mag_Array_In [][NCOMP_MAG][ FLU_NXT_P1*SQR(FLU_NXT) ],
         real   g_FC_Var       [][6][NCOMP_TOTAL_PLUS_MAG][ CUBE(N_FC_VAR)    ],
   const LR_Limiter_t LR_Limiter, const real MinMod_Coeff,
   const real dt, const real dh,
   const real MinDens, const real MinPres, const real MinEint,
   const long PassiveFloor, const bool FracPassive, const int NFrac,
   const bool JeansMinPres, const real JeansMinPres_Coeff,
   const EoS_t EoS )
{
   const int P = blockIdx.x;
   const real (*const g_Flu_Array_In_1PG)[CUBE(FLU_NXT)] = g_Flu_Array_In[P];
         real (*const g_FC_Var_1PG   )[NCOMP_TOTAL_PLUS_MAG][CUBE(N_FC_VAR)] = g_FC_Var[P];

// 5 z-slabs of PriVar: covers (k_cc-2..k_cc+2) for each k_fc plane
// Con2Pri computed inline from g_Flu_Array_In; eliminates g_PriVar + g_Slope_PPM global roundtrip
   __shared__ real s_PriVar[5][NCOMP_LR][ SQR(FLU_NXT) ];

   const int NIn           = FLU_NXT;
   const int NGhost        = LR_GHOST_SIZE;
   const int N_FC_VAR2     = SQR( N_FC_VAR );
   const int NIn2          = SQR( NIn );
   const real dt_dh2       = (real)0.5*dt/dh;
   const int idx_wave[NWAVE] = { 0, 1, 2, 3, 4 };

#  if (  ( RSOLVER == HLLE || RSOLVER == HLLC || RSOLVER == HLLD )  &&  defined HLL_NO_REF_STATE  )
#  ifdef HLL_INCLUDE_ALL_WAVES
   const bool HLL_Include_All_Waves = true;
#  else
   const bool HLL_Include_All_Waves = false;
#  endif
#  endif

   real EigenVal[3][NWAVE];
   real LEigenVec[NWAVE][NWAVE] = { { 0.0, NULL_REAL, 0.0, 0.0, NULL_REAL },
                                    { 1.0,       0.0, 0.0, 0.0, NULL_REAL },
                                    { 0.0,       0.0, 1.0, 0.0,       0.0 },
                                    { 0.0,       0.0, 0.0, 1.0,       0.0 },
                                    { 0.0, NULL_REAL, 0.0, 0.0, NULL_REAL } };
   real REigenVec[NWAVE][NWAVE] = { { 1.0, NULL_REAL, 0.0, 0.0, NULL_REAL },
                                    { 1.0,       0.0, 0.0, 0.0,       0.0 },
                                    { 0.0,       0.0, 1.0, 0.0,       0.0 },
                                    { 0.0,       0.0, 0.0, 1.0,       0.0 },
                                    { 1.0, NULL_REAL, 0.0, 0.0, NULL_REAL } };

   for (int k_fc = 0; k_fc < N_FC_VAR; k_fc++)
   {
//    cooperatively load 5 z-slabs into shared memory, computing Con2Pri inline from g_Flu_Array_In
      const int k_cc = k_fc + NGhost;
      for (int slab = 0; slab < 5; slab++)
      {
         const int k_slab = k_cc - 2 + slab;
         CGPU_LOOP( idx_xy, NIn2 )
         {
            real ConVar[NCOMP_TOTAL], PriVar_tmp[NCOMP_LR];
            for (int v = 0; v < NCOMP_TOTAL; v++)
               ConVar[v] = g_Flu_Array_In_1PG[v][ k_slab*NIn2 + idx_xy ];
            Hydro_Con2Pri( ConVar, PriVar_tmp, MinPres, PassiveFloor, FracPassive, NFrac, c_FracIdx,
                           JeansMinPres, JeansMinPres_Coeff,
                           EoS.DensEint2Pres_FuncPtr, EoS.DensPres2Eint_FuncPtr,
                           EoS.GuessHTilde_FuncPtr, EoS.HTilde2Temp_FuncPtr,
                           EoS.AuxArrayDevPtr_Flt, EoS.AuxArrayDevPtr_Int, EoS.Table, NULL, NULL );
            for (int v = 0; v < NCOMP_LR; v++)
               s_PriVar[slab][v][idx_xy] = PriVar_tmp[v];
         }
      }
      __syncthreads();

//    process all (i_fc, j_fc) in this k-plane
      CGPU_LOOP( idx_ij, N_FC_VAR2 )
      {
         const int i_fc      = idx_ij % N_FC_VAR;
         const int j_fc      = idx_ij / N_FC_VAR;
         const int i_cc      = i_fc + NGhost;
         const int j_cc      = j_fc + NGhost;
         const int idx_fc_out = k_fc*N_FC_VAR2 + idx_ij;
         const int xy_C      = j_cc*NIn + i_cc;

         real cc_C_ncomp[NCOMP_LR], fcCon[6][NCOMP_LR], fcPri[6][NCOMP_LR], dfc[NCOMP_LR], dfc6[NCOMP_LR];

         for (int v = 0; v < NCOMP_LR; v++)   cc_C_ncomp[v] = s_PriVar[2][v][xy_C];

         Hydro_GetEigenSystem( cc_C_ncomp, EigenVal, LEigenVec, REigenVec, &EoS );

         for (int d = 0; d < 3; d++)
         {
            const int faceL      = 2*d;
            const int faceR      = faceL + 1;

//          slab/xy indices for +-1 neighbors; z-dir uses slabs 1/3, x/y-dir uses center slab 2
            const int slab_L = (d == 2) ? 1 : 2;
            const int slab_R = (d == 2) ? 3 : 2;
            const int  xy_L  = (d == 0) ? xy_C - 1     : (d == 1) ? xy_C - NIn   : xy_C;
            const int  xy_R  = (d == 0) ? xy_C + 1     : (d == 1) ? xy_C + NIn   : xy_C;
            const int  xy_LL = (d == 0) ? xy_C - 2     : (d == 1) ? xy_C - 2*NIn : -1;
            const int  xy_RR = (d == 0) ? xy_C + 2     : (d == 1) ? xy_C + 2*NIn : -1;

#           ifdef FLOAT8
            const real round_err = 1.e-12;
#           else
            const real round_err = 1.e-6;
#           endif
            const real C_factor  = 1.25;

            for (int v = 0; v < NCOMP_LR; v++)
            {
               real fc_L, fc_R;

               if ( LR_Limiter == LR_LIMITER_ATHENA )
               {
//                all +-1 and +-2 neighbors from shmem (5-slab window covers k_cc-2..k_cc+2)
                  const real cc_LL = (d < 2) ? s_PriVar[2    ][v][xy_LL]
                                             : s_PriVar[0    ][v][xy_C  ];
                  const real cc_L  = s_PriVar[slab_L][v][xy_L];
                  const real cc_C  = cc_C_ncomp[v];
                  const real cc_R  = s_PriVar[slab_R][v][xy_R];
                  const real cc_RR = (d < 2) ? s_PriVar[2    ][v][xy_RR]
                                             : s_PriVar[4    ][v][xy_C  ];

                  real tmp, rho, cc_abs_max;
                  real d_L, d_R, dd_L, dd_C, dd_R;
                  real dh_LL, dh_L, dh_R, dh_RR, ddh_L, ddh_C, ddh_R;

                  d_L  = cc_C  - cc_L;
                  d_R  = cc_R  - cc_C;
                  dd_L = cc_LL - (real)2.*cc_L + cc_C;
                  dd_C = cc_L  - (real)2.*cc_C + cc_R;
                  dd_R = cc_C  - (real)2.*cc_R + cc_RR;

                  cc_abs_max = FMAX( FABS(cc_LL), FMAX( FABS(cc_L), FMAX( FABS(cc_C), FMAX( FABS(cc_R), FABS(cc_RR) ) ) ) );

                  fc_L = ( -cc_LL + (real)7.*cc_L + (real)7.*cc_C - cc_R  ) / (real)12.0;
                  fc_R = ( -cc_L  + (real)7.*cc_C + (real)7.*cc_R - cc_RR ) / (real)12.0;

                  dh_LL = fc_L - cc_L;
                  dh_L  = cc_C - fc_L;
                  dh_R  = fc_R - cc_C;
                  dh_RR = cc_R - fc_R;
                  ddh_L = cc_L - (real)2.*fc_L + cc_C;
                  ddh_C = fc_L - (real)2.*cc_C + fc_R;
                  ddh_R = cc_C - (real)2.*fc_R + cc_R;

                  if ( dh_LL*dh_L < (real)0.0 ) {
                     if ( SIGN(dd_L) == SIGN(ddh_L)  &&  SIGN(ddh_L) == SIGN(dd_C) )
                        tmp = SIGN(dd_C)*FMIN( C_factor*FABS(dd_L), FMIN( (real)3.*FABS(ddh_L), C_factor*FABS(dd_C) ) );
                     else
                        tmp = (real)0.0;
                     fc_L  = (real)0.5*(cc_L+cc_C) - tmp/(real)6.0;
                     dh_L  = cc_C - fc_L;
                     ddh_C = fc_L - (real)2.*cc_C + fc_R;
                  }

                  if ( dh_R*dh_RR < (real)0.0 ) {
                     if ( SIGN(dd_C) == SIGN(ddh_R)  &&  SIGN(ddh_R) == SIGN(dd_R) )
                        tmp = SIGN(dd_C)*FMIN( C_factor*FABS(dd_C), FMIN( (real)3.*FABS(ddh_R), C_factor*FABS(dd_R) ) );
                     else
                        tmp = (real)0.0;
                     fc_R  = (real)0.5*(cc_C+cc_R) - tmp/(real)6.0;
                     dh_R  = fc_R - cc_C;
                     ddh_C = fc_L - (real)2.*cc_C + fc_R;
                  }

                  if ( SIGN(dd_L) == SIGN(dd_C)  &&  SIGN(dd_C) == SIGN(dd_R)  &&  SIGN(dd_R) == SIGN(ddh_C) )
                     tmp = SIGN(dd_C)*FMIN( C_factor*FABS(dd_L), FMIN( C_factor*FABS(dd_C), FMIN( C_factor*FABS(dd_R), (real)6.0*FABS(ddh_C) ) ) );
                  else
                     tmp = (real)0.0;

                  rho = ( (real)6.*FABS(ddh_C) > round_err*cc_abs_max ) ? tmp/ddh_C/(real)6. : (real)0.0;

                  if ( dh_L*dh_R < (real)0.0  ||  d_L*d_R < (real)0.0 ) {
                     if ( rho < (real)1.-round_err ) { fc_L = cc_C - rho*dh_L;  fc_R = cc_C + rho*dh_R; }
                  } else {
                     if ( FABS(dh_L) >= (real)2.*FABS(dh_R) )   fc_L = cc_C - (real)2.*dh_R;
                     if ( FABS(dh_R) >= (real)2.*FABS(dh_L) )   fc_R = cc_C + (real)2.*dh_L;
                  }

               } else { // non-ATHENA: all reads from shmem; slopes computed inline
                  const real cc_L_val = s_PriVar[slab_L][v][xy_L];
                  const real cc_R_val = s_PriVar[slab_R][v][xy_R];
                  const real cc_C_val = cc_C_ncomp[v];
//                +-2 cell values for slope stencil at C-1 and C+1
                  const real cc_LL_val = (d == 0) ? s_PriVar[2][v][xy_C - 2    ]
                                       : (d == 1) ? s_PriVar[2][v][xy_C - 2*NIn]
                                       :            s_PriVar[0][v][xy_C        ];
                  const real cc_RR_val = (d == 0) ? s_PriVar[2][v][xy_C + 2    ]
                                       : (d == 1) ? s_PriVar[2][v][xy_C + 2*NIn]
                                       :            s_PriVar[4][v][xy_C        ];
                  const real dcc_L = Hydro_LimitSlope_Scalar( cc_LL_val, cc_L_val, cc_C_val, LR_Limiter, MinMod_Coeff );
                  const real dcc_C = Hydro_LimitSlope_Scalar( cc_L_val,  cc_C_val, cc_R_val, LR_Limiter, MinMod_Coeff );
                  const real dcc_R = Hydro_LimitSlope_Scalar( cc_C_val,  cc_R_val, cc_RR_val, LR_Limiter, MinMod_Coeff );

                  fc_L = (real)0.5*( cc_C_val + cc_L_val ) - (real)1.0/(real)6.0*( dcc_C - dcc_L );
                  fc_R = (real)0.5*( cc_C_val + cc_R_val ) - (real)1.0/(real)6.0*( dcc_R - dcc_C );

                  if ( LR_Limiter == LR_LIMITER_CENTRAL )
                  {
                     if ( (cc_C_val-fc_L)*(fc_L-cc_L_val) < (real)0.0 )   fc_L = (real)0.5*(cc_C_val+cc_L_val);
                     if ( (cc_R_val-fc_R)*(fc_R-cc_C_val) < (real)0.0 )   fc_R = (real)0.5*(cc_C_val+cc_R_val);
                  }

                  dfc [v] = fc_R - fc_L;
                  dfc6[v] = (real)6.0*(  cc_C_val - (real)0.5*( fc_L + fc_R )  );

                  if (  ( fc_R - cc_C_val )*( cc_C_val - fc_L ) <= (real)0.0  )
                  {
                     fc_L = cc_C_val;
                     fc_R = cc_C_val;
                  }
                  else if ( dfc[v]*dfc6[v] > +dfc[v]*dfc[v] )
                     fc_L = (real)3.0*cc_C_val - (real)2.0*fc_R;
                  else if ( dfc[v]*dfc6[v] < -dfc[v]*dfc[v] )
                     fc_R = (real)3.0*cc_C_val - (real)2.0*fc_L;

                  real Min, Max;
                  Min  = ( cc_C_val < cc_L_val ) ? cc_C_val : cc_L_val;
                  Max  = ( cc_C_val > cc_L_val ) ? cc_C_val : cc_L_val;
                  fc_L = ( fc_L > Min ) ? fc_L : Min;
                  fc_L = ( fc_L < Max ) ? fc_L : Max;

                  Min  = ( cc_C_val < cc_R_val ) ? cc_C_val : cc_R_val;
                  Max  = ( cc_C_val > cc_R_val ) ? cc_C_val : cc_R_val;
                  fc_R = ( fc_R > Min ) ? fc_R : Min;
                  fc_R = ( fc_R < Max ) ? fc_R : Max;
               } // if ( LR_Limiter == LR_LIMITER_ATHENA ) ... else ...

               fcPri[faceL][v] = fc_L;
               fcPri[faceR][v] = fc_R;

            } // for (int v=0; v<NCOMP_LR; v++)


//          4. advance by half time-step (CTU)
            real Coeff_L, Coeff_R;
            real Correct_L[NCOMP_LR], Correct_R[NCOMP_LR];

            for (int v = 0; v < NCOMP_LR; v++)
            {
               dfc [v] = fcPri[faceR][v] - fcPri[faceL][v];
               dfc6[v] = (real)6.0*(  cc_C_ncomp[v] - (real)0.5*( fcPri[faceL][v] + fcPri[faceR][v] )  );
            }

            Hydro_Rotate3D( dfc,  d, true, MAG_OFFSET );
            Hydro_Rotate3D( dfc6, d, true, MAG_OFFSET );

#           if (  ( RSOLVER == HLLE || RSOLVER == HLLC || RSOLVER == HLLD )  &&  defined HLL_NO_REF_STATE  )
            for (int v = 0; v < NWAVE; v++)
            {
               Correct_L[ idx_wave[v] ] = (real)0.0;
               Correct_R[ idx_wave[v] ] = (real)0.0;
            }

            for (int Mode = 0; Mode < NWAVE; Mode++)
            {
               Coeff_L = (real)0.0;
               Coeff_R = (real)0.0;

               if ( HLL_Include_All_Waves  ||  EigenVal[d][Mode] <= (real)0.0 )
               {
                  const real Coeff_C = -dt_dh2*EigenVal[d][Mode];
                  const real Coeff_D = real(-4.0/3.0)*SQR(Coeff_C);

                  for (int v = 0; v < NWAVE; v++)
                     Coeff_L += LEigenVec[Mode][v]*(  Coeff_C*( dfc[ idx_wave[v] ] + dfc6[ idx_wave[v] ] ) +
                                                      Coeff_D*( dfc6[ idx_wave[v] ]                      )  );
                  for (int v = 0; v < NWAVE; v++)
                     Correct_L[ idx_wave[v] ] += Coeff_L*REigenVec[Mode][v];
               }

               if ( HLL_Include_All_Waves  ||  EigenVal[d][Mode] >= (real)0.0 )
               {
                  const real Coeff_A = -dt_dh2*EigenVal[d][Mode];
                  const real Coeff_B = real(-4.0/3.0)*SQR(Coeff_A);

                  for (int v = 0; v < NWAVE; v++)
                     Coeff_R += LEigenVec[Mode][v]*(  Coeff_A*( dfc[ idx_wave[v] ] - dfc6[ idx_wave[v] ] ) +
                                                      Coeff_B*( dfc6[ idx_wave[v] ]                      )  );
                  for (int v = 0; v < NWAVE; v++)
                     Correct_R[ idx_wave[v] ] += Coeff_R*REigenVec[Mode][v];
               }
            } // for (int Mode=0; Mode<NWAVE; Mode++)

#           else // Roe / exact solver
            Coeff_L = -dt_dh2*FMIN( EigenVal[d][       0 ], (real)0.0 );
            Coeff_R = -dt_dh2*FMAX( EigenVal[d][ NWAVE-1 ], (real)0.0 );

            for (int v = 0; v < NWAVE; v++)
            {
               Correct_L[ idx_wave[v] ] = Coeff_L*(  dfc[ idx_wave[v] ] + ( (real)1.0 - real(4.0/3.0)*Coeff_L )*dfc6[ idx_wave[v] ]  );
               Correct_R[ idx_wave[v] ] = Coeff_R*(  dfc[ idx_wave[v] ] - ( (real)1.0 + real(4.0/3.0)*Coeff_R )*dfc6[ idx_wave[v] ]  );
            }

            for (int Mode = 0; Mode < NWAVE; Mode++)
            {
               Coeff_L = (real)0.0;
               Coeff_R = (real)0.0;

               if ( EigenVal[d][Mode] <= (real)0.0 )
               {
                  const real Coeff_C = dt_dh2*( EigenVal[d][0] - EigenVal[d][Mode] );
                  const real Coeff_D = real(4.0/3.0)*dt_dh2*Coeff_C*( EigenVal[d][0] + EigenVal[d][Mode] );

                  for (int v = 0; v < NWAVE; v++)
                     Coeff_L += LEigenVec[Mode][v]*(  Coeff_C*( dfc[ idx_wave[v] ] + dfc6[ idx_wave[v] ] ) +
                                                      Coeff_D*( dfc6[ idx_wave[v] ]                      )  );
                  for (int v = 0; v < NWAVE; v++)
                     Correct_L[ idx_wave[v] ] += Coeff_L*REigenVec[Mode][v];
               }

               if ( EigenVal[d][Mode] >= (real)0.0 )
               {
                  const real Coeff_A = dt_dh2*( EigenVal[d][ NWAVE-1 ] - EigenVal[d][Mode] );
                  const real Coeff_B = real(4.0/3.0)*dt_dh2*Coeff_A*( EigenVal[d][ NWAVE-1 ] + EigenVal[d][Mode] );

                  for (int v = 0; v < NWAVE; v++)
                     Coeff_R += LEigenVec[Mode][v]*(  Coeff_A*( dfc[ idx_wave[v] ] - dfc6[ idx_wave[v] ] ) +
                                                      Coeff_B*( dfc6[ idx_wave[v] ]                      )  );
                  for (int v = 0; v < NWAVE; v++)
                     Correct_R[ idx_wave[v] ] += Coeff_R*REigenVec[Mode][v];
               }
            } // for (int Mode=0; Mode<NWAVE; Mode++)

#           endif // if (  ( RSOLVER == HLLE || RSOLVER == HLLC || RSOLVER == HLLD )  &&  defined HLL_NO_REF_STATE  ) ... else ...

#           if ( NCOMP_PASSIVE > 0 )
            Coeff_L = -dt_dh2*FMIN( EigenVal[d][1], (real)0.0 );
            Coeff_R = -dt_dh2*FMAX( EigenVal[d][1], (real)0.0 );
            for (int v = NCOMP_FLUID; v < NCOMP_TOTAL; v++)
            {
               Correct_L[v] = Coeff_L*(  dfc[v] + ( (real)1.0 - real(4.0/3.0)*Coeff_L )*dfc6[v]  );
               Correct_R[v] = Coeff_R*(  dfc[v] - ( (real)1.0 + real(4.0/3.0)*Coeff_R )*dfc6[v]  );
            }
#           endif

            Hydro_Rotate3D( Correct_L, d, false, MAG_OFFSET );
            Hydro_Rotate3D( Correct_R, d, false, MAG_OFFSET );

            for (int v = 0; v < NCOMP_LR; v++)
            {
               fcPri[faceL][v] += Correct_L[v];
               fcPri[faceR][v] += Correct_R[v];
            }

            fcPri[faceL][0] = FMAX( fcPri[faceL][0], MinDens );
            fcPri[faceR][0] = FMAX( fcPri[faceR][0], MinDens );
            fcPri[faceL][4] = Hydro_CheckMinPres( fcPri[faceL][4], MinPres );
            fcPri[faceR][4] = Hydro_CheckMinPres( fcPri[faceR][4], MinPres );

#           if ( NCOMP_PASSIVE > 0 )
            for (int v = NCOMP_FLUID; v < NCOMP_TOTAL; v++)
               if ( PassiveFloor & BIDX(v) ) {
               fcPri[faceL][v] = FMAX( fcPri[faceL][v], TINY_NUMBER );
               fcPri[faceR][v] = FMAX( fcPri[faceR][v], TINY_NUMBER ); }
#           endif

            Hydro_Pri2Con( fcPri[faceL], fcCon[faceL], FracPassive, NFrac, c_FracIdx,
                           EoS.DensPres2Eint_FuncPtr, EoS.Temp2HTilde_FuncPtr, EoS.HTilde2Temp_FuncPtr,
                           EoS.AuxArrayDevPtr_Flt, EoS.AuxArrayDevPtr_Int, EoS.Table, NULL );
            Hydro_Pri2Con( fcPri[faceR], fcCon[faceR], FracPassive, NFrac, c_FracIdx,
                           EoS.DensPres2Eint_FuncPtr, EoS.Temp2HTilde_FuncPtr, EoS.HTilde2Temp_FuncPtr,
                           EoS.AuxArrayDevPtr_Flt, EoS.AuxArrayDevPtr_Int, EoS.Table, NULL );

         } // for (int d=0; d<3; d++)

         for (int f = 0; f < 6; f++)
         for (int v = 0; v < NCOMP_TOTAL_PLUS_MAG; v++)
            g_FC_Var_1PG[f][v][idx_fc_out] = fcCon[f][v];

      } // CGPU_LOOP( idx_ij, N_FC_VAR2 )
      __syncthreads();
   } // for (int k_fc=0; k_fc<N_FC_VAR; k_fc++)
}
#endif // #ifndef MHD
#endif // __CUDACC__

#ifdef __CUDACC__
__global__
void CUFLU_FluidSolver_CTU(
   const real   g_Flu_Array_In [][NCOMP_TOTAL][ CUBE(FLU_NXT) ],
         real   g_Flu_Array_Out[][NCOMP_TOTAL][ CUBE(PS2) ],
   const real   g_Mag_Array_In [][NCOMP_MAG][ FLU_NXT_P1*SQR(FLU_NXT) ],
         real   g_Mag_Array_Out[][NCOMP_MAG][ PS2P1*SQR(PS2) ],
         char   g_DE_Array_Out [][ CUBE(PS2) ],
         real   g_Flux_Array   [][9][NCOMP_TOTAL][ SQR(PS2) ],
         real   g_Ele_Array    [][9][NCOMP_ELE][ PS2P1*PS2 ],
   const double g_Corner_Array [][3],
   const real   g_Pot_Array_USG[][ CUBE(USG_NXT_F) ],
         real   g_PriVar       []   [NCOMP_LR            ][ CUBE(FLU_NXT) ],
         real   g_Slope_PPM    [][3][NCOMP_LR            ][ CUBE(N_SLOPE_PPM) ],
         real   g_FC_Var       [][6][NCOMP_TOTAL_PLUS_MAG][ CUBE(N_FC_VAR) ],
         real   g_FC_Flux      [][3][NCOMP_TOTAL_PLUS_MAG][ CUBE(N_FC_FLUX) ],
         real   g_FC_Mag_Half  [][NCOMP_MAG][ FLU_NXT_P1*SQR(FLU_NXT) ],
         real   g_EC_Ele       [][NCOMP_MAG][ CUBE(N_EC_ELE) ],
   const real dt, const real dh,
   const bool StoreFlux, const bool StoreElectric,
   const LR_Limiter_t LR_Limiter, const real MinMod_Coeff, const double Time,
   const bool UsePot, const OptExtAcc_t ExtAcc, const ExtAcc_t ExtAcc_Func,
   const real MinDens, const real MinPres, const real MinEint,
   const real DualEnergySwitch,
   const long PassiveFloor,
   const bool NormPassive, const int NNorm,
   const bool FracPassive, const int NFrac,
   const bool JeansMinPres, const real JeansMinPres_Coeff,
   const EoS_t EoS )
#else
void CPU_FluidSolver_CTU(
   const real   g_Flu_Array_In [][NCOMP_TOTAL][ CUBE(FLU_NXT) ],
         real   g_Flu_Array_Out[][NCOMP_TOTAL][ CUBE(PS2) ],
   const real   g_Mag_Array_In [][NCOMP_MAG][ FLU_NXT_P1*SQR(FLU_NXT) ],
         real   g_Mag_Array_Out[][NCOMP_MAG][ PS2P1*SQR(PS2) ],
         char   g_DE_Array_Out [][ CUBE(PS2) ],
         real   g_Flux_Array   [][9][NCOMP_TOTAL][ SQR(PS2) ],
         real   g_Ele_Array    [][9][NCOMP_ELE][ PS2P1*PS2 ],
   const double g_Corner_Array [][3],
   const real   g_Pot_Array_USG[][ CUBE(USG_NXT_F) ],
         real   g_PriVar       []   [NCOMP_LR            ][ CUBE(FLU_NXT) ],
         real   g_Slope_PPM    [][3][NCOMP_LR            ][ CUBE(N_SLOPE_PPM) ],
         real   g_FC_Var       [][6][NCOMP_TOTAL_PLUS_MAG][ CUBE(N_FC_VAR) ],
         real   g_FC_Flux      [][3][NCOMP_TOTAL_PLUS_MAG][ CUBE(N_FC_FLUX) ],
         real   g_FC_Mag_Half  [][NCOMP_MAG][ FLU_NXT_P1*SQR(FLU_NXT) ],
         real   g_EC_Ele       [][NCOMP_MAG][ CUBE(N_EC_ELE) ],
   const int NPatchGroup,
   const real dt, const real dh,
   const bool StoreFlux, const bool StoreElectric,
   const LR_Limiter_t LR_Limiter, const real MinMod_Coeff, const double Time,
   const bool UsePot, const OptExtAcc_t ExtAcc, const ExtAcc_t ExtAcc_Func,
   const double c_ExtAcc_AuxArray[],
   const real MinDens, const real MinPres, const real MinEint,
   const real DualEnergySwitch,
   const long PassiveFloor,
   const bool NormPassive, const int NNorm, const int c_NormIdx[],
   const bool FracPassive, const int NFrac, const int c_FracIdx[],
   const bool JeansMinPres, const real JeansMinPres_Coeff,
   const EoS_t EoS )
#endif // #ifdef __CUDACC__ ... else ...
{

#  ifdef UNSPLIT_GRAVITY
   const bool CorrHalfVel          = true;
#  else
   const bool CorrHalfVel          = false;
#  endif
#  if !( defined __CUDACC__  &&  !defined MHD )
   const bool CorrHalfVel_No       = false;
#  endif
#  if !( defined __CUDACC__  &&  !defined MHD )
   const bool Con2Pri_Yes          = true;
#  endif
#  ifdef MHD
   const bool StoreElectric_No     = false;
#  endif
#  if ( defined __CUDACC__  &&  !defined GRAVITY )
   const double *c_ExtAcc_AuxArray = NULL;
#  endif


// openmp pragma for the CPU solver
#  ifndef __CUDACC__
#  pragma omp parallel
#  endif
   {
//    loop over all patch groups
//    --> CPU/GPU solver: use different (OpenMP threads) / (CUDA thread blocks)
//        to work on different patch groups
#     ifdef __CUDACC__
      const int P = blockIdx.x;
#     else
#     pragma omp for schedule( runtime )
      for (int P=0; P<NPatchGroup; P++)
#     endif
      {
//       0. point to the arrays associated with different patch groups
//          --> necessary because different patch groups are computed by different OpenMP threads or CUDA blocks in parallel
         real (*const g_FC_Var_1PG   )[NCOMP_TOTAL_PLUS_MAG][ CUBE(N_FC_VAR)    ] = g_FC_Var   [P];
         real (*const g_FC_Flux_1PG  )[NCOMP_TOTAL_PLUS_MAG][ CUBE(N_FC_FLUX)   ] = g_FC_Flux  [P];
#        if !( defined __CUDACC__  &&  !defined MHD )
         real (*const g_PriVar_1PG   )                      [ CUBE(FLU_NXT)     ] = g_PriVar   [P];
#        endif
#        if !( defined __CUDACC__  &&  !defined MHD )
         real (*const g_Slope_PPM_1PG)[NCOMP_LR            ][ CUBE(N_SLOPE_PPM) ] = g_Slope_PPM[P];
#        endif

#        ifdef MHD
         real (*const g_FC_Mag_Half_1PG)[ FLU_NXT_P1*SQR(FLU_NXT) ] = g_FC_Mag_Half[P];
         real (*const g_EC_Ele_1PG     )[ CUBE(N_EC_ELE)          ] = g_EC_Ele     [P];
         real (*const g_PriVar_Half_1PG)[ CUBE(FLU_NXT) ]           = g_PriVar_1PG;
#        elif !defined __CUDACC__
         real (*const g_FC_Mag_Half_1PG)[ FLU_NXT_P1*SQR(FLU_NXT) ] = NULL;
         real (*const g_EC_Ele_1PG     )[ CUBE(N_EC_ELE)          ] = NULL;
#        endif


//       1. evaluate the face-centered values at the half time-step
//          GPU non-MHD: CUFLU_DR_FCVar has already filled g_FC_Var; skip here
//          CPU or GPU+MHD: full inline
#        if !( defined __CUDACC__  &&  !defined MHD )
         Hydro_DataReconstruction( g_Flu_Array_In[P], g_Mag_Array_In[P], g_PriVar_1PG, g_FC_Var_1PG, g_Slope_PPM_1PG,
                                   NULL, Con2Pri_Yes, LR_Limiter, MinMod_Coeff, dt, dh,
                                   MinDens, MinPres, MinEint, PassiveFloor, FracPassive, NFrac, c_FracIdx,
                                   JeansMinPres, JeansMinPres_Coeff, &EoS );
#        endif


//       2-4. GPU non-MHD: k-slab pipeline — half-step Riemanns into shmem ping-pong buffers,
//            inline TGrad from shmem; eliminates ~50 MB DRAM roundtrip of g_FC_Flux_1PG[half-step]
#        if ( defined __CUDACC__  &&  !defined MHD )
         {
            const int  N_FC_VAR2  = SQR(N_FC_VAR);
            const real dt_dh2     = (real)0.5*dt/dh;
            const int  didx_FC[3] = { 1, N_FC_VAR, N_FC_VAR2 };

//          ping-pong shmem for half-step fluxes; pp = current k-plane, pq = previous
            __shared__ real s_FluxX[2][NCOMP_TOTAL][SQR(N_FC_VAR)];
            __shared__ real s_FluxY[2][NCOMP_TOTAL][SQR(N_FC_VAR)];
            __shared__ real s_FluxZ[2][NCOMP_TOTAL][SQR(N_FC_VAR)];

            for (int k_fc = 0; k_fc < N_FC_VAR; k_fc++)
            {
               const int pp = k_fc & 1;
               const int pq = 1 - pp;

//             step 1: half-step Riemann for d=0 (x), d=1 (y), d=2 (z, skipped at k_fc=0)
               for (int d_r = 0; d_r < 3; d_r++)
               {
                  if ( d_r == 2  &&  k_fc == 0 )   continue;

                  const int faceR    = 2*d_r + 1;    // right face of L cell
                  const int faceL    = 2*d_r;         // left  face of R cell
                  const int n_flux_i = (d_r == 0) ? N_FC_VAR-1 : N_FC_VAR;
                  const int n_flux_j = (d_r == 1) ? N_FC_VAR-1 : N_FC_VAR;

                  CGPU_LOOP( idx, n_flux_i*n_flux_j )
                  {
                     const int i_flux = idx % n_flux_i;
                     const int j_flux = idx / n_flux_i;

                     int idx_Lfc, idx_Rfc;
                     if ( d_r < 2 ) {
                        idx_Lfc = IDX321( i_flux, j_flux, k_fc,   N_FC_VAR, N_FC_VAR );
                        idx_Rfc = idx_Lfc + didx_FC[d_r];
                     } else {
                        idx_Lfc = IDX321( i_flux, j_flux, k_fc-1, N_FC_VAR, N_FC_VAR );
                        idx_Rfc = IDX321( i_flux, j_flux, k_fc,   N_FC_VAR, N_FC_VAR );
                     }

                     real ConVar_L[NCOMP_TOTAL_PLUS_MAG], ConVar_R[NCOMP_TOTAL_PLUS_MAG], Flux_1F[NCOMP_TOTAL_PLUS_MAG];
                     for (int v=0; v<NCOMP_TOTAL_PLUS_MAG; v++) {
                        ConVar_L[v] = g_FC_Var_1PG[faceR][v][idx_Lfc];
                        ConVar_R[v] = g_FC_Var_1PG[faceL][v][idx_Rfc];
                     }

#                    if   ( RSOLVER == EXACT )
                     Hydro_RiemannSolver_Exact( d_r, Flux_1F, ConVar_L, ConVar_R, MinDens, MinPres, PassiveFloor,
                                                EoS.DensEint2Pres_FuncPtr, EoS.DensPres2CSqr_FuncPtr,
                                                EoS.AuxArrayDevPtr_Flt, EoS.AuxArrayDevPtr_Int, EoS.Table );
#                    elif ( RSOLVER == ROE )
                     Hydro_RiemannSolver_Roe( d_r, Flux_1F, ConVar_L, ConVar_R, MinDens, MinPres, PassiveFloor,
                                              EoS.DensEint2Pres_FuncPtr, EoS.DensPres2CSqr_FuncPtr,
                                              EoS.AuxArrayDevPtr_Flt, EoS.AuxArrayDevPtr_Int, EoS.Table );
#                    elif ( RSOLVER == HLLE )
                     Hydro_RiemannSolver_HLLE( d_r, Flux_1F, ConVar_L, ConVar_R, MinDens, MinPres, PassiveFloor,
                                               EoS.DensEint2Pres_FuncPtr, EoS.DensPres2CSqr_FuncPtr,
                                               EoS.GuessHTilde_FuncPtr, EoS.HTilde2Temp_FuncPtr,
                                               EoS.AuxArrayDevPtr_Flt, EoS.AuxArrayDevPtr_Int, EoS.Table );
#                    elif ( RSOLVER == HLLC )
                     Hydro_RiemannSolver_HLLC( d_r, Flux_1F, ConVar_L, ConVar_R, MinDens, MinPres, PassiveFloor,
                                               EoS.DensEint2Pres_FuncPtr, EoS.DensPres2CSqr_FuncPtr,
                                               EoS.GuessHTilde_FuncPtr, EoS.HTilde2Temp_FuncPtr,
                                               EoS.AuxArrayDevPtr_Flt, EoS.AuxArrayDevPtr_Int, EoS.Table );
#                    else
#                    error : unsupported Riemann solver in k-slab path
#                    endif

#                    if ( RSOLVER_RESCUE != OPTION_NONE )
                     for (int v=0; v<NCOMP_TOTAL_PLUS_MAG; v++) {
                        if ( Flux_1F[v] != Flux_1F[v] ) {
#                          if   ( RSOLVER_RESCUE == EXACT )
                           Hydro_RiemannSolver_Exact( d_r, Flux_1F, ConVar_L, ConVar_R, MinDens, MinPres, PassiveFloor,
                                                      EoS.DensEint2Pres_FuncPtr, EoS.DensPres2CSqr_FuncPtr,
                                                      EoS.AuxArrayDevPtr_Flt, EoS.AuxArrayDevPtr_Int, EoS.Table );
#                          elif ( RSOLVER_RESCUE == ROE )
                           Hydro_RiemannSolver_Roe( d_r, Flux_1F, ConVar_L, ConVar_R, MinDens, MinPres, PassiveFloor,
                                                    EoS.DensEint2Pres_FuncPtr, EoS.DensPres2CSqr_FuncPtr,
                                                    EoS.AuxArrayDevPtr_Flt, EoS.AuxArrayDevPtr_Int, EoS.Table );
#                          elif ( RSOLVER_RESCUE == HLLE )
                           Hydro_RiemannSolver_HLLE( d_r, Flux_1F, ConVar_L, ConVar_R, MinDens, MinPres, PassiveFloor,
                                                     EoS.DensEint2Pres_FuncPtr, EoS.DensPres2CSqr_FuncPtr,
                                                     EoS.GuessHTilde_FuncPtr, EoS.HTilde2Temp_FuncPtr,
                                                     EoS.AuxArrayDevPtr_Flt, EoS.AuxArrayDevPtr_Int, EoS.Table );
#                          elif ( RSOLVER_RESCUE == HLLC )
                           Hydro_RiemannSolver_HLLC( d_r, Flux_1F, ConVar_L, ConVar_R, MinDens, MinPres, PassiveFloor,
                                                     EoS.DensEint2Pres_FuncPtr, EoS.DensPres2CSqr_FuncPtr,
                                                     EoS.GuessHTilde_FuncPtr, EoS.HTilde2Temp_FuncPtr,
                                                     EoS.AuxArrayDevPtr_Flt, EoS.AuxArrayDevPtr_Int, EoS.Table );
#                          endif
                           break;
                        }
                     }
#                    endif // RSOLVER_RESCUE

                     const int ij = j_flux*N_FC_VAR + i_flux;
                     for (int v=0; v<NCOMP_TOTAL; v++) {
                        if      ( d_r == 0 ) s_FluxX[pp][v][ij] = Flux_1F[v];
                        else if ( d_r == 1 ) s_FluxY[pp][v][ij] = Flux_1F[v];
                        else                 s_FluxZ[pp][v][ij] = Flux_1F[v];
                     }
                  } // CGPU_LOOP Riemann d_r
               } // for d_r = 0..2

               __syncthreads();

//             step 2: TGrad at k_out = k_fc-1 (deferred by 1 to allow z-flux at k_fc-1/k_fc boundary)
//             convention: s_FluxX/Y[pq] = x/y-fluxes at k_out; s_FluxZ[pp]/[pq] = z-fluxes at interfaces k_out/k_out+1 and k_out-1/k_out
               if ( k_fc >= 1 )
               {
                  const int k_out = k_fc - 1;

                  for (int d=0; d<3; d++)
                  {
                     const int faceL = 2*d;
                     const int faceR = faceL + 1;

                     int nskip_i, nskip_j, nskip_k;
                     switch (d) {
                        case 0:  nskip_i=0; nskip_j=1; nskip_k=1;  break;
                        case 1:  nskip_i=1; nskip_j=0; nskip_k=1;  break;
                        default: nskip_i=1; nskip_j=1; nskip_k=0;  break;
                     }

                     if ( k_out < nskip_k  ||  k_out >= N_FC_VAR - nskip_k )   continue;

                     const int size_i = N_FC_VAR - 2*nskip_i;
                     const int size_j = N_FC_VAR - 2*nskip_j;

                     CGPU_LOOP( idx_ij, size_i*size_j )
                     {
                        const int i_fc   = idx_ij % size_i + nskip_i;
                        const int j_fc   = idx_ij / size_i + nskip_j;
                        const int idx_fc = IDX321( i_fc, j_fc, k_out, N_FC_VAR, N_FC_VAR );
                        const int ij     = j_fc*N_FC_VAR + i_fc;

                        real fc[2][NCOMP_TOTAL_PLUS_MAG];
                        for (int v=0; v<NCOMP_TOTAL_PLUS_MAG; v++) {
                           fc[0][v] = g_FC_Var_1PG[faceL][v][idx_fc];
                           fc[1][v] = g_FC_Var_1PG[faceR][v][idx_fc];
                        }

                        for (int v=0; v<NCOMP_TOTAL; v++) {
                           real TGrad1, TGrad2;
                           switch (d) {
                           case 0:   // TDir1=y, TDir2=z
                              TGrad1 = s_FluxY[pq][v][ij]                         - s_FluxY[pq][v][(j_fc-1)*N_FC_VAR+i_fc];
                              TGrad2 = s_FluxZ[pp][v][ij]                         - s_FluxZ[pq][v][ij];
                              break;
                           case 1:   // TDir1=z, TDir2=x
                              TGrad1 = s_FluxZ[pp][v][ij]                         - s_FluxZ[pq][v][ij];
                              TGrad2 = s_FluxX[pq][v][ij]                         - s_FluxX[pq][v][j_fc*N_FC_VAR+(i_fc-1)];
                              break;
                           default:  // TDir1=x, TDir2=y
                              TGrad1 = s_FluxX[pq][v][ij]                         - s_FluxX[pq][v][j_fc*N_FC_VAR+(i_fc-1)];
                              TGrad2 = s_FluxY[pq][v][ij]                         - s_FluxY[pq][v][(j_fc-1)*N_FC_VAR+i_fc];
                              break;
                           }
                           const real Correct = -dt_dh2 * ( TGrad1 + TGrad2 );
                           fc[0][v] += Correct;
                           fc[1][v] += Correct;
                        }

                        for (int f=0; f<2; f++) {
                           fc[f][0] = FMAX( fc[f][0], MinDens );
                           fc[f][4] = Hydro_CheckMinEintInEngy( fc[f][0], fc[f][1], fc[f][2], fc[f][3],
                                                                 fc[f][4], MinEint, PassiveFloor, NULL_REAL );
#                          if ( NCOMP_PASSIVE > 0 )
                           for (int v=NCOMP_FLUID; v<NCOMP_TOTAL; v++)
                              if ( PassiveFloor & BIDX(v) )
                                 fc[f][v] = FMAX( fc[f][v], TINY_NUMBER );
#                          endif
                        }

                        for (int v=0; v<NCOMP_TOTAL_PLUS_MAG; v++) {
                           g_FC_Var_1PG[faceL][v][idx_fc] = fc[0][v];
                           g_FC_Var_1PG[faceR][v][idx_fc] = fc[1][v];
                        }
                     } // CGPU_LOOP TGrad
                  } // for d = 0..2
               } // if k_fc >= 1

               __syncthreads();

            } // for k_fc = 0..N_FC_VAR-1

//          post-loop: TGrad d=2 at k_out = N_FC_VAR-1 (d=0,1 skip due to nskip_k=1)
            {
               const int pp_last = (N_FC_VAR-1) & 1;   // = 1 for N_FC_VAR=18
               const int k_out   = N_FC_VAR - 1;

               CGPU_LOOP( idx_ij, SQR(N_FC_VAR-2) )
               {
                  const int i_fc   = idx_ij % (N_FC_VAR-2) + 1;
                  const int j_fc   = idx_ij / (N_FC_VAR-2) + 1;
                  const int idx_fc = IDX321( i_fc, j_fc, k_out, N_FC_VAR, N_FC_VAR );
                  const int ij     = j_fc*N_FC_VAR + i_fc;

                  real fc[2][NCOMP_TOTAL_PLUS_MAG];
                  for (int v=0; v<NCOMP_TOTAL_PLUS_MAG; v++) {
                     fc[0][v] = g_FC_Var_1PG[4][v][idx_fc];
                     fc[1][v] = g_FC_Var_1PG[5][v][idx_fc];
                  }

                  for (int v=0; v<NCOMP_TOTAL; v++) {
                     const real TGrad1 = s_FluxX[pp_last][v][ij] - s_FluxX[pp_last][v][j_fc*N_FC_VAR+(i_fc-1)];
                     const real TGrad2 = s_FluxY[pp_last][v][ij] - s_FluxY[pp_last][v][(j_fc-1)*N_FC_VAR+i_fc];
                     const real Correct = -dt_dh2 * ( TGrad1 + TGrad2 );
                     fc[0][v] += Correct;
                     fc[1][v] += Correct;
                  }

                  for (int f=0; f<2; f++) {
                     fc[f][0] = FMAX( fc[f][0], MinDens );
                     fc[f][4] = Hydro_CheckMinEintInEngy( fc[f][0], fc[f][1], fc[f][2], fc[f][3],
                                                           fc[f][4], MinEint, PassiveFloor, NULL_REAL );
#                    if ( NCOMP_PASSIVE > 0 )
                     for (int v=NCOMP_FLUID; v<NCOMP_TOTAL; v++)
                        if ( PassiveFloor & BIDX(v) )
                           fc[f][v] = FMAX( fc[f][v], TINY_NUMBER );
#                    endif
                  }

                  for (int v=0; v<NCOMP_TOTAL_PLUS_MAG; v++) {
                     g_FC_Var_1PG[4][v][idx_fc] = fc[0][v];
                     g_FC_Var_1PG[5][v][idx_fc] = fc[1][v];
                  }
               }
               __syncthreads();
            } // post-loop TGrad d=2

         } // k-slab block

#        else // !( defined __CUDACC__  &&  !defined MHD ): CPU or MHD — standard path

//       2. evaluate the face-centered half-step fluxes by solving the Riemann problem
         Hydro_ComputeFlux( g_FC_Var_1PG, g_FC_Flux_1PG, N_HF_FLUX, 0, 0, CorrHalfVel_No,
                            NULL, NULL, NULL_REAL, NULL_REAL, NULL_REAL,
                            EXT_POT_NONE, EXT_ACC_NONE, NULL, NULL,
                            MinDens, MinPres, PassiveFloor, &EoS );


//       3. evaluate electric field and update B field at the half time-step
#        ifdef MHD
         MHD_ComputeElectric( g_EC_Ele_1PG, g_FC_Flux_1PG, g_PriVar_1PG, N_HF_ELE, N_HF_FLUX,
                              FLU_NXT, LR_GHOST_SIZE, dt, dh, StoreElectric_No, NULL,
                              CorrHalfVel_No, NULL, NULL, NULL_REAL,
                              EXT_POT_NONE, EXT_ACC_NONE, NULL, NULL );

         MHD_UpdateMagnetic( g_FC_Mag_Half_1PG[0], g_FC_Mag_Half_1PG[1], g_FC_Mag_Half_1PG[2],
                             g_Mag_Array_In[P], g_EC_Ele_1PG, (real)0.5*dt, dh, N_HF_VAR, N_HF_ELE, FLU_GHOST_SIZE-1 );
#        endif


//       4. correct the face-centered variables by the transverse flux gradients
         Hydro_TGradientCorrection( g_FC_Var_1PG, g_FC_Flux_1PG, g_Mag_Array_In[P], g_FC_Mag_Half_1PG, g_EC_Ele_1PG, g_PriVar_1PG,
                                    dt, dh, MinDens, MinEint, PassiveFloor );

#        endif // #if ( defined __CUDACC__  &&  !defined MHD )


//       5. evaluate the cell-centered primitive variables at the half time-step
//          --> for computing CT electric field later
#        ifdef MHD
         MHD_HalfStepPrimitive( g_Flu_Array_In[P], g_FC_Mag_Half_1PG, g_PriVar_Half_1PG, g_FC_Flux_1PG, dt, dh, MinDens );
#        endif


//       6. evaluate the face-centered full-step fluxes by solving the Riemann problem with the corrected data
#        ifdef MHD
         const int NSkip_N = 1;
         const int NSkip_T = 1;
#        else
         const int NSkip_N = 0;
         const int NSkip_T = 1;
#        endif
         Hydro_ComputeFlux( g_FC_Var_1PG, g_FC_Flux_1PG, N_FL_FLUX, NSkip_N, NSkip_T, CorrHalfVel,
                            g_Pot_Array_USG[P], g_Corner_Array[P], dt, dh, Time,
                            UsePot, ExtAcc, ExtAcc_Func, c_ExtAcc_AuxArray,
                            MinDens, MinPres, PassiveFloor, &EoS );

         if ( StoreFlux )
            Hydro_StoreIntFlux( g_FC_Flux_1PG, g_Flux_Array[P], N_FL_FLUX );


//       7. evaluate electric field and update B field at the full time-step
//          --> must update B field before Hydro_FullStepUpdate() since the latter requires
//              the updated magnetic energy when adopting the dual-energy formalism
#        ifdef MHD
         MHD_ComputeElectric( g_EC_Ele_1PG, g_FC_Flux_1PG, g_PriVar_Half_1PG, N_FL_ELE, N_FL_FLUX,
                              N_HF_VAR, 0, dt, dh, StoreElectric, g_Ele_Array[P],
                              CorrHalfVel, g_Pot_Array_USG[P], g_Corner_Array[P], Time,
                              UsePot, ExtAcc, ExtAcc_Func, c_ExtAcc_AuxArray );

         MHD_UpdateMagnetic( g_Mag_Array_Out[P][0], g_Mag_Array_Out[P][1], g_Mag_Array_Out[P][2],
                             g_Mag_Array_In[P], g_EC_Ele_1PG, dt, dh, PS2, N_FL_ELE, FLU_GHOST_SIZE );
#        endif


//       8. full-step evolution of the fluid data
//          --> CTU does not support reducing the min-mod coefficient
         Hydro_FullStepUpdate( g_Flu_Array_In[P], g_Flu_Array_Out[P], g_DE_Array_Out[P], g_Mag_Array_Out[P],
                               g_FC_Flux_1PG, dt, dh, MinDens, MinEint, DualEnergySwitch,
                               PassiveFloor, NormPassive, NNorm, c_NormIdx, &EoS, NULL, NULL_INT, NULL_INT );

      } // loop over all patch groups
   } // OpenMP parallel region

} // FUNCTION : CPU_FluidSolver_CTU



//-------------------------------------------------------------------------------------------------------
// Function    :  Hydro_TGradientCorrection
// Description :  Correct the face-centered variables by the transverse flux gradients
//
// Note        :  1. Ref: (a) Stone et al., ApJS, 178, 137 (2008)
//                        (b) Gardiner & Stone, J. Comput. Phys., 227, 4123 (2008)
//                2. Assuming "N_FC_VAR == N_HF_FLUX"
//
// Parameter   :  g_FC_Var     : Array to store the input and output face-centered conserved variables
//                               --> Accessed with the stride N_FC_VAR
//                g_FC_Flux    : Array storing the input face-centered fluxes
//                               --> Accessed with the stride N_HF_FLUX
//                g_FC_B_In    : Array storing the input initial   face-centered B field
//                g_FC_B_Half  : Array storing the input half-step face-centered B field
//                g_EC_Ele     : Array storing the input edge-centered electric field
//                g_PriVar     : Array storing the input cell-centered primitive variables
//                dt           : Time interval to advance solution
//                dh           : Cell size
//                MinDens/Eint : Density and internal energy floors
//                PassiveFloor : Bitwise flag to specify the passive scalars to be floored
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE
void Hydro_TGradientCorrection(       real g_FC_Var   [][NCOMP_TOTAL_PLUS_MAG][ CUBE(N_FC_VAR)  ],
                                const real g_FC_Flux  [][NCOMP_TOTAL_PLUS_MAG][ CUBE(N_FC_FLUX) ],
                                const real g_FC_B_In  [][ FLU_NXT_P1*SQR(FLU_NXT) ],
                                const real g_FC_B_Half[][ FLU_NXT_P1*SQR(FLU_NXT) ],
                                const real g_EC_Ele   [][ CUBE(N_EC_ELE) ],
                                const real g_PriVar   [][ CUBE(FLU_NXT) ],
                                const real dt, const real dh, const real MinDens, const real MinEint,
                                const long PassiveFloor )
{

   const int  didx_flux[3]   = { 1, N_HF_FLUX, SQR(N_HF_FLUX) };
   const real dt_dh2         = (real)0.5*dt/dh;
#  ifdef MHD
   const int  didx_b_in  [3] = { 1, FLU_NXT,  SQR(FLU_NXT)  };
   const int  didx_b_half[3] = { 1, N_HF_VAR, SQR(N_HF_VAR) };
   const int  didx_ele   [3] = { 1, N_HF_ELE, SQR(N_HF_ELE) };
   const real _dh            = (real)1.0/dh;
   const real dt_dh4         = (real)0.25*dt*_dh;
   const real dt_2           = (real)0.5*dt;

   real PriVar_1Cell[NCOMP_FLUID+NCOMP_MAG], B_Face[NCOMP_MAG][2];   // [2]=left/right faces
#  endif

// loop over different spatial directions
   for (int d=0; d<3; d++)
   {
      const int faceL = 2*d;
      const int faceR = faceL+1;
      const int TDir1 = (d+1)%3;    // transverse direction 1
      const int TDir2 = (d+2)%3;    // transverse direction 2

      real fc_var[2][NCOMP_TOTAL_PLUS_MAG];  // [2]=left/right faces

#     ifdef MHD
      const int nskip[3] = { 1, 1, 1 };
#     else
      int nskip[3];
      switch ( d )
      {
         case 0 : nskip[0] = 0;  nskip[1] = 1;  nskip[2] = 1;  break;
         case 1 : nskip[0] = 1;  nskip[1] = 0;  nskip[2] = 1;  break;
         case 2 : nskip[0] = 1;  nskip[1] = 1;  nskip[2] = 0;  break;
      }
#     endif

      const int size_i  = ( N_FC_VAR - 2*nskip[0] );
      const int size_j  = ( N_FC_VAR - 2*nskip[1] );
      const int size_k  = ( N_FC_VAR - 2*nskip[2] );
      const int size_ij = size_i*size_j;

      CGPU_LOOP( idx0, size_i*size_j*size_k )
      {
//       i/j/k0 start from zero
         const int i0         = idx0 % size_i;
         const int j0         = idx0 % size_ij / size_i;
         const int k0         = idx0 / size_ij;

         const int i_fc_var   = i0 + nskip[0];
         const int j_fc_var   = j0 + nskip[1];
         const int k_fc_var   = k0 + nskip[2];
         const int idx_fc_var = IDX321( i_fc_var, j_fc_var, k_fc_var, N_FC_VAR, N_FC_VAR );

         const int idx_fluxR  = idx_fc_var;  // assuming N_FC_VAR == N_HF_FLUX
         const int idx_fluxL1 = idx_fluxR - didx_flux[TDir1];
         const int idx_fluxL2 = idx_fluxR - didx_flux[TDir2];

//       0. load g_FC_Var[] to the local variable fc[] to reduce the GPU global memory access
         for (int v=0; v<NCOMP_TOTAL_PLUS_MAG; v++)
         {
            fc_var[0][v] = g_FC_Var[faceL][v][idx_fc_var];
            fc_var[1][v] = g_FC_Var[faceR][v][idx_fc_var];
         }


//       1. calculate the transverse fluid flux gradients and update the corresponding face-centered fluid variables
         for (int v=0; v<NCOMP_TOTAL; v++)
         {
            real Correct, TGrad1, TGrad2;

            TGrad1  = g_FC_Flux[TDir1][v][idx_fluxR] - g_FC_Flux[TDir1][v][idx_fluxL1];
            TGrad2  = g_FC_Flux[TDir2][v][idx_fluxR] - g_FC_Flux[TDir2][v][idx_fluxL2];
            Correct = -dt_dh2*( TGrad1 + TGrad2 );

            fc_var[0][v] += Correct;
            fc_var[1][v] += Correct;
         }


#        ifdef MHD
//       2. correct the transverse B field
         const int idx_ele = IDX321( i0, j0, k0, N_HF_ELE, N_HF_ELE );

         for (int v=1; v<NCOMP_MAG; v++)
         {
            real Correct, TGrad1, TGrad2, Sign;

            const int TD1 = (d+v)%3;            // transverse direction 1
            const int TD2 = (d+2*v)%3;          // transverse direction 2
            const int TB  = TD1 + MAG_OFFSET;   // target transverse B field

            TGrad1  = g_EC_Ele[d][ idx_ele + didx_ele[TD2]                 ] - g_EC_Ele[d][ idx_ele                 ];
            TGrad2  = g_EC_Ele[d][ idx_ele + didx_ele[TD2] + didx_ele[TD1] ] - g_EC_Ele[d][ idx_ele + didx_ele[TD1] ];
            Sign    = (real)2.0*v - (real)3.0;  // v=1/2 --> sign=-1/+1
            Correct = Sign*dt_dh4*( TGrad2 + TGrad1 );

            fc_var[0][TB] += Correct;
            fc_var[1][TB] += Correct;
         } // for (int v=1; v<NCOMP_MAG; v++)


//       3. add the divergence(B) source terms
         Hydro_Rotate3D( fc_var[0], d, true, MAG_OFFSET );
         Hydro_Rotate3D( fc_var[1], d, true, MAG_OFFSET );

//       3-1. get the initial cell-centered primitive variables
         const int i_pri   = i_fc_var + LR_GHOST_SIZE;
         const int j_pri   = j_fc_var + LR_GHOST_SIZE;
         const int k_pri   = k_fc_var + LR_GHOST_SIZE;
         const int idx_pri = IDX321( i_pri, j_pri, k_pri, FLU_NXT, FLU_NXT );

//       skip passive scalars
         for (int v=0; v<NCOMP_FLUID; v++)   PriVar_1Cell[ v               ] = g_PriVar[ v              ][idx_pri];
         for (int v=0; v<NCOMP_MAG;   v++)   PriVar_1Cell[ v + NCOMP_FLUID ] = g_PriVar[ v + MAG_OFFSET ][idx_pri];

         Hydro_Rotate3D( PriVar_1Cell, d, true, NCOMP_FLUID );

//       3-2. get the initial face-centered B field
         const int idx_b_in[3] = { IDX321( i_pri, j_pri, k_pri, FLU_NXT_P1, FLU_NXT    ),
                                   IDX321( i_pri, j_pri, k_pri, FLU_NXT,    FLU_NXT_P1 ),
                                   IDX321( i_pri, j_pri, k_pri, FLU_NXT,    FLU_NXT    ) };

         B_Face[0][0] = g_FC_B_In[  d  ][ idx_b_in[  d  ]                    ];
         B_Face[0][1] = g_FC_B_In[  d  ][ idx_b_in[  d  ] + didx_b_in[  d  ] ];

         B_Face[1][0] = g_FC_B_In[TDir1][ idx_b_in[TDir1]                    ];
         B_Face[1][1] = g_FC_B_In[TDir1][ idx_b_in[TDir1] + didx_b_in[TDir1] ];

         B_Face[2][0] = g_FC_B_In[TDir2][ idx_b_in[TDir2]                    ];
         B_Face[2][1] = g_FC_B_In[TDir2][ idx_b_in[TDir2] + didx_b_in[TDir2] ];


//       3-3. add the divergence(B) source term
         const real Vy = PriVar_1Cell[ 2 ];
         const real Vz = PriVar_1Cell[ 3 ];
         const real Bx = PriVar_1Cell[ 0 + NCOMP_FLUID ];
         const real By = PriVar_1Cell[ 1 + NCOMP_FLUID ];
         const real Bz = PriVar_1Cell[ 2 + NCOMP_FLUID ];

         real dB[NCOMP_MAG], SrcFlu[NCOMP_FLUID-1], SrcMag[2], Vy_MinModBxz, Vz_MinModBxy;

         for (int v=0; v<NCOMP_MAG; v++)  dB[v] = ( B_Face[v][1] - B_Face[v][0] )*_dh;

#        define MINMOD( a , b )  (  ( (a)*(b)>(real)0.0 ) ? ( SIGN(a)*FMIN(FABS(a),FABS(b)) ) : (real)0.0  )
         Vy_MinModBxz = Vy*MINMOD( -dB[2], dB[0] );
         Vz_MinModBxy = Vz*MINMOD( -dB[1], dB[0] );
#        undef MINMOD

         SrcFlu[0] = dt_2*Bx*dB[0];
         SrcFlu[1] = dt_2*By*dB[0];
         SrcFlu[2] = dt_2*Bz*dB[0];
         SrcFlu[3] = dt_2*( By*Vy_MinModBxz + Bz*Vz_MinModBxy );

         SrcMag[0] = dt_2*Vy_MinModBxz;
         SrcMag[1] = dt_2*Vz_MinModBxy;

         for (int f=0; f<2; f++)
         {
            fc_var[f][ 1 ] += SrcFlu[0];
            fc_var[f][ 2 ] += SrcFlu[1];
            fc_var[f][ 3 ] += SrcFlu[2];
            fc_var[f][ 4 ] += SrcFlu[3];

            fc_var[f][ 1 + MAG_OFFSET ] += SrcMag[0];
            fc_var[f][ 2 + MAG_OFFSET ] += SrcMag[1];
         }


//       4. set the longitudinal B field to the half-step values updated by MHD_UpdateMagnetic()
         int idx_b_half;

         switch ( d )
         {
            case 0 : idx_b_half = IDX321( i0, j0, k0, N_HF_VAR+1, N_HF_VAR   );  break;
            case 1 : idx_b_half = IDX321( i0, j0, k0, N_HF_VAR,   N_HF_VAR+1 );  break;
            case 2 : idx_b_half = IDX321( i0, j0, k0, N_HF_VAR,   N_HF_VAR   );  break;
         }

         fc_var[0][MAG_OFFSET] = g_FC_B_Half[d][ idx_b_half                  ];
         fc_var[1][MAG_OFFSET] = g_FC_B_Half[d][ idx_b_half + didx_b_half[d] ];

         Hydro_Rotate3D( fc_var[0], d, false, MAG_OFFSET );
         Hydro_Rotate3D( fc_var[1], d, false, MAG_OFFSET );
#        endif // #ifdef MHD


//       5. apply density and internal energy floors
         for (int f=0; f<2; f++)
         {
#           ifdef MHD
            const real Bx   = fc_var[f][ MAG_OFFSET + 0 ];
            const real By   = fc_var[f][ MAG_OFFSET + 1 ];
            const real Bz   = fc_var[f][ MAG_OFFSET + 2 ];
            const real Emag= (real)0.5*( SQR(Bx) + SQR(By) + SQR(Bz) );
#           else
            const real Emag = NULL_REAL;
#           endif
            fc_var[f][0] = FMAX( fc_var[f][0], MinDens );
            fc_var[f][4] = Hydro_CheckMinEintInEngy( fc_var[f][0], fc_var[f][1], fc_var[f][2], fc_var[f][3], fc_var[f][4],
                                                     MinEint, PassiveFloor, Emag );
#           if ( NCOMP_PASSIVE > 0 )
            for (int v=NCOMP_FLUID; v<NCOMP_TOTAL; v++)
            if ( PassiveFloor & BIDX(v) )
            fc_var[f][v] = FMAX( fc_var[f][v], TINY_NUMBER );
#           endif
         }

//       store the results to g_FC_Var[]
         for (int v=0; v<NCOMP_TOTAL_PLUS_MAG; v++)
         {
            g_FC_Var[faceL][v][idx_fc_var] = fc_var[0][v];
            g_FC_Var[faceR][v][idx_fc_var] = fc_var[1][v];
         }

      } // CGPU_LOOP( idx0, size_i*size_j*size_k )
   } // for (int d=0; d<3; d++)


#  ifdef __CUDACC__
   __syncthreads();
#  endif

} // FUNCTION : Hydro_TGradientCorrection



#endif // #if ( MODEL == HYDRO  &&  FLU_SCHEME == CTU )
