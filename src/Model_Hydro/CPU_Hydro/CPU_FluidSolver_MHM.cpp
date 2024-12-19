#include "CUFLU.h"

#if (  MODEL == HYDRO  &&  ( FLU_SCHEME == MHM || FLU_SCHEME == MHM_RP )  )



// external functions
#ifdef __CUDACC__

#include "CUFLU_Shared_FluUtility.cu"
#include "CUFLU_Shared_DataReconstruction.cu"
#include "CUFLU_Shared_ComputeFlux.cu"
#include "CUFLU_Shared_FullStepUpdate.cu"
#ifdef MHD
#include "CUFLU_Shared_ConstrainedTransport.cu"
#endif

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

#include "CUDA_ConstMemory.h"

#ifdef COSMIC_RAY
# include "CUFLU_CosmicRay.cu"
#ifdef CR_DIFFUSION
# include "../../Microphysics/CosmicRayDiffusion/CUFLU_CR_AddDiffuseFlux.cu"
#endif
#endif // #ifdef COSMIC_RAY

#else // #ifdef __CUDACC__

void Hydro_DataReconstruction( const real g_ConVar   [][ CUBE(FLU_NXT) ],
                               const real g_FC_B     [][ SQR(FLU_NXT)*FLU_NXT_P1 ],
                                     real g_PriVar   [][ CUBE(FLU_NXT) ],
                                     real g_FC_Var   [][NCOMP_TOTAL_PLUS_MAG][ CUBE(N_FC_VAR) ],
                                     real g_Slope_PPM[][NCOMP_LR            ][ CUBE(N_SLOPE_PPM) ],
                                     real g_EC_Ele   [][ CUBE(N_EC_ELE) ],
                               const bool Con2Pri, const LR_Limiter_t LR_Limiter, const real MinMod_Coeff,
                               const real dt, const real dh,
                               const real MinDens, const real MinPres, const real MinEint,
                               const bool FracPassive, const int NFrac, const int FracIdx[],
                               const bool JeansMinPres, const real JeansMinPres_Coeff,
                               const EoS_t *EoS );
void Hydro_ComputeFlux( const real g_FC_Var [][NCOMP_TOTAL_PLUS_MAG][ CUBE(N_FC_VAR) ],
                              real g_FC_Flux[][NCOMP_TOTAL_PLUS_MAG][ CUBE(N_FC_FLUX) ],
                        const int NFlux, const int NSkip_N, const int NSkip_T,
                        const bool CorrHalfVel, const real g_Pot_USG[], const double g_Corner[],
                        const real dt, const real dh, const double Time, const bool UsePot,
                        const OptExtAcc_t ExtAcc, const ExtAcc_t ExtAcc_Func, const double ExtAcc_AuxArray[],
                        const real MinDens, const real MinPres, const EoS_t *EoS );
void Hydro_StoreIntFlux( const real g_FC_Flux[][NCOMP_TOTAL_PLUS_MAG][ CUBE(N_FC_FLUX) ],
                               real g_IntFlux[][NCOMP_TOTAL][ SQR(PS2) ],
                         const int NFlux );
void Hydro_FullStepUpdate( const real g_Input[][ CUBE(FLU_NXT) ], real g_Output[][ CUBE(PS2) ], char g_DE_Status[],
                           const real g_FC_B[][ PS2P1*SQR(PS2) ], const real g_Flux[][NCOMP_TOTAL_PLUS_MAG][ CUBE(N_FC_FLUX) ],
                           const real dt, const real dh, const real MinDens, const real MinEint,
                           const real DualEnergySwitch, const bool NormPassive, const int NNorm, const int NormIdx[],
                           const EoS_t *EoS, int *s_FullStepFailure, const int Iteration, const int MinMod_MaxIter );
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
                               const double EoS_AuxArray_Flt[],
                               const int EoS_AuxArray_Int[], const real* const EoS_Table[EOS_NTABLE_MAX] );
#endif
#if ( RSOLVER == HLLC   ||  RSOLVER_RESCUE == HLLC  )
void Hydro_RiemannSolver_HLLC( const int XYZ, real Flux_Out[], const real L_In[], const real R_In[],
                               const real MinDens, const real MinPres, const EoS_DE2P_t EoS_DensEint2Pres,
                               const EoS_DP2C_t EoS_DensPres2CSqr, const EoS_GUESS_t EoS_GuessHTilde,
                               const EoS_H2TEM_t EoS_HTilde2Temp,
                               const double EoS_AuxArray_Flt[],
                               const int EoS_AuxArray_Int[], const real* const EoS_Table[EOS_NTABLE_MAX] );
#endif
#if ( RSOLVER == HLLD   ||  RSOLVER_RESCUE == HLLD  )
void Hydro_RiemannSolver_HLLD( const int XYZ, real Flux_Out[], const real L_In[], const real R_In[],
                               const real MinDens, const real MinPres, const EoS_DE2P_t EoS_DensEint2Pres,
                               const EoS_DP2C_t EoS_DensPres2CSqr, const double EoS_AuxArray_Flt[],
                               const int EoS_AuxArray_Int[], const real* const EoS_Table[EOS_NTABLE_MAX] );
#endif
#if ( FLU_SCHEME == MHM_RP )
void Hydro_Con2Pri( const real In[], real Out[], const real MinPres,
                    const bool FracPassive, const int NFrac, const int FracIdx[],
                    const bool JeansMinPres, const real JeansMinPres_Coeff,
                    const EoS_DE2P_t EoS_DensEint2Pres, const EoS_DP2E_t EoS_DensPres2Eint,
                    const EoS_GUESS_t EoS_GuessHTilde, const EoS_H2TEM_t EoS_HTilde2Temp,
                    const double EoS_AuxArray_Flt[], const int EoS_AuxArray_Int[],
                    const real *const EoS_Table[EOS_NTABLE_MAX], real* const EintOut, real* LorentzFactorPtr );
#endif // #if ( FLU_SCHEME == MHM_RP )
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
#endif // #ifdef MHD

#ifdef COSMIC_RAY
void CR_AdiabaticWork_HalfStep_MHM_RP( real OneCell[NCOMP_TOTAL_PLUS_MAG],
                                       const real g_ConVar_In[][ CUBE(FLU_NXT) ],
                                       const real g_Flux_Half[][NCOMP_TOTAL_PLUS_MAG][ CUBE(N_FC_FLUX) ],
                                       const int idx_fc, const int didx_fc[3],
                                       const int idx_flux, const int didx_flux[3],
                                       const real dt_dh2, const EoS_t *EoS );
void CR_AdiabaticWork_FullStep( const real g_PriVar_Half[][ CUBE(FLU_NXT) ],
                                      real g_Output[][ CUBE(PS2) ],
                                const real g_Flux[][NCOMP_TOTAL_PLUS_MAG][ CUBE(N_FC_FLUX) ],
                                const real g_FC_Var[][NCOMP_TOTAL_PLUS_MAG][ CUBE(N_FC_VAR) ],
                                const real dt, const real dh, const EoS_t *EoS );
#ifdef CR_DIFFUSION
void CR_AddDiffuseFlux_HalfStep( const real g_ConVar[][ CUBE(FLU_NXT) ],
                                       real g_Flux_Half[][NCOMP_TOTAL_PLUS_MAG][ CUBE(N_FC_FLUX) ],
                                 const real g_FC_B[][ SQR(FLU_NXT)*FLU_NXT_P1 ],
                                 const real g_CC_B[][ CUBE(FLU_NXT) ],
                                 const real dh, const MicroPhy_t *MicroPhy );
void CR_AddDiffuseFlux_FullStep( const real g_PriVar_Half[][ CUBE(FLU_NXT) ],
                                       real g_FC_Flux[][NCOMP_TOTAL_PLUS_MAG][ CUBE(N_FC_FLUX) ],
                                 const real g_FC_B_Half[][ FLU_NXT_P1*SQR(FLU_NXT) ],
                                 const int NFlux, const real dh, const MicroPhy_t *MicroPhy );
#endif // #ifdef CR_DIFFUSION
#endif // #ifdef COSMIC_RAY

#endif // #ifdef __CUDACC__ ... else ...


// internal functions
#if ( FLU_SCHEME == MHM_RP )
GPU_DEVICE
static void Hydro_RiemannPredict_Flux( const real g_ConVar[][ CUBE(FLU_NXT) ],
                                             real g_Flux_Half[][NCOMP_TOTAL_PLUS_MAG][ CUBE(N_FC_FLUX) ],
                                       const real g_FC_B[][ SQR(FLU_NXT)*FLU_NXT_P1 ],
                                       const real g_CC_B[][ CUBE(FLU_NXT) ],
                                       const real MinDens, const real MinPres,
                                       const EoS_t *EoS );
GPU_DEVICE
static void Hydro_RiemannPredict( const real g_ConVar_In[][ CUBE(FLU_NXT) ],
                                  const real g_FC_B_Half[][ FLU_NXT_P1*SQR(FLU_NXT) ],
                                  const real g_Flux_Half[][NCOMP_TOTAL_PLUS_MAG][ CUBE(N_FC_FLUX) ],
                                        real g_PriVar_Half[][ CUBE(FLU_NXT) ],
                                  const real dt, const real dh,
                                  const real MinDens, const real MinPres, const real MinEint,
                                  const bool FracPassive, const int NFrac, const int FracIdx[],
                                  const bool JeansMinPres, const real JeansMinPres_Coeff,
                                  const EoS_t *EoS );
#endif




//-------------------------------------------------------------------------------------------------------
// Function    :  CPU/CUFLU_FluidSolver_MHM
// Description :  CPU/GPU fluid solver based on the MUSCL-Hancock scheme
//
// Note        :  1. The three-dimensional evolution is achieved by using the unsplit method
//                2. Two half-step prediction schemes are supported, including "MHM" and "MHM_RP"
//                   MHM    : use interpolated face-centered values to calculate the half-step fluxes
//                   MHM_RP : use Riemann solver to calculate the half-step fluxes
//                3. Ref :
//                   MHM    : "Riemann Solvers and Numerical Methods for Fluid Dynamics
//                             - A Practical Introduction ~ by Eleuterio F. Toro"
//                   MHM_RP : Stone & Gardiner, NewA, 14, 139 (2009)
//                4. See include/CUFLU.h for the values and description of different symbolic constants
//                   such as N_FC_VAR, N_FC_FLUX, N_SLOPE_PPM, N_FL_FLUX, N_HF_VAR
//                5. Arrays with a prefix "g_" are stored in the global memory of GPU
//                6. If an unphysical result occurs in the full-step update, we redo data reconstruction by
//                   reducing the original minmod coefficient repeatedly until either the unphysical result is
//                   solved or the reduced minmod coefficient equals zero. Note that interpolating with a
//                   vanished minmod coefficient is equivalent to the piecewise constant spatial reconstruction.
//
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
//                MinMod_MaxIter     : Maximum number of iterations to reduce MinMod_Coeff
//                Time               : Current physical time                                 (for UNSPLIT_GRAVITY only)
//                UsePot             : Add self-gravity and/or external potential            (for UNSPLIT_GRAVITY only)
//                ExtAcc             : Add external acceleration                             (for UNSPLIT_GRAVITY only)
//                ExtAcc_Func        : Function pointer to the external acceleration routine (for UNSPLIT_GRAVITY only)
//                c_ExtAcc_AuxArray  : Auxiliary array for adding external acceleration      (for UNSPLIT_GRAVITY and CPU only)
//                                     --> When using GPU, this array is stored in the constant memory header
//                                         CUDA_ConstMemory.h and does not need to be passed as a function argument
//                MinDens/Pres/Eint  : Density, pressure, and internal energy floors
//                DualEnergySwitch   : Use the dual-energy formalism if E_int/E_kin < DualEnergySwitch
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
//                MicroPhy           : Microphysics object
//-------------------------------------------------------------------------------------------------------
#ifdef __CUDACC__
__global__
void CUFLU_FluidSolver_MHM(
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
   const LR_Limiter_t LR_Limiter, const real MinMod_Coeff, const int MinMod_MaxIter, const double Time,
   const bool UsePot, const OptExtAcc_t ExtAcc, const ExtAcc_t ExtAcc_Func,
   const real MinDens, const real MinPres, const real MinEint,
   const real DualEnergySwitch,
   const bool NormPassive, const int NNorm,
   const bool FracPassive, const int NFrac,
   const bool JeansMinPres, const real JeansMinPres_Coeff,
   const EoS_t EoS, const MicroPhy_t MicroPhy )
#else
void CPU_FluidSolver_MHM(
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
   const LR_Limiter_t LR_Limiter, const real MinMod_Coeff, const int MinMod_MaxIter, const double Time,
   const bool UsePot, const OptExtAcc_t ExtAcc, const ExtAcc_t ExtAcc_Func,
   const double c_ExtAcc_AuxArray[],
   const real MinDens, const real MinPres, const real MinEint,
   const real DualEnergySwitch,
   const bool NormPassive, const int NNorm, const int c_NormIdx[],
   const bool FracPassive, const int NFrac, const int c_FracIdx[],
   const bool JeansMinPres, const real JeansMinPres_Coeff,
   const EoS_t EoS, const MicroPhy_t MicroPhy )
#endif // #ifdef __CUDACC__ ... else ...
{

#  ifdef UNSPLIT_GRAVITY
   const bool CorrHalfVel          = true;
#  else
   const bool CorrHalfVel          = false;
#  endif
#  if   ( FLU_SCHEME == MHM )
   const bool Con2Pri_Yes          = true;
#  elif ( FLU_SCHEME == MHM_RP )
   const bool Con2Pri_No           = false;
#  endif
#  if ( defined MHD  &&  FLU_SCHEME == MHM_RP )
   const bool CorrHalfVel_No       = false;
   const bool StoreElectric_No     = false;
#  endif
#  if ( defined __CUDACC__  &&  !defined GRAVITY )
   const double *c_ExtAcc_AuxArray = NULL;
#  endif

   int Iteration;
#  ifdef __CUDACC__
   __shared__ int s_FullStepFailure;
#  else
   int s_FullStepFailure;
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
#     pragma omp for schedule( runtime ) private ( Iteration, s_FullStepFailure )
      for (int P=0; P<NPatchGroup; P++)
#     endif
      {
         Iteration = 0;

//       0. point to the arrays associated with different patch groups
//          --> necessary because different patch groups are computed by different OpenMP threads or CUDA blocks in parallel
         real (*const g_FC_Var_1PG   )[NCOMP_TOTAL_PLUS_MAG][ CUBE(N_FC_VAR)    ] = g_FC_Var   [P];
         real (*const g_FC_Flux_1PG  )[NCOMP_TOTAL_PLUS_MAG][ CUBE(N_FC_FLUX)   ] = g_FC_Flux  [P];
         real (*const g_PriVar_1PG   )                      [ CUBE(FLU_NXT)     ] = g_PriVar   [P];
         real (*const g_Slope_PPM_1PG)[NCOMP_LR            ][ CUBE(N_SLOPE_PPM) ] = g_Slope_PPM[P];

#        if ( FLU_SCHEME == MHM_RP )
         real (*const g_Flux_Half_1PG)[NCOMP_TOTAL_PLUS_MAG][ CUBE(N_FC_FLUX) ] = g_FC_Flux_1PG;
         real (*const g_PriVar_Half_1PG )                   [ CUBE(FLU_NXT)   ] = g_PriVar_1PG;

#        ifdef MHD
         real (*const g_FC_Mag_Half_1PG)[ FLU_NXT_P1*SQR(FLU_NXT) ] = g_FC_Mag_Half[P];
         real (*const g_EC_Ele_1PG     )[ CUBE(N_EC_ELE)          ] = g_EC_Ele     [P];
#        else
         real (*const g_FC_Mag_Half_1PG)[ FLU_NXT_P1*SQR(FLU_NXT) ] = NULL;
#        endif
#        endif // if ( FLU_SCHEME == MHM_RP )

#        if ( FLU_SCHEME == MHM )
//       half-step array size is over-estimated here as it only needs CUBE(N_FC_VAR)
//       --> we use the same array size as the half-step variables of MHM_RP to avoid
//           changing the MHM_RP full-step MHD_ComputeElectric()
         real (*const g_PriVar_Half_1PG )[ CUBE(FLU_NXT)  ] = g_PriVar_1PG;
#        ifdef MHD
         real (*const g_EC_Ele_1PG      )[ CUBE(N_EC_ELE) ] = g_EC_Ele[P];
#        else
         real (*const g_EC_Ele_1PG      )[ CUBE(N_EC_ELE) ] = NULL;
#        endif
#        endif // #if ( FLU_SCHEME == MHM )


//       1. half-step prediction
//       1-a. MHM_RP: use Riemann solver to calculate the half-step fluxes
#        if ( FLU_SCHEME == MHM_RP )

#        ifdef MHD
//       1-a-1. evaluate the cell-centered B field and store in g_PriVar[]
//              --> also copy density and compute velocity for MHD_ComputeElectric()
         real CC_B[NCOMP_MAG];

         CGPU_LOOP( idx, CUBE(FLU_NXT) )
         {
            const int size_ij = SQR( FLU_NXT );
            const int i       = idx % FLU_NXT;
            const int j       = idx % size_ij / FLU_NXT;
            const int k       = idx / size_ij;

//          density and velocity
            const real  Dens     = g_Flu_Array_In[P][0][idx];
            const real _Dens     = (real)1.0/Dens;

            g_PriVar_1PG[0][idx] = Dens;
            g_PriVar_1PG[1][idx] = g_Flu_Array_In[P][1][idx]*_Dens;
            g_PriVar_1PG[2][idx] = g_Flu_Array_In[P][2][idx]*_Dens;
            g_PriVar_1PG[3][idx] = g_Flu_Array_In[P][3][idx]*_Dens;

//          magnetic field
            MHD_GetCellCenteredBField( CC_B, g_Mag_Array_In[P][0], g_Mag_Array_In[P][1], g_Mag_Array_In[P][2],
                                       FLU_NXT, FLU_NXT, FLU_NXT, i, j, k );

            for (int v=0; v<NCOMP_MAG; v++)  g_PriVar_1PG[ MAG_OFFSET + v ][idx] = CC_B[v];
         }

#        ifdef __CUDACC__
         __syncthreads();
#        endif
#        endif // #ifdef MHD


//       1-a-2. evaluate the half-step first-order fluxes by Riemann solver
//       hydrodynamic fluxes
         Hydro_RiemannPredict_Flux( g_Flu_Array_In[P], g_Flux_Half_1PG, g_Mag_Array_In[P], g_PriVar_1PG+MAG_OFFSET,
                                    MinDens, MinPres, &EoS );

//       add cosmic-ray fluxes
#        ifdef CR_DIFFUSION
         CR_AddDiffuseFlux_HalfStep( g_Flu_Array_In[P], g_Flux_Half_1PG, g_Mag_Array_In[P], g_PriVar_1PG+MAG_OFFSET, dh, &MicroPhy );
#        endif


//       1-a-3. evaluate electric field and update B field at the half time-step
#        ifdef MHD
         MHD_ComputeElectric( g_EC_Ele_1PG, g_Flux_Half_1PG, g_PriVar_1PG, N_HF_ELE, N_HF_FLUX,
                              FLU_NXT, 0, dt, dh, StoreElectric_No, NULL,
                              CorrHalfVel_No, NULL, NULL, NULL_REAL,
                              EXT_POT_NONE, EXT_ACC_NONE, NULL, NULL );

         MHD_UpdateMagnetic( g_FC_Mag_Half_1PG[0], g_FC_Mag_Half_1PG[1], g_FC_Mag_Half_1PG[2],
                             g_Mag_Array_In[P], g_EC_Ele_1PG, (real)0.5*dt, dh, N_HF_VAR, N_HF_ELE, 1 );
#        endif


//       1-a-4. evaluate the half-step solutions
         Hydro_RiemannPredict( g_Flu_Array_In[P], g_FC_Mag_Half_1PG, g_Flux_Half_1PG, g_PriVar_Half_1PG,
                               dt, dh, MinDens, MinPres, MinEint, FracPassive, NFrac, c_FracIdx,
                               JeansMinPres, JeansMinPres_Coeff, &EoS );


         do {

#           ifdef __CUDACC__
            __syncthreads();
#           endif
            s_FullStepFailure = 0;
#           ifdef __CUDACC__
            __syncthreads();
#           endif

            real AdaptiveMinModCoeff = ( MinMod_MaxIter == 0 ) ? MinMod_Coeff :
            MinMod_Coeff - (real)Iteration * MinMod_Coeff / (real)MinMod_MaxIter;


//          ensure adaptive MinMod_Coeff is non-negative
            AdaptiveMinModCoeff = FMAX( AdaptiveMinModCoeff, (real)0.0 );


//          1-a-5. evaluate the face-centered values by data reconstruction
//                 --> note that g_PriVar_Half_1PG[] returned by Hydro_RiemannPredict() stores the primitive variables
            Hydro_DataReconstruction( NULL, g_FC_Mag_Half_1PG, g_PriVar_Half_1PG, g_FC_Var_1PG, g_Slope_PPM_1PG,
                                      NULL, Con2Pri_No, LR_Limiter, AdaptiveMinModCoeff, dt, dh,
                                      MinDens, MinPres, MinEint, FracPassive, NFrac, c_FracIdx,
                                      JeansMinPres, JeansMinPres_Coeff, &EoS );


//       1-b. MHM: use interpolated face-centered values to calculate the half-step fluxes
#        elif ( FLU_SCHEME == MHM )

         do {

#           ifdef __CUDACC__
            __syncthreads();
#           endif
            s_FullStepFailure = 0;
#           ifdef __CUDACC__
            __syncthreads();
#           endif

            real AdaptiveMinModCoeff = ( MinMod_MaxIter == 0 ) ? MinMod_Coeff :
            MinMod_Coeff - (real)Iteration * MinMod_Coeff / (real)MinMod_MaxIter;


//          ensure adaptive MinMod_Coeff is non-negative
            AdaptiveMinModCoeff = FMAX( AdaptiveMinModCoeff, (real)0.0 );


//          evaluate the face-centered values by data reconstruction
            Hydro_DataReconstruction( g_Flu_Array_In[P], g_Mag_Array_In[P], g_PriVar_Half_1PG, g_FC_Var_1PG, g_Slope_PPM_1PG,
                                      g_EC_Ele_1PG, Con2Pri_Yes, LR_Limiter, AdaptiveMinModCoeff, dt, dh,
                                      MinDens, MinPres, MinEint, FracPassive, NFrac, c_FracIdx,
                                      JeansMinPres, JeansMinPres_Coeff, &EoS );

#        endif // #if ( FLU_SCHEME == MHM_RP ) ... else ...


//          2. evaluate the full-step fluxes
#           ifdef MHD
            const int NSkip_N = 0;
            const int NSkip_T = 0;
#           else
            const int NSkip_N = 0;
            const int NSkip_T = 1;
#           endif

//          hydrodynamic fluxes
            Hydro_ComputeFlux( g_FC_Var_1PG, g_FC_Flux_1PG, N_FL_FLUX, NSkip_N, NSkip_T,
                               CorrHalfVel, g_Pot_Array_USG[P], g_Corner_Array[P],
                               dt, dh, Time, UsePot, ExtAcc, ExtAcc_Func, c_ExtAcc_AuxArray,
                               MinDens, MinPres, &EoS );

//          add cosmic-ray fluxes
#           ifdef CR_DIFFUSION
            CR_AddDiffuseFlux_FullStep( g_PriVar_Half_1PG, g_FC_Flux_1PG, g_FC_Mag_Half_1PG, N_FL_FLUX, dh, &MicroPhy );
#           endif


            if ( StoreFlux )
               Hydro_StoreIntFlux( g_FC_Flux_1PG, g_Flux_Array[P], N_FL_FLUX );


//          3. evaluate electric field and update B field at the full time-step
//             --> must update B field before Hydro_FullStepUpdate() since the latter requires
//                 the updated magnetic energy when adopting the dual-energy formalism
#           ifdef MHD
#           if ( FLU_SCHEME == MHM )
            const int OffsetPri = 0;
#           else
            const int OffsetPri = LR_GHOST_SIZE;
#           endif

            MHD_ComputeElectric( g_EC_Ele_1PG, g_FC_Flux_1PG, g_PriVar_Half_1PG, N_FL_ELE, N_FL_FLUX,
                                 N_HF_VAR, OffsetPri, dt, dh, StoreElectric, g_Ele_Array[P],
                                 CorrHalfVel, g_Pot_Array_USG[P], g_Corner_Array[P], Time,
                                 UsePot, ExtAcc, ExtAcc_Func, c_ExtAcc_AuxArray );

            MHD_UpdateMagnetic( g_Mag_Array_Out[P][0], g_Mag_Array_Out[P][1], g_Mag_Array_Out[P][2],
                                g_Mag_Array_In[P], g_EC_Ele_1PG, dt, dh, PS2, N_FL_ELE, FLU_GHOST_SIZE );
#           endif // #ifdef MHD


//          4. full-step evolution
            Hydro_FullStepUpdate( g_Flu_Array_In[P], g_Flu_Array_Out[P], g_DE_Array_Out[P], g_Mag_Array_Out[P],
                                  g_FC_Flux_1PG, dt, dh, MinDens, MinEint, DualEnergySwitch,
                                  NormPassive, NNorm, c_NormIdx, &EoS, &s_FullStepFailure, Iteration, MinMod_MaxIter );

//          add the cosmic-ray source term of adiabatic work
#           ifdef COSMIC_RAY
            CR_AdiabaticWork_FullStep( g_PriVar_Half_1PG, g_Flu_Array_Out[P], g_FC_Flux_1PG, g_FC_Var_1PG,
                                       dt, dh, &EoS );
#           endif


//          5. counter increment
            Iteration++;

         } while ( s_FullStepFailure  &&  Iteration <= MinMod_MaxIter );

      } // loop over all patch groups
   } // OpenMP parallel region

} // FUNCTION : CPU_FluidSolver_MHM



#if ( FLU_SCHEME == MHM_RP )
//-------------------------------------------------------------------------------------------------------
// Function    :  Hydro_RiemannPredict_Flux
// Description :  Evaluate the half-step face-centered fluxes by Riemann solver
//
// Note        :  1. Work for the MHM_RP scheme
//                2. Currently support the exact, HLLC, HLLE, HLLD, and Roe solvers
//                3. g_Flux_Half[] is accessed with a stride N_HF_FLUX
//                   --> Fluxes on the **left** face of the (i+1,j+1,k+1) element in g_ConVar[] will
//                       be stored in the (i,j,k) element of g_Flux_Half[]
//
// Parameter   :  g_ConVar     : Array storing the input conserved variables
//                g_Flux_Half  : Array to store the output face-centered fluxes
//                g_FC_B       : Array storing the input face-centered magnetic field (for MHD only)
//                               --> Accessed with strides FLU_NXT/FLU_NXT+1 along the
//                                   transverse/longitudinal directions
//                g_CC_B       : Array storing the input cell-centered magnetic field (for MHD only)
//                               --> Accessed with a stride FLU_NXT
//                MinDens/Pres : Density and pressure floors
//                EoS          : EoS object
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE
void Hydro_RiemannPredict_Flux( const real g_ConVar[][ CUBE(FLU_NXT) ],
                                      real g_Flux_Half[][NCOMP_TOTAL_PLUS_MAG][ CUBE(N_FC_FLUX) ],
                                const real g_FC_B[][ SQR(FLU_NXT)*FLU_NXT_P1 ],
                                const real g_CC_B[][ CUBE(FLU_NXT) ],
                                const real MinDens, const real MinPres,
                                const EoS_t *EoS )
{

   const int didx_cvar[3] = { 1, FLU_NXT, SQR(FLU_NXT) };
   real ConVar_L[NCOMP_TOTAL_PLUS_MAG], ConVar_R[NCOMP_TOTAL_PLUS_MAG], Flux_1Face[NCOMP_TOTAL_PLUS_MAG];

// loop over different spatial directions
   for (int d=0; d<3; d++)
   {
#     ifdef MHD
      const int TDir1          = (d+1)%3;    // transverse direction 1
      const int TDir2          = (d+2)%3;    // transverse direction 2
      const int stride_fc_B[3] = { 1, FLU_NXT, SQR(FLU_NXT) };

      int sizeB_i, sizeB_j;
#     endif

      int i_cvar_s=0, j_cvar_s=0, k_cvar_s=0, size_i, size_j, size_k;

      switch ( d )
      {
#        ifdef MHD
         case 0 : size_i  = N_HF_FLUX-1;  size_j  = N_HF_FLUX-0;  size_k = N_HF_FLUX-0;
                  sizeB_i = FLU_NXT_P1;   sizeB_j = FLU_NXT;
                  break;

         case 1 : size_i  = N_HF_FLUX-0;  size_j  = N_HF_FLUX-1;  size_k = N_HF_FLUX-0;
                  sizeB_i = FLU_NXT;      sizeB_j = FLU_NXT_P1;
                  break;

         case 2 : size_i  = N_HF_FLUX-0;  size_j  = N_HF_FLUX-0;  size_k = N_HF_FLUX-1;
                  sizeB_i = FLU_NXT;      sizeB_j = FLU_NXT;
                  break;

#        else // #ifdef MHD
         case 0 : i_cvar_s = 0;            j_cvar_s = 1;            k_cvar_s = 1;
                  size_i   = N_HF_FLUX-0;  size_j   = N_HF_FLUX-1;  size_k   = N_HF_FLUX-1;
                  break;

         case 1 : i_cvar_s = 1;            j_cvar_s = 0;            k_cvar_s = 1;
                  size_i   = N_HF_FLUX-1;  size_j   = N_HF_FLUX-0;  size_k   = N_HF_FLUX-1;
                  break;

         case 2 : i_cvar_s = 1;            j_cvar_s = 1;            k_cvar_s = 0;
                  size_i   = N_HF_FLUX-1;  size_j   = N_HF_FLUX-1;  size_k   = N_HF_FLUX-0;
                  break;
#        endif // #ifdef MHD ... else ...
      } // switch ( d )

      const int size_ij = size_i*size_j;

      CGPU_LOOP( idx, size_i*size_j*size_k )
      {
         const int i_flux   = idx % size_i;
         const int j_flux   = idx % size_ij / size_i;
         const int k_flux   = idx / size_ij;
         const int idx_flux = IDX321( i_flux, j_flux, k_flux, N_HF_FLUX, N_HF_FLUX );

         const int i_cvar   = i_flux + i_cvar_s;
         const int j_cvar   = j_flux + j_cvar_s;
         const int k_cvar   = k_flux + k_cvar_s;
         const int idx_cvar = IDX321( i_cvar, j_cvar, k_cvar, FLU_NXT, FLU_NXT );

//       get the left and right fluid variables
         for (int v=0; v<NCOMP_TOTAL; v++)
         {
            ConVar_L[v] = g_ConVar[v][ idx_cvar                ];
            ConVar_R[v] = g_ConVar[v][ idx_cvar + didx_cvar[d] ];
         }

//       get the left and right B field
#        ifdef MHD
//       longitudinal component is face-centered
         const int idx_fc_B = IDX321( i_cvar, j_cvar, k_cvar, sizeB_i, sizeB_j ) + stride_fc_B[d];
         ConVar_L[ MAG_OFFSET + d     ] = g_FC_B[d][idx_fc_B];
         ConVar_R[ MAG_OFFSET + d     ] = g_FC_B[d][idx_fc_B];

//       transverse components are cell-centered
         ConVar_L[ MAG_OFFSET + TDir1 ] = g_CC_B[TDir1][ idx_cvar                ];
         ConVar_L[ MAG_OFFSET + TDir2 ] = g_CC_B[TDir2][ idx_cvar                ];
         ConVar_R[ MAG_OFFSET + TDir1 ] = g_CC_B[TDir1][ idx_cvar + didx_cvar[d] ];
         ConVar_R[ MAG_OFFSET + TDir2 ] = g_CC_B[TDir2][ idx_cvar + didx_cvar[d] ];

//       correct total energy by the difference between the face- and cell-centered longitudinal B field
         ConVar_L[4] += (real)0.5*(  SQR( ConVar_L[MAG_OFFSET + d] ) - SQR( g_CC_B[d][idx_cvar               ] )  );
         ConVar_R[4] += (real)0.5*(  SQR( ConVar_R[MAG_OFFSET + d] ) - SQR( g_CC_B[d][idx_cvar + didx_cvar[d]] )  );
#        endif

//       invoke the Riemann solver
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

//       switch to a different Riemann solver if the default one fails
#        if ( RSOLVER_RESCUE != NONE )
         for (int v=0; v<NCOMP_TOTAL_PLUS_MAG; v++)
         {
//          only check NaN for now
            if ( Flux_1Face[v] != Flux_1Face[v] )
            {
#              ifdef CHECK_UNPHYSICAL_IN_FLUID
               printf( "WARNING : default Riemann solver failed in Hydro_RiemannPredict_Flux() --> switch to RSOLVER_RESCUE (%d) !!\n", RSOLVER_RESCUE );
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

//       store the results in g_Flux_Half[]
         for (int v=0; v<NCOMP_TOTAL_PLUS_MAG; v++)   g_Flux_Half[d][v][idx_flux] = Flux_1Face[v];
      } // CGPU_LOOP( idx, N_HF_FLUX*SQR(N_HF_FLUX-1) )
   } // for (int d=0; d<3; d++)


#  ifdef __CUDACC__
   __syncthreads();
#  endif

} // FUNCTION : Hydro_RiemannPredict_Flux



//-------------------------------------------------------------------------------------------------------
// Function    :  Hydro_RiemannPredict
// Description :  Evolve the cell-centered variables by half time-step using the fluxes returned
//                by Hydro_RiemannPredict_Flux()
//
// Note        :  1. Work for the MHM_RP scheme
//                2. For the performance consideration, the output data are converted to primitive variables
//                   --> Reducing the global memory access on GPU
//                3. Cell-centered B field is simply obtained by averaging the half-step face-centered B field
//
// Parameter   :  g_ConVar_In        : Array storing the input conserved variables
//                g_FC_B_Half        : Array storing the input half-step face-centered B field
//                g_Flux_Half        : Array storing the input face-centered fluxes
//                                     --> Accessed with the stride N_HF_FLUX
//                g_PriVar_Half      : Array to store the output primitive variables
//                                     --> Accessed with the stride N_HF_VAR
//                                     --> Although its actually allocated size is FLU_NXT^3 since it points to g_PriVar_1PG[]
//                dt                 : Time interval to advance solution
//                dh                 : Cell size
//                MinDens/Pres/Eint  : Density, pressure, and internal energy floors
//                FracPassive        : true --> convert passive scalars to mass fraction during data reconstruction
//                NFrac              : Number of passive scalars for the option "FracPassive"
//                FracIdx            : Target variable indices for the option "FracPassive"
//                JeansMinPres       : Apply minimum pressure estimated from the Jeans length
//                JeansMinPres_Coeff : Coefficient used by JeansMinPres = G*(Jeans_NCell*Jeans_dh)^2/(Gamma*pi);
//                EoS                : EoS object
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE
void Hydro_RiemannPredict( const real g_ConVar_In[][ CUBE(FLU_NXT) ],
                           const real g_FC_B_Half[][ FLU_NXT_P1*SQR(FLU_NXT) ],
                           const real g_Flux_Half[][NCOMP_TOTAL_PLUS_MAG][ CUBE(N_FC_FLUX) ],
                                 real g_PriVar_Half[][ CUBE(FLU_NXT) ],
                           const real dt, const real dh,
                           const real MinDens, const real MinPres, const real MinEint,
                           const bool FracPassive, const int NFrac, const int FracIdx[],
                           const bool JeansMinPres, const real JeansMinPres_Coeff,
                           const EoS_t *EoS )
{

   const int  didx_flux[3] = { 1, N_HF_FLUX, SQR(N_HF_FLUX) };
#  ifdef COSMIC_RAY
   const int  didx_in[3]   = { 1, FLU_NXT, SQR(FLU_NXT) };
#  endif
   const real dt_dh2       = (real)0.5*dt/dh;

   const int N_HF_VAR2 = SQR(N_HF_VAR);
   CGPU_LOOP( idx_out, CUBE(N_HF_VAR) )
   {
      const int i_out    = idx_out % N_HF_VAR;
      const int j_out    = idx_out % N_HF_VAR2 / N_HF_VAR;
      const int k_out    = idx_out / N_HF_VAR2;

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
      const int idx_flux = IDX321( i_flux, j_flux, k_flux, N_HF_FLUX, N_HF_FLUX );

      const int i_in     = i_out + 1;
      const int j_in     = j_out + 1;
      const int k_in     = k_out + 1;
      const int idx_in   = IDX321( i_in, j_in, k_in, FLU_NXT, FLU_NXT );

      real out_con[NCOMP_TOTAL_PLUS_MAG], out_pri[NCOMP_TOTAL_PLUS_MAG], dflux[3][NCOMP_TOTAL];
#     ifdef LR_EINT
      real Eint;
      real* const EintPtr = &Eint;
#     else
      real* const EintPtr = NULL;
#     endif

//    calculate the flux differences of the fluid variables
      for (int d=0; d<3; d++)
      for (int v=0; v<NCOMP_TOTAL; v++)
      {
#        ifdef MHD
         dflux[d][v] = g_Flux_Half[d][v][idx_flux] - g_Flux_Half[d][v][ idx_flux - didx_flux[d] ];
#        else
         dflux[d][v] = g_Flux_Half[d][v][ idx_flux + didx_flux[d] ] - g_Flux_Half[d][v][idx_flux];
#        endif
      }

//    update the input cell-centered conserved variables with the flux differences
      for (int v=0; v<NCOMP_TOTAL; v++)
         out_con[v] = g_ConVar_In[v][idx_in] - dt_dh2*( dflux[0][v] + dflux[1][v] + dflux[2][v] );


//    add the cosmic-ray source term of adiabatic work
#     ifdef COSMIC_RAY
      CR_AdiabaticWork_HalfStep_MHM_RP( out_con, g_ConVar_In, g_Flux_Half, idx_in, didx_in,
                                        idx_flux, didx_flux, dt_dh2, EoS );
#     endif


//    compute the cell-centered half-step B field
#     ifdef MHD
      MHD_GetCellCenteredBField( out_con+MAG_OFFSET, g_FC_B_Half[0], g_FC_B_Half[1], g_FC_B_Half[2],
                                 N_HF_VAR, N_HF_VAR, N_HF_VAR, i_out, j_out, k_out );
#     endif

//    apply density and internal energy floors
      out_con[0] = FMAX( out_con[0], MinDens );
#     ifndef SRHD
#     ifndef BAROTROPIC_EOS
#     ifdef MHD
      const real Emag = (real)0.5*( SQR(out_con[MAG_OFFSET+0]) + SQR(out_con[MAG_OFFSET+1]) + SQR(out_con[MAG_OFFSET+2]) );
#     else
      const real Emag = NULL_REAL;
#     endif
      out_con[4] = Hydro_CheckMinEintInEngy( out_con[0], out_con[1], out_con[2], out_con[3], out_con[4], MinEint, Emag );
#     endif // #ifndef BAROTROPIC_EOS
#     endif // #ifndef SRHD
#     if ( NCOMP_PASSIVE > 0 )
      for (int v=NCOMP_FLUID; v<NCOMP_TOTAL; v++)
      out_con[v] = FMAX( out_con[v], TINY_NUMBER );
#     endif

//    conserved --> primitive variables
      Hydro_Con2Pri( out_con, out_pri, MinPres, FracPassive, NFrac, FracIdx, JeansMinPres, JeansMinPres_Coeff,
                     EoS->DensEint2Pres_FuncPtr, EoS->DensPres2Eint_FuncPtr,
                     EoS->GuessHTilde_FuncPtr, EoS->HTilde2Temp_FuncPtr,
                     EoS->AuxArrayDevPtr_Flt, EoS->AuxArrayDevPtr_Int,
                     EoS->Table, EintPtr, NULL );

//    store the results in g_PriVar_Half[]
      for (int v=0; v<NCOMP_TOTAL_PLUS_MAG; v++)   g_PriVar_Half[v][idx_out] = out_pri[v];

//    store Eint in the last variable for LR_EINT
#     ifdef LR_EINT
      g_PriVar_Half[NCOMP_TOTAL_PLUS_MAG][idx_out] = Hydro_CheckMinEint( Eint, MinEint );
#     endif
   } // i,j,k


#  ifdef __CUDACC__
   __syncthreads();
#  endif

} // FUNCTION : Hydro_RiemannPredict



#endif // #if ( FLU_SCHEME == MHM_RP )



#endif // #if (  MODEL == HYDRO  &&  ( FLU_SCHEME == MHM || FLU_SCHEME == MHM_RP )  )
