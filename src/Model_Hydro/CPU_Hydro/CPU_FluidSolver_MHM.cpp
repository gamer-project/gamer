#include "CUFLU.h"

#if (  MODEL == HYDRO  &&  ( FLU_SCHEME == MHM || FLU_SCHEME == MHM_RP )  )



// external functions
#ifdef __CUDACC__

#include "CUFLU_Shared_FluUtility.cu"
#include "CUFLU_Shared_DataReconstruction.cu"
#include "CUFLU_Shared_ComputeFlux.cu"
#include "CUFLU_Shared_FullStepUpdate.cu"
#include "CUFLU_SetConstMem_FluidSolver.cu"
#ifdef MHD
#include "CUFLU_Shared_ConstrainedTransport.cu"
#endif

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

#else // #ifdef __CUDACC__

void Hydro_DataReconstruction( const real g_ConVar   [][ CUBE(FLU_NXT) ],
                               const real g_FC_B     [][ SQR(FLU_NXT)*FLU_NXT_P1 ],
                                     real g_PriVar   [][ CUBE(FLU_NXT) ],
                                     real g_FC_Var   [][NCOMP_TOTAL_PLUS_MAG][ CUBE(N_FC_VAR) ],
                                     real g_Slope_PPM[][NCOMP_TOTAL_PLUS_MAG][ CUBE(N_SLOPE_PPM) ],
                               const bool Con2Pri, const int NIn, const int NGhost, const real Gamma,
                               const LR_Limiter_t LR_Limiter, const real MinMod_Coeff,
                               const real dt, const real dh, const real MinDens, const real MinPres,
                               const bool NormPassive, const int NNorm, const int NormIdx[],
                               const bool JeansMinPres, const real JeansMinPres_Coeff );
void Hydro_ComputeFlux( const real g_FC_Var [][NCOMP_TOTAL_PLUS_MAG][ CUBE(N_FC_VAR) ],
                              real g_FC_Flux[][NCOMP_TOTAL_PLUS_MAG][ CUBE(N_FC_FLUX) ],
                        const int NFlux, const int NSkip_N, const int NSkip_T, const real Gamma,
                        const bool CorrHalfVel, const real g_Pot_USG[], const double g_Corner[],
                        const real dt, const real dh, const double Time,
                        const OptGravityType_t GravityType, const double ExtAcc_AuxArray[],
                        const real MinPres, const bool DumpIntFlux, real g_IntFlux[][NCOMP_TOTAL][ SQR(PS2) ] );
void Hydro_FullStepUpdate( const real g_Input[][ CUBE(FLU_NXT) ], real g_Output[][ CUBE(PS2) ], char g_DE_Status[],
                           const real g_FC_B[][ PS2P1*SQR(PS2) ], const real g_Flux[][NCOMP_TOTAL_PLUS_MAG][ CUBE(N_FC_FLUX) ],
                           const real dt, const real dh, const real Gamma, const real MinDens, const real MinPres,
                           const real DualEnergySwitch, const bool NormPassive, const int NNorm, const int NormIdx[] );
#if   ( RSOLVER == EXACT )
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
#if ( FLU_SCHEME == MHM_RP )
void Hydro_Con2Pri( const real In[], real Out[], const real Gamma_m1, const real MinPres,
                    const bool NormPassive, const int NNorm, const int NormIdx[],
                    const bool JeansMinPres, const real JeansMinPres_Coeff );
real Hydro_CheckMinPresInEngy( const real Dens, const real MomX, const real MomY, const real MomZ, const real Engy,
                               const real Gamma_m1, const real _Gamma_m1, const real MinPres );
#ifdef MHD
void MHD_ComputeElectric(       real g_EC_Ele[][ CUBE(N_EC_ELE) ],
                          const real g_FC_Flux[][NCOMP_TOTAL_PLUS_MAG][ CUBE(N_FC_FLUX) ],
                          const real g_PriVar[][ CUBE(FLU_NXT) ],
                          const int NEle, const int NFlux, const int NPri, const int OffsetPri,
                          const real dt, const real dh,
                          const bool DumpIntEle, real g_IntEle[][NCOMP_ELE][ PS2P1*PS2 ],
                          const bool CorrHalfVel, const real g_Pot_USG[], const double g_Corner[],
                          const double Time, const OptGravityType_t GravityType, const double ExtAcc_AuxArray[] );
void MHD_UpdateMagnetic( real *g_FC_Bx_Out, real *g_FC_By_Out, real *g_FC_Bz_Out,
                         const real g_FC_B_In[][ FLU_NXT_P1*SQR(FLU_NXT) ],
                         const real g_EC_Ele[][ CUBE(N_EC_ELE) ],
                         const real dt, const real dh, const int NOut, const int NEle, const int Offset_B_In );
#endif // #ifdef MHD
#endif // #if ( FLU_SCHEME == MHM_RP )

#endif // #ifdef __CUDACC__ ... else ...


// internal functions
#if ( FLU_SCHEME == MHM_RP )
GPU_DEVICE
static void Hydro_RiemannPredict_Flux( const real g_ConVar[][ CUBE(FLU_NXT) ],
                                             real g_Flux_Half[][NCOMP_TOTAL_PLUS_MAG][ CUBE(N_FC_FLUX) ],
                                       const real g_FC_B[][ SQR(FLU_NXT)*FLU_NXT_P1 ],
                                       const real g_CC_B[][ CUBE(FLU_NXT) ],
                                       const real Gamma, const real MinPres );
GPU_DEVICE
static void Hydro_RiemannPredict( const real g_ConVar_In[][ CUBE(FLU_NXT) ],
                                  const real g_FC_B_Half[][ FLU_NXT_P1*SQR(FLU_NXT) ],
                                  const real g_Flux_Half[][NCOMP_TOTAL_PLUS_MAG][ CUBE(N_FC_FLUX) ],
                                        real g_PriVar_Half[][ CUBE(FLU_NXT) ],
                                  const real dt, const real dh, const real Gamma, const real MinDens, const real MinPres,
                                  const bool NormPassive, const int NNorm, const int NormIdx[],
                                  const bool JeansMinPres, const real JeansMinPres_Coeff );
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
//                Gamma              : Ratio of specific heats
//                StoreFlux          : true --> store the coarse-fine fluxes
//                StoreElectric      : true --> store the coarse-fine electric field
//                LR_Limiter         : Slope limiter for the data reconstruction in the MHM/MHM_RP/CTU schemes
//                                     (0/1/2/3/4) = (vanLeer/generalized MinMod/vanAlbada/
//                                                    vanLeer + generalized MinMod/extrema-preserving) limiter
//                MinMod_Coeff       : Coefficient of the generalized MinMod limiter
//                Time               : Current physical time                                     (for UNSPLIT_GRAVITY only)
//                GravityType        : Types of gravity --> self-gravity, external gravity, both (for UNSPLIT_GRAVITY only)
//                c_ExtAcc_AuxArray  : Auxiliary array for adding external acceleration          (for UNSPLIT_GRAVITY only)
//                                     --> When using GPU, this array is stored in the constant memory and does
//                                         not need to be passed as a function argument
//                                         --> Declared in CUFLU_SetConstMem_FluidSolver.cu with the prefix "c_" to
//                                             highlight that this is a constant variable on GPU
//                MinDens/Pres       : Minimum allowed density and pressure
//                DualEnergySwitch   : Use the dual-energy formalism if E_int/E_kin < DualEnergySwitch
//                NormPassive        : true --> normalize passive scalars so that the sum of their mass density
//                                              is equal to the gas mass density
//                NNorm              : Number of passive scalars to be normalized
//                                     --> Should be set to the global variable "PassiveNorm_NVar"
//                c_NormIdx          : Target variable indices to be normalized
//                                     --> Should be set to the global variable "PassiveNorm_VarIdx"
//                                     --> When using GPU, this array is stored in the constant memory and does
//                                         not need to be passed as a function argument
//                                         --> Declared in CUFLU_SetConstMem_FluidSolver.cu with the prefix "c_" to
//                                             highlight that this is a constant variable on GPU
//                JeansMinPres       : Apply minimum pressure estimated from the Jeans length
//                JeansMinPres_Coeff : Coefficient used by JeansMinPres = G*(Jeans_NCell*Jeans_dh)^2/(Gamma*pi);
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
         real   g_PriVar       []   [NCOMP_TOTAL_PLUS_MAG][ CUBE(FLU_NXT) ],
         real   g_Slope_PPM    [][3][NCOMP_TOTAL_PLUS_MAG][ CUBE(N_SLOPE_PPM) ],
         real   g_FC_Var       [][6][NCOMP_TOTAL_PLUS_MAG][ CUBE(N_FC_VAR) ],
         real   g_FC_Flux      [][3][NCOMP_TOTAL_PLUS_MAG][ CUBE(N_FC_FLUX) ],
         real   g_FC_Mag_Half  [][NCOMP_MAG][ FLU_NXT_P1*SQR(FLU_NXT) ],
         real   g_EC_Ele       [][NCOMP_MAG][ CUBE(N_EC_ELE) ],
   const real dt, const real dh, const real Gamma,
   const bool StoreFlux, const bool StoreElectric,
   const LR_Limiter_t LR_Limiter, const real MinMod_Coeff,
   const double Time, const OptGravityType_t GravityType,
   const real MinDens, const real MinPres, const real DualEnergySwitch,
   const bool NormPassive, const int NNorm,
   const bool JeansMinPres, const real JeansMinPres_Coeff )
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
         real   g_PriVar       []   [NCOMP_TOTAL_PLUS_MAG][ CUBE(FLU_NXT) ],
         real   g_Slope_PPM    [][3][NCOMP_TOTAL_PLUS_MAG][ CUBE(N_SLOPE_PPM) ],
         real   g_FC_Var       [][6][NCOMP_TOTAL_PLUS_MAG][ CUBE(N_FC_VAR) ],
         real   g_FC_Flux      [][3][NCOMP_TOTAL_PLUS_MAG][ CUBE(N_FC_FLUX) ],
         real   g_FC_Mag_Half  [][NCOMP_MAG][ FLU_NXT_P1*SQR(FLU_NXT) ],
         real   g_EC_Ele       [][NCOMP_MAG][ CUBE(N_EC_ELE) ],
   const int NPatchGroup, const real dt, const real dh, const real Gamma,
   const bool StoreFlux, const bool StoreElectric,
   const LR_Limiter_t LR_Limiter, const real MinMod_Coeff,
   const double Time, const OptGravityType_t GravityType,
   const double c_ExtAcc_AuxArray[], const real MinDens, const real MinPres,
   const real DualEnergySwitch, const bool NormPassive, const int NNorm, const int c_NormIdx[],
   const bool JeansMinPres, const real JeansMinPres_Coeff )
#endif // #ifdef __CUDACC__ ... else ...
{

#  ifdef UNSPLIT_GRAVITY
   const bool CorrHalfVel          = true;
#  else
   const bool CorrHalfVel          = false;
#  endif
   const bool CorrHalfVel_No       = false;
#  if   ( FLU_SCHEME == MHM )
   const bool Con2Pri_Yes          = true;
#  elif ( FLU_SCHEME == MHM_RP )
   const bool Con2Pri_No           = false;
#  endif
#  ifdef MHD
   const bool StoreElectric_No     = false;
#  endif
#  if ( defined __CUDACC__  &&  !defined UNSPLIT_GRAVITY )
   const double *c_ExtAcc_AuxArray = NULL;
#  endif


// openmp pragma for the CPU solver
#  ifndef __CUDACC__
#  pragma omp parallel
#  endif
   {
//    point to the arrays associated with different OpenMP threads (for CPU) or CUDA thread blocks (for GPU)
#     ifdef __CUDACC__
      const int array_idx = blockIdx.x;
#     else
#     ifdef OPENMP
      const int array_idx = omp_get_thread_num();
#     else
      const int array_idx = 0;
#     endif
#     endif // #ifdef __CUDACC__ ... else ...

      real (*const g_FC_Var_1PG   )[NCOMP_TOTAL_PLUS_MAG][ CUBE(N_FC_VAR)    ] = g_FC_Var   [array_idx];
      real (*const g_FC_Flux_1PG  )[NCOMP_TOTAL_PLUS_MAG][ CUBE(N_FC_FLUX)   ] = g_FC_Flux  [array_idx];
      real (*const g_PriVar_1PG   )                      [ CUBE(FLU_NXT)     ] = g_PriVar   [array_idx];
      real (*const g_Slope_PPM_1PG)[NCOMP_TOTAL_PLUS_MAG][ CUBE(N_SLOPE_PPM) ] = g_Slope_PPM[array_idx];

#     ifdef MHD
      real (*const g_FC_Mag_Half_1PG)[ FLU_NXT_P1*SQR(FLU_NXT) ] = g_FC_Mag_Half[array_idx];
      real (*const g_EC_Ele_1PG     )[ CUBE(N_EC_ELE)          ] = g_EC_Ele     [array_idx];
#     else
      real (*const g_FC_Mag_Half_1PG)[ FLU_NXT_P1*SQR(FLU_NXT) ] = NULL;
      real (*const g_EC_Ele_1PG     )[ CUBE(N_EC_ELE)          ] = NULL;
#     endif

#     if ( FLU_SCHEME == MHM_RP )
      real (*const g_Flux_Half_1PG)[NCOMP_TOTAL_PLUS_MAG][ CUBE(N_FC_FLUX) ] = g_FC_Flux_1PG;
      real (*const g_PriVar_Half_1PG )                   [ CUBE(FLU_NXT)   ] = g_PriVar_1PG;
#     endif


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
         Hydro_RiemannPredict_Flux( g_Flu_Array_In[P], g_Flux_Half_1PG, g_Mag_Array_In[P], g_PriVar_1PG+MAG_OFFSET,
                                    Gamma, MinPres );


//       1-a-3. evaluate electric field and update B field at the half time-step
#        ifdef MHD
         MHD_ComputeElectric( g_EC_Ele_1PG, g_Flux_Half_1PG, g_PriVar_1PG, N_HF_ELE, N_HF_FLUX,
                              FLU_NXT, 0, dt, dh, StoreElectric_No, NULL,
                              CorrHalfVel_No, NULL, NULL, NULL_REAL, GRAVITY_NONE, NULL );

         MHD_UpdateMagnetic( g_FC_Mag_Half_1PG[0], g_FC_Mag_Half_1PG[1], g_FC_Mag_Half_1PG[2],
                             g_Mag_Array_In[P], g_EC_Ele_1PG, (real)0.5*dt, dh, N_HF_VAR, N_HF_ELE, 1 );
#        endif


//       1-a-4. evaluate the half-step solutions
         Hydro_RiemannPredict( g_Flu_Array_In[P], g_FC_Mag_Half_1PG, g_Flux_Half_1PG, g_PriVar_Half_1PG,
                               dt, dh, Gamma, MinDens, MinPres, NormPassive, NNorm, c_NormIdx,
                               JeansMinPres, JeansMinPres_Coeff );


//       1-a-5. evaluate the face-centered values by data reconstruction
//              --> note that g_PriVar_Half_1PG[] returned by Hydro_RiemannPredict() stores the primitive variables
         Hydro_DataReconstruction( NULL, g_FC_Mag_Half_1PG, g_PriVar_Half_1PG, g_FC_Var_1PG, g_Slope_PPM_1PG,
                                   Con2Pri_No, N_HF_VAR, LR_GHOST_SIZE, Gamma, LR_Limiter, MinMod_Coeff, dt, dh,
                                   MinDens, MinPres, NormPassive, NNorm, c_NormIdx, JeansMinPres, JeansMinPres_Coeff );


//       1-b. MHM: use interpolated face-centered values to calculate the half-step fluxes
#        elif ( FLU_SCHEME == MHM )

//       evaluate the face-centered values by data reconstruction
         Hydro_DataReconstruction( g_Flu_Array_In[P], NULL, g_PriVar_1PG, g_FC_Var_1PG, g_Slope_PPM_1PG,
                                   Con2Pri_Yes, FLU_NXT, LR_GHOST_SIZE, Gamma, LR_Limiter, MinMod_Coeff, dt, dh,
                                   MinDens, MinPres, NormPassive, NNorm, c_NormIdx, JeansMinPres, JeansMinPres_Coeff );

#        endif // #if ( FLU_SCHEME == MHM_RP ) ... else ...


//       2. evaluate the full-step fluxes
#        ifdef MHD
         const int NSkip_N = 0;
         const int NSkip_T = 0;
#        else
         const int NSkip_N = 0;
         const int NSkip_T = 1;
#        endif
         Hydro_ComputeFlux( g_FC_Var_1PG, g_FC_Flux_1PG, N_FL_FLUX, NSkip_N, NSkip_T, Gamma,
                            CorrHalfVel, g_Pot_Array_USG[P], g_Corner_Array[P],
                            dt, dh, Time, GravityType, c_ExtAcc_AuxArray, MinPres,
                            StoreFlux, g_Flux_Array[P] );


//       3. evaluate electric field and update B field at the full time-step
//          --> must update B field before Hydro_FullStepUpdate() since the latter requires
//              the updated magnetic energy when adopting the dual-energy formalism
#        ifdef MHD
         MHD_ComputeElectric( g_EC_Ele_1PG, g_FC_Flux_1PG, g_PriVar_Half_1PG, N_FL_ELE, N_FL_FLUX,
                              N_HF_VAR, LR_GHOST_SIZE, dt, dh, StoreElectric, g_Ele_Array[P],
                              CorrHalfVel, g_Pot_Array_USG[P], g_Corner_Array[P],
                              Time, GravityType, c_ExtAcc_AuxArray );

         MHD_UpdateMagnetic( g_Mag_Array_Out[P][0], g_Mag_Array_Out[P][1], g_Mag_Array_Out[P][2],
                             g_Mag_Array_In[P], g_EC_Ele_1PG, dt, dh, PS2, N_FL_ELE, FLU_GHOST_SIZE );
#        endif


//       4. full-step evolution
         Hydro_FullStepUpdate( g_Flu_Array_In[P], g_Flu_Array_Out[P], g_DE_Array_Out[P], g_Mag_Array_Out[P],
                               g_FC_Flux_1PG, dt, dh, Gamma, MinDens, MinPres, DualEnergySwitch,
                               NormPassive, NNorm, c_NormIdx );

      } // loop over all patch groups
   } // OpenMP parallel region

} // FUNCTION : CPU_FluidSolver_MHM



#if ( FLU_SCHEME == MHM_RP )
//-------------------------------------------------------------------------------------------------------
// Function    :  Hydro_RiemannPredict_Flux
// Description :  Evaluate the half-step face-centered fluxes by Riemann solver
//
// Note        :  1. Work for the MHM_RP scheme
//                2. Currently support the exact, Roe, HLLE, and HLLC solvers
//                3. g_Flux_Half[] is accessed with a stride N_HF_FLUX
//                   --> Fluxes on the **left** face of the (i+1,j+1,k+1) element in g_ConVar[] will
//                       be stored in the (i,j,k) element of g_Flux_Half[]
//
// Parameter   :  g_ConVar    : Array storing the input conserved variables
//                g_Flux_Half : Array to store the output face-centered fluxes
//                g_FC_B      : Array storing the input face-centered magnetic field (for MHD only)
//                              --> Accessed with strides FLU_NXT/FLU_NXT+1 along the
//                                  transverse/longitudinal directions
//                g_CC_B      : Array storing the input cell-centered magnetic field (for MHD only)
//                              --> Accessed with a stride FLU_NXT
//                Gamma       : Ratio of specific heats
//                MinPres     : Minimum allowed pressure
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE
void Hydro_RiemannPredict_Flux( const real g_ConVar[][ CUBE(FLU_NXT) ],
                                      real g_Flux_Half[][NCOMP_TOTAL_PLUS_MAG][ CUBE(N_FC_FLUX) ],
                                const real g_FC_B[][ SQR(FLU_NXT)*FLU_NXT_P1 ],
                                const real g_CC_B[][ CUBE(FLU_NXT) ],
                                const real Gamma, const real MinPres )
{

   const int didx_cvar[3] = { 1, FLU_NXT, SQR(FLU_NXT) };
   real ConVar_L[NCOMP_TOTAL_PLUS_MAG], ConVar_R[NCOMP_TOTAL_PLUS_MAG], Flux_1Face[NCOMP_TOTAL_PLUS_MAG];

#  if ( RSOLVER == EXACT )
   const real Gamma_m1 = Gamma - (real)1.0;
   real PriVar_L[NCOMP_TOTAL], PriVar_R[NCOMP_TOTAL];
#  endif


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

//       correct total energy by the difference between the face- and cell-centered logitudinal B field
         ConVar_L[4] += (real)0.5*(  SQR( ConVar_L[MAG_OFFSET + d] ) - SQR( g_CC_B[d][idx_cvar               ] )  );
         ConVar_R[4] += (real)0.5*(  SQR( ConVar_R[MAG_OFFSET + d] ) - SQR( g_CC_B[d][idx_cvar + didx_cvar[d]] )  );
#        endif

//       invoke the Riemann solver
#        if   ( RSOLVER == EXACT  &&  !defined MHD )
         const bool NormPassive_No  = false;  // do NOT convert any passive variable to mass fraction for the Riemann solvers
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
//                Gamma              : Ratio of specific heats
//                MinDens/Pres       : Minimum allowed density and pressure
//                NormPassive        : true --> convert passive scalars to mass fraction
//                NNorm              : Number of passive scalars for the option "NormPassive"
//                                     --> Should be set to the global variable "PassiveNorm_NVar"
//                NormIdx            : Target variable indices for the option "NormPassive"
//                                     --> Should be set to the global variable "PassiveNorm_VarIdx"
//                JeansMinPres       : Apply minimum pressure estimated from the Jeans length
//                JeansMinPres_Coeff : Coefficient used by JeansMinPres = G*(Jeans_NCell*Jeans_dh)^2/(Gamma*pi);
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE
void Hydro_RiemannPredict( const real g_ConVar_In[][ CUBE(FLU_NXT) ],
                           const real g_FC_B_Half[][ FLU_NXT_P1*SQR(FLU_NXT) ],
                           const real g_Flux_Half[][NCOMP_TOTAL_PLUS_MAG][ CUBE(N_FC_FLUX) ],
                                 real g_PriVar_Half[][ CUBE(FLU_NXT) ],
                           const real dt, const real dh, const real Gamma, const real MinDens, const real MinPres,
                           const bool NormPassive, const int NNorm, const int NormIdx[],
                           const bool JeansMinPres, const real JeansMinPres_Coeff )
{

   const int  didx_flux[3] = { 1, N_HF_FLUX, SQR(N_HF_FLUX) };
   const real dt_dh2       = (real)0.5*dt/dh;
   const real  Gamma_m1    = Gamma - (real)1.0;
   const real _Gamma_m1    = (real)1.0 / Gamma_m1;

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

//    compute the cell-centered half-step B field
#     ifdef MHD
      MHD_GetCellCenteredBField( out_con+MAG_OFFSET, g_FC_B_Half[0], g_FC_B_Half[1], g_FC_B_Half[2],
                                 N_HF_VAR, N_HF_VAR, N_HF_VAR, i_out, j_out, k_out );
#     endif

//    ensure positive density and pressure
#     ifdef MHD
      const real EngyB = (real)0.5*( SQR(out_con[MAG_OFFSET+0]) + SQR(out_con[MAG_OFFSET+1]) + SQR(out_con[MAG_OFFSET+2]) );
#     else
      const real EngyB = NULL_REAL;
#     endif
      out_con[0] = FMAX( out_con[0], MinDens );
      out_con[4] = Hydro_CheckMinPresInEngy( out_con[0], out_con[1], out_con[2], out_con[3], out_con[4], Gamma_m1, _Gamma_m1, MinPres, EngyB );
#     if ( NCOMP_PASSIVE > 0 )
      for (int v=NCOMP_FLUID; v<NCOMP_TOTAL; v++)
      out_con[v] = FMAX( out_con[v], TINY_NUMBER );
#     endif

//    conserved --> primitive variables
      Hydro_Con2Pri( out_con, out_pri, Gamma_m1, MinPres, NormPassive, NNorm, NormIdx, JeansMinPres, JeansMinPres_Coeff );

//    store the results in g_PriVar_Half[]
      for (int v=0; v<NCOMP_TOTAL_PLUS_MAG; v++)   g_PriVar_Half[v][idx_out] = out_pri[v];
   } // i,j,k


#  ifdef __CUDACC__
   __syncthreads();
#  endif

} // FUNCTION : Hydro_RiemannPredict
#endif // #if ( FLU_SCHEME == MHM_RP )



#endif // #if (  MODEL == HYDRO  &&  ( FLU_SCHEME == MHM || FLU_SCHEME == MHM_RP )  )
