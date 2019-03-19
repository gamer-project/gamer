#include "CUFLU.h"

#if ( MODEL == HYDRO  &&  FLU_SCHEME == CTU )



// external functions
#ifdef __CUDACC__

#include "CUFLU_Shared_FluUtility.cu"
#include "CUFLU_Shared_DataReconstruction.cu"
#include "CUFLU_Shared_ComputeFlux.cu"
#include "CUFLU_Shared_FullStepUpdate.cu"
#include "CUFLU_SetConstMem_FluidSolver.cu"

#else // #ifdef __CUDACC__

void Hydro_DataReconstruction( const real g_ConVar   [][ CUBE(FLU_NXT) ],
                                     real g_PriVar   [][ CUBE(FLU_NXT) ],
                                     real g_FC_Var   [][NCOMP_TOTAL][ CUBE(N_FC_VAR) ],
                                     real g_Slope_PPM[][NCOMP_TOTAL][ CUBE(N_SLOPE_PPM) ],
                               const bool Con2Pri, const int NIn, const int NGhost, const real Gamma,
                               const LR_Limiter_t LR_Limiter, const real MinMod_Coeff,
                               const real dt, const real dh, const real MinDens, const real MinPres,
                               const bool NormPassive, const int NNorm, const int NormIdx[],
                               const bool JeansMinPres, const real JeansMinPres_Coeff );
void Hydro_ComputeFlux( const real g_FC_Var [][NCOMP_TOTAL][ CUBE(N_FC_VAR) ],
                              real g_FC_Flux[][NCOMP_TOTAL][ CUBE(N_FC_FLUX) ],
                        const int Gap, const real Gamma, const bool CorrHalfVel, const real g_Pot_USG[],
                        const double g_Corner[], const real dt, const real dh, const double Time,
                        const OptGravityType_t GravityType, const double ExtAcc_AuxArray[], const real MinPres,
                        const bool DumpIntFlux, real g_IntFlux[][NCOMP_TOTAL][ SQR(PS2) ] );
void Hydro_FullStepUpdate( const real g_Input[][ CUBE(FLU_NXT) ], real g_Output[][ CUBE(PS2) ], char g_DE_Status[],
                           const real g_Flux[][NCOMP_TOTAL][ CUBE(N_FC_FLUX) ], const real dt, const real dh,
                           const real Gamma, const real MinDens, const real MinPres, const real DualEnergySwitch,
                           const bool NormPassive, const int NNorm, const int NormIdx[] );
real Hydro_CheckMinPresInEngy( const real Dens, const real MomX, const real MomY, const real MomZ, const real Engy,
                               const real Gamma_m1, const real _Gamma_m1, const real MinPres );

#endif // #ifdef __CUDACC__ ... else ...


// internal functions
GPU_DEVICE
static void Hydro_TGradient_Correction(       real g_FC_Var [][NCOMP_TOTAL][ CUBE(N_FC_VAR)  ],
                                        const real g_FC_Flux[][NCOMP_TOTAL][ CUBE(N_FC_FLUX) ],
                                        const real dt, const real dh, const real Gamma,
                                        const real MinDens, const real MinPres );




//-------------------------------------------------------------------------------------------------------
// Function    :  CPU/CUFLU_FluidSolver_CTU
// Description :  CPU/GPU fluid solver based on the Corner-Transport-Upwind (CTU) scheme
//
// Note        :  1. Ref: Stone et al., ApJS, 178, 137 (2008)
//                2. See include/CUFLU.h for the values and description of different symbolic constants
//                   such as N_FC_VAR, N_FC_FLUX, N_SLOPE_PPM, N_FL_FLUX, N_HF_VAR
//                3. Arrays with a prefix "g_" are stored in the global memory of GPU
//
// Parameter   :  g_Flu_Array_In     : Array storing the input fluid variables
//                g_Flu_Array_Out    : Array to store the output fluid variables
//                g_DE_Array_Out     : Array to store the dual-energy status
//                g_Flux_Array       : Array to store the output fluxes
//                g_Corner_Array     : Array storing the physical corner coordinates of each patch group (for UNSPLIT_GRAVITY)
//                g_Pot_Array_USG    : Array storing the input potential for UNSPLIT_GRAVITY
//                g_PriVar           : Array to store the primitive variables
//                g_Slope_PPM        : Array to store the slope for the PPM reconstruction
//                g_FC_Var           : Array to store the half-step variables
//                g_FC_Flux          : Array to store the face-centered fluxes
//                NPatchGroup        : Number of patch groups to be evaluated
//                dt                 : Time interval to advance solution
//                dh                 : Cell size
//                Gamma              : Ratio of specific heats
//                StoreFlux          : true --> store the coarse-fine fluxes
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
void CUFLU_FluidSolver_CTU(
   const real   g_Flu_Array_In [][NCOMP_TOTAL][ CUBE(FLU_NXT) ],
         real   g_Flu_Array_Out[][NCOMP_TOTAL][ CUBE(PS2) ],
         char   g_DE_Array_Out [][ CUBE(PS2) ],
         real   g_Flux_Array   [][9][NCOMP_TOTAL][ SQR(PS2) ],
   const double g_Corner_Array [][3],
   const real   g_Pot_Array_USG[][ CUBE(USG_NXT_F) ],
         real   g_PriVar       [][NCOMP_TOTAL][ CUBE(FLU_NXT) ],
         real   g_Slope_PPM    [][3][NCOMP_TOTAL][ CUBE(N_SLOPE_PPM) ],
         real   g_FC_Var       [][6][NCOMP_TOTAL][ CUBE(N_FC_VAR) ],
         real   g_FC_Flux      [][3][NCOMP_TOTAL][ CUBE(N_FC_FLUX) ],
   const real dt, const real dh, const real Gamma, const bool StoreFlux,
   const LR_Limiter_t LR_Limiter, const real MinMod_Coeff,
   const double Time, const OptGravityType_t GravityType,
   const real MinDens, const real MinPres, const real DualEnergySwitch,
   const bool NormPassive, const int NNorm,
   const bool JeansMinPres, const real JeansMinPres_Coeff )
#else
void CPU_FluidSolver_CTU(
   const real   g_Flu_Array_In [][NCOMP_TOTAL][ CUBE(FLU_NXT) ],
         real   g_Flu_Array_Out[][NCOMP_TOTAL][ CUBE(PS2) ],
         char   g_DE_Array_Out [][ CUBE(PS2) ],
         real   g_Flux_Array   [][9][NCOMP_TOTAL][ SQR(PS2) ],
   const double g_Corner_Array [][3],
   const real   g_Pot_Array_USG[][ CUBE(USG_NXT_F) ],
         real   g_PriVar       [][NCOMP_TOTAL][ CUBE(FLU_NXT) ],
         real   g_Slope_PPM    [][3][NCOMP_TOTAL][ CUBE(N_SLOPE_PPM) ],
         real   g_FC_Var       [][6][NCOMP_TOTAL][ CUBE(N_FC_VAR) ],
         real   g_FC_Flux      [][3][NCOMP_TOTAL][ CUBE(N_FC_FLUX) ],
   const int NPatchGroup, const real dt, const real dh, const real Gamma,
   const bool StoreFlux, const LR_Limiter_t LR_Limiter, const real MinMod_Coeff,
   const double Time, const OptGravityType_t GravityType,
   const double c_ExtAcc_AuxArray[], const real MinDens, const real MinPres,
   const real DualEnergySwitch, const bool NormPassive, const int NNorm, const int c_NormIdx[],
   const bool JeansMinPres, const real JeansMinPres_Coeff )
#endif // #ifdef __CUDACC__ ... else ...
{

#  ifdef UNSPLIT_GRAVITY
   const bool CorrHalfVel_Yes = true;
#  endif
   const bool CorrHalfVel_No  = false;
   const bool StoreFlux_No    = false;
   const bool Con2Pri_Yes     = true;


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

      real (*const g_FC_Var_1PG   )[NCOMP_TOTAL][ CUBE(N_FC_VAR)    ] = g_FC_Var   [array_idx];
      real (*const g_FC_Flux_1PG  )[NCOMP_TOTAL][ CUBE(N_FC_FLUX)   ] = g_FC_Flux  [array_idx];
      real (*const g_PriVar_1PG   )             [ CUBE(FLU_NXT)     ] = g_PriVar   [array_idx];
      real (*const g_Slope_PPM_1PG)[NCOMP_TOTAL][ CUBE(N_SLOPE_PPM) ] = g_Slope_PPM[array_idx];


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
//       1. evaluate the face-centered values at the half time-step
         Hydro_DataReconstruction( g_Flu_Array_In[P], g_PriVar_1PG, g_FC_Var_1PG, g_Slope_PPM_1PG,
                                   Con2Pri_Yes, FLU_NXT, FLU_GHOST_SIZE-1,
                                   Gamma, LR_Limiter, MinMod_Coeff, dt, dh, MinDens, MinPres,
                                   NormPassive, NNorm, c_NormIdx, JeansMinPres, JeansMinPres_Coeff );


//       2. evaluate the face-centered half-step fluxes by solving the Riemann problem
         Hydro_ComputeFlux( g_FC_Var_1PG, g_FC_Flux_1PG, 0, Gamma, CorrHalfVel_No, NULL, NULL,
                            NULL_REAL, NULL_REAL, NULL_REAL, GRAVITY_NONE, NULL, MinPres,
                            StoreFlux_No, NULL );


//       3. correct the face-centered variables by the transverse flux gradients
         Hydro_TGradient_Correction( g_FC_Var_1PG, g_FC_Flux_1PG, dt, dh, Gamma, MinDens, MinPres );


//       4. evaluate the face-centered full-step fluxes by solving the Riemann problem with the corrected data
#        ifdef UNSPLIT_GRAVITY
         Hydro_ComputeFlux( g_FC_Var_1PG, g_FC_Flux_1PG, 1, Gamma, CorrHalfVel_Yes,
                            g_Pot_Array_USG[P], g_Corner_Array[P],
                            dt, dh, Time, GravityType, c_ExtAcc_AuxArray, MinPres,
                            StoreFlux, g_Flux_Array[P] );
#        else
         Hydro_ComputeFlux( g_FC_Var_1PG, g_FC_Flux_1PG, 1, Gamma, CorrHalfVel_No,
                            NULL, NULL,
                            NULL_REAL, NULL_REAL, NULL_REAL, GRAVITY_NONE, NULL, MinPres,
                            StoreFlux, g_Flux_Array[P] );
#        endif


//       5. full-step evolution
         Hydro_FullStepUpdate( g_Flu_Array_In[P], g_Flu_Array_Out[P], g_DE_Array_Out[P],
                               g_FC_Flux_1PG, dt, dh, Gamma, MinDens, MinPres, DualEnergySwitch,
                               NormPassive, NNorm, c_NormIdx );

      } // loop over all patch groups
   } // OpenMP parallel region

} // FUNCTION : CPU_FluidSolver_CTU



//-------------------------------------------------------------------------------------------------------
// Function    :  Hydro_TGradient_Correction
// Description :  1. Correct the face-centered variables by the transverse flux gradients
//                2. Assuming "N_FC_VAR == N_FC_FLUX"
//
// Parameter   :  g_FC_Var     : Array to store the input and output face-centered conserved variables
//                               --> Accessed with the stride N_FC_VAR
//                g_FC_Flux    : Array storing the input face-centered fluxes
//                               --> Accessed with the stride N_FC_FLUX
//                dt           : Time interval to advance solution
//                dh           : Cell size
//                Gamma        : Ratio of specific heats
//                MinDens/Pres : Minimum allowed density and pressure
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE
void Hydro_TGradient_Correction(       real g_FC_Var [][NCOMP_TOTAL][ CUBE(N_FC_VAR)  ],
                                 const real g_FC_Flux[][NCOMP_TOTAL][ CUBE(N_FC_FLUX) ],
                                 const real dt, const real dh, const real Gamma,
                                 const real MinDens, const real MinPres )
{

   const int  NCell     = N_FC_VAR;    // size of g_FC_Var[] and g_FC_Flux[] in each direction
   const int  didx[3]   = { 1, NCell, SQR(NCell) };
   const real dt_dh2    = (real)0.5*dt/dh;
   const real  Gamma_m1 = Gamma - (real)1.0;
   const real _Gamma_m1 = (real)1.0 / Gamma_m1;

// loop over different spatial directions
   for (int d=0; d<3; d++)
   {
      const int faceL = 2*d;
      const int faceR = faceL+1;
      const int TDir1 = (d+1)%3;    // transverse direction 1
      const int TDir2 = (d+2)%3;    // transverse direction 2

      real fc_var[2][NCOMP_TOTAL];  // 0/1 = left/right faces
      int  gap[3];

      switch ( d )
      {
         case 0 : gap[0] = 0;   gap[1] = 1;   gap[2] = 1;   break;
         case 1 : gap[0] = 1;   gap[1] = 0;   gap[2] = 1;   break;
         case 2 : gap[0] = 1;   gap[1] = 1;   gap[2] = 0;   break;
      }

      const int size_i  = ( NCell - 2*gap[0] );
      const int size_ij = ( NCell - 2*gap[1] )*size_i;

      CGPU_LOOP( idx, NCell*SQR(NCell-2) )
      {
         const int i_var   = gap[0] + idx % size_i;
         const int j_var   = gap[1] + idx % size_ij / size_i;
         const int k_var   = gap[2] + idx / size_ij;
         const int idx_var = IDX321( i_var, j_var, k_var, NCell, NCell );

         const int idx_fluxR  = idx_var;
         const int idx_fluxL1 = idx_fluxR - didx[TDir1];
         const int idx_fluxL2 = idx_fluxR - didx[TDir2];

//       calculate the transverse flux gradients and update the corresponding face-centered variables
         for (int v=0; v<NCOMP_TOTAL; v++)
         {
            real Correct, TGrad1, TGrad2;

//          load g_FC_Var[] to a local variable fc[] to reduce the GPU global memory access
            fc_var[0][v] = g_FC_Var[faceL][v][idx_var];
            fc_var[1][v] = g_FC_Var[faceR][v][idx_var];

            TGrad1  = g_FC_Flux[TDir1][v][idx_fluxR] - g_FC_Flux[TDir1][v][idx_fluxL1];
            TGrad2  = g_FC_Flux[TDir2][v][idx_fluxR] - g_FC_Flux[TDir2][v][idx_fluxL2];
            Correct = -dt_dh2*( TGrad1 + TGrad2 );

            fc_var[0][v] += Correct;
            fc_var[1][v] += Correct;
         }

//       ensure positive density and pressure
         for (int f=0; f<2; f++)
         {
            fc_var[f][0] = FMAX( fc_var[f][0], MinDens );
            fc_var[f][4] = Hydro_CheckMinPresInEngy( fc_var[f][0], fc_var[f][1], fc_var[f][2], fc_var[f][3], fc_var[f][4],
                                                     Gamma_m1, _Gamma_m1, MinPres );
#           if ( NCOMP_PASSIVE > 0 )
            for (int v=NCOMP_FLUID; v<NCOMP_TOTAL; v++)
            fc_var[f][v] = FMAX( fc_var[f][v], TINY_NUMBER );
#           endif
         }

//       store the results to g_FC_Var[]
         for (int v=0; v<NCOMP_TOTAL; v++)
         {
            g_FC_Var[faceL][v][idx_var] = fc_var[0][v];
            g_FC_Var[faceR][v][idx_var] = fc_var[1][v];
         }

      } // CGPU_LOOP( idx, NCell*SQR(NCell-2) )
   } // for (int d=0; d<3; d++)


#  ifdef __CUDACC__
   __syncthreads();
#  endif

} // FUNCTION : Hydro_TGradient_Correction



#endif // #if ( !defined GPU  &&  MODEL == HYDRO  &&  FLU_SCHEME == CTU )
