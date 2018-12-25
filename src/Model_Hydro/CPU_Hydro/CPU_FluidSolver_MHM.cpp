#include "CUFLU.h"

#if (  MODEL == HYDRO  &&  ( FLU_SCHEME == MHM || FLU_SCHEME == MHM_RP )  )



// external functions
#ifdef __CUDACC__

#include "CUFLU_Shared_FluUtility.cu"
#include "CUFLU_Shared_DataReconstruction.cu"
#include "CUFLU_Shared_ComputeFlux.cu"
#include "CUFLU_Shared_FullStepUpdate.cu"

#if   ( RSOLVER == EXACT )
# include "CUFLU_Shared_RiemannSolver_Exact.cu"
#elif ( RSOLVER == ROE )
# include "CUFLU_Shared_RiemannSolver_Roe.cu"
#elif ( RSOLVER == HLLE )
# include "CUFLU_Shared_RiemannSolver_HLLE.cu"
#elif ( RSOLVER == HLLC )
# include "CUFLU_Shared_RiemannSolver_HLLC.cu"
#endif

#include "CUFLU_SetConstMem_FluidSolver.cu"

#else // #ifdef __CUDACC__

void Hydro_DataReconstruction( const real ConVar   [][ CUBE(FLU_NXT) ],
                                     real PriVar   [][ CUBE(FLU_NXT) ],
                                     real FC_Var   [][NCOMP_TOTAL][ CUBE(N_FC_VAR) ],
                                     real Slope_PPM[][NCOMP_TOTAL][ CUBE(N_SLOPE_PPM) ],
                               const bool Con2Pri, const int NIn, const int NGhost, const real Gamma,
                               const LR_Limiter_t LR_Limiter, const real MinMod_Coeff,
                               const real dt, const real dh, const real MinDens, const real MinPres,
                               const bool NormPassive, const int NNorm, const int NormIdx[],
                               const bool JeansMinPres, const real JeansMinPres_Coeff );
void Hydro_ComputeFlux( const real FC_Var [][NCOMP_TOTAL][ CUBE(N_FC_VAR)],
                              real FC_Flux[][NCOMP_TOTAL][ CUBE(N_FC_FLUX) ],
                        const int Gap, const real Gamma, const bool CorrHalfVel, const real Pot_USG[],
                        const double Corner[], const real dt, const real dh, const double Time,
                        const OptGravityType_t GravityType, const double ExtAcc_AuxArray[], const real MinPres,
                        const bool DumpIntFlux, real IntFlux[][NCOMP_TOTAL][ SQR(PS2) ] );
void Hydro_FullStepUpdate( const real Input[][ CUBE(FLU_NXT) ], real Output[][ CUBE(PS2) ], char DE_Status[],
                           const real Flux[][NCOMP_TOTAL][ CUBE(N_FC_FLUX) ], const real dt, const real dh,
                           const real Gamma_m1, const real _Gamma_m1, const real MinDens, const real MinPres,
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
#endif
#if ( FLU_SCHEME == MHM_RP )
void Hydro_Con2Pri( const real In[], real Out[], const real Gamma_m1, const real MinPres,
                    const bool NormPassive, const int NNorm, const int NormIdx[],
                    const bool JeansMinPres, const real JeansMinPres_Coeff );
real Hydro_CheckMinPresInEngy( const real Dens, const real MomX, const real MomY, const real MomZ, const real Engy,
                               const real Gamma_m1, const real _Gamma_m1, const real MinPres );
#endif

#endif // #ifdef __CUDACC__ ... else ...


// internal functions
#if ( FLU_SCHEME == MHM_RP )
GPU_DEVICE
static void Hydro_RiemannPredict_Flux( const real ConVar[][ CUBE(FLU_NXT) ], real Half_Flux[][NCOMP_TOTAL][ CUBE(N_FC_FLUX) ],
                                       const real Gamma, const real MinPres );
GPU_DEVICE
static void Hydro_RiemannPredict( const real ConVar_In[][ CUBE(FLU_NXT) ],
                                  const real Half_Flux[][NCOMP_TOTAL][ CUBE(N_FC_FLUX) ],
                                        real Half_Var [][ CUBE(FLU_NXT) ],
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
//
// Parameter   :  Flu_Array_In       : Array storing the input fluid variables
//                Flu_Array_Out      : Array to store the output fluid variables
//                DE_Array_Out       : Array to store the dual-energy status
//                Flux_Array         : Array to store the output fluxes
//                Corner_Array       : Array storing the physical corner coordinates of each patch group (for UNSPLIT_GRAVITY)
//                Pot_Array_USG      : Array storing the input potential for UNSPLIT_GRAVITY
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
void CUFLU_FluidSolver_MHM(
   const real   Flu_Array_In [][NCOMP_TOTAL][ CUBE(FLU_NXT) ],
         real   Flu_Array_Out[][NCOMP_TOTAL][ CUBE(PS2) ],
         char   DE_Array_Out [][ CUBE(PS2) ],
         real   Flux_Array   [][9][NCOMP_TOTAL][ SQR(PS2) ],
   const double Corner_Array [][3],
   const real   Pot_Array_USG[][ CUBE(USG_NXT_F) ],
         real   PriVar       [][NCOMP_TOTAL][ CUBE(FLU_NXT) ],
         real   Slope_PPM    [][3][NCOMP_TOTAL][ CUBE(N_SLOPE_PPM) ],
         real   FC_Var       [][6][NCOMP_TOTAL][ CUBE(N_FC_VAR) ],
         real   FC_Flux      [][3][NCOMP_TOTAL][ CUBE(N_FC_FLUX) ],
   const real dt, const real dh, const real Gamma, const bool StoreFlux,
   const LR_Limiter_t LR_Limiter, const real MinMod_Coeff,
   const double Time, const OptGravityType_t GravityType,
   const real MinDens, const real MinPres, const real DualEnergySwitch,
   const bool NormPassive, const int NNorm,
   const bool JeansMinPres, const real JeansMinPres_Coeff )
#else
void CPU_FluidSolver_MHM(
   const real   Flu_Array_In [][NCOMP_TOTAL][ CUBE(FLU_NXT) ],
         real   Flu_Array_Out[][NCOMP_TOTAL][ CUBE(PS2) ],
         char   DE_Array_Out [][ CUBE(PS2) ],
         real   Flux_Array   [][9][NCOMP_TOTAL][ SQR(PS2) ],
   const double Corner_Array [][3],
   const real   Pot_Array_USG[][ CUBE(USG_NXT_F) ],
   const int NPatchGroup, const real dt, const real dh, const real Gamma,
   const bool StoreFlux, const LR_Limiter_t LR_Limiter, const real MinMod_Coeff,
   const double Time, const OptGravityType_t GravityType,
   const double c_ExtAcc_AuxArray[], const real MinDens, const real MinPres,
   const real DualEnergySwitch, const bool NormPassive, const int NNorm, const int c_NormIdx[],
   const bool JeansMinPres, const real JeansMinPres_Coeff )
#endif // #ifdef __CUDACC__ ... else ...
{

   const real  Gamma_m1       = Gamma - (real)1.0;
   const real _Gamma_m1       = (real)1.0 / Gamma_m1;
#  ifdef UNSPLIT_GRAVITY
   const bool CorrHalfVel_Yes = true;
#  else
   const bool CorrHalfVel_No  = false;
#  endif
#  if   ( FLU_SCHEME == MHM )
   const bool Con2Pri_Yes     = true;
#  elif ( FLU_SCHEME == MHM_RP )
   const bool Con2Pri_No      = false;
#  endif


// openmp pragma for the CPU solver
#  ifndef __CUDACC__
#  pragma omp parallel
#  endif
   {
//    FC: Face-Centered variables/fluxes
//    --> CPU solver: allocate data
//        GPU solver: link to the input global memory array
#     ifdef __CUDACC__
      real (*FC_Var_1PG   )[NCOMP_TOTAL][ CUBE(N_FC_VAR)    ] = FC_Var   [blockIdx.x];
      real (*FC_Flux_1PG  )[NCOMP_TOTAL][ CUBE(N_FC_FLUX)   ] = FC_Flux  [blockIdx.x];
      real (*PriVar_1PG   )             [ CUBE(FLU_NXT)     ] = PriVar   [blockIdx.x];
#     if ( LR_SCHEME == PPM )
      real (*Slope_PPM_1PG)[NCOMP_TOTAL][ CUBE(N_SLOPE_PPM) ] = Slope_PPM[blockIdx.x];
#     else
      real (*Slope_PPM_1PG)[NCOMP_TOTAL][ CUBE(N_SLOPE_PPM) ] = NULL;
#     endif

#     else // #ifdef __CUDACC__

      real (*FC_Var_1PG   )[NCOMP_TOTAL][ CUBE(N_FC_VAR)    ] = new real [6][NCOMP_TOTAL][ CUBE(N_FC_VAR)    ];
      real (*FC_Flux_1PG  )[NCOMP_TOTAL][ CUBE(N_FC_FLUX)   ] = new real [3][NCOMP_TOTAL][ CUBE(N_FC_FLUX)   ];   // also used by "Half_Flux"
      real (*PriVar_1PG   )             [ CUBE(FLU_NXT)     ] = new real    [NCOMP_TOTAL][ CUBE(FLU_NXT)     ];   // also used by "Half_Var"
#     if ( LR_SCHEME == PPM )
      real (*Slope_PPM_1PG)[NCOMP_TOTAL][ CUBE(N_SLOPE_PPM) ] = new real [3][NCOMP_TOTAL][ CUBE(N_SLOPE_PPM) ];
#     else
      real (*Slope_PPM_1PG)[NCOMP_TOTAL][ CUBE(N_SLOPE_PPM) ] = NULL;
#     endif
#     endif // #ifdef __CUDACC__ ... else ...

#     if ( FLU_SCHEME == MHM_RP )
      real (*const Half_Flux_1PG)[NCOMP_TOTAL][ CUBE(N_FC_FLUX) ] = FC_Flux_1PG;
      real (*const Half_Var_1PG )             [ CUBE(FLU_NXT)   ] = PriVar_1PG;
#     endif


//    loop over all patch groups
//    --> CPU/GPU solver: use different (OpenMP threads) / (CUDA thread blocks)
//                        to work on different patch groups
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

//       1-a-1. evaluate the half-step first-order fluxes by Riemann solver
         Hydro_RiemannPredict_Flux( Flu_Array_In[P], Half_Flux_1PG, Gamma, MinPres );


//       1-a-2. evaluate the half-step solutions
         Hydro_RiemannPredict( Flu_Array_In[P], Half_Flux_1PG, Half_Var_1PG, dt, dh, Gamma, MinDens, MinPres,
                               NormPassive, NNorm, c_NormIdx, JeansMinPres, JeansMinPres_Coeff );


//       1-a-3. evaluate the face-centered values by data reconstruction
//              --> note that Half_Var_1PG[] returned by Hydro_RiemannPredict() stores the primitive variables
         Hydro_DataReconstruction( NULL, Half_Var_1PG, FC_Var_1PG, Slope_PPM_1PG,
                                   Con2Pri_No, N_HF_VAR, FLU_GHOST_SIZE-2,
                                   Gamma, LR_Limiter, MinMod_Coeff, dt, dh, MinDens, MinPres,
                                   NormPassive, NNorm, c_NormIdx, JeansMinPres, JeansMinPres_Coeff );


//       1-b. MHM: use interpolated face-centered values to calculate the half-step fluxes
#        elif ( FLU_SCHEME == MHM )

//       evaluate the face-centered values by data reconstruction
         Hydro_DataReconstruction( Flu_Array_In[P], PriVar_1PG, FC_Var_1PG, Slope_PPM_1PG,
                                   Con2Pri_Yes, FLU_NXT, FLU_GHOST_SIZE-1,
                                   Gamma, LR_Limiter, MinMod_Coeff, dt, dh, MinDens, MinPres,
                                   NormPassive, NNorm, c_NormIdx, JeansMinPres, JeansMinPres_Coeff );
#        endif // #if ( FLU_SCHEME == MHM_RP ) ... else ...


//       2. evaluate the full-step fluxes
#        ifdef UNSPLIT_GRAVITY
         Hydro_ComputeFlux( FC_Var_1PG, FC_Flux_1PG, 1, Gamma, CorrHalfVel_Yes,
                            Pot_Array_USG[P], Corner_Array[P],
                            dt, dh, Time, GravityType, c_ExtAcc_AuxArray, MinPres,
                            StoreFlux, Flux_Array[P] );
#        else
         Hydro_ComputeFlux( FC_Var_1PG, FC_Flux_1PG, 1, Gamma, CorrHalfVel_No,
                            NULL, NULL,
                            NULL_REAL, NULL_REAL, NULL_REAL, GRAVITY_NONE, NULL, MinPres,
                            StoreFlux, Flux_Array[P] );
#        endif


//       3. full-step evolution
         Hydro_FullStepUpdate( Flu_Array_In[P], Flu_Array_Out[P], DE_Array_Out[P],
                               FC_Flux_1PG, dt, dh, Gamma_m1, _Gamma_m1, MinDens, MinPres, DualEnergySwitch,
                               NormPassive, NNorm, c_NormIdx );

      } // loop over all patch groups

#     ifndef __CUDACC__
      delete [] FC_Var_1PG;
      delete [] FC_Flux_1PG;
      delete [] PriVar_1PG;
      delete [] Slope_PPM_1PG;
#     endif

   } // OpenMP parallel region

} // FUNCTION : CPU_FluidSolver_MHM



#if ( FLU_SCHEME == MHM_RP )
//-------------------------------------------------------------------------------------------------------
// Function    :  Hydro_RiemannPredict_Flux
// Description :  Evaluate the half-step face-centered fluxes by Riemann solver
//
// Note        :  1. Work for the MHM_RP scheme
//                2. Currently support the exact, Roe, HLLE, and HLLC solvers
//                3. Half_Flux[] is accessed with the stride N_FC_FLUX
//                   --> Fluxes on the **left** face of the (i+1,j+1,k+1) element in ConVar[] will
//                       be stored in the (i,j,k) element of Half_Flux[]
//
// Parameter   :  ConVar    : Array storing the input conserved variables
//                Half_Flux : Array to store the output face-centered fluxes
//                Gamma     : Ratio of specific heats
//                MinPres   : Minimum allowed pressure
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE
void Hydro_RiemannPredict_Flux( const real ConVar[][ CUBE(FLU_NXT) ], real Half_Flux[][NCOMP_TOTAL][ CUBE(N_FC_FLUX) ],
                                const real Gamma, const real MinPres )
{

   const int didx_cvar[3] = { 1, FLU_NXT, SQR(FLU_NXT) };
   real ConVar_L[NCOMP_TOTAL], ConVar_R[NCOMP_TOTAL], Flux_1Face[NCOMP_TOTAL];

#  if ( RSOLVER == EXACT )
   const real Gamma_m1 = Gamma - (real)1.0;
   real PriVar_L[NCOMP_TOTAL], PriVar_R[NCOMP_TOTAL];
#  endif


// loop over different spatial directions
   for (int d=0; d<3; d++)
   {
      int gap[3];

      switch ( d )
      {
         case 0 : gap[0] = 0;  gap[1] = 1;  gap[2] = 1;  break;
         case 1 : gap[0] = 1;  gap[1] = 0;  gap[2] = 1;  break;
         case 2 : gap[0] = 1;  gap[1] = 1;  gap[2] = 0;  break;
      }

      const int size_i  = ( N_FC_FLUX - gap[0] );
      const int size_ij = ( N_FC_FLUX - gap[1] )*size_i;

      CGPU_LOOP( idx, N_FC_FLUX*SQR(N_FC_FLUX-1) )
      {
         const int i_flux   = idx % size_i;
         const int j_flux   = idx % size_ij / size_i;
         const int k_flux   = idx / size_ij;
         const int idx_flux = IDX321( i_flux, j_flux, k_flux, N_FC_FLUX, N_FC_FLUX );

         const int i_cvar   = i_flux + gap[0];
         const int j_cvar   = j_flux + gap[1];
         const int k_cvar   = k_flux + gap[2];
         const int idx_cvar = IDX321( i_cvar, j_cvar, k_cvar, FLU_NXT, FLU_NXT );

//       get the left and right states
         for (int v=0; v<NCOMP_TOTAL; v++)
         {
            ConVar_L[v] = ConVar[v][ idx_cvar              ];
            ConVar_R[v] = ConVar[v][ idx_cvar+didx_cvar[d] ];
         }

//       invoke the Riemann solver
#        if   ( RSOLVER == EXACT )
         const bool NormPassive_No  = false;  // do NOT convert any passive variable to mass fraction for the Riemann solvers
         const bool JeansMinPres_No = false;

         Hydro_Con2Pri( ConVar_L, PriVar_L, Gamma_m1, MinPres, NormPassive_No, NULL_INT, NULL, JeansMinPres_No, NULL_REAL );
         Hydro_Con2Pri( ConVar_R, PriVar_R, Gamma_m1, MinPres, NormPassive_No, NULL_INT, NULL, JeansMinPres_No, NULL_REAL );

         Hydro_RiemannSolver_Exact( d, Flux_1Face, PriVar_L, PriVar_R, Gamma );
#        elif ( RSOLVER == ROE )
         Hydro_RiemannSolver_Roe  ( d, Flux_1Face, ConVar_L, ConVar_R, Gamma, MinPres );
#        elif ( RSOLVER == HLLE )
         Hydro_RiemannSolver_HLLE ( d, Flux_1Face, ConVar_L, ConVar_R, Gamma, MinPres );
#        elif ( RSOLVER == HLLC )
         Hydro_RiemannSolver_HLLC ( d, Flux_1Face, ConVar_L, ConVar_R, Gamma, MinPres );
#        else
#        error : ERROR : unsupported Riemann solver (EXACT/ROE) !!
#        endif

//       store the results in Half_Flux[]
         for (int v=0; v<NCOMP_TOTAL; v++)   Half_Flux[d][v][idx_flux] = Flux_1Face[v];
      } // CGPU_LOOP( idx, N_FC_FLUX*SQR(N_FC_FLUX-1) )
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
//
// Parameter   :  ConVar_In          : Array storing the input conserved variables
//                Half_Flux          : Array storing the input face-centered fluxes
//                                     --> Accessed with the stride N_FC_FLUX
//                Half_Var           : Array to store the output primitive variables
//                                     --> Accessed with the stride N_HF_VAR
//                                     --> Although its actually allocated size is FLU_NXT^3 since it points to PriVar_1PG[]
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
void Hydro_RiemannPredict( const real ConVar_In[][ CUBE(FLU_NXT) ],
                           const real Half_Flux[][NCOMP_TOTAL][ CUBE(N_FC_FLUX) ],
                                 real Half_Var [][ CUBE(FLU_NXT) ],
                           const real dt, const real dh, const real Gamma, const real MinDens, const real MinPres,
                           const bool NormPassive, const int NNorm, const int NormIdx[],
                           const bool JeansMinPres, const real JeansMinPres_Coeff )
{

   const int  didx_flux[3] = { 1, N_FC_FLUX, SQR(N_FC_FLUX) };
   const real dt_dh2       = (real)0.5*dt/dh;
   const real  Gamma_m1    = Gamma - (real)1.0;
   const real _Gamma_m1    = (real)1.0 / Gamma_m1;

   const int N_HF_VAR2 = SQR(N_HF_VAR);
   CGPU_LOOP( idx_out, CUBE(N_HF_VAR) )
   {
      const int i_flux   = idx_out % N_HF_VAR;
      const int j_flux   = idx_out % N_HF_VAR2 / N_HF_VAR;
      const int k_flux   = idx_out / N_HF_VAR2;
      const int idx_flux = IDX321( i_flux, j_flux, k_flux, N_FC_FLUX, N_FC_FLUX );

      const int i_in     = i_flux + 1;
      const int j_in     = j_flux + 1;
      const int k_in     = k_flux + 1;
      const int idx_in   = IDX321( i_in, j_in, k_in, FLU_NXT, FLU_NXT );

      real out_con[NCOMP_TOTAL], out_pri[NCOMP_TOTAL], dflux[3][NCOMP_TOTAL];

//    calculate the flux differences
      for (int d=0; d<3; d++)
      for (int v=0; v<NCOMP_TOTAL; v++)    dflux[d][v] = Half_Flux[d][v][ idx_flux+didx_flux[d] ] - Half_Flux[d][v][idx_flux];

//    update the input cell-centered conserved variables with the flux differences
      for (int v=0; v<NCOMP_TOTAL; v++)
         out_con[v] = ConVar_In[v][idx_in] - dt_dh2*( dflux[0][v] + dflux[1][v] + dflux[2][v] );

//    ensure positive density and pressure
      out_con[0] = FMAX( out_con[0], MinDens );
      out_con[4] = Hydro_CheckMinPresInEngy( out_con[0], out_con[1], out_con[2], out_con[3], out_con[4], Gamma_m1, _Gamma_m1, MinPres );
#     if ( NCOMP_PASSIVE > 0 )
      for (int v=NCOMP_FLUID; v<NCOMP_TOTAL; v++)
      out_con[v] = FMAX( out_con[v], TINY_NUMBER );
#     endif

//    conserved --> primitive variables
      Hydro_Con2Pri( out_con, out_pri, Gamma_m1, MinPres, NormPassive, NNorm, NormIdx, JeansMinPres, JeansMinPres_Coeff );

//    store the results to Half_Var[]
      for (int v=0; v<NCOMP_TOTAL; v++)   Half_Var[v][idx_out] = out_pri[v];
   } // i,j,k


#  ifdef __CUDACC__
   __syncthreads();
#  endif

} // FUNCTION : Hydro_RiemannPredict
#endif // #if ( FLU_SCHEME == MHM_RP )



#endif // #if (  MODEL == HYDRO  &&  ( FLU_SCHEME == MHM || FLU_SCHEME == MHM_RP )  )
