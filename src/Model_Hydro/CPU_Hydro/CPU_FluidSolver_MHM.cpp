#include "GAMER.h"
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

#else // #ifdef __CUDACC__

extern void Hydro_DataReconstruction( const real PriVar[][ FLU_NXT*FLU_NXT*FLU_NXT    ],
                                      const real ConVar[][ FLU_NXT*FLU_NXT*FLU_NXT    ],
                                            real FC_Var[][NCOMP_TOTAL][ N_FC_VAR*N_FC_VAR*N_FC_VAR ],
                                      const int NIn, const int NGhost, const real Gamma,
                                      const LR_Limiter_t LR_Limiter, const real MinMod_Coeff,
                                      const real EP_Coeff, const real dt, const real dh,
                                      const real MinDens, const real MinPres,
                                      const bool NormPassive, const int NNorm, const int NormIdx[] );
extern void Hydro_Con2Pri_AllPatch( const real ConVar[][ FLU_NXT*FLU_NXT*FLU_NXT ],
                                          real PriVar[][ FLU_NXT*FLU_NXT*FLU_NXT ],
                                    const real Gamma_m1, const real MinPres,
                                    const bool NormPassive, const int NNorm, const int NormIdx[],
                                    const bool JeansMinPres, const real JeansMinPres_Coeff );
extern void Hydro_ComputeFlux( const real FC_Var [][NCOMP_TOTAL][ N_FC_VAR *N_FC_VAR *N_FC_VAR  ],
                                     real FC_Flux[][NCOMP_TOTAL][ N_FC_FLUX*N_FC_FLUX*N_FC_FLUX ],
                               const int Gap, const real Gamma, const bool CorrHalfVel, const real Pot_USG[],
                               const double Corner[], const real dt, const real dh, const double Time,
                               const OptGravityType_t GravityType, const double ExtAcc_AuxArray[], const real MinPres,
                               const bool DumpIntFlux, real IntFlux[][NCOMP_TOTAL][ PS2*PS2 ] );
extern void Hydro_FullStepUpdate( const real Input[][ FLU_NXT*FLU_NXT*FLU_NXT ], real Output[][ PS2*PS2*PS2 ], char DE_Status[],
                                  const real Flux[][NCOMP_TOTAL][ N_FC_FLUX*N_FC_FLUX*N_FC_FLUX ], const real dt, const real dh,
                                  const real Gamma_m1, const real _Gamma_m1, const real MinDens, const real MinPres, const real DualEnergySwitch,
                                  const bool NormPassive, const int NNorm, const int NormIdx[] );
#if   ( RSOLVER == EXACT )
extern void Hydro_RiemannSolver_Exact( const int XYZ, real eival_out[], real L_star_out[], real R_star_out[],
                                       real Flux_Out[], const real L_In[], const real R_In[], const real Gamma );
#elif ( RSOLVER == ROE )
extern void Hydro_RiemannSolver_Roe( const int XYZ, real Flux_Out[], const real L_In[], const real R_In[],
                                     const real Gamma, const real MinPres );
#elif ( RSOLVER == HLLE )
extern void Hydro_RiemannSolver_HLLE( const int XYZ, real Flux_Out[], const real L_In[], const real R_In[],
                                      const real Gamma, const real MinPres );
#elif ( RSOLVER == HLLC )
extern void Hydro_RiemannSolver_HLLC( const int XYZ, real Flux_Out[], const real L_In[], const real R_In[],
                                      const real Gamma, const real MinPres );
#endif
extern real Hydro_CheckMinPres( const real InPres, const real MinPres );

#if   ( FLU_SCHEME == MHM_RP )
extern void Hydro_Con2Pri( const real In[], real Out[], const real Gamma_m1, const real MinPres,
                           const bool NormPassive, const int NNorm, const int NormIdx[],
                           const bool JeansMinPres, const real JeansMinPres_Coeff );
#endif

#endif // #ifdef __CUDACC__ ... else ...


// internal functions
#if   ( FLU_SCHEME == MHM_RP )
static void Hydro_RiemannPredict( const real Flu_Array_In[][ FLU_NXT*FLU_NXT*FLU_NXT ],
                                  const real Half_Flux[][3][NCOMP_TOTAL], real Half_Var[][NCOMP_TOTAL], const real dt,
                                const real dh, const real Gamma, const real MinDens, const real MinPres );
static void Hydro_RiemannPredict_Flux( const real Flu_Array_In[][ FLU_NXT*FLU_NXT*FLU_NXT ], real Half_Flux[][3][NCOMP_TOTAL],
                                       const real Gamma, const real MinPres );
#endif



#ifdef __CUDACC__

#ifdef UNSPLIT_GRAVITY
#include "CUPOT.h"
__constant__ double ExtAcc_AuxArray_d_Flu[EXT_ACC_NAUX_MAX];

//-------------------------------------------------------------------------------------------------------
// Function    :  CUFLU_FluidSolver_SetConstMem_ExtAcc
// Description :  Set the constant memory of ExtAcc_AuxArray_d_Flu used by CUFLU_FluidSolver_CTU/MHM
//
// Note        :  Adopt the suggested approach for CUDA version >= 5.0
//
// Parameter   :  None
//
// Return      :  0/-1 : successful/failed
//---------------------------------------------------------------------------------------------------
int CUFLU_FluidSolver_SetConstMem_ExtAcc( double ExtAcc_AuxArray_h[] )
{

   if (  cudaSuccess != cudaMemcpyToSymbol( ExtAcc_AuxArray_d_Flu, ExtAcc_AuxArray_h, EXT_ACC_NAUX_MAX*sizeof(double),
                                            0, cudaMemcpyHostToDevice)  )
      return -1;

   else
      return 0;

} // FUNCTION : CUFLU_FluidSolver_SetConstMem_ExtAcc
#endif // #ifdef UNSPLIT_GRAVITY


#if ( NCOMP_PASSIVE > 0 )
__constant__ int NormIdx[NCOMP_PASSIVE];

//-------------------------------------------------------------------------------------------------------
// Function    :  CUFLU_FluidSolver_SetConstMem_NormIdx
// Description :  Set the constant memory of NormIdx[] used by CUFLU_FluidSolver_CTU/MHM
//
// Note        :  Adopt the suggested approach for CUDA version >= 5.0
//
// Parameter   :  None
//
// Return      :  0/-1 : successful/failed
//---------------------------------------------------------------------------------------------------
int CUFLU_FluidSolver_SetConstMem_NormIdx( int NormIdx_h[] )
{

   if (  cudaSuccess != cudaMemcpyToSymbol( NormIdx, NormIdx_h, NCOMP_PASSIVE*sizeof(int),
                                            0, cudaMemcpyHostToDevice)  )
      return -1;

   else
      return 0;

} // FUNCTION : CUFLU_FluidSolver_SetConstMem_NormIdx

#else
__constant__ int *NormIdx = NULL;

#endif // #if ( NCOMP_PASSIVE > 0 ) ... else ...

#endif // #ifdef __CUDACC__




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
//
// Parameter   :  Flu_Array_In       : Array storing the input fluid variables
//                Flu_Array_Out      : Array to store the output fluid variables
//                DE_Array_Out       : Array to store the dual-energy status
//                Flux_Array         : Array to store the output fluxes
//                Corner_Array       : Array storing the physical corner coordinates of each patch group (for UNSPLIT_GRAVITY)
//                Pot_Array_USG      : Array storing the input potential for UNSPLIT_GRAVITY
//                NPatchGroup        : Number of patch groups to be evaluated
//                dt                 : Time interval to advance solution
//                dh                 : Grid size
//                Gamma              : Ratio of specific heats
//                StoreFlux          : true --> store the coarse-fine fluxes
//                LR_Limiter         : Slope limiter for the data reconstruction in the MHM/MHM_RP/CTU schemes
//                                     (0/1/2/3/4) = (vanLeer/generalized MinMod/vanAlbada/
//                                                    vanLeer + generalized MinMod/extrema-preserving) limiter
//                MinMod_Coeff       : Coefficient of the generalized MinMod limiter
//                EP_Coeff           : Coefficient of the extrema-preserving limiter
//                Time               : Current physical time                                     (for UNSPLIT_GRAVITY only)
//                GravityType        : Types of gravity --> self-gravity, external gravity, both (for UNSPLIT_GRAVITY only)
//                ExtAcc_AuxArray    : Auxiliary array for adding external acceleration          (for UNSPLIT_GRAVITY only)
//                MinDens/Pres       : Minimum allowed density and pressure
//                DualEnergySwitch   : Use the dual-energy formalism if E_int/E_kin < DualEnergySwitch
//                NormPassive        : true --> normalize passive scalars so that the sum of their mass density
//                                              is equal to the gas mass density
//                NNorm              : Number of passive scalars to be normalized
//                                     --> Should be set to the global variable "PassiveNorm_NVar"
//                NormIdx            : Target variable indices to be normalized
//                                     --> Should be set to the global variable "PassiveNorm_VarIdx"
//                JeansMinPres       : Apply minimum pressure estimated from the Jeans length
//                JeansMinPres_Coeff : Coefficient used by JeansMinPres = G*(Jeans_NCell*Jeans_dh)^2/(Gamma*pi);
//-------------------------------------------------------------------------------------------------------
#ifdef __CUDACC__
__global__ void CUFLU_FluidSolver_MHM(
   const real   Flu_Array_In []   [NCOMP_TOTAL][ FLU_NXT*FLU_NXT*FLU_NXT ],
         real   Flu_Array_Out[]   [NCOMP_TOTAL][ PS2*PS2*PS2 ],
         char   DE_Array_Out []                [ PS2*PS2*PS2 ],
         real   Flux_Array   [][9][NCOMP_TOTAL][ PS2*PS2 ],
   const double Corner_Array [][3],
   const real   Pot_Array_USG[]                [ USG_NXT_F*USG_NXT_F*USG_NXT_F ],
   real PriVar   []   [NCOMP_TOTAL][ FLU_NXT*FLU_NXT*FLU_NXT ],
   real Slope_PPM[][3][NCOMP_TOTAL][ N_SLOPE_PPM*N_SLOPE_PPM*N_SLOPE_PPM],
   real FC_Var   [][6][NCOMP_TOTAL][ N_FC_VAR*N_FC_VAR*N_FC_VAR ],
   real FC_Flux  [][3][NCOMP_TOTAL][ N_FC_FLUX*N_FC_FLUX*N_FC_FLUX ],
   const real dt, const real dh, const real Gamma, const bool StoreFlux,
   const LR_Limiter_t LR_Limiter, const real MinMod_Coeff,
   const real EP_Coeff, const double Time, const OptGravityType_t GravityType,
   const real MinDens, const real MinPres, const real DualEnergySwitch,
   const bool NormPassive, const int NNorm,
   const bool JeansMinPres, const real JeansMinPres_Coeff )
#else
void CPU_FluidSolver_MHM(
   const real   Flu_Array_In []   [NCOMP_TOTAL][ FLU_NXT*FLU_NXT*FLU_NXT ],
         real   Flu_Array_Out[]   [NCOMP_TOTAL][ PS2*PS2*PS2 ],
         char   DE_Array_Out []                [ PS2*PS2*PS2 ],
         real   Flux_Array   [][9][NCOMP_TOTAL][ PS2*PS2 ],
   const double Corner_Array [][3],
   const real   Pot_Array_USG[]                [ USG_NXT_F*USG_NXT_F*USG_NXT_F ],
   const int NPatchGroup, const real dt, const real dh, const real Gamma,
   const bool StoreFlux, const LR_Limiter_t LR_Limiter, const real MinMod_Coeff,
   const real EP_Coeff, const double Time, const OptGravityType_t GravityType,
   const double ExtAcc_AuxArray[], const real MinDens, const real MinPres,
   const real DualEnergySwitch, const bool NormPassive, const int NNorm, const int NormIdx[],
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


// openmp pragma for the CPU solver
#  ifndef __CUDACC__
#  pragma omp parallel
#  endif
   {
//    FC: Face-Centered variables/fluxes
//    --> CPU solver: allocate data
//        GPU solver: link to the input global memory array
#     ifdef __CUDACC__
      real (*FC_Var_1PG )[NCOMP_TOTAL][ N_FC_VAR*N_FC_VAR*N_FC_VAR    ] = FC_Var [blockIdx.x];
      real (*FC_Flux_1PG)[NCOMP_TOTAL][ N_FC_FLUX*N_FC_FLUX*N_FC_FLUX ] = FC_Flux[blockIdx.x];
      real (*PriVar_1PG )             [ FLU_NXT*FLU_NXT*FLU_NXT       ] = PriVar [blockIdx.x];
#     else
      real (*FC_Var_1PG )[NCOMP_TOTAL][ N_FC_VAR*N_FC_VAR*N_FC_VAR    ] = new real [6][NCOMP_TOTAL][ N_FC_VAR*N_FC_VAR*N_FC_VAR    ];
      real (*FC_Flux_1PG)[NCOMP_TOTAL][ N_FC_FLUX*N_FC_FLUX*N_FC_FLUX ] = new real [3][NCOMP_TOTAL][ N_FC_FLUX*N_FC_FLUX*N_FC_FLUX ];   // also used by "Half_Flux"
      real (*PriVar_1PG )             [ FLU_NXT*FLU_NXT*FLU_NXT       ] = new real    [NCOMP_TOTAL][ FLU_NXT*FLU_NXT*FLU_NXT       ];   // also used by "Half_Var"

#     if ( FLU_SCHEME == MHM_RP )
      real (*const Half_Flux_1PG)[NCOMP_TOTAL][ N_FC_FLUX*N_FC_FLUX*N_FC_FLUX ] = FC_Flux_1PG;
      real (*const Half_Var_1PG)              [ N_FC_FLUX*N_FC_FLUX*N_FC_FLUX ] = PriVar_1PG;
#     endif
#     endif // #ifdef __CUDACC__ ... else ...


//    loop over all patch groups
//    --> CPU solver: use different OpenMP threads    to work on different patch groups
//        GPU solver: ...           CUDA thread block ...
#     ifdef __CUDACC__
      const int P = blockIdx.x;
#     else
#     pragma omp for schedule( runtime )
      for (int P=0; P<NPatchGroup; P++)
#     endif
      {

//       1. half-step prediction
//       (1-a) MHM_RP: use Riemann solver to calculate the half-step fluxes
#        if ( FLU_SCHEME == MHM_RP )

         /*
//       (1-a-1) evaluate the half-step first-order fluxes by Riemann solver
         Hydro_RiemannPredict_Flux( Flu_Array_In[P], Half_Flux, Gamma, MinPres );


//       (1-a-2) evaluate the half-step solutions
         Hydro_RiemannPredict( Flu_Array_In[P], Half_Flux, Half_Var, dt, dh, Gamma, MinDens, MinPres );


//       (1-a-3) conserved variables --> primitive variables
         for (int k=0; k<N_HF_VAR; k++)
         for (int j=0; j<N_HF_VAR; j++)
         for (int i=0; i<N_HF_VAR; i++)
         {
            ID1 = (k*N_HF_VAR + j)*N_HF_VAR + i;

            for (int v=0; v<NCOMP_TOTAL; v++)   Input[v] = Half_Var[ID1][v];

            Hydro_Con2Pri( Input, Half_Var[ID1], Gamma_m1, MinPres, NormPassive, NNorm, NormIdx,
                           JeansMinPres, JeansMinPres_Coeff );
         }


//       (1-a-4) evaluate the face-centered values by data reconstruction
         Hydro_DataReconstruction( Half_Var, FC_Var, N_HF_VAR, FLU_GHOST_SIZE-2, Gamma, LR_Limiter,
                                   MinMod_Coeff, EP_Coeff, NULL_REAL, NULL_INT, MinDens, MinPres );
         */


//       (1-b) MHM: use interpolated face-centered values to calculate the half-step fluxes
#        elif ( FLU_SCHEME == MHM )

//       (1-b-1) conserved variables --> primitive variables
         Hydro_Con2Pri_AllPatch( Flu_Array_In[P], PriVar_1PG, Gamma_m1, MinPres, NormPassive, NNorm, NormIdx,
                                 JeansMinPres, JeansMinPres_Coeff );
#        ifdef __CUDACC__
         __syncthreads();
#        endif


//       (1-b-2) evaluate the face-centered values by data reconstruction
         Hydro_DataReconstruction( PriVar_1PG, Flu_Array_In[P], FC_Var_1PG, FLU_NXT, FLU_GHOST_SIZE-1,
                                   Gamma, LR_Limiter, MinMod_Coeff, EP_Coeff, NULL_REAL, NULL_INT,
                                   MinDens, MinPres, NormPassive, NNorm, NormIdx );
#        ifdef __CUDACC__
         __syncthreads();
#        endif

#        endif // #if ( FLU_SCHEME == MHM_RP ) ... else ...


//       2. evaluate the full-step fluxes
#        ifdef UNSPLIT_GRAVITY
         Hydro_ComputeFlux( FC_Var_1PG, FC_Flux_1PG, 1, Gamma, CorrHalfVel_Yes,
                            Pot_Array_USG[P], Corner_Array[P],
                            dt, dh, Time, GravityType, ExtAcc_AuxArray, MinPres,
                            StoreFlux, Flux_Array[P] );
#        else
         Hydro_ComputeFlux( FC_Var_1PG, FC_Flux_1PG, 1, Gamma, CorrHalfVel_No,
                            NULL, NULL,
                            NULL_REAL, NULL_REAL, NULL_REAL, GRAVITY_NONE, NULL, MinPres,
                            StoreFlux, Flux_Array[P] );
#        endif

#        ifdef __CUDACC__
         __syncthreads();
#        endif


//       3. full-step evolution
         Hydro_FullStepUpdate( Flu_Array_In[P], Flu_Array_Out[P], DE_Array_Out[P],
                               FC_Flux_1PG, dt, dh, Gamma_m1, _Gamma_m1, MinDens, MinPres, DualEnergySwitch,
                               NormPassive, NNorm, NormIdx );

      } // loop over all patch groups

#     ifndef __CUDACC__
      delete [] FC_Var_1PG;
      delete [] FC_Flux_1PG;
      delete [] PriVar_1PG;
#     endif

   } // OpenMP parallel region

} // FUNCTION : CPU_FluidSolver_MHM



#if ( FLU_SCHEME == MHM_RP )
//-------------------------------------------------------------------------------------------------------
// Function    :  Hydro_RiemannPredict_Flux
// Description :  Evaluate the half-step face-centered fluxes by Riemann solver
//
// Note        :  1. Work for the MUSCL-Hancock method + Riemann-prediction (MHM_RP)
//                2. Currently support the exact, Roe, HLLE, and HLLC solvers
//
// Parameter   :  Flu_Array_In : Array storing the input conserved variables
//                Half_Flux    : Array to store the output face-centered fluxes
//                               --> The size is assumed to be N_HF_FLUX^3
//                Gamma        : Ratio of specific heats
//                MinPres      : Minimum allowed pressure
//-------------------------------------------------------------------------------------------------------
void Hydro_RiemannPredict_Flux( const real Flu_Array_In[][ FLU_NXT*FLU_NXT*FLU_NXT ], real Half_Flux[][3][NCOMP_TOTAL],
                                const real Gamma, const real MinPres )
{

   const int dr[3] = { 1, FLU_NXT, FLU_NXT*FLU_NXT };
   int ID1, ID2, dN[3]={ 0 };
   real ConVar_L[NCOMP_TOTAL], ConVar_R[NCOMP_TOTAL];

#  if ( RSOLVER == EXACT )
   const real Gamma_m1 = Gamma - (real)1.0;
   real PriVar_L[NCOMP_TOTAL], PriVar_R[NCOMP_TOTAL];
#  endif


// loop over different spatial directions
   for (int d=0; d<3; d++)
   {
      switch ( d )
      {
         case 0 : dN[0] = 0;  dN[1] = 1;  dN[2] = 1;  break;
         case 1 : dN[0] = 1;  dN[1] = 0;  dN[2] = 1;  break;
         case 2 : dN[0] = 1;  dN[1] = 1;  dN[2] = 0;  break;
      }

      for (int k1=0, k2=dN[2];  k1<N_HF_FLUX-dN[2];  k1++, k2++)
      for (int j1=0, j2=dN[1];  j1<N_HF_FLUX-dN[1];  j1++, j2++)
      for (int i1=0, i2=dN[0];  i1<N_HF_FLUX-dN[0];  i1++, i2++)
      {
         ID1 = (k1*N_HF_FLUX + j1)*N_HF_FLUX + i1;
         ID2 = (k2*FLU_NXT   + j2)*FLU_NXT   + i2;

//       get the left and right states
         for (int v=0; v<NCOMP_TOTAL; v++)
         {
            ConVar_L[v] = Flu_Array_In[v][ ID2       ];
            ConVar_R[v] = Flu_Array_In[v][ ID2+dr[d] ];
         }

//       invoke the Riemann solver
#        if   ( RSOLVER == EXACT )
         const bool NormPassive_No  = false;  // do NOT convert any passive variable to mass fraction for the Riemann solvers
         const bool JeansMinPres_No = false;

         Hydro_Con2Pri( ConVar_L, PriVar_L, Gamma_m1, MinPres, NormPassive_No, NULL_INT, NULL, JeansMinPres_No, NULL_REAL );
         Hydro_Con2Pri( ConVar_R, PriVar_R, Gamma_m1, MinPres, NormPassive_No, NULL_INT, NULL, JeansMinPres_No, NULL_REAL );

         Hydro_RiemannSolver_Exact( d, NULL, NULL, NULL, Half_Flux[ID1][d], PriVar_L, PriVar_R, Gamma );
#        elif ( RSOLVER == ROE )
         Hydro_RiemannSolver_Roe ( d, Half_Flux[ID1][d], ConVar_L, ConVar_R, Gamma, MinPres );
#        elif ( RSOLVER == HLLE )
         Hydro_RiemannSolver_HLLE( d, Half_Flux[ID1][d], ConVar_L, ConVar_R, Gamma, MinPres );
#        elif ( RSOLVER == HLLC )
         Hydro_RiemannSolver_HLLC( d, Half_Flux[ID1][d], ConVar_L, ConVar_R, Gamma, MinPres );
#        else
#        error : ERROR : unsupported Riemann solver (EXACT/ROE) !!
#        endif
      }
   } // for (int d=0; d<3; d++)

} // FUNCTION : Hydro_RiemannPredict_Flux



//-------------------------------------------------------------------------------------------------------
// Function    :  Hydro_RiemannPredict
// Description :  Evolve the cell-centered variables by half time-step by using the Riemann solvers
//
// Note        :  Work for the MUSCL-Hancock method + Riemann-prediction (MHM_RP)
//
// Parameter   :  Flu_Array_In : Array storing the input conserved variables
//                Half_Flux    : Array storing the input face-centered fluxes
//                               --> The size is assumed to be N_HF_FLUX^3
//                Half_Var     : Array to store the output conserved variables
//                               --> The size is assumed to be N_HF_VAR^3
//                dt           : Time interval to advance solution
//                dh           : Grid size
//                Gamma        : Ratio of specific heats
//                MinDens/Pres : Minimum allowed density and pressure
//-------------------------------------------------------------------------------------------------------
void Hydro_RiemannPredict( const real Flu_Array_In[][ FLU_NXT*FLU_NXT*FLU_NXT ], const real Half_Flux[][3][NCOMP_TOTAL],
                           real Half_Var[][NCOMP_TOTAL], const real dt, const real dh, const real Gamma,
                           const real MinDens, const real MinPres )
{

   const int  dID3[3]   = { 1, N_HF_FLUX, N_HF_FLUX*N_HF_FLUX };
   const real dt_dh2    = (real)0.5*dt/dh;
   const real  Gamma_m1 = Gamma - (real)1.0;
   const real _Gamma_m1 = (real)1.0 / Gamma_m1;

   real dF[3][NCOMP_TOTAL];
   int ID1, ID2, ID3;


   for (int k1=0, k2=1;  k1<N_HF_VAR;  k1++, k2++)
   for (int j1=0, j2=1;  j1<N_HF_VAR;  j1++, j2++)
   for (int i1=0, i2=1;  i1<N_HF_VAR;  i1++, i2++)
   {
      ID1 = (k1*N_HF_VAR  + j1)*N_HF_VAR  + i1;
      ID2 = (k2*FLU_NXT   + j2)*FLU_NXT   + i2;
      ID3 = (k1*N_HF_FLUX + j1)*N_HF_FLUX + i1;

      for (int d=0; d<3; d++)
      for (int v=0; v<NCOMP_TOTAL; v++)    dF[d][v] = Half_Flux[ ID3+dID3[d] ][d][v] - Half_Flux[ID3][d][v];

      for (int v=0; v<NCOMP_TOTAL; v++)
         Half_Var[ID1][v] = Flu_Array_In[v][ID2] - dt_dh2*( dF[0][v] + dF[1][v] + dF[2][v] );

//    ensure positive density and pressure
      Half_Var[ID1][0] = FMAX( Half_Var[ID1][0], MinDens );
      Half_Var[ID1][4] = Hydro_CheckMinPresInEngy( Half_Var[ID1][0], Half_Var[ID1][1], Half_Var[ID1][2],
                                                   Half_Var[ID1][3], Half_Var[ID1][4], Gamma_m1, _Gamma_m1, MinPres );
#     if ( NCOMP_PASSIVE > 0 )
      for (int v=NCOMP_FLUID; v<NCOMP_TOTAL; v++)
      Half_Var[ID1][v] = FMAX( Half_Var[ID1][v], TINY_NUMBER );
#     endif
   } // i,j,k

} // FUNCTION : Hydro_RiemannPredict
#endif // #if ( FLU_SCHEME == MHM_RP )



#endif // #if (  MODEL == HYDRO  &&  ( FLU_SCHEME == MHM || FLU_SCHEME == MHM_RP )  )
