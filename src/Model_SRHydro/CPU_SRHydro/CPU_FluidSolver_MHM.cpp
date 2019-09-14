#include "CUFLU.h"

#if (  MODEL == SR_HYDRO  &&  ( FLU_SCHEME == MHM || FLU_SCHEME == MHM_RP )  )



// external functions
#ifdef __CUDACC__

#include "CUFLU_Shared_FluUtility.cu"
#include "CUFLU_Shared_DataReconstruction.cu"
#include "CUFLU_Shared_ComputeFlux.cu"
#include "CUFLU_Shared_FullStepUpdate.cu"

#if ( RSOLVER == HLLE )
# include "CUFLU_Shared_RiemannSolver_HLLE.cu"
#elif ( RSOLVER == HLLC )
# include "CUFLU_Shared_RiemannSolver_HLLC.cu"
#endif

#include "CUFLU_SetConstMem_FluidSolver.cu"

#else // #ifdef __CUDACC__
# include "../../../include/SRHydroPrototypes.h"

#endif // #ifdef __CUDACC__ ... else ...


// internal functions
#if ( FLU_SCHEME == MHM_RP )
GPU_DEVICE
static void SRHydro_RiemannPredict_Flux( const real g_ConVar[][ CUBE(FLU_NXT) ],
                                               real g_Half_Flux[][NCOMP_TOTAL][ CUBE(N_FC_FLUX) ],
                                         const real Gamma, const real MinTemp );
GPU_DEVICE
static void SRHydro_RiemannPredict( const real g_ConVar_In[][ CUBE(FLU_NXT) ],
                                    const real g_Half_Flux[][NCOMP_TOTAL][ CUBE(N_FC_FLUX) ],
                                          real g_Half_Var [][ CUBE(FLU_NXT) ],
                                    const real dt, const real dh, const real Gamma, const real MinDens, const real MinTemp );
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
// Parameter   :  [ 1] g_Flu_Array_In     : Array storing the input fluid variables
//                [ 2] g_Flu_Array_Out    : Array to store the output fluid variables
//                [ 3] g_Flux_Array       : Array to store the output fluxes
//                [ 4] g_PriVar           : Array to store the primitive variables
//                [ 5] g_Slope_PPM        : Array to store the slope for the PPM reconstruction
//                [ 6] g_FC_Var           : Array to store the half-step variables
//                [ 7] g_FC_Flux          : Array to store the face-centered fluxes
//                [ 8] NPatchGroup        : Number of patch groups to be evaluated
//                [ 9] dt                 : Time interval to advance solution
//                [10] dh                 : Cell size
//                [11] Gamma              : Ratio of specific heats
//                [12] StoreFlux          : true --> store the coarse-fine fluxes
//                [13] LR_Limiter         : Slope limiter for the data reconstruction in the MHM/MHM_RP/CTU schemes
//                                          (0/1/2/3/4) = (vanLeer/generalized MinMod/vanAlbada/
//                                                         vanLeer + generalized MinMod/extrema-preserving) limiter
//                [14] MinMod_Coeff       : Coefficient of the generalized MinMod limiter
//             [15/16] MinDens/Temp       : Minimum allowed density and temperature
//-------------------------------------------------------------------------------------------------------
#ifdef __CUDACC__
__global__
void CUFLU_FluidSolver_MHM(
   const real   g_Flu_Array_In [][NCOMP_TOTAL][ CUBE(FLU_NXT) ],
         real   g_Flu_Array_Out[][NCOMP_TOTAL][ CUBE(PS2) ],
         real   g_Flux_Array   [][9][NCOMP_TOTAL][ SQR(PS2) ],
         real   g_PriVar       [][NCOMP_TOTAL][ CUBE(FLU_NXT) ],
         real   g_Slope_PPM    [][3][NCOMP_TOTAL][ CUBE(N_SLOPE_PPM) ],
         real   g_FC_Var       [][6][NCOMP_TOTAL][ CUBE(N_FC_VAR) ],
         real   g_FC_Flux      [][3][NCOMP_TOTAL][ CUBE(N_FC_FLUX) ],
   const real dt, const real dh, const real Gamma, const bool StoreFlux,
   const LR_Limiter_t LR_Limiter, const real MinMod_Coeff,
   const real MinDens, const real MinTemp )
#else
void CPU_FluidSolver_MHM(
   const real   g_Flu_Array_In [][NCOMP_TOTAL][ CUBE(FLU_NXT) ],
         real   g_Flu_Array_Out[][NCOMP_TOTAL][ CUBE(PS2) ],
         real   g_Flux_Array   [][9][NCOMP_TOTAL][ SQR(PS2) ],
         real   g_PriVar       [][NCOMP_TOTAL][ CUBE(FLU_NXT) ],
         real   g_Slope_PPM    [][3][NCOMP_TOTAL][ CUBE(N_SLOPE_PPM) ],
         real   g_FC_Var       [][6][NCOMP_TOTAL][ CUBE(N_FC_VAR) ],
         real   g_FC_Flux      [][3][NCOMP_TOTAL][ CUBE(N_FC_FLUX) ],
   const int NPatchGroup, const real dt, const real dh, const real Gamma,
   const bool StoreFlux, const LR_Limiter_t LR_Limiter, const real MinMod_Coeff,
   const real MinDens, const real MinTemp )
#endif // #ifdef __CUDACC__ ... else ...
{

   const char Max = 4;
   char iteration;
   real AdaptiveMinModCoeff;

#  ifdef __CUDACC__
   __shared__ char state;
#  else
   char state;
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

      real (*const g_FC_Var_1PG   )[NCOMP_TOTAL][ CUBE(N_FC_VAR)    ] = g_FC_Var   [array_idx];
      real (*const g_FC_Flux_1PG  )[NCOMP_TOTAL][ CUBE(N_FC_FLUX)   ] = g_FC_Flux  [array_idx];
      real (*const g_PriVar_1PG   )             [ CUBE(FLU_NXT)     ] = g_PriVar   [array_idx];
      real (*const g_Slope_PPM_1PG)[NCOMP_TOTAL][ CUBE(N_SLOPE_PPM) ] = g_Slope_PPM[array_idx];

#     if ( FLU_SCHEME == MHM_RP )
      real (*const g_Half_Flux_1PG)[NCOMP_TOTAL][ CUBE(N_FC_FLUX) ] = g_FC_Flux_1PG;
      real (*const g_Half_Var_1PG )             [ CUBE(FLU_NXT)   ] = g_PriVar_1PG;
#     endif


//    loop over all patch groups
//    --> CPU/GPU solver: use different (OpenMP threads) / (CUDA thread blocks)
//        to work on different patch groups
#     ifdef __CUDACC__
      const int P = blockIdx.x;
#     else
#     pragma omp for schedule( runtime ) private( iteration, AdaptiveMinModCoeff, state )
      for (int P=0; P<NPatchGroup; P++)
#     endif
      {
         iteration = 0;

//       1. half-step prediction
//       1-a. MHM_RP: use Riemann solver to calculate the half-step fluxes
#        if ( FLU_SCHEME == MHM_RP )

//       1-a-1. evaluate the half-step first-order fluxes by Riemann solver
//              --> check unphysical cells in g_Flu_Array_In[] before computing flux
         SRHydro_RiemannPredict_Flux( g_Flu_Array_In[P], g_Half_Flux_1PG, Gamma, MinTemp );


//       1-a-2. evaluate the half-step solutions
//              --> check unphysical cells in g_Half_Var_1PG[] after prediction
         SRHydro_RiemannPredict( g_Flu_Array_In[P], g_Half_Flux_1PG, g_Half_Var_1PG, dt, dh, Gamma, MinDens, MinTemp );

         do {
               state = 0;

//             adaptive minmod coefficient
               AdaptiveMinModCoeff = ( Max - iteration ) * ( MinMod_Coeff / (real) Max );


//             1-a-3. evaluate the face-centered values by data reconstruction
//                    --> note that g_Half_Var_1PG[] returned by SRHydro_RiemannPredict() stores the primitive variables
               SRHydro_DataReconstruction( NULL, g_Half_Var_1PG, g_FC_Var_1PG, g_Slope_PPM_1PG,
                                           N_HF_VAR, FLU_GHOST_SIZE-2,
                                           Gamma, LR_Limiter, AdaptiveMinModCoeff, dt, dh, MinDens, MinTemp );


//       1-b. MHM: use interpolated face-centered values to calculate the half-step fluxes
#        elif ( FLU_SCHEME == MHM )

         do {
               state = 0;

//             adaptive minmod coefficient         
               AdaptiveMinModCoeff = ( Max - iteration ) * ( MinMod_Coeff / (real) Max );


//             evaluate the face-centered values by data reconstruction
//             --> check unphysical cells in g_Flu_Array_In[] before data reconstruction
               SRHydro_DataReconstruction( g_Flu_Array_In[P], g_PriVar_1PG, g_FC_Var_1PG, g_Slope_PPM_1PG,
                                           FLU_NXT, FLU_GHOST_SIZE-1,
                                           Gamma, LR_Limiter, AdaptiveMinModCoeff, dt, dh, MinDens, MinTemp );
#        endif // #if ( FLU_SCHEME == MHM_RP ) ... else ...



//             2. evaluate the full-step fluxes
//                --> check unphysical cells in g_FC_Var_1PG[] before computing flux
               SRHydro_ComputeFlux( g_FC_Var_1PG, g_FC_Flux_1PG, 1, Gamma,
                                    MinTemp, StoreFlux, g_Flux_Array[P] );


//             3. full-step evolution
//                --> check unphysical cells in g_Flu_Array_Out[] after full update
               SRHydro_FullStepUpdate( g_Flu_Array_In[P], g_Flu_Array_Out[P], NULL,
                                       g_FC_Flux_1PG, dt, dh, Gamma, MinDens, MinTemp, &state );

#              ifdef __CUDACC__
               if ( threadIdx.x == 0 && state == 1 )
#              else
               if ( state == 1 ) 
#              endif
                 printf("iteration=%d, AdaptiveMinModCoeff=%13.10f\n", iteration, AdaptiveMinModCoeff );

               iteration++;


            } while( state && iteration <= Max );

      } // loop over all patch groups
   } // OpenMP parallel region

} // FUNCTION : CPU_FluidSolver_MHM



#if ( FLU_SCHEME == MHM_RP )
//-------------------------------------------------------------------------------------------------------
// Function    :  SRHydro_RiemannPredict_Flux
// Description :  Evaluate the half-step face-centered fluxes by Riemann solver
//
// Note        :  1. Work for the MHM_RP scheme
//                2. Currently support the exact, Roe, HLLE, and HLLC solvers
//                3. g_Half_Flux[] is accessed with the stride N_FC_FLUX
//                   --> Fluxes on the **left** face of the (i+1,j+1,k+1) element in g_ConVar[] will
//                       be stored in the (i,j,k) element of g_Half_Flux[]
//
// Parameter   :  [1] g_ConVar    : Array storing the input conserved variables
//                [2] g_Half_Flux : Array to store the output face-centered fluxes
//                [3] Gamma       : Ratio of specific heats
//                [4] MinTemp     : Minimum allowed pressure
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE
void SRHydro_RiemannPredict_Flux( const real g_ConVar[][ CUBE(FLU_NXT) ],
                                        real g_Half_Flux[][NCOMP_TOTAL][ CUBE(N_FC_FLUX) ],
                                  const real Gamma, const real MinTemp )
{
   const int didx_cvar[3] = { 1, FLU_NXT, SQR(FLU_NXT) };
   real ConVar_L[NCOMP_TOTAL], ConVar_R[NCOMP_TOTAL], Flux_1Face[NCOMP_TOTAL];


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
            ConVar_L[v] = g_ConVar[v][ idx_cvar              ];
            ConVar_R[v] = g_ConVar[v][ idx_cvar+didx_cvar[d] ];
         }

//       check unphysical cells before computing flux
#        ifdef CHECK_NEGATIVE_IN_FLUID
         SRHydro_CheckUnphysical(ConVar_L, NULL, Gamma, MinTemp, __FUNCTION__, __LINE__, true);
         SRHydro_CheckUnphysical(ConVar_R, NULL, Gamma, MinTemp, __FUNCTION__, __LINE__, true);
#        endif

#        ifdef CHECK_MIN_TEMP
         ConVar_L[ENGY] = SRHydro_CheckMinTempInEngy(ConVar_L, MinTemp, Gamma);
         ConVar_R[ENGY] = SRHydro_CheckMinTempInEngy(ConVar_R, MinTemp, Gamma);
#        endif

//       invoke the Riemann solver
#        if ( RSOLVER == HLLE )
         SRHydro_RiemannSolver_HLLE ( d, Flux_1Face, ConVar_L, ConVar_R, Gamma, MinTemp );
#        elif ( RSOLVER == HLLC )
         SRHydro_RiemannSolver_HLLC ( d, Flux_1Face, ConVar_L, ConVar_R, Gamma, MinTemp );
#        else
#        error : ERROR : unsupported Riemann solver !!
#        endif

//       store the results in g_Half_Flux[]
         for (int v=0; v<NCOMP_TOTAL; v++)   g_Half_Flux[d][v][idx_flux] = Flux_1Face[v];
      } // CGPU_LOOP( idx, N_FC_FLUX*SQR(N_FC_FLUX-1) )
   } // for (int d=0; d<3; d++)


#  ifdef __CUDACC__
   __syncthreads();
#  endif

} // FUNCTION : SRHydro_RiemannPredict_Flux



//-------------------------------------------------------------------------------------------------------
// Function    :  SRHydro_RiemannPredict
// Description :  Evolve the cell-centered variables by half time-step using the fluxes returned
//                by SRHydro_RiemannPredict_Flux()
//
// Note        :  1. Work for the MHM_RP scheme
//                2. For the performance consideration, the output data are converted to primitive variables
//                   --> Reducing the global memory access on GPU
//
// Parameter   :  [1] g_ConVar_In        : Array storing the input conserved variables
//                [2] g_Half_Flux        : Array storing the input face-centered fluxes
//                                         --> Accessed with the stride N_FC_FLUX
//                [3] g_Half_Var         : Array to store the output primitive variables
//                                         --> Accessed with the stride N_HF_VAR
//                                         --> Although its actually allocated size is FLU_NXT^3 since it points to g_PriVar_1PG[]
//                [4] dt                 : Time interval to advance solution
//                [5] dh                 : Cell size
//                [6] Gamma              : Ratio of specific heats
//              [7/8]MinDens/Temp        : Minimum allowed density and temperature
//-------------------------------------------------------------------------------------------------------
GPU_DEVICE
void SRHydro_RiemannPredict( const real g_ConVar_In[][ CUBE(FLU_NXT) ],
                             const real g_Half_Flux[][NCOMP_TOTAL][ CUBE(N_FC_FLUX) ],
                                   real g_Half_Var [][ CUBE(FLU_NXT) ],
                             const real dt, const real dh, const real Gamma, const real MinDens, const real MinTemp )
{

   const int  didx_flux[3] = { 1, N_FC_FLUX, SQR(N_FC_FLUX) };
   const real dt_dh2       = (real)0.5*dt/dh;

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
      for (int v=0; v<NCOMP_TOTAL; v++)    dflux[d][v] = g_Half_Flux[d][v][ idx_flux+didx_flux[d] ] - g_Half_Flux[d][v][idx_flux];

//    update the input cell-centered conserved variables with the flux differences
      for (int v=0; v<NCOMP_TOTAL; v++)
         out_con[v] = g_ConVar_In[v][idx_in] - dt_dh2*( dflux[0][v] + dflux[1][v] + dflux[2][v] );

#     ifdef CHECK_NEGATIVE_IN_FLUID
      SRHydro_CheckUnphysical(out_con, NULL, Gamma, MinTemp, __FUNCTION__, __LINE__, true);
#     endif

#     ifdef CHECK_MIN_TEMP
      out_con[ENGY] = SRHydro_CheckMinTempInEngy( out_con, MinTemp, Gamma);
#     endif

//    conserved --> primitive variables
      SRHydro_Con2Pri( out_con, out_pri, Gamma, MinTemp );
	  SRHydro_3Velto4Vel( out_pri, out_pri );

//    store the results to g_Half_Var[]
      for (int v=0; v<NCOMP_TOTAL; v++)   g_Half_Var[v][idx_out] = out_pri[v];
   } // i,j,k


#  ifdef __CUDACC__
   __syncthreads();
#  endif

} // FUNCTION : SRHydro_RiemannPredict
#endif // #if ( FLU_SCHEME == MHM_RP )



#endif // #if (  MODEL == SR_HYDRO  &&  ( FLU_SCHEME == MHM || FLU_SCHEME == MHM_RP )  )
